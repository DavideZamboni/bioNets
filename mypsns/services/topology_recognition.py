from os import path, makedirs
import math
import pandas as pd
from pathlib import Path
import os
import MDAnalysis as mda
import timeit
import logging
import sys
script_path = path.abspath(__file__)
script_dir = path.dirname(script_path)
import warnings

warnings.filterwarnings(
    "ignore",
    message="Invalid atom serials were present, bonds will not be parsed",
    category=UserWarning,
    module="MDAnalysis.topology.PDBParser"
)

warnings.filterwarnings(
    "ignore",
     message="Serial numbers went over 100,000.  Higher serials have been guessed",
     category=UserWarning,
     module="MDAnalysis.topology.PDBParser"
 )

warnings.filterwarnings(
    "ignore",
    message="Unknown element  found for some atoms. These have been given an empty element record. If needed they can be guessed using universe.guess_TopologyAttrs(context='default', to_guess=['elements']).",
    category=UserWarning,
    module="MDAnalysis.topology.PDBParser"
)

sys.path.insert(0, script_dir + '/scripts')

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from mypsns.utils.topology_recognition_utilities import (
    associate_ff_params, get_residue_types, parse_initial_topology_info,
    COMPLETE_TABLE, FF_ATOM_COLUMNS, UNMATCHED_TABLE,
    candidate_files, candidate_index, dangling_bonds, pair_candidate_atoms,
    read_str_csv, read_topology_table, topology_key, with_structure_names,
)
from mypsns.utils.general_utilities import get_endtime
pd.options.mode.chained_assignment = None
log = logging.getLogger(__name__)


class TopologyAssociationError(RuntimeError):
    """Raised when topology recognition cannot associate every distinct residue
    topology in the input structure to a forcefield entry. Halts the pipeline
    before the interactors stage (which would otherwise fail with a confusing
    'file not found' for the unmatched residues). The fix is either a corrected
    forcefield / structure or `mypsns accept-topology` (accept_topologies)."""

class TopologyAcceptError(RuntimeError):
    """Raised by accept_topologies when a requested choice cannot be applied
    (unknown topology, missing candidate, atoms that cannot be paired). Nothing
    is written when it is raised."""

class Topology_recognition():
    def __init__(
            self, 
            topology: str,
            outpath: Path, 
            ff_path: Path,
            na_values: list = ["", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NULL", "NaN", "n/a", "nan", "null"]                ):
        '''
         - topology: topology file (pdb, gro, ...)
         - outpath: output folder
         - na_values: specifies strings recognized as pandas.NA values.
                      feature added to avoid problems with sodium

        Obtains useful information from input topology, producing .csv file with all the possible following info:

        - chainID: chain index
        - resindex: internal (unique and 0-based) residue index
        - resid: residue index information (from topology)
        - resname: residue name
        - index: internal (unique and 0-based) atom index
        - id: atom index information (from topology)
        - name: atom name
        - type: atom type (could be element or forcefield atomtype, depending on input topology file)
        - element: atomic element
        - mass: atomic mass 
        - charge: atom's partial charge 

        NW: obtained info depend on info defined in the topology file provided. For more information, visit https://userguide.mdanalysis.org/stable/topology_system.html#format-specific-attributes

        Defines:
         - self.u: reference Universe
         - self.out: output folder, exactly outpath
         - self.nas: na_values list
        '''

        self.u = mda.Universe(topology)

        # Internalize output folder, create log 
        self.out = Path(outpath)

        # Internalize ff folder
        self.ff = ff_path


        # 1 - Obtain information from current topology
        parse_initial_topology_info(self.u, self.out)
    
        if not hasattr(self, 'start'):
            self.start = timeit.default_timer()

        self.nas = na_values
        log.info('Topology_recognition run starting; initial topology info written to %s',
                 (self.out / 'initial_topology_info.csv').absolute())
      
    def get_restypes(self):
        '''
        Produces a .csv file containing a table containing all the possible topologies found for each of the residue names in universe. E.g.:

         resname    atnames             rixs
             ALA    ALA_atoms.csv       ALA_rixs.txt
             ARG    ARG_atoms.csv       ARG_rixs.txt
             ASN    ASN_atoms.csv       ASN_rixs.txt
             ASP    ASP_atoms.csv       ASP_rixs.txt

        where:
         - resname: resname found in the topology file
         - atnames: atom names in the residue(s) indicated in the corresponding index list
         - rixs: list of residue indices containing the residues constituted exactly by the corresponding list of atom names

        Defines
        
         - self.res_topology_folder: self.out / resid_topologies, 2 relevant file types resulting:
            - residuename(_topId)_rsix.txt -> with all residue indices for the respective topology
            - residuename(_topId)_atoms.csv -> csv file with atom names and masses
            NW: (_topId): numer of different topologies (different number/definition of atoms for each residue name)
                          present in the filename only if it is > 0

         - self.restypes: self.out / different_residue_topologies.csv
       '''

        start = timeit.default_timer()
        log.info('Detecting residue types')

        self.res_topology_folder = self.out / "resid_topologies"
        makedirs(self.res_topology_folder, exist_ok = True)

        # Find all different types of residues associated with the same resname:
        restypes = pd.DataFrame.from_dict({'resname':[], 'atnames':[], 'rixs':[], 'masses':[]})
        top_df = pd.read_csv(self.out / 'initial_topology_info.csv', header = 0, index_col = False,
                             na_values = self.nas, keep_default_na = False)
       
        restypes = get_residue_types(top_df, self.res_topology_folder)

        self.restypes = restypes
        self.restypes.to_csv(self.out / 'different_residue_topologies.csv', header = True, index = False)
        proc_time, tot_time = get_endtime(start, self.start)
        log.info('Residue-type detection complete; get_restypes took %s (h:m:s), total %s (h:m:s)',
                 proc_time, tot_time)

    def get_ff_parameters(self, 
                          maxToCheck: int = -1
                         ):
        '''
        Function that associate the unique residue topologies found to specific ff entries

         - verbose: whether to write report file of whole process
         - maxToCheck: max number of residues to check from ff, default (-1) = check all.
                       residue names in forcefield are sorted by similarity to topology residue name

        '''
        # To associate unique restypes to ff entries
        start = timeit.default_timer()

        complete, unmatched = associate_ff_params(self.ff, self.out, self.nas, maxToCheck)

        proc_time, tot_time = get_endtime(start, self.start)
        log.info('Residue-to-forcefield association complete; get_ff_parameters took %s (h:m:s), total %s (h:m:s)',
                 proc_time, tot_time)

        if not complete:
            # Build a readable list of the residues that could not be matched.
            lines = []
            for u in unmatched:
                matched, total = u.get('matched_atoms'), u.get('total_atoms')
                n_cand = u.get('n_candidates', 0)
                is_count = (isinstance(matched, (int, float)) and total is not None
                            and math.isfinite(matched) and matched >= 0)
                if is_count and n_cand:
                    lines.append(
                        f"  - {u['resname']} (topology {u['topology']}): no exact match; "
                        f"the best of {n_cand} candidate(s) shares {int(matched)}/{int(total)} "
                        f"atom names")
                else:
                    lines.append(
                        f"  - {u['resname']} (topology {u['topology']}): no candidate "
                        f"(none of the forcefield entries inspected, with or without "
                        f"terminal patches, has the same number of atoms; raise "
                        f"--max-to-check, currently {maxToCheck}, to inspect more)")
            detail = "\n".join(lines)
            # accept-topology works on the run root, the parent of Topology_recognition/
            run_root = self.out.parent
            example = next((u['topology'] for u in unmatched if u.get('n_candidates')),
                           unmatched[0]['topology'])
            msg = (
                f"Topology recognition could not associate {len(unmatched)} residue "
                f"topology/ies to the forcefield:\n{detail}\n"
                f"No forcefield entry exactly matches these atom names, so the "
                f"per-residue topology files needed downstream were not written and "
                f"the pipeline stops before the interactors stage. To continue, either:\n"
                f"  1. promote a candidate. List the candidates and how their atoms "
                f"pair with the structure's:\n"
                f"       mypsns accept-topology -o {run_root}\n"
                f"     then accept one per topology (by index, or a hand-edited CSV path):\n"
                f"       mypsns accept-topology -o {run_root} {example}=0\n"
                f"  2. add the missing residue(s) to the forcefield, or correct the "
                f"structure's atom names, and re-run `mypsns topology --overwrite`.\n"
                f"Details: '{self.out / UNMATCHED_TABLE}', candidates in "
                f"'{self.out / 'possible_topologies'}'.")
            log.error(msg)
            raise TopologyAssociationError(msg)


def describe_unresolved_topologies(topology_out: Path, na_values: list) -> list:
    '''
    What `mypsns accept-topology` lists when called without choices: one entry
    per residue topology that recognition could not match (including the ones
    already accepted by hand), with every proposed candidate and how its atoms
    pair with the structure's.

     - topology_out: <outpath>/Topology_recognition

    Each entry is a dict with key, resname, n_residues, accepted_from ('' if
    not accepted yet) and candidates: a list of dicts with index, file,
    ff_resnames (the forcefield residue names the candidate is built from,
    e.g. ['NH3+', 'MET'] for a MET with its N-terminal patch), pairs (the
    pair_candidate_atoms table), error (why it cannot be paired, or None) and
    same_as (index of an earlier candidate with identical content, or None).

    same_as matters because recognition often proposes the same topology
    several times: check_ters keeps an unpatched copy of the residue for every
    terminal patch that does not apply, so e.g. a GLY gets one GLY candidate
    per patch tried. Marking the copies leaves the real alternatives visible.
    Returns [] when there is no unmatching table.
    '''
    topology_out = Path(topology_out)
    table_path = topology_out / UNMATCHED_TABLE
    if not table_path.exists():
        return []
    table = read_topology_table(table_path, na_values)

    entries = []
    for _, row in table.iterrows():
        if row['ff_atoms'] and not row['accepted_from']:
            continue                     # matched exactly by recognition
        key = topology_key(row['atoms'])
        top_atoms = read_str_csv(topology_out / 'resid_topologies' / row['atoms'], na_values)
        with open(topology_out / 'resid_topologies' / row['rixs']) as f:
            n_residues = len(f.read().split())

        candidates = []
        first_with_content: dict = {}    # file content -> first candidate index
        for fname in candidate_files(row):
            path = topology_out / 'possible_topologies' / fname
            cand = read_str_csv(path, na_values)
            try:
                pairs, error = pair_candidate_atoms(top_atoms, cand), None
            except ValueError as e:
                pairs, error = None, str(e)
            content = path.read_bytes()
            index = candidate_index(fname)
            candidates.append({
                'index':       index,
                'file':        fname,
                'ff_resnames': list(dict.fromkeys(cand['resname'].astype(str))),
                'pairs':       pairs,
                'error':       error,
                'same_as':     first_with_content.get(content),
            })
            first_with_content.setdefault(content, index)
        entries.append({
            'key':           key,
            'resname':       row['resname'],
            'n_residues':    n_residues,
            'accepted_from': row['accepted_from'],
            'candidates':    candidates,
        })
    return entries


def accept_topologies(
        topology_out: Path,
        choices: dict,
        na_values: list,
        dry_run: bool = False,
        ) -> dict:
    '''
    Promote chosen candidates to the forcefield topology of residues that
    recognition could not match, so the interactors stage can run. This is the
    manual possible_topologies -> found_topologies fix, done safely.

     - topology_out: <outpath>/Topology_recognition
     - choices: {topology key: choice}. The key is the residue topology name
                (e.g. MET_1, as in the TopologyAssociationError message). The
                choice is either a candidate index n, meaning
                possible_topologies/<key>_uncertain_atoms_<n>.csv, or the path
                of a CSV with the same columns (e.g. a candidate edited by hand).
     - dry_run: check and report, but write nothing.

    For every choice the candidate's atoms are paired with the structure's
    (pair_candidate_atoms), and the structure's name of every atom is written
    to the structure_atom_name column next to atom_name, which (like bonds)
    keeps the forcefield's names. The interactors stage renames the atoms to
    the structure's names (structure_named), so the topology joins with the
    structure in the contacts stage. A pair that renames an atom of a
    different mass means the positional pairing is wrong: that choice is
    refused and the user is asked for a hand-edited CSV instead.

    Writes found_topologies/<key>_ff_atoms.csv, fills the row's ff_atoms (and
    accepted_from, the provenance) in unmatching_residue_topologies.csv, and
    once every row has ff_atoms writes complete_residue_topologies.csv, the file
    the interactors stage reads. All choices are checked before anything is
    written; any failure raises TopologyAcceptError listing every problem.
    Bond partners that are not atoms of the residue (dangling_bonds) only
    produce a warning.

    Returns {'accepted': [{key, source, target, renames}], 'pending': [keys
    still unresolved], 'complete': bool, 'dry_run': bool, 'warnings': [str]}.
    '''
    topology_out = Path(topology_out)
    table_path = topology_out / UNMATCHED_TABLE
    if not table_path.exists():
        if (topology_out / COMPLETE_TABLE).exists():
            raise TopologyAcceptError(
                f'{topology_out / COMPLETE_TABLE} exists and there is no {UNMATCHED_TABLE}: '
                f'recognition matched every residue topology, there is nothing to accept.')
        raise TopologyAcceptError(
            f'{table_path} not found: run `mypsns topology` first.')
    table = read_topology_table(table_path, na_values)

    row_of = {topology_key(a): ix for ix, a in table['atoms'].items()}
    pending = [k for k, ix in row_of.items() if not table.loc[ix, 'ff_atoms']]

    errors, warnings_, plans = [], [], []
    for key, choice in choices.items():
        choice = str(choice).strip()
        if key not in row_of:
            errors.append(f'{key}: no such residue topology. Unresolved: '
                          f'{", ".join(pending) or "none"}.')
            continue
        row = table.loc[row_of[key]]
        if row['ff_atoms'] and not row['accepted_from']:
            errors.append(f'{key}: already matched exactly by recognition '
                          f'({row["ff_atoms"]}), there is nothing to accept.')
            continue

        # resolve the choice to a source CSV
        files = candidate_files(row)
        if choice.isdigit():
            fname = f'{key}_uncertain_atoms_{int(choice)}.csv'
            if fname not in files:
                available = ', '.join(str(candidate_index(f)) for f in files)
                errors.append(f'{key}={choice}: no such candidate. Available: '
                              f'{available or "none, pass the path of a CSV instead"}.')
                continue
            source, label = topology_out / 'possible_topologies' / fname, fname
        else:
            source = Path(choice).expanduser()
            if not source.is_file():
                errors.append(f'{key}={choice}: file not found.')
                continue
            source = source.resolve()
            label = str(source)

        cand = read_str_csv(source, na_values)
        missing = [c for c in FF_ATOM_COLUMNS if c not in cand.columns]
        if missing:
            errors.append(f'{key}: {label} lacks the column(s) {missing}.')
            continue

        top_atoms = read_str_csv(topology_out / 'resid_topologies' / row['atoms'], na_values)
        try:
            pairs = pair_candidate_atoms(top_atoms, cand)
        except ValueError as e:
            errors.append(f'{key}: cannot pair {label} with the structure residue: {e}.')
            continue
        bad = pairs[pairs['renamed'] & ~pairs['mass_ok']]
        if not bad.empty:
            detail = ', '.join(f'{r.ff_name} ({r.ff_mass:g}) -> {r.top_name} ({r.top_mass:g})'
                               for r in bad.itertuples())
            errors.append(
                f'{key}: pairing {label} with the structure by position would rename '
                f'atoms of a different mass: {detail}. Copy the candidate, set its '
                f'atom_name column to the structure\'s names, and pass the copy '
                f'instead ({key}=path/to/file.csv).')
            continue

        renamed = pairs[pairs['renamed']]
        renames = dict(zip(renamed['ff_name'], renamed['top_name']))
        # pairs has one row per candidate atom, in candidate order
        ff_atoms = with_structure_names(cand, pairs['top_name'])
        dangling = dangling_bonds(ff_atoms)
        if dangling:
            detail = '; '.join(f'{a}: {" ".join(t)}' for a, t in dangling.items())
            warnings_.append(f'{key}: the bonds column names atoms the residue does not '
                             f'have ({detail}). Check {label}.')
        plans.append({
            'key':      key,
            'ix':       row_of[key],
            'source':   label,
            'target':   f'{key}_ff_atoms.csv',
            'renames':  renames,
            'ff_atoms': ff_atoms,
        })

    if errors:
        msg = ('accept-topology could not apply these choices, so nothing was written:\n'
               + '\n'.join(f'  - {e}' for e in errors))
        log.error(msg)
        raise TopologyAcceptError(msg)

    found_dir = topology_out / 'found_topologies'
    for p in plans:
        table.loc[p['ix'], 'ff_atoms'] = p['target']
        table.loc[p['ix'], 'accepted_from'] = p['source']
        log.info('accept-topology: %s <- %s%s', p['key'], p['source'],
                 ' (structure names: ' + ', '.join(f'{a}->{b}' for a, b in p['renames'].items()) + ')'
                 if p['renames'] else '')
        if not dry_run:
            p['ff_atoms'].to_csv(found_dir / p['target'], header = True, index = False)
    for w in warnings_:
        log.warning('accept-topology: %s', w)

    pending = [k for k, ix in row_of.items() if not table.loc[ix, 'ff_atoms']]
    complete = not pending
    if not dry_run:
        table.to_csv(table_path, header = True, index = False)
        if complete:
            # same layout save_results gives a complete match (no possible_atoms);
            # accepted_from stays as the record of what was accepted by hand
            table.drop(columns = ['possible_atoms']).to_csv(
                topology_out / COMPLETE_TABLE, header = True, index = False)
            log.info('accept-topology: every residue topology resolved, wrote %s',
                     topology_out / COMPLETE_TABLE)

    return {
        'accepted': [{k: p[k] for k in ('key', 'source', 'target', 'renames')} for p in plans],
        'pending':  pending,
        'complete': complete,
        'dry_run':  dry_run,
        'warnings': warnings_,
    }