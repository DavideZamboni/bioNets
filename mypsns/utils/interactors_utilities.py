"""
interactors_utilities.py – helper functions for the Interactors service.
 
Public API (called from interactors.py)
----------------------------------------
reconstruct_whole_topology  – orchestrator; calls five private helpers below
get_charge_based_interactors         – assign PIE interaction roles to topology atoms
get_named_interactors  – assign PyInteraph interaction roles to topology atoms
 
Private helpers (called only by reconstruct_whole_topology)
------------------------------------------------------------
 extract_chain_bounds     – parse chain start/end residue indices from Universe
 fill_bond_info           – complete bidirectional bond lists for each residue type
 build_wholetop           – concatenate per-residue DataFrames into one table
 detect_chains            – detect chain boundaries from ter entries / index jumps
 assign_chains_and_entity – stamp Chain and Entity columns; fix polymer bonds;
                             adjust ter chargegroups; assign atom indices
 
Supporting utilities
--------------------
calc_counts          – contact-index of residue atoms toward ter entry atoms
adjust_ters          – re-assign chargegroups of ter atoms
checkif_polymer      – classify a chain as Polymer / Monomer / Unknown
add_poly_bonds       – add -X back-bonds for polymer chains
get_nextbonds        – collect all +X forward-bond entries across the topology
below_max_num_bonds  – BFS reachability (used by add_hc_atoms)
cluster              – greedy bond-contiguous grouping of HC atoms
get_1st_letter       – index of first non-digit character in an atom name
add_hc_atoms         – assign HC clusters and q_interactors = 'HC'
add_sb_atoms         – assign Acidic / Basic q_interactors
add_hbonds           – assign HAcc / HDon / HBond q_interactors
"""
from __future__ import annotations
 
import numpy as np
import pandas as pd
from os import makedirs
from pathlib import Path
from typing import Dict, List, Literal, Optional, Set, Tuple, Union
import copy
import networkx as nx
import MDAnalysis as mda
import logging

from mypsns.config import NAMED_INTERACTOR_DEFAULTS

log = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════════
# Named-interactor (PyInteraph-style) role definitions
# ═══════════════════════════════════════════════════════════════════════════════
# NAMED_INTERACTOR_DEFAULTS lives in mypsns/config.py (single source of truth;
# it used to be duplicated here and kept in sync by hand). Each role maps to a
# list of RULES. A rule is a dict {attribute: values}; all conditions in one
# dict are AND-combined, and the rules in a role are OR-combined (an atom
# matches the role if it satisfies ANY rule). A condition key prefixed with
# "not " negates it. This is the YAML-native form recorded in run_config.yaml
# and overridable from config; it is converted to the internal (attr, values)
# grammar by `_rules_from_config`. get_named_interactors deep-copies it before
# use, so the shared constant is never mutated. Imported above and
# re-exported from this module for existing callers.

# Roles evaluated as selection rules (everything except the scalar hb_masses).
NAMED_INTERACTOR_ROLES: List[str] = ["HC", "HDon", "HAcc", "Basic", "Acidic"]


def _rules_from_config(role_rules: list) -> List[list]:
    """Convert a role's config rules (list of ``{attribute: values}`` dicts) into
    the internal grammar `apply_rule` expects: a list of rules, each a list of
    ``(attribute, values)`` tuples (AND-combined). Tolerates rules already given
    in tuple/list grammar (e.g. from code) so both forms work."""
    out: List[list] = []
    for rule in role_rules:
        if isinstance(rule, dict):
            conds = [(str(k), (v if isinstance(v, list) else [v]))
                     for k, v in rule.items()]
            if conds:
                out.append(conds)
        elif isinstance(rule, (list, tuple)):
            # already grammar: a single (attr, vals) or a list of them
            if rule and isinstance(rule[0], str):
                out.append([(rule[0], rule[1])])           # (attr, vals)
            else:
                conds = []
                for cond in rule:
                    if isinstance(cond, (list, tuple)) and len(cond) == 2:
                        conds.append((cond[0], cond[1]))
                if conds:
                    out.append(conds)
    return out


# ═══════════════════════════════════════════════════════════════════════════════
# Low-level chargegroup helpers
# ═══════════════════════════════════════════════════════════════════════════════
 
def calc_counts(
    res_df: pd.DataFrame,
    bonds: Dict[str, List[str]],
) -> Tuple[Dict[str, float], float, Optional[str]]:
    """
    Compute the *contact index* of each non-ter atom toward the ter entry.
 
    The contact index of atom *a* is::
 
        contact_index(a) = |{bonds(a) that are in the ter entry}| / |bonds(a)|
 
    Parameters
    ----------
    res_df:
        DataFrame for a single ff residue containing a ter entry.
        Rows with ``chargegroup == -1`` are the ter atoms.
    bonds:
        ``{atom_name: [bound_atom_name, ...]}`` for all atoms in the residue.
 
    Returns
    -------
    counts:   ``{atom_name: contact_index}`` for non-ter atoms.
    maxc:     Highest contact index strictly less than 1.
    maxc_an:  Atom name achieving *maxc*; ``None`` if none found.
    """
    ter_atoms   = set(res_df.loc[res_df['chargegroup'] == -1, 'atom_name'])
    other_atoms = res_df.loc[res_df['chargegroup'] != -1, 'atom_name'].tolist()
 
    counts: Dict[str, float] = {}
    maxc:    float             = 0.0
    maxc_an: Optional[str]    = None
 
    for an in other_atoms:
        bound = bonds.get(an, [])
        if not bound:
            continue
        in_ter = sum(1 for b in bound if b in ter_atoms)
        c = in_ter / len(bound)
        counts[an] = c
        if 0.0 < c < 1.0 and c > maxc:
            maxc    = c
            maxc_an = an
 
    return counts, maxc, maxc_an
 
def get_all_bonds(data: pd.DataFrame) -> Dict[str, Set[str]]:
    '''
    to obtain allbonds: Bidirectional bond dict for this residue.

    '''
    allbonds: Dict[str, Set[str]] = {}
    for _, row in data.iterrows():
        an = row['atom_name']
        bs = [b for b in row['bonds'].split() if b != '']

        if an not in allbonds.keys():
            allbonds[an] = set()
        
        for b in bs:
            allbonds[an].add(b)
            if b not in allbonds.keys():
                allbonds[b] = set([an])
            else:
                allbonds[b].add(an)
    return allbonds


def adjust_ters(
    data: pd.DataFrame,
) -> pd.DataFrame:
    """
    Re-assign the ``chargegroup`` of ter atoms (default value ``-1``).
 
    **Strategy**
 
    Pass 1 – single-bond bridges:
        For each non-ter atom *a* with ``contact_index == 1`` and only one
        bond: if *a*'s chargegroup is large enough, the ter atom inherits it;
        otherwise *a* itself is marked ``-1``.
 
    Pass 2 – remaining unassigned atoms:
        Handled by the atom with the highest contact index.  If the bridge
        atom has only one bond into the ter entry:
        * ter atoms with ``|ter| > 2`` form their own new chargegroup;
        * smaller ter groups inherit the bridging atom's chargegroup.
 
    Parameters
    ----------
    data:     DataFrame of the ff residue including ter atoms.
    s:        Running log string.
 
    Returns
    -------
    trial: Updated DataFrame with corrected chargegroups.
    s:     Updated log string.
    """
    trial = data.copy().reset_index(drop=True)
    log.debug('adjust_ters: adjusting chargegroups for ter entry (%s atoms)', trial.shape[0])
    allbonds = get_all_bonds(data)

    # ── Pass 1: single-bond bridges ───────────────────────────────────────────
    counts, maxc, maxc_an = calc_counts(trial, allbonds)
    log.debug('adjust_ters: pass-1 max contact index %s (atom %s)', maxc, maxc_an)

    for an, c in counts.items():
        log.debug('adjust_ters:   %s contact index %s', an, c)
        if c == 1.0 and len(allbonds.get(an, [])) == 1:
            res_at  = trial[trial['atom_name'] == an]
            ter_at  = trial[trial['atom_name'].isin(allbonds[an])]
            at_chgg = trial[
                trial['chargegroup'].isin(res_at['chargegroup'].tolist()) &
                (round(trial['mass']) > 1)
            ]
            if at_chgg.shape[0] > 1:
                trial.loc[ter_at.index, 'chargegroup'] = (
                    trial.loc[trial['atom_name'] == an, 'chargegroup'].values
                )
            else:
                trial.loc[res_at.index, 'chargegroup'] = -1
 
    # ── Pass 2: remaining unassigned atoms ────────────────────────────────────
    ter_df = trial[trial['chargegroup'] == -1]
    counts, maxc, maxc_an = calc_counts(trial, allbonds)
    log.debug('adjust_ters: pass-2 max contact index %s (atom %s)', maxc, maxc_an)


    if maxc_an is not None:
        if maxc * len(allbonds.get(maxc_an, [])) <= 1:
            if ter_df.shape[0] > 2:
                next_cg = int(trial['chargegroup'].max()) + 1
                trial.loc[trial['chargegroup'] == -1, 'chargegroup'] = next_cg
            else:
                my_cg = trial.loc[trial['atom_name'] == maxc_an, 'chargegroup'].values
                trial.loc[trial['chargegroup'] == -1, 'chargegroup'] = my_cg
        else:
            log.warning('adjust_ters: bridging atom has multiple bonds to ter - unhandled')
    elif ter_df.shape[0] != 0:
        if ter_df.shape[0] > 2:
            next_cg = int(trial['chargegroup'].max()) + 1
            trial.loc[trial['chargegroup'] == -1, 'chargegroup'] = next_cg
        else:
            log.warning('adjust_ters: ter entry too small for own chargegroup - left at -1')

    return trial
 
 
# ═══════════════════════════════════════════════════════════════════════════════
# Chain / polymer helpers
# ═══════════════════════════════════════════════════════════════════════════════
 
def add_poly_bonds(chain_dataframe: pd.DataFrame) -> List[Tuple[int, str]]:
    """
    Add ``-prevAtom`` back-bonds to atoms in residue *i+1* that are targeted
    by ``+X`` forward bonds in residue *i*.
 
    Parameters
    ----------
    chain_dataframe:
        Rows of the whole-topology table for a single polymer chain.
 
    Returns
    -------
    List of ``(row_index, new_bond_string)`` pairs.
    """
    rixs = sorted(np.unique(chain_dataframe['Resindex'].values).tolist())
    newbonds: List[Tuple[int, str]] = []
 
    for count in range(1, len(rixs)):
        prev = chain_dataframe[chain_dataframe['Resindex'] == rixs[count - 1]]
        prev_fwd = [
            row for _, row in prev.iterrows()
            if isinstance(row['bonds'], str) and '+' in row['bonds']
        ]
        cur = chain_dataframe[chain_dataframe['Resindex'] == rixs[count]]
 
        for ix, row in cur.iterrows():
            an     = row['atom_name']
            bs     = row['bonds'] if isinstance(row['bonds'], str) else ''
            truebs = (
                [b for b in bs.split() if b and not b.startswith('+')]
                if count == len(rixs) - 1
                else [b for b in bs.split() if b]
            )
            for prevat in prev_fwd:
                targets = [b[1:] for b in prevat['bonds'].split() if b.startswith('+')]
                if an in targets:
                    truebs.append(f'-{prevat["atom_name"]}')
            newbonds.append((ix, ' '.join(truebs)))
 
    return newbonds

def get_nextbonds(
    wholetop: pd.DataFrame,
) -> Dict[str, Dict[str, Dict[str, Set[str]]]]:
    """
    Collect all ``+X`` forward-bond entries from the topology, grouped by
    class → residue name → atom name.
 
    Parameters
    ----------
    wholetop: Full concatenated topology DataFrame.
    s:        Running log string.
 
    Returns
    -------
    bonds_to_next: ``{rtp_file: {resname: {atom_name: {next_atom_names}}}}``.
    s:             Updated log string.
    """
    bonds_to_next: Dict[str, Dict[str, Dict[str, Set[str]]]] = {}
    #print(wholetop.head())
    for c, cdf in wholetop.groupby('rtp_file'):
        for resn, resdf in cdf.groupby('resname'):
            for atn, atdf in resdf.groupby('atom_name'):
                for _, row in atdf.iterrows():
                    if not isinstance(row['bonds'], str):
                        continue
                    fwd = [b[1:] for b in row['bonds'].split() if b.startswith('+')]
                    if not fwd:
                        continue
                    (bonds_to_next
                     .setdefault(c, {})
                     .setdefault(resn, {})
                     .setdefault(atn, set())
                     .update(fwd))
 
    for rtp_file in bonds_to_next:
        for resn in bonds_to_next[rtp_file]:
            for atn, nexts in bonds_to_next[rtp_file][resn].items():
                log.debug('forward bond +X: rtp_file=%s res=%s atom=%s -> %s', rtp_file, resn, atn, nexts)
    return bonds_to_next

# ═══════════════════════════════════════════════════════════════════════════════
# helpers for reconstruct_whole_topology
# ═══════════════════════════════════════════════════════════════════════════════
 
def extract_chain_bounds(
    universe: mda.Universe,
) -> Dict[Optional[str], List[int]]:
    """
    Parse chain start/end residue indices from a MDAnalysis Universe.
 
    Returns
    -------
    top_chains: ``{chainID: [first_resid, last_resid]}``.
    s:          Updated log string.
    """
    # get defined chains
    try:
        chids = list({(at.resnum, at.chainID) for at in universe.atoms})
    except mda.exceptions.NoDataError:
        chids = list({(at.resnum, None) for at in universe.atoms})
 
    chids.sort(key=lambda x: x[0])
    # dict to store start/end of each chain
    top_chains: Dict[str, Tuple[int, int]] = {}
 
    for ri, c in chids:
        if c not in top_chains:
            top_chains[c] = [ri, ri]
        else:
            if ri > top_chains[c][1]:
                top_chains[c][1] = ri
 
    for c, (first, last) in top_chains.items():
        log.debug('topology chain %s: residues %s -> %s', c, first, last)

    return top_chains
 
def fill_bond_info(
    topologies: List[Path],
    nas: List[str]
) -> Dict[str, Path]:
    """
    Build a bidirectional bond dict for each per-residue CSV, write the
    completed version to *tmp_dir*, and record bond info for ter residues.
 
    Parameters
    ----------
    topologies: Paths to per-residue ``<TopName>.csv`` files.
    s:          Running log string.
 
    Returns
    -------
    top_df_dict: completed_csv_path
    s:          Updated log string.
    """
    top_df_dict: Dict[str, Path] = {}
 
    for res_df_path in topologies:
        rdf_name = Path(res_df_path).name
        res_df = pd.read_csv(res_df_path, header=0, index_col=False, na_values = nas)
 
        if res_df.shape[0] <= 1:
            top_df_dict[rdf_name] = res_df_path
            log.debug('bond info completed: %s', rdf_name)
            continue
 
        # Build symmetric bond dict
        allbonds: Dict[str, List[str]] = {}
        for _, row in res_df.iterrows():
            an    = row['atom_name']
            bonds = row['bonds']
            allbonds.setdefault(an, [])
            if isinstance(bonds, str): # na float
                for b in bonds.split():
                    if b and b not in allbonds[an]:
                        allbonds[an].append(b)
                    allbonds.setdefault(b, [])
                    if an not in allbonds[b]:
                        allbonds[b].append(an)
 
        for ix, row in res_df.iterrows():
            res_df.loc[ix, 'bonds'] = ' '.join(allbonds.get(row['atom_name'], []))

        res_df.to_csv(res_df_path, header = True, index = False)
        top_df_dict[rdf_name] = res_df_path
        log.debug('bond info completed: %s', rdf_name)

    return top_df_dict

def build_wholetop(
    top_df_dict: Dict[str, Path],
    resid_to_top: pd.DataFrame,
    nas: List[str]
) -> Tuple[pd.DataFrame, List[int]]:
    """
    Concatenate per-residue DataFrames into a single whole-topology table,
    stamping each row with its ``Resindex``.
 
    Also re-keys *terbonds* from type-name strings to numeric Resindex values.
 
    Parameters
    ----------
    top_df_ass:   ``{stem: completed_csv_path}``.
    terbonds:     ``{stem: bond_dict}`` (updated in-place with Resindex keys).
    resid_to_top: Table with columns ``['Resindex', 'topology']``.
    s:            Running log string.
 
    Returns
    -------
    wholetop: Concatenated DataFrame (index reset).
    terbonds: Updated with ``{Resindex: bond_dict}`` entries.
    allres:   Sorted list of unique residue indices.
    s:        Updated log string.
    """
    chunks: List[pd.DataFrame] = []
    resid_to_top = resid_to_top.sort_values('Resindex')
 
    for _, row in resid_to_top.iterrows():
        rdf_name = row['topology']
        tmp = pd.read_csv(top_df_dict[rdf_name], header=0, index_col=False, na_values=nas)
        tmp['Resindex'] = row['Resindex']
        chunks.append(tmp)
 
    wholetop = pd.concat(chunks, ignore_index=True)
    allres   = sorted(wholetop['Resindex'].unique().tolist())
    log.debug('whole topology assembled: %s atoms, %s residues', wholetop.shape[0], len(allres))
    return wholetop, allres
 
def detect_chains(
    wholetop: pd.DataFrame,
    allres: List[int],
) -> Tuple[Dict[str, Tuple[int, int, bool]], pd.DataFrame]:
    """
    Two-pass chain detection:
      1) Walk +bond connectivity to find polymer stretches.
      2) Assign remaining residues as non-polymer chains, grouped by
         contiguous runs sharing the same rtp_file.
    """
    LETTERS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    # ── Pass 1: find polymer chains via +bond connectivity ────────────
    polymer_resids: Set[int] = set()
    poly_chains: List[List[int]] = []  # each entry is a list of contiguous resids
    current_chain: List[int] = []

    for i, resid in enumerate(allres):
        res_df = wholetop[wholetop['Resindex'] == resid]

        # collect atoms this residue expects in the next residue
        to_next = []
        for _, row in res_df.iterrows():
            bonds_str = row.get('bonds', '')
            if not isinstance(bonds_str, str):
                continue
            for b in bonds_str.split():
                if b.startswith('+'):
                    to_next.append(b.lstrip('+'))

        linked_forward = False
        if to_next and i + 1 < len(allres):
            next_rix = allres[i + 1]
            next_ats = wholetop[
                (wholetop['Resindex'] == next_rix)
                & (wholetop['atom_name'].isin(to_next))
            ]
            if next_ats.shape[0] == len(to_next):
                linked_forward = True

        if linked_forward:
            if not current_chain:
                current_chain.append(resid)
            elif current_chain[-1] != resid:
                current_chain.append(resid)
            # also mark the next residue (it's the target of the + bond)
        else:
            if current_chain:
                # this residue is the last in the polymer stretch
                current_chain.append(resid)
                poly_chains.append(current_chain)
                polymer_resids.update(current_chain)
                current_chain = []

    # close any open chain
    if current_chain:
        poly_chains.append(current_chain)
        polymer_resids.update(current_chain)

    # ── Pass 2: group leftover residues by contiguous rtp_file runs ───
    nonpoly_chains: List[List[int]] = []
    current_group: List[int] = []
    current_rtp: Optional[str] = None

    for resid in allres:
        if resid in polymer_resids:
            # flush any pending non-polymer group
            if current_group:
                nonpoly_chains.append(current_group)
                current_group = []
                current_rtp = None
            continue

        rtp = wholetop[wholetop['Resindex'] == resid]['rtp_file'].iloc[0]
        if current_rtp is None:
            current_rtp = rtp
            current_group = [resid]
        elif rtp == current_rtp:
            current_group.append(resid)
        else:
            nonpoly_chains.append(current_group)
            current_group = [resid]
            current_rtp = rtp

    if current_group:
        nonpoly_chains.append(current_group)

    # ── Assign chain letters ──────────────────────────────────────────
    # polymers first (in residue order), then non-polymers
    all_segments = []
    for pc in poly_chains:
        all_segments.append((pc[0], pc[-1], True))
    for nc in nonpoly_chains:
        all_segments.append((nc[0], nc[-1], False))

    # sort by first resindex so letters follow residue order
    all_segments.sort(key=lambda x: x[0])

    chains: Dict[str, Tuple[int, int, bool]] = {}
    for idx, (first, last, is_poly) in enumerate(all_segments):
        if idx < len(LETTERS):
            cid = LETTERS[idx]
        else:
            cid = f'Z{idx - len(LETTERS) + 1}'
        chains[cid] = (first, last, is_poly)

    # ── Log ───────────────────────────────────────────────────────────
    for c in sorted(chains.keys()):
        first, last, is_poly = chains[c]
        log.debug('chain %s: %s -> %s, polymer=%s', c, first, last, is_poly)

    return chains, wholetop

def get_resname(rdf: pd.DataFrame) -> str:
    maxc = 0
    resname = ''
    for r, rn_df in rdf.groupby(by = 'resname'):
        if rn_df.shape[0] > maxc:
            resname = r
            maxc = rn_df.shape[0]
    return resname

## NEW METHOD
def build_residue_variants(
    wholetop:   pd.DataFrame,
    mychains:   Dict[str, Tuple[int, int, bool]],
    universe:   mda.Universe,
    outpath:    Path,
    topology_resname_map: Optional[Dict[int, str]] = None,
) -> None:
    """
    Per-residue reconstruction → deduplication → write three normalised tables.

    Steps
    -----
    1. Stamp Chain and Entity columns on wholetop.
    2. Add polymer back-bonds (-X) via add_poly_bonds.
    3. Adjust ter chargegroups via adjust_ters on terminal polymer residues.
    4. Iterate residues; group by hash of intrinsic columns; first time a
       hash is seen, save its rows to found_topologies/<variant_id>.csv.
    5. Write residue_index.csv (one row per residue, including the input-topology
       resname for selection building) and atom_index_map.csv (one row per atom).

    Parameters
    ----------
    topology_resname_map:
        ``{Resindex: input-topology resname}`` (e.g. HIS), written to
        residue_index.csv as ``topology_resname``. This is the name the
        structure/PDB uses and the correct basis for MDAnalysis selections,
        unlike variant_id which carries the rtp/forcefield name (e.g. HSE).
    """
    topology_resname_map = topology_resname_map or {}
    found_dir = outpath / 'found_topologies'
    makedirs(found_dir, exist_ok=True)

    wholetop = wholetop.reset_index(drop=True).copy()
    wholetop['bonds'] = wholetop['bonds'].fillna('')

    # ── Chain + Entity ──────────────────────────────────────────────────
    for c, (first, last, poly) in mychains.items():
        mask = (wholetop['Resindex'] >= first) & (wholetop['Resindex'] <= last)
        wholetop.loc[mask, 'Chain']  = c
        wholetop.loc[mask, 'Entity'] = 'P' if poly else 'M'

    # ── Polymer back-bonds ──────────────────────────────────────────────
    for c, (first, last, poly) in mychains.items():
        if not poly:
            continue
        cdf = wholetop[wholetop['Chain'] == c]
        for ix, bond_str in add_poly_bonds(cdf):
            wholetop.loc[ix, 'bonds'] = bond_str

    # ── Ter chargegroup adjustments ─────────────────────────────────────
    for c, (first, last, poly) in mychains.items():
        if not poly:
            continue
        for terminal_resix in (first, last):
            tres = wholetop[wholetop['Resindex'] == terminal_resix]
            if tres[tres['chargegroup'] == -1].shape[0] == 0:
                continue
            tres_adj = adjust_ters(tres)
            tres_adj.index = tres.index
            wholetop.loc[tres.index] = tres_adj

    # ── Variant deduplication ───────────────────────────────────────────
    intrinsic_cols = [
        'rtp_file', 'resname', 'atom_name', 'atom_type',
        'charge', 'chargegroup', 'bonds', 'mass',
    ]
    intrinsic_cols = [c for c in intrinsic_cols if c in wholetop.columns]

    variants_seen:  Dict[str, str] = {}    # hash → variant_id
    variant_counts: Dict[str, int] = {}    # resname → counter
    res_to_variant: Dict[int, str] = {}

    allres = sorted(wholetop['Resindex'].unique().tolist())
    for resix in allres:
        rdf = (wholetop[wholetop['Resindex'] == resix][intrinsic_cols]
               .sort_values('atom_name')
               .reset_index(drop=True))
        h = str(int(pd.util.hash_pandas_object(rdf, index=False).sum()))

        if h in variants_seen:
            res_to_variant[resix] = variants_seen[h]
            continue

        resname = get_resname(rdf)
        n = variant_counts.get(resname, 0) + 1
        variant_counts[resname] = n
        variant_id = resname if n == 1 else f'{resname}_{n}'

        variants_seen[h]        = variant_id
        res_to_variant[resix]   = variant_id
        rdf.to_csv(found_dir / f'{variant_id}.csv',
                   header=True, index=False)
        log.debug('variant %s (first instance: Resindex=%s)', variant_id, resix)

    log.info('Unique residue variants: %s', len(variants_seen))

    # ── residue_index.csv ───────────────────────────────────────────────
    chain_map  = {r: wholetop.loc[wholetop['Resindex'] == r, 'Chain'].iloc[0]
                  for r in allres}
    entity_map = {r: wholetop.loc[wholetop['Resindex'] == r, 'Entity'].iloc[0]
                  for r in allres}
    resnum_map = {int(res.resindex): int(res.resnum)
                  for res in universe.residues}

    res_idx = pd.DataFrame([{
        'Resindex':         int(r),
        'resnum':           resnum_map.get(int(r), int(r)),
        'Chain':            chain_map[r],
        'Entity':           entity_map[r],
        'variant_id':       res_to_variant[r],
        'topology_resname': topology_resname_map.get(int(r), ''),
    } for r in allres])
    res_idx.to_csv(outpath / 'residue_index.csv', header=True, index=False)
    log.debug('residue_index.csv written (%s residues)', res_idx.shape[0])

    # ── atom_index_map.csv ──────────────────────────────────────────────
    atom_rows = [
        {'Resindex':   int(res.resindex),
         'atom_name':  at.name,
         'atom_index': int(at.id)}
        for res in universe.residues for at in res.atoms
    ]
    aim = (pd.DataFrame(atom_rows)
           .sort_values(['Resindex', 'atom_index'])
           .reset_index(drop=True))
    aim.to_csv(outpath / 'atom_index_map.csv', header=True, index=False)
    log.debug('atom_index_map.csv written (%s atoms)', aim.shape[0])

def get_resid_to_top(resid_df: pd.DataFrame, outpath: Path) -> pd.DataFrame:
    """Map every Resindex to its forcefield topology file and to the residue
    name as found in the INPUT topology (``resname`` column of
    complete_residue_topologies.csv).

    The input-topology resname (e.g. HIS) is what the structure/PDB actually
    uses, unlike the rtp/forcefield name baked into variant_id (e.g. HSE), so it
    is the correct name to build MDAnalysis selections from. Returned columns:
    ``['Resindex', 'topology', 'topology_resname']``.
    """
    has_resname = 'resname' in resid_df.columns
    data = {'Resindex': [], 'topology': [], 'topology_resname': []}
    for ix, row in resid_df.iterrows():
        ix_file = outpath / 'resid_topologies' / row['rixs']
        with open(ix_file, 'r') as f:
            resixs = [int(i) for i in f.readline().split(' ') if i != '']
            for myix in resixs:
                data['Resindex'].append(myix)
                data['topology'].append(row['ff_atoms'])
                data['topology_resname'].append(str(row['resname']) if has_resname else '')
    return pd.DataFrame.from_dict(data)

# ═══════════════════════════════════════════════════════════════════════════════
# Charge-based interaction-type assignment helpers
# ═══════════════════════════════════════════════════════════════════════════════

def get_1st_letter(bondlist: List[str]) -> List[Tuple[str, int]]:
    """
    Return ``(entry, index_of_first_alphabetic_char)`` for each entry.
    Handles atom names beginning with digits (common in carbohydrates).
    """
    result: List[Tuple[str, int]] = []
    for b in bondlist:
        idx = next((i for i, ch in enumerate(b) if not ch.isdigit()), len(b))
        result.append((b, idx))
    return result

def below_max_num_bonds(
    ref_atom:       Tuple[str, int],
    atoms_to_reach: List[Tuple[str, int]],
    df:             pd.DataFrame,
    max_bonds:      int  = 2,
    verbose:        bool = False,
) -> int:
    """
    BFS reachability: minimum bond distance from *ref_atom* to any atom in
    *atoms_to_reach*, capped at ``max_bonds + 1``.

    Parameters
    ----------
    ref_atom:       ``(atom_name, Resindex)`` of the start atom.
    atoms_to_reach: ``[(atom_name, Resindex)]`` targets.
    df:             Residue (+ neighbour) DataFrame with ``bonds`` column.
    max_bonds:      Search depth ceiling.
    verbose:        Print intermediate BFS state (debug).

    Returns
    -------
    int – minimum distance found, or ``max_bonds + 1`` if not reachable.
    """
    rixs      = sorted(df['Resindex'].unique().tolist())
    reach_set = set(atoms_to_reach)
    visited:  Set[Tuple[str, int]] = {ref_atom}
    frontier: List[Tuple[str, int]] = [ref_atom]
    bonds = 0

    while bonds <= max_bonds and frontier:
        next_frontier: List[Tuple[str, int]] = []
        for atn, rix in frontier:
            rows = df[(df['atom_name'] == atn) & (df['Resindex'] == rix)]
            for _, at_row in rows.iterrows():
                if not isinstance(at_row['bonds'], str):
                    continue
                for b in at_row['bonds'].split():
                    if b.startswith('H'):
                        continue

                    if b.startswith('-'):
                        if rix in rixs[1:]: # cant be first
                            brix      = rixs[rixs.index(rix) - 1]
                            neighbour = (b[1:], brix)
                        else:
                            continue   # no previous residue
                    elif b.startswith('+'):
                        if rix in rixs[:-1]: # cant be last
                            brix      = rixs[rixs.index(rix) + 1]
                            neighbour = (b[1:], brix)
                        else:
                            continue   # no next residue
                    else:
                        neighbour = (b, rix)

                    if neighbour in reach_set:
                        return bonds + 1
                    if neighbour not in visited:
                        visited.add(neighbour)
                        next_frontier.append(neighbour)
        frontier = next_frontier
        bonds   += 1
    return bonds

def get_terminal_tree(resdf: pd.DataFrame) -> Dict[int, List[int]]:
    """
    Build a map: root_index → [all recursively terminal atom indices].
    
    An atom is "terminal" if it has exactly one intra-residue heavy-atom
    bond (i.e. it's a leaf in the heavy-atom bond graph).
    Terminal-of-terminal atoms propagate up to the same root.
    """
    res_names = set(resdf['atom_name'])
    
    # For each atom, find its intra-residue heavy-atom neighbors
    heavy_neighbors: Dict[int, List[int]] = {}
    for ix, row in resdf.iterrows():
        bonds_str = row.get('bonds', '')
        if not isinstance(bonds_str, str) or bonds_str.strip() == '':
            heavy_neighbors[ix] = []
            continue
        # only intra-residue bonds (no +/- prefixes)
        bond_names = [
            b for b in bonds_str.split()
            if not b.startswith('+') and not b.startswith('-')
            and b in res_names
        ]
        # heavy neighbors only
        neighbors = resdf[
            resdf['atom_name'].isin(bond_names)
            & (round(resdf['mass']) != 1)
        ].index.tolist()
        heavy_neighbors[ix] = neighbors

    is_heavy = {ix: round(resdf.at[ix, 'mass']) != 1 for ix in resdf.index}

    # Find terminal atoms: those with exactly 1 heavy-atom neighbor
    # H atoms are always terminal (0 or 1 heavy neighbor)
    terminal_parent: Dict[int, int] = {}  # terminal_ix → parent_ix
    for ix in resdf.index:
        hn = heavy_neighbors[ix]
        if is_heavy[ix]:
            if len(hn) == 1:
                terminal_parent[ix] = hn[0]
        else:
            # H atom: find its heavy parent from bonds
            bonds_str = resdf.at[ix, 'bonds']
            if not isinstance(bonds_str, str):
                continue
            bond_names = [
                b for b in bonds_str.split()
                if not b.startswith('+') and not b.startswith('-')
                and b in res_names
            ]
            h_heavy = resdf[
                resdf['atom_name'].isin(bond_names)
                & (round(resdf['mass']) != 1)
            ].index.tolist()
            if len(h_heavy) == 1:
                terminal_parent[ix] = h_heavy[0]

    # Propagate: if parent is itself terminal, walk up to root
    def find_root(ix: int) -> int:
        visited = {ix}
        current = ix
        while current in terminal_parent:
            parent = terminal_parent[current]
            if parent in visited:
                break  # cycle guard
            if parent in terminal_parent:
                visited.add(parent)
                current = parent
            else:
                return parent
        return terminal_parent.get(current, current)

    # Build root → [all terminals] map
    root_terminals: Dict[int, List[int]] = {}
    for tix in terminal_parent:
        root = find_root(tix)
        root_terminals.setdefault(root, []).append(tix)

    return root_terminals

def sb_charged_atoms(
    resdf:          pd.DataFrame,
    min_abs_charge: float,
) -> List[Tuple[str, int]]:
    '''
    Heavy root atoms of any SB-forming group in `resdf`, using the same recursive
    local-charge rule as add_sb_atoms: local_q = own_q + sum(charges of recursive
    terminals); group is SB-active when |local_q| >= min_abs_charge.

    Returns (atom_name, Resindex) tuples — the form below_max_num_bonds consumes.
    '''
    root_terminals = get_terminal_tree(resdf)
    all_terminals: Set[int] = set()
    for terms in root_terminals.values():
        all_terminals.update(terms)

    heavy = resdf[round(resdf['mass']) != 1]
    out: List[Tuple[str, int]] = []
    for ix, row in heavy.iterrows():
        if ix in all_terminals:
            continue
        term_ixs = root_terminals.get(ix, [])
        term_q   = resdf.loc[term_ixs, 'charge'].sum() if term_ixs else 0.0
        if abs(row['charge'] + term_q) >= min_abs_charge:
            out.append((row['atom_name'], int(row['Resindex'])))
    return out

def score(
        sizes: List[int], 
        min_num_heavy_ats: int, 
        max_num_heavy_ats: int
    ) -> Tuple[int, int]:
    viol = sum(
        max(0, s - max_num_heavy_ats) + max(0, min_num_heavy_ats - s)
        for s in sizes
    )
    return (viol, -min(sizes, default=0))

def split(
        G: nx.Graph,
        nodes: List[str],
        min_num_heavy_ats: int,
        max_num_heavy_ats: int
        ) -> List[List[str]]:
    
    if len(nodes) <= max_num_heavy_ats:
        return [nodes]
    
    sub = G.subgraph(nodes)
    edges = list(sub.edges())
    if not edges:
        return [nodes]
    best = (float('inf'), float('inf')), None, None
    for u, v in edges:
        G.remove_edge(u, v)
        comps = [list(c) for c in nx.connected_components(G.subgraph(nodes))]
        sc = score([len(c) for c in comps], min_num_heavy_ats, max_num_heavy_ats)
        if sc < best[0]:
            best = sc, (u, v), comps
        G.add_edge(u, v)
    
    _, (u, v), comps = best
    G.remove_edge(u, v)

    out: List[List[str]] = []
    for c in comps:
        out.extend(split(G, c, min_num_heavy_ats, max_num_heavy_ats))
    return out

def get_constrained_conncomps(
    adj: Dict[str, List[str]],
    min_num_heavy_ats: int,
    max_num_heavy_ats: int,
) -> List[List[str]]:
    '''
    Builds an undirected graph from `adj` and greedily removes edges so
    every connected component has size in [min_num_heavy_ats,
    max_num_heavy_ats]. At each step every internal edge of an oversized
    component is probed; the cut minimizing (total size violations,
    -smallest sub-component) is committed, then both sides recurse.
    '''
    G = nx.Graph()
    G.add_nodes_from(adj.keys())
    for u, nbrs in adj.items():
        for v in nbrs:
            if u != v:
                G.add_edge(u, v)

    result: List[List[str]] = []
    for c in nx.connected_components(G):
        result.extend(split(G, list(c), min_num_heavy_ats, max_num_heavy_ats))
    return result

def cluster(
    hcs:               List[str],
    df:                pd.DataFrame,
    max_num_heavy_ats: int = 3,
    min_num_heavy_ats: int = 2
) -> List[List[str]]:
    """
    Bond-contiguous grouping of HC-eligible heavy atoms into clusters.
    Uses BFS on the bond graph restricted to eligible atoms to find
    connected components, then splits oversized components.

    NW wrong: componnts split in non connected subcomponents, i need to remove edge and assess wheter 
    conn comps respect constraints number heavy <= max_num_heavy_ats and > 1
    """
    hcs_set = set(hcs)

    # Build adjacency restricted to hcs atoms
    #print(f'Starting HC clustering of residue::\n{df.head(50)}')
    #print(f'Detectd HC atoms: {hcs}')
    adj: Dict[str, List[str]] = {a: [] for a in hcs}
    for a in hcs:
        at = df[df['atom_name'] == a]
        if at.empty:
            log.warning('cluster: atom %s not found in residue', a)
            continue
        bonds_str = at['bonds'].values[0]
        if not isinstance(bonds_str, str):
            continue
        for b in bonds_str.split():
            if b in hcs_set and b != a:
                if b not in adj[a]:
                    adj[a].append(b)
                # symmetric
                if a not in adj.get(b, []):
                    adj.setdefault(b, []).append(a)

    components = get_constrained_conncomps(
            adj,
            min_num_heavy_ats=min_num_heavy_ats,
            max_num_heavy_ats=max_num_heavy_ats,
        )
    #print(components)
    return components

def add_hc_atoms(
    df:                   pd.DataFrame,
    max_abs_charge:       float = 0.27,
    max_num_heavy_ats:    int   = 3,
    min_bonds_to_charged: int   = 3,
    min_bonds_to_polar:   int   = 1,
    min_abs_charge:       float = 0.75,   # SB-rule threshold (charged)
    max_DA_charge:        float = -0.3,   # HBond donor/acceptor thr → polar (heavy)
    min_H_charge:         float =  0.3,   # Hbond H charge
) -> pd.DataFrame:
    """
    Identify hydrophobic-cluster atoms.  A heavy atom is HC iff:
      * |its own charge| <= max_abs_charge (not itself polar/charged), and
      * BFS distance to any SB-active charged anchor > min_bonds_to_charged
        (anchors = same atoms add_sb_atoms would label, via local-q rule), and
      * BFS distance to any polar heavy atom        > min_bonds_to_polar
        (polar = |charge| >= min(|max_DA_charge|, min_H_charge), i.e. anything
         that could donate/accept in add_hbonds).

    Polar set is broader and uses the smaller bond-distance threshold; SB anchors
    are sparser and use the stricter (larger) one.
    """
    log.debug('HC identification: max_abs_charge=%s max_num_heavy_ats=%s '
              'min_bonds_to_charged=%s min_bonds_to_polar=%s min_abs_charge=%s '
              'max_DA_charge=%s min_H_charge=%s',
              max_abs_charge, max_num_heavy_ats,
              min_bonds_to_charged, min_bonds_to_polar,
              min_abs_charge, max_DA_charge, min_H_charge)

    polar_thr = min(abs(max_DA_charge), abs(min_H_charge))
    df['HC_cluster'] = ''
    heavy = df[(round(df['mass']) != 1) & (~df['bonds'].isna())].copy()
    myresids = sorted(heavy['Resindex'].unique().tolist())

    for count, resix in enumerate(myresids):
        indf = heavy[heavy['Resindex'] == resix]
        ent  = np.unique(indf['Entity'].values)
        log.debug('HC residue %s (Entity %s)', resix, ent)

        # context (residue + neighbours for polymers)
        if 'P' in ent:
            lo = myresids[max(count - 1, 0)]
            hi = myresids[min(count + 1, len(myresids) - 1)]
            heavy_ctx = heavy[(heavy['Resindex'] >= lo) & (heavy['Resindex'] <= hi)].copy()
            ctx_full  = df[(df['Resindex']    >= lo) & (df['Resindex']    <= hi)].copy()
        else:
            heavy_ctx = heavy[heavy['Resindex'] == resix].copy()
            ctx_full  = df[df['Resindex']      == resix].copy()

        # polar atoms (heavy, |q| >= polar_thr) in context
        polar_mask  = abs(heavy_ctx['charge']) >= polar_thr
        polar_atoms = [(heavy_ctx.at[k, 'atom_name'], int(heavy_ctx.at[k, 'Resindex']))
                       for k in heavy_ctx[polar_mask].index]

        # SB-active charged anchors per residue (uses full atoms incl. H so
        # terminal H charges enter local_q correctly via get_terminal_tree)
        charged_atoms: List[Tuple[str, int]] = []
        for _ri, sub_full in ctx_full.groupby('Resindex'):
            charged_atoms.extend(sb_charged_atoms(sub_full, min_abs_charge))

        log.debug('HC   ctx: %d polar, %d charged anchors',
                  len(polar_atoms), len(charged_atoms))

        # candidates
        true_hcs: List[Tuple[str, int]] = []
        for _, row in indf.iterrows():
            if abs(row['charge']) > max_abs_charge:
                continue
            at = (row['atom_name'], int(row['Resindex']))
            ok = True
            if charged_atoms:
                d = below_max_num_bonds(at, charged_atoms, heavy_ctx,
                                        max_bonds=min_bonds_to_charged)
                if d <= min_bonds_to_charged:
                    log.debug('HC   %s: %d bonds to charged → reject', at[0], d)
                    ok = False
            if ok and polar_atoms:
                d = below_max_num_bonds(at, polar_atoms, heavy_ctx,
                                        max_bonds=min_bonds_to_polar)
                if d <= min_bonds_to_polar:
                    log.debug('HC   %s: %d bonds to polar → reject', at[0], d)
                    ok = False
            if ok:
                true_hcs.append(at)

        if not true_hcs:
            continue

        hcs    = list({a[0] for a in true_hcs})
        act_df = df[df['Resindex'] == resix]
        hc_clusters = (
            cluster(hcs, indf, max_num_heavy_ats=max_num_heavy_ats)
            if len(hcs) > 1 else [hcs]
        )
        log.debug('HC   clusters: %s', hc_clusters)

        for ci, hc_group in enumerate(hc_clusters, start=1):
            extended = list(hc_group)
            for an in hc_group:
                an_row = act_df[act_df['atom_name'] == an]
                if an_row.empty:
                    continue
                bonds_str = an_row['bonds'].values[0]
                if not isinstance(bonds_str, str):
                    continue
                hs = act_df[
                    act_df['atom_name'].isin(bonds_str.split())
                    & (round(act_df['mass']) == 1)
                ]['atom_name'].tolist()
                extended.extend(hs)

            cluster_atoms = act_df[act_df['atom_name'].isin(extended)]
            df.loc[cluster_atoms.index, 'HC_cluster'] = f'HC{ci}'
            for aix in cluster_atoms.index:
                if 'HC' not in str(df.at[aix, 'q_interactors']):
                    df.loc[aix, 'q_interactors'] += ' HC'

    return df

def add_sb_atoms(
    df:             pd.DataFrame,
    min_abs_charge: float = 0.75,
) -> pd.DataFrame:
    """
    Identify salt-bridge atoms using recursive local charge.
    
    For each non-terminal heavy atom, local_charge = its own charge +
    charges of all recursively terminal atoms (heavy or H) hanging off it.
    
    If |local_charge| >= min_abs_charge:
      - positive → Basic (label H atoms with max charge, 
        or heavy atom if no H)
      - negative → Acidic (label heavy atoms with min charge 
        among the group)
    """
    log.debug('salt-bridge identification: min_abs_charge=%s', min_abs_charge)

    for resix, resdf in df.groupby('Resindex'):
        log.debug('SB residue %s', resix)
        root_terminals = get_terminal_tree(resdf)

        # Also handle single-atom residues (ions) with no terminals
        heavy_in_res = resdf[round(resdf['mass']) != 1]

        # Atoms that are terminal somewhere — skip as roots
        all_terminals = set()
        for terms in root_terminals.values():
            all_terminals.update(terms)

        for ix, row in heavy_in_res.iterrows():
            if ix in all_terminals:
                continue  # this atom's charge is counted in its root

            aname = row['atom_name']
            own_q = row['charge']

            term_ixs = root_terminals.get(ix, [])
            term_q = df.loc[term_ixs, 'charge'].sum() if term_ixs else 0.0
            local_q = own_q + term_q

            if abs(local_q) < min_abs_charge:
                continue

            # Collect the full group: root + all terminals
            group_ixs = [ix] + term_ixs
            group_df = df.loc[group_ixs]

            if local_q > 0:
                # Basic: prefer H atoms with max charge
                h_atoms = group_df[round(group_df['mass']) == 1]
                if h_atoms.empty:
                    if 'Basic' not in df.at[ix, 'q_interactors']:
                        df.loc[ix, 'q_interactors'] += ' Basic'
                    log.debug('SB   %s: local_q=%.3f -> Basic (no H)', aname, local_q)
                else:
                    max_hq = h_atoms['charge'].max()
                    best_h = h_atoms[h_atoms['charge'] == max_hq]
                    for hix in best_h.index:
                        if 'Basic' not in df.at[hix, 'q_interactors']:
                            df.loc[hix, 'q_interactors'] += ' Basic'
                    log.debug('SB   %s: local_q=%.3f -> Basic H %s', aname, local_q, best_h["atom_name"].tolist())
            else:
                # Acidic: label heavy atoms with most negative charge
                heavy_group = group_df[round(group_df['mass']) != 1]
                min_q = heavy_group['charge'].min()
                most_neg = heavy_group[heavy_group['charge'] == min_q]
                for aix in most_neg.index:
                    if 'Acidic' not in df.at[aix, 'q_interactors']:
                        df.loc[aix, 'q_interactors'] += ' Acidic'
                log.debug('SB   %s: local_q=%.3f -> Acidic %s', aname, local_q, most_neg["atom_name"].tolist())

    return df

def add_hbonds(
    df:            pd.DataFrame,
    max_DA_charge: float = -0.3,
    min_H_charge:  float =  0.3,
    DA_masses:     Optional[List[int]] = None,
) -> pd.DataFrame:
    """
    Identify and label H-bond acceptor (HAcc), donor (HDon), and hydrogen
    (HBond) atoms.

    * **HAcc**: N/O/F atom with charge ≤ *max_DA_charge*.
      N excluded if it forms ≥ 4 bonds.
    * **HDon + HBond**: same heavy atom with at least one bound H atom with
      charge ≥ *min_H_charge*.

    Parameters
    ----------
    df:            Full topology DataFrame.
    out:           Log string (``None`` suppresses logging).
    max_DA_charge: Maximum (least negative) charge for D/A heavy atoms.
    min_H_charge:  Minimum H-atom charge for H-bonding participation.
    DA_masses:     Rounded masses of D/A elements (default [14, 16, 19] = N, O, F).

    Returns
    -------
    df:  Updated DataFrame.
    out: Updated log string (or ``None``).
    """
    if DA_masses is None:
        DA_masses = [14, 16, 19]

    log.debug('H-bond identification: max_DA_charge=%s min_H_charge=%s DA_masses=%s',
              max_DA_charge, min_H_charge, DA_masses)

    if 'q_interactors' not in df.columns:
        df['q_interactors'] = ''

    cols  = ['rtp_file', 'resname', 'atom_name', 'atom_type',
             'charge', 'chargegroup', 'bonds', 'mass',
             'Resindex', 'q_interactors']
    usedf = df[cols]

    for resix, resdf in usedf.groupby('Resindex'):
        log.debug('HBond residue %s', resix)
        DA_candidates = resdf[round(resdf['mass']).isin(DA_masses)]
        if DA_candidates.empty:
            continue

        for row in DA_candidates.itertuples():
            ix     = row.Index
            an     = row.atom_name
            charge = row.charge
            bonds  = row.bonds
            inter  = row.q_interactors
            atype  = row.atom_type if hasattr(row, 'atom_type') else ''

            first_letter = get_1st_letter([atype])[0][0][get_1st_letter([atype])[0][1]].upper() if atype else ''

            # ── Acceptor ───────────────────────────────────────────────
            if charge <= max_DA_charge and 'HAcc' not in inter:
                if first_letter != 'N':
                    df.loc[ix, 'q_interactors'] += ' HAcc'
                    log.debug('HBond   %s: HAcc', an)
                elif isinstance(bonds, str) and len(bonds.split()) < 4:
                    df.loc[ix, 'q_interactors'] += ' HAcc'
                    log.debug('HBond   %s: HAcc (N <4 bonds)', an)

            # ── Donor — search bound H across the whole residue ─────────
            if isinstance(bonds, str) and 'HDon' not in inter:
                bound_H = resdf[
                    resdf['atom_name'].isin(bonds.split()) &
                    (round(resdf['mass']) == 1) &
                    (resdf['charge'] >= min_H_charge)
                ]
                if not bound_H.empty:
                    df.loc[ix, 'q_interactors']            += ' HDon'
                    df.loc[bound_H.index, 'q_interactors'] += ' HBond'
                    log.debug('HBond   %s: HDon; H atoms %s', an, bound_H["atom_name"].tolist())
    return df

# ═══════════════════════════════════════════════════════════════════════════════
# Name-based interaction-type assignment helpers
# ═══════════════════════════════════════════════════════════════════════════════

def apply_rule(rule: List[Tuple[str, str]], inttype: str, df: pd.DataFrame,
               resid_net_charge: Optional[pd.Series] = None,
               charge_sign: int = 0, charge_thr: float = 0.0) -> pd.DataFrame:
    """Tag atoms matching ``rule`` with ``inttype``.

    If ``charge_sign`` is non-zero, the match is additionally gated by the net
    charge of the atom's residue (``resid_net_charge``, a per-atom Series of the
    residue's summed partial charge): only atoms whose residue carries net charge
    of the requested sign and magnitude are tagged. ``charge_sign=+1`` keeps
    residues with net charge > +``charge_thr`` (basic); ``charge_sign=-1`` keeps
    residues with net charge < -``charge_thr`` (acidic). This makes the
    salt-bridge roles protonation-state aware (neutral His ≈ 0 → not basic;
    protonated His ≈ +1 → basic) regardless of resname.
    """
    seldf = df.copy()
    log.debug('apply_rule: %s -> %s', rule, inttype)
    for r in rule:
        attr = r[0]
        vals = r[1]
        if isinstance(vals, str):
            vals = [vals]

        if attr.startswith('not '):
            seldf = seldf[~seldf[attr[4:]].isin(vals)]
        else:
            seldf = seldf[seldf[attr].isin(vals)]

    if charge_sign and resid_net_charge is not None:
        gate = resid_net_charge.reindex(seldf.index)
        if charge_sign > 0:
            seldf = seldf[gate > charge_thr]
        else:
            seldf = seldf[gate < -charge_thr]

    df, seldf = df.align(seldf, axis=1)
    df.loc[seldf.index, 'named_interactors'] += f' {inttype}'
    log.debug('apply_rule: matched %s atoms for %s', seldf.shape[0], inttype)
    return df

def fill_hydrogen_selections(
    top_df: pd.DataFrame,
    rounded_HB_masses: List[int] = [1],
) -> pd.DataFrame:
    """
    For named_interactors:in residues where HDon exists but HBond (hydrogens) are not present,
    mark H atoms bound to donor atom as HBond.

    Parameters
    ----------
    top_df: pd.DataFrame
        topology detailed with all information ( top_df of Interactors.add_inttypes())
    rounded_HB_masses: List[int] = [1]:
        List of all possibl masses (rounded to integers) of atoms that, if found to be bound to any Donor atom in the same
        residue, will be considered as H candidates for HBond  

    Returns
    -------
    top_df: pd.DataFrame
        Same DataFrame with 'HBond' entries filled where missing.
    """
    top_df['named_interactors'] = top_df['named_interactors'].fillna('')
    top_df['bonds'] = top_df['bonds'].fillna('')

    hdon_mask = top_df['named_interactors'].str.contains(r'\bHDon\b', regex=True)
    hdon_rows = top_df[hdon_mask][['Resindex', 'bonds']].copy()

    new_hbond_indices = []

    for _, row in hdon_rows.iterrows():
        bonded_names = [b for b in row['bonds'].split(' ') if b != '']
        if not bonded_names:
            continue

        rix = int(row['Resindex'])
        candidates = top_df[
            (top_df['Resindex'] == rix) &
            (top_df['atom_name'].isin(bonded_names)) &
            (top_df['mass'].round().astype(int).isin(rounded_HB_masses)) &
            (~top_df['named_interactors'].str.contains(r'\bHBond\b', regex=True))
        ]
        new_hbond_indices.extend(candidates.index.tolist())

    if new_hbond_indices:
        new_hbond_indices = list(set(new_hbond_indices))
        for ix in new_hbond_indices:
            # Tag as an HBond hydrogen but KEEP any HC tag: a donor-bound H such
            # as the backbone amide HN, or Trp HE1, is both an HBond hydrogen and
            # part of the hydrophobic side-chain COM in PyInteraph. Stripping HC
            # here removed HN from every hydrophobic residue and shifted every
            # side-chain COM (~11% HC error vs PyInteraph). Each contact type
            # uses its own tag, so carrying both is correct.
            top_df.at[ix, 'named_interactors'] = top_df.at[ix, 'named_interactors'] + ' HBond'

        log.debug('fill_hydrogen_selections: tagged %s H atoms as HBond (HC tag preserved)', len(new_hbond_indices))

    return top_df

# ═══════════════════════════════════════════════════════════════════════════════
# Called methods
# ═══════════════════════════════════════════════════════════════════════════════
 
def reconstruct_whole_topology(
    universe: mda.Universe,
    topologies: List[Path],
    resid_df: pd.DataFrame,
    inpath: Path,
    outpath: Path,
    nas: List[str]
) -> None:
    """
    Perform full topology reconstruction.
 
    Steps
    -----
    1. Fill bidirectional bond info for each residue type.
    2. Concatenate residue DataFrames into a single whole-topology table.
    3. Detect chain boundaries from ter entries.
    4. Stamp chains, classify polymer/monomer, add back-bonds, adjust ter
       chargegroups, assign atom indices, and write ``whole_topology.csv``.
 
    Parameters
    ----------
    universe:
        MDAnalysis Universe from the original topology file.
    topologies:
        Paths to per-residue ``<TopName>.csv`` files from
        :meth:`Topology_recognition.get_ff_parameters`.
    resid_to_top:
        DataFrame with columns ``['Resindex', 'topology']``.
    outpath:
        Destination folder; ``whole_topology.csv`` is written here.
 
    Returns
    -------
    None
    """
    log.info('Reconstructing whole topology')

    top_chains  = extract_chain_bounds(universe)
    top_df_dict = fill_bond_info(topologies, nas)

    resid_to_top = get_resid_to_top(resid_df, inpath)
    # Resindex → input-topology resname (HIS, …) for building selections later
    topology_resname_map = {}
    if 'topology_resname' in resid_to_top.columns:
        topology_resname_map = dict(zip(
            resid_to_top['Resindex'].astype(int),
            resid_to_top['topology_resname'].astype(str),
        ))
    wholetop, allres = build_wholetop(top_df_dict, resid_to_top, nas)

    wholetop['bonds'] = wholetop['bonds'].fillna('')

    mychains, wholetop = detect_chains(wholetop, allres)

    build_residue_variants(wholetop, mychains, universe, outpath,
                           topology_resname_map=topology_resname_map)

    tmp = wholetop.dropna(ignore_index=True)
    if wholetop.shape != tmp.shape:
        log.warning('NaNs introduced after chain detection: topology %s vs cleaned %s',
                    wholetop.shape, tmp.shape)

def get_charge_based_interactors(
    df:  pd.DataFrame,
    **kwargs,
) -> pd.DataFrame:
    """
    Assign interaction roles to all topology atoms based on their partial charge and input chrge thresholds.

    Calls :func:`add_hbonds` → :func:`add_sb_atoms` → :func:`add_hc_atoms`.

    Parameters
    ----------
    df:       Full topology DataFrame with a ``q_interactors`` column.
    out:      Log string (``None`` suppresses logging).
    **kwargs: 
    Override any default:
     - H-bond: 
         - ``max_DA_charge``(-0.3): maximum donor/acceptor charge
         - ``min_H_charge`` (0.3): minimum H charge
         - ``DA_masses``: list of rounded to int donor/acceptor masses (default [14, 16, 19] for N, O, F)
     - Salt-bridge: 
         - ``min_abs_charge`` (0.75): min absolute charge of a chargegroup 
     - HC:
         - ``max_abs_charge`` (0.27): max possible charge for HC atoms
         - ``max_num_heavy_ats`` (3): Bond-contiguous HC groups are split into sub-clusters of 3 heavy atoms each 
         - ``min_bonds_to_charged`` (3): minimum nuber of bonds separing a charged etom (w |charge| > max_abs_charge) 
                                         and an HC atom to include it in HC clusters 

    Returns
    -------
    df:  Updated DataFrame.
    out: Updated log string (or ``None``).
    """
    defaults = dict(
        max_DA_charge        = -0.3,
        min_H_charge         =  0.3,
        DA_masses            = [14, 16, 19],
        min_abs_charge       =  0.75,
        max_abs_charge       =  0.27,
        max_num_heavy_ats    =  3,
        min_bonds_to_charged =  3,
    )
    defaults.update(kwargs)

    log.debug('Running get_charge_based_interactors')

    df = add_hbonds(
        df,
        max_DA_charge = defaults['max_DA_charge'],
        min_H_charge  = defaults['min_H_charge'],
        DA_masses     = defaults['DA_masses'],
    )
    df = add_sb_atoms(df, min_abs_charge=defaults['min_abs_charge'])
    df = add_hc_atoms(
        df,
        max_abs_charge       = defaults['max_abs_charge'],
        max_num_heavy_ats    = defaults['max_num_heavy_ats'],
        min_bonds_to_charged = defaults['min_bonds_to_charged'],
        min_bonds_to_polar   = defaults['min_bonds_to_polar'],
        # shared semantics with the other two:
        min_abs_charge       = defaults['min_abs_charge'],
        max_DA_charge        = defaults['max_DA_charge'],
        min_H_charge         = defaults['min_H_charge'],
    )
    return df

def get_named_interactors(
    df:  pd.DataFrame,
    **kwargs,
) -> pd.DataFrame:
    """
    Assign PyInteraph-compatible interaction roles (protein-centric).

    Each role is evaluated against a list of selection rules.  A rule is
    either a single ``(attribute, values)`` tuple **or** a list of such tuples
    (all conditions AND-combined).  Negative conditions use attribute names
    prefixed with ``'not '``.

    Default selections mirror the original PyInteraph atom definitions for
    HDon, HAcc, HC, Basic, and Acidic.

    Parameters
    ----------
    df:       Full topology DataFrame with a ``named_interactors`` column.
    out:      Log string (``None`` suppresses logging).
    **kwargs: Override any default selection key.
    In particular: 
    'OG', 'NE2', 'OG1', 'NE', 'N', 'ND1', 'NZ',
                'NH1', 'NH2', 'OW', 'ND2', 'SG', 'OH2', 'NE1',
        (selection refined to exclude TRP NE)
        - H-Bonds:
         - HDon = [('atom_name', [
                'OH', 'N', 'OG', 'NE2', 'OG1', 'NE', 'N', 'ND1', 'NZ',
                'NH1', 'NH2', 'OW', 'ND2', 'SG', 'OH2']
                   ), 
                
                ],
         - HAcc = [('atom_name', [
                'OH', 'OG', 'OD1', 'OD2', 'OG1', 'O', 'ND1', 'NE2', 'OE2', 'OW', 'SG', 'OE1',
                'OH2', 'SD', 'O12', 'O13', 'O14', 'O11', 'O21', 'O22', 'O31', 'O32', 'O13A',
                'O13B'
            ])],
         - hb_masses = [1]
        - Hydrophobic: HC = [[
                ('resname',       ['ALA', 'VAL', 'LEU', 'ILE',
                                   'PHE', 'PRO', 'TRP', 'MET', 'IIL', 'GLM']),
                ('not atom_name', ['CA', 'C', 'O', 'N', 'H', 'H1', 'H2', 'H3',
                                   'O1', 'O2', 'OXT', 'OT1', 'OT2', 
                                   'O12', 'O13', 'O14', 'O11', 'O21', 'O22', 'O31', 'O32', 
                                   'O13A', 'O13B'
                               ]),       
         - Salt-bridges: 
            - Basic = [
                    [('resname',  ['GLY-NH3+', 'NH3+', 'NHT3+', 'PRO-NH2+']),
                         ('atom_name', ['N', 'H1', 'H2', 'H3', 'HT1', 'HT2', 'HT3'])],
                    [('resname',  ['HIS', 'HSD', 'HSE', 'HSP', 'HID', 'HIE', 'HIP']),
                     ('atom_name', ['CG', 'ND1', 'HD1', 'CD2', 'CE1', 'NE2', 'HE2'])],
                    [('resname',  ['LYS']),
                     ('atom_name', ['CE', ' NZ', 'HZ1', 'HZ2', 'HZ3'])],
                    [('resname',  ['ARG']),
                     ('atom_name', ['CD', 'NE', 'CZ', 'NH1', 'HH11', 'HH12',
                                    'NH2', 'HH21', 'HH22'])],
                    [('atom_name', 'CAL')],
                    [('atom_name', 'MG')]
                ]
            - Acidic = [
                    [('resname',  ['ASP']),
                     ('atom_name', ['CG', 'OD1', 'OD2'])],
                    [('resname',  ['GLU']),
                     ('atom_name', ['CD', 'OE1', 'OE2'])],
                    [('resname',  ['COO-']),
                         ('atom_name', ['C', 'O1', 'O2', 'O1T', 'O2T', 'OT1', 'OT2'])],
                
                ]
    Returns
    -------
    df:  Updated DataFrame.
    out: Updated log string (or ``None``).
    """
    # Named-role rule set. The caller passes the FULL definition via
    # `named_rules` (a {role: [rule, ...]} dict, e.g. from config) which
    # completely replaces the built-in defaults; when absent, the built-in
    # defaults are used. Each role's rules are converted to the internal grammar
    # and applied (AND within a rule, OR across a role's rules).
    effective: Dict[str, list] = kwargs.get('named_rules') or copy.deepcopy(
        NAMED_INTERACTOR_DEFAULTS
    )

    log.debug('Running get_named_interactors')

    df['named_interactors'] = ''

    # Per-charged-group net charge gate for the salt-bridge roles (Basic/Acidic).
    # Atoms are gated by the summed partial charge of their (Resindex, resname)
    # GROUP — not the whole residue. The forcefield/ter resname distinguishes a
    # terminus patch (NH3+, COO-) from the residue body, so each charged group is
    # gated independently: e.g. an N-terminal Asp splits into NH3+ (≈ +1 → basic)
    # and ASP (≈ -1 → acidic) instead of cancelling to ~0. Independent of
    # protonation/resname (neutral His ≈ 0 → not basic; protonated His ≈ +1 →
    # basic). Threshold via `min_residue_charge` (default 0.8); ≤ 0 disables.
    charge_thr = float(kwargs.get('min_residue_charge', 0.8) or 0.0)
    resid_net_charge = None
    if charge_thr > 0 and 'charge' in df.columns and 'Resindex' in df.columns:
        group_keys = ['Resindex', 'resname'] if 'resname' in df.columns else ['Resindex']
        resid_net_charge = df.groupby(group_keys)['charge'].transform('sum')

    for inttype in NAMED_INTERACTOR_ROLES:
        role_rules = effective.get(inttype, [])
        log.debug('named role: %s (%d rule(s))', inttype, len(role_rules))
        sign = {'Basic': 1, 'Acidic': -1}.get(inttype, 0)
        for myrule in _rules_from_config(role_rules):
            if sign and resid_net_charge is not None:
                df = apply_rule(myrule, inttype, df,
                                resid_net_charge=resid_net_charge,
                                charge_sign=sign, charge_thr=charge_thr)
            else:
                df = apply_rule(myrule, inttype, df)

    # NOTE: previously the HC tag was stripped from atoms that also got
    # HDon/HAcc/Basic/Acidic. That removed Met SD and Trp NE1 (acceptor/donor)
    # from the hydrophobic side-chain COM, which PyInteraph KEEPS — and it
    # accounted for a ~5% HC discrepancy on the validation set. The whole
    # side chain now contributes to the HC COM (an atom may carry several roles;
    # each contact type uses its own tag), matching PyInteraph (HC ~0.1% diff).

    df = fill_hydrogen_selections(df, effective.get('hb_masses', [1]))
    return df