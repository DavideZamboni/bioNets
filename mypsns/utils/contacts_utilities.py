## contacts_utilities.py
## Helper functions for the Contacts service module.
##
## Sections
## --------
## INIT          – pre-flight checks, input loading, defaults
## CORE          – distance evaluators, position getters, RAM guard
## WORKERS       – core_sb / core_hb / core_hc worker functions + dispatcher
## WRAPPERS      – salt_bridges / Hbonds / hydrophobics orchestrators
## POST          – merge_results, get_unique_indices
## SELECTION     – get_selection, generate_combinations

import datetime
import logging
import multiprocessing as mp
import re
import timeit
from logging.handlers import QueueHandler
from os import makedirs
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union
import math

import MDAnalysis as mda
import numpy as np
import pandas as pd
import psutil
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis as HBA
from MDAnalysis.lib.distances import capped_distance, self_capped_distance, distance_array, calc_bonds
import networkx as nx
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from itertools import combinations
from typing import Tuple, Set
from mypsns.utils.matrix_utilities import write_matrix, build_label_map, build_chain_map, resindex_from_id

log = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════════
# INIT – checks, input loading, defaults
# ═══════════════════════════════════════════════════════════════════════════════

def assemble_full_interactors(
    found_topologies_dir: Path,
    residue_index_path:   Path,
    atom_index_map_path:  Path,
    na_values: Optional[List[str]] = ["", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NULL", "NaN", "n/a", "nan", "null"]
) -> pd.DataFrame:
    """
    Reconstruct the flat per-atom interactors table from the three normalised
    files written by Interactors:

      found_topologies/<variant_id>.csv  – intrinsic + interactor columns
      residue_index.csv                  – Resindex, resnum, Chain, Entity, variant_id
      atom_index_map.csv                 – Resindex, atom_name, atom_index

    Each variant CSV is read once and reused for every Resindex that maps
    to it, so I/O scales with unique chemistries, not residue instances.
    The output is sorted by atom_index so it can be assigned positionally
    to a Universe via add_TopologyAttr (as check_input_top does).

    Returned columns mirror OLD whole_topology_interactors.csv plus
    the new resnum / variant_id columns.
    """
    res_idx = pd.read_csv(residue_index_path, header=0, index_col=False, na_values=na_values)
    aim     = pd.read_csv(atom_index_map_path, header=0, index_col=False, na_values=na_values)

    variant_cache: Dict[str, pd.DataFrame] = {}
    chunks:        List[pd.DataFrame]      = []

    for _, row in res_idx.iterrows():
        vid = str(row['variant_id'])
        if vid not in variant_cache:
            vpath = found_topologies_dir / f'{vid}.csv'
            if not vpath.exists():
                raise FileNotFoundError(
                    f'Variant CSV missing: {vpath} '
                    f'(referenced by Resindex {int(row["Resindex"])})'
                )
            variant_cache[vid] = pd.read_csv(vpath, header=0, index_col=False, na_values=na_values)

        vdf = variant_cache[vid].copy()
        vdf['Resindex']   = int(row['Resindex'])
        vdf['resnum']     = int(row['resnum'])
        vdf['Chain']      = row['Chain']
        vdf['Entity']     = row['Entity']
        vdf['variant_id'] = vid
        chunks.append(vdf)

    if not chunks:
        raise ValueError(
            f'residue_index.csv at {residue_index_path} is empty'
        )

    full = pd.concat(chunks, ignore_index=True)
    full = full.merge(aim, on=['Resindex', 'atom_name'], how='left')

    if full['atom_index'].isna().any():
        missing = full[full['atom_index'].isna()][['Resindex', 'atom_name']]
        raise ValueError(
            f'atom_index_map.csv missing entries for '
            f'{missing.shape[0]} (Resindex, atom_name) pairs:\n{missing.head(20)}'
        )

    full['atom_index'] = full['atom_index'].astype(int)
    full = full.sort_values('atom_index', kind='mergesort').reset_index(drop=True)
    return full


def startContacts(
    topology: str,
    inp: Path,
) -> Tuple[bool, str, str]:
    '''
    Validate Topology_recognition output.

    Returns (stop, message, topology_stem) where stop is True when a required
    file or folder is missing.
    '''
    name = Path(topology).stem
    log.info(
        f'\n\nContacts run starting at:\n'
        f'{datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")}\n'
        f'Topology: {name}{Path(topology).suffix}\n'
        f'Topology_recognition folder: {inp}\n'
    )

    if not inp.exists():
        return True, f'Topology_recognition folder not found: {inp}\n', name

    required = ['residue_index.csv', 'atom_index_map.csv']
    present  = {f.name for f in inp.iterdir() if f.is_file()}
    missing  = [r for r in required if r not in present]
    if missing:
        return (
            True,
            f'Missing required files in {inp.as_posix()}: {missing}\n'
            f'Present CSVs: {sorted(present)}\n',
            name,
        )

    found_dir = inp / 'found_topologies'
    if not found_dir.exists():
        return (
            True,
            f'found_topologies/ folder missing in {inp.as_posix()}.\n',
            name,
        )
    if not any(p.suffix == '.csv' for p in found_dir.iterdir() if p.is_file()):
        return (
            True,
            f'found_topologies/ contains no variant CSVs in {inp.as_posix()}.\n',
            name,
        )

    return False, '', name

def check_input_top(topology: str, inters: pd.DataFrame) -> mda.Universe:
    u = mda.Universe(topology)
    
    # Build lookup: for each atom, find its row in inters
    # Match by position (atom order should be identical)
    if len(u.atoms) != len(inters):
        log.warning(
            f'Atom count mismatch: Universe={len(u.atoms)}, '
            f'topology CSV={len(inters)}'
        )
    
    # Assign resids from Resindex column
    # Group by Resindex to get one resid per residue
    res_map = inters.groupby('Resindex').first().reset_index()
    new_resids = inters['Resindex'].values.tolist()
    
    # Assign chain IDs per atom
    new_chains = inters['Chain'].values.tolist()
    
    # Assign atom IDs
    new_ids = inters['atom_index'].astype(int).values.tolist()
    
    # Assign charges
    new_charges = inters['charge'].astype(float).values.tolist()
    
    # Stamp everything onto the Universe
    # per-atom attributes
    u.add_TopologyAttr('chainID', new_chains)
    
    # per-residue attributes
    residue_resids = []
    for res in u.residues:
        # Use the Resindex from the first atom of this residue
        first_atom_idx = res.atoms[0].index
        if first_atom_idx < len(inters):
            residue_resids.append(int(inters.iloc[first_atom_idx]['Resindex']))
        else:
            residue_resids.append(res.resid)
    
    u.add_TopologyAttr('resid', residue_resids)
    u.add_TopologyAttr('id', new_ids)
    u.add_TopologyAttr('charge', new_charges)
    
    return u

def resync_segid_to_chain(path: Path) -> None:
    out = []
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                chain = line[21]
                line  = line[:72] + chain.ljust(4) + line[76:]
            out.append(line)
    with open(path, "w") as f:
        f.writelines(out)


def load_inputs(
    inp: Path,
    outp: Path,
    name: str,
    topology: str,
    trajectory: str,
    na_values: Optional[List[str]] = ["", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NULL", "NaN", "n/a", "nan", "null"]
) -> Tuple[pd.DataFrame, Path, Path]:
    '''
    Load the interaction-topology table and ensure a corrected topology PDB exists.
    '''
    inters = assemble_full_interactors(
        found_topologies_dir = inp / 'found_topologies',
        residue_index_path   = inp / 'residue_index.csv',
        atom_index_map_path  = inp / 'atom_index_map.csv',
        na_values=na_values
    )
    top = outp / f'{name}_correct.pdb'

    existing_pdbs = [f.name for f in outp.iterdir() if f.suffix == '.pdb' and f.is_file()]
    log.info(f'Checking for {top.name} in {outp}; found: {existing_pdbs}')

    if top.name not in existing_pdbs:
        log.info(f'Generating corrected topology: {top.name} → {outp}\n')
        u = check_input_top(topology, inters)
        u.atoms.write(str(top))
        resync_segid_to_chain(top)
        

    log.info('Topology loaded.\n')

    trj = Path(trajectory)
    if not (trj.exists()):
        log.error("trajectory file provided does not exists")
        exit
    
    if not (trj.is_file()):
        log.error("trajectory file provided is not a file")
        exit
    return inters, top, trj

def init_defaults(poly_chains: List[str], found_chains: List[str], log_queue: mp.Queue) -> dict:
    '''
    Return the default parameter dictionary for a Contacts run.
    '''
    return dict(
        nc                = 2,
        ram_thresh        = 80,
        nframes           = 500,
        beg_frm           = 0,
        end_frm           = -1,
        interactors_type  = ['q_interactors', 'named_interactors'], #
        targets           = ['HBonds', 'Salt_bridges', 'HCs'], # , ],#
        subject_chains    = poly_chains,
        surrounding_chains= found_chains,
        d_a_dist          = 4.5, # should be: a_h_dist = 3.5 NW PyInteraph2: 3
        angle             = 150, # angle = 120 NW PyInteraph2: 150
        sb_cutoff         = 3.5, # 4.5, at least 1 pair ()
        hc_cutoff         = 5.5, # for GLM, POPC
        log_queue         = log_queue,
        matrix_type       = ['persistency', 'energy'],
    )

# ═══════════════════════════════════════════════════════════════════════════════
# SELECTION – get_selection, generate_combinations
# ═══════════════════════════════════════════════════════════════════════════════

def _ids_to_sel(sub_df: pd.DataFrame, chainID: str) -> Tuple[Optional[str], int]:
    """Build a compact, readable selection grouped by (resname, atom-name set):
        (chainID X) and (
            (resname A and resid R1 R2 ... and name N1 N2)
            or (resname B and resid R3 R4 ... and name M1 M2)
            ...
        )
    Residues sharing both their resname and the exact set of tagged atom
    names collapse into a single clause. Assumes atom names are unique
    within a residue (standard protein/lipid topologies).

    Returns (selection_string, expected_atom_count). The string is None and
    the count is 0 when there are no atoms to select.
    """
    if 'atom_name' not in sub_df.columns or 'Resindex' not in sub_df.columns:
        return None, 0
    sub_df = sub_df.dropna(subset=['atom_name', 'Resindex'])
    if sub_df.empty:
        return None, 0
    has_resname = 'resname' in sub_df.columns

    # Group residues by (resname, frozenset of atom names)
    groups: Dict[Tuple[str, frozenset], List[int]] = {}
    for resindex, gdf in sub_df.groupby('Resindex'):
        names = frozenset(gdf['atom_name'].astype(str).unique())
        if not names:
            continue
        rn = ''
        if has_resname:
            rn_series = gdf['resname'].dropna().astype(str)
            if not rn_series.empty:
                rn = rn_series.iloc[0]
        groups.setdefault((rn, names), []).append(int(resindex))

    clauses: List[str] = []
    for (rn, names), resids in groups.items():
        names_str  = ' '.join(sorted(names))
        resids_str = ' '.join(str(r) for r in sorted(resids))
        if rn:
            clauses.append(
                f'(resname {rn} and resid {resids_str} and name {names_str})'
            )
        else:
            clauses.append(f'(resid {resids_str} and name {names_str})')

    if not clauses:
        return None, 0
    return f'(chainID {chainID}) and ({" or ".join(clauses)})', len(sub_df)

def get_selection(
    chain_df: pd.DataFrame,
    inttype:  str,
    role:     str,
    chainID:  str,
) -> Optional[Union[Tuple[str, int], Dict[str, Tuple[str, int]]]]:
    """
    Build MDAnalysis selection string(s) for atoms matching a role.

    HC role → {cluster_label: (sel_str, expected_count)}
    Other roles → (sel_str, expected_count)
    Returns None when no atoms match (so build_selections' `if sel is not None`
    drops the role cleanly and the wrappers' .get(role, (None, 0)) default
    kicks in).
    """
    if inttype not in chain_df.columns:
        return None

    mask = chain_df[inttype].astype(str).str.contains(role, na=False)
    matched = chain_df[mask]
    if matched.empty:
        return None

    # HC q_interactors → per-cluster, per-resid selections
    if role == 'HC' and inttype == 'q_interactors' and 'HC_cluster' in matched.columns:
        cluster_seles: Dict[str, Tuple[str, int]] = {}
        for ri, resdf in matched.groupby(by='Resindex'):
            for clust_label, clust_df in resdf.groupby('HC_cluster'):
                if pd.isna(clust_label):
                    continue
                sel_str, expected = _ids_to_sel(clust_df, chainID=chainID)
                if sel_str is not None:
                    cluster_seles[f'{clust_label}_{ri}'] = (sel_str, expected)
        return cluster_seles or None

    if role == 'HC' and inttype == 'named_interactors':
        cluster_seles = {}
        for ri, resdf in matched.groupby(by='Resindex'):
            sel_str, expected = _ids_to_sel(resdf, chainID=chainID)
            if sel_str is not None:
                cluster_seles[f'HC_{ri}'] = (sel_str, expected)
        return cluster_seles or None

    # Single selection
    sel_str, expected = _ids_to_sel(matched, chainID=chainID)
    return (sel_str, expected) if sel_str is not None else None

def generate_combinations(
    subjects: List[str],
    others:   List[str],
) -> List[List[str]]:
    '''
    Build chain-pair combinations for contact computation.

    For each subject chain, produce [subject, *remaining_others].
    '''
    #print(subjects, others)
    for s in subjects:
        if s not in others:
            others.append(s)

    indices        = [others.index(s) for s in subjects]
    sorted_indices = sorted(set(indices))
    return [others[i:] for i in sorted_indices]

# ═══════════════════════════════════════════════════════════════════════════════
# CORE – selection checker, distance evaluators, position getters, RAM guard
# ═══════════════════════════════════════════════════════════════════════════════

def check_selection(u: mda.Universe, selection: str,
                    expected: Optional[int] = None) -> bool:
    """Verify a selection returns atoms (optionally a specific count)."""
    n = len(u.select_atoms(selection))
    if n == 0:
        log.warning('Selection matched 0 atoms: %s', selection[:200])
        return False
    if expected is not None and n != expected:
        log.warning('Selection matched %d atoms, expected %d: %s',
                    n, expected, selection[:200])
        return False
    return True

def split_frames(u: mda.Universe, **kwargs) -> List[List[int]]:
    '''
    Partition trajectory frames into nc roughly equal sublists.
    '''
    defaults: dict = {}
    defaults.update(kwargs)

    allframes  = [ts.frame for ts in u.trajectory]
    myend      = allframes[-1] if defaults['end_frm'] == -1 else defaults['end_frm']
    to_analyze = [fr for fr in allframes if defaults['beg_frm'] <= fr <= myend]
    nc         = defaults['nc']

    if len(to_analyze) <= nc:
        return [
            [to_analyze[i]] if i < len(to_analyze) else []
            for i in range(nc)
        ]

    split_size = len(to_analyze) // nc
    remainder  = len(to_analyze) % nc
    split_lists: List[List[int]] = []
    last = 0

    for i in range(nc):
        if remainder > 0:
            split_lists.append(to_analyze[last:last + split_size + 1])
            remainder -= 1
            last      += split_size + 1
        else:
            split_lists.append(to_analyze[last:last + split_size])
            last += split_size

    return split_lists

def check_ram_usage(threshold: int = 80) -> bool:
    '''Return True when current RAM usage meets or exceeds the threshold.'''
    return psutil.virtual_memory().percent >= threshold

def atom_distance_evaluator(
    mydict:     dict,
    frame:      int,
    sel1:       mda.AtomGroup,
    cutoff:     float,
    dimensions: np.ndarray,
    sel2:       Optional[mda.AtomGroup] = None,
) -> dict:
    '''
    Evaluate pairwise distances between atom selections and append results.
    Pairs whose two atoms belong to the same residue are discarded.
    '''
    if sel2 is None:
        pairs, distances = self_capped_distance(
            sel1.positions,
            max_cutoff=cutoff,
            return_distances=True,
            box=dimensions,
        )
        for k, (i, j) in enumerate(pairs):
            r1, r2 = sel1[i].resid, sel1[j].resid
            if r1 != r2:
                mydict['frame'].append(frame)
                mydict['id_1'].append(f'{sel1[i].index}_r{r1}')
                mydict['id_2'].append(f'{sel1[j].index}_r{r2}')
                mydict['distance'].append(distances[k])
    else:
        pairs, distances = capped_distance(
            sel1.positions,
            sel2.positions,
            max_cutoff=cutoff,
            return_distances=True,
            box=dimensions,
        )
        for k, (i, j) in enumerate(pairs):
            r1, r2 = sel1[i].resid, sel2[j].resid
            if r1 != r2:
                mydict['frame'].append(frame)
                mydict['id_1'].append(f'{sel1[i].index}_r{r1}')
                mydict['id_2'].append(f'{sel2[j].index}_r{r2}')
                mydict['distance'].append(distances[k])

    return mydict

def points_distance_evaluator(
    mydict:     dict,
    frame:      int,
    pos1:       np.ndarray,
    cutoff:     float,
    dimensions: np.ndarray,
    names:      List[str],
    pos2:       Optional[np.ndarray] = None,
    names2:     Optional[List[str]]  = None,
) -> dict:
    '''
    Evaluate distances between sets of representative positions (e.g. CoMs).
    '''
    if pos2 is None:
        pairs, distances = self_capped_distance(
            np.array(pos1),
            max_cutoff=cutoff,
            return_distances=True,
            box=dimensions,
        )
        for k, (i, j) in enumerate(pairs):
            if names[i] != names[j]:
                mydict['frame'].append(frame)
                mydict['id_1'].append(names[i])
                mydict['id_2'].append(names[j])
                mydict['distance'].append(distances[k])
    else:
        pairs, distances = capped_distance(
            np.array(pos1),
            np.array(pos2),
            max_cutoff=cutoff,
            return_distances=True,
            box=dimensions,
        )
        for k, (i, j) in enumerate(pairs):
            if names[i] != names2[j]:
                mydict['frame'].append(frame)
                mydict['id_1'].append(names[i])
                mydict['id_2'].append(names2[j])
                mydict['distance'].append(distances[k])

    return mydict

def collect_result(result: Optional[Tuple[str, int, int]]) -> None:
    '''Callback for pool.apply_async; logs completion or skips empty results.'''
    if result is None:
        log.warning('Worker returned None – result skipped.')
    else:
        log.info(
            f'{result[0]} worker finished: frames {result[1]}–{result[2]}.'
        )

# ═══════════════════════════════════════════════════════════════════════════════
# WORKERS – per-interaction-type worker functions + dispatcher
# ═══════════════════════════════════════════════════════════════════════════════

def worker_init(log_queue: mp.Queue) -> None:
    '''Per-worker initialiser: route all log records to the shared queue.'''
    qh   = QueueHandler(log_queue)
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    root.addHandler(qh)

def core_sb(
    top:              str,
    trj:              str,
    frames:           List[int],
    acidic_selection: str,
    basic_selection:  str,
    cutoff:           float,
    outpath:          Path,
    outname:          str,
    threshold:        int,
) -> Tuple[str, int, int]:
    '''
    Worker: compute salt-bridge distances for a subset of frames.
    '''
    log.info(f'SB worker: frames {frames[0]}–{frames[-1]}')

    u = mda.Universe(top, trj)
    acidic_ats  = u.select_atoms(acidic_selection)
    basic_ats   = u.select_atoms(basic_selection)
        
    salt_dict   = {'frame': [], 'id_1': [], 'id_2': [], 'distance': []}
    last_frame  = frames[0]

    for ts in u.trajectory[frames]:
        if check_ram_usage(threshold):
            pd.DataFrame(salt_dict).to_csv(
                outpath / f'{outname}_upto_{last_frame}.csv', header=True, index=False
            )
            salt_dict = {'frame': [], 'id_1': [], 'id_2': [], 'distance': []}

        salt_dict  = atom_distance_evaluator(
            salt_dict, ts.frame, acidic_ats, cutoff, u.dimensions, sel2=basic_ats
        )
        last_frame = ts.frame

    pd.DataFrame(salt_dict).to_csv(
        outpath / f'{outname}_upto_{last_frame}.csv', header=True, index=False
    )
    log.info(f'SB worker done: frames {frames[0]}–{frames[-1]}')
    return 'Salt bridges', frames[0], frames[-1]

def estimate_hb_frame_window(
    u, donors: str, acceptors: str, hydrogens: Optional[str],
    distc: float, angc: float, nframes: int, threshold: int, nc: int,
) -> int:
    '''
    Estimate a per-worker HBA batch size (frames) that keeps `nc` concurrent
    workers under `threshold`% system RAM.

    Runs ONE small HBA batch to measure result memory per frame, fits that to
    the current RAM headroom, then splits the budget across nc workers. Computed
    once (in Hbonds) instead of re-probed inside every worker.
    '''
    n_traj  = len(u.trajectory)
    probe_n = max(1, min(nframes, n_traj))
    probe   = run_batch_hba(u, donors, acceptors, hydrogens, distc, angc, 0, probe_n)
    bytes_per_frame = max(1.0, probe.memory_usage(deep=True).sum() / probe_n)

    vm     = psutil.virtual_memory()
    used   = vm.total - vm.available
    budget = vm.total * (threshold / 100.0) - used          # headroom before threshold
    if budget <= bytes_per_frame * nframes:
        log.warning('HB window: little RAM headroom – falling back to nframes=%d', nframes)
        return nframes

    per_worker = int((budget / max(1, nc)) / bytes_per_frame)
    window = max(nframes, min(per_worker, n_traj))
    log.info('HB window estimate: %.1f KB/frame, budget=%.2f GB, nc=%d -> %d frames/worker',
             bytes_per_frame / 1024, budget / 1e9, nc, window)
    return window

def run_batch_hba(u,
                  donors,
                  acceptors,
                  hydrogens,
                  distc,
                  angc,
                  start,
                  stop,                  
                  ) -> pd.DataFrame:
        
    log.info(f"run_batch_hba running with frames from {start} to {stop}")

    myhbs = HBA(u,             
        donors_sel         = donors,
        acceptors_sel      = acceptors,
        hydrogens_sel      = hydrogens,
        d_a_cutoff         = distc,
        d_h_a_angle_cutoff = angc,
        update_selections  = True
    )
    myhbs.run(start=start, stop=stop, verbose=False)

    cols = ['Frame', 'Donor_ix', 'Hydrogen_ix', 'Acceptor_ix', 'DA_Distance', 'HA_distance']
    if myhbs.results.hbonds.size == 0:
        log.warning("No results found in HBonds Analysis")
        return pd.DataFrame(columns=cols)
    
    pad_df = pd.DataFrame(
        myhbs.results.hbonds[:, :4].astype(int),
        columns=['Frame', 'Donor_ix', 'Hydrogen_ix', 'Acceptor_ix'],
    )
    pad_df['DA_Distance'] = myhbs.results.hbonds[:, 4]

    # hydrogen–acceptor distance (for energy); HBA doesn't report it → per-frame recompute
    ha = np.empty(len(pad_df), dtype=float)
    for f, sub in pad_df.groupby('Frame'):
        u.trajectory[int(f)]
        hpos = u.atoms[sub['Hydrogen_ix'].to_numpy()].positions
        apos = u.atoms[sub['Acceptor_ix'].to_numpy()].positions
        ha[sub.index.to_numpy()] = calc_bonds(hpos, apos, box=u.dimensions)
    pad_df['HA_distance'] = ha
    return pad_df

def core_hb(
    top:          str,
    trj:          str,
    frames:       List[int],
    acceptors:    str,
    donors:       str,
    distc:        float,
    angc:         float,
    outpath:      Path,
    outname:      str,
    hydrogens:    Optional[str],
    frame_window: int,
) -> Tuple[str, int, int]:
    '''
    Worker: compute H-bonds for a subset of frames using a FIXED batch window
    (sized once in Hbonds). No per-worker RAM probing.
    '''
    log.info(f'HB worker: frames {frames[0]}–{frames[-1]} (window={frame_window})')
    log.debug('donors=%s acceptors=%s hydrogens=%s d_a=%s angle=%s',
              donors, acceptors, hydrogens, distc, angc)

    u    = mda.Universe(top, trj)
    last = 0
    while last < len(frames):
        batch = frames[last: last + frame_window]
        if not batch:
            break
        start, stop = batch[0], batch[-1] + 1        # run() stop is exclusive
        log.info(f'HB batch frames {batch[0]}–{batch[-1]} (window={frame_window})')
        hb_df = run_batch_hba(u, donors, acceptors, hydrogens, distc, angc, start, stop)
        hb_df.to_csv(outpath / f'{outname}_upto_{batch[-1]}.csv', header=True, index=False)
        last += len(batch)

    log.info(f'HB worker done: frames {frames[0]}–{frames[-1]}')
    return 'H-bonds', frames[0], frames[-1]

def core_hc(
    top:       Union[str, Path],
    trj:       str,
    frames:    List[int],
    sel_dict:  Dict[str, str],
    cutoff:    float,
    outpath:   Path,
    outname:   str,
    threshold: int,
    sel2_dict: Optional[Dict[str, str]],
) -> Tuple[str, int, int]:
    '''
    Worker: compute hydrophobic-cluster distances for a subset of frames.
 
    Each entry in sel_dict / sel2_dict is one CoM unit:
      key   – label  (e.g. 'HC1_0', 'HC_760')
      value – MDAnalysis selection string
 
    For intrachain (sel2_dict is None):
        all-vs-all CoM distances within sel_dict.
    For interchain:
        sel_dict CoMs  vs  sel2_dict CoMs.
    '''
    log.info(f'HC worker: frames {frames[0]}–{frames[-1]}')
 
    u = mda.Universe(str(top), trj)
 
    # ── Pre-select atom groups once ──────────────────────────────────────────
    groups1: Dict[str, mda.AtomGroup] = {}
    for label, sel in sel_dict.items():
        grp = u.select_atoms(sel)
        if len(grp) > 0:
            groups1[label] = grp
 
    groups2: Optional[Dict[str, mda.AtomGroup]] = None
    if sel2_dict is not None:
        groups2 = {}
        for label, sel in sel2_dict.items():
            grp = u.select_atoms(sel)
            if len(grp) > 0:
                groups2[label] = grp
 
    if not groups1:
        log.warning(f'HC worker: no atoms matched sel_dict – writing empty CSV')
        pd.DataFrame({'frame': [], 'id_1': [], 'id_2': [], 'distance': []}).to_csv(
            outpath / f'{outname}_upto_{frames[-1]}.csv', header=True, index=False
        )
        return 'HC', frames[0], frames[-1]
 
    # ── Frame loop ───────────────────────────────────────────────────────────
    hc_dict: dict = {'frame': [], 'id_1': [], 'id_2': [], 'distance': []}
    last_frame = frames[0]

    for frame in frames:
        if check_ram_usage(threshold):
            pd.DataFrame(hc_dict).to_csv(
                outpath / f'{outname}_upto_{last_frame}.csv', header=True, index=False
            )
            hc_dict = {'frame': [], 'id_1': [], 'id_2': [], 'distance': []}
 
        u.trajectory[frame]
        dims = u.dimensions
 
        # Compute CoMs for set 1
        pos1: List[np.ndarray] = []
        names1: List[str] = []
        for label, grp in groups1.items():
            pos1.append(grp.center_of_mass())
            names1.append(label)
 
        if not pos1:
            continue
 
        if groups2 is None:
            # Intrachain: all-vs-all within set 1
            hc_dict = points_distance_evaluator(
                hc_dict, frame, np.array(pos1), cutoff, dims, names1
            )
        else:
            # Interchain: set 1 vs set 2
            pos2: List[np.ndarray] = []
            names2: List[str] = []
            for label, grp in groups2.items():
                pos2.append(grp.center_of_mass())
                names2.append(label)
            if pos2:
                hc_dict = points_distance_evaluator(
                    hc_dict, frame, np.array(pos1), cutoff, dims, names1,
                    np.array(pos2), names2
                )
        last_frame = frame

    pd.DataFrame(hc_dict).to_csv(
        outpath / f'{outname}_upto_{last_frame}.csv', header=True, index=False
    )
    log.info(f'HC worker done: frames {frames[0]}–{frames[-1]}')
    return 'HC', frames[0], frames[-1]

def multi_run_wrapper(args: tuple) -> Tuple:
    '''
    Dispatcher for parallel worker calls.
    Last element of args is a tag: 'sb' → core_sb, 'hb' → core_hb, 'hc' → core_hc.
    '''
    tag = args[-1]
    if tag == 'sb':
        log.info('starting core_sb')
        return core_sb(*args[:-1])
    elif tag == 'hb':
        log.info('starting core_hb')
        return core_hb(*args[:-1])
    elif tag == 'hc':
        log.info('starting core_hc')
        return core_hc(*args[:-1])

# ═══════════════════════════════════════════════════════════════════════════════
# WRAPPERS – salt_bridges / Hbonds / hydrophobics orchestrators
# ═══════════════════════════════════════════════════════════════════════════════

def salt_bridges(
    top:        str,
    trj:        str,
    selections: Dict[str, Dict[str, str]],
    outpath:    Path,
    frames:     List[List[int]],
    **kwargs,
) -> None:
    '''Orchestrate parallel salt-bridge computation across all chain combinations.'''
    tot_start = timeit.default_timer()
    defaults: dict = {}
    defaults.update(kwargs)

    cutoff    = float(defaults['sb_cutoff'])
    ram_thr   = int(defaults['ram_thresh'])
    nc        = defaults['nc']
    log_queue = defaults['log_queue']
    combos    = defaults['chain_combinations']
    outnames: List[str] = []

    log.info(f'salt_bridges: cutoff={cutoff} Å\n')

    for cs in combos:
        c1 = cs[0]
        sel_1_ac, ac1_expected = selections.get(c1, {}).get('Acidic', (None, 0))
        sel_1_bas, bas1_expected = selections.get(c1, {}).get('Basic', (None, 0))

        if sel_1_ac is None or sel_1_bas is None:
            log.info(f'  Chain {c1}: empty Acidic or Basic selection – intrachain skipped.\n')
        else:
            u = mda.Universe(top, trj)
            keep_running = check_selection(u, sel_1_ac, expected = ac1_expected)
            if not keep_running:
                break

            keep_running = check_selection(u, sel_1_bas, expected = bas1_expected)
            if not keep_running:
                break

            fname   = f'{c1}_intrachain_sb'
            outnames.append(fname)
            allargs = [
                (top, trj, mf, sel_1_ac, sel_1_bas, cutoff, outpath, fname, ram_thr, 'sb')
                for mf in frames
            ]
            t0   = timeit.default_timer()
            pool = mp.Pool(processes=nc, initializer=worker_init, initargs=(log_queue,))
            for args in allargs:
                pool.apply_async(multi_run_wrapper, (args,), callback=collect_result, error_callback=lambda e: log.error(f'HB worker crashed: {e}'))
            pool.close()
            pool.join()
            log.info(f'  Chain {c1} intrachain SB: {timeit.default_timer()-t0:.2f}s\n')

        for c2 in cs[1:]:
            sel_2_ac, ac2_expected  = selections.get(c2, {}).get('Acidic', (None, 0))
            sel_2_bas, bas2_expected = selections.get(c2, {}).get('Basic', (None, 0))
            t0 = timeit.default_timer()

            if sel_1_bas is not None and sel_2_ac is not None:
                keep_running = check_selection(u, sel_2_ac, expected = ac2_expected)
                if not keep_running:
                    break

                fname   = f'{c1}_bas_{c2}_ac_sb'
                outnames.append(fname)
                allargs = [
                    (top, trj, mf, sel_2_ac, sel_1_bas, cutoff, outpath, fname, ram_thr, 'sb')
                    for mf in frames
                ]
                pool = mp.Pool(processes=nc, initializer=worker_init, initargs=(log_queue,))
                for args in allargs:
                    pool.apply_async(multi_run_wrapper, (args,), callback=collect_result, error_callback=lambda e: log.error(f'SB worker crashed: {e}'))
                pool.close()
                pool.join()
            else:
                log.info(f'  Chains {c1}/{c2}: empty selection – {c1}bas×{c2}ac skipped.\n')

            if sel_1_ac is not None and sel_2_bas is not None:
                keep_running = check_selection(u, sel_2_bas,expected = bas2_expected)
                if not keep_running:
                    break
                
                fname   = f'{c1}_ac_{c2}_bas_sb'
                outnames.append(fname)
                allargs = [
                    (top, trj, mf, sel_1_ac, sel_2_bas, cutoff, outpath, fname, ram_thr, 'sb')
                    for mf in frames
                ]
                pool = mp.Pool(processes=nc, initializer=worker_init, initargs=(log_queue,))
                for args in allargs:
                    pool.apply_async(multi_run_wrapper, (args,), callback=collect_result, error_callback=lambda e: log.error(f'SB worker crashed: {e}'))
                pool.close()
                pool.join()
            else:
                log.info(f'  Chains {c1}/{c2}: empty selection – {c1}ac×{c2}bas skipped.\n')

            log.info(f'  Chains {c1}/{c2} interchain SB: {timeit.default_timer()-t0:.2f}s\n')

    merge_results(outpath, outnames)
    log.info(f'Total SB time: {timeit.default_timer()-tot_start:.2f}s\n')

def Hbonds(
    top:        str,
    trj:        str,
    selections: Dict[str, Dict[str, Tuple[str, int]]],
    outpath:    Path,
    frames:     List[List[int]],
    **kwargs,
) -> str:
    '''
    Orchestrate parallel H-bond computation across all chain combinations.

    '''
    log.info('starting HBonds computation')
    tot_start = timeit.default_timer()
    defaults: dict = {}
    defaults.update(kwargs)

    distc     = float(defaults['d_a_dist'])
    angc      = float(defaults['angle'])
    nframes   = int(defaults['nframes'])
    nc        = defaults['nc']
    log_queue = defaults['log_queue']
    combos    = defaults['chain_combinations']
    threshold = defaults['ram_thresh']
    outnames: List[str] = []
    log.info(f'Hbonds: d_a_dist={distc} Å, angle={angc}°, nframes={nframes}\n')

    # ── Adaptive frame window computation ──
    per_worker_window = nframes
    # ── Per-worker batch window, computed once on the LARGEST selection ──────
    per_worker_window = nframes
    u_probe = mda.Universe(top, trj)
    best = None     # (n_atoms, donors, acceptors, hydrogens, chain)
    for ch, roles in selections.items():
        d, _ = roles.get('HDon',  (None, 0))
        a, _ = roles.get('HAcc',  (None, 0))
        if not (d and a):
            continue
        h, _ = roles.get('HBond', (None, 0))    # may be None — fine, HBA accepts that
        n_at = len(u_probe.select_atoms(d)) + len(u_probe.select_atoms(a))
        if best is None or n_at > best[0]:
            best = (n_at, d, a, h, ch)

    if best is not None:
        _, d, a, h, ch = best
        per_worker_window = estimate_hb_frame_window(
            u_probe, d, a, h, distc, angc, nframes, threshold, nc)
        log.info('HB per-worker frame window (probed on largest selection: '
                 'chain %s, %d donor+acceptor atoms): %d', ch, best[0], per_worker_window)
    del u_probe
    
    for cs in combos:
        c1         = cs[0]
        sel_1_hacc, hacc1_expected = selections.get(c1, {}).get('HAcc', (None, 0))
        sel_1_hdon, hdon1_expected = selections.get(c1, {}).get('HDon', (None, 0))
        sel_1_hs, hs1_expected   = selections.get(c1, {}).get('HBond', (None, 0))

        log.info(f'Chain {c1}:\n\tAcceptors: {sel_1_hacc}\n\tDonors: {sel_1_hdon}')

        if sel_1_hacc is None or sel_1_hdon is None:
            log.info('skipped')
            log.info(f'  Chain {c1}: empty HAcc or HDon – intrachain skipped.\n')
        else:
            u = mda.Universe(top, trj)
            keep_running = check_selection(u, sel_1_hacc, expected = hacc1_expected)
            if not keep_running:
                break
            keep_running = check_selection(u, sel_1_hdon, expected = hdon1_expected)
            if not keep_running:
                break
            if (sel_1_hs is not None):
                keep_running = check_selection(u, sel_1_hs, expected = hs1_expected)
                if not keep_running:
                    break

            fname   = f'{c1}_intrachain_Hbonds'
            outnames.append(fname)

            allargs = [
                (top, trj, mf, sel_1_hacc, sel_1_hdon, distc, angc,
                 outpath, fname, sel_1_hs, per_worker_window, 'hb')
                for mf in frames
            ]

            t0   = timeit.default_timer()
            pool = mp.Pool(processes=nc, initializer=worker_init, initargs=(log_queue,))
            for args in allargs:
                pool.apply_async(multi_run_wrapper, (args,), callback=collect_result, error_callback=lambda e: log.error(f'HB worker crashed: {e}'))
            pool.close()
            pool.join()
            log.info(f'  Chain {c1} intrachain HB: {timeit.default_timer()-t0:.2f}s\n')

        for c2 in cs[1:]:
            sel_2_hacc, hacc2_expected = selections.get(c2, {}).get('HAcc', (None, 0))
            sel_2_hdon, hdon2_expected = selections.get(c2, {}).get('HDon', (None, 0))
            sel_2_hs, hs2_expected   = selections.get(c2, {}).get('HBond', (None, 0))
            t0 = timeit.default_timer()

            if sel_1_hdon is not None and sel_2_hacc is not None:
                keep_running = check_selection(u, sel_2_hacc, expected = hacc2_expected)
                if not keep_running:
                    break

                fname   = f'{c1}_HDon_{c2}_HAcc'
                outnames.append(fname)
                allargs = [
                    (top, trj, mf, sel_2_hacc, sel_1_hdon, distc, angc,
                     outpath, fname, sel_1_hs, per_worker_window, 'hb')
                    for mf in frames
                ]
                pool = mp.Pool(processes=nc, initializer=worker_init, initargs=(log_queue,))
                for args in allargs:
                    pool.apply_async(multi_run_wrapper, (args,), callback=collect_result, error_callback=lambda e: log.error(f'HB worker crashed: {e}'))
                pool.close()
                pool.join()
            else:
                log.info(f'  Chains {c1}/{c2}: empty selection – {c1}Don×{c2}Acc skipped.\n')

            if sel_1_hacc is not None and sel_2_hdon is not None:
                keep_running = check_selection(u, sel_2_hdon, expected = hdon2_expected)
                if not keep_running:
                    break
                if (sel_2_hs is not None):
                    keep_running = check_selection(u, sel_2_hs, expected = hs2_expected)
                    if not keep_running:
                        break

                fname   = f'{c1}_HAcc_{c2}_HDon'
                outnames.append(fname)
                
                allargs = [
                    (top, trj, mf, sel_1_hacc, sel_2_hdon, distc, angc,
                     outpath, fname, sel_2_hs, per_worker_window, 'hb')
                    for mf in frames
                ]
                pool = mp.Pool(processes=nc, initializer=worker_init, initargs=(log_queue,))
                for args in allargs:
                    pool.apply_async(multi_run_wrapper, (args,), callback=collect_result, error_callback=lambda e: log.error(f'HB worker crashed: {e}'))
                pool.close()
                pool.join()
            else:
                log.info(f'  Chains {c1}/{c2}: empty selection – {c1}Acc×{c2}Don skipped.\n')

            log.info(f'  Chains {c1}/{c2} interchain HB: {timeit.default_timer()-t0:.2f}s\n')

    merge_results(outpath, outnames)
    log.info(f'Total HB time: {timeit.default_timer()-tot_start:.2f}s\n')

def hydrophobics(
    top:        str,
    trj:        str,
    selections: Dict[str, Dict[str, Union[str, Dict[str, str]]]],
    outpath:    Path,
    frames:     List[List[int]],
    **kwargs,
) -> str:
    '''
    Orchestrate parallel hydrophobic-cluster computation.
 
    HC selections from get_selection are dicts:
      q_interactors  → {cluster_resindex_label: selection_string}
      named_interactors → {HC_resindex_label: selection_string}
 
    Each dict entry is one CoM unit. dict passed
    to core_hc.
    '''
    log.info('starting HC computation')

    tot_start = timeit.default_timer()
    defaults: dict = {}
    defaults.update(kwargs)

    cutoff       = float(defaults['hc_cutoff'])
    ram_thr      = int(defaults['ram_thresh'])
    nc           = int(defaults['nc'])
    log_queue    = defaults['log_queue']
    combos       = defaults['chain_combinations']
    outnames: List[str] = []

    log.info(f'hydrophobics: cutoff={cutoff} Å\n')

    for cs in combos:
        c1 = cs[0]
        sel_1_hc = selections.get(c1, {}).get('HC', {})        # dict default
        log.info(f'Chain {c1}:\n\tHCs: {list(sel_1_hc)}')

        if not sel_1_hc:
            log.info(f'  Chain {c1}: no HC selection – skipped.\n')
            continue

        u = mda.Universe(top, trj)

        # Validate each cluster, build a string-only dict for core_hc
        sel_1_hc_str: Dict[str, str] = {}
        n_atoms_total = 0
        keep_running  = True
        for label, (sel_str, expected) in sel_1_hc.items():
            if sel_str is None:
                continue
            if not check_selection(u, sel_str, expected=expected):
                keep_running = False
                break
            sel_1_hc_str[label] = sel_str
            n_atoms_total += expected
        if not keep_running:
            break
        log.info(f'Chain {c1}: {len(sel_1_hc_str)} HC groups, {n_atoms_total} atoms total')

        # Intrachain — pass the string-only dict
        fname = f'{c1}_intrachain_HC'
        outnames.append(fname)
        allargs = [
            (top, trj, mf, sel_1_hc_str, cutoff, outpath, fname, ram_thr, None, 'hc')
            for mf in frames
        ]

        t0   = timeit.default_timer()
        pool = mp.Pool(processes=nc, initializer=worker_init, initargs=(log_queue,))
        for args in allargs:
            pool.apply_async(
                multi_run_wrapper, (args,),
                callback=collect_result,
                error_callback=lambda e: log.error(f'HC worker crashed: {e}')
            )
        pool.close()
        pool.join()
        log.info(f'  Chain {c1} intrachain HC: {timeit.default_timer()-t0:.2f}s\n')
 
        # ── Interchain ───────────────────────────────────────────────────────
        t0 = timeit.default_timer()
        for c2 in cs[1:]:
            sel_2_hc = selections.get(c2, {}).get('HC', {})
            if not sel_2_hc:
                log.info(f'  Chains {c1}/{c2}: no HC selection for {c2} – skipped.\n')
                continue

            sel_2_hc_str: Dict[str, str] = {}
            n_atoms2 = 0
            keep_running = True
            for label, (sel_str, expected) in sel_2_hc.items():
                if sel_str is None:
                    continue
                if not check_selection(u, sel_str, expected=expected):
                    keep_running = False
                    break
                sel_2_hc_str[label] = sel_str
                n_atoms2 += expected
            if not keep_running:
                break
            log.info(f'Chain {c2}: {len(sel_2_hc_str)} HC groups, {n_atoms2} atoms total')

            fname = f'{c1}_{c2}_HC'
            outnames.append(fname)
            allargs = [
                (top, trj, mf, sel_1_hc_str, cutoff, outpath, fname,
                 ram_thr, sel_2_hc_str, 'hc')
                for mf in frames
            ]
            pool = mp.Pool(processes=nc, initializer=worker_init, initargs=(log_queue,))
            for args in allargs:
                pool.apply_async(
                    multi_run_wrapper, (args,),
                    callback=collect_result,
                    error_callback=lambda e: log.error(f'HC worker crashed: {e}')
                )
            pool.close()
            pool.join()
        log.info(f'  Chain {c1} interchain HC: {timeit.default_timer()-t0:.2f}s\n')
 
    merge_results(outpath, outnames)
    log.info(f'Total HC time: {timeit.default_timer()-tot_start:.2f}s\n')

# ═══════════════════════════════════════════════════════════════════════════════
# POST – merge_results, get_unique_indices
# ═══════════════════════════════════════════════════════════════════════════════

def merge_results(
        outpath: Path, 
        outnames: List[str],
        na_values: Optional[List[str]] = ["", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NULL", "NaN", "n/a", "nan", "null"]
) -> None:
    '''Concatenate per-frame partial CSVs into a single merged CSV per output name.'''
    all_files = [f for f in outpath.iterdir() if f.is_file()]

    log.debug('Merging results in %s', outpath)

    if not outnames:
        log.debug('merge_results: no output names provided')
        return

    for outname in outnames:
        partials = [
            f for f in all_files
            if f.name.startswith(outname) and '_upto_' in f.name and not f.name.endswith('merged.csv')
        ]
        partials.sort(key=lambda x: int(x.stem.split('_')[-1]))

        if not partials:
            log.debug('merge_results: %s - no partial files found', outname)
            continue

        merged = pd.concat(
            [pd.read_csv(f, header=0, index_col=False, na_values=na_values) for f in partials],
            ignore_index=True,
        )

        if merged.shape[0] > 0:
            merged_name = outname + '_merged.csv'
            merged.to_csv(outpath / merged_name, header=True, index=False)
            log.debug('merge_results: %s -> %s', outname, merged_name)
        else:
            log.debug('merge_results: %s - empty result (no interactions found)', outname)

def get_unique_indices(
    outpath: Union[str, Path],
    targets: Optional[List[str]] = None,
) -> None:
    '''Parse merged contact CSVs and logging.info the unique chain-pair identifiers found.'''
    if targets is None:
        targets = ['HBonds', 'Salt_bridges', 'HCs']

    folder_map  = {'Salt_bridges': 'Salt_bridges', 'HBonds': 'H_bonds', 'HCs': 'HCs'}
    role_tokens = {
        'Salt_bridges': {'bas', 'ac'},
        'HBonds':       {'HAcc', 'HDon'},
        'HCs':          {'HC'},
    }
    for target in targets:
        folder = folder_map.get(target)
        if folder is None:
            continue
        p = Path(outpath) / folder
        if not p.exists():
            log.warning(f'get_unique_indices: folder not found: {p}')
            continue

        ignore  = role_tokens.get(target, set()) | {''}
        myfiles = [f for f in p.iterdir() if f.is_file() and f.name.endswith('merged.csv')]

        for fpath in myfiles:
            parts      = fpath.stem.split('_')
            env_chains = [t for t in parts[1:-2] if t not in ignore]
            log.info(f'{target} | {fpath.name} | surrounding chains: {env_chains}')

# ═══════════════════════════════════════════════════════════════════════════════
# POST 2 –  adjacency matrix builder and initial persistency analyses for threshold detection
# ═══════════════════════════════════════════════════════════════════════════════
# ═══════════════════════════════════════════════════════════════════════════════
# helpers
# ═══════════════════════════════════════════════════════════════════════════════

# Folder names under each Adjacency/inttype_folder/ that hold merged CSVs
_SUBFOLDER: Dict[str, str] = {
    'Salt_bridges': 'Salt_bridges',
    'HBonds':       'H_bonds',
    'HCs':          'HCs',
}

# Column names for the two CSV schemas produced by Contacts
_SCHEMA_PAIR  = {'frame', 'id_1', 'id_2', 'distance'}          # SB and HC
_SCHEMA_HBOND = {'Frame', 'Donor_ix', 'Hydrogen_ix', 'Acceptor_ix', 'HA_distance'}  # HB

# contact loaders
# for HC/SB contacts
def load_pair_contacts(
    csv_paths:  List[Path],
    int_df:     pd.DataFrame,
    outpath:    Path,
    inttype:    str,
    chunksize:  int = 100_000,
    na_values: Optional[List[str]] = ["", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NULL", "NaN", "n/a", "nan", "null"]
    ) -> Path:
    '''
    Load SB or HC merged CSVs and return a per-residue-pair frame count table.

    Each CSV has columns: frame, id_1, id_2, distance.
    The id labels encode the residue index as the suffix _r{Resindex}.

    - csv_paths: list of *_merged.csv paths (one per chain-pair computation)
    - int_df:    whole_topology_interactors DataFrame (for node-label look-up)

    Returns a DataFrame with columns:
        res_1, res_2, chain_1, chain_2, label_1, label_2, n_frames_contact, mean_dist
    where res_1 < res_2 (canonical ordering).
    '''
    op = outpath / f'{inttype}_final.csv'

    #print(f'Obtaining merged results for {inttype}, files: ', csv_paths)
    if not csv_paths:
        return op

    # ── per-chunk pre-aggregation → small DataFrames only ────────────────
    chunks = []
    for p in csv_paths:
        for chunk in pd.read_csv(p, header=0, index_col=False, chunksize=chunksize, na_values=na_values):
            if not _SCHEMA_PAIR.issubset(chunk.columns):
                log.warning(f'load_pair_contacts: unexpected columns in {p.name} – skipped')
                break
            chunk = chunk.dropna(subset=['id_1', 'id_2'])
            if chunk.empty:
                continue
            chunk['res_1'] = chunk['id_1'].astype(str).apply(resindex_from_id)
            chunk['res_2'] = chunk['id_2'].astype(str).apply(resindex_from_id)

            # canonical ordering
            swap = chunk['res_1'] > chunk['res_2']
            chunk.loc[swap, ['res_1', 'res_2']] = chunk.loc[swap, ['res_2', 'res_1']].values

            # pre-aggregate this chunk → collapses rows fast
            agg = (chunk.groupby(['res_1', 'res_2'], as_index=False)
                        .agg(frames=('frame', lambda x: set(x.unique())),
                             dsum=('distance', 'sum'),
                             dcnt=('distance', 'size')))
            chunks.append(agg)

    if not chunks:
        return op

    # ── merge per-chunk aggregates ───────────────────────────────────────
    merged = pd.concat(chunks, ignore_index=True)
    final = (merged.groupby(['res_1', 'res_2'], as_index=False)
                   .agg(frames=('frames', lambda s: set().union(*s)),
                        dsum=('dsum', 'sum'),
                        dcnt=('dcnt', 'sum')))
    final['n_frames_contact'] = final['frames'].apply(len)
    final['mean_dist'] = final['dsum'] / final['dcnt']
    # ── visited frames CSV (exploded from set) ───────────────────────────
    visited = (final[['res_1', 'res_2', 'frames']]
               .explode('frames')
               .rename(columns={'frames': 'frame'}))
    visited.to_csv(op.with_name(op.stem + '_visited_frames.csv'),
                   header=True, index=False)

    # ── final aggregated table ───────────────────────────────────────────
    label_map = build_label_map(int_df)
    chain_map = build_chain_map(int_df)
    final['chain_1'] = final['res_1'].map(chain_map)
    final['chain_2'] = final['res_2'].map(chain_map)
    final['label_1'] = final['res_1'].map(label_map)
    final['label_2'] = final['res_2'].map(label_map)

    final[['res_1', 'res_2', 'n_frames_contact', 'mean_dist',
           'chain_1', 'chain_2', 'label_1', 'label_2']].to_csv(
        op, header=True, index=False
    )
    return op

# for Hbonds
def load_hbond_contacts(
    csv_paths: List[Path],
    int_df:    pd.DataFrame,
    outpath: Path,
    chunksize: int = 100_000,
    na_values: Optional[List[str]] = ["", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NULL", "NaN", "n/a", "nan", "null"]
) -> Path:
    '''
    Load H-bond merged CSVs and return a per-residue-pair # frames table.

    Each CSV has columns: Frame, Donor_ix, Hydrogen_ix, Acceptor_ix, Distance.
    Donor_ix and Acceptor_ix are MDAnalysis atom indices that match the
    'atom_index' column of whole_topology_interactors.  The H atom
    (Hydrogen_ix) is excluded from the residue-pair identification.

    Contact frequency is defined as:

        frequency(u, v) = (# distinct frames with ≥1 u–v H-bond) / n_frames

    - csv_paths: list of *_merged.csv paths
    - int_df:    whole_topology_interactors DataFrame

    Returns a DataFrame with the same columns as load_pair_contacts.
    '''
    op = outpath / 'Hbonds_final.csv'

    #print('Obtaining merged results for HBonds, files: ', csv_paths)
    if not csv_paths:
        return op

    # Build atom_index → Resindex map from whole_topology_interactors
    pos_df    = int_df.drop_duplicates('atom_index').sort_values('atom_index').reset_index(drop=True)
    res_by_ix = pos_df['Resindex'].astype(int).to_numpy()

    chunks = []
    for p in csv_paths:
        for chunk in pd.read_csv(p, header=0, index_col=False, chunksize=chunksize, na_values=na_values):
            if not _SCHEMA_HBOND.issubset(chunk.columns):
                log.warning(f'load_hbond_contacts: unexpected columns in {p.name} – skipped')
                break
            chunk = chunk.dropna(subset=['Donor_ix', 'Acceptor_ix'])
            if chunk.empty:
                continue
            chunk = chunk.rename(columns={'Frame': 'frame', 'Distance': 'distance'})
            chunk['res_1'] = res_by_ix[chunk['Donor_ix'].astype(int).to_numpy()]
            chunk['res_2'] = res_by_ix[chunk['Acceptor_ix'].astype(int).to_numpy()]            
            chunk = chunk.dropna(subset=['res_1', 'res_2'])
            if chunk.empty:
                continue
            chunk['res_1'] = chunk['res_1'].astype(int)
            chunk['res_2'] = chunk['res_2'].astype(int)

            # canonical ordering
            swap = chunk['res_1'] > chunk['res_2']
            chunk.loc[swap, ['res_1', 'res_2']] = chunk.loc[swap, ['res_2', 'res_1']].values

            # drop self-contacts
            chunk = chunk[chunk['res_1'] != chunk['res_2']]
            if chunk.empty:
                continue

            agg = (chunk.groupby(['res_1', 'res_2'], as_index=False)
                        .agg(frames=('frame', lambda x: set(x.unique())),
                             dsum=('HA_distance', 'sum'),
                             dcnt=('HA_distance', 'size')))
            chunks.append(agg)

    if not chunks:
        return op

    merged = pd.concat(chunks, ignore_index=True)
    final = (merged.groupby(['res_1', 'res_2'], as_index=False)
                   .agg(frames=('frames', lambda s: set().union(*s)),
                        dsum=('dsum', 'sum'),
                        dcnt=('dcnt', 'sum')))
    final['n_frames_contact'] = final['frames'].apply(len)
    final['mean_dist'] = final['dsum'] / final['dcnt']

    visited = (final[['res_1', 'res_2', 'frames']]
               .explode('frames')
               .rename(columns={'frames': 'frame'}))
    visited.to_csv(op.with_name(op.stem + '_visited_frames.csv'),
                   header=True, index=False)

    label_map = build_label_map(int_df)
    chain_map = build_chain_map(int_df)
    final['chain_1'] = final['res_1'].map(chain_map)
    final['chain_2'] = final['res_2'].map(chain_map)
    final['label_1'] = final['res_1'].map(label_map)
    final['label_2'] = final['res_2'].map(label_map)

    final[['res_1', 'res_2', 'n_frames_contact', 'mean_dist',
           'chain_1', 'chain_2', 'label_1', 'label_2']].to_csv(
        op, header=True, index=False
    )
    return op

def all_subset_combos(lst: List[str]) -> List[Tuple[str, ...]]:
    out = []
    for r in range(1, len(lst) + 1):
        out.extend(combinations(lst, r))
    return out

def threshold_worker(args):
    name, mats, csvdir, plotdir, step, label, xlabel = args
    vals     = np.arange(0.0, 100.0, step)
    n        = next(iter(mats.values())).shape[0]
    maxclust = np.zeros(len(vals), dtype=float)
    deg_distrib: Dict[str, List[Tuple[int, int]]] = {}
    if len(vals) <= 10:
        sample_every = len(vals)
    else:    
        sample_every = max(1, len(vals) // 10)

    for k, v in enumerate(vals):
        merged = np.zeros((n, n), dtype=bool)
        for m in mats.values():
            merged |= (m > v)                     # <- single, uniform rule
        np.fill_diagonal(merged, False)
        G = nx.from_numpy_array(merged.astype(np.int8))
        if G.number_of_edges() == 0:
            maxclust[k] = 0.0
        else:
            biggest = max(len(c) for c in nx.connected_components(G))
            maxclust[k] = biggest * 100.0 / n
        if k % sample_every == 0:
            degs = [d for _, d in G.degree()]
            uniq = sorted(set(degs))
            deg_distrib[f'{v:.1f}'] = [(d, degs.count(d)) for d in uniq]

    pd.DataFrame({'threshold': vals, 'biggest_cluster': maxclust}).to_csv(
        csvdir / f'{name}_{label}_elbow.csv', index=False)
    rows = [{'threshold': thr, 'degree': d, 'count': c}
            for thr, dist in deg_distrib.items() for d, c in dist]
    pd.DataFrame(rows).to_csv(
        csvdir / f'{name}_{label}_degree_distribution.csv', index=False)

    plt.figure()
    plt.scatter(vals, maxclust, s=3)
    plt.xlabel(xlabel); plt.ylabel('Largest cluster size (% of nodes)')
    plt.ylim(0, 100); plt.title(name); plt.tight_layout()
    plt.savefig(plotdir / f'{name}_{label}_elbow.png', dpi=300); plt.close()

    if deg_distrib:
        plt.figure()
        for thr, dist in deg_distrib.items():
            if not dist:
                continue
            xs, ys = zip(*dist)
            plt.plot(xs, ys, marker='o', linewidth=0.6, markersize=2, label=thr)
        plt.xlabel('Degree'); plt.ylabel('# nodes')
        plt.legend(loc='best', fontsize=7, title=xlabel)
        plt.title(f'{name}: degree distribution'); plt.tight_layout()
        plt.savefig(plotdir / f'{name}_{label}_degree_distribution.png', dpi=300)
        plt.close()
    return name

# ═══════════════════════════════════════════════════════════════════════════════
# energy mats methods
# ═══════════════════════════════════════════════════════════════════════════════

# ── energy constants
_NM_PER_ANGSTROM = 0.1
_KE = 138.935458  # kJ*mol^-1*nm*e^-2  (GROMACS Coulomb constant, eps_r = 1)

def load_lj_params(ff_dir: Path, comb_rule: int = 2, na_values=None):
    '''lj_lookup(type_i, type_j) -> combined (sigma_ij, epsilon_ij) in nm, kJ/mol.'''
    at  = pd.read_csv(ff_dir / 'atomtypes_se.csv', na_values=na_values)
    sig = dict(zip(at['atom_type'].astype(str), at['sigma'].astype(float)))
    eps = dict(zip(at['atom_type'].astype(str), at['epsilon'].astype(float)))
    cache: Dict[Tuple[str, str], Tuple[float, float]] = {}
    def lj_lookup(t1: str, t2: str) -> Tuple[float, float]:
        key = (t1, t2) if t1 <= t2 else (t2, t1)
        hit = cache.get(key)
        if hit is not None:
            return hit
        if comb_rule == 2:                       # Lorentz-Berthelot (CHARMM)
            val = ((sig[t1] + sig[t2]) / 2.0, math.sqrt(eps[t1] * eps[t2]))
        else:                                    # comb-rule 1/3
            val = (math.sqrt(sig[t1] * sig[t2]), math.sqrt(eps[t1] * eps[t2]))
        cache[key] = val
        return val
    return lj_lookup

# ─────────────────────────────────────────────────────────────────────────────
# Internal helper: aggregate one buffer to (frame, res_1, res_2) and append.
# ─────────────────────────────────────────────────────────────────────────────
 
def _aggregate_and_append(
    rows_df:      Optional[pd.DataFrame],
    out_path:     Path,
    wrote_header: bool,
) -> bool:
    """Groupby (frame, res_1, res_2) sum across atom-pair rows and append."""
    if rows_df is None or rows_df.empty:
        return wrote_header
    agg = (rows_df.groupby(['frame', 'res_1', 'res_2'], as_index=False)
                  .agg(n_contacts     = ('energy_LJ', 'size'),
                       energy_LJ      = ('energy_LJ', 'sum'),
                       energy_Coulomb = ('energy_Coulomb', 'sum')))
    agg['energy_total'] = agg['energy_LJ'] + agg['energy_Coulomb']
    agg = agg[['frame', 'res_1', 'res_2', 'n_contacts',
               'energy_LJ', 'energy_Coulomb', 'energy_total']]
    mode = 'w' if not wrote_header else 'a'
    agg.to_csv(out_path, header=not wrote_header, index=False, mode=mode)
    return True

def plot_energy_hist(vals: np.ndarray, name: str, pdir: Path) -> None:
    '''Histogram of off-diagonal energies with a fitted normal overlay.'''
    from scipy import stats
    mu, sigma = float(np.mean(vals)), float(np.std(vals))
    skew = float(stats.skew(vals)) if vals.size > 2 else 0.0
    plt.figure()
    plt.hist(vals, bins='auto', density=True, alpha=0.6)
    if sigma > 0:
        xs = np.linspace(float(vals.min()), float(vals.max()), 200)
        plt.plot(xs, stats.norm.pdf(xs, mu, sigma), linewidth=1.2)
    plt.xlabel('Interaction energy (kJ/mol)'); plt.ylabel('density')
    plt.title(f'{name}: energy distribution (μ={mu:.1f}, σ={sigma:.1f}, skew={skew:.2f})')
    plt.tight_layout()
    plt.savefig(pdir / f'{name}_energy_hist.png', dpi=300); plt.close()
 
def get_hc_energy(
    hc_files:  List[Path],
    int_df:    pd.DataFrame,
    label_sel: Dict[str, str],
    top:       Union[str, Path],
    trj:       Union[str, Path],
    n_frames:  int,
    outpath:   Path,
    lj_lookup,
    chunksize: int  = 100_000,
    na_values: Optional[List[str]] = None,
    write_per_frame_csv: bool = True,
) -> Dict[str, Path]:
    '''
    HC energy (kJ/mol, mean per frame), streaming.
 
    For each merged CSV: process one frame at a time. When the frame number
    changes in the chunk stream, the previous frame's HC-unit pairs are
    energy-evaluated (universe.trajectory[frame]) and their (res_1, res_2)
    aggregates flushed to <stem>_energy.csv (append mode).
 
    HC-unit selections are loaded lazily on first use. Bounded memory: one
    frame's HC-unit-pair rows + atom-group cache.
    '''
    from mypsns.utils.matrix_utilities import write_matrix, resindex_from_id
 
    pos = int_df.drop_duplicates('atom_index').sort_values('atom_index').reset_index(drop=True)
    pos_type = pos['atom_type'].astype(str).to_numpy()
    pos_chg  = pos['charge'].astype(float).to_numpy()
 
    all_res = sorted(int_df['Resindex'].unique().astype(int).tolist())
    n   = len(all_res)
    idx = {r: i for i, r in enumerate(all_res)}
 
    if not hc_files:
        return {}
 
    u = mda.Universe(str(top), str(trj))
 
    # lazy atom-group cache
    groups: Dict[str, "mda.AtomGroup"] = {}
    gidx:   Dict[str, np.ndarray]      = {}
 
    def _get_group(lab: str):
        if lab in groups:
            return groups[lab]
        sel = label_sel.get(lab)
        if sel is None:
            return None
        g = u.select_atoms(sel)
        if not len(g):
            return None
        groups[lab] = g
        gidx[lab]   = g.indices
        return g
 
    lj_mat   = np.zeros((n, n))
    coul_mat = np.zeros((n, n))
 
    def _process_frame(frame_int: int, rows: List[Dict]) -> Optional[pd.DataFrame]:
        if not rows:
            return None
        u.trajectory[int(frame_int)]
        dims = u.dimensions
        out = []
        for r in rows:
            l1, l2 = r['id_1'], r['id_2']
            g1, g2 = _get_group(l1), _get_group(l2)
            if g1 is None or g2 is None:
                continue
            ix1, ix2 = gidx[l1], gidx[l2]
            d    = distance_array(g1.positions, g2.positions, box=dims)
            r_nm = np.maximum(d * _NM_PER_ANGSTROM, 1e-6)
            t1, t2 = pos_type[ix1], pos_type[ix2]
            q1, q2 = pos_chg[ix1],  pos_chg[ix2]
            lut = {(a, b): lj_lookup(a, b) for a in set(t1) for b in set(t2)}
            sig = np.array([[lut[(t1[a], t2[b])][0] for b in range(len(t2))] for a in range(len(t1))])
            eps = np.array([[lut[(t1[a], t2[b])][1] for b in range(len(t2))] for a in range(len(t1))])
            sr6  = (sig / r_nm) ** 6
            e_lj = float(np.sum(4.0 * eps * (sr6 * sr6 - sr6)))
            e_co = float(np.sum(_KE * (q1[:, None] * q2[None, :]) / r_nm))
            r1_, r2_ = sorted((int(r['res_1']), int(r['res_2'])))
            i, j = idx[r1_], idx[r2_]
            lj_mat[i, j]   += e_lj; lj_mat[j, i]   += e_lj
            coul_mat[i, j] += e_co; coul_mat[j, i] += e_co
            if write_per_frame_csv:
                out.append({
                    'frame':          frame_int,
                    'res_1':          r1_,
                    'res_2':          r2_,
                    'energy_LJ':      e_lj,
                    'energy_Coulomb': e_co,
                })
        return pd.DataFrame(out) if out else None
 
    written: Dict[str, Path] = {}
 
    for p in hc_files:
        per_frame_path = p.with_name(
            p.stem.replace('_merged', '') + '_energy.csv'
            if '_merged' in p.stem else p.stem + '_energy.csv'
        )
        wrote_header = False
        carry_frame: Optional[int] = None
        carry_rows: List[Dict]     = []
 
        for chunk in pd.read_csv(p, chunksize=chunksize, na_values=na_values):
            chunk = chunk.dropna(subset=['id_1', 'id_2'])
            if chunk.empty:
                continue
            chunk['id_1']  = chunk['id_1'].astype(str)
            chunk['id_2']  = chunk['id_2'].astype(str)
            chunk['res_1'] = chunk['id_1'].apply(resindex_from_id)
            chunk['res_2'] = chunk['id_2'].apply(resindex_from_id)
            chunk = chunk[chunk['res_1'] != chunk['res_2']]
            if chunk.empty:
                continue
 
            for frame, fdf in chunk.groupby('frame', sort=True):
                frame_int = int(frame)
                recs = fdf[['id_1', 'id_2', 'res_1', 'res_2']].to_dict('records')
                if carry_frame is None or carry_frame == frame_int:
                    carry_frame = frame_int
                    carry_rows.extend(recs)
                else:
                    fdf_out = _process_frame(carry_frame, carry_rows)
                    wrote_header = _aggregate_and_append(fdf_out, per_frame_path, wrote_header)
                    carry_frame, carry_rows = frame_int, recs
 
        if carry_frame is not None and carry_rows:
            fdf_out = _process_frame(carry_frame, carry_rows)
            wrote_header = _aggregate_and_append(fdf_out, per_frame_path, wrote_header)
 
        if write_per_frame_csv and wrote_header:
            log.info('  per-frame HC energy: %s', per_frame_path.name)
 
    if n_frames:
        lj_mat   /= n_frames
        coul_mat /= n_frames
    tot = lj_mat + coul_mat
    for tag, m in (('LJ', lj_mat), ('Coulomb', coul_mat), ('total', tot)):
        op = outpath / f'HC_{tag}_energy.dat'
        write_matrix(m, op, fmt='%.6e')
        written[f'HC_{tag}'] = op
        log.info('  HC %s: %s (%dx%d, %d non-zero pairs)',
                 tag, op.name, n, n, np.count_nonzero(m) // 2)
 
    return written

# ═══════════════════════════════════════════════════════════════════════════════
# called methods
# ═══════════════════════════════════════════════════════════════════════════════

def build_results(
        merged: List[Path], 
        interactors: pd.DataFrame,
        frames: Set[int], 
        outp: Path,
        na_values: Optional[List[str]] = ["", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NULL", "NaN", "n/a", "nan", "null"]
):
    log.info(f'Building results, files: {merged}')
    hbonds = [f for f in merged if f.name.endswith('_Hbonds_merged.csv')]
    hcs = [f for f in merged if f.name.endswith('_HC_merged.csv')]
    sbs = [f for f in merged if f.name.endswith('_sb_merged.csv')]
    
    hbond_res = load_hbond_contacts(hbonds, interactors, outp, na_values=na_values)
    log.info(f'Results for Hbonds found at: {hbond_res.resolve()}')

    hcs_res = load_pair_contacts(hcs, interactors, outp, 'HC', na_values=na_values)
    log.info(f'Results for HCs: {hcs_res.resolve()}')

    sbs_res = load_pair_contacts(sbs, interactors, outp, 'SB', na_values=na_values)
    log.info(f'Results for sbs: {sbs_res.resolve()}')      

def build_persistency_matrices(
    final_csvs:   Dict[str, Path],   # {'SB': .../SB_final.csv, ...}
    visited_csvs: Dict[str, Path],   # {'SB': .../SB_final_visited_frames.csv, ...}
    int_df:       pd.DataFrame,
    n_frames:     int,
    outpath:      Path,
    fmt:          str = '%.2f',
    na_values: Optional[List[str]] = ["", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A", "NULL", "NaN", "n/a", "nan", "null"]
) -> Tuple[Dict[str, Path], Dict[str, np.ndarray], Path]:
    """
    Build per-type adjacency matrices (persistency %) + union-combined matrix.
    per-type:  100 * n_frames_contact / n_frames        (from *_final.csv)
    combined:  100 * |union(frames over types)| / n_frames   (from *_visited_frames.csv)
    Returns (paths_dict, matrices_dict).
    """
    all_resindices = sorted(int_df['Resindex'].unique().astype(int).tolist())
    n = len(all_resindices)
    idx_map = {r: i for i, r in enumerate(all_resindices)}

    label_map = build_label_map(int_df)
    labels_path = outpath / 'adjacency_labels.dat'
    with open(labels_path, 'w') as fh:
        for r in all_resindices:
            fh.write(f'{r}\t{label_map.get(r, str(r))}\n')

    written: Dict[str, Path]    = {'labels': labels_path}
    matrices: Dict[str, np.ndarray] = {}

    # per-type
    for inttype, csv_p in final_csvs.items():
        if not csv_p.exists():
            log.warning(f'build_adjacency_matrices: {csv_p} missing → skip')
            continue
        df  = pd.read_csv(csv_p, header=0, index_col=False, na_values=na_values)
        mat = np.zeros((n, n), dtype=float)

        for _, row in df.iterrows():
            r1, r2 = int(row['res_1']), int(row['res_2'])
            if r1 not in idx_map or r2 not in idx_map:
                continue
            i, j = idx_map[r1], idx_map[r2]
            val  = float(row['n_frames_contact']) / n_frames * 100.0
            mat[i, j] = val
            mat[j, i] = val
        matrices[inttype] = mat
        p = outpath / f'{inttype}_adjacency.dat'
        np.savetxt(p, mat, fmt=fmt, delimiter=' ')
        written[inttype] = p
        log.info(f'  {inttype}: {p.name} ({n}×{n}, '
                     f'{np.count_nonzero(mat)//2} non-zero pairs)')

    # combined → union of frames per pair
    pair_frames: Dict[Tuple[int, int], Set[int]] = {}
    
    for inttype, vis_p in visited_csvs.items():
        if not vis_p.exists():
            log.warning(f'build_adjacency_matrices: {vis_p} missing → skip union')
            continue
        df = pd.read_csv(vis_p, header=0, index_col=False, na_values=na_values)
        for _, row in df.iterrows():
            r1, r2 = int(row['res_1']), int(row['res_2'])
            if r1 > r2:
                r1, r2 = r2, r1
            pair_frames.setdefault((r1, r2), set()).add(int(row['frame']))

    combined = np.zeros((n, n), dtype=float)
    for (r1, r2), frames in pair_frames.items():
        if r1 not in idx_map or r2 not in idx_map:
            continue
        i, j = idx_map[r1], idx_map[r2]
        val  = len(frames) / n_frames * 100.0
        combined[i, j] = val
        combined[j, i] = val
    matrices['combined'] = combined
    p = outpath / 'combined_adjacency.dat'
    np.savetxt(p, combined, fmt=fmt, delimiter=' ')
    written['combined'] = p
    log.info(f'  combined (union): {p.name} '
                 f'({n}×{n}, {np.count_nonzero(combined)//2} non-zero pairs)')
    return written, matrices, labels_path

def get_interaction_energy(
    sb_files:  List[Path],
    hb_files:  List[Path],
    int_df:    pd.DataFrame,
    n_frames:  int,
    outpath:   Path,
    lj_lookup,
    chunksize: int  = 100_000,
    na_values: Optional[List[str]] = None,
    write_per_frame_csv: bool = True,
) -> Dict[str, Path]:
    '''
    SB + H-bond LJ + Coulomb energy (kJ/mol). Bounded memory.
 
    Matrix outputs (unchanged):
      {SB,Hbonds}_LJ_energy.dat
      {SB,Hbonds}_Coulomb_energy.dat
      {SB,Hbonds}_total_energy.dat
 
    Per-merged-file companion (when write_per_frame_csv=True):
      <stem>_energy.csv with columns
        frame, res_1, res_2, n_contacts, energy_LJ, energy_Coulomb, energy_total
    '''
    from mypsns.utils.matrix_utilities import write_matrix
 
    pos = int_df.drop_duplicates('atom_index').sort_values('atom_index').reset_index(drop=True)
    pos_type = pos['atom_type'].astype(str).to_numpy()
    pos_res  = pos['Resindex'].astype(int).to_numpy()
    pos_chg  = pos['charge'].astype(float).to_numpy()
 
    all_res = sorted(int_df['Resindex'].unique().astype(int).tolist())
    n   = len(all_res)
    idx = {r: i for i, r in enumerate(all_res)}
 
    written: Dict[str, Path] = {}
 
    for kind, files in (('SB', sb_files), ('Hbonds', hb_files)):
        if not files:
            continue
        lj_mat   = np.zeros((n, n))
        coul_mat = np.zeros((n, n))
 
        for p in files:
            per_frame_path = p.with_name(
                p.stem.replace('_merged', '') + '_energy.csv'
                if '_merged' in p.stem else p.stem + '_energy.csv'
            )
            wrote_header = False
            carry: Optional[pd.DataFrame] = None     # rows whose frame may continue
 
            for chunk in pd.read_csv(p, chunksize=chunksize, na_values=na_values):
                # ── decode chunk by kind ─────────────────────────────────────
                if kind == 'SB':
                    chunk = chunk.dropna(subset=['id_1', 'id_2'])
                    if chunk.empty:
                        continue
                    ix1  = chunk['id_1'].str.split('_r').str[0].astype(int).to_numpy()
                    ix2  = chunk['id_2'].str.split('_r').str[0].astype(int).to_numpy()
                    dist = chunk['distance'].to_numpy(dtype=float)
                    frames = chunk['frame'].astype(int).to_numpy()
                else:  # H-bonds — H···Acceptor at HA_distance
                    chunk = chunk.dropna(subset=['Hydrogen_ix', 'Acceptor_ix'])
                    if chunk.empty:
                        continue
                    ix1  = chunk['Hydrogen_ix'].astype(int).to_numpy()
                    ix2  = chunk['Acceptor_ix'].astype(int).to_numpy()
                    dist = chunk['HA_distance'].to_numpy(dtype=float)
                    fc = 'Frame' if 'Frame' in chunk.columns else 'frame'
                    frames = chunk[fc].astype(int).to_numpy()
 
                # ── exclude within-residue ───────────────────────────────────
                res1 = pos_res[ix1]
                res2 = pos_res[ix2]
                keep = res1 != res2
                if not keep.any():
                    continue
                ix1, ix2, dist = ix1[keep], ix2[keep], dist[keep]
                res1, res2 = res1[keep], res2[keep]
                frames = frames[keep]
 
                # ── LJ params per pair (dedup → tiny lookup table) ───────────
                df = pd.DataFrame({'t1': pos_type[ix1], 't2': pos_type[ix2]})
                up = df.drop_duplicates().copy()
                up[['sig', 'eps']] = up.apply(
                    lambda x: pd.Series(lj_lookup(x['t1'], x['t2'])), axis=1)
                df = df.merge(up, on=['t1', 't2'], how='left')
 
                # ── energies ─────────────────────────────────────────────────
                r_nm = dist * _NM_PER_ANGSTROM
                sr6  = (df['sig'].to_numpy() / r_nm) ** 6
                e_lj = 4.0 * df['eps'].to_numpy() * (sr6 * sr6 - sr6)
                e_co = _KE * pos_chg[ix1] * pos_chg[ix2] / r_nm
                lo   = np.minimum(res1, res2)
                hi   = np.maximum(res1, res2)
 
                # ── accumulate matrices (sum over all atom pairs of all frames)
                mat_agg = (pd.DataFrame({'lo': lo, 'hi': hi, 'e_lj': e_lj, 'e_co': e_co})
                             .groupby(['lo', 'hi'], as_index=False).sum())
                for _, row in mat_agg.iterrows():
                    r1, r2 = int(row['lo']), int(row['hi'])
                    if r1 in idx and r2 in idx:
                        i, j = idx[r1], idx[r2]
                        lj_mat[i, j]   += row['e_lj']
                        lj_mat[j, i]   += row['e_lj']
                        coul_mat[i, j] += row['e_co']
                        coul_mat[j, i] += row['e_co']
 
                # ── stream per-frame CSV ─────────────────────────────────────
                if not write_per_frame_csv:
                    continue
 
                chunk_df = pd.DataFrame({
                    'frame':          frames,
                    'res_1':          lo,
                    'res_2':          hi,
                    'energy_LJ':      e_lj,
                    'energy_Coulomb': e_co,
                })
                if carry is not None:
                    chunk_df = pd.concat([carry, chunk_df], ignore_index=True)
                    carry = None
 
                # Hold back the trailing frame — it may continue in next chunk
                last_frame = int(chunk_df['frame'].iloc[-1])
                last_mask  = chunk_df['frame'].values == last_frame
                carry      = chunk_df.loc[last_mask].copy()
                to_write   = chunk_df.loc[~last_mask]
 
                wrote_header = _aggregate_and_append(to_write, per_frame_path, wrote_header)
 
            # Flush trailing carry
            if write_per_frame_csv and carry is not None:
                wrote_header = _aggregate_and_append(carry, per_frame_path, wrote_header)
 
            if write_per_frame_csv and wrote_header:
                log.info('  per-frame energy: %s', per_frame_path.name)
 
        # Write matrices
        if n_frames:
            lj_mat   /= n_frames
            coul_mat /= n_frames
        tot = lj_mat + coul_mat
        for tag, m in (('LJ', lj_mat), ('Coulomb', coul_mat), ('total', tot)):
            op = outpath / f'{kind}_{tag}_energy.dat'
            write_matrix(m, op, fmt='%.6e')
            written[f'{kind}_{tag}'] = op
            log.info('  %s %s: %s (%dx%d, %d non-zero pairs)',
                     kind, tag, op.name, n, n, np.count_nonzero(m) // 2)
 
    return written
 
def energy_to_percentage(
    mats:    Dict[str, np.ndarray],
    outdir:  Path,
    use_abs: bool = False,
    plot:    bool = True,
) -> Dict[str, np.ndarray]:
    '''
    Map raw energy matrices onto a 0..100 "interaction strength %" so the
    interaction-agnostic threshold_analysis can scan them with the same
    'edge if value > threshold' rule as persistency.

        strength(E) = |E|              if use_abs   (keeps strong repulsive too)
                    = max(0, -E)       otherwise    (attractive only)
        pct         = 100 * strength / strength.max()   over nonzero off-diagonal
        zeros (no contact) stay 0.

    With plot=True, each type's off-diagonal energy distribution is written to
    outdir/energy_distribution/<type>_energy_hist.png (normal fit + μ/σ/skew) so
    you can judge whether |E| (≈Gaussian, symmetric) or the signed mapping fits.
    '''
    if plot:
        pdir = outdir / 'energy_distribution'
        pdir.mkdir(parents=True, exist_ok=True)

    pct_mats: Dict[str, np.ndarray] = {}
    for name, m in mats.items():
        nz = m[np.triu_indices_from(m, k=1)]
        nz = nz[nz != 0.0]
        if plot and nz.size:
            plot_energy_hist(nz, name, pdir)

        strength = np.abs(m) if use_abs else np.clip(-m, 0.0, None)
        np.fill_diagonal(strength, 0.0)
        smax = strength[np.triu_indices_from(strength, k=1)].max() if nz.size else 0.0
        pct_mats[name] = (100.0 * strength / smax) if smax > 0 else np.zeros_like(strength)
    return pct_mats

def threshold_analysis(
    per_type_mats: Dict[str, np.ndarray],
    outdir:        Path,
    step:          float = 0.1,
    nc:            int   = 1,
    label:         str   = 'persistency',
    xlabel:        Optional[str] = None,
) -> None:
    
    if xlabel is None:
        xlabel = r'$p_{min}$ (%)' if label == 'persistency' else r'strength (%)'
    csvdir  = outdir / 'csvs';  csvdir.mkdir(parents=True, exist_ok=True)
    plotdir = outdir / 'plots'; plotdir.mkdir(parents=True, exist_ok=True)

    types = [t for t in per_type_mats if t != 'combined']
    if not types:
        log.warning('threshold_analysis: no per-type matrices → abort')
        return

    args = [
        ('_'.join(combo), {k: per_type_mats[k] for k in combo},
         csvdir, plotdir, step, label, xlabel)
        for combo in all_subset_combos(types)
    ]
    start = timeit.default_timer()
    if nc > 1:
        with mp.Pool(processes=nc) as pool:
            pool.map(threshold_worker, args)
    else:
        for a in args:
            threshold_worker(a)
    log.info('threshold_analysis (%s): done in %.2f s on %d proc(s)',
             label, timeit.default_timer() - start, nc)
 
 
# ═════════════════════════════════════════════════════════════════════════════
# get_interaction_energy — streaming
# ═════════════════════════════════════════════════════════════════════════════
 
 
