"""
matrix_utilities.py — shared matrix I/O and residue-pair helpers.

Used by:
- mypsns/utils/contacts_utilities.py    (writes matrices + per-frame energy CSVs)
- mypsns/utils/network_utilities.py     (reads matrices for the static network)
- mypsns/services/network.py            (per-frame energy CSV → time-resolved networks)

Layout of an `Adjacency/<inttype>/` folder:
  adjacency_labels.dat              ← write_labels / load_labels
  <kind>_adjacency.dat              ← persistency matrix (kind: SB | HC | Hbonds)
  combined_adjacency.dat            ← persistency, union across types
  <kind>_{LJ,Coulomb,total}_energy.dat
  <kind>_final.csv                  ← per-residue-pair static summary
  <kind>_final_visited_frames.csv   ← per (frame, res_1, res_2): which frames had contact
  <type-folder>/                    ← e.g. HCs/  – per-chain-pair merged CSVs
       <pair>_<kind>_merged.csv      (atom/cluster pair × frame, distance)
       <pair>_<kind>_energy.csv      (frame × res_1 × res_2, energies)   ← NEW
"""
from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

# ═════════════════════════════════════════════════════════════════════════════
# Residue label / chain / id helpers
# ═════════════════════════════════════════════════════════════════════════════

def build_label_map(int_df: pd.DataFrame) -> Dict[int, str]:
    """{Resindex: 'Chain:resname:Resindex'} from whole_topology_interactors."""
    unique_res = int_df.drop_duplicates(subset='Resindex')
    return {
        int(row['Resindex']): f"{row['Chain']}:{row['resname']}:{int(row['Resindex'])}"
        for _, row in unique_res.iterrows()
    }

def build_chain_map(int_df: pd.DataFrame) -> Dict[int, str]:
    """{Resindex: chain_id} from whole_topology_interactors."""
    unique_res = int_df.drop_duplicates(subset='Resindex')
    return dict(zip(
        unique_res['Resindex'].astype(int),
        unique_res['Chain'],
    ))

_RESINDEX_RE = re.compile(r'_r(\d+)(?:_c.+)?$')

def resindex_from_id(atom_id_str: str) -> int:
    """
    Extract Resindex from atom-id strings written by the contact evaluators.
    Supported forms:
      '{ix}_r{Resindex}'              (SB, HB atom_distance_evaluator)
      '{ids}_r{Resindex}'             (HC residue-mode)
      '{ids}_r{Resindex}_c{cluster}'  (HC cluster-mode)
      'HC{n}_{Resindex}' / 'HC_{Resindex}'  (simplified HC dict-keys)
    """
    m = _RESINDEX_RE.search(atom_id_str)
    if m is not None:
        return int(m.group(1))
    m = re.search(r'_(\d+)$', atom_id_str)
    if m is not None:
        return int(m.group(1))
    raise ValueError(f'Cannot parse Resindex from atom label: {atom_id_str!r}')

# ═════════════════════════════════════════════════════════════════════════════
# Matrix I/O
# ═════════════════════════════════════════════════════════════════════════════

def write_labels(int_df: pd.DataFrame, path: Path) -> List[int]:
    """
    Write 'adjacency_labels.dat' (one '<resindex>\\t<label>' per line) at `path`
    and return the canonical ordered Resindex list used for rows/cols.
    """
    path = Path(path)
    label_map = build_label_map(int_df)
    all_resindices = sorted(int_df['Resindex'].unique().astype(int).tolist())
    with open(path, 'w') as f:
        for r in all_resindices:
            f.write(f'{r}\t{label_map.get(r, str(r))}\n')
    return all_resindices

def load_labels(path: Path) -> List[int]:
    """Read 'adjacency_labels.dat' → ordered list of Resindex (int)."""
    out: List[int] = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split('\t', 1)
        try:
            out.append(int(parts[0]))
        except ValueError:
            continue
    return out

def write_matrix(mat: np.ndarray, path: Path, fmt: str = '%.6e') -> Path:
    """Square matrix → whitespace-separated .dat (np.savetxt)."""
    np.savetxt(path, mat, fmt=fmt, delimiter=' ')
    return Path(path)

def load_matrix(path: Path) -> np.ndarray:
    """Load whitespace-separated square matrix; pure np.loadtxt."""
    return np.loadtxt(path)

# ═════════════════════════════════════════════════════════════════════════════
# Pairs ↔ matrix pivots
# ═════════════════════════════════════════════════════════════════════════════

def pairs_to_matrix(
    pairs_df:  pd.DataFrame,
    labels:    List[int],
    value_col: str   = 'value',
    scale:     float = 1.0,
    symmetric: bool  = True,
) -> np.ndarray:
    """
    Pivot (res_1, res_2, value) → NxN matrix on `labels` order.

    Pairs whose endpoints aren't in `labels` are silently dropped.
    """
    n = len(labels)
    idx = {r: i for i, r in enumerate(labels)}
    mat = np.zeros((n, n), dtype=float)
    for _, row in pairs_df.iterrows():
        r1, r2 = int(row['res_1']), int(row['res_2'])
        if r1 not in idx or r2 not in idx:
            continue
        i, j = idx[r1], idx[r2]
        v = float(row[value_col]) * scale
        mat[i, j] = v
        if symmetric and i != j:
            mat[j, i] = v
    return mat

def matrix_to_pairs(
    mat:        np.ndarray,
    labels:     List[int],
    value_name: str   = 'value',
    min_value:  float = 0.0,
    upper_only: bool  = True,
) -> pd.DataFrame:
    """
    Inverse of pairs_to_matrix. By default returns only upper-triangle entries
    strictly above `min_value`.
    """
    rows = []
    n = mat.shape[0]
    for i in range(n):
        for j in range((i + 1) if upper_only else 0, n):
            if i == j:
                continue
            v = float(mat[i, j])
            if v > min_value:
                rows.append({'res_1': labels[i], 'res_2': labels[j], value_name: v})
    return pd.DataFrame(rows)

# ═════════════════════════════════════════════════════════════════════════════
# Time-resolved CSV ↔ matrix
# ═════════════════════════════════════════════════════════════════════════════

def per_frame_to_matrix(
    per_frame_df: pd.DataFrame,
    labels:       List[int],
    value_col:    str   = 'energy_total',
    n_frames:     int   = 1,
    agg:          str   = 'mean',
) -> np.ndarray:
    """
    Collapse a per-(frame, res_1, res_2) table to a single residue×residue
    matrix.

    agg='mean'  → sum over frames, divide by n_frames (mean per frame; matches
                  what get_interaction_energy stores in *_total_energy.dat).
    agg='sum'   → sum over frames only.
    """
    if per_frame_df.empty:
        n = len(labels)
        return np.zeros((n, n), dtype=float)
    pair_sum = (per_frame_df
                .groupby(['res_1', 'res_2'], as_index=False)[value_col]
                .sum())
    if agg == 'mean':
        pair_sum[value_col] = pair_sum[value_col] / max(n_frames, 1)
    elif agg != 'sum':
        raise ValueError(f"per_frame_to_matrix: unknown agg={agg!r}")
    return pairs_to_matrix(pair_sum, labels, value_col=value_col)

# ═════════════════════════════════════════════════════════════════════════════
# Contacts → Network matrix loaders
# ═════════════════════════════════════════════════════════════════════════════
 
# Maps the public target name (used in NetworkConfig.targets) to the filename
# prefix Contacts uses for matrices and merged CSVs.
_MATRIX_PREFIX: Dict[str, str] = {
    'Salt_bridges': 'SB',
    'HBonds':       'Hbonds',
    'HCs':          'HC',
}
_MATRIX_SUFFIX: Dict[str, str] = {
    'persistency': '_adjacency.dat',
    'energy':      '_total_energy.dat',
}
_PER_FRAME_SUBFOLDER: Dict[str, str] = {
    'Salt_bridges': 'Salt_bridges',
    'HBonds':       'H_bonds',
    'HCs':          'HCs',
}
_PER_FRAME_SUFFIX: Dict[str, str] = {
    'Salt_bridges': '_sb_energy.csv',
    'HBonds':       '_Hbonds_energy.csv',
    'HCs':          '_HC_energy.csv',
}
 
 
def energy_transform(
    mat:     np.ndarray,
    use_abs: bool = False,
) -> np.ndarray:
    """
    Convert a signed energy matrix to a non-negative weight matrix.
 
    use_abs=False (default): only attractive contributions → weight = max(0, -E)
    use_abs=True            : both directions count       → weight = |E|
 
    Matches the convention in get_persistence_analysis / energy_to_percentage.
    """
    return np.abs(mat) if use_abs else np.maximum(0.0, -mat)


def strength_to_percentage(mat: np.ndarray) -> np.ndarray:
    """
    Map a non-negative strength matrix onto 0..100 "% of the strongest edge".

        pct[i,j] = 100 * mat[i,j] / max(off-diagonal mat)

    The diagonal is zeroed. If every off-diagonal entry is 0, returns zeros.
    Matches the normalisation in contacts_utilities.energy_to_percentage, so a
    percentage threshold here means the same thing the energy elbow/degree
    plots are scanning.
    """
    m = np.array(mat, dtype=float, copy=True)
    np.fill_diagonal(m, 0.0)
    iu = np.triu_indices_from(m, k=1)
    smax = m[iu].max() if m[iu].size else 0.0
    if smax <= 0:
        return np.zeros_like(m)
    return 100.0 * m / smax
 
 
def load_per_type_matrix(
    folder: Path,
    target: str,
    source: str,
) -> Optional[np.ndarray]:
    """
    Read one per-type matrix from a Contacts output folder.
 
    `target` is the public name (Salt_bridges / HBonds / HCs).
    `source` is 'persistency' or 'energy'. Returns None if missing.
    """
    folder = Path(folder)
    prefix = _MATRIX_PREFIX.get(target)
    suffix = _MATRIX_SUFFIX.get(source)
    if prefix is None or suffix is None:
        raise ValueError(f'unknown target={target!r} or source={source!r}')
    p = folder / f'{prefix}{suffix}'
    if not p.exists():
        log.warning('load_per_type_matrix: %s not found', p)
        return None
    return load_matrix(p) 
 
def load_combined_per_type(
    folder:        Path,
    targets:       List[str],
    source:        str,
    use_abs:       bool = False,
) -> Tuple[Optional[np.ndarray], List[int]]:
    """
    Sum per-type matrices across the requested `targets` for one `source`.
 
    For source='energy', each per-type matrix is first transformed by
    `energy_transform(..., use_abs)` so the result is a non-negative weight.
    For source='persistency', matrices are added directly (entries are %).
 
    Returns (matrix or None, labels). Missing per-type files are skipped with
    a warning; a None matrix means no requested target had a file present.
    """
    labels_path = Path(folder) / 'adjacency_labels.dat'
    if not labels_path.exists():
        log.warning('load_combined_per_type: %s missing', labels_path)
        return None, []
    labels = load_labels(labels_path)
 
    summed: Optional[np.ndarray] = None
    for target in targets:
        m = load_per_type_matrix(folder, target, source)
        if m is None:
            continue
        if source == 'energy':
            m = energy_transform(m, use_abs=use_abs)
        summed = m.copy() if summed is None else summed + m
    return summed, labels

def load_per_frame_energy(
    folder:   Path,
    targets:  List[str],
    use_abs:  bool = False,
) -> pd.DataFrame:
    """
    Concatenate the per-frame energy CSVs for the requested `targets`.
 
    Each requested target's CSVs live in `folder/<subfolder>/*<suffix>` where
    subfolder/suffix come from _PER_FRAME_{SUBFOLDER,SUFFIX}.
 
    Per-row weight is derived from `energy_total` via :func:`energy_transform`.
    Returns columns: frame, res_1, res_2, n_contacts, energy_total, weight.
    Empty DataFrame if no matching files.
    """
    folder = Path(folder)
    parts: List[pd.DataFrame] = []
    for target in targets:
        sub    = _PER_FRAME_SUBFOLDER.get(target)
        suffix = _PER_FRAME_SUFFIX.get(target)
        if sub is None or suffix is None:
            log.warning('load_per_frame_energy: unknown target %s', target)
            continue
        target_dir = folder / sub
        if not target_dir.exists():
            log.warning('load_per_frame_energy: %s missing', target_dir)
            continue
        files = sorted(target_dir.glob(f'*{suffix}'))
        if not files:
            log.warning('load_per_frame_energy: no %s files in %s', suffix, target_dir)
            continue
        for f in files:
            df = pd.read_csv(f, usecols=['frame', 'res_1', 'res_2',
                                         'n_contacts', 'energy_total'])
            df['target'] = target
            parts.append(df)
            log.debug('  %s: %d rows', f.name, len(df))
    if not parts:
        return pd.DataFrame(columns=['frame', 'res_1', 'res_2',
                                     'n_contacts', 'energy_total', 'weight'])
    out = pd.concat(parts, ignore_index=True)
    out['res_1'] = out['res_1'].astype(int)
    out['res_2'] = out['res_2'].astype(int)
    out['frame'] = out['frame'].astype(int)
    # Apply attractive-only / abs transform to per-row energy
    e = out['energy_total'].to_numpy(dtype=float)
    out['weight'] = np.abs(e) if use_abs else np.maximum(0.0, -e)
    return out

def window_per_frame(
    per_frame_df: pd.DataFrame,
    window_size:  int,
    value_cols:   List[str],
    stride:       Optional[int] = None,
) -> pd.DataFrame:
    """
    Time-windowed aggregation of a per-(frame, res_1, res_2) table.

    For every window [t, t+window_size) (stride defaults to window_size, i.e.
    non-overlapping), sums all listed `value_cols` per (res_1, res_2).

    Returns a DataFrame indexed by (window_start, res_1, res_2) with summed
    columns, suitable for one graph per window.
    """
    if per_frame_df.empty:
        return pd.DataFrame(columns=['window_start', 'res_1', 'res_2', *value_cols])
    stride = stride or window_size
    df = per_frame_df.copy()
    df['window_start'] = (df['frame'].astype(int) // stride) * stride
    agg = (df.groupby(['window_start', 'res_1', 'res_2'], as_index=False)
             [value_cols].sum())
    return agg