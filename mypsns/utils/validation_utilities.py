"""
validation_utilities.py — compare myPSNs residue matrices with reference graphs
(e.g. PyInteraph 1 ``hb-graph.dat`` / ``hc-graph.dat`` / ``sb-graph.dat``).

Replaces the ad-hoc ``mypsns/prove.py``. Used by:
- mypsns/commands/compare_cmd.py          (`mypsns compare`)
- tests/test_pyinteraph_regression.py     (SB/HC ≤ 0.2 %, H-bonds ≤ 5 %)

Conventions — all inherited from prove.py, so numbers and CSVs produced before
the rewrite stay comparable (the CSVs are byte-identical):

* Only the strict lower triangle (row > col) is compared. The diagonal is
  ignored: PyInteraph's ``hb-graph.dat`` carries intra-residue H-bonds there
  (24 on the test case), which myPSNs never writes, and PyInteraph's own
  ``filter_graph`` zeroes it as well.
* An *edge* is an off-diagonal pair with value > 0. "shared" edges are > 0 in
  both matrices.
* "% weight difference" = 100 · Σ|computed − truth| / Σ truth over that
  triangle. It is the number quoted in CLAUDE.md / README.
* Two entries *differ* when ``not np.isclose(truth, computed, atol=atol)``
  (default atol 1e-2; numpy's default rtol). Differing entries are split into
    - selection mismatches (``possible_sele_mismatch.csv``): exactly one of the
      two is 0, i.e. the pair is an edge in one matrix only. These point at a
      role/selection difference (which atoms may interact at all);
    - parameter mismatches (``possible_params_mismatch.csv``): both non-zero
      but different. These point at cutoffs, geometry or normalisation.
  ``row``/``col`` are 0-based matrix indices (the ``adjacency_labels.dat``
  order); pass a labels file to also get residue labels.

Matrices are matched by a free-form *type* key. The built-in names
(``HBonds``/``HC``/``SB``) map to PyInteraph and myPSNs file names; any other
key works when the paths are given explicitly (e.g. a macro-IIN).
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

# ═════════════════════════════════════════════════════════════════════════════
# File-name conventions and legacy datasets
# ═════════════════════════════════════════════════════════════════════════════

# Type keys in the order prove.py reported them (and wrote the CSV rows).
DEFAULT_TYPES: Tuple[str, ...] = ('HBonds', 'HC', 'SB')

# PyInteraph 1 output names: <stem><suffix>.dat (suffix e.g. '_1micros').
PYINTERAPH_GRAPH_STEMS: Dict[str, str] = {
    'HBonds': 'hb-graph',
    'HC':     'hc-graph',
    'SB':     'sb-graph',
}
# myPSNs contacts output names inside Adjacency/<inttype>/.
MYPSNS_ADJACENCY_NAMES: Dict[str, str] = {
    'HBonds': 'Hbonds_adjacency.dat',
    'HC':     'HC_adjacency.dat',
    'SB':     'SB_adjacency.dat',
}

# The two datasets prove.py compared, relative to a data root. The root is
# $MYPSNS_VALIDATION_ROOT if set, else the lab folder they were produced in.
VALIDATION_ROOT_ENV = 'MYPSNS_VALIDATION_ROOT'
DEFAULT_VALIDATION_ROOT = Path('/home/davide/storage/lab_dati')
PRESETS: Dict[str, Dict[str, str]] = {
    # CHARMM test case shipped with PyInteraph 1 (149 residues). This is the
    # validation quoted in CLAUDE.md and asserted by the regression test.
    'pyinteraph_testcase': {
        'truth_dir':    'risultati_paragone/pyinter_test',
        'truth_suffix': '',
        'computed_dir': 'psns_out/cli_outs/pyinteraph_testcase/Adjacency/named_interactors',
    },
    # Recoverin–GRK1 (219 residues), PyInteraph run on 1 µs.
    'rec_grk1': {
        'truth_dir':    'risultati_paragone/rec',
        'truth_suffix': '_1micros',
        'computed_dir': 'psns_out/cli_outs/rec_grk1/Adjacency/named_interactors',
    },
}

# Column layout of the mismatch CSVs (unchanged from prove.py).
SELE_COLUMNS   = ['row', 'col', 'intertype', 'present_in', 'non-0 val']
PARAMS_COLUMNS = ['row', 'col', 'intertype', 'difference']


def validation_root(data_root: Optional[Path] = None) -> Path:
    """Data root for :data:`PRESETS`: explicit > $MYPSNS_VALIDATION_ROOT > default."""
    if data_root is not None:
        return Path(data_root)
    env = os.environ.get(VALIDATION_ROOT_ENV)
    return Path(env) if env else DEFAULT_VALIDATION_ROOT


def reference_paths(
    truth_dir: Path,
    suffix:    str = '',
    types:     Sequence[str] = DEFAULT_TYPES,
) -> Dict[str, Path]:
    """{type: <truth_dir>/<pyinteraph stem><suffix>.dat} for the built-in types."""
    return {t: Path(truth_dir) / f'{PYINTERAPH_GRAPH_STEMS[t]}{suffix}.dat'
            for t in types if t in PYINTERAPH_GRAPH_STEMS}


def adjacency_paths(
    computed_dir: Path,
    types:        Sequence[str] = DEFAULT_TYPES,
) -> Dict[str, Path]:
    """{type: <computed_dir>/<myPSNs adjacency name>} for the built-in types."""
    return {t: Path(computed_dir) / MYPSNS_ADJACENCY_NAMES[t]
            for t in types if t in MYPSNS_ADJACENCY_NAMES}


def preset_paths(
    name:      str,
    data_root: Optional[Path] = None,
) -> Tuple[Dict[str, Path], Dict[str, Path]]:
    """(truth paths, computed paths) for a legacy dataset in :data:`PRESETS`."""
    if name not in PRESETS:
        raise ValueError(f'unknown preset {name!r}; choose from {sorted(PRESETS)}')
    root = validation_root(data_root)
    p = PRESETS[name]
    return (reference_paths(root / p['truth_dir'], p['truth_suffix']),
            adjacency_paths(root / p['computed_dir']))


def read_label_names(path: Path) -> List[str]:
    """``adjacency_labels.dat`` → residue labels ('Chain:resname:Resindex') in
    matrix order. Lines without a label fall back to the Resindex."""
    out: List[str] = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split('\t', 1)
        out.append(parts[1].strip() if len(parts) > 1 else parts[0].strip())
    return out

# ═════════════════════════════════════════════════════════════════════════════
# Comparison
# ═════════════════════════════════════════════════════════════════════════════

@dataclass
class MatrixComparison:
    """Result of comparing one reference matrix (truth) with one computed matrix."""

    inttype:        str
    truth_path:     Optional[Path] = None
    computed_path:  Optional[Path] = None
    shape:          Optional[Tuple[int, ...]] = None
    edges_truth:    int = 0
    edges_computed: int = 0
    shared:         int = 0
    only_truth:     int = 0
    only_computed:  int = 0
    value_mismatch: int = 0        # shared edges whose values differ (> atol)
    truth_total:    float = 0.0
    computed_total: float = 0.0
    abs_diff:       float = 0.0
    sele_rows:      List[dict] = field(default_factory=list)
    params_rows:    List[dict] = field(default_factory=list)
    error:          Optional[str] = None

    @property
    def pct_diff(self) -> float:
        """100 · Σ|computed − truth| / Σ truth (nan when the truth is empty)."""
        if self.truth_total == 0:
            return float('nan')
        return 100.0 * self.abs_diff / self.truth_total

    def summary_row(self) -> dict:
        """One row of ``comparison_summary.csv``."""
        return {
            'intertype':      self.inttype,
            'edges_truth':    self.edges_truth,
            'edges_computed': self.edges_computed,
            'shared':         self.shared,
            'only_truth':     self.only_truth,
            'only_computed':  self.only_computed,
            'value_mismatch': self.value_mismatch,
            'truth_total':    self.truth_total,
            'computed_total': self.computed_total,
            'abs_diff':       self.abs_diff,
            'pct_diff':       self.pct_diff,
            'error':          self.error or '',
            'truth_path':     str(self.truth_path or ''),
            'computed_path':  str(self.computed_path or ''),
        }


def compare_matrices(
    truth:    np.ndarray,
    computed: np.ndarray,
    inttype:  str,
    atol:     float = 1e-2,
) -> MatrixComparison:
    """
    Compare two square residue matrices on the strict lower triangle (see the
    module docstring for every definition). A shape mismatch is reported in
    ``error`` with all counts left at 0, as prove.py skipped such a type.
    """
    res = MatrixComparison(inttype=inttype)
    truth = np.asarray(truth, dtype=float)
    computed = np.asarray(computed, dtype=float)
    if truth.shape != computed.shape:
        res.error = f'shape mismatch: truth {truth.shape} vs computed {computed.shape}'
        return res
    if truth.ndim != 2 or truth.shape[0] != truth.shape[1]:
        res.error = f'not a square matrix: {truth.shape}'
        return res
    res.shape = truth.shape

    lower = np.tril(np.ones(truth.shape, dtype=bool), k=-1)
    t, k = truth[lower], computed[lower]
    t_edge, k_edge = t > 0, k > 0
    res.edges_truth    = int(t_edge.sum())
    res.edges_computed = int(k_edge.sum())
    res.shared         = int((t_edge & k_edge).sum())
    res.only_truth     = res.edges_truth - res.shared
    res.only_computed  = res.edges_computed - res.shared
    res.truth_total    = float(t.sum())
    res.computed_total = float(k.sum())
    res.abs_diff       = float(np.abs(k - t).sum())

    # Row-major scan of the lower triangle, so rows come out in the same order
    # prove.py wrote them. Same arithmetic too (|k - t| and t - k on numpy
    # scalars), so the CSVs are byte-identical to the old ones.
    differ = lower & ~np.isclose(truth, computed, atol=atol)
    for r, c in zip(*np.where(differ)):
        tv, kv = truth[r, c], computed[r, c]
        t_zero, k_zero = np.equal(tv, 0.0), np.equal(kv, 0.0)
        if t_zero ^ k_zero:                    # an edge in one matrix only
            res.sele_rows.append({
                'row': r, 'col': c, 'intertype': inttype,
                'present_in': 'computed' if t_zero else 'true',
                'non-0 val': abs(kv - tv),
            })
        elif not t_zero and not k_zero:        # both edges, values differ
            res.params_rows.append({
                'row': r, 'col': c, 'intertype': inttype,
                'difference': tv - kv,
            })
    res.value_mismatch = len(res.params_rows)
    return res


def compare_files(
    pairs: Dict[str, Tuple[Path, Path]],
    atol:  float = 1e-2,
) -> List[MatrixComparison]:
    """Load and compare every ``{type: (truth_path, computed_path)}`` pair, in
    insertion order. Unreadable files are reported in ``error``."""
    results: List[MatrixComparison] = []
    for inttype, (truth_p, comp_p) in pairs.items():
        try:
            truth = np.loadtxt(truth_p)
            computed = np.loadtxt(comp_p)
        except (OSError, ValueError) as e:
            res = MatrixComparison(inttype=inttype, error=f'cannot load: {e}')
        else:
            res = compare_matrices(truth, computed, inttype, atol=atol)
        res.truth_path, res.computed_path = Path(truth_p), Path(comp_p)
        log.info('%s: edges truth/computed/shared %d/%d/%d, %.4g %% weight diff%s',
                 inttype, res.edges_truth, res.edges_computed, res.shared,
                 res.pct_diff, f' ({res.error})' if res.error else '')
        results.append(res)
    return results


def write_comparison_csvs(
    results: List[MatrixComparison],
    out_dir: Path,
    labels:  Optional[List[str]] = None,
) -> List[Path]:
    """
    Write ``possible_sele_mismatch.csv`` and ``possible_params_mismatch.csv``
    (prove.py's files and columns, every type in one file) plus
    ``comparison_summary.csv`` (one row per type) to ``out_dir``.

    With ``labels`` (from :func:`read_label_names`), the mismatch CSVs gain
    ``label_row`` / ``label_col`` columns after the original ones.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    def _table(rows: List[dict], columns: List[str]) -> pd.DataFrame:
        df = pd.DataFrame(rows, columns=columns)
        if labels is not None:
            df['label_row'] = [labels[int(r)] if int(r) < len(labels) else ''
                               for r in df['row']]
            df['label_col'] = [labels[int(c)] if int(c) < len(labels) else ''
                               for c in df['col']]
        return df

    sele   = [row for res in results for row in res.sele_rows]
    params = [row for res in results for row in res.params_rows]
    written = [out_dir / 'possible_sele_mismatch.csv',
               out_dir / 'possible_params_mismatch.csv',
               out_dir / 'comparison_summary.csv']
    _table(sele, SELE_COLUMNS).to_csv(written[0], index=False)
    _table(params, PARAMS_COLUMNS).to_csv(written[1], index=False)
    pd.DataFrame([res.summary_row() for res in results]).to_csv(written[2], index=False)
    return written
