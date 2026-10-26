"""Unit tests for mypsns/utils/validation_utilities.py and `mypsns compare`
(the replacement of the old ad-hoc prove.py), on small synthetic matrices."""
import math

import numpy as np
import pandas as pd
import pytest
from typer.testing import CliRunner

from mypsns.cli import app
from mypsns.utils.validation_utilities import (
    PARAMS_COLUMNS, SELE_COLUMNS, VALIDATION_ROOT_ENV,
    compare_files, compare_matrices, preset_paths, read_label_names,
    write_comparison_csvs,
)


def _sym(n, entries, diag=None):
    m = np.zeros((n, n))
    for (i, j), v in entries.items():
        m[i, j] = m[j, i] = v
    if diag is not None:
        np.fill_diagonal(m, diag)
    return m


# Pairs (row > col, lower triangle):
#   (1,0) shared, equal          (2,0) shared, differs by 5
#   (3,0) truth only             (3,2) computed only
#   (2,1) shared, within atol    diagonal differs but must be ignored
TRUTH    = _sym(4, {(1, 0): 50.0, (2, 0): 20.0, (3, 0): 10.0, (2, 1): 30.0}, diag=7.0)
COMPUTED = _sym(4, {(1, 0): 50.0, (2, 0): 25.0, (3, 2): 4.0,  (2, 1): 30.005})


def test_edge_counts_ignore_the_diagonal():
    res = compare_matrices(TRUTH, COMPUTED, 'X')
    assert res.error is None
    assert (res.edges_truth, res.edges_computed, res.shared) == (4, 4, 3)
    assert (res.only_truth, res.only_computed) == (1, 1)
    assert res.value_mismatch == 1


def test_weight_difference_is_relative_to_truth_total():
    res = compare_matrices(TRUTH, COMPUTED, 'X')
    assert res.truth_total == pytest.approx(110.0)
    assert res.computed_total == pytest.approx(109.005)
    # |25-20| + |0-10| + |4-0| + |30.005-30|
    assert res.abs_diff == pytest.approx(19.005)
    assert res.pct_diff == pytest.approx(100 * 19.005 / 110.0)


def test_mismatch_rows_match_prove_py_layout():
    res = compare_matrices(TRUTH, COMPUTED, 'X')
    assert res.sele_rows == [
        {'row': 3, 'col': 0, 'intertype': 'X', 'present_in': 'true',     'non-0 val': 10.0},
        {'row': 3, 'col': 2, 'intertype': 'X', 'present_in': 'computed', 'non-0 val': 4.0},
    ]
    assert res.params_rows == [
        {'row': 2, 'col': 0, 'intertype': 'X', 'difference': -5.0},
    ]


def test_atol_controls_what_counts_as_different():
    res = compare_matrices(TRUTH, COMPUTED, 'X', atol=1e-4)
    assert res.value_mismatch == 2                 # now (2,1) differs too


def test_identical_matrices():
    res = compare_matrices(TRUTH, TRUTH, 'X')
    assert res.abs_diff == 0.0 and res.pct_diff == 0.0
    assert res.sele_rows == [] and res.params_rows == []


def test_shape_mismatch_is_reported_not_raised():
    res = compare_matrices(np.zeros((3, 3)), np.zeros((4, 4)), 'X')
    assert 'shape mismatch' in res.error
    assert res.edges_truth == 0


def test_empty_truth_gives_nan_pct():
    assert math.isnan(compare_matrices(np.zeros((2, 2)), np.eye(2), 'X').pct_diff)


def test_compare_files_reports_unreadable_input(tmp_path):
    good = tmp_path / 'good.dat'
    np.savetxt(good, TRUTH)
    res, = compare_files({'X': (tmp_path / 'missing.dat', good)})
    assert res.error.startswith('cannot load')


def test_write_csvs_columns_and_labels(tmp_path):
    res = compare_matrices(TRUTH, COMPUTED, 'X')
    labels_file = tmp_path / 'adjacency_labels.dat'
    labels_file.write_text('0\tA:ALA:0\n1\tA:GLY:1\n2\tA:LYS:2\n3\tA:ASP:3\n')

    plain = write_comparison_csvs([res], tmp_path / 'plain')
    assert list(pd.read_csv(plain[0]).columns) == SELE_COLUMNS
    assert list(pd.read_csv(plain[1]).columns) == PARAMS_COLUMNS
    summary = pd.read_csv(plain[2])
    assert summary.loc[0, 'shared'] == 3

    labelled = write_comparison_csvs([res], tmp_path / 'lab',
                                     labels=read_label_names(labels_file))
    sele = pd.read_csv(labelled[0])
    assert list(sele.columns) == SELE_COLUMNS + ['label_row', 'label_col']
    assert sele.loc[0, 'label_row'] == 'A:ASP:3' and sele.loc[0, 'label_col'] == 'A:ALA:0'


def test_preset_paths_use_env_root(tmp_path, monkeypatch):
    monkeypatch.setenv(VALIDATION_ROOT_ENV, str(tmp_path))
    truth, computed = preset_paths('rec_grk1')
    assert truth['HBonds'] == tmp_path / 'risultati_paragone/rec/hb-graph_1micros.dat'
    assert computed['SB'] == (tmp_path / 'psns_out/cli_outs/rec_grk1/Adjacency/'
                              'named_interactors/SB_adjacency.dat')
    with pytest.raises(ValueError, match='unknown preset'):
        preset_paths('nope')

# ═════════════════════════════════════════════════════════════════════════════
# CLI wiring
# ═════════════════════════════════════════════════════════════════════════════

@pytest.fixture
def mats(tmp_path):
    t, c = tmp_path / 't.dat', tmp_path / 'c.dat'
    np.savetxt(t, TRUTH)
    np.savetxt(c, COMPUTED)
    return t, c


def test_cli_writes_outputs_and_passes(mats, tmp_path):
    t, c = mats
    out = tmp_path / 'out'
    result = CliRunner().invoke(app, [
        'compare', '--truth', f'X={t}', '--computed', f'X={c}',
        '--out-dir', str(out), '--max-pct-diff', 'X=20',
    ])
    assert result.exit_code == 0, result.output
    assert (out / 'possible_sele_mismatch.csv').exists()
    assert (out / 'comparison_summary.csv').exists()


def test_cli_fails_above_threshold(mats):
    t, c = mats
    result = CliRunner().invoke(app, [
        'compare', '--truth', f'X={t}', '--computed', f'X={c}',
        '--max-pct-diff', 'X=1',
    ])
    assert result.exit_code == 1


def test_cli_rejects_missing_inputs(tmp_path):
    result = CliRunner().invoke(app, [
        'compare', '--truth', f'X={tmp_path / "a.dat"}',
        '--computed', f'X={tmp_path / "b.dat"}',
    ])
    assert result.exit_code == 2
