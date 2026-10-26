"""Regression test against PyInteraph 1 on its CHARMM test case (149 residues).

Compares PyInteraph's hb/hc/sb-graph.dat with the myPSNs contacts output of the
same run (named_interactors, default cutoffs) — the validation recorded in
CLAUDE.md. The data are not in the repo: the test is skipped when they are
missing. Point $MYPSNS_VALIDATION_ROOT at another copy of the lab folder to
run it elsewhere (layout in validation_utilities.PRESETS).

Note this checks the matrices on disk, so re-run `contacts` on the test case
(into that run folder) after changing contact code, then run the test.
"""
import pytest

from mypsns.utils.validation_utilities import compare_files, preset_paths

# Tolerances on the % weight difference. SB/HC sit at 0.07 % / 0.11 %; H-bonds
# plateau at ~4.56 % because of the detection engine (MDAnalysis
# HydrogenBondAnalysis vs PyInteraph's geometric test), not a parameter.
MAX_PCT_DIFF = {'SB': 0.2, 'HC': 0.2, 'HBonds': 5.0}

_truth, _computed = preset_paths('pyinteraph_testcase')
PAIRS = {t: (_truth[t], _computed[t]) for t in MAX_PCT_DIFF}
MISSING = [str(p) for pair in PAIRS.values() for p in pair if not p.exists()]

pytestmark = pytest.mark.skipif(
    bool(MISSING), reason=f'PyInteraph test-case data not found: {MISSING}')


@pytest.fixture(scope='module')
def results():
    return {res.inttype: res for res in compare_files(PAIRS)}


@pytest.mark.parametrize('inttype', list(MAX_PCT_DIFF))
def test_weight_difference_within_tolerance(results, inttype):
    res = results[inttype]
    assert res.error is None, res.error
    assert res.pct_diff <= MAX_PCT_DIFF[inttype], (
        f'{inttype}: {res.pct_diff:.4f} % weight difference vs PyInteraph '
        f'(limit {MAX_PCT_DIFF[inttype]} %); edges truth/computed/shared '
        f'{res.edges_truth}/{res.edges_computed}/{res.shared}. '
        f'Inspect with `mypsns compare --preset pyinteraph_testcase -o <dir>`.')
