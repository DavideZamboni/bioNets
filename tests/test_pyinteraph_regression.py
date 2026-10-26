"""Regression test against PyInteraph 1 on its CHARMM test case (149 residues).

Runs the whole pipeline (ff-reader → topology → interactors → contacts) on the
inputs shipped in tests/data/pyinteraph_testcase/ into a temporary folder, with
that preset's config.yaml (named_interactors, default cutoffs), then compares
the adjacency matrices with PyInteraph's own hb/hc/sb-graph.dat through the
`mypsns compare` logic. Nothing outside the repository is read.

Slow (~3 min, most of it ff-reader parsing CHARMM36): deselect with
`pytest -m "not slow"`.
"""
import os

import pytest
import yaml
from typer.testing import CliRunner

from mypsns.cli import app
from mypsns.utils.validation_utilities import (
    TEST_DATA_DIR, adjacency_paths, compare_files, get_preset, preset_paths,
)

pytestmark = pytest.mark.slow

DATA   = TEST_DATA_DIR / 'pyinteraph_testcase'
PRESET = get_preset('pyinteraph_testcase')


@pytest.fixture(scope='module')
def results(tmp_path_factory):
    out = tmp_path_factory.mktemp('pyinteraph_testcase')
    cfg = yaml.safe_load(PRESET.config.read_text())
    cfg.update(
        topology=str(DATA / 'structures' / 'sim.prot.gro'),
        trajectory=str(DATA / 'structures' / 'traj.xtc'),
        forcefield=str(DATA / 'forcefields' / 'charmm36-apr2019_vm.ff'),
        outpath=str(out),
    )
    cfg.setdefault('contacts', {})['nc'] = min(4, os.cpu_count() or 1)
    cfg_file = out / 'config.yaml'
    cfg_file.write_text(yaml.safe_dump(cfg))

    run = CliRunner().invoke(app, ['run', '--config', str(cfg_file)])
    assert run.exit_code == 0, run.output

    truth    = preset_paths('pyinteraph_testcase')
    computed = adjacency_paths(out / 'Adjacency' / 'named_interactors')
    return {r.inttype: r for r in compare_files({t: (truth[t], computed[t]) for t in truth})}


@pytest.mark.parametrize('inttype', sorted(PRESET.max_pct_diff))
def test_weight_difference_within_tolerance(results, inttype):
    # SB/HC sit at 0.07 % / 0.11 %. H-bonds at 3.93 %: the detection engine
    # (MDAnalysis HydrogenBondAnalysis vs PyInteraph's geometric test), not a
    # parameter. Limits are the preset's, shared with `mypsns compare`.
    res, limit = results[inttype], PRESET.max_pct_diff[inttype]
    assert res.error is None, res.error
    assert res.pct_diff <= limit, (
        f'{inttype}: {res.pct_diff:.4f} % weight difference vs PyInteraph '
        f'(limit {limit} %); edges truth/computed/shared '
        f'{res.edges_truth}/{res.edges_computed}/{res.shared}.')


@pytest.mark.parametrize('inttype', sorted(PRESET.max_pct_diff))
def test_every_pyinteraph_edge_is_found(results, inttype):
    # No PyInteraph-only edge. For H-bonds this needs d_h_cutoff ≥ 1.32 Å: with
    # MDAnalysis' 1.2 Å default the Cys86 thiol (S–H 1.32 Å) never donates and
    # Cys86–Asp118/Pro119 and one more pair go missing.
    res = results[inttype]
    assert res.only_truth == 0, (
        f'{inttype}: {res.only_truth} PyInteraph edge(s) missing from myPSNs')


def test_salt_bridges_and_hydrophobic_edges_are_identical(results):
    for inttype in ('SB', 'HC'):
        res = results[inttype]
        assert res.edges_truth == res.edges_computed == res.shared, inttype
