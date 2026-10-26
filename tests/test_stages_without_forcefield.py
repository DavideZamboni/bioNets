"""`interactors` and `contacts` run from flags alone, without a forcefield.

Neither reads the forcefield path (interactors works from
Topology_recognition/, contacts from <outpath>/forcefield/ and Interactors/),
but they used to inherit the `forcefield: Field required` error, and neither has
a --forcefield flag. That only showed when `topology` was the first command
of a folder and halted, which leaves no run_config.yaml behind.
"""
from typer.testing import CliRunner

from mypsns.cli import app

runner = CliRunner()


def run_folder(tmp_path):
    (tmp_path / 'Topology_recognition').mkdir()
    (tmp_path / 'Interactors').mkdir()
    (tmp_path / 's.pdb').write_text('END\n')
    (tmp_path / 't.xtc').write_text('')
    return tmp_path


def test_interactors_from_flags(tmp_path):
    run = run_folder(tmp_path)
    res = runner.invoke(app, ['interactors', '-t', str(run / 's.pdb'), '-o', str(run), '--dry-run'])
    assert res.exit_code == 0, res.output
    assert 'Field required' not in res.output


def test_contacts_from_flags(tmp_path):
    run = run_folder(tmp_path)
    res = runner.invoke(app, ['contacts', '-t', str(run / 's.pdb'), '--trajectory', str(run / 't.xtc'),
                              '-o', str(run), '--dry-run'])
    assert res.exit_code == 0, res.output
    assert 'Field required' not in res.output


def test_topology_still_needs_the_forcefield(tmp_path):
    run = run_folder(tmp_path)
    res = runner.invoke(app, ['topology', '-t', str(run / 's.pdb'), '-o', str(run), '--dry-run'])
    assert res.exit_code != 0
