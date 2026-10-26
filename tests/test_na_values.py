"""Atoms and residues named "NA" survive every CSV round trip.

DEFAULT_NA_VALUES leaves "NA" out so that sodium (and the heme pyrrole
nitrogen NA, both in the CHARMM forcefields) is not read as missing. That only
works with keep_default_na=False: otherwise pandas adds its own list, which
has "NA", to na_values. 18 of the 22 read_csv calls lacked it, and the very
first one (initial_topology_info.csv) passed no na_values at all, so a
sodium residue NA was dropped from recognition by groupby('resname').
"""
import ast
from pathlib import Path

import pandas as pd

import mypsns
from mypsns.config import DEFAULT_NA_VALUES
from mypsns.services.interactors import Interactors
from mypsns.services.topology_recognition import Topology_recognition
from mypsns.utils.contacts_utilities import assemble_full_interactors

PACKAGE = Path(mypsns.__file__).parent


def read_csv_calls():
    for path in sorted(PACKAGE.rglob('*.py')):
        for node in ast.walk(ast.parse(path.read_text())):
            if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                    and node.func.attr == 'read_csv'):
                yield path.relative_to(PACKAGE.parent), node


def test_every_read_csv_keeps_na_names():
    bad = []
    for where, call in read_csv_calls():
        kw = {k.arg: k.value for k in call.keywords}
        keep = kw.get('keep_default_na')
        if 'na_values' not in kw or not (isinstance(keep, ast.Constant) and keep.value is False):
            bad.append(f'{where}:{call.lineno}')
    assert not bad, f'read_csv without na_values + keep_default_na=False: {bad}'


# Two residues whose names pandas would read as missing: a sodium ion (resname
# and atom NA) and a heme-like residue with a nitrogen NA.
PDB = """\
HETATM    1 NA    NA A   1       0.000   0.000   0.000  1.00  0.00          NA
HETATM    2 FE   HEM A   2       5.000   0.000   0.000  1.00  0.00          FE
HETATM    3 NA   HEM A   2       6.900   0.000   0.000  1.00  0.00           N
END
"""
FF_RESIDUES = """\
rtp_file,resname,atom_name,atom_type,charge,chargegroup,bonds,sigma,epsilon,mass
ions,NA,NA,SOD,1.0,1,,0.25,0.20,22.98977
heme,HEM,FE,FE,0.24,1,NA,0.05,0.0,55.847
heme,HEM,NA,NPH,-0.24,1,FE,0.33,0.84,14.007
"""
FF_TERS = "tdb_file,resname,atom_name,atom_type,mass,charge,chargegroup,bonds,sigma,epsilon\n"


def test_na_names_through_recognition_and_interactors(tmp_path):
    pdb = tmp_path / 'na.pdb'
    pdb.write_text(PDB)
    ff = tmp_path / 'forcefield'
    ff.mkdir()
    (ff / 'final_ff_residues.csv').write_text(FF_RESIDUES)
    (ff / 'final_ff_ters.csv').write_text(FF_TERS)
    topo = tmp_path / 'Topology_recognition'
    topo.mkdir()

    tr = Topology_recognition(topology=str(pdb), outpath=topo, ff_path=ff,
                              na_values=DEFAULT_NA_VALUES)
    tr.get_restypes()
    assert sorted(tr.restypes['resname']) == ['HEM', 'NA']
    tr.get_ff_parameters(maxToCheck=-1)       # raises if a residue is unmatched
    sodium = pd.read_csv(topo / 'found_topologies' / 'NA_ff_atoms.csv',
                         na_values=DEFAULT_NA_VALUES, keep_default_na=False)
    assert sodium['atom_name'].tolist() == ['NA']

    stage = Interactors(topology=str(pdb), general_outpath=tmp_path, na_values=DEFAULT_NA_VALUES)
    stage.get_complete_topology()
    stage.add_inttypes(na_values=DEFAULT_NA_VALUES)
    out = tmp_path / 'Interactors'
    full = assemble_full_interactors(out / 'found_topologies', out / 'residue_index.csv',
                                     out / 'atom_index_map.csv', na_values=DEFAULT_NA_VALUES)
    assert full['atom_name'].tolist() == ['NA', 'FE', 'NA']
    assert full['topology_resname'].tolist() == ['NA', 'HEM', 'HEM']
    heme_n = full[(full['topology_resname'] == 'HEM') & (full['atom_name'] == 'NA')]
    assert heme_n['bonds'].item() == 'FE'
    assert 'Basic' in full.loc[0, 'q_interactors']        # the ion's +1 charge
