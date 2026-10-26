"""The default HC rule takes a terminal residue's side chain the way PyInteraph
does: every atom of the residue except CA C O N H H1 H2 H3 O1 O2 OXT OT1 OT2.

Terminal patch atoms carry the patch resname (COOH, NH2, …), so the rule
selects by `residue`, the residue each atom belongs to, set by
Interactors.add_inttypes through `residue_name`. The tables below are the
C-terminal Leu (COOH) and N-terminal Met (NH2) of the rec_grk1 run, as
interactors stores them before adding roles.
"""
import io

import pandas as pd

from mypsns.utils.interactors_utilities import get_named_interactors, residue_name

COLUMNS = 'rtp_file,resname,atom_name,atom_type,charge,chargegroup,bonds,mass\n'

LEU_COOH = COLUMNS + """\
merged,COOH,C,CD,0.72,17,CA OT1 OT2,12.011
merged,LEU,CA,CT1,0.07,2,N CB C HA,12.011
merged,LEU,CB,CT2,-0.18,4,CA CG HB1 HB2,12.011
merged,LEU,CD1,CT3,-0.27,9,CG HD11 HD12 HD13,12.011
merged,LEU,CD2,CT3,-0.27,13,CG HD21 HD22 HD23,12.011
merged,LEU,CG,CT1,-0.09,7,CB CD1 CD2 HG,12.011
merged,LEU,HA,HB1,0.09,3,CA,1.008
merged,LEU,HB1,HA2,0.09,5,CB,1.008
merged,LEU,HB2,HA2,0.09,6,CB,1.008
merged,LEU,HD11,HA3,0.09,10,CD1,1.008
merged,LEU,HD12,HA3,0.09,11,CD1,1.008
merged,LEU,HD13,HA3,0.09,12,CD1,1.008
merged,LEU,HD21,HA3,0.09,14,CD2,1.008
merged,LEU,HD22,HA3,0.09,15,CD2,1.008
merged,LEU,HD23,HA3,0.09,16,CD2,1.008
merged,LEU,HG,HA1,0.09,8,CG,1.008
merged,LEU,HN,H,0.31,1,N,1.008
merged,COOH,HT2,H,0.44,17,OT2,1.008
merged,LEU,N,NH1,-0.47,0,HN CA -C,14.007
merged,COOH,OT1,OB,-0.55,17,C,15.9994
merged,COOH,OT2,OH1,-0.61,17,C HT2,15.9994
"""

MET_NH2 = COLUMNS + """\
merged,MET,C,C,0.51,15,CA +N O,12.011
merged,NH2,CA,CT1,0.19,17,N CB C HA,12.011
merged,MET,CB,CT2,-0.18,4,CA CG HB1 HB2,12.011
merged,MET,CE,CT3,-0.22,11,SD HE1 HE2 HE3,12.011
merged,MET,CG,CT2,-0.14,7,CB SD HG1 HG2,12.011
merged,NH2,H1,H,0.34,17,N,1.008
merged,NH2,H2,H,0.34,17,N,1.008
merged,NH2,HA,HB1,0.09,17,CA,1.008
merged,MET,HB1,HA2,0.09,5,CB,1.008
merged,MET,HB2,HA2,0.09,6,CB,1.008
merged,MET,HE1,HA3,0.09,12,CE,1.008
merged,MET,HE2,HA3,0.09,13,CE,1.008
merged,MET,HE3,HA3,0.09,14,CE,1.008
merged,MET,HG1,HA2,0.09,8,CG,1.008
merged,MET,HG2,HA2,0.09,9,CG,1.008
merged,NH2,N,NH2,-0.96,17,CA H1 H2,14.007
merged,MET,O,O,-0.51,16,C,15.999
merged,MET,SD,S,-0.09,10,CG CE,32.06
"""

PATCHES = {'COO-', 'COOH', 'NH2', 'NH3+', 'GLY-NH3+', 'PRO-NH2+'}
# PyInteraph 1.1 generate_sc_identifiers: the side chain is everything else
PYINTERAPH_BACKBONE = {'CA', 'C', 'O', 'N', 'H', 'H1', 'H2', 'H3',
                       'O1', 'O2', 'OXT', 'OT1', 'OT2'}


def _hc_atoms(text):
    df = pd.read_csv(io.StringIO(text), keep_default_na=False)
    df['Resindex'] = 0
    df['residue'] = residue_name(df, PATCHES)
    out = get_named_interactors(df)
    hc = out[out['named_interactors'].str.contains('HC')]
    return set(hc['atom_name']), set(df['atom_name'])


def test_terminal_side_chain_matches_pyinteraph():
    for text in (LEU_COOH, MET_NH2):
        hc, every = _hc_atoms(text)
        assert hc == every - PYINTERAPH_BACKBONE


def test_patch_atoms_are_included():
    assert 'HT2' in _hc_atoms(LEU_COOH)[0]        # COOH hydrogen
    assert 'HA' in _hc_atoms(MET_NH2)[0]          # HA of the NH2 patch


def test_residue_name_ignores_patch_rows():
    # A Gly whose NH3+ patch outnumbers its own rows is still GLY.
    gly = pd.DataFrame({'resname': ['GLY-NH3+'] * 6 + ['GLY'] * 2})
    assert residue_name(gly, PATCHES) == 'GLY'
    met = pd.read_csv(io.StringIO(MET_NH2), keep_default_na=False)
    assert residue_name(met, PATCHES) == 'MET'
