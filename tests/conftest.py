"""Shared synthetic inputs for the interactors tests.

PHE_VARIANT is a CHARMM36 phenylalanine as the interactors stage stores it in
Interactors/found_topologies/<variant>.csv before roles are added: one residue,
bonds already two-sided, previous/next residue written as -C / +N. Its ring is
the side chain on which the q_interactors HC split used to depend on
PYTHONHASHSEED.
"""
import io

import pandas as pd
import pytest

PHE_VARIANT = """\
rtp_file,resname,atom_name,atom_type,charge,chargegroup,bonds,mass
aminoacids,PHE,N,NH1,-0.47,1,HN CA -C,14.007
aminoacids,PHE,HN,H,0.31,1,N,1.008
aminoacids,PHE,CA,CT1,0.07,1,N CB C HA,12.011
aminoacids,PHE,HA,HB1,0.09,1,CA,1.008
aminoacids,PHE,CB,CT2,-0.18,2,CA CG HB1 HB2,12.011
aminoacids,PHE,HB1,HA2,0.09,2,CB,1.008
aminoacids,PHE,HB2,HA2,0.09,2,CB,1.008
aminoacids,PHE,CG,CA,0.0,3,CB CD1 CD2,12.011
aminoacids,PHE,CD1,CA,-0.115,4,CG HD1 CE1,12.011
aminoacids,PHE,HD1,HP,0.115,4,CD1,1.008
aminoacids,PHE,CE1,CA,-0.115,5,CD1 HE1 CZ,12.011
aminoacids,PHE,HE1,HP,0.115,5,CE1,1.008
aminoacids,PHE,CZ,CA,-0.115,6,CE1 HZ CE2,12.011
aminoacids,PHE,HZ,HP,0.115,6,CZ,1.008
aminoacids,PHE,CD2,CA,-0.115,7,CG HD2 CE2,12.011
aminoacids,PHE,HD2,HP,0.115,7,CD2,1.008
aminoacids,PHE,CE2,CA,-0.115,8,CD2 HE2 CZ,12.011
aminoacids,PHE,HE2,HP,0.115,8,CE2,1.008
aminoacids,PHE,C,C,0.51,9,CA O +N,12.011
aminoacids,PHE,O,O,-0.51,9,C,15.9994
"""


def variant_frame(csv_text: str = PHE_VARIANT) -> pd.DataFrame:
    """A variant table with the scaffolding columns add_inttypes adds before
    typing (single residue, no neighbours)."""
    df = pd.read_csv(io.StringIO(csv_text), keep_default_na=False)
    df['Resindex'] = 0
    df['Entity'] = 'M'
    return df


@pytest.fixture
def phe_variant() -> pd.DataFrame:
    return variant_frame()
