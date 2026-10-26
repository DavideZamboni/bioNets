"""Terminal patches (.tdb entries) applied by topology recognition, and the
structure_atom_name column of found_topologies/.

get_ter_ffres used to edit only the replaced atoms' own rows: renamed and
deleted names survived in the other rows' bonds, an added atom whose control
atom did not list it kept the patch's "-X" delete instructions as bonds, and
the hand-made NHT3+ of psns_out/{wt,M654V}_noh2o reached found_topologies with
N bonded to H1 H2 H3 (atoms it does not have) and HT1-HT3 bonded to nothing.
"""
import io

import pandas as pd
import pytest

from mypsns.config import DEFAULT_NA_VALUES
from mypsns.services.topology_recognition import accept_topologies
from mypsns.utils import ff_reader_utilities as ffr
from mypsns.utils.interactors_utilities import build_wholetop
from mypsns.utils.topology_recognition_utilities import (
    STRUCTURE_NAME_COLUMN, bond_tokens, dangling_bonds, get_res_w_ter,
    save_results, structure_named, unbonded_patch_atoms,
)

NA = DEFAULT_NA_VALUES

# A CHARMM-like Ala as ff_reader writes it to final_ff_residues.csv.
ALA = """\
rtp_file,resname,atom_name,atom_type,charge,chargegroup,bonds,sigma,epsilon,mass
aa,ALA,N,NH1,-0.47,1,HN CA,0.33,0.84,14.007
aa,ALA,HN,H,0.31,1,N,0.04,0.19,1.008
aa,ALA,CA,CT1,0.07,1,N CB C HA,0.36,0.13,12.011
aa,ALA,HA,HB1,0.09,1,CA,0.24,0.09,1.008
aa,ALA,CB,CT3,-0.27,2,CA HB1 HB2 HB3,0.36,0.33,12.011
aa,ALA,HB1,HA3,0.09,2,CB,0.24,0.10,1.008
aa,ALA,HB2,HA3,0.09,2,CB,0.24,0.10,1.008
aa,ALA,HB3,HA3,0.09,2,CB,0.24,0.10,1.008
aa,ALA,C,C,0.51,3,CA +N O,0.36,0.46,12.011
aa,ALA,O,O,-0.51,3,C,0.30,0.50,15.9994
"""

TER_HEADER = "tdb_file,resname,atom_name,atom_type,mass,charge,chargegroup,bonds,sigma,epsilon\n"


def frame(text):
    df = pd.read_csv(io.StringIO(text), na_values=NA, keep_default_na=False)
    df['bonds'] = df['bonds'].fillna('')
    return df


def ter(rows):
    return frame(TER_HEADER + rows)


def bonds_of(df):
    return {a: set(bond_tokens(b)) for a, b in zip(df['atom_name'], df['bonds'])}


def assert_consistent(df):
    '''No patch instruction left, every bond token names an atom (or the
    previous / next residue), and every intra-residue bond is two-sided.'''
    bonds = bonds_of(df)
    for atom, toks in bonds.items():
        for t in toks:
            assert not t.startswith('r'), (atom, t)
            if t[0] in '+-':
                assert atom in ('C', 'N'), f'{atom} lists {t}: a patch instruction left as a bond'
                continue
            assert t in bonds, (atom, t)
            assert atom in bonds[t], f'{atom}->{t} is one-sided'


# ── patches as the existing final_ff_ters.csv files store them ────────────────

# jul2022 NH3+: the added H1-H3 carry only the delete instructions, because
# ff_reader's substring check found "N" in "-HN" and never wrote their bond.
NH3P_OLD_FORMAT = ter("""\
aa,NH3+,CA,CT1,12.011,0.21,-1,-HN -HN1 rCA,0.36,0.13
aa,NH3+,H1,HC,1.008,0.33,-1,-HN -HN1,0.04,0.19
aa,NH3+,H2,HC,1.008,0.33,-1,-HN -HN1,0.04,0.19
aa,NH3+,H3,HC,1.008,0.33,-1,-HN -HN1,0.04,0.19
aa,NH3+,N,NH3,14.007,-0.3,-1,-HN -HN1 CA H1 H2 H3 rN,0.33,0.84
""")


def test_nterm_patch_bonds_added_atoms_and_drops_instructions():
    out = get_res_w_ter(NH3P_OLD_FORMAT, frame(ALA))
    bonds = bonds_of(out)
    assert set(out['atom_name']) == {'N', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'HB3',
                                     'C', 'O', 'H1', 'H2', 'H3'}
    assert bonds['N'] == {'CA', 'H1', 'H2', 'H3'}
    assert bonds['H1'] == bonds['H2'] == bonds['H3'] == {'N'}
    assert_consistent(out)
    # rows the patch does not change keep their cell byte for byte
    ala = frame(ALA).set_index('atom_name')['bonds']
    kept = out.set_index('atom_name')['bonds']
    for atom in ('CA', 'CB', 'C', 'O'):
        assert kept[atom] == ala[atom]
    assert unbonded_patch_atoms(out) == []


def test_rename_and_delete_reach_every_bonds_list():
    # The patch renames O -> OT1 and deletes HB3, but replaces neither C nor
    # CB: their rows must still follow the rename / deletion.
    patch = ter("""\
aa,XTER,OT1,OC,15.9994,-0.67,-1,-HB3 rO,0.30,0.50
aa,XTER,OT2,OC,15.9994,-0.67,-1,-HB3 C,0.30,0.50
""")
    out = get_res_w_ter(patch, frame(ALA))
    bonds = bonds_of(out)
    assert 'O' not in bonds and 'HB3' not in bonds
    assert bonds['C'] == {'CA', '+N', 'OT1', 'OT2'}   # C not replaced: keeps +N
    assert bonds['CB'] == {'CA', 'HB1', 'HB2'}
    assert bonds['OT1'] == bonds['OT2'] == {'C'}
    assert_consistent(out)


def test_cterm_patch_bond_listed_by_the_added_atom_only():
    # apr2019 COOH: OT2 lists C only, HT2 lists OT2. The old code overwrote
    # OT2's bonds with "C" and lost OT2-HT2 (interactors' symmetrisation
    # restored it later); the patched topology itself must have it.
    patch = ter("""\
aa,COOH,C,CD,12.011,0.72,-1,CA OT1 OT2 rC,0.36,0.29
aa,COOH,HT2,H,1.008,0.44,-1,OT2,0.04,0.19
aa,COOH,OT1,OB,15.9994,-0.55,-1,C rO,0.30,0.50
aa,COOH,OT2,OH1,15.9994,-0.61,-1,C,0.32,0.64
""")
    out = get_res_w_ter(patch, frame(ALA))
    bonds = bonds_of(out)
    assert bonds['C'] == {'CA', 'OT1', 'OT2'}        # replaced: +N dropped
    assert bonds['OT2'] == {'C', 'HT2'}
    assert bonds['HT2'] == {'OT2'}
    assert_consistent(out)


def test_hand_made_nht3_patch_is_reported_not_kept():
    # psns_out/{wt,M654V}_noh2o: NH3+ copied with the added atoms renamed
    # HT1-HT3 but N still listing H1 H2 H3.
    nht3 = NH3P_OLD_FORMAT.replace({'NH3+': 'NHT3+', 'H1': 'HT1', 'H2': 'HT2', 'H3': 'HT3'})
    out = get_res_w_ter(nht3, frame(ALA))
    bonds = bonds_of(out)
    assert bonds['N'] == {'CA'}                      # H1-H3 are not atoms: dropped
    assert bonds['HT1'] == bonds['HT2'] == bonds['HT3'] == set()   # no "-HN"
    assert dangling_bonds(out) == {}
    assert unbonded_patch_atoms(out) == ['HT1', 'HT2', 'HT3']

    # the same entry with its bonds fixed (as fix_nht3 in the scratch notes)
    fixed = nht3.copy()
    fixed.loc[fixed['atom_name'] == 'N', 'bonds'] = '-HN -HN1 CA HT1 HT2 HT3 rN'
    out = get_res_w_ter(fixed, frame(ALA))
    assert bonds_of(out)['N'] == {'CA', 'HT1', 'HT2', 'HT3'}
    assert bonds_of(out)['HT1'] == {'N'}
    assert unbonded_patch_atoms(out) == []


# ── a .tdb file, parsed by ff_reader and applied by recognition ───────────────

N_TDB = """\
; test termini database
[ None ]

[ NH3+ ]
[ delete ]
 HN
 HT1
[ replace ]
N      NH3       14.007000  -0.3000
CA     CT1       12.011000   0.2100
[ add ]
3   4   H     N      CA     C
  HC         1.008000   0.3300  -1
"""

C_TDB = """\
[ COO- ]
[ replace ]
C      CC        12.011000   0.3400
O      OT1  OC   15.999400  -0.6700
[ add ]
1   8   OT2   C      CA     N
  OC        15.999400  -0.6700  -1
"""


def parse_tdb(tmp_path, name, text):
    '''ff_reader's path for one .tdb file, round-tripped through CSV like
    final_ff_ters.csv (sigma/epsilon come from the atom types there).'''
    (tmp_path / name).write_text(text)
    residues = {k: [] for k in ('tdb_file', 'resname', 'atom_name', 'atom_type',
                                'mass', 'charge', 'chargegroup', 'bonds')}
    parsed = ffr.merge_duplicate_atoms(pd.DataFrame(ffr.read_tbd_file(tmp_path, residues, name)))
    parsed['sigma'] = parsed['epsilon'] = 0.1
    return frame(parsed.to_csv(index=False))


def test_tdb_patch_parsed_and_applied(tmp_path):
    nter = parse_tdb(tmp_path, 'test.n.tdb', N_TDB)
    # ff_reader now writes the [ add ] bond in the added atom's row too
    h1 = nter.loc[nter['atom_name'] == 'H1', 'bonds'].item().split()
    assert 'N' in h1 and {'-HN', '-HT1'} <= set(h1)

    out = get_res_w_ter(nter[nter['resname'] == 'NH3+'], frame(ALA))
    bonds = bonds_of(out)
    assert 'HN' not in bonds
    assert bonds['N'] == {'CA', 'H1', 'H2', 'H3'}
    assert bonds['H1'] == {'N'}
    assert_consistent(out)

    cter = parse_tdb(tmp_path, 'test.c.tdb', C_TDB)
    out = get_res_w_ter(cter[cter['resname'] == 'COO-'], frame(ALA))
    bonds = bonds_of(out)
    assert 'O' not in bonds                          # renamed by [ replace ]
    assert bonds['C'] == {'CA', 'OT1', 'OT2'}
    assert bonds['OT1'] == bonds['OT2'] == {'C'}
    assert out.loc[out['atom_name'] == 'OT1', 'charge'].item() == pytest.approx(-0.67)
    assert_consistent(out)


def test_ff_reader_bond_checks_are_per_token():
    assert ffr.has_bond_token('-HN -HN1 rN CA', 'CA')
    assert ffr.has_bond_token('-HN -HN1 rN', '-HN1')
    assert not ffr.has_bond_token('-HN -HN1 rN', 'N')    # was True as a substring
    assert not ffr.has_bond_token('-HT1 -HT2', 'HT1')


# ── structure_atom_name ──────────────────────────────────────────────────────

def recognition_folder(tmp_path, key, atoms):
    out = tmp_path / 'Topology_recognition'
    (out / 'resid_topologies').mkdir(parents=True)
    pd.DataFrame({'resname': [key], 'atoms': [f'{key}_atoms.csv'], 'rixs': [f'{key}_rixs.txt']}
                 ).to_csv(out / 'different_residue_topologies.csv', index=False)
    atoms[['atom_name', 'mass']].rename(columns={'atom_name': 'name'}).to_csv(
        out / 'resid_topologies' / f'{key}_atoms.csv', index=False)
    (out / 'resid_topologies' / f'{key}_rixs.txt').write_text('0 3')
    return out


def test_found_topologies_have_structure_names(tmp_path):
    ala = frame(ALA)
    out = recognition_folder(tmp_path, 'ALA', ala)
    complete, _ = save_results(out, {'ALA': ala}, {}, NA)
    assert complete
    found = pd.read_csv(out / 'found_topologies' / 'ALA_ff_atoms.csv', keep_default_na=False)
    cols = found.columns.tolist()
    assert cols[cols.index('atom_name') + 1] == STRUCTURE_NAME_COLUMN
    assert (found[STRUCTURE_NAME_COLUMN] == found['atom_name']).all()


def test_recognition_warns_about_unbonded_patch_atoms(tmp_path, caplog):
    nht3 = NH3P_OLD_FORMAT.replace({'NH3+': 'NHT3+', 'H1': 'HT1', 'H2': 'HT2', 'H3': 'HT3'})
    patched = get_res_w_ter(nht3, frame(ALA))
    out = recognition_folder(tmp_path, 'ALA', patched)
    with caplog.at_level('WARNING'):
        save_results(out, {'ALA': patched}, {}, NA)
    assert 'NHT3+' in caplog.text and 'HT1 HT2 HT3' in caplog.text


GLY_FF = """\
rtp_file,resname,atom_name,atom_type,charge,chargegroup,bonds,sigma,epsilon,mass
aa,GLY,N,NH1,-0.47,1,HN CA,0.33,0.84,14.007
aa,GLY,HN,H,0.31,1,N,0.04,0.19,1.008
aa,GLY,CA,CT2,-0.02,2,N C HA1 HA2,0.36,0.23,12.011
aa,GLY,HA1,HB2,0.09,2,CA,0.24,0.12,1.008
aa,GLY,HA2,HB2,0.09,2,CA,0.24,0.12,1.008
aa,GLY,C,C,0.51,3,CA O +N,0.36,0.46,12.011
aa,GLY,O,O,-0.51,3,C,0.30,0.50,15.9994
"""


def test_accept_topology_keeps_ff_names_and_records_structure_names(tmp_path):
    # the structure calls HA1 "HA3"; recognition halted and proposed GLY_FF
    out = tmp_path / 'Topology_recognition'
    for d in ('resid_topologies', 'possible_topologies', 'found_topologies'):
        (out / d).mkdir(parents=True)
    pd.DataFrame({'name': ['N', 'HN', 'CA', 'HA3', 'HA2', 'C', 'O'],
                  'mass': [14.007, 1.008, 12.011, 1.008, 1.008, 12.011, 15.9994]}
                 ).to_csv(out / 'resid_topologies' / 'GLY_atoms.csv', index=False)
    (out / 'resid_topologies' / 'GLY_rixs.txt').write_text('2 5')
    (out / 'possible_topologies' / 'GLY_uncertain_atoms_0.csv').write_text(GLY_FF)
    pd.DataFrame({'resname': ['GLY'], 'atoms': ['GLY_atoms.csv'], 'rixs': ['GLY_rixs.txt'],
                  'ff_atoms': [''], 'possible_atoms': ['GLY_uncertain_atoms_0.csv '],
                  'unmatching_masses': [0]}
                 ).to_csv(out / 'unmatching_residue_topologies.csv', index=False)

    result = accept_topologies(out, {'GLY': '0'}, NA)
    assert result['complete'] and result['accepted'][0]['renames'] == {'HA1': 'HA3'}

    found = pd.read_csv(out / 'found_topologies' / 'GLY_ff_atoms.csv', keep_default_na=False)
    assert found['atom_name'].tolist() == ['N', 'HN', 'CA', 'HA1', 'HA2', 'C', 'O']
    assert found[STRUCTURE_NAME_COLUMN].tolist() == ['N', 'HN', 'CA', 'HA3', 'HA2', 'C', 'O']
    assert found.loc[found['atom_name'] == 'CA', 'bonds'].item() == 'N C HA1 HA2'

    # what the interactors stage builds from it: the structure's names
    named = structure_named(found)
    assert STRUCTURE_NAME_COLUMN not in named.columns
    assert named['atom_name'].tolist() == ['N', 'HN', 'CA', 'HA3', 'HA2', 'C', 'O']
    assert named.loc[named['atom_name'] == 'CA', 'bonds'].item() == 'N C HA3 HA2'

    resid_to_top = pd.DataFrame({'Resindex': [2, 5], 'topology': ['GLY_ff_atoms.csv'] * 2})
    wholetop, allres = build_wholetop(
        {'GLY_ff_atoms.csv': out / 'found_topologies' / 'GLY_ff_atoms.csv'}, resid_to_top, NA)
    assert allres == [2, 5]
    assert STRUCTURE_NAME_COLUMN not in wholetop.columns
    assert (wholetop['atom_name'] == 'HA3').sum() == 2 and 'HA1' not in set(wholetop['atom_name'])


def test_structure_named_leaves_old_files_alone():
    # found_topologies written before the column existed are used as they are
    ala = frame(ALA)
    assert structure_named(ala) is ala
