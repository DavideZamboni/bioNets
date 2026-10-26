"""build_results must hand every merged contact file of a type to the
aggregation, inter-chain ones included. It used to pick H-bond files by the
'_Hbonds_merged.csv' suffix, which only the intra-chain files carry, so on a
multi-chain system (rec_grk1: GRK1 peptide ↔ recoverin) every inter-chain
H-bond was missing from Hbonds_final.csv and the adjacency matrices."""
from mypsns.utils import contacts_utilities as cu

# File names exactly as the Salt_bridges / Hbonds / hydrophobics wrappers
# write them, for subject chain A with partner chain C.
WRITTEN = {
    'H_bonds':      ['A_intrachain_Hbonds', 'A_HDon_C_HAcc', 'A_HAcc_C_HDon'],
    'Salt_bridges': ['A_intrachain_sb', 'A_bas_C_ac_sb', 'A_ac_C_bas_sb'],
    'HCs':          ['A_intrachain_HC', 'A_C_HC'],
}


def test_build_results_keeps_inter_chain_files(tmp_path, monkeypatch):
    merged = []
    for folder, stems in WRITTEN.items():
        (tmp_path / folder).mkdir()
        for stem in stems:
            f = tmp_path / folder / f'{stem}_merged.csv'
            f.write_text('')
            merged.append(f)

    seen = {}

    def fake_hbonds(files, *args, **kwargs):
        seen['HBonds'] = sorted(f.stem for f in files)
        return tmp_path

    def fake_pairs(files, interactors, outp, label, **kwargs):
        seen[label] = sorted(f.stem for f in files)
        return tmp_path

    monkeypatch.setattr(cu, 'load_hbond_contacts', fake_hbonds)
    monkeypatch.setattr(cu, 'load_pair_contacts', fake_pairs)
    cu.build_results(merged, interactors=None, frames=set(), outp=tmp_path)

    expected = {key: sorted(f'{s}_merged' for s in WRITTEN[folder])
                for key, folder in [('HBonds', 'H_bonds'), ('SB', 'Salt_bridges'),
                                    ('HC', 'HCs')]}
    assert seen == expected
