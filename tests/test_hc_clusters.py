"""q_interactors hydrophobic clusters (HC_cluster) are deterministic, avoid
size violations where the bonds allow it (exact search, greedy fallback), are
stored by the interactors stage, and are only read back by contacts.

The split used to follow the iteration order of Python sets of atom names,
which changes with PYTHONHASHSEED: over 12 seeds this PHE got 7 different
partitions.
"""
import inspect
import io
import json
import os
import subprocess
import sys

import pandas as pd

import mypsns.services.contacts as contacts_service
import mypsns.utils.contacts_utilities as contacts_utilities
from mypsns.services.interactors import Interactors
from mypsns.utils import interactors_utilities
from mypsns.utils.contacts_utilities import get_selection
from mypsns.utils.interactors_utilities import (
    get_charge_based_interactors, get_constrained_conncomps,
)

from conftest import PHE_VARIANT, variant_frame

HEAVY_RING = ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']

CLUSTERS_SCRIPT = """
import json, sys
import pandas as pd
from mypsns.utils.interactors_utilities import get_charge_based_interactors
df = pd.read_csv(sys.argv[1], keep_default_na=False)
df['Resindex'] = 0
df['Entity'] = 'M'
df = get_charge_based_interactors(df)
print(json.dumps(dict(zip(df['atom_name'], df['HC_cluster']))))
"""


def clusters(df):
    '''{label: set of atom names} of a typed variant.'''
    hc = df[df['HC_cluster'] != '']
    return {label: set(g['atom_name']) for label, g in hc.groupby('HC_cluster')}


def test_same_clusters_under_any_hash_seed(tmp_path):
    csv = tmp_path / 'PHE.csv'
    csv.write_text(PHE_VARIANT)
    runs = set()
    for seed in range(8):
        res = subprocess.run([sys.executable, '-c', CLUSTERS_SCRIPT, str(csv)],
                             env={**os.environ, 'PYTHONHASHSEED': str(seed)},
                             capture_output=True, text=True, check=True)
        runs.add(res.stdout.strip().splitlines()[-1])
    assert len(runs) == 1, runs
    labels = json.loads(runs.pop())
    by_label = {}
    for atom, label in labels.items():
        if label:
            by_label.setdefault(label, set()).add(atom)
    # the clusters partition the ring + CB (with their H) and are numbered in
    # sorted order of their atoms
    heavy = [a for a in HEAVY_RING if labels[a]]
    assert sorted(heavy) == sorted(HEAVY_RING)
    assert sorted(by_label) == [f'HC{i}' for i in range(1, len(by_label) + 1)]
    firsts = [min(a for a in by_label[l] if a in HEAVY_RING) for l in sorted(by_label)]
    assert firsts == sorted(firsts)


def test_phe_split_has_no_size_violation(phe_variant):
    # The greedy split gave CB alone + CD1 CD2 CG + CE1 CE2 CZ (one cluster
    # under min_num_heavy_ats=2); the exact search finds 3 + 2 + 2.
    got = clusters(get_charge_based_interactors(phe_variant))
    assert got == {'HC1': {'CB', 'CD1', 'CG', 'HB1', 'HB2', 'HD1'},
                   'HC2': {'CD2', 'CE2', 'HD2', 'HE2'},
                   'HC3': {'CE1', 'CZ', 'HE1', 'HZ'}}


def test_chain_split_prefers_fewest_violations_then_fewest_clusters():
    chain = {f'C{i}': [f'C{j}' for j in (i - 1, i + 1) if 1 <= j <= 7] for i in range(1, 8)}
    parts = get_constrained_conncomps(chain, min_num_heavy_ats=2, max_num_heavy_ats=3)
    assert sorted(len(p) for p in parts) == [2, 2, 3]      # not 3 + 3 + 1, not 2 + 2 + 2 + 1
    long_chain = {f'C{i:02d}': [f'C{j:02d}' for j in (i - 1, i + 1) if 1 <= j <= 16]
                  for i in range(1, 17)}
    parts = get_constrained_conncomps(long_chain, min_num_heavy_ats=2, max_num_heavy_ats=3)
    assert sorted(len(p) for p in parts) == [2, 2, 3, 3, 3, 3]   # 6 clusters, the fewest possible


def test_falls_back_to_greedy_past_the_budget(monkeypatch, caplog):
    ring = {'CB': ['CG'], 'CG': ['CB', 'CD1', 'CD2'], 'CD1': ['CG', 'CE1'],
            'CE1': ['CD1', 'CZ'], 'CZ': ['CE1', 'CE2'], 'CE2': ['CZ', 'CD2'],
            'CD2': ['CE2', 'CG']}
    monkeypatch.setattr(interactors_utilities, 'MAX_PARTITION_WORK', 1)
    with caplog.at_level('WARNING'):
        parts = get_constrained_conncomps(ring, min_num_heavy_ats=2, max_num_heavy_ats=3)
    assert parts == [['CB'], ['CD1', 'CD2', 'CG'], ['CE1', 'CE2', 'CZ']]   # the greedy split
    assert 'exact search budget' in caplog.text


def test_split_ignores_input_order():
    ring = {'CB': ['CG'], 'CG': ['CB', 'CD1', 'CD2'], 'CD1': ['CG', 'CE1'],
            'CE1': ['CD1', 'CZ'], 'CZ': ['CE1', 'CE2'], 'CE2': ['CZ', 'CD2'],
            'CD2': ['CE2', 'CG']}
    reordered = {k: list(reversed(ring[k])) for k in reversed(list(ring))}
    a = get_constrained_conncomps(ring, min_num_heavy_ats=2, max_num_heavy_ats=3)
    b = get_constrained_conncomps(reordered, min_num_heavy_ats=2, max_num_heavy_ats=3)
    assert a == b
    assert a == sorted(sorted(c) for c in a)
    assert sorted(x for c in a for x in c) == sorted(ring)


def test_interactors_stores_the_clusters(tmp_path):
    found = tmp_path / 'Interactors' / 'found_topologies'
    found.mkdir(parents=True)
    pd.read_csv(io.StringIO(PHE_VARIANT), keep_default_na=False).to_csv(
        found / 'PHE.csv', index=False)
    stage = object.__new__(Interactors)      # add_inttypes only needs out/start
    stage.out, stage.start = tmp_path / 'Interactors', 0.0
    stage.add_inttypes()
    stored = pd.read_csv(found / 'PHE.csv', keep_default_na=False)
    assert clusters(stored) == clusters(get_charge_based_interactors(variant_frame()))
    assert {'q_interactors', 'named_interactors', 'HC_cluster'} <= set(stored.columns)


def test_contacts_uses_the_stored_labels():
    # Labels that the clustering would never produce (CB with CZ, not bonded):
    # contacts must take them as they are.
    chain = variant_frame()
    chain['Resindex'] = 7
    chain['q_interactors'] = ''
    chain['HC_cluster'] = pd.Series([None] * len(chain), dtype=object)   # NaN when read back
    stored = {'HC1': ['CB', 'CZ', 'HZ'], 'HC2': ['CG', 'CD1', 'CD2', 'CE1', 'CE2']}
    for label, atoms in stored.items():
        rows = chain['atom_name'].isin(atoms)
        chain.loc[rows, 'q_interactors'] = ' HC'
        chain.loc[rows, 'HC_cluster'] = label
    sel = get_selection(chain, 'q_interactors', 'HC', 'A')
    assert set(sel) == {'HC1_7', 'HC2_7'}
    for label, atoms in stored.items():
        sel_str, n = sel[f'{label}_7']
        assert n == len(atoms)
        assert sel_str.endswith(f'name {" ".join(sorted(atoms))}))')


def test_contacts_never_recomputes_the_clusters():
    for module in (contacts_utilities, contacts_service):
        src = inspect.getsource(module)
        for name in ('interactors_utilities', 'add_hc_atoms', 'get_constrained_conncomps',
                     'get_charge_based_interactors'):
            assert name not in src, f'{module.__name__} refers to {name}'
