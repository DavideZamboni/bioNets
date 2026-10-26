"""Unit tests for the network-matrix helpers in mypsns/utils/matrix_utilities.py:
binarize_by_threshold, graph_to_matrix and load_combined_per_type."""
import networkx as nx
import numpy as np
import pytest

from mypsns.utils.matrix_utilities import (
    binarize_by_threshold, graph_to_matrix, load_combined_per_type,
    matrix_to_pairs, write_matrix,
)

# ═════════════════════════════════════════════════════════════════════════════
# binarize_by_threshold
# ═════════════════════════════════════════════════════════════════════════════

PERS = np.array([
    [50.0, 20.0,  0.0, 35.0],
    [20.0,  0.0, 10.0,  0.0],
    [ 0.0, 10.0,  5.0, 19.9],
    [35.0,  0.0, 19.9,  0.0],
])


def test_binarize_keeps_values_at_or_above_threshold():
    out = binarize_by_threshold(PERS, 20.0)
    expected = np.zeros_like(PERS)
    for i, j in [(0, 1), (0, 3)]:                  # 20.0 (== threshold) and 35.0
        expected[i, j] = expected[j, i] = 1.0
    np.testing.assert_array_equal(out, expected)


def test_binarize_zeroes_the_diagonal():
    out = binarize_by_threshold(PERS, 1.0)
    assert not np.any(np.diag(out))


def test_binarize_threshold_zero_keeps_only_positive_values():
    # A value of exactly 0 is a non-contact, not an edge, even at threshold 0.
    out = binarize_by_threshold(PERS, 0.0)
    off = ~np.eye(4, dtype=bool)
    np.testing.assert_array_equal(out[off], (PERS[off] > 0).astype(float))


def test_binarize_negative_threshold_still_needs_a_positive_value():
    mat = np.array([[0.0, -3.0], [-3.0, 0.0]])
    assert not np.any(binarize_by_threshold(mat, -5.0))


def test_binarize_returns_float_and_leaves_input_untouched():
    before = PERS.copy()
    out = binarize_by_threshold(PERS.astype(int), 20)
    assert out.dtype == float
    np.testing.assert_array_equal(PERS, before)


def test_threshold_semantics_weighted_strict_unweighted_inclusive():
    # The documented asymmetry: a pair exactly at the threshold is an edge of
    # the unweighted graph (>=) but not of the weighted one (>).
    labels = [0, 1, 2, 3]
    weighted = matrix_to_pairs(PERS, labels, value_name='weight', min_value=20.0)
    weighted_pairs = set(zip(weighted['res_1'], weighted['res_2']))
    binary = binarize_by_threshold(PERS, 20.0)
    assert (0, 1) not in weighted_pairs
    assert binary[0, 1] == 1.0
    assert (0, 3) in weighted_pairs and binary[0, 3] == 1.0

# ═════════════════════════════════════════════════════════════════════════════
# graph_to_matrix
# ═════════════════════════════════════════════════════════════════════════════

def _graph(edges, nodes=(10, 20, 30)):
    G = nx.Graph()
    G.add_nodes_from(nodes)
    for u, v, w in edges:
        G.add_edge(u, v, weight=w)
    return G


def test_graph_to_matrix_follows_label_order_and_is_symmetric():
    G = _graph([(10, 20, 2.5), (30, 20, 1.0)])
    mat = graph_to_matrix(G, labels=[30, 10, 20])
    expected = np.array([
        [0.0, 0.0, 1.0],    # 30
        [0.0, 0.0, 2.5],    # 10
        [1.0, 2.5, 0.0],    # 20
    ])
    np.testing.assert_array_equal(mat, expected)
    np.testing.assert_array_equal(mat, mat.T)


def test_graph_to_matrix_skips_unknown_endpoints_and_self_loops():
    G = _graph([(10, 20, 2.0), (10, 99, 7.0), (30, 30, 4.0)], nodes=(10, 20, 30, 99))
    mat = graph_to_matrix(G, labels=[10, 20, 30])
    assert np.count_nonzero(mat) == 2              # only 10–20, both triangles
    assert mat[0, 1] == mat[1, 0] == 2.0
    assert not np.any(np.diag(mat))


def test_graph_to_matrix_custom_weight_column():
    G = nx.Graph()
    G.add_edge(1, 2, weight=9.0, energy=3.0)
    mat = graph_to_matrix(G, labels=[1, 2], weight_col='energy')
    assert mat[0, 1] == mat[1, 0] == 3.0


def test_graph_to_matrix_empty_graph_gives_zeros():
    mat = graph_to_matrix(nx.Graph(), labels=[1, 2, 3])
    np.testing.assert_array_equal(mat, np.zeros((3, 3)))


def test_graph_to_matrix_round_trips_matrix_to_pairs():
    # What Network.write_matrices relies on: matrix → pairs → graph → matrix
    # reproduces the thresholded matrix.
    labels = [4, 7, 9, 12]
    sym = np.triu(PERS, 1)
    sym = sym + sym.T
    pairs = matrix_to_pairs(sym, labels, value_name='weight', min_value=15.0)
    G = nx.Graph()
    G.add_nodes_from(labels)
    for _, row in pairs.iterrows():
        G.add_edge(int(row['res_1']), int(row['res_2']), weight=row['weight'])
    expected = np.where(sym > 15.0, sym, 0.0)
    np.testing.assert_array_equal(graph_to_matrix(G, labels), expected)

# ═════════════════════════════════════════════════════════════════════════════
# load_combined_per_type
# ═════════════════════════════════════════════════════════════════════════════

ALL_TARGETS = ['Salt_bridges', 'HBonds', 'HCs']


def _sym(upper):
    """3×3 symmetric matrix from its upper-triangle entries (01, 02, 12)."""
    a01, a02, a12 = upper
    return np.array([[0.0, a01, a02], [a01, 0.0, a12], [a02, a12, 0.0]])


@pytest.fixture
def adjacency(tmp_path):
    """A minimal Adjacency/<inttype>/ folder in the format Contacts writes."""
    (tmp_path / 'adjacency_labels.dat').write_text(
        '5\tA:ARG:5\n6\tA:GLU:6\n8\tA:LEU:8\n')
    files = {
        # persistency (%)
        'SB_adjacency.dat':      _sym((40.0,  0.0,  0.0)),
        'Hbonds_adjacency.dat':  _sym((60.0, 10.0,  0.0)),
        'HC_adjacency.dat':      _sym(( 0.0, 30.0, 25.0)),
        # signed energies (kJ/mol); positive = repulsive
        'SB_total_energy.dat':     _sym((-8.0,  0.0,  0.0)),
        'Hbonds_total_energy.dat': _sym((-2.0,  3.0,  0.0)),
        'HC_total_energy.dat':     _sym(( 0.0, -1.0, -0.5)),
    }
    for name, mat in files.items():
        write_matrix(mat, tmp_path / name)
    return tmp_path


def test_persistency_defaults_to_elementwise_max(adjacency):
    mat, labels = load_combined_per_type(adjacency, ALL_TARGETS, 'persistency')
    assert labels == [5, 6, 8]
    np.testing.assert_allclose(mat, _sym((60.0, 30.0, 25.0)))


def test_persistency_sum_override(adjacency):
    mat, _ = load_combined_per_type(adjacency, ALL_TARGETS, 'persistency', combine='sum')
    np.testing.assert_allclose(mat, _sym((100.0, 40.0, 25.0)))


def test_energy_defaults_to_sum_of_attractive_weights(adjacency):
    # max(0, -E) per type, then summed: repulsive +3.0 (Hbonds 0–2) is dropped.
    mat, _ = load_combined_per_type(adjacency, ALL_TARGETS, 'energy')
    np.testing.assert_allclose(mat, _sym((10.0, 1.0, 0.5)))


def test_energy_use_abs_counts_repulsive_terms(adjacency):
    mat, _ = load_combined_per_type(adjacency, ALL_TARGETS, 'energy', use_abs=True)
    np.testing.assert_allclose(mat, _sym((10.0, 4.0, 0.5)))


def test_energy_max_override(adjacency):
    mat, _ = load_combined_per_type(adjacency, ALL_TARGETS, 'energy', combine='max')
    np.testing.assert_allclose(mat, _sym((8.0, 1.0, 0.5)))


def test_only_requested_targets_are_combined(adjacency):
    mat, _ = load_combined_per_type(adjacency, ['Salt_bridges', 'HCs'], 'persistency')
    np.testing.assert_allclose(mat, _sym((40.0, 30.0, 25.0)))


def test_missing_target_file_is_skipped(adjacency):
    (adjacency / 'HC_adjacency.dat').unlink()
    mat, _ = load_combined_per_type(adjacency, ALL_TARGETS, 'persistency')
    np.testing.assert_allclose(mat, _sym((60.0, 10.0, 0.0)))


def test_no_target_file_gives_none_with_labels(adjacency):
    for name in ('SB_adjacency.dat', 'Hbonds_adjacency.dat', 'HC_adjacency.dat'):
        (adjacency / name).unlink()
    mat, labels = load_combined_per_type(adjacency, ALL_TARGETS, 'persistency')
    assert mat is None
    assert labels == [5, 6, 8]


def test_missing_labels_gives_none_and_empty_labels(adjacency):
    (adjacency / 'adjacency_labels.dat').unlink()
    assert load_combined_per_type(adjacency, ALL_TARGETS, 'persistency') == (None, [])


def test_invalid_combine_raises(adjacency):
    with pytest.raises(ValueError, match='combine'):
        load_combined_per_type(adjacency, ALL_TARGETS, 'persistency', combine='mean')


def test_inputs_on_disk_are_not_modified(adjacency):
    before = np.loadtxt(adjacency / 'Hbonds_adjacency.dat')
    load_combined_per_type(adjacency, ALL_TARGETS, 'persistency')
    load_combined_per_type(adjacency, ALL_TARGETS, 'energy')
    np.testing.assert_array_equal(np.loadtxt(adjacency / 'Hbonds_adjacency.dat'), before)
