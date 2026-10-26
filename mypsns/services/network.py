"""
network.py — Network service for myPSNs (energy + time-resolved support).

Workflow
--------
1. :meth:`load_static`              – build one residue-interaction graph per
   (interactors_type × edge_weight_source) by reading the per-type matrices
   that ``Contacts`` writes (e.g. ``SB_total_energy.dat`` + ``Hbonds_total_energy.dat``
   + ``HC_total_energy.dat``), transformed (attractive-only or absolute) and
   summed across the requested targets.

2. :meth:`run_analysis`             – centrality + community detection per graph.

3. :meth:`write_pymol`              – one .pml per
   (interactors_type × source × node_metric × community_method).

4. :meth:`load_time_resolved_energy` (optional) – read the per-frame energy
   CSVs that the patched ``get_interaction_energy`` / ``get_hc_energy``
   produce, aggregate into windows of ``window_size`` frames, and build one
   graph per (interactors_type × window_start).

5. :meth:`track_communities_over_time` – run community detection on every
   windowed graph and dump ``community_evolution_*.csv`` with one row per
   (window_start, residue) — pivot externally for the "community vs time"
   visualisation.

6. :meth:`write_outputs` / :meth:`get_shortest_paths` / :meth:`apply_filter`
   — utilities that operate on a chosen graph key.

Graph identifiers
-----------------
Static graphs use a 2-tuple key ``(interactors_type, source)``.
Time-resolved graphs use a 3-tuple ``(interactors_type, source, window_start)``.
Currently `source` is always ``'energy'`` for time-resolved (only the per-frame
energy CSVs carry per-frame information).
"""
from __future__ import annotations

import logging
import timeit
from os import makedirs
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import networkx as nx
import numpy as np
import pandas as pd

from mypsns.utils.matrix_utilities import (
    build_chain_map, build_label_map,
    binarize_by_threshold, energy_transform, graph_to_matrix,
    load_combined_per_type, load_per_frame_energy, load_per_type_matrix,
    matrix_to_pairs, pairs_to_matrix, strength_to_percentage,
    write_labels, write_matrix, window_per_frame,
)
from mypsns.utils.network_utilities import (
    build_resindex_to_pdb_map,
    compute_centrality, detect_communities,
    filter_by_degree, filter_by_weight,
    write_edge_csv, write_graphml, write_node_csv,
    write_pymol_nodes, write_pymol_communities,
)
from mypsns.utils.general_utilities import get_endtime
from mypsns.utils.contacts_utilities import assemble_full_interactors

pd.options.mode.chained_assignment = None
log = logging.getLogger(__name__)

_NA_VALUES: List[str] = [
    "", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN",
    "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A",
    "NULL", "NaN", "n/a", "nan", "null",
]


# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────

def _threshold_for_source(
    source:          str,
    min_persistency: float,
    min_energy:      float,
) -> float:
    """Pick the per-source edge threshold.

    persistency graphs are on a 0..100 % scale → use min_persistency.
    energy graphs are on a kJ/mol (post-transform) scale → use min_energy.
    Any other source falls back to 0.0 (no filtering).
    """
    if source == 'persistency':
        return min_persistency
    if source == 'energy':
        return min_energy
    return 0.0


def _resolve_weightings(weighting: str) -> List[str]:
    """Expand a weighting choice into the list of edge modes to build.

    'weighted'   → weighted edges only (edge weight = the value)
    'unweighted' → binary edges only   (edge weight = 1 above threshold)
    'both'       → build both (default)
    """
    w = (weighting or 'both').lower()
    if w == 'both':
        return ['weighted', 'unweighted']
    if w in ('weighted', 'unweighted'):
        return [w]
    raise ValueError(
        f"weighting must be 'weighted', 'unweighted', or 'both' — got {weighting!r}")


def _source_label(source: str, mode: str) -> str:
    """Graph/output label for a (source, mode) pair.

    Weighted keeps the bare source name (backward-compatible file names);
    unweighted gets an explicit suffix. Both variants therefore coexist under
    distinct keys/files when weighting='both'.
    """
    return source if mode == 'weighted' else f'{source}_unweighted'


def _build_graph_from_pairs(
    pairs:    pd.DataFrame,
    int_df:   pd.DataFrame,
    weight_col: str = 'weight',
) -> nx.Graph:
    """Build a nx.Graph from a (res_1, res_2, weight) table with node attrs
    pulled from `int_df` (label, chain, resname, resindex)."""
    G = nx.Graph()
    label_map = build_label_map(int_df)
    chain_map = build_chain_map(int_df)
    uniq = int_df.drop_duplicates('Resindex')
    resname_map = dict(zip(uniq['Resindex'].astype(int), uniq['resname']))
    for resix, label in label_map.items():
        G.add_node(
            int(resix),
            label    = label,
            chain    = chain_map.get(int(resix), ''),
            resname  = resname_map.get(int(resix), ''),
            resindex = int(resix),
        )
    for _, row in pairs.iterrows():
        G.add_edge(
            int(row['res_1']), int(row['res_2']),
            weight=float(row[weight_col]),
        )
    return G


# ═════════════════════════════════════════════════════════════════════════════
# Network service
# ═════════════════════════════════════════════════════════════════════════════

class Network:
    """
    Residue-interaction network analysis from Contacts output.

    Parameters
    ----------
    topology_pdb : str or Path
        Used to build the Resindex → (chain, PDB resnum) map for correct PyMOL
        selectors, and as the default reference structure for the PyMOL sessions.
    general_outpath : str or Path
        Run root. Must contain ``Interactors/`` and ``Adjacency/``.
    n_frames : int, optional
        Total trajectory frames analysed during Contacts. Informational only
        (logged); analysis derives frame counts directly from the data on disk.
    na_values : list[str], optional
        Pandas NA strings.
    """

    def __init__(
        self,
        topology_pdb:    Union[str, Path],
        general_outpath: Union[str, Path],
        n_frames:        Optional[int] = None,
        na_values:       Optional[List[str]] = None,
    ) -> None:
        general_outpath = Path(general_outpath)
        self.topology_pdb     = Path(topology_pdb)
        self.adjacency_root   = general_outpath / 'Adjacency'
        self.interactors_root = general_outpath / 'Interactors'
        self.out              = general_outpath / 'Network'
        self.nas              = na_values if na_values is not None else list(_NA_VALUES)
        self.n_frames         = int(n_frames) if n_frames else None
        self.start            = timeit.default_timer()

        # static
        self.graphs:         Dict[Tuple[str, str], nx.Graph]        = {}
        self.centrality_dfs: Dict[Tuple[str, str], pd.DataFrame]    = {}
        self.community_dfs:  Dict[Tuple[str, str],
                                  Dict[str, pd.DataFrame]]          = {}
        # pre-threshold matrices captured during load_static, so write_matrices()
        # can dump the combined and per-target contributions alongside the final
        # (thresholded) network matrix. {(inttype, source): {'combined': ndarray,
        # 'per_target': {target: ndarray, ...}}}.
        self._source_matrices: Dict[Tuple[str, str], Dict[str, object]] = {}
        # time-resolved
        self.windowed_graphs: Dict[Tuple[str, str, int], nx.Graph]  = {}

        makedirs(self.out, exist_ok=True)
        log.info('Network run | adj=%s | interactors=%s | out=%s | n_frames=%s',
                 self.adjacency_root, self.interactors_root, self.out,
                 self.n_frames if self.n_frames is not None else 'all')

        # ── Reference structure for PyMOL selectors ──────────────────────────
        # The graph's residue indices / chains were built against the CORRECTED
        # topology Contacts wrote (Adjacency/<name>_correct.pdb): Resindex stamped
        # as resid, chains re-stamped to match. Selectors (chain X and resi N)
        # only resolve against that structure, so prefer it over whatever the
        # user passed in. Fall back to the supplied topology if none is found.
        self.reference_pdb = self._find_corrected_pdb()
        if self.reference_pdb is not None:
            log.info('Using corrected topology for PyMOL: %s', self.reference_pdb)
        else:
            self.reference_pdb = self.topology_pdb
            log.warning('No *_correct.pdb found in %s — falling back to the '
                        'supplied topology (%s). PyMOL edges/nodes may not '
                        'resolve if its chains/resids differ from the network.',
                        self.adjacency_root, self.topology_pdb)

        # ── Interactors table ────────────────────────────────────────────────
        # Reconstructed in memory from the three normalised files Interactors
        # writes (found_topologies/<variant>.csv, residue_index.csv,
        # atom_index_map.csv) — exactly as Contacts does. There is no monolithic
        # whole_topology_interactors.csv on disk (kept normalised to stay small).
        residue_index_path  = self.interactors_root / 'residue_index.csv'
        atom_index_map_path = self.interactors_root / 'atom_index_map.csv'
        found_topologies_dir = self.interactors_root / 'found_topologies'
        missing = [p for p in (residue_index_path, atom_index_map_path,
                               found_topologies_dir) if not p.exists()]
        if missing:
            raise FileNotFoundError(
                f'Interactors output incomplete in {self.interactors_root.absolute()} '
                f'(missing: {[p.name for p in missing]}). '
                'Run Topology_recognition / Interactors first.'
            )
        self.int_top = assemble_full_interactors(
            found_topologies_dir = found_topologies_dir,
            residue_index_path   = residue_index_path,
            atom_index_map_path  = atom_index_map_path,
            na_values            = self.nas,
        )
        log.info('interactors table: %d atoms, %d residues',
                 self.int_top.shape[0], self.int_top['Resindex'].nunique())

        # ── Resindex → (chain, PDB resid) for PyMOL ──────────────────────────
        try:
            self.resindex_to_pdb = build_resindex_to_pdb_map(self.reference_pdb)
            log.debug('Resindex→(chain,resid): %d residues', len(self.resindex_to_pdb))
        except Exception as e:
            log.warning('build_resindex_to_pdb_map failed (%s) — PyMOL selectors '
                        'will fall back to Resindex.', e)
            self.resindex_to_pdb = {}

    def _find_corrected_pdb(self) -> Optional[Path]:
        """Locate the corrected topology Contacts wrote, ``<name>_correct.pdb``
        under ``Adjacency/``.

        Preference order:
          1. ``Adjacency/<supplied_topology_stem>_correct.pdb`` (the exact match
             for the structure this pipeline was run on);
          2. any single ``Adjacency/*_correct.pdb`` if exactly one exists;
          3. None (caller falls back to the supplied topology).
        """
        if not self.adjacency_root.exists():
            return None
        # 1. exact stem match
        stem = self.topology_pdb.stem
        exact = self.adjacency_root / f'{stem}_correct.pdb'
        if exact.exists():
            return exact
        # 2. unique *_correct.pdb
        candidates = sorted(self.adjacency_root.glob('*_correct.pdb'))
        if len(candidates) == 1:
            return candidates[0]
        if len(candidates) > 1:
            log.warning('Multiple *_correct.pdb in %s (%s) — cannot pick '
                        'unambiguously; using supplied topology.',
                        self.adjacency_root, [c.name for c in candidates])
        return None

    # ═════════════════════════════════════════════════════════════════════════
    # Static network from contacts matrices
    # ═════════════════════════════════════════════════════════════════════════

    def load_static(
        self,
        interactors_type:    List[str],
        targets:             List[str],
        edge_weight_sources: List[str],
        energy_use_abs:      bool  = False,
        min_persistency:     float = 0.0,
        min_energy:          float = 0.0,
        weighting:           str   = 'both',
        frame_start:         Optional[int] = None,
        frame_stop:          Optional[int] = None,
    ) -> None:
        """
        For each (inttype × source × weighting mode), build a graph stored at
        ``self.graphs[(inttype, label)]`` where ``label`` is the source name for
        weighted graphs and ``<source>_unweighted`` for unweighted ones.

        Two paths
        ---------
        * Fast (no frame range): reads the per-type matrices Contacts wrote
          (e.g. ``SB_total_energy.dat``) and combines across ``targets``
          (persistency by MAX/union, energy by SUM).
        * Per-frame (frame_start and/or frame_stop set): re-aggregates from
          the per-frame energy CSVs restricted to ``[frame_start, frame_stop)``.

        Weighting
        ---------
        Both persistency and energy can be built weighted, unweighted, or both:
          * weighted   — edge weight is the combined value (persistency percent,
            or energy as "% of the strongest edge"); edges with value > the
            threshold are kept.
          * unweighted — edge weight is 1 wherever the combined value ≥ the
            threshold (and > 0); this is the binary macro-IIN analogue for
            persistency, and the same idea applied to energy strength.

        Parameters
        ----------
        interactors_type : list of {'q_interactors', 'named_interactors'}
            Sub-folders of Adjacency/ to consume.
        targets : list of {'Salt_bridges', 'HBonds', 'HCs'}
            Contact types contributing to each graph.
        edge_weight_sources : list of {'persistency', 'energy'}
            One graph (per weighting mode) is built per source.
        energy_use_abs : bool
            False (default) → attractive-only: weight = max(0, -E).
            True            → absolute:        weight = |E|.
        min_persistency : float
            Persistency-graph edge threshold (percent, 0..100). Weighted keeps
            values > this; unweighted keeps an edge (weight 1) where value ≥ this.
        min_energy : float
            Energy-graph edge threshold (percent of the strongest edge, 0..100).
            Energy strengths are normalised to 0..100 before thresholding.
        weighting : {'both', 'weighted', 'unweighted'}
            Which edge weighting(s) to build for every source. Default 'both'.
        frame_start, frame_stop : int, optional
            If either is set, re-aggregate from per-frame CSVs over
            ``[frame_start, frame_stop)``.
        """
        modes = _resolve_weightings(weighting)

        if frame_start is not None or frame_stop is not None:
            return self._load_static_from_per_frame(
                interactors_type, targets, edge_weight_sources,
                energy_use_abs, min_persistency, min_energy, modes,
                frame_start, frame_stop,
            )

        start = timeit.default_timer()

        for inttype in interactors_type:
            folder = self.adjacency_root / inttype
            if not folder.exists():
                log.warning('load_static: %s not found — skipped', folder)
                continue
            for source in edge_weight_sources:
                # combine per type: persistency by MAX (union → macro-IIN when
                # binarised), energy by SUM of attractive weights. That rule is
                # load_combined_per_type's default for each source.
                mat, labels = load_combined_per_type(
                    folder, targets=targets, source=source,
                    use_abs=energy_use_abs,
                )
                if mat is None:
                    log.warning('load_static: no matrices for %s/%s — skipped',
                                inttype, source)
                    continue
                # energy is normalised to 0..100 "% of strongest edge" so its
                # threshold and edge weights are on the same scale as the energy
                # elbow/degree plots.
                if source == 'energy':
                    mat = strength_to_percentage(mat)
                # capture per-target contributions (energy → attractive-only/abs
                # weight; persistency → raw %) and the combined pre-threshold
                # matrix, shared by both weighting modes, for write_matrices().
                per_target: Dict[str, np.ndarray] = {}
                for target in targets:
                    pm = load_per_type_matrix(folder, target, source)
                    if pm is None:
                        continue
                    per_target[target] = (energy_transform(pm, use_abs=energy_use_abs)
                                          if source == 'energy' else pm)
                threshold = _threshold_for_source(source, min_persistency, min_energy)

                for mode in modes:
                    label = _source_label(source, mode)
                    self._source_matrices[(inttype, label)] = {
                        'combined':   mat.copy(),
                        'per_target': per_target,
                    }
                    if mode == 'unweighted':
                        binary = binarize_by_threshold(mat, threshold)
                        pairs  = matrix_to_pairs(binary, labels, value_name='weight',
                                                 min_value=0.0)
                        G = _build_graph_from_pairs(pairs, self.int_top, weight_col='weight')
                        self.graphs[(inttype, label)] = G
                        log.info('  %s/%s graph (unweighted, value>=%g): '
                                 '%d nodes, %d edges', inttype, label, threshold,
                                 G.number_of_nodes(), G.number_of_edges())
                    else:
                        pairs = matrix_to_pairs(mat, labels, value_name='weight',
                                                min_value=threshold)
                        G = _build_graph_from_pairs(pairs, self.int_top, weight_col='weight')
                        self.graphs[(inttype, label)] = G
                        log.info('  %s/%s graph (weighted, value>%g): '
                                 '%d nodes, %d edges', inttype, label, threshold,
                                 G.number_of_nodes(), G.number_of_edges())

        proc_time, total_time = get_endtime(start, self.start)
        log.info('load_static complete | took %s (h:m:s), total %s (h:m:s)',
                 proc_time, total_time)

    # ─────────────────────────────────────────────────────────────────────────

    def _load_static_from_per_frame(
        self,
        interactors_type:    List[str],
        targets:             List[str],
        edge_weight_sources: List[str],
        energy_use_abs:      bool,
        min_persistency:     float,
        min_energy:          float,
        modes:               List[str],
        frame_start:         Optional[int],
        frame_stop:          Optional[int],
    ) -> None:
        """Frame-range path for load_static — uses per-frame energy CSVs.

        Persistency is scaled to percent (0..100) so the value and the
        ``min_persistency`` threshold mean the same thing as in the
        matrix-based path. Energy is mean-per-frame weight (kJ/mol).
        """
        start = timeit.default_timer()

        for inttype in interactors_type:
            folder = self.adjacency_root / inttype
            if not folder.exists():
                log.warning('load_static (per-frame): %s missing — skipped', folder)
                continue
            df = load_per_frame_energy(folder, targets=targets, use_abs=energy_use_abs)
            if df.empty:
                log.warning('load_static (per-frame): no per-frame energy CSVs '
                            'under %s — skipped (re-run contacts with the '
                            'streaming-energy patch?)', folder)
                continue

            if frame_start is not None:
                df = df[df['frame'] >= frame_start]
            if frame_stop is not None:
                df = df[df['frame'] <  frame_stop]
            if df.empty:
                log.warning('load_static (per-frame): frame range [%s, %s) '
                            'excluded all data under %s — skipped',
                            frame_start, frame_stop, folder)
                continue

            # canonical order
            swap = df['res_1'] > df['res_2']
            df.loc[swap, ['res_1', 'res_2']] = df.loc[swap, ['res_2', 'res_1']].values

            # effective frame count for normalisation
            eff_start = frame_start if frame_start is not None else int(df['frame'].min())
            eff_stop  = frame_stop  if frame_stop  is not None else int(df['frame'].max()) + 1
            n_frames_eff = max(1, eff_stop - eff_start)
            log.info('  %s: per-frame aggregation, frames [%d, %d)  (n=%d)',
                     inttype, eff_start, eff_stop, n_frames_eff)

            for source in edge_weight_sources:
                if source == 'energy':
                    agg = df.groupby(['res_1', 'res_2'], as_index=False)['weight'].sum()
                    agg['weight'] = agg['weight'] / n_frames_eff
                    # normalise to 0..100 "% of strongest edge" (matches matrix path)
                    wmax = agg['weight'].max() if len(agg) else 0.0
                    if wmax > 0:
                        agg['weight'] = 100.0 * agg['weight'] / wmax
                    else:
                        agg['weight'] = 0.0
                elif source == 'persistency':
                    pers = (df.groupby(['target', 'res_1', 'res_2'])['frame']
                              .nunique().reset_index(name='hit'))
                    # percent, to match the matrix-based persistency convention
                    pers['weight'] = pers['hit'] / n_frames_eff * 100.0
                    # combine targets by MAX (union), matching the matrix path
                    agg = pers.groupby(['res_1', 'res_2'], as_index=False)['weight'].max()
                else:
                    log.warning('load_static (per-frame): unknown source %r — skipped',
                                source)
                    continue

                # capture combined pre-threshold matrix for write_matrices()
                # (per-target contributions aren't separable in the per-frame
                # path — those files are only written from the matrix path).
                labels_all = sorted(self.int_top['Resindex'].unique().astype(int).tolist())
                combined_mat = pairs_to_matrix(agg, labels_all, value_col='weight')
                threshold = _threshold_for_source(source, min_persistency, min_energy)

                for mode in modes:
                    label = _source_label(source, mode)
                    self._source_matrices[(inttype, label)] = {
                        'combined':   combined_mat.copy(),
                        'per_target': {},
                    }
                    if mode == 'unweighted':
                        sub = agg[(agg['weight'] >= threshold) &
                                  (agg['weight'] > 0.0)].copy()
                        if sub.empty:
                            log.warning('  %s/%s: no edges with value>=%g',
                                        inttype, label, threshold)
                            continue
                        sub['weight'] = 1.0
                        G = _build_graph_from_pairs(sub, self.int_top, weight_col='weight')
                        self.graphs[(inttype, label)] = G
                        log.info('  %s/%s graph (per-frame, unweighted, value>=%g): '
                                 '%d nodes, %d edges', inttype, label, threshold,
                                 G.number_of_nodes(), G.number_of_edges())
                    else:
                        sub = agg[agg['weight'] > threshold]
                        if sub.empty:
                            log.warning('  %s/%s: no edges above value>%g',
                                        inttype, label, threshold)
                            continue
                        G = _build_graph_from_pairs(sub, self.int_top, weight_col='weight')
                        self.graphs[(inttype, label)] = G
                        log.info('  %s/%s graph (per-frame, weighted, value>%g): '
                                 '%d nodes, %d edges', inttype, label, threshold,
                                 G.number_of_nodes(), G.number_of_edges())

        proc_time, total_time = get_endtime(start, self.start)
        log.info('load_static (per-frame) complete | took %s (h:m:s), total %s (h:m:s)',
                 proc_time, total_time)

    # ─────────────────────────────────────────────────────────────────────────

    def write_matrices(
        self,
        subdir: str = 'matrices',
    ) -> List[Path]:
        """
        Serialise every static network's adjacency matrices to ``.dat`` plus a
        shared label list, in the exact format the Adjacency stage uses
        (``numpy.savetxt`` whitespace-separated square matrix; ``adjacency_labels.dat``
        as ``<resindex>\\t<Chain:resname:Resindex>`` per line).

        Layout under ``Network/<inttype>/<subdir>/``::

            adjacency_labels.dat                 (once per inttype)
            <label>_network_adjacency.dat        final graph (post threshold/degree filter)
            <label>_combined_adjacency.dat       combined-across-targets, pre-threshold
            <label>_<target>_adjacency.dat       per-target contribution (matrix path only)

        where ``<label>`` is ``<source>`` for a weighted graph and
        ``<source>_unweighted`` for an unweighted one (so both coexist when
        ``weighting='both'``).

        Scales (documented so the numbers are unambiguous):
          * persistency  — per-type files (``*_<target>``) carry the raw
            per-type persistency percent (0..100); ``*_combined`` is the per-type
            MAX (union) percent, pre-threshold. The weighted ``*_network`` keeps
            those percents; the unweighted ``*_network`` is **binary** (1 where
            combined persistency ≥ min_persistency) — the residue-network
            analogue of PyInteraph's macro-IIN.
          * energy       — ``*_combined`` and the weighted ``*_network`` are on
            the network's 0..100 "% of the strongest edge" scale; per-target
            ``*_<target>`` files are the pre-normalisation attractive-only /
            absolute weight (kJ/mol). The unweighted ``*_network`` is binary
            (1 where that strength ≥ min_energy).

        The ``*_network_adjacency.dat`` matrix is reconstructed straight from the
        built graph, so re-loading it reproduces exactly the thresholded edge set
        the rest of the pipeline analyses. Returns the list of written paths.
        """
        if not self.graphs:
            raise RuntimeError('write_matrices() called before load_static()')

        start = timeit.default_timer()
        written: List[Path] = []

        # group graph keys by interactors_type so the label file is shared
        by_inttype: Dict[str, List[Tuple[str, str]]] = {}
        for key in self.graphs:
            by_inttype.setdefault(key[0], []).append(key)

        for inttype, keys in by_inttype.items():
            mdir = self.out / inttype / subdir
            makedirs(mdir, exist_ok=True)

            # canonical label order — identical content/format to Adjacency's
            # adjacency_labels.dat; also the row/col order for every matrix here.
            labels = write_labels(self.int_top, mdir / 'adjacency_labels.dat')
            written.append(mdir / 'adjacency_labels.dat')

            for key in keys:
                _inttype, source = key
                G = self.graphs[key]

                # final (thresholded) network matrix — round-trips with the graph.
                # (+ 0.0 normalises IEEE negative zeros so the .dat has no
                #  confusing '-0.000000e+00' entries; values are unchanged.)
                final_mat = graph_to_matrix(G, labels, weight_col='weight') + 0.0
                fp = mdir / f'{source}_network_adjacency.dat'
                write_matrix(final_mat, fp, fmt='%.6e')
                written.append(fp)
                nz = int(np.count_nonzero(final_mat) // 2)
                if nz != G.number_of_edges():
                    log.warning('  %s/%s: matrix has %d pairs but graph has %d '
                                'edges (endpoints outside label set?)',
                                inttype, source, nz, G.number_of_edges())
                log.info('  %s/%s: %s (%d×%d, %d non-zero pairs)',
                         inttype, source, fp.name, len(labels), len(labels), nz)

                cap = self._source_matrices.get(key, {})

                # combined-across-targets, pre-threshold
                combined = cap.get('combined')
                if isinstance(combined, np.ndarray) and combined.size:
                    cp = mdir / f'{source}_combined_adjacency.dat'
                    write_matrix(combined + 0.0, cp, fmt='%.6e')
                    written.append(cp)

                # per-target contributions (matrix path only)
                for target, pm in cap.get('per_target', {}).items():
                    if not (isinstance(pm, np.ndarray) and pm.size):
                        continue
                    pp = mdir / f'{source}_{target}_adjacency.dat'
                    write_matrix(pm + 0.0, pp, fmt='%.6e')
                    written.append(pp)

        proc_time, total_time = get_endtime(start, self.start)
        log.info('write_matrices: %d files across %d interactors_type(s) | '
                 'took %s (h:m:s), total %s (h:m:s)',
                 len(written), len(by_inttype), proc_time, total_time)
        return written

    # ─────────────────────────────────────────────────────────────────────────

    def run_analysis(
        self,
        community_methods: Optional[List[str]] = None,
        community_seed:    int                  = 42,
        community_min_members: int              = 2,
    ) -> None:
        """Centrality + community detection, per static graph."""
        if not self.graphs:
            raise RuntimeError('run_analysis() called before load_static()')
        if community_methods is None:
            community_methods = ['louvain']

        start = timeit.default_timer()

        for key, G in self.graphs.items():
            inttype, source = key
            cdf_centrality = compute_centrality(G)
            self.centrality_dfs[key] = cdf_centrality
            self.community_dfs[key] = {}

            csv_dir = self.out / inttype / 'csv'
            makedirs(csv_dir, exist_ok=True)

            for method in community_methods:
                cdf, n_comm = detect_communities(
                    G, method=method, seed=community_seed,
                    min_members=community_min_members,
                )
                self.community_dfs[key][method] = cdf
                log.info('  %s/%s: %s → %d communities',
                         inttype, source, method, n_comm)
                nout = write_node_csv(
                    G, cdf_centrality, cdf, csv_dir,
                    filename=f'nodes_{source}_{method}.csv',
                )
                log.info('    nodes CSV: %s', nout)
            eout = write_edge_csv(
                G, csv_dir, filename=f'edges_{source}.csv',
            )
            log.info('    edges CSV: %s', eout)

        proc_time, total_time = get_endtime(start, self.start)
        log.info('run_analysis complete | took %s (h:m:s), total %s (h:m:s)',
                 proc_time, total_time)

    # ─────────────────────────────────────────────────────────────────────────

    def write_pymol(
        self,
        node_metrics:      List[str],
        community_methods: Optional[List[str]] = None,
        sphere_scale:      float = 2.0,
        stick_scale:       float = 0.5,
        reference_pdb:     Optional[Union[str, Path]] = None,
    ) -> List[Path]:
        """
        Write PyMOL scripts, separating node-centrality and community renders.

        Layout under ``Network/``::

            <inttype>/node_metrics/<source>_<metric>.pml
            <inttype>/community_metrics/<source>_<method>.pml

        Each node-metric file shows CA spheres scaled by that one metric, with
        gray network edges for connectivity. Each community file shows CA
        spheres coloured by that one partition; its intra-community edges take
        the community colour and cross-community edges stay gray. Node and
        community renders are independent, so the count is
        (metrics + methods) per (inttype × source), not their product.
        """
        if not self.graphs:
            raise RuntimeError('write_pymol() called before load_static()')
        if not self.centrality_dfs:
            raise RuntimeError('write_pymol() called before run_analysis()')

        ref_pdb = Path(reference_pdb) if reference_pdb is not None else self.reference_pdb
        if community_methods is None:
            first = next(iter(self.community_dfs.values()), {})
            community_methods = list(first.keys())

        out_files: List[Path] = []
        for key, G in self.graphs.items():
            inttype, source = key
            cdf_centrality = self.centrality_dfs.get(key)
            if cdf_centrality is None:
                continue

            node_dir = self.out / inttype / 'node_metrics'
            comm_dir = self.out / inttype / 'community_metrics'
            makedirs(node_dir, exist_ok=True)
            makedirs(comm_dir, exist_ok=True)

            # ── one file per node metric ─────────────────────────────────────
            for nm in node_metrics:
                if nm not in cdf_centrality.columns:
                    log.warning('node_metric %r not in centrality for %s/%s — skipped',
                                nm, inttype, source)
                    continue
                p = write_pymol_nodes(
                    G               = G,
                    pdb_path        = ref_pdb,
                    outpath         = node_dir,
                    node_metric     = nm,
                    centrality_df   = cdf_centrality,
                    resindex_to_pdb = self.resindex_to_pdb,
                    sphere_scale    = sphere_scale,
                    stick_scale     = stick_scale,
                    filename        = f'{source}_{nm}.pml',
                )
                out_files.append(p)
                log.info('PyMOL node metric: %s/%s', inttype, p.name)

            # ── one file per community method ────────────────────────────────
            for cm in community_methods:
                cdf_comm = self.community_dfs.get(key, {}).get(cm)
                if cdf_comm is None:
                    log.warning('community %r not run for %s/%s — skipped',
                                cm, inttype, source)
                    continue
                p = write_pymol_communities(
                    G                = G,
                    pdb_path         = ref_pdb,
                    outpath          = comm_dir,
                    community_df     = cdf_comm,
                    community_method = cm,
                    resindex_to_pdb  = self.resindex_to_pdb,
                    sphere_scale     = sphere_scale,
                    stick_scale      = stick_scale,
                    filename         = f'{source}_{cm}.pml',
                )
                out_files.append(p)
                log.info('PyMOL community: %s/%s', inttype, p.name)
        return out_files

    # ═════════════════════════════════════════════════════════════════════════
    # Time-resolved
    # ═════════════════════════════════════════════════════════════════════════

    def load_time_resolved_energy(
        self,
        interactors_type: List[str],
        targets:          List[str],
        window_size:      int,
        stride:           Optional[int] = None,
        energy_use_abs:   bool  = False,
        min_energy:       float = 0.0,
        frame_start:      Optional[int] = None,
        frame_stop:       Optional[int] = None,
    ) -> None:
        """
        Read the per-frame energy CSVs (one per merged file, sitting in the
        target sub-folder), window-aggregate, and build one graph per
        ``(inttype, 'energy', window_start)``.

        Parameters
        ----------
        window_size : int
            Frames per window.
        stride : int, optional
            Step between windows (defaults to ``window_size`` — non-overlapping).
        energy_use_abs : bool
            Per-row weight: ``|energy_total|`` if True, else ``max(0, -energy_total)``.
        min_energy : float
            Drop edges (per window) with weight ≤ this (percent of the strongest
            edge across all windows, 0..100). Energy strengths are normalised to
            0..100 using a global max so the threshold is comparable between
            windows. Windowed graphs are always energy-sourced.
        frame_start, frame_stop : int, optional
            Restrict the per-frame data to ``[frame_start, frame_stop)``. Either
            or both may be None — None means "no bound on that side". When both
            are None (default), all frames in the CSVs are analysed.
        """
        start = timeit.default_timer()
        stride = stride or window_size

        for inttype in interactors_type:
            folder = self.adjacency_root / inttype
            if not folder.exists():
                log.warning('load_time_resolved_energy: %s missing — skipped', folder)
                continue
            df = load_per_frame_energy(folder, targets=targets, use_abs=energy_use_abs)
            if df.empty:
                log.warning('load_time_resolved_energy: no per-frame energy CSVs '
                            'under %s — skipped (did you re-run contacts after '
                            'the streaming-energy patch?)', folder)
                continue

            # ── frame-range filter ──────────────────────────────────────────
            n_before = df['frame'].nunique()
            if frame_start is not None:
                df = df[df['frame'] >= frame_start]
            if frame_stop is not None:
                df = df[df['frame'] <  frame_stop]
            if df.empty:
                log.warning('load_time_resolved_energy: frame range '
                            '[%s, %s) excluded all data under %s — skipped',
                            frame_start, frame_stop, folder)
                continue
            n_after = df['frame'].nunique()
            if frame_start is not None or frame_stop is not None:
                log.info('  %s: frame range [%s, %s) → %d/%d frames retained',
                         inttype, frame_start, frame_stop, n_after, n_before)

            # canonical order: res_1 < res_2 (already done by writer, but be safe)
            swap = df['res_1'] > df['res_2']
            df.loc[swap, ['res_1', 'res_2']] = df.loc[swap, ['res_2', 'res_1']].values

            # window aggregation (sum weight + n_contacts within each window/pair)
            wdf = window_per_frame(
                df, window_size=window_size, stride=stride,
                value_cols=['weight', 'n_contacts'],
            )
            log.info('  %s: %d frames → %d (window,pair) rows',
                     inttype, df['frame'].nunique(), len(wdf))

            # normalise to 0..100 "% of strongest edge" using a GLOBAL max across
            # all windows, so min_energy means the same % in every window and the
            # windows stay comparable for community-evolution tracking.
            wmax = wdf['weight'].max() if len(wdf) else 0.0
            if wmax > 0:
                wdf = wdf.copy()
                wdf['weight'] = 100.0 * wdf['weight'] / wmax
            else:
                wdf = wdf.copy()
                wdf['weight'] = 0.0

            for window_start, gdf in wdf.groupby('window_start'):
                pairs = gdf[['res_1', 'res_2', 'weight']]
                pairs = pairs[pairs['weight'] > min_energy]
                if pairs.empty:
                    continue
                G = _build_graph_from_pairs(pairs, self.int_top, weight_col='weight')
                self.windowed_graphs[(inttype, 'energy', int(window_start))] = G

        proc_time, total_time = get_endtime(start, self.start)
        log.info('load_time_resolved_energy: %d windowed graphs | '
                 'took %s (h:m:s), total %s (h:m:s)',
                 len(self.windowed_graphs), proc_time, total_time)

    # ─────────────────────────────────────────────────────────────────────────

    def track_communities_over_time(
        self,
        community_methods: Optional[List[str]] = None,
        community_seed:    int                  = 42,
        min_window_edges:  int                  = 0,
        community_min_members: int              = 2,
    ) -> Dict[Tuple[str, str, str], Path]:
        """
        Detect communities on every windowed graph and write a tidy CSV per
        (inttype × source × community_method):

            community_evolution_<inttype>_<source>_<method>.csv
              window_start, resindex, label, community

        Use this as the input for the "community vs time" pivot/heatmap.
        """
        if not self.windowed_graphs:
            raise RuntimeError('track_communities_over_time() called before '
                               'load_time_resolved_energy()')
        if community_methods is None:
            community_methods = ['louvain']

        start = timeit.default_timer()
        # group windowed_graphs by (inttype, source)
        by_run: Dict[Tuple[str, str], List[Tuple[int, nx.Graph]]] = {}
        for (inttype, source, ws), G in self.windowed_graphs.items():
            by_run.setdefault((inttype, source), []).append((ws, G))

        written: Dict[Tuple[str, str, str], Path] = {}
        label_map = build_label_map(self.int_top)

        for (inttype, source), window_list in by_run.items():
            window_list.sort(key=lambda x: x[0])
            for method in community_methods:
                rows: List[Dict] = []
                for ws, G in window_list:
                    if G.number_of_edges() < min_window_edges:
                        continue
                    cdf, n_comm = detect_communities(
                        G, method=method, seed=community_seed,
                        min_members=community_min_members,
                    )
                    log.debug('  %s/%s window=%d (%s): %d communities',
                              inttype, source, ws, method, n_comm)
                    for _, row in cdf.iterrows():
                        rows.append({
                            'window_start': int(ws),
                            'resindex':     int(row['resindex']),
                            'label':        label_map.get(int(row['resindex']),
                                                          str(int(row['resindex']))),
                            'community':    int(row['community']),
                        })
                if not rows:
                    continue
                out_df = pd.DataFrame(rows)
                ev_dir = self.out / inttype / 'time_resolved'
                makedirs(ev_dir, exist_ok=True)
                out_path = ev_dir / f'community_evolution_{source}_{method}.csv'
                out_df.to_csv(out_path, index=False)
                written[(inttype, source, method)] = out_path
                log.info('community evolution: %s (%d (window,residue) rows)',
                         out_path.name, len(out_df))

        proc_time, total_time = get_endtime(start, self.start)
        log.info('track_communities_over_time complete | '
                 'took %s (h:m:s), total %s (h:m:s)',
                 proc_time, total_time)
        return written

    # ═════════════════════════════════════════════════════════════════════════
    # Auxiliary outputs
    # ═════════════════════════════════════════════════════════════════════════

    def write_outputs(
        self,
        graphml_per_source: bool = True,
    ) -> List[Path]:
        """One GraphML per static graph (or only the first if disabled)."""
        if not self.graphs:
            raise RuntimeError('write_outputs() called before load_static()')
        out_paths: List[Path] = []

        for key, G in self.graphs.items():
            inttype, source = key
            H = G.copy()
            if key in self.centrality_dfs:
                for _, row in self.centrality_dfs[key].iterrows():
                    n = int(row['resindex'])
                    if H.has_node(n):
                        H.nodes[n]['betweenness_centrality'] = row['betweenness_centrality']
                        H.nodes[n]['closeness_centrality']   = row['closeness_centrality']
                        H.nodes[n]['degree_centrality']      = row['degree_centrality']
                        H.nodes[n]['weighted_degree']        = row['weighted_degree']
            if key in self.community_dfs and self.community_dfs[key]:
                first_method = next(iter(self.community_dfs[key]))
                comm_map = dict(zip(
                    self.community_dfs[key][first_method]['resindex'],
                    self.community_dfs[key][first_method]['community'],
                ))
                for n in H.nodes():
                    H.nodes[n]['community'] = int(comm_map.get(int(n), -1))

            gml_dir = self.out / inttype
            makedirs(gml_dir, exist_ok=True)
            out_path = write_graphml(
                H, gml_dir, filename=f'network_{source}.graphml',
            )
            out_paths.append(out_path)
            log.info('GraphML: %s', out_path)
            if not graphml_per_source:
                break
        return out_paths

    # ─────────────────────────────────────────────────────────────────────────

    def get_shortest_paths(
        self,
        key:    Tuple[str, str],
        source: Optional[Union[int, List[int]]] = None,
        target: Optional[Union[int, List[int]]] = None,
        filename: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Weighted shortest paths on one chosen static graph.
        ``key`` is ``(interactors_type, edge_weight_source)``.
        """
        from mypsns.utils.network_utilities import shortest_paths as _sp
        if key not in self.graphs:
            raise KeyError(f'no static graph for {key!r}. Available: {list(self.graphs)}')
        df = _sp(self.graphs[key], source=source, target=target)
        inttype, src = key
        csv_dir = self.out / inttype / 'csv'
        makedirs(csv_dir, exist_ok=True)
        if filename is None:
            filename = f'shortest_paths_{src}.csv'
        out = csv_dir / filename
        df.to_csv(out, header=True, index=False)
        log.info('shortest_paths: %d pairs → %s', len(df), out)
        return df

    # ─────────────────────────────────────────────────────────────────────────

    def apply_filter(
        self,
        key:        Tuple[str, str],
        min_weight: Optional[float] = None,
        min_degree: Optional[int]   = None,
    ) -> nx.Graph:
        """Filtered copy of one static graph; ``self.graphs[key]`` untouched.

        ``min_weight`` is interpreted on the chosen graph's own scale
        (persistency % or energy kJ/mol), since ``key`` already identifies the
        source.
        """
        if key not in self.graphs:
            raise KeyError(f'no static graph for {key!r}')
        H = self.graphs[key].copy()
        if min_weight is not None:
            H = filter_by_weight(H, min_weight)
        if min_degree is not None:
            H = filter_by_degree(H, min_degree)
        log.info('apply_filter %s: min_w=%s min_d=%s → %d nodes, %d edges',
                 key, min_weight, min_degree,
                 H.number_of_nodes(), H.number_of_edges())
        return H