"""
compare_variants.py
───────────────────
`CompareVariants(Network)` — overlay ClinVar pathogenicity onto a myPSNs network
and rank node/edge metrics by how well they separate pathogenic missense sites.

It is a **subclass of `Network`**, so it inherits the loaded interactor table
(`self.int_top`, incl. the `Entity` P/M column), the `Resindex → (chain, resid)`
map, and — once you call the inherited `load_static()` — the in-memory
`self.graphs`. No reloading, no CSV round-trip.

Design notes addressed here:
  • Protein-chain residues are selected via the `Entity` column ('P' = polymeric),
    NOT a hand-maintained ion/water blocklist.
  • 3→1 letter coding is an inlined lookup table (the standard 20 + IUPAC
    ambiguity codes + forcefield variants); Biopython is used only for the
    pairwise alignment, where a maintained library is worth it.
  • Positions are mapped to nodes by SEQUENCE ALIGNMENT (robust to offsets, gaps,
    insertion codes, and heterodimers — auto-picks the chain that is the gene).
  • ClinVar is queried purely through NCBI E-utilities (esearch/esummary; efetch
    for the reference protein).
  • Only residues aligned to the reference are scored. Partner chains, ions and
    ligands carry no ClinVar information for this gene, so counting them as
    negatives (negative='rest') would be unfounded and would inflate AUROC.
  • Centralities come from `run_analysis()` if it has been run; otherwise they
    are computed in memory, so ranking has no CSV/community side effects on the
    `network` outputs.

Labels (per reference position, from missense_truth_table): 'non_benign' when
ClinVar holds at least one pathogenic / likely-pathogenic missense record there,
'benign' when only benign / likely-benign ones, 'vus' otherwise. Positives are
the 'non_benign' sites; negatives are every other aligned residue
(negative='rest') or only the 'benign' sites (negative='benign').

Driven by `mypsns compare-variants` (commands/compare_variants_cmd.py), run
after `network`: it needs internet, a RefSeq protein accession, and the optional
`variants` extra (requests + biopython).
"""
from __future__ import annotations

import logging
import re
from importlib.util import find_spec
from os import makedirs
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

# requests (E-utilities) and biopython (alignment) are the optional `variants`
# extra. Check both here, before importing the helpers: otherwise a missing
# requests is a bare ModuleNotFoundError, and a missing Bio (imported lazily by
# the aligner) only fails after the network has been rebuilt and ClinVar queried.
_missing = [pkg for mod, pkg in (("requests", "requests"), ("Bio", "biopython"))
            if find_spec(mod) is None]
if _missing:
    raise ImportError(
        f"compare-variants needs the optional 'variants' extra (missing: "
        f"{', '.join(_missing)}). Install it with `pip install -e '.[variants]'` "
        f"from the myPSNs checkout, or `pip install requests biopython`."
    )

from mypsns.services.network import Network
from mypsns.utils.clinvar_query_utilities import fetch_clinvar_variants, fetch_reference_protein
from mypsns.utils.metric_comparison_utilities import (
    map_reference_positions, missense_truth_table, node_metric_table, score,
)
from mypsns.utils.network_utilities import compute_centrality

log = logging.getLogger(__name__)

_PROTEIN_SEQ = re.compile(r"[ACDEFGHIKLMNPQRSTVWYX]+", re.I)


def variant_file_prefix(gene: str) -> str:
    """File-name prefix for one gene's outputs. Prefixing (rather than one
    file per inttype/label) lets several genes of a complex — e.g. RCVRN and
    GRK1 — share Network/<inttype>/variants/ without overwriting each other."""
    return re.sub(r"[^A-Za-z0-9_.-]", "_", gene.strip())


class CompareVariants(Network):
    """A `Network` that can also rank its metrics against ClinVar pathogenicity.

    Usage:
        cv = CompareVariants(topology_pdb=..., general_outpath=..., ...)
        cv.load_static(...)         # inherited → populates self.graphs
        ranked = cv.rank_against_clinvar(gene="RCVRN", reference="NP_002894.1")
        cv.write_rankings(ranked, gene="RCVRN")
    """

    NODE_METRICS = ["degree", "weighted_degree", "degree_centrality",
                    "betweenness_centrality", "closeness_centrality"]
    EDGE_FEATURES = ["edge_strength", "edge_count", "edge_weight_mean", "edge_weight_max"]

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        # Filled by rank_against_clinvar(), written by write_rankings().
        self.clinvar_variants: pd.DataFrame = pd.DataFrame()   # raw esummary rows
        self.clinvar_sites:    pd.DataFrame = pd.DataFrame()   # per missense site
        self.site_tables: Dict[Tuple[str, str], pd.DataFrame] = {}  # per graph
        self.alignment:   Dict[str, object] = {}

    # ── orchestration ────────────────────────────────────────────────────────
    def rank_against_clinvar(self, gene: str, reference: str,
                             chrom: Optional[str] = None, start: Optional[int] = None,
                             stop: Optional[int] = None, email: Optional[str] = None,
                             api_key: Optional[str] = None, negative: str = "rest",
                             min_identity: float = 70.0,
                             ) -> Dict[Tuple[str, str], pd.DataFrame]:
        """Query ClinVar for `gene`, map its missense sites onto the structure by
        aligning to `reference`, and score every graph in ``self.graphs``.

        `reference` is a RefSeq protein accession (efetched; e.g. NP_002894.1)
        or a raw protein sequence. chrom/start/stop is an optional GRCh38 window
        on top of the gene query.

        Returns ``{(inttype, label): ranking}``; a ranking is empty (header only)
        when no pathogenic site falls on an aligned residue.
        """
        if not self.graphs:
            raise RuntimeError("Call load_static() first.")
        if negative not in ("rest", "benign"):
            raise ValueError(f"negative must be 'rest' or 'benign', got {negative!r}")

        self.clinvar_variants = fetch_clinvar_variants(gene, chrom, start, stop,
                                                       email=email, api_key=api_key)
        truth = missense_truth_table(self.clinvar_variants)
        if truth.empty:
            raise ValueError(f"ClinVar returned {len(self.clinvar_variants)} records "
                             f"for {gene} but no missense variant with a parseable "
                             "protein position.")

        ref_seq = (reference if _PROTEIN_SEQ.fullmatch(reference or "")
                   else fetch_reference_protein(reference, email=email, api_key=api_key))
        if not ref_seq:
            raise ValueError(f"Empty reference protein sequence for {reference!r}.")
        chain, mapping, ident = map_reference_positions(
            self.int_top, self.resindex_to_pdb, ref_seq, min_identity=min_identity)

        # One row per ClinVar missense site, with the structure residue it
        # aligns to (NaN when the reference position is not modelled).
        aligned = pd.DataFrame(
            [(pos, *res) for pos, res in mapping.items()],
            columns=["ref_pos", "resindex", "resnum", "resname", "structure_aa"])
        aligned["ref_aa"] = [ref_seq[p - 1] for p in aligned["ref_pos"]]
        sites = truth.reset_index().merge(aligned, on="ref_pos", how="left")
        sites["ref_aa"] = [ref_seq[p - 1] if 0 < p <= len(ref_seq) else ""
                           for p in sites["ref_pos"]]
        for col in ("resindex", "resnum"):
            sites[col] = sites[col].astype("Int64")
        sites.insert(sites.columns.get_loc("resindex"), "chain",
                     sites["resindex"].notna().map({True: chain, False: ""}))
        self.clinvar_sites = sites

        mapped = sites.dropna(subset=["resindex"])
        n_pos = int((mapped["label"] == "non_benign").sum())
        self.alignment = {"gene": gene, "reference": reference, "chain": chain,
                          "identity": round(float(ident), 1),
                          "aligned_residues": len(mapping)}
        log.info("[clinvar] %s: %d records → %d missense sites (%d non-benign, "
                 "%d benign, %d VUS)", gene, len(self.clinvar_variants), len(sites),
                 int((sites["label"] == "non_benign").sum()),
                 int((sites["label"] == "benign").sum()),
                 int((sites["label"] == "vus").sum()))
        log.info("[align] %s → chain %s @ %.1f%% identity; %d/%d sites on modelled "
                 "residues (%d non-benign)", reference, chain, ident, len(mapped),
                 len(sites), n_pos)
        if n_pos == 0:
            log.warning("No pathogenic / likely-pathogenic missense site of %s falls "
                        "on a modelled residue: rankings will be empty.", gene)

        # Every residue aligned to the reference, with its ClinVar class
        # ('rest' where ClinVar has no missense record). The inner merge below
        # restricts scoring to these residues.
        labels = (aligned[["resindex", "ref_pos", "ref_aa"]]
                  .merge(sites[["ref_pos", "label"]], on="ref_pos", how="left")
                  .astype({"resindex": int}))
        labels["label"] = labels["label"].fillna("rest")

        results: Dict[Tuple[str, str], pd.DataFrame] = {}
        for key, G in self.graphs.items():
            cdf = self.centrality_dfs.get(key)
            cdf = compute_centrality(G) if cdf is None else cdf.copy()
            # centrality's 'label' is the residue label; score() reads the
            # ClinVar class from 'label', so rename before merging.
            metrics = node_metric_table(G, cdf).rename(columns={"label": "residue"})
            merged = metrics.merge(labels, on="resindex", how="inner")
            self.site_tables[key] = merged
            results[key] = score(merged, negative, self.NODE_METRICS, self.EDGE_FEATURES)
            log.info("=== %s/%s (negative='%s') ===\n%s", key[0], key[1], negative,
                     results[key].to_string(index=False) if not results[key].empty
                     else "(nothing to rank)")
        return results

    # ── output ───────────────────────────────────────────────────────────────
    def write_rankings(self, results: Dict[Tuple[str, str], pd.DataFrame],
                       gene: str) -> List[Path]:
        """Write the tables of one rank_against_clinvar() call under
        ``Network/<inttype>/variants/``, every file prefixed with the gene:

          <gene>_ranking_<label>.csv   one row per metric (AUROC, AUPRC, MWU p)
          <gene>_nodes_<label>.csv     per aligned residue: metrics + ClinVar label
          <gene>_clinvar_sites.csv     per missense site: counts, label, mapping
          <gene>_clinvar_variants.csv  raw ClinVar records (a dated snapshot —
                                       ClinVar classifications change over time)

        The two clinvar tables are the same for every inttype; they are repeated
        so each variants/ folder is self-contained.
        """
        prefix = variant_file_prefix(gene)
        written: List[Path] = []
        for inttype in dict.fromkeys(k[0] for k in results):
            vdir = self.out / inttype / "variants"
            makedirs(vdir, exist_ok=True)
            for name, df in (("clinvar_variants", self.clinvar_variants),
                             ("clinvar_sites", self.clinvar_sites)):
                path = vdir / f"{prefix}_{name}.csv"
                df.to_csv(path, index=False)
                written.append(path)
        for (inttype, label), ranking in results.items():
            vdir = self.out / inttype / "variants"
            path = vdir / f"{prefix}_ranking_{label}.csv"
            ranking.to_csv(path, index=False)
            written.append(path)
            path = vdir / f"{prefix}_nodes_{label}.csv"
            self.site_tables[(inttype, label)].to_csv(path, index=False)
            written.append(path)
        log.info("write_rankings: %d files under %s/<inttype>/variants/",
                 len(written), self.out)
        return written
