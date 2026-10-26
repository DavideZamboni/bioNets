"""
compare_variants.py
───────────────────
`CompareVariants(Network)` — overlay ClinVar pathogenicity onto a myPSNs network
and rank node/edge metrics by how well they separate non-benign missense sites.

It is a **subclass of `Network`**, so it inherits the loaded interactor table
(`self.int_top`, incl. the `Entity` P/M column), the `Resindex → (chain, resid)`
map, and — once you call the inherited `load_static()` / `run_analysis()` — the
in-memory `self.graphs` and `self.centrality_dfs`. No reloading, no CSV round-trip.

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

Suggested placement: a separate `mypsns compare-variants` subcommand run *after*
`network` (it needs internet + a RefSeq protein accession and is an optional
validation layer) — a thin Typer wrapper around `rank_against_clinvar()`.

Deps: requests pandas numpy scipy biopython  (MDAnalysis via Network)
"""
from __future__ import annotations
import re
import numpy as np
import pandas as pd
from mypsns.services.network import Network
from mypsns.utils.metric_comparison_utilities import score, node_metric_table, missense_truth_table, _auroc, _auprc
from mypsns.utils.clinvar_query_utilities import fetch_clinvar_variants, fetch_reference_protein

class CompareVariants(Network):
    """A `Network` that can also rank its metrics against ClinVar pathogenicity.

    Usage:
        cv = CompareVariants(topology_pdb=..., general_outpath=..., ...)
        cv.load_static(...)         # inherited
        cv.run_analysis(...)        # inherited → populates graphs + metrics
        ranked = cv.rank_against_clinvar(gene="RCVRN", reference="NP_002894.1",
                                         chrom="chr17", start=9898189, stop=9905271)
    """

    NODE_METRICS = ["degree", "weighted_degree", "degree_centrality",
                    "betweenness_centrality", "closeness_centrality"]
    EDGE_FEATURES = ["edge_strength", "edge_count", "edge_weight_mean", "edge_weight_max"]

    # ── orchestration ────────────────────────────────────────────────────────
    def rank_against_clinvar(self, gene, reference, chrom=None, start=None,
                             stop=None, assembly="GRCh38", email=None,
                             api_key=None, negative="rest", verbose=True):
        if not self.graphs:
            raise RuntimeError("Call load_static()+run_analysis() first.")

        variants = fetch_clinvar_variants(gene, chrom, start, stop, assembly,
                                          email=email, api_key=api_key)
        truth = missense_truth_table(variants)
        if truth.empty:
            raise ValueError("No missense variants with parseable positions.")

        ref_seq = (reference if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYX]+", reference or "", re.I)
                   else fetch_reference_protein(reference, email=email, api_key=api_key))
        chain, mapping, ident = self.map_reference_positions(ref_seq)
        pos2ri = {p: rr[0] for p, rr in mapping.items()}

        truth = truth.reset_index()
        truth["resindex"] = truth["ref_pos"].map(pos2ri)
        truth = truth.dropna(subset=["resindex"])
        truth["resindex"] = truth["resindex"].astype(int)
        if verbose:
            print(f"[clinvar] {len(variants)} variants → "
                  f"{(truth.label=='non_benign').sum()} non-benign / "
                  f"{(truth.label=='benign').sum()} benign sites")
            print(f"[align] {gene} → chain {chain} @ {ident:.1f}% identity")

        results = {}
        for key in self.graphs:
            metrics = node_metric_table(self.graphs[key], self.centrality_dfs[key].copy())
            merged = metrics.merge(truth[["resindex", "label"]], on="resindex", how="left")
            merged["label"] = merged["label"].fillna("rest")
            results[key] = score(merged, negative, self.NODE_METRICS, self.EDGE_FEATURES)
            if verbose:
                print(f"\n=== {key} (negative='{negative}') ===")
                print(results[key].to_string(index=False))
        return results