"""`mypsns compare-variants` — rank network metrics against ClinVar pathogenicity
(experimental; needs internet and the optional `variants` extra)."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional

import typer
from rich.markup import escape

from ..logging_utils import (
    print_banner, print_error, print_info, print_success, print_summary,
    print_warning, stage_log_file,
)
from ._helpers import (
    collect_overrides, dump_run_config, load_or_die, require_path, warn_skip,
)


log = logging.getLogger(__name__)

PANEL_GENERAL = "General"
PANEL_CLINVAR = "ClinVar query"
PANEL_SCORE   = "Scoring"
PANEL_NETWORK = "Network rebuild (overrides the network: section)"
PANEL_BEHAVE  = "Behaviour"


def compare_variants(
    config: Optional[Path] = typer.Option(
        None, "--config",
        help="YAML config file. Point it at the run_config.yaml that `network` "
             "wrote to rank exactly the networks it built.",
        rich_help_panel=PANEL_GENERAL,
    ),
    topology: Optional[Path] = typer.Option(
        None, "--topology", "-t",
        help="Topology file (used to build the Resindex→PDB-resid map).",
        rich_help_panel=PANEL_GENERAL,
    ),
    outpath: Optional[Path] = typer.Option(
        None, "--outpath", "-o",
        help="Root output folder containing Interactors/ and Adjacency/.",
        rich_help_panel=PANEL_GENERAL,
    ),
    gene: Optional[str] = typer.Option(
        None, "--gene",
        help="Gene symbol queried in ClinVar (e.g. RCVRN). Required.",
        rich_help_panel=PANEL_CLINVAR,
    ),
    reference: Optional[str] = typer.Option(
        None, "--reference",
        help=("RefSeq protein accession the ClinVar protein positions refer to "
              "(e.g. NP_002894.1), or a raw protein sequence. Required."),
        rich_help_panel=PANEL_CLINVAR,
    ),
    chrom: Optional[str] = typer.Option(
        None, "--chrom",
        help="Optional GRCh38 window: chromosome (17 or chr17). Needs --start/--stop.",
        rich_help_panel=PANEL_CLINVAR,
    ),
    start: Optional[int] = typer.Option(
        None, "--start",
        help="Optional GRCh38 window: first base (inclusive).",
        rich_help_panel=PANEL_CLINVAR,
    ),
    stop: Optional[int] = typer.Option(
        None, "--stop",
        help="Optional GRCh38 window: last base (inclusive).",
        rich_help_panel=PANEL_CLINVAR,
    ),
    email: Optional[str] = typer.Option(
        None, "--email", envvar="NCBI_EMAIL",
        help="Contact email sent to NCBI E-utilities (NCBI asks for one).",
        rich_help_panel=PANEL_CLINVAR,
    ),
    api_key: Optional[str] = typer.Option(
        None, "--api-key", envvar="NCBI_API_KEY",
        help="NCBI API key (10 requests/s instead of 3). Never written to "
             "run_config.yaml.",
        rich_help_panel=PANEL_CLINVAR,
    ),
    negative: Optional[str] = typer.Option(
        None, "--negative",
        help=("Negative class: 'rest' (default; every other residue aligned to "
              "the reference) or 'benign' (only sites with a benign ClinVar record)."),
        rich_help_panel=PANEL_SCORE,
    ),
    min_identity: Optional[float] = typer.Option(
        None, "--min-identity",
        help="Minimum percent identity for a chain to be taken as the gene's "
             "chain (default 70).",
        rich_help_panel=PANEL_SCORE,
    ),
    interactors_type: Optional[List[str]] = typer.Option(
        None, "--interactors-type",
        help="Interactor schema(s). Repeat for multiple.",
        rich_help_panel=PANEL_NETWORK,
    ),
    targets: Optional[List[str]] = typer.Option(
        None, "--target",
        help="Contact type(s) to include: HBonds, Salt_bridges, HCs.",
        rich_help_panel=PANEL_NETWORK,
    ),
    edge_weight_sources: Optional[List[str]] = typer.Option(
        None, "--edge-weight-source",
        help="'persistency' and/or 'energy'. Repeat for multiple.",
        rich_help_panel=PANEL_NETWORK,
    ),
    energy_use_abs: Optional[bool] = typer.Option(
        None, "--energy-use-abs/--no-energy-use-abs",
        help="Energy → weight: |E| instead of attractive-only max(0,-E).",
        rich_help_panel=PANEL_NETWORK,
    ),
    weighting: Optional[str] = typer.Option(
        None, "--weighting",
        help="'weighted', 'unweighted' or 'both' (default).",
        rich_help_panel=PANEL_NETWORK,
    ),
    min_persistency: Optional[float] = typer.Option(
        None, "--min-persistency",
        help="Persistency edge threshold (percent, 0..100).",
        rich_help_panel=PANEL_NETWORK,
    ),
    min_energy: Optional[float] = typer.Option(
        None, "--min-energy",
        help="Energy edge threshold (percent of strongest edge, 0..100).",
        rich_help_panel=PANEL_NETWORK,
    ),
    min_degree: Optional[int] = typer.Option(
        None, "--min-degree",
        help="Drop nodes with degree below this value.",
        rich_help_panel=PANEL_NETWORK,
    ),
    frame_start: Optional[int] = typer.Option(
        None, "--frame-start",
        help="First frame (inclusive); re-aggregates from per-frame energy CSVs.",
        rich_help_panel=PANEL_NETWORK,
    ),
    frame_stop: Optional[int] = typer.Option(
        None, "--frame-stop",
        help="Last frame (exclusive).",
        rich_help_panel=PANEL_NETWORK,
    ),
    overwrite: bool = typer.Option(
        False, "--overwrite",
        help="Overwrite this gene's tables if they already exist.",
        rich_help_panel=PANEL_BEHAVE,
    ),
    dry_run: bool = typer.Option(
        False, "--dry-run",
        help="Validate inputs and print the plan; do not query ClinVar.",
        rich_help_panel=PANEL_BEHAVE,
    ),
) -> None:
    """
    Rank network node/edge metrics by how well they separate pathogenic ClinVar
    missense sites from the rest of the protein (AUROC / AUPRC / Mann–Whitney).

    Rebuilds the static graphs from ``<outpath>/Adjacency/`` with the
    ``network:`` settings, maps ClinVar positions onto the structure by aligning
    to the RefSeq protein, and writes ``<gene>_*.csv`` tables under
    ``<outpath>/Network/<inttype>/variants/``.

    Experimental. Needs internet and ``pip install -e '.[variants]'``
    (requests + biopython). Not part of `mypsns run` — invoke after `network`.
    """
    # Two stage namespaces: the ClinVar/scoring options are this stage's own,
    # the graph-rebuild options override `network:` so a single config drives
    # both stages identically.
    overrides = collect_overrides({
        "topology":     topology,
        "outpath":      outpath,
        "gene":         gene,
        "reference":    reference,
        "chrom":        chrom,
        "start":        start,
        "stop":         stop,
        "email":        email,
        "api_key":      api_key,
        "negative":     negative,
        "min_identity": min_identity,
    }, stage_key="compare_variants")
    network_overrides = collect_overrides({
        "interactors_type":    interactors_type,
        "targets":             targets,
        "edge_weight_sources": edge_weight_sources,
        "energy_use_abs":      energy_use_abs,
        "weighting":           weighting,
        "min_persistency":     min_persistency,
        "min_energy":          min_energy,
        "min_degree":          min_degree,
        "frame_start":         frame_start,
        "frame_stop":          frame_stop,
    }, stage_key="network")
    overrides.update(network_overrides)

    cfg = load_or_die(config, overrides)
    cv_cfg, net_cfg = cfg.compare_variants, cfg.network

    require_path("topology", cfg.topology, must_exist=True)
    require_path("outpath",  cfg.outpath,  must_exist=True)
    for flag, value in (("gene", cv_cfg.gene), ("reference", cv_cfg.reference)):
        if not value:
            print_error(f"--{flag} is required (pass via CLI or set "
                        f"compare_variants.{flag} in --config).")
            raise typer.Exit(code=2)

    net_out = cfg.outpath / "Network"
    adj_in  = cfg.outpath / "Adjacency"
    int_in  = cfg.outpath / "Interactors"
    if not adj_in.exists():
        print_error(f"Adjacency folder not found in {cfg.outpath}. "
                    "Run `mypsns contacts` first.")
        raise typer.Exit(code=2)
    if not int_in.exists():
        print_error(f"Interactors folder not found in {cfg.outpath}. "
                    "Run `mypsns interactors` first.")
        raise typer.Exit(code=2)

    # Import before the dry-run exit so a missing optional extra is reported by
    # --dry-run too. The service raises an ImportError naming the extra.
    try:
        from ..services.compare_variants import CompareVariants, variant_file_prefix
    except ImportError as e:
        # escape(): the message names the '.[variants]' extra, which Rich
        # would otherwise swallow as a markup tag.
        print_error(escape(str(e)))
        raise typer.Exit(code=3)

    window = (f"chr{str(cv_cfg.chrom).lower().replace('chr', '')}:"
              f"{cv_cfg.start}-{cv_cfg.stop} (GRCh38)" if cv_cfg.chrom else "whole gene")
    print_banner("compare-variants", subtitle=f"{cv_cfg.gene} → {net_out}/<inttype>/variants")
    print_summary("Inputs", {
        "topology":    cfg.topology,
        "outpath":     cfg.outpath,
        "Adjacency":   adj_in,
        "Interactors": int_in,
    })
    print_summary("ClinVar query", {
        "gene":      cv_cfg.gene,
        "reference": cv_cfg.reference,
        "window":    window,
        "email":     cv_cfg.email or "not set",
        "api_key":   "set" if cv_cfg.api_key else "not set",
    })
    print_summary("Scoring", {
        "negative":     cv_cfg.negative,
        "min_identity": cv_cfg.min_identity,
    })
    print_summary("Network rebuild", {
        "interactors_type":    net_cfg.interactors_type,
        "targets":             net_cfg.targets,
        "edge_weight_sources": net_cfg.edge_weight_sources,
        "energy_use_abs":      net_cfg.energy_use_abs,
        "weighting":           net_cfg.weighting,
        "min_persistency":     net_cfg.min_persistency,
        "min_energy":          net_cfg.min_energy,
        "min_degree":          net_cfg.min_degree,
        "frame_start":         net_cfg.frame_start if net_cfg.frame_start is not None else "start",
        "frame_stop":          net_cfg.frame_stop  if net_cfg.frame_stop  is not None else "end",
    })

    if dry_run:
        print_info("Dry run — nothing executed.")
        raise typer.Exit()

    prefix = variant_file_prefix(cv_cfg.gene)
    for inttype in net_cfg.interactors_type:
        marker = net_out / inttype / "variants" / f"{prefix}_clinvar_sites.csv"
        if warn_skip(marker, "compare-variants", overwrite):
            raise typer.Exit()

    try:
        with stage_log_file(net_out / "compare_variants.log"):
            cv = CompareVariants(
                topology_pdb    = cfg.topology,
                general_outpath = cfg.outpath,
                n_frames        = net_cfg.n_frames,
                na_values       = net_cfg.na_values,
            )

            # Same graph construction as `mypsns network` (static path only).
            cv.load_static(
                interactors_type    = net_cfg.interactors_type,
                targets             = net_cfg.targets,
                edge_weight_sources = net_cfg.edge_weight_sources,
                energy_use_abs      = net_cfg.energy_use_abs,
                min_persistency     = net_cfg.min_persistency,
                min_energy          = net_cfg.min_energy,
                weighting           = net_cfg.weighting,
                frame_start         = net_cfg.frame_start,
                frame_stop          = net_cfg.frame_stop,
            )
            if not cv.graphs:
                raise RuntimeError("no graphs were built — check --interactors-type "
                                   "and the Adjacency/ sub-folders.")
            if net_cfg.min_degree and net_cfg.min_degree > 0:
                for key in list(cv.graphs.keys()):
                    cv.graphs[key] = cv.apply_filter(key, min_degree=net_cfg.min_degree)

            results = cv.rank_against_clinvar(
                gene         = cv_cfg.gene,
                reference    = cv_cfg.reference,
                chrom        = cv_cfg.chrom,
                start        = cv_cfg.start,
                stop         = cv_cfg.stop,
                email        = cv_cfg.email,
                api_key      = cv_cfg.api_key,
                negative     = cv_cfg.negative,
                min_identity = cv_cfg.min_identity,
            )
            written = cv.write_rankings(results, gene=cv_cfg.gene)

    except Exception as e:
        log.exception("compare-variants failed")
        print_error(f"compare-variants failed: {escape(str(e))}")
        raise typer.Exit(code=1)

    dump_run_config(cfg, cfg.outpath)

    sites = cv.clinvar_sites
    mapped = sites[sites["resindex"].notna()]
    print_summary("ClinVar sites (missense)", {
        "records":                    len(cv.clinvar_variants),
        "missense sites":             len(sites),
        "on modelled residues":       len(mapped),
        "non-benign / benign / VUS":  " / ".join(
            str(int((mapped["label"] == lab).sum()))
            for lab in ("non_benign", "benign", "vus")),
        "aligned chain":              f"{cv.alignment['chain']} "
                                      f"({cv.alignment['identity']}% identity)",
    })
    best = {}
    for (inttype, label), ranking in results.items():
        if ranking.empty:
            best[f"{inttype}/{label}"] = "nothing to rank"
            continue
        top = ranking.iloc[0]
        best[f"{inttype}/{label}"] = (
            f"{top['metric']}  AUROC {top['auroc_oriented']:.3f} ({top['direction']}), "
            f"AUPRC {top['auprc']:.3f}, n={top['n_pos']}/{top['n_neg']}")
    print_summary("Best metric per network", best)
    if all(r.empty for r in results.values()):
        print_warning(f"No pathogenic / likely-pathogenic missense site of "
                      f"{cv_cfg.gene} falls on a modelled residue, so every "
                      "ranking is empty. The site and ClinVar tables are still written.")

    print_success(f"compare-variants complete → {len(written)} files under "
                  f"{net_out}/<inttype>/variants/ (prefix '{prefix}_')")
