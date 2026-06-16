"""`mypsns network` — build & visualise residue-interaction networks
(persistency- and/or energy-weighted; optionally time-resolved)."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional

import typer

from ..logging_utils import (
    print_banner, print_error, print_info, print_success, print_summary,
    stage_log_file,
)
from ._helpers import (
    collect_overrides, dump_run_config, ensure_dir, load_or_die,
    require_path, warn_skip,
)


log = logging.getLogger(__name__)

PANEL_GENERAL = "General"
PANEL_SELECT  = "Selectors"
PANEL_WEIGHT  = "Edge weight"
PANEL_FRAMES  = "Frame range"
PANEL_TIME    = "Time-resolved"
PANEL_METRICS = "Output metrics"
PANEL_FILTER  = "Filters"
PANEL_PYMOL   = "PyMOL rendering"
PANEL_BEHAVE  = "Behaviour"


def network(
    config: Optional[Path] = typer.Option(
        None, "--config", help="YAML config file.", rich_help_panel=PANEL_GENERAL,
    ),
    topology: Optional[Path] = typer.Option(
        None, "--topology", "-t",
        help="Topology file (used to build Resindex→PDB-resid map for PyMOL).",
        rich_help_panel=PANEL_GENERAL,
    ),
    outpath: Optional[Path] = typer.Option(
        None, "--outpath", "-o",
        help="Root output folder containing Interactors/ and Adjacency/.",
        rich_help_panel=PANEL_GENERAL,
    ),
    n_frames: Optional[int] = typer.Option(
        None, "--n-frames", "-n",
        help="Total trajectory frames analysed during Contacts (informational; "
             "logged only).",
        rich_help_panel=PANEL_GENERAL,
    ),
    interactors_type: Optional[List[str]] = typer.Option(
        None, "--interactors-type",
        help="Interactor schema(s). Repeat for multiple.",
        rich_help_panel=PANEL_SELECT,
    ),
    targets: Optional[List[str]] = typer.Option(
        None, "--target",
        help="Contact type(s) to include: HBonds, Salt_bridges, HCs.",
        rich_help_panel=PANEL_SELECT,
    ),
    edge_weight_sources: Optional[List[str]] = typer.Option(
        None, "--edge-weight-source",
        help=("Build one network per source: 'persistency' (frequency) "
              "or 'energy' (LJ+Coulomb). Repeat for multiple."),
        rich_help_panel=PANEL_WEIGHT,
    ),
    energy_use_abs: Optional[bool] = typer.Option(
        None, "--energy-use-abs/--no-energy-use-abs",
        help=("Energy → weight transform. Off (default): only attractive "
              "(weight=max(0,-E)). On: |E|, including repulsive."),
        rich_help_panel=PANEL_WEIGHT,
    ),
    window_size: Optional[int] = typer.Option(
        None, "--window-size",
        help=("Also build per-window energy networks of this many frames "
              "and write community_evolution_*.csv."),
        rich_help_panel=PANEL_TIME,
    ),
    window_stride: Optional[int] = typer.Option(
        None, "--window-stride",
        help="Step between windows. Defaults to window_size (non-overlapping).",
        rich_help_panel=PANEL_TIME,
    ),
    frame_start: Optional[int] = typer.Option(
        None, "--frame-start",
        help=("First frame (inclusive) for both static and time-resolved "
              "analysis. When set, static networks are re-aggregated from "
              "the per-frame energy CSVs instead of the pre-aggregated "
              "matrices. Default: start of trajectory."),
        rich_help_panel=PANEL_FRAMES,
    ),
    frame_stop: Optional[int] = typer.Option(
        None, "--frame-stop",
        help=("Last frame (exclusive) for both static and time-resolved "
              "analysis. Default: end of trajectory."),
        rich_help_panel=PANEL_FRAMES,
    ),
    node_metrics: Optional[List[str]] = typer.Option(
        None, "--node-metric",
        help=("Node centrality metric(s) for CA-sphere scaling. Repeat. "
              "One .pml is produced per (source × node_metric × community)."),
        rich_help_panel=PANEL_METRICS,
    ),
    community_methods: Optional[List[str]] = typer.Option(
        None, "--community",
        help=("Community detection algorithm(s). Repeat. Choices: louvain, "
              "greedy_modularity, label_propagation, connected_components."),
        rich_help_panel=PANEL_METRICS,
    ),
    community_seed: Optional[int] = typer.Option(
        None, "--community-seed",
        help="Seed for stochastic community algorithms.",
        rich_help_panel=PANEL_METRICS,
    ),
    min_persistency: Optional[float] = typer.Option(
        None, "--min-persistency",
        help="Drop persistency edges with summed weight \u2264 this (percent, 0..100).",
        rich_help_panel=PANEL_FILTER,
    ),
    min_energy: Optional[float] = typer.Option(
        None, "--min-energy",
        help="Drop energy edges with summed weight \u2264 this (kJ/mol, post-transform).",
        rich_help_panel=PANEL_FILTER,
    ),
    min_degree: Optional[int] = typer.Option(
        None, "--min-degree",
        help="Drop nodes with degree below this value.",
        rich_help_panel=PANEL_FILTER,
    ),
    reference_pdb: Optional[Path] = typer.Option(
        None, "--reference-pdb",
        help="PDB to load in PyMOL (default: --topology).",
        rich_help_panel=PANEL_PYMOL,
    ),
    sphere_scale: Optional[float] = typer.Option(
        None, "--sphere-scale",
        help="Max CA sphere vdw radius [Å].",
        rich_help_panel=PANEL_PYMOL,
    ),
    stick_scale: Optional[float] = typer.Option(
        None, "--stick-scale",
        help="Max edge stick (dash) radius [Å].",
        rich_help_panel=PANEL_PYMOL,
    ),
    write_graphml: Optional[bool] = typer.Option(
        None, "--write-graphml/--no-write-graphml",
        help="Also write a GraphML per static network.",
        rich_help_panel=PANEL_PYMOL,
    ),
    shortest_paths: Optional[bool] = typer.Option(
        None, "--shortest-paths/--no-shortest-paths",
        help="Compute all-pairs shortest paths on each static graph.",
        rich_help_panel=PANEL_PYMOL,
    ),
    resume: bool = typer.Option(
        False, "--resume",
        help="Skip if Network/ already populated.",
        rich_help_panel=PANEL_BEHAVE,
    ),
    overwrite: bool = typer.Option(
        False, "--overwrite",
        help="Force recomputation.",
        rich_help_panel=PANEL_BEHAVE,
    ),
    dry_run: bool = typer.Option(
        False, "--dry-run",
        help="Validate inputs and print the plan; do not run.",
        rich_help_panel=PANEL_BEHAVE,
    ),
) -> None:
    """
    Build residue-interaction networks from Contacts output, render PyMOL
    sessions, and (optionally) track community evolution in time.

    Reads ``<outpath>/Adjacency/`` + ``<outpath>/Interactors/`` and writes
    node / edge CSVs, GraphML, and one .pml per
    (interactors_type × edge_weight_source × node_metric × community_method)
    under ``<outpath>/Network/``. With ``--window-size``, also writes
    ``community_evolution_*.csv``.

    Not part of `mypsns run` — invoke explicitly.
    """
    overrides = collect_overrides({
        "topology":            topology,
        "outpath":             outpath,
        "n_frames":            n_frames,
        "interactors_type":    interactors_type,
        "targets":             targets,
        "edge_weight_sources": edge_weight_sources,
        "energy_use_abs":      energy_use_abs,
        "window_size":         window_size,
        "window_stride":       window_stride,
        "frame_start":         frame_start,
        "frame_stop":          frame_stop,
        "node_metrics":        node_metrics,
        "community_methods":   community_methods,
        "community_seed":      community_seed,
        "min_persistency":     min_persistency,
        "min_energy":          min_energy,
        "min_degree":          min_degree,
        "reference_pdb":       reference_pdb,
        "sphere_scale":        sphere_scale,
        "stick_scale":         stick_scale,
        "write_graphml":       write_graphml,
        "shortest_paths":      shortest_paths,
    }, stage_key="network")

    cfg = load_or_die(config, overrides)

    require_path("topology", cfg.topology, must_exist=True)
    require_path("outpath",  cfg.outpath,  must_exist=True)

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

    ref_pdb = cfg.network.reference_pdb or cfg.topology

    print_banner("network", subtitle=str(net_out))
    print_summary("Inputs", {
        "topology":      cfg.topology,
        "outpath":       cfg.outpath,
        "Adjacency":     adj_in,
        "Interactors":   int_in,
        "n_frames":      cfg.network.n_frames if cfg.network.n_frames else "all",
        "reference_pdb": ref_pdb,
    })
    print_summary("Selectors", {
        "interactors_type":     cfg.network.interactors_type,
        "targets":              cfg.network.targets,
        "edge_weight_sources":  cfg.network.edge_weight_sources,
        "energy_use_abs":       cfg.network.energy_use_abs,
    })
    if cfg.network.frame_start is not None or cfg.network.frame_stop is not None:
        print_summary("Frame range", {
            "frame_start": cfg.network.frame_start if cfg.network.frame_start is not None else "start",
            "frame_stop":  cfg.network.frame_stop  if cfg.network.frame_stop  is not None else "end",
            "note":        "Static path re-aggregates from per-frame energy CSVs.",
        })
    if cfg.network.window_size:
        print_summary("Time-resolved", {
            "window_size":   cfg.network.window_size,
            "window_stride": cfg.network.window_stride or cfg.network.window_size,
            "frame_start":   cfg.network.frame_start if cfg.network.frame_start is not None else "start",
            "frame_stop":    cfg.network.frame_stop  if cfg.network.frame_stop  is not None else "end",
        })
    print_summary("Output metrics", {
        "node_metrics":      cfg.network.node_metrics,
        "community_methods": cfg.network.community_methods,
        "community_seed":    cfg.network.community_seed,
    })
    print_summary("Filters", {
        "min_persistency": cfg.network.min_persistency,
        "min_energy":      cfg.network.min_energy,
        "min_degree":      cfg.network.min_degree,
    })
    print_summary("PyMOL", {
        "sphere_scale":   cfg.network.sphere_scale,
        "stick_scale":    cfg.network.stick_scale,
        "write_graphml":  cfg.network.write_graphml,
        "shortest_paths": cfg.network.shortest_paths,
    })

    if dry_run:
        print_info("Dry run — nothing executed.")
        raise typer.Exit()

    if resume and net_out.exists() and any(net_out.glob("*.pml")) and not overwrite:
        print_info(f"network: .pml outputs already present in {net_out} — skipping (resume).")
        raise typer.Exit()
    if warn_skip(net_out, "network", overwrite):
        raise typer.Exit()

    ensure_dir(net_out)

    try:
        from ..services.network import Network  # type: ignore
    except ImportError as e:
        print_error(
            "Could not import `mypsns.services.network.Network`. "
            f"Underlying error: {e}"
        )
        raise typer.Exit(code=3)

    try:
        with stage_log_file(net_out / "network.log"):
            net = Network(
                topology_pdb    = cfg.topology,
                general_outpath = cfg.outpath,
                n_frames        = cfg.network.n_frames,
                na_values       = cfg.network.na_values,
            )

            # ── Static networks (one per inttype × source) ───────────────────
            net.load_static(
                interactors_type    = cfg.network.interactors_type,
                targets             = cfg.network.targets,
                edge_weight_sources = cfg.network.edge_weight_sources,
                energy_use_abs      = cfg.network.energy_use_abs,
                min_persistency     = cfg.network.min_persistency,
                min_energy          = cfg.network.min_energy,
                frame_start         = cfg.network.frame_start,
                frame_stop          = cfg.network.frame_stop,
            )

            # Optional per-graph degree filter
            if cfg.network.min_degree and cfg.network.min_degree > 0:
                for key in list(net.graphs.keys()):
                    net.graphs[key] = net.apply_filter(
                        key, min_degree=cfg.network.min_degree
                    )

            net.run_analysis(
                community_methods = cfg.network.community_methods,
                community_seed    = cfg.network.community_seed,
            )

            written_pml = net.write_pymol(
                node_metrics      = cfg.network.node_metrics,
                community_methods = cfg.network.community_methods,
                sphere_scale      = cfg.network.sphere_scale,
                stick_scale       = cfg.network.stick_scale,
                reference_pdb     = cfg.network.reference_pdb,
            )

            if cfg.network.write_graphml:
                net.write_outputs(graphml_per_source=True)

            if cfg.network.shortest_paths:
                for key in net.graphs.keys():
                    net.get_shortest_paths(key=key)

            # ── Time-resolved (optional) ─────────────────────────────────────
            written_evolution = {}
            if cfg.network.window_size:
                net.load_time_resolved_energy(
                    interactors_type = cfg.network.interactors_type,
                    targets          = cfg.network.targets,
                    window_size      = cfg.network.window_size,
                    stride           = cfg.network.window_stride,
                    energy_use_abs   = cfg.network.energy_use_abs,
                    min_energy       = cfg.network.min_energy,
                    frame_start      = cfg.network.frame_start,
                    frame_stop       = cfg.network.frame_stop,
                )
                written_evolution = net.track_communities_over_time(
                    community_methods = cfg.network.community_methods,
                    community_seed    = cfg.network.community_seed,
                )

    except Exception as e:
        log.exception("network failed")
        print_error(f"network failed: {e}")
        raise typer.Exit(code=1)

    dump_run_config(cfg, cfg.outpath)

    pml_list = "\n  ".join(str(p) for p in written_pml) if written_pml else "(none)"
    msg = f"network complete → {net_out}\n  PyMOL sessions:\n  {pml_list}"
    if written_evolution:
        ev_list = "\n  ".join(str(p) for p in written_evolution.values())
        msg += f"\n  Community evolution:\n  {ev_list}"
    print_success(msg)