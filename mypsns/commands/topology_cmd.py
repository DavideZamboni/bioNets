"""`mypsns topology` — match topology residues to forcefield entries."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import typer

from ..logging_utils import (
    print_banner, print_error, print_info, print_success, print_summary, stage_log_file,
)
from ._helpers import (
    collect_overrides, dump_run_config, ensure_dir,
    load_or_die, require_path, warn_skip,
)


log = logging.getLogger(__name__)

PANEL_GENERAL = "General"
PANEL_TOPOPTS = "Recognition options"
PANEL_BEHAVE  = "Behaviour"


def topology_recognition(
    config: Optional[Path] = typer.Option(
        None, "--config",
        help="YAML config file. CLI args override config values.",
        rich_help_panel=PANEL_GENERAL,
    ),
    topology: Optional[Path] = typer.Option(
        None, "--topology", "-t",
        help="Topology file (PDB, GRO, …).",
        rich_help_panel=PANEL_GENERAL,
    ),
    forcefield: Optional[Path] = typer.Option(
        None, "--forcefield", "-ff",
        help=(
            "Forcefield folder — either the original (will be parsed via "
            "ff_reader) or the parsed copy containing final_ff_residues.csv."
        ),
        rich_help_panel=PANEL_GENERAL,
    ),
    outpath: Optional[Path] = typer.Option(
        None, "--outpath", "-o",
        help="Root output folder. Creates Topology_recognition/ inside.",
        rich_help_panel=PANEL_GENERAL,
    ),
    max_to_check: Optional[int] = typer.Option(
        None, "--max-to-check",
        help=(
            "Maximum number of forcefield entries to inspect per residue. "
            "-1 = all (sorted by name similarity)."
        ),
        rich_help_panel=PANEL_TOPOPTS,
    ),
    resume: bool = typer.Option(
        False, "--resume",
        help="Skip if topology recognition outputs already exist.",
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
    Match each residue in the topology to a forcefield entry.

    Produces (under `<outpath>/Topology_recognition/`):
    `initial_topology_info.csv`, `different_residue_topologies.csv`,
    and `complete_residue_topologies.csv` (or `unmatching_…` when partial).
    """
    overrides = collect_overrides({
        "topology":     topology,
        "forcefield":   forcefield,
        "outpath":      outpath,
        "max_to_check": max_to_check,
    }, stage_key="topology_recognition")

    cfg = load_or_die(config, overrides)

    require_path("topology",   cfg.topology,   must_exist=True)
    require_path("forcefield", cfg.forcefield, must_exist=True)
    require_path("outpath",    cfg.outpath,    must_exist=False)

    topo_out = cfg.outpath / "Topology_recognition"
    print_banner("topology", subtitle=cfg.topology.name)
    print_summary("Inputs", {
        "topology":     cfg.topology,
        "forcefield":   cfg.forcefield,
        "outpath":      topo_out,
        "max_to_check": cfg.topology_recognition.max_to_check,
    })

    if dry_run:
        print_info("Dry run — nothing executed.")
        raise typer.Exit()

    marker = topo_out / "complete_residue_topologies.csv"
    if resume and warn_skip(marker, "topology", overwrite=False):
        raise typer.Exit()
    if marker.exists() and not overwrite and not resume:
        if warn_skip(marker, "topology", overwrite):
            raise typer.Exit()

    ensure_dir(topo_out)

    # Resolve forcefield path: prefer the parsed copy if present
    ff_path = _resolve_ff_path(cfg.forcefield, cfg.outpath)

    try:
        from ..services.topology_recognition import (  # type: ignore
            Topology_recognition, TopologyAssociationError,
        )
    except ImportError as e:
        print_error(
            "Could not import `mypsns.services.topology_recognition.Topology_recognition`. "
            f"Underlying error: {e}"
        )
        raise typer.Exit(code=3)

    try:
        with stage_log_file(topo_out / 'topology_recognition.log'):

            tr = Topology_recognition(
                topology=str(cfg.topology),
                outpath=topo_out,
                ff_path=ff_path,
                na_values=cfg.topology_recognition.na_values,
            )
            tr.get_restypes()
            tr.get_ff_parameters(
                maxToCheck=cfg.topology_recognition.max_to_check,
            )
    except TopologyAssociationError:
        # Expected halt, not a crash: get_ff_parameters already logged the full
        # message (unmatched residues + the accept-topology commands), so no
        # traceback that would bury it.
        print_error("topology halted: some residue topologies are unmatched (see above).")
        raise typer.Exit(code=1)
    except Exception as e:
        log.exception("topology recognition failed")
        print_error(f"topology failed: {e}")
        raise typer.Exit(code=1)

    dump_run_config(cfg, cfg.outpath)
    print_success(f"topology recognition complete → {topo_out}")


def _resolve_ff_path(ff_arg: Path, outpath: Path) -> Path:
    """
    Prefer a parsed forcefield copy in <outpath>/<ff_name>/ when present;
    otherwise fall back to the path the user gave.
    """
    candidate = outpath / 'forcefield'
    if (candidate / "final_ff_residues.csv").exists():
        return candidate
    return ff_arg