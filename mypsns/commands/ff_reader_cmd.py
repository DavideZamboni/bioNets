"""`mypsns ff-reader` — parse a GROMACS-style forcefield."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import typer

from ..logging_utils import (
    print_banner, print_error, print_info, print_success, print_summary,
    stage_log_file,
)
from ._helpers import (
    collect_overrides, dump_run_config, ensure_dir,
    load_or_die, require_path, warn_skip,
)

log = logging.getLogger(__name__)

PANEL_GENERAL = "General"
PANEL_FFOPTS  = "Forcefield options"
PANEL_BEHAVE  = "Behaviour"


def ff_reader(
    config: Optional[Path] = typer.Option(
        None, "--config",
        help="YAML config file. CLI args override config values.",
        rich_help_panel=PANEL_GENERAL,
    ),
    forcefield: Optional[Path] = typer.Option(
        None, "--forcefield", "-ff",
        help="Input forcefield folder (e.g. amber99sb.ff/).",
        rich_help_panel=PANEL_GENERAL,
    ),
    outpath: Optional[Path] = typer.Option(
        None, "--outpath", "-o",
        help="Root output folder.",
        rich_help_panel=PANEL_GENERAL,
    ),
    nonbonded_name: Optional[str] = typer.Option(
        None, "--nonbonded-name",
        help="Name of the .itp file with nonbonded parameters.",
        rich_help_panel=PANEL_FFOPTS,
    ),
    no_additional_residues: bool = typer.Option(
        False, "--no-additional-residues",
        help="Skip merging extra residues from resources/.",
        rich_help_panel=PANEL_FFOPTS,
    ),
    overwrite: bool = typer.Option(
        False, "--overwrite",
        help="Force recomputation if outputs already exist.",
        rich_help_panel=PANEL_BEHAVE,
    ),
    dry_run: bool = typer.Option(
        False, "--dry-run",
        help="Validate inputs and print the plan; do not run.",
        rich_help_panel=PANEL_BEHAVE,
    ),
) -> None:
    """
    Parse a GROMACS forcefield folder.

    Reads all .rtp / .tdb files plus the nonbonded .itp file and writes
    `final_ff_residues.csv` and `final_ff_ters.csv` under
    `<outpath>/<sanitised_ff_name>/`.
    """
    overrides = collect_overrides({
        "forcefield":          forcefield,
        "outpath":             outpath,
        "nonbonded_name":      nonbonded_name,
        "additional_residues": (False if no_additional_residues else None),
    }, stage_key="ff_reader")

    cfg = load_or_die(config, overrides)

    # Validate required paths
    require_path("forcefield", cfg.forcefield, must_exist=True)
    require_path("outpath",    cfg.outpath,    must_exist=False)

    ff_out = cfg.outpath / "forcefield"

    print_banner("ff-reader", subtitle=f"{cfg.forcefield.name} → {ff_out.name}")
    print_summary("Inputs", {
        "forcefield folder":   cfg.forcefield,
        "output folder":       ff_out,
        "nonbonded file":      cfg.ff_reader.nonbonded_name,
        "additional residues": cfg.ff_reader.additional_residues,
    })

    if dry_run:
        print_info("Dry run — nothing executed.")
        raise typer.Exit()

    marker = ff_out / "final_ff_residues.csv"
    if warn_skip(marker, "ff-reader", overwrite):
        raise typer.Exit()

    ensure_dir(cfg.outpath)

    try:
        from ..services.ff_reader import ff_reader as ff_reader_service
    except ImportError as e:
        print_error(
            "Could not import `mypsns.services.ff_reader.ff_reader`. "
            f"Underlying error: {e}"
        )
        raise typer.Exit(code=3)

    try:
        with stage_log_file(ff_out / 'ff_reader.log'):
            ff_reader_service(
                inpath=str(cfg.forcefield),
                outpath=str(cfg.outpath),
                na_values=cfg.ff_reader.na_values,
                nonbonded_name=cfg.ff_reader.nonbonded_name,
                additional_residues=cfg.ff_reader.additional_residues,
            )
    except FileNotFoundError as e:
        print_error(f"Missing forcefield file: {e}")
        raise typer.Exit(code=1)
    except Exception as e:
        log.exception("ff-reader failed")
        print_error(f"ff-reader failed: {e}")
        raise typer.Exit(code=1)

    dump_run_config(cfg, cfg.outpath)
    print_success(f"ff-reader complete → {ff_out}")
