"""`mypsns interactors` — assign interaction roles to topology atoms."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional

import typer

from ..logging_utils import (
    print_banner, print_error, print_info, print_success, print_summary, stage_log_file,
)
from ._helpers import (
    collect_overrides, dump_run_config, ensure_dir, load_or_die,
    require_path, warn_skip,
)


log = logging.getLogger(__name__)

PANEL_GENERAL = "General"
PANEL_MODES   = "Modes"
PANEL_HBOND   = "HBond tuning"
PANEL_SB      = "Salt bridge tuning"
PANEL_HC      = "Hydrophobic cluster tuning"
PANEL_BEHAVE  = "Behaviour"


def interactors(
    config: Optional[Path] = typer.Option(
        None, "--config", help="YAML config file.", rich_help_panel=PANEL_GENERAL,
    ),
    topology: Optional[Path] = typer.Option(
        None, "--topology", "-t",
        help="Topology file (same as topology stage).",
        rich_help_panel=PANEL_GENERAL,
    ),
    outpath: Optional[Path] = typer.Option(
        None, "--outpath", "-o",
        help="Root output folder containing Topology_recognition/.",
        rich_help_panel=PANEL_GENERAL,
    ),
    modes: Optional[List[str]] = typer.Option(
        None, "--mode", "-m",
        help=(
            "Interactor-typing scheme(s) to run. Repeat for multiple. "
            "Choices: q_interactors, named_interactors."
        ),
        rich_help_panel=PANEL_MODES,
    ),
    max_DA_charge: Optional[float] = typer.Option(
        None, "--max-da-charge",
        help="Max charge for HBond donor/acceptor.",
        rich_help_panel=PANEL_HBOND,
    ),
    min_H_charge: Optional[float] = typer.Option(
        None, "--min-h-charge",
        help="Min H charge for HBond participation.",
        rich_help_panel=PANEL_HBOND,
    ),
    min_abs_charge: Optional[float] = typer.Option(
        None, "--min-abs-charge",
        help="Min |local charge| for Acidic/Basic classification.",
        rich_help_panel=PANEL_SB,
    ),
    max_abs_charge: Optional[float] = typer.Option(
        None, "--max-abs-charge",
        help="Max |charge| for HC eligibility.",
        rich_help_panel=PANEL_HC,
    ),
    max_num_heavy_ats: Optional[int] = typer.Option(
        None, "--max-num-heavy-ats",
        help="Max heavy atoms per HC sub-cluster.",
        rich_help_panel=PANEL_HC,
    ),
    min_bonds_to_charged: Optional[int] = typer.Option(
        None, "--min-bonds-to-charged",
        help="Min bond-distance from any charged atom for HC eligibility.",
        rich_help_panel=PANEL_HC,
    ),
    min_bonds_to_polar: Optional[int] = typer.Option(
        None, "--min-bonds-to-polar",
        help="Min bond-distance from any polar atom for HC eligibility.",
        rich_help_panel=PANEL_HC,
    ),
    resume: bool = typer.Option(
        False, "--resume",
        help="Skip if interactor outputs already exist.",
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
    Assign interaction roles to each atom in the topology.

    Reads `<outpath>/Topology_recognition/` and writes per-variant CSVs to
    `<outpath>/Interactors/found_topologies/`, plus `residue_index.csv`
    and `atom_index_map.csv`.
    """
    overrides = collect_overrides({
        "topology":              topology,
        "outpath":               outpath,
        "modes":                 modes,
        "max_DA_charge":         max_DA_charge,
        "min_H_charge":          min_H_charge,
        "min_abs_charge":        min_abs_charge,
        "max_abs_charge":        max_abs_charge,
        "max_num_heavy_ats":     max_num_heavy_ats,
        "min_bonds_to_charged":  min_bonds_to_charged,
        "min_bonds_to_polar":    min_bonds_to_polar,
    }, stage_key="interactors")

    cfg = load_or_die(config, overrides)

    require_path("topology", cfg.topology, must_exist=True)
    require_path("outpath",  cfg.outpath,  must_exist=False)

    int_out = cfg.outpath / "Interactors"
    topo_in = cfg.outpath / "Topology_recognition"
    if not topo_in.exists():
        print_error(
            f"Topology_recognition folder not found in {cfg.outpath}. "
            "Run `mypsns topology` first."
        )
        raise typer.Exit(code=2)

    print_banner("interactors", subtitle=", ".join(cfg.interactors.modes))
    print_summary("Inputs", {
        "topology":   cfg.topology,
        "input dir":  topo_in,
        "output dir": int_out,
        "modes":      cfg.interactors.modes,
    })
    print_summary("Tuning", {
        "max_DA_charge":        cfg.interactors.max_DA_charge,
        "min_H_charge":         cfg.interactors.min_H_charge,
        "min_abs_charge":       cfg.interactors.min_abs_charge,
        "max_abs_charge":       cfg.interactors.max_abs_charge,
        "max_num_heavy_ats":    cfg.interactors.max_num_heavy_ats,
        "min_bonds_to_charged": cfg.interactors.min_bonds_to_charged,
        "min_bonds_to_polar":   cfg.interactors.min_bonds_to_polar,
    })

    if dry_run:
        print_info("Dry run — nothing executed.")
        raise typer.Exit()

    marker = int_out / "residue_index.csv"
    if resume and stage_done_quietly(marker):
        print_info(f"interactors: outputs found at {marker} — skipping (resume).")
        raise typer.Exit()
    if warn_skip(marker, "interactors", overwrite):
        raise typer.Exit()

    ensure_dir(int_out)

    try:
        from ..services.interactors import Interactors  # type: ignore
    except ImportError as e:
        print_error(
            "Could not import `mypsns.services.interactors.Interactors`. "
            f"Underlying error: {e}"
        )
        raise typer.Exit(code=3)

    try:
        with stage_log_file(int_out / 'interactors_recognition.log'):
            I = Interactors(
                topology=str(cfg.topology),
                general_outpath=cfg.outpath,
                na_values=cfg.interactors.na_values,
            )
            I.get_complete_topology()
            I.add_inttypes(
                max_DA_charge=cfg.interactors.max_DA_charge,
                min_H_charge=cfg.interactors.min_H_charge,
                DA_masses=cfg.interactors.DA_masses,
                min_abs_charge=cfg.interactors.min_abs_charge,
                max_abs_charge=cfg.interactors.max_abs_charge,
                max_num_heavy_ats=cfg.interactors.max_num_heavy_ats,
                min_bonds_to_charged=cfg.interactors.min_bonds_to_charged,
                min_bonds_to_polar=cfg.interactors.min_bonds_to_polar,
                named_rules=cfg.interactors.named_interactors_rules,
                min_residue_charge=cfg.interactors.min_residue_charge,
            )
    except Exception as e:
        log.exception("interactors failed")
        print_error(f"interactors failed: {e}")
        raise typer.Exit(code=1)

    dump_run_config(cfg, cfg.outpath)
    print_success(f"interactors complete → {int_out}")


def stage_done_quietly(marker: Path) -> bool:
    return Path(marker).exists()