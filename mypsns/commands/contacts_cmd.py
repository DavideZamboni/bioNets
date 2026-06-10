"""`mypsns contacts` — compute MD-trajectory contact statistics."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional

import typer

from ..logging_utils import (
    print_banner, print_error, print_info, print_success, print_summary,
)
from ._helpers import (
    collect_overrides, dump_run_config, ensure_dir, load_or_die,
    require_path, warn_skip,
)


log = logging.getLogger(__name__)

PANEL_GENERAL = "General"
PANEL_PARA    = "Parallelization"
PANEL_WINDOW  = "Trajectory window"
PANEL_SELECT  = "Selectors"
PANEL_CUTOFFS = "Geometric cutoffs"
PANEL_BEHAVE  = "Behaviour"


def contacts(
    config: Optional[Path] = typer.Option(
        None, "--config", help="YAML config file.", rich_help_panel=PANEL_GENERAL,
    ),
    topology: Optional[Path] = typer.Option(
        None, "--topology", "-t",
        help="Topology file (same as previous stages).",
        rich_help_panel=PANEL_GENERAL,
    ),
    trajectory: Optional[Path] = typer.Option(
        None, "--trajectory", "-x",
        help="MD trajectory file (xtc/dcd/trr/…).",
        rich_help_panel=PANEL_GENERAL,
    ),
    outpath: Optional[Path] = typer.Option(
        None, "--outpath", "-o",
        help="Root output folder containing Interactors/.",
        rich_help_panel=PANEL_GENERAL,
    ),
    nc: Optional[int] = typer.Option(
        None, "--nc",
        help="Number of parallel workers.",
        rich_help_panel=PANEL_PARA,
    ),
    ram_thresh: Optional[int] = typer.Option(
        None, "--ram-thresh",
        help="RAM-usage %% at which HBond batch growth freezes.",
        rich_help_panel=PANEL_PARA,
    ),
    nframes: Optional[int] = typer.Option(
        None, "--nframes",
        help="Initial HBond frame batch size.",
        rich_help_panel=PANEL_PARA,
    ),
    beg_frm: Optional[int] = typer.Option(
        None, "--beg-frm",
        help="First frame to analyse.",
        rich_help_panel=PANEL_WINDOW,
    ),
    end_frm: Optional[int] = typer.Option(
        None, "--end-frm",
        help="Last frame to analyse (-1 = end of trajectory).",
        rich_help_panel=PANEL_WINDOW,
    ),
    targets: Optional[List[str]] = typer.Option(
        None, "--target",
        help=(
            "Contact type(s) to compute. Repeat for multiple. "
            "Choices: HBonds, Salt_bridges, HCs."
        ),
        rich_help_panel=PANEL_SELECT,
    ),
    interactors_type: Optional[List[str]] = typer.Option(
        None, "--interactors-type",
        help="Interactor schema(s) to consult. Repeat for multiple.",
        rich_help_panel=PANEL_SELECT,
    ),
    subject_chains: Optional[List[str]] = typer.Option(
        None, "--subject-chain",
        help="Subject chain(s). Repeat for multiple. Auto-detect if omitted.",
        rich_help_panel=PANEL_SELECT,
    ),
    surrounding_chains: Optional[List[str]] = typer.Option(
        None, "--surrounding-chain",
        help="Surrounding chain(s). Repeat for multiple. Auto-detect if omitted.",
        rich_help_panel=PANEL_SELECT,
    ),
    d_a_dist: Optional[float] = typer.Option(
        None, "--d-a-dist",
        help="Donor-acceptor distance cutoff [Å].",
        rich_help_panel=PANEL_CUTOFFS,
    ),
    angle: Optional[float] = typer.Option(
        None, "--angle",
        help="Donor-H-acceptor angle cutoff [°].",
        rich_help_panel=PANEL_CUTOFFS,
    ),
    sb_cutoff: Optional[float] = typer.Option(
        None, "--sb-cutoff",
        help="Salt-bridge atom-pair cutoff [Å].",
        rich_help_panel=PANEL_CUTOFFS,
    ),
    hc_cutoff: Optional[float] = typer.Option(
        None, "--hc-cutoff",
        help="Hydrophobic-cluster CoM cutoff [Å].",
        rich_help_panel=PANEL_CUTOFFS,
    ),
    matrix_type: Optional[List[str]] = typer.Option(
        None, "--matrix-type",
        help="Type of matrix to compute (persistency, energy). Repeat for multiple.",
        rich_help_panel=PANEL_SELECT,
    ),
    energy_use_abs: Optional[bool] = typer.Option(
        None, "--energy-use-abs",
        help="Whether to retain repulsive forces present in energy network analysis",
        rich_help_panel=PANEL_SELECT,
    ),
    thr_step: Optional[float] = typer.Option(
        None, "--thr-step",
        help="Step of % persistency or energy used when computing network statistics.",
        rich_help_panel=PANEL_SELECT,
    ),
    resume: bool = typer.Option(
        False, "--resume",
        help="Skip contact computation if final outputs already exist.",
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
    Compute contact statistics across the trajectory.

    Reads `<outpath>/Interactors/` and writes per-frame partial CSVs,
    merged CSVs, final per-residue-pair tables, adjacency matrices, and
    persistency-threshold analyses to `<outpath>/Contacts/`.
    """
    overrides = collect_overrides({
        "topology":           topology,
        "trajectory":         trajectory,
        "outpath":            outpath,
        "nc":                 nc,
        "ram_thresh":         ram_thresh,
        "nframes":            nframes,
        "beg_frm":            beg_frm,
        "end_frm":            end_frm,
        "targets":            targets,
        "interactors_type":   interactors_type,
        "subject_chains":     subject_chains,
        "surrounding_chains": surrounding_chains,
        "d_a_dist":           d_a_dist,
        "angle":              angle,
        "sb_cutoff":          sb_cutoff,
        "hc_cutoff":          hc_cutoff,
        "matrix_type":        matrix_type,
        "energy_use_abs":     energy_use_abs,
        "thr_step":           thr_step,
    }, stage_key="contacts")

    cfg = load_or_die(config, overrides)

    require_path("topology",   cfg.topology,   must_exist=True)
    require_path("trajectory", cfg.trajectory, must_exist=True)
    require_path("outpath",    cfg.outpath,    must_exist=True)

    cnt_out = cfg.outpath / "Adjacency"
    int_in  = cfg.outpath / "Interactors"
    if not int_in.exists():
        print_error(
            f"Interactors folder not found in {cfg.outpath}. "
            "Run `mypsns interactors` first."
        )
        raise typer.Exit(code=2)

    print_banner("contacts", subtitle=cfg.trajectory.name)
    print_summary("Inputs", {
        "topology":    cfg.topology,
        "trajectory":  cfg.trajectory,
        "interactors": int_in,
        "outpath":     cnt_out,
    })
    print_summary("Parallelization", {
        "nc":         cfg.contacts.nc,
        "ram_thresh": cfg.contacts.ram_thresh,
        "nframes":    cfg.contacts.nframes,
    })
    print_summary("Selectors", {
        "targets":            cfg.contacts.targets,
        "interactors_type":   cfg.contacts.interactors_type,
        "subject_chains":     cfg.contacts.subject_chains,
        "surrounding_chains": cfg.contacts.surrounding_chains,
        "frame_window":       f"{cfg.contacts.beg_frm}..{cfg.contacts.end_frm}",
    })
    print_summary("Cutoffs", {
        "d_a_dist":  f"{cfg.contacts.d_a_dist} Å",
        "angle":     f"{cfg.contacts.angle}°",
        "sb_cutoff": f"{cfg.contacts.sb_cutoff} Å",
        "hc_cutoff": f"{cfg.contacts.hc_cutoff} Å",
    })

    print_summary("Network analysis parameters", {
        "matrix_type":  f"{cfg.contacts.matrix_type} Å",
        "energy_use_abs":     f"{cfg.contacts.energy_use_abs}°",
        "thr_step": f"{cfg.contacts.thr_step} Å",
    })

    if dry_run:
        print_info("Dry run — nothing executed.")
        raise typer.Exit()

    # Resume / overwrite logic — check for the canonical final marker
    marker = cnt_out / 'graph_report' / 'plot'   # final preliminary data output
    if resume and marker.exists() and not overwrite:
        print_info(f"contacts: final output found at {marker} — skipping (resume).")
        raise typer.Exit()
    if warn_skip(marker, "contacts", overwrite):
        raise typer.Exit()

    ensure_dir(cnt_out)

    try:
        from ..services.contacts import Contacts  # type: ignore
    except ImportError as e:
        print_error(
            "Could not import `mypsns.services.contacts.Contacts`. "
            f"Underlying error: {e}"
        )
        raise typer.Exit(code=3)

    try:
        # Runtime params go to the CONSTRUCTOR (**kwargs). Leaving
        # subject_chains / surrounding_chains as None lets Contacts auto-infer
        # them (subject = polymeric chains, surrounding = all chains).
        c = Contacts(
            topology=str(cfg.topology),
            trj_file=str(cfg.trajectory),
            outpath=cfg.outpath,
            nc=cfg.contacts.nc,
            ram_thresh=cfg.contacts.ram_thresh,
            nframes=cfg.contacts.nframes,
            beg_frm=cfg.contacts.beg_frm,
            end_frm=cfg.contacts.end_frm,
            interactors_type=cfg.contacts.interactors_type,
            targets=cfg.contacts.targets,
            subject_chains=cfg.contacts.subject_chains,
            surrounding_chains=cfg.contacts.surrounding_chains,
            d_a_dist=cfg.contacts.d_a_dist,
            angle=cfg.contacts.angle,
            sb_cutoff=cfg.contacts.sb_cutoff,
            hc_cutoff=cfg.contacts.hc_cutoff,
            na_values=cfg.contacts.na_values,
            matrix_type=cfg.contacts.matrix_type,
            energy_use_abs=cfg.contacts.energy_use_abs,
            thr_step=cfg.contacts.thr_step,
        )
        c.run()
        c.get_matrices()
    except Exception as e:
        log.exception("contacts failed")
        print_error(f"contacts failed: {e}")
        raise typer.Exit(code=1)

    dump_run_config(cfg, cfg.outpath)
    print_success(f"contacts complete → {cnt_out}")