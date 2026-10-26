"""`mypsns run` — orchestrate the full pipeline."""
from __future__ import annotations

import logging
import timeit
from pathlib import Path
from typing import List, Optional
import pandas as pd
import typer

from ..logging_utils import (
    print_banner, print_env_info, print_error, print_info, print_success,
    print_summary, stage_log_file,
)

from ._helpers import (
    contacts_done, dump_run_config, ensure_dir,
    load_or_die, require_path,
)


log = logging.getLogger(__name__)

PANEL_GENERAL = "General"
PANEL_PIPE    = "Pipeline control"
PANEL_QUICK   = "Common overrides"


def run_pipeline(
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
    trajectory: Optional[Path] = typer.Option(
        None, "--trajectory", "-x",
        help="MD trajectory file.",
        rich_help_panel=PANEL_GENERAL,
    ),
    forcefield: Optional[Path] = typer.Option(
        None, "--forcefield", "-ff",
        help="Forcefield folder.",
        rich_help_panel=PANEL_GENERAL,
    ),
    outpath: Optional[Path] = typer.Option(
        None, "--outpath", "-o",
        help="Root output folder for the run.",
        rich_help_panel=PANEL_GENERAL,
    ),
    nc: Optional[int] = typer.Option(
        None, "--nc",
        help="Number of parallel workers (contacts stage).",
        rich_help_panel=PANEL_QUICK,
    ),
    targets: Optional[List[str]] = typer.Option(
        None, "--target",
        help="Contact types to compute. Repeat for multiple.",
        rich_help_panel=PANEL_QUICK,
    ),
    skip: Optional[List[str]] = typer.Option(
        None, "--skip",
        help=(
            "Stage(s) to skip. Repeat for multiple. "
            "Choices: ff-reader, topology, interactors, contacts."
        ),
        rich_help_panel=PANEL_PIPE,
    ),
    only: Optional[List[str]] = typer.Option(
        None, "--only",
        help="Run only the listed stage(s). Repeat for multiple.",
        rich_help_panel=PANEL_PIPE,
    ),
    resume: bool = typer.Option(
        False, "--resume",
        help="Skip stages whose outputs already exist.",
        rich_help_panel=PANEL_PIPE,
    ),
    overwrite: bool = typer.Option(
        False, "--overwrite",
        help="Force recomputation of every stage.",
        rich_help_panel=PANEL_PIPE,
    ),
    dry_run: bool = typer.Option(
        False, "--dry-run",
        help="Validate everything; do not run any stage.",
        rich_help_panel=PANEL_PIPE,
    ),
) -> None:
    """
    Run the full pipeline:
    [yellow]ff-reader → topology → interactors → contacts[/yellow].
    """
    # Quick CLI overrides go into the right nested namespaces.
    overrides: dict = {}
    if topology   is not None: overrides["topology"]   = topology
    if trajectory is not None: overrides["trajectory"] = trajectory
    if forcefield is not None: overrides["forcefield"] = forcefield
    if outpath    is not None: overrides["outpath"]    = outpath
    if nc         is not None: overrides.setdefault("contacts", {})["nc"]       = nc
    if targets:                 overrides.setdefault("contacts", {})["targets"] = targets

    cfg = load_or_die(config, overrides)

    require_path("topology",   cfg.topology,   must_exist=True)
    require_path("forcefield", cfg.forcefield, must_exist=True)
    require_path("outpath",    cfg.outpath,    must_exist=False)
    if cfg.trajectory is not None:
        require_path("trajectory", cfg.trajectory, must_exist=True)

    stages_all   = ["ff-reader", "topology", "interactors", "contacts"]
    stages_to_run = _resolve_stage_set(stages_all, only=only, skip=skip)

    if "contacts" in stages_to_run and cfg.trajectory is None:
        print_error("contacts stage requires --trajectory (or a trajectory in --config).")
        raise typer.Exit(code=2)

    print_env_info()
    print_banner("myPSNs run", subtitle=str(cfg.outpath))
    print_summary("Plan", {
        "topology":   cfg.topology,
        "trajectory": cfg.trajectory,
        "forcefield": cfg.forcefield,
        "outpath":    cfg.outpath,
        "stages":     stages_to_run,
        "resume":     resume,
        "overwrite":  overwrite,
    })

    if dry_run:
        print_info("Dry run — nothing executed.")
        raise typer.Exit()

    ensure_dir(cfg.outpath)
    # the four stages `run` chains; network/compare_variants keep what they recorded
    dump_run_config(cfg, cfg.outpath, sections=[
        "ff_reader", "topology_recognition", "interactors", "contacts"])

    t0 = timeit.default_timer()

    # imported here, like the stage services, so `mypsns --help` stays fast
    from ..services.topology_recognition import TopologyAssociationError  # type: ignore

    for stage in stages_to_run:
        try:
            _run_stage(stage, cfg, resume=resume, overwrite=overwrite)
        except typer.Exit as e:
            # Sub-stage signalled exit (success or skip). Continue iff code 0.
            if e.exit_code not in (0, None):
                raise
        except TopologyAssociationError:
            # Expected halt, not a crash: the message (unmatched residues +
            # accept-topology commands) is already logged; a traceback would bury it.
            print_error(f"Stage {stage} halted: some residue topologies are unmatched (see above).")
            raise typer.Exit(code=1)
        except Exception as e:
            log.exception("Pipeline stage %s failed", stage)
            print_error(f"Stage {stage} failed: {e}")
            raise typer.Exit(code=1)

    dt = timeit.default_timer() - t0
    print_success(f"Pipeline complete in {dt:.1f}s → {cfg.outpath}")


# ─────────────────────────────────────────────────────────────────────────────
# Internals
# ─────────────────────────────────────────────────────────────────────────────

def _resolve_stage_set(
    all_stages: list[str],
    *,
    only: Optional[List[str]],
    skip: Optional[List[str]],
) -> list[str]:
    """Determine which stages to run respecting --only / --skip."""
    if only:
        return [s for s in all_stages if s in set(only)]
    skip_set = set(skip or [])
    return [s for s in all_stages if s not in skip_set]


def _run_stage(
    stage: str,
    cfg,
    *,
    resume: bool,
    overwrite: bool,
) -> None:
    """Dispatch a single stage."""
    if stage == "ff-reader":
        _run_ff_reader(cfg, resume=resume, overwrite=overwrite)
    elif stage == "topology":
        _run_topology(cfg, resume=resume, overwrite=overwrite)
    elif stage == "interactors":
        _run_interactors(cfg, resume=resume, overwrite=overwrite)
    elif stage == "contacts":
        _run_contacts(cfg, resume=resume, overwrite=overwrite)
    else:
        print_error(f"Unknown stage: {stage}")
        raise typer.Exit(code=2)


def _run_ff_reader(cfg, *, resume: bool, overwrite: bool) -> None:
    ff_out = cfg.outpath
    marker = ff_out / "final_ff_residues.csv"
    print_banner("Stage 1/4 — ff-reader")

    if marker.exists() and (resume or not overwrite):
        print_info(f"ff-reader: skipped ({marker} already present).")
        return

    ensure_dir(ff_out)
    from ..services.ff_reader import ff_reader
    if not ((ff_out / 'forcefield' / 'final_ff_ters.csv').exists() and (ff_out / 'forcefield' / 'final_ff_residues.csv').exists()):
        with stage_log_file(ff_out / 'forcefield' / 'ff_reader.log'):
            ff_reader(
                inpath=cfg.forcefield,
                outpath=ff_out,
                na_values=cfg.ff_reader.na_values,
                nonbonded_name=cfg.ff_reader.nonbonded_name,
                additional_residues=cfg.ff_reader.additional_residues,
            )
        print_success(f"ff-reader → {ff_out}")
    else:
        print_success(f"forcefield-derived files already present in {ff_out}, skipping")

def _run_topology(cfg, *, resume: bool, overwrite: bool) -> None:
    topo_out = cfg.outpath / "Topology_recognition"
    marker   = topo_out / "complete_residue_topologies.csv"
    print_banner("Stage 2/4 — topology")

    if marker.exists() and (resume or not overwrite):
        print_info(f"topology: skipped ({marker} already present).")
        return

    ensure_dir(topo_out)
    from ..services.topology_recognition import Topology_recognition  # type: ignore
    candidate = cfg.outpath / 'forcefield'
    #ff_path = candidate if (candidate / "final_ff_residues.csv").exists() else cfg.forcefield
    
    with stage_log_file(topo_out / 'topology_recognition.log'):
        tr = Topology_recognition(
            topology=str(cfg.topology),
            outpath=topo_out,
            ff_path=candidate,
            na_values=cfg.topology_recognition.na_values,
        )
        tr.get_restypes()
        tr.get_ff_parameters(
            maxToCheck=cfg.topology_recognition.max_to_check,
        )
        print_success(f"topology → {topo_out}")


def _run_interactors(cfg, *, resume: bool, overwrite: bool) -> None:
    int_out = cfg.outpath / "Interactors"
    skip = True
    gen_out  = int_out / "found_topologies"
    print_banner("Stage 3/4 — interactors")

    if gen_out.exists():
        already_foud = [f for f in gen_out.iterdir() if f.is_file() and f.suffix == '.csv']
        for f in already_foud:
            tmp = pd.read_csv(f, header=0, index_col=False,
                              na_values=cfg.interactors.na_values, keep_default_na=False)
            cols = tmp.columns.to_list()
            if not ('named_interactors' in cols and 'q_interactors' in cols):
                skip = False
    else:
        skip = False

    if skip and (resume or not overwrite):
        print_info(f"interactors: skipped (interactors already assigned).")
        return

    ensure_dir(int_out)
    from ..services.interactors import Interactors  # type: ignore
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
            min_num_heavy_ats=cfg.interactors.min_num_heavy_ats,
            min_bonds_to_charged=cfg.interactors.min_bonds_to_charged,
            min_bonds_to_polar=cfg.interactors.min_bonds_to_polar,
            # same set as `mypsns interactors`: without these two, `run` typed
            # named_interactors with the built-in rules and 0.8 e gate whatever
            # the config said
            named_rules=cfg.interactors.named_interactors_rules,
            min_residue_charge=cfg.interactors.min_residue_charge,
        )
        print_success(f"interactors → {int_out}")

def _run_contacts(cfg, *, resume: bool, overwrite: bool) -> None:
    cnt_out = cfg.outpath / "Adjacency"
    print_banner("Stage 4/4 — contacts")
 
    if contacts_done(cfg) and (resume or not overwrite):
        print_info(f"contacts: skipped (adjacency output already present in {cnt_out}).")
        return
 
    from ..services.contacts import Contacts  # type: ignore
    # Runtime params go to the CONSTRUCTOR (**kwargs); None chain values are
    # auto-inferred inside Contacts (subject=polymers, surrounding=all).
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
        d_h_cutoff=cfg.contacts.d_h_cutoff,
        sb_cutoff=cfg.contacts.sb_cutoff,
        hc_cutoff=cfg.contacts.hc_cutoff,
        na_values=cfg.contacts.na_values,
        matrix_type=cfg.contacts.matrix_type,
        energy_use_abs=cfg.contacts.energy_use_abs,
        thr_step=cfg.contacts.thr_step
    )
    c.run()
    c.get_matrices()
    print_success(f"contacts → {cnt_out}")
 