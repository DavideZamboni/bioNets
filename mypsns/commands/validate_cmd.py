"""`mypsns validate` — validate a config and its referenced files."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import typer

from ..config import ConfigError, build_pipeline_config
from ..logging_utils import (
    print_banner, print_error, print_info, print_kv_table, print_success,
    print_warning,
)


log = logging.getLogger(__name__)

def validate(
    config: Path = typer.Option(
        ..., "--config",
        help="YAML config file to validate.",
    ),
    deep: bool = typer.Option(
        False, "--deep",
        help="Also try to open the topology and trajectory with MDAnalysis.",
    ),
) -> None:
    """
    Validate a config file and (optionally) the inputs it points to.

    Without --deep: validates schema + checks that referenced paths exist.
    With --deep: additionally tries to load the topology/trajectory and
    reports atom/residue/frame counts.
    """
    print_banner("validate", subtitle=str(config))

    rows: list[tuple[str, str, str]] = []

    # ── Schema ──────────────────────────────────────────────────────────
    try:
        cfg = build_pipeline_config(config_path=config)
    except ConfigError as e:
        rows.append(("Config schema", "[red]✗[/red]", str(e).replace("\n", " | ")))
        print_kv_table("Validation", rows)
        raise typer.Exit(code=2)
    rows.append(("Config schema", "[green]✓[/green]", "all fields valid"))

    # ── Path existence ─────────────────────────────────────────────────
    rows.append(_check_path("topology", cfg.topology, must_exist=True))
    rows.append(_check_path("forcefield", cfg.forcefield, must_exist=True))
    rows.append(_check_path("trajectory", cfg.trajectory, must_exist=True))
    rows.append(_check_path("outpath", cfg.outpath, must_exist=False))

    # ── Forcefield internals ───────────────────────────────────────────
    if cfg.forcefield and cfg.forcefield.exists():
        rows.append(_check_ff(cfg.forcefield, cfg.ff_reader.nonbonded_name))

    # ── Output directory writability ───────────────────────────────────
    rows.append(_check_writable(cfg.outpath))

    # ── Deep checks ────────────────────────────────────────────────────
    if deep:
        rows.extend(_deep_checks(cfg))

    print_kv_table("Validation", rows)

    failures = [r for r in rows if "✗" in r[1]]
    if failures:
        print_error(f"{len(failures)} validation check(s) failed.")
        raise typer.Exit(code=1)

    print_success("All validation checks passed.")


# ─────────────────────────────────────────────────────────────────────────────
# Individual checks
# ─────────────────────────────────────────────────────────────────────────────

def _check_path(label: str, p: Optional[Path], *, must_exist: bool) -> tuple[str, str, str]:
    if p is None:
        return (label, "[yellow]—[/yellow]", "not set")
    p = Path(p).expanduser()
    if must_exist and not p.exists():
        return (label, "[red]✗[/red]", f"missing: {p}")
    return (label, "[green]✓[/green]", str(p))


def _check_writable(outpath: Path) -> tuple[str, str, str]:
    try:
        outpath = Path(outpath).expanduser()
        outpath.mkdir(parents=True, exist_ok=True)
        probe = outpath / ".mypsns_validate_probe"
        probe.write_text("ok")
        probe.unlink()
        return ("outpath writable", "[green]✓[/green]", str(outpath))
    except Exception as e:
        return ("outpath writable", "[red]✗[/red]", f"{outpath}: {e}")


def _check_ff(ff: Path, nonbonded_name: str) -> tuple[str, str, str]:
    rtps = list(ff.glob("*.rtp"))
    itp  = ff / nonbonded_name

    if not rtps:
        return ("forcefield content", "[red]✗[/red]", f"no .rtp files found in {ff}")
    if not itp.exists():
        return ("forcefield content", "[red]✗[/red]",
                f"nonbonded file missing: {itp.name}")
    return ("forcefield content", "[green]✓[/green]",
            f"{len(rtps)} .rtp file(s); {itp.name} present")


def _deep_checks(cfg) -> list[tuple[str, str, str]]:
    out: list[tuple[str, str, str]] = []

    try:
        import MDAnalysis as mda  # type: ignore
    except ImportError as e:
        out.append(("MDAnalysis", "[red]✗[/red]", f"import failed: {e}"))
        return out

    # Topology
    try:
        u = mda.Universe(str(cfg.topology))
        out.append((
            "topology loadable",
            "[green]✓[/green]",
            f"{u.atoms.n_atoms} atoms, {u.residues.n_residues} residues",
        ))
    except Exception as e:
        out.append(("topology loadable", "[red]✗[/red]", str(e)))
        return out

    # Trajectory
    if cfg.trajectory is not None:
        try:
            u = mda.Universe(str(cfg.topology), str(cfg.trajectory))
            out.append((
                "trajectory loadable",
                "[green]✓[/green]",
                f"{u.trajectory.n_frames} frames",
            ))
        except Exception as e:
            out.append(("trajectory loadable", "[red]✗[/red]", str(e)))

    return out
