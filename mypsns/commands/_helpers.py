"""Shared helpers for command modules."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Optional

import typer

from ..config import (
    ConfigError,
    PipelineConfig,
    build_pipeline_config,
    save_pipeline_config,
)
from ..logging_utils import print_error, print_info, print_summary, print_warning


log = logging.getLogger(__name__)


def collect_overrides(
    values: dict[str, Any],
    *,
    stage_key: Optional[str] = None,
) -> dict:
    """
    Build a nested override dict from CLI option values.

    Values that are None are stripped (they correspond to unset CLI options).
    If stage_key is provided, the cleaned values are nested under that key.
    Top-level keys ("topology", "trajectory", "forcefield", "outpath") are
    extracted before nesting so they can sit at the PipelineConfig root.
    """
    top_keys = {"topology", "trajectory", "forcefield", "outpath"}
    stage_values: dict[str, Any] = {}
    top_values: dict[str, Any] = {}

    for k, v in values.items():
        if v is None:
            continue
        if k in top_keys:
            top_values[k] = v
        else:
            stage_values[k] = v

    result: dict = dict(top_values)
    if stage_values and stage_key is not None:
        result[stage_key] = stage_values
    elif stage_values:
        # No stage namespace requested → merge flat (rarely used).
        result.update(stage_values)
    return result


def load_or_die(
    config_path: Optional[Path],
    overrides: Optional[dict] = None,
) -> PipelineConfig:
    """Build the pipeline config or exit cleanly with a clear error."""
    try:
        return build_pipeline_config(
            config_path=config_path, overrides=overrides or {},
        )
    except ConfigError as e:
        print_error(str(e))
        raise typer.Exit(code=2)


def ensure_dir(path: Path) -> None:
    """Create a directory (parents inclusive) idempotently."""
    Path(path).mkdir(parents=True, exist_ok=True)


def dump_run_config(cfg: PipelineConfig, outpath: Path) -> Path:
    """Save the resolved config to <outpath>/run_config.yaml for reproducibility."""
    ensure_dir(outpath)
    dump_path = outpath / "run_config.yaml"
    save_pipeline_config(cfg, dump_path)
    return dump_path


def stage_done(marker: Path) -> bool:
    """True if the given marker file/folder exists (used by --resume)."""
    return Path(marker).exists()

def show_paths(cfg: PipelineConfig, *, include_trajectory: bool = True) -> None:
    """Print a compact paths summary."""
    data: dict[str, Any] = {
        "topology":   cfg.topology,
        "forcefield": cfg.forcefield,
        "outpath":    cfg.outpath,
    }
    if include_trajectory:
        data["trajectory"] = cfg.trajectory
    print_summary("Paths", data)


def require_path(label: str, p: Optional[Path], must_exist: bool = True) -> Path:
    """Validate a path argument; print an error and exit if missing/non-existent."""
    if p is None:
        print_error(f"--{label} is required (pass via CLI or include in --config).")
        raise typer.Exit(code=2)
    p = Path(p).expanduser()
    if must_exist and not p.exists():
        print_error(f"{label} does not exist: {p}")
        raise typer.Exit(code=2)
    return p


def warn_skip(marker: Path, stage: str, overwrite: bool) -> bool:
    """
    Return True if the stage should be skipped because output already exists
    and --overwrite was not passed.
    """
    if stage_done(marker) and not overwrite:
        print_warning(
            f"{stage}: output already present at {marker} — skipping. "
            f"Use --overwrite to force."
        )
        return True
    return False