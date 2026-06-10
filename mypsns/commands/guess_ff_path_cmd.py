"""`mypsns ff-reader` — parse a GROMACS-style forcefield."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, List
import typer

from ..logging_utils import (
    print_banner, print_error, print_success, print_summary,
    stage_log_file,
)

log = logging.getLogger(__name__)

PANEL_GENERAL = "General"
PANEL_FFGUESS  = "Search options"

def guess_ff_path(
    parent_dir: Optional[str] = typer.Option(
        None, "--parent", "-parent",
        help="Parent folder of possible forcefield subfolders.",
        rich_help_panel=PANEL_GENERAL,
    ),
    nonbonded_name: str = typer.Option(
        "ffnonbonded.itp", "--nonbonded-name",
        help="Name of the .itp file with nonbonded parameters.",
        rich_help_panel=PANEL_FFGUESS,
    ),
    recursive: bool = typer.Option(
        False, "--recursive-search", "-recursive",
        help="Look recursively in subdirectories up to max_depth (False).",
        rich_help_panel=PANEL_FFGUESS,
    ),
    max_depth: int = typer.Option(
        10, "--maximum-depth", "-maxd",
        help="Maximum depth for recursive search.",
        rich_help_panel=PANEL_FFGUESS,
    ),
) -> None:
    """Scan a directory tree for GROMACS forcefield layouts and print matches."""
    if parent_dir is None:
        print_error("--parent / -parent is required.")
        raise typer.Exit(code=2)

    pdir = Path(parent_dir).expanduser().resolve()
    if not pdir.is_dir():
        print_error(f"{pdir} does not exist or is not a directory.")
        raise typer.Exit(code=2)

    print_banner("guess-ff-path", subtitle=str(pdir))
    print_summary("Inputs", {
        "Parent folder":   str(pdir),
        "Nonbonded file":  nonbonded_name,
        "Recursive":       recursive,
        "Max depth":       max_depth,
    })

    try:
        from ..services.guess_ff_path import guess_ff_path as guess_ff_path_service
    except ImportError as e:
        print_error(
            "Could not import `mypsns.services.guess_ff_path.guess_ff_path`. "
            f"Underlying error: {e}"
        )
        raise typer.Exit(code=3)

    try:
        with stage_log_file(pdir / "guess_ff_path.log"):
            paths: List[Path] = guess_ff_path_service(
                parent_dir=str(pdir),
                recursive=recursive,
                max_depth=max_depth,
                nonbonded_name=nonbonded_name,
            )
    except Exception as e:
        log.exception("guess-ff-path failed")
        print_error(f"guess-ff-path failed: {e}")
        raise typer.Exit(code=1)

    if not paths:
        print_success(f"guess-ff-path complete → no forcefield layouts under {pdir}")
        return

    str_paths = "\n".join(str(p.resolve()) for p in paths)
    print_success(f"guess-ff-path complete → forcefield directories identified:\n{str_paths}")
