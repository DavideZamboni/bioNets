"""
cli.py
myPSNs command-line interface.

Entry point exposed as `mypsns` via pyproject.toml console_scripts.
Subcommands live in ``mypsns.commands.*``; this module just wires the Typer
app together and configures global options (verbosity, logging).
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from .logging_utils import get_console, setup_logging
from .commands import (
    accept_topology_cmd,
    compare_cmd,
    compare_variants_cmd,
    contacts_cmd,
    ff_reader_cmd,
    interactors_cmd,
    report_cmd,
    run_cmd,
    topology_cmd,
    validate_cmd,
    guess_ff_path_cmd,
    network_cmd
)

__version__ = "0.1.0"


HELP = (
    "[bold cyan]myPSNs[/bold cyan]: Protein Structure Network analysis from "
    "molecular dynamics trajectories.\n\n"
    "Pipeline: [yellow]ff-reader → topology → interactors → contacts[/yellow]\n\n"
    "Use [bold]mypsns run[/bold] for the full pipeline, or invoke individual "
    "stages directly. Pass a [bold]--config[/bold] YAML file for reproducible runs."
)

app = typer.Typer(
    name="mypsns",
    help=HELP,
    rich_markup_mode="rich",
    no_args_is_help=True,
    pretty_exceptions_show_locals=False,
    add_completion=False,
)


def _version_cb(value: bool) -> None:
    if value:
        get_console().print(f"myPSNs version [bold]{__version__}[/bold]")
        raise typer.Exit()


@app.callback()
def main_callback(
    version: Optional[bool] = typer.Option(
        None, "--version",
        callback=_version_cb, is_eager=True,
        help="Print version and exit.",
    ),
    verbose: int = typer.Option(
        0, "-v", "--verbose", count=True,
        help="Increase verbosity (-v info, -vv debug). Default is INFO.",
    ),
    quiet: bool = typer.Option(
        False, "-q", "--quiet",
        help="Suppress everything below WARNING.",
    ),
    log_file: Optional[Path] = typer.Option(
        None, "--log-file",
        help="Write a DEBUG-level log to this file (alongside stdout).",
    ),
) -> None:
    """Configure global logging before any subcommand runs."""
    setup_logging(verbosity=verbose, quiet=quiet, log_file=log_file)


# ── Subcommand registration ─────────────────────────────────────────────────
app.command(name="run",         help="Run the full pipeline.")(run_cmd.run_pipeline)
app.command(name="ff-reader",   help="Parse a GROMACS forcefield folder.")(ff_reader_cmd.ff_reader)
app.command(name="topology",    help="Match topology residues against the forcefield.")(topology_cmd.topology_recognition)
app.command(name="accept-topology", help="Promote a proposed candidate for residues `topology` could not match.")(accept_topology_cmd.accept_topology)
app.command(name="interactors", help="Assign interaction roles to topology atoms.")(interactors_cmd.interactors)
app.command(name="contacts",    help="Compute contacts across an MD trajectory.")(contacts_cmd.contacts)
app.command(name="validate",    help="Validate a config and its referenced files.")(validate_cmd.validate)
app.command(name="report",      help="Summarise outputs in a run folder.")(report_cmd.report)
app.command(name="guess_ff",    help="Guesses possible forcefield subfolders.")(guess_ff_path_cmd.guess_ff_path)
app.command(name="network",    help="Network construction and subsequent analyses.")(network_cmd.network)
app.command(name="compare",    help="Compare residue matrices with reference (e.g. PyInteraph) graphs.")(compare_cmd.compare)
app.command(name="compare-variants", help="Rank network metrics against ClinVar pathogenicity (experimental).")(compare_variants_cmd.compare_variants)

def main() -> None:
    """Console-script entry point."""
    app()


if __name__ == "__main__":
    main()
