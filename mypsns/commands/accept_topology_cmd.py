"""`mypsns accept-topology` — promote a proposed candidate for residues that
topology recognition could not match, so the pipeline can continue."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import List, Optional

import typer
from rich.table import Table

from ..config import DEFAULT_NA_VALUES
from ..logging_utils import (
    get_console, print_banner, print_error, print_info, print_success,
    print_summary, print_warning, stage_log_file,
)
from ._helpers import collect_overrides, load_or_die, require_path


log = logging.getLogger(__name__)

PANEL_GENERAL = "General"
PANEL_BEHAVE  = "Behaviour"


def accept_topology(
    choices: Optional[List[str]] = typer.Argument(
        None, metavar="[TOPOLOGY=CHOICE]...",
        help=(
            "Residue topology to resolve and what to use for it. CHOICE is a "
            "candidate index n (possible_topologies/<TOPOLOGY>_uncertain_atoms_<n>.csv) "
            "or the path of a CSV with the same columns, e.g. a hand-edited "
            "candidate. Without choices, list the unresolved topologies and "
            "their candidates."
        ),
    ),
    config: Optional[Path] = typer.Option(
        None, "--config",
        help="YAML config file (only outpath and na_values are used).",
        rich_help_panel=PANEL_GENERAL,
    ),
    outpath: Optional[Path] = typer.Option(
        None, "--outpath", "-o",
        help="Root output folder of the run (contains Topology_recognition/).",
        rich_help_panel=PANEL_GENERAL,
    ),
    dry_run: bool = typer.Option(
        False, "--dry-run",
        help="Check the choices and show the renames; write nothing.",
        rich_help_panel=PANEL_BEHAVE,
    ),
) -> None:
    """
    Resolve residues that `topology` could not match exactly.

    Copies the chosen candidate to `found_topologies/<TOPOLOGY>_ff_atoms.csv`
    with each atom's structure name in `structure_atom_name` (`interactors`
    renames the atoms to those names), and fills
    `ff_atoms` in `unmatching_residue_topologies.csv`. Once every residue
    topology is resolved it writes `complete_residue_topologies.csv`, and
    `interactors` can run.
    """
    # A config is optional: accepting only needs the run folder. When one is
    # given it supplies outpath (unless -o overrides it) and na_values.
    if config is not None:
        cfg = load_or_die(config, collect_overrides({"outpath": outpath}))
        outpath, na_values = cfg.outpath, cfg.topology_recognition.na_values
    else:
        na_values = DEFAULT_NA_VALUES.copy()
    outpath = require_path("outpath", outpath, must_exist=True)
    topo_out = outpath / "Topology_recognition"
    if not topo_out.is_dir():
        print_error(f"{topo_out} not found: run `mypsns topology` first.")
        raise typer.Exit(code=2)

    parsed = _parse_choices(choices or [])

    from ..services.topology_recognition import (  # type: ignore
        TopologyAcceptError, accept_topologies, describe_unresolved_topologies,
    )

    print_banner("accept-topology", subtitle=str(topo_out))

    if not parsed:
        _print_unresolved(describe_unresolved_topologies(topo_out, na_values), topo_out, outpath)
        return

    print_summary("Inputs", {
        "outpath": topo_out,
        "choices": " ".join(f"{k}={v}" for k, v in parsed.items()),
        "dry_run": dry_run,
    })

    try:
        # its own log: stage_log_file truncates, and topology_recognition.log
        # must keep the recognition run that halted
        with stage_log_file(topo_out / "accept_topology.log"):
            result = accept_topologies(topo_out, parsed, na_values, dry_run=dry_run)
    except TopologyAcceptError:
        # accept_topologies already logged every problem; don't repeat them
        print_error("accept-topology: nothing was written (see the problems above).")
        raise typer.Exit(code=1)

    _print_accepted(result["accepted"])
    if dry_run:
        print_info("Dry run — nothing written.")
    if result["pending"]:
        print_warning(
            f"Still unresolved: {', '.join(result['pending'])}. Accept a candidate "
            f"for each (run `mypsns accept-topology -o {outpath}` to list them).")
        return
    if dry_run:
        print_info("These choices would resolve every residue topology.")
        return

    print_success(f"every residue topology resolved → {topo_out / 'complete_residue_topologies.csv'}")
    if (outpath / "Interactors" / "residue_index.csv").exists():
        print_warning(
            "Interactors/ holds outputs built before this change: re-run "
            "`mypsns interactors --overwrite`, then `contacts`.")
    else:
        print_info("Next: `mypsns interactors` (or `mypsns run`, which skips the finished stages).")


def _parse_choices(raw: List[str]) -> dict:
    """['MET_1=0', 'LIG=lig.csv'] -> {'MET_1': '0', 'LIG': 'lig.csv'}."""
    parsed: dict = {}
    for item in raw:
        key, sep, value = item.partition("=")
        key, value = key.strip(), value.strip()
        if not sep or not key or not value:
            print_error(f"Bad choice {item!r}: expected TOPOLOGY=CHOICE, e.g. MET_1=0.")
            raise typer.Exit(code=2)
        if key in parsed:
            print_error(f"{key} is given more than once.")
            raise typer.Exit(code=2)
        parsed[key] = value
    return parsed


def _format_renames(pairs) -> str:
    """'H1→HT1, H2→HT2' from a pair_candidate_atoms table; '!' marks a pair
    whose masses differ (that candidate cannot be accepted as is)."""
    renamed = pairs[pairs["renamed"]]
    if renamed.empty:
        return "none (same atom names)"
    parts = []
    for r in renamed.itertuples():
        mark = "" if r.mass_ok else f" [red]! {r.ff_mass:g}≠{r.top_mass:g}[/red]"
        parts.append(f"{r.ff_name}→{r.top_name}{mark}")
    return ", ".join(parts)


def _print_unresolved(entries: list, topo_out: Path, outpath: Path) -> None:
    if not entries:
        if (topo_out / "complete_residue_topologies.csv").exists():
            print_success("Every residue topology is matched: nothing to accept.")
        else:
            print_info("No unmatching_residue_topologies.csv: run `mypsns topology` first.")
        return

    console = get_console()
    for e in entries:
        status = (f"accepted from {e['accepted_from']}" if e["accepted_from"]
                  else "[yellow]unresolved[/yellow]")
        tbl = Table(
            title=f"{e['key']}  ({e['resname']}, {e['n_residues']} residue(s)) — {status}",
            show_header=True, header_style="bold magenta",
        )
        tbl.add_column("#", style="cyan", no_wrap=True)
        tbl.add_column("forcefield residue(s)", style="white")
        tbl.add_column("renames (forcefield → structure)", style="white")
        for c in e["candidates"]:
            if c["same_as"] is not None:
                renames = f"[dim]identical to #{c['same_as']}[/dim]"
            elif c["error"]:
                renames = f"[red]{c['error']}[/red]"
            else:
                renames = _format_renames(c["pairs"])
            tbl.add_row(str(c["index"]), ", ".join(c["ff_resnames"]), renames)
        if not e["candidates"]:
            tbl.add_row("—", "no candidate", "write the topology by hand and pass its path")
        console.print(tbl)

    first = entries[0]["key"]
    print_info(
        f"Accept with `mypsns accept-topology -o {outpath} {first}=<#>` "
        f"(one TOPOLOGY=CHOICE per topology; CHOICE may also be a CSV path). "
        f"A [red]![/red] pair renames atoms of different mass: that candidate "
        f"needs a hand-edited copy.")


def _print_accepted(accepted: list) -> None:
    tbl = Table(title="Accepted", show_header=True, header_style="bold magenta")
    tbl.add_column("topology", style="cyan", no_wrap=True)
    tbl.add_column("from", style="white")
    tbl.add_column("written as", style="white")
    tbl.add_column("renames (forcefield → structure)", style="white")
    for a in accepted:
        renames = ", ".join(f"{k}→{v}" for k, v in a["renames"].items()) or "none"
        tbl.add_row(a["key"], a["source"], f"found_topologies/{a['target']}", renames)
    get_console().print(tbl)
