"""`mypsns compare` — compare residue matrices with reference graphs (e.g. the
PyInteraph 1 validation). Logic lives in utils/validation_utilities.py."""
from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Dict, List, Optional

import typer
from rich.table import Table

from ..logging_utils import (
    get_console, print_banner, print_error, print_info, print_success,
    print_warning,
)
from ..utils.validation_utilities import (
    DEFAULT_TYPES, PRESETS, VALIDATION_ROOT_ENV,
    adjacency_paths, compare_files, preset_paths, read_label_names,
    reference_paths, write_comparison_csvs,
)


log = logging.getLogger(__name__)

PANEL_INPUT  = "Inputs"
PANEL_OUTPUT = "Output"
PANEL_CHECK  = "Checks"


def _parse_kv(items: Optional[List[str]], flag: str) -> Dict[str, str]:
    """['TYPE=VALUE', ...] → {TYPE: VALUE}; exits with a clear error otherwise."""
    out: Dict[str, str] = {}
    for item in items or []:
        key, sep, value = item.partition("=")
        if not sep or not key.strip() or not value.strip():
            print_error(f"{flag} expects TYPE=VALUE, got {item!r}")
            raise typer.Exit(code=2)
        out[key.strip()] = value.strip()
    return out


def compare(
    preset: Optional[str] = typer.Option(
        None, "--preset",
        help=(f"Legacy dataset from the PyInteraph validation: "
              f"{', '.join(sorted(PRESETS))}. Paths are taken under --data-root."),
        rich_help_panel=PANEL_INPUT,
    ),
    data_root: Optional[Path] = typer.Option(
        None, "--data-root",
        help=(f"Root for --preset paths (default: ${VALIDATION_ROOT_ENV}, else "
              "/home/davide/storage/lab_dati)."),
        rich_help_panel=PANEL_INPUT,
    ),
    truth_dir: Optional[Path] = typer.Option(
        None, "--truth-dir",
        help="Folder with reference graphs named like PyInteraph: "
             "{hb,hc,sb}-graph<suffix>.dat.",
        rich_help_panel=PANEL_INPUT,
    ),
    truth_suffix: str = typer.Option(
        "", "--truth-suffix",
        help="Text between the PyInteraph stem and '.dat' (e.g. _1micros).",
        rich_help_panel=PANEL_INPUT,
    ),
    computed_dir: Optional[Path] = typer.Option(
        None, "--computed-dir",
        help="myPSNs Adjacency/<inttype>/ folder ({Hbonds,HC,SB}_adjacency.dat).",
        rich_help_panel=PANEL_INPUT,
    ),
    truth_files: Optional[List[str]] = typer.Option(
        None, "--truth",
        help="TYPE=PATH reference matrix; overrides --truth-dir/--preset for "
             "TYPE. Any TYPE name works (e.g. macroIIN=...). Repeatable.",
        rich_help_panel=PANEL_INPUT,
    ),
    computed_files: Optional[List[str]] = typer.Option(
        None, "--computed",
        help="TYPE=PATH computed matrix, paired with --truth by TYPE. Repeatable.",
        rich_help_panel=PANEL_INPUT,
    ),
    types: Optional[List[str]] = typer.Option(
        None, "--type",
        help=f"Only compare these types (default: every type with both files; "
             f"built-in: {', '.join(DEFAULT_TYPES)}). Repeatable.",
        rich_help_panel=PANEL_INPUT,
    ),
    out_dir: Optional[Path] = typer.Option(
        None, "--out-dir", "-o",
        help="Write possible_sele_mismatch.csv, possible_params_mismatch.csv and "
             "comparison_summary.csv here. Without it nothing is written.",
        rich_help_panel=PANEL_OUTPUT,
    ),
    labels: Optional[Path] = typer.Option(
        None, "--labels",
        help="adjacency_labels.dat: adds residue labels (label_row/label_col) "
             "to the mismatch CSVs.",
        rich_help_panel=PANEL_OUTPUT,
    ),
    atol: float = typer.Option(
        1e-2, "--atol",
        help="Entries closer than this (np.isclose) count as equal.",
        rich_help_panel=PANEL_CHECK,
    ),
    max_pct_diff: Optional[List[str]] = typer.Option(
        None, "--max-pct-diff",
        help="TYPE=PCT: exit 1 if TYPE's % weight difference exceeds PCT "
             "(e.g. SB=0.2). Repeatable.",
        rich_help_panel=PANEL_CHECK,
    ),
) -> None:
    """
    Compare residue × residue matrices with reference ones, type by type.

    Reports edges (truth / computed / shared), entries whose values differ, and
    the % weight difference 100·Σ|computed−truth|/Σtruth on the strict lower
    triangle (diagonal ignored). With --out-dir, also writes the mismatch CSVs
    the old prove.py produced.

    Examples:

      mypsns compare --preset pyinteraph_testcase -o cmp/

      mypsns compare --truth-dir ref/ --truth-suffix _1micros
        --computed-dir run/Adjacency/named_interactors

      mypsns compare --truth macroIIN=macro.dat
        --computed macroIIN=run/Network/named_interactors/matrices/persistency_unweighted_network_adjacency.dat
    """
    truth_paths: Dict[str, Path] = {}
    comp_paths:  Dict[str, Path] = {}
    if preset is not None:
        try:
            t, c = preset_paths(preset, data_root)
        except ValueError as e:
            print_error(str(e))
            raise typer.Exit(code=2)
        truth_paths.update(t)
        comp_paths.update(c)
    if truth_dir is not None:
        truth_paths.update(reference_paths(truth_dir, truth_suffix))
    if computed_dir is not None:
        comp_paths.update(adjacency_paths(computed_dir))
    truth_paths.update({k: Path(v) for k, v in _parse_kv(truth_files, "--truth").items()})
    comp_paths.update({k: Path(v) for k, v in _parse_kv(computed_files, "--computed").items()})
    thresholds = {k: float(v) for k, v in _parse_kv(max_pct_diff, "--max-pct-diff").items()}

    wanted = list(types) if types else list(dict.fromkeys([*truth_paths, *comp_paths]))
    unpaired = [t for t in wanted if t not in truth_paths or t not in comp_paths]
    pairs = {t: (truth_paths[t], comp_paths[t]) for t in wanted if t not in unpaired}
    if unpaired:
        print_warning(f"no truth/computed pair for: {', '.join(unpaired)} — skipped")
    if not pairs:
        print_error("Nothing to compare. Give --preset, --truth-dir + --computed-dir, "
                    "or --truth/--computed TYPE=PATH pairs.")
        raise typer.Exit(code=2)
    missing = [str(p) for tp in pairs.values() for p in tp if not p.exists()]
    if missing:
        print_error("Missing input file(s):\n  " + "\n  ".join(missing))
        raise typer.Exit(code=2)

    print_banner("compare", subtitle=preset or "custom paths")
    results = compare_files(pairs, atol=atol)

    tbl = Table(title=f"Truth (T) vs computed (C) — lower triangle, atol={atol:g}",
                show_header=True, header_style="bold magenta")
    # T = truth, C = computed; "≠ value" = shared edges whose values differ.
    for col in ("Type", "Edges T", "Edges C", "Shared", "Only T", "Only C",
                "≠ value", "% diff", "Status"):
        tbl.add_column(col, justify="left" if col in ("Type", "Status") else "right")
    failed = False
    for res in results:
        limit = thresholds.get(res.inttype)
        if res.error:
            status, failed = f"[red]{res.error}[/red]", True
        elif limit is not None and not res.pct_diff <= limit:   # nan fails too
            status, failed = f"[red]FAIL > {limit:g} %[/red]", True
        elif limit is not None:
            status = f"[green]ok ≤ {limit:g} %[/green]"
        else:
            status = ""
        pct = "n/a" if math.isnan(res.pct_diff) else f"{res.pct_diff:.4g}"
        tbl.add_row(res.inttype, str(res.edges_truth), str(res.edges_computed),
                    str(res.shared), str(res.only_truth), str(res.only_computed),
                    str(res.value_mismatch), pct, status)
    get_console().print(tbl)
    unchecked = sorted(set(thresholds) - set(pairs))
    if unchecked:
        print_warning(f"--max-pct-diff given for types not compared: {', '.join(unchecked)}")

    if out_dir is not None:
        names = read_label_names(labels) if labels is not None else None
        written = write_comparison_csvs(results, out_dir, labels=names)
        print_info("Wrote:\n  " + "\n  ".join(str(p) for p in written))

    if failed:
        print_error("compare: at least one type failed.")
        raise typer.Exit(code=1)
    print_success("compare complete")
