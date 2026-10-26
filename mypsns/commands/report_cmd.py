"""`mypsns report` — summarise outputs in a run folder."""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional
import typer

from ..logging_utils import (
    get_console, print_banner, print_error, print_info, print_kv_table,
    print_summary, print_warning,
)

log = logging.getLogger(__name__)


def report(
    outpath: Path = typer.Argument(
        ..., help="Run output folder to inspect.",
    ),
    json_output: bool = typer.Option(
        False, "--json",
        help="Emit the report as JSON instead of formatted tables.",
    ),
) -> None:
    """
    Summarise outputs found in a run folder.

    Walks the standard pipeline subfolders (forcefield/, Topology_recognition/,
    Interactors/, Adjacency/<inttype>/, Network/<inttype>/) and reports which
    stages have produced final outputs.
    """
    outpath = Path(outpath).expanduser()
    if not outpath.exists():
        print_error(f"outpath does not exist: {outpath}")
        raise typer.Exit(code=2)

    payload = {
        "outpath":              str(outpath),
        "ff_reader":            _ff_reader_status(outpath),
        "topology_recognition": _topology_status(outpath),
        "interactors":          _interactors_status(outpath),
        "contacts":             _contacts_status(outpath),
        "network":              _network_status(outpath),
    }

    if json_output:
        get_console().print_json(json.dumps(payload, default=str))
        return

    print_banner("report", subtitle=str(outpath))

    # Stage table
    stage_rows: list[tuple[str, str, str]] = []
    for stage in ("ff_reader", "topology_recognition", "interactors", "contacts", "network"):
        s = payload[stage]
        status = "[green]✓[/green]" if s["complete"] else "[yellow]—[/yellow]"
        stage_rows.append((stage, status, s.get("summary", "")))
    print_kv_table("Stage status", stage_rows)

    # Interactors details
    int_info = payload["interactors"]
    if int_info["complete"]:
        print_summary("Interactors", {
            "variants":     int_info.get("n_variants", "—"),
            "residues":     int_info.get("n_residues", "—"),
            "atoms":        int_info.get("n_atoms", "—"),
        })

    # Contacts details
    cn_info = payload["contacts"]
    for inttype, info in cn_info.get("schemes", {}).items():
        print_summary(f"Contacts — {inttype}", {
            "final tables":    info["finals"],
            "adjacency files": info["adjacencies"],
            "energy matrices": info["energies"],
            "elbow CSVs":      info["elbow_csvs"],
        })

    nw_info = payload["network"]
    for inttype, info in nw_info.get("schemes", {}).items():
        print_summary(f"Network — {inttype}", {
            "node tables":   info["node_tables"],
            "PyMOL scripts": info["pml"],
            "matrices":      info["matrices"],
        })

# ─────────────────────────────────────────────────────────────────────────────
# Stage-specific status collectors
# ─────────────────────────────────────────────────────────────────────────────

def _ff_reader_status(outpath: Path) -> dict:
    candidates = [
        d for d in outpath.iterdir()
        if d.is_dir() and (d / "final_ff_residues.csv").exists()
    ]
    if not candidates:
        return {"complete": False, "summary": "no parsed forcefield found"}
    names = [c.name for c in candidates]
    return {
        "complete": True,
        "summary":  f"{len(candidates)} forcefield(s): {', '.join(names)}",
        "folders":  [str(c) for c in candidates],
    }


def _topology_status(outpath: Path) -> dict:
    topo = outpath / "Topology_recognition"
    if not topo.exists():
        return {"complete": False, "summary": "stage not run"}

    complete  = topo / "complete_residue_topologies.csv"
    unmatched = topo / "unmatching_residue_topologies.csv"
    if complete.exists():
        return {"complete": True, "summary": "complete_residue_topologies.csv present"}
    if unmatched.exists():
        return {"complete": False,
                "summary": "stage ran but some residues unmatched"}
    return {"complete": False, "summary": "stage partial / not finished"}


def _interactors_status(outpath: Path) -> dict:
    intd = outpath / "Interactors"
    if not intd.exists():
        return {"complete": False, "summary": "stage not run"}

    res_idx = intd / "residue_index.csv"
    aim     = intd / "atom_index_map.csv"
    variants_dir = intd / "found_topologies"
    if not res_idx.exists() or not aim.exists():
        return {"complete": False, "summary": "outputs incomplete"}

    n_variants = (
        len(list(variants_dir.glob("*.csv"))) if variants_dir.exists() else 0
    )
    n_residues = _csv_row_count(res_idx)
    n_atoms    = _csv_row_count(aim)

    return {
        "complete":   True,
        "summary":    f"{n_variants} variants, {n_residues} residues",
        "n_variants": n_variants,
        "n_residues": n_residues,
        "n_atoms":    n_atoms,
    }


def _scheme_dirs(root: Path) -> list:
    """<root>/<inttype>/ folders (named_interactors, q_interactors), sorted."""
    if not root.exists():
        return []
    return sorted(d for d in root.iterdir()
                  if d.is_dir() and d.name.endswith("_interactors"))


def _contacts_status(outpath: Path) -> dict:
    # contacts writes one Adjacency/<inttype>/ folder per interactor scheme
    schemes = {}
    for d in _scheme_dirs(outpath / "Adjacency"):
        finals = sorted(p.name for p in d.glob("*_final.csv"))
        if not finals:
            continue
        schemes[d.name] = {
            "finals":      finals,
            "adjacencies": sorted(p.name for p in d.glob("*_adjacency.dat")),
            "energies":    len(list(d.glob("*_energy.dat"))),
            "elbow_csvs":  len(list(d.glob("*graph_report/csvs/*_elbow.csv"))),
        }
    if not schemes:
        return {"complete": False, "summary": "stage not run"}
    return {
        "complete": True,
        "summary":  ", ".join(f"{k}: {len(v['finals'])} final table(s)"
                              for k, v in schemes.items()),
        "schemes":  schemes,
    }


def _network_status(outpath: Path) -> dict:
    schemes = {}
    for d in _scheme_dirs(outpath / "Network"):
        nodes = list((d / "csv").glob("nodes_*.csv"))
        if not nodes:
            continue
        schemes[d.name] = {
            "node_tables": len(nodes),
            "pml":         len(list(d.glob("*_metrics/*.pml"))),
            "matrices":    len(list(d.glob("matrices/*_adjacency.dat"))),
        }
    if not schemes:
        return {"complete": False, "summary": "stage not run"}
    return {
        "complete": True,
        "summary":  ", ".join(f"{k}: {v['node_tables']} node table(s)"
                              for k, v in schemes.items()),
        "schemes":  schemes,
    }


def _csv_row_count(path: Path) -> int:
    """Count data rows (header excluded) without loading the whole file."""
    try:
        with path.open("r") as f:
            return max(sum(1 for _ in f) - 1, 0)
    except Exception:
        return 0
