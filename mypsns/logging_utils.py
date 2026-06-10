"""
Rich-based logging utilities for the myPSNs CLI.

Sets up a colourful, scientific-user-friendly log/stdout stream and provides
helpers for banners, summary tables, and standardised status messages.
"""
from __future__ import annotations

import logging
import platform
from datetime import datetime
from pathlib import Path
from typing import Any, Optional
from contextlib import contextmanager

from rich.console import Console
from rich.logging import RichHandler
from rich.panel import Panel
from rich.table import Table
from rich.text import Text


_console: Optional[Console] = None

def get_console() -> Console:
    """Lazy singleton for the Rich console."""
    global _console
    if _console is None:
        _console = Console(stderr=False, soft_wrap=False)
    return _console

def setup_logging(
    verbosity: int = 0,
    quiet: bool = False,
    log_file: Optional[Path] = None,
) -> None:
    """
    Configure root logging.

    Parameters
    ----------
    verbosity:
        0 = INFO (default), 1 = INFO with caller info, 2 = DEBUG.
    quiet:
        Overrides verbosity; only WARNING and above are shown on stdout.
    log_file:
        Optional file to receive a DEBUG-level full log alongside stdout.
    """
    if quiet:
        console_level = logging.WARNING
    elif verbosity >= 2:
        console_level = logging.DEBUG
    else:
        console_level = logging.INFO

    root = logging.getLogger()
    for h in list(root.handlers):
        root.removeHandler(h)

    rich = RichHandler(
        console=get_console(),
        rich_tracebacks=True,
        show_path=verbosity >= 2,
        show_time=True,
        markup=True,
        log_time_format="%H:%M:%S",
    )
    rich.setLevel(console_level)
    root.addHandler(rich)

    file_attached = False
    file_error: Optional[Exception] = None
    if log_file is not None:
        try:
            log_file = Path(log_file).expanduser()
            log_file.parent.mkdir(parents=True, exist_ok=True)
            fh = logging.FileHandler(log_file, mode="a", encoding="utf-8")
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(logging.Formatter(
                "%(asctime)s %(levelname)-7s %(name)s | %(message)s"
            ))
            root.addHandler(fh)
            file_attached = True
        except OSError as e:          # permissions, read-only FS, bad path, …
            file_error = e

    root.setLevel(logging.DEBUG if file_attached else console_level)
    
    # Silence third-party loggers
    for noisy in ("MDAnalysis", "matplotlib", "numexpr", "h5py", "PIL", "filelock"):
        logging.getLogger(noisy).setLevel(logging.WARNING)

    if log_file is not None and not file_attached:
        print_warning(
            f"Could not open log file '{log_file}': {file_error}. "
            f"Continuing with console logging only."
        )

# ─────────────────────────────────────────────────────────────────────────────
# Status helpers
# ─────────────────────────────────────────────────────────────────────────────

def print_banner(title: str, subtitle: Optional[str] = None) -> None:
    """Print a section banner."""
    body = Text(title, style="bold cyan")
    if subtitle:
        body.append(f"\n{subtitle}", style="dim")
    get_console().print(
        Panel(body, expand=False, border_style="cyan", padding=(0, 2))
    )

def print_summary(
    title: str,
    data: dict[str, Any],
    style_key: str = "cyan",
) -> None:
    """Print a 2-column (parameter, value) summary."""
    tbl = Table(
        title=title, show_header=True,
        header_style="bold magenta", show_lines=False,
    )
    tbl.add_column("Parameter", style=style_key, no_wrap=True)
    tbl.add_column("Value", style="white")
    for k, v in data.items():
        tbl.add_row(str(k), _fmt_value(v))
    get_console().print(tbl)

def print_kv_table(title: str, rows: list[tuple[str, str, str]]) -> None:
    """Print a 3-column status table (label, status, detail)."""
    tbl = Table(title=title, show_header=True, header_style="bold magenta")
    tbl.add_column("Check", style="cyan", no_wrap=True)
    tbl.add_column("Status", style="white", no_wrap=True)
    tbl.add_column("Detail", style="white")
    for label, status, detail in rows:
        tbl.add_row(label, status, detail)
    get_console().print(tbl)

def print_env_info() -> None:
    """Log a small banner with environment info (useful for reproducibility)."""
    info = {
        "started": datetime.now().isoformat(timespec="seconds"),
        "python":  platform.python_version(),
        "system":  f"{platform.system()} {platform.release()}",
        "machine": platform.machine(),
    }
    print_summary("Run environment", info, style_key="dim")

def _fmt_value(v: Any) -> str:
    if isinstance(v, (list, tuple)):
        return ", ".join(str(x) for x in v) if v else "—"
    if v is None:
        return "[dim]auto[/dim]"
    if isinstance(v, bool):
        return "[green]on[/green]" if v else "[dim]off[/dim]"
    return str(v)

def print_error(msg: str) -> None:
    get_console().print(f"[bold red]✗ Error:[/bold red] {msg}")


def print_warning(msg: str) -> None:
    get_console().print(f"[bold yellow]⚠ Warning:[/bold yellow] {msg}")


def print_success(msg: str) -> None:
    get_console().print(f"[bold green]✓[/bold green] {msg}")


def print_info(msg: str) -> None:
    get_console().print(f"[cyan]→[/cyan] {msg}")

@contextmanager
def stage_log_file(
    path: Path,
    level: int = logging.DEBUG,
    fmt:   str = "%(asctime)s %(name)s %(levelname)s: %(message)s",
):
    """
    Attach a per-stage FileHandler to the root logger for the duration of a
    pipeline stage, then remove and close it.

    Mirrors how the Contacts service writes Contacts.log, but driven from the
    caller so the simpler services (ff_reader, topology_recognition, interactors)
    don't need to add handlers themselves.

    The handler is set to DEBUG; the root level (configured by setup_logging)
    decides what actually reaches it — use `-vv` if you want DEBUG details in
    the stage log too.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(path, mode="w")
    fh.setLevel(level)
    fh.setFormatter(logging.Formatter(fmt))
    root = logging.getLogger()
    root.addHandler(fh)
    try:
        yield fh
    finally:
        root.removeHandler(fh)
        fh.close()