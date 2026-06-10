# myPSNs CLI

A production-grade command-line interface for the **myPSNs** pipeline, built on **Typer + Rich + Pydantic + PyYAML**.

## Directory layout

The files in this drop-in are designed to live inside the existing `mypsns/` package:

```
mypsns/
├── cli.py                          ← Typer app: subcommand registration + global options
├── config.py                       ← Pydantic models + YAML loading/merging
├── logging_utils.py                ← Rich console + log handler + banners/tables
├── commands/
│   ├── __init__.py
│   ├── _helpers.py                 ← shared CLI helpers 
│   ├── run_cmd.py                  ← full-pipeline orchestrator
│   ├── ff_reader_cmd.py
│   ├── topology_cmd.py
│   ├── interactors_cmd.py
│   ├── contacts_cmd.py
│   ├── validate_cmd.py
│   └── report_cmd.py
├── services/                       ← existing pipeline modules (untouched)
│   ├── ff_reader.py
│   ├── topology_recognition.py
│   ├── interactors.py
│   └── contacts.py
└── utils/
```

Plus, alongside the package:

```
example_config.yaml                 ← sample YAML config (drop in your repo's docs/ or examples/)
pyproject_snippet.toml              ← entry-point + dependency snippet for pyproject.toml
```

## Architecture at a glance

The CLI is layered so that each concern owns exactly one file.

**`cli.py`** is the Typer app. It owns nothing but registration:
- global options (`--verbose`, `--quiet`, `--log-file`, `--version`)
- subcommand wiring (`run`, `ff-reader`, `topology`, `interactors`, `contacts`, `validate`, `report`)
- a single entry point `main()` for the `mypsns` console script

**`config.py`** is the source of truth for parameters. Every stage has a Pydantic model (`FFReaderConfig`, `TopologyConfig`, `InteractorsConfig`, `ContactsConfig`); `PipelineConfig` composes them. Values flow through three layers, lowest to highest priority:
1. Pydantic field defaults
2. YAML file (`--config`)
3. CLI option overrides

Pydantic does the schema validation; the helpers do the merging. This means typos in YAML or out-of-range numerics produce structured `ConfigError`s rather than mysterious `KeyError`s deep in the pipeline.

**`logging_utils.py`** centralises Rich. `setup_logging()` is called once by `cli.py`'s callback before any subcommand runs. Stage banners, summary tables, and status messages all flow through helpers here so visual style stays consistent.

**`commands/`** has one file per subcommand. Each module exposes a single function the Typer app registers as a command. They share helpers (`collect_overrides`, `load_or_die`, `dump_run_config`, `warn_skip`, …) through `_helpers.py`. Each command:
1. Collects its CLI arguments into a nested overrides dict
2. Loads the YAML config (if any) and applies overrides
3. Validates required paths
4. Prints inputs and (if `--dry-run`) exits
5. Honours `--resume` / `--overwrite` against per-stage marker files
6. Calls the underlying service module
7. Dumps the resolved config back as `run_config.yaml` for reproducibility

## Quick start

Minimal beginner-friendly invocation:

```bash
mypsns run \
  --topology protein.pdb \
  --trajectory traj.xtc \
  --forcefield amber99sb.ff \
  --outpath results/
```

Reproducible HPC-style invocation:

```bash
mypsns run --config config.yaml --log-file results/run.log
```

CLI overrides on top of a config (e.g. bump worker count for a fat node):

```bash
mypsns run --config config.yaml --nc 32 --target HBonds --target Salt_bridges
```

Individual stages:

```bash
mypsns ff-reader   --forcefield amber99sb.ff --outpath results/
mypsns topology    --topology protein.pdb --forcefield amber99sb.ff --outpath results/
mypsns interactors --topology protein.pdb --outpath results/
mypsns contacts    --topology protein.pdb --trajectory traj.xtc --outpath results/ --nc 16
```

Pre-flight check (no calculations):

```bash
mypsns validate --config config.yaml --deep
```

Post-run summary:

```bash
mypsns report results/
mypsns report results/ --json | jq '.contacts'
```

Resume an interrupted run (skips stages whose output marker is already present):

```bash
mypsns run --config config.yaml --resume
```

Run only some stages:

```bash
mypsns run --config config.yaml --only interactors --only contacts
mypsns run --config config.yaml --skip ff-reader
```

Verbosity:

```bash
mypsns run --config config.yaml          # INFO
mypsns -v run --config config.yaml       # INFO with caller paths
mypsns -vv run --config config.yaml      # DEBUG
mypsns -q run --config config.yaml       # only WARNING+
```

## YAML config format

See `example_config.yaml`. Top-level keys (`topology`, `trajectory`, `forcefield`, `outpath`, `log_file`) sit at the root; everything else nests under stage keys (`ff_reader:`, `topology_recognition:`, `interactors:`, `contacts:`). Unknown keys produce a hard error (`extra="forbid"` on every Pydantic model) so config drift is caught immediately.

## Design rationale

**Why Typer and not argparse/Click directly?** Type-hint–driven option declaration keeps each command's signature literally readable; `rich_help_panel` groups related flags in `--help` without scaffolding; Typer is a thin layer over Click so escape hatches stay available.

**Why Pydantic for config?** Three reasons. First, structured validation errors (`max_to_check: must be ≥ 0`) instead of pipeline crashes ten minutes in. Second, the same model that validates a config file also serialises the resolved run back to YAML for reproducibility (`run_config.yaml` is dumped after each successful stage). Third, `extra="forbid"` catches typos in YAML keys immediately.

**Why one file per subcommand?** Adding a new stage means adding one module and one `app.command()` line in `cli.py`. No god-file to navigate. Shared logic stays in `_helpers.py`; the command files contain only their own Typer option declarations and the call into the underlying service.

**Resume/overwrite semantics**: each stage has a canonical "completion marker" file (e.g. `Contacts/Hbonds_final.csv`). `--resume` skips a stage iff its marker exists; `--overwrite` forces re-execution regardless. Default behaviour is "skip with warning if marker exists" — protects against accidental overwrites of long runs.

**HPC friendliness**: every run dumps `run_config.yaml` to `<outpath>/`. Pair with `--log-file` and you get a reproducible artefact you can ship between machines. The `run` command also logs environment info (Python version, OS, machine arch, timestamp) before any work starts.

**Error handling**: every command wraps the service call in try/except that prints a clean error line via `print_error()` and exits with a non-zero code. Stack traces only surface at `-vv`. `ConfigError` becomes "exit code 2 + Rich error"; missing service modules become "exit code 3"; pipeline crashes become "exit code 1".

## Adding a new subcommand

1. Create `mypsns/commands/<new>_cmd.py` with a single function whose parameters are `typer.Option(...)` calls.
2. Re-export it from `mypsns/commands/__init__.py`.
3. Register it in `cli.py`:

   ```python
   app.command(name="my-new", help="…")(new_cmd.my_new)
   ```

That's it — no parser plumbing, no manual help-text wrangling.

## Adding a new config field

1. Add the field to the appropriate Pydantic model in `config.py` (with default + description).
2. Add a matching `typer.Option(None, "--…", help="…")` to the relevant command.
3. Add it to the `collect_overrides({...})` dict in that command.

The validation, YAML loading, override merging, and `--help` formatting come for free.

## Service module assumptions

The CLI wraps the existing service classes with these expected signatures:

| Service | Class | Method calls |
|---|---|---|
| ff_reader.py | `FFReader(inpath, outpath, na_values, nonbonded_name, additional_residues)` | `.run()` |
| topology_recognition.py | `Topology_recognition(topology, outpath, ff_path, na_values, outlog)` | `.get_restypes(verbose=…)`, `.get_ff_parameters(verbose=…, maxToCheck=…)` |
| interactors.py | `Interactors(topology, general_outpath, na_values, outlog)` | `.get_complete_topology(verbose=…)`, `.add_inttypes(verbose=…, **kwargs)` |
| contacts.py | `Contacts(topology, trajectory, general_outpath)` | `.run(**kwargs)` |

The existing `Topology_recognition` and `Interactors` classes already match this. The CLI assumes `FFReader` and `Contacts` expose `.run()` orchestrator methods; if your current modules differ, the only change needed is in the `try: from ..services.X import Y` block of the corresponding `*_cmd.py` and the call site. The CLI itself stays untouched.
