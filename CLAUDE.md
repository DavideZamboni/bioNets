# CLAUDE.md

myPSNs: Protein Structure Network analysis from MD trajectories (GROMACS forcefield + topology + trajectory →
H-bond / salt-bridge / hydrophobic-cluster contacts → persistency/energy matrices → residue networks + PyMOL).
It is a re-implementation of **PyInteraph 1** and is validated against it (see below).
`README.md` is the user-facing reference (flags, defaults, output layouts). Read it before changing behaviour.

## Environment

```bash
conda env create -f environment.yml && conda activate mypsns   # python 3.11, installs `pip install -e .`
mypsns --help                                                  # console script → mypsns/cli.py:main
```

There is effectively **no test suite**: `tests/test_basic.py` is a 2-line placeholder. Check changes by running the
CLI on a real run folder (below), not with pytest. (The README says a test keeps the two `NAMED_INTERACTOR_DEFAULTS`
in sync. That test does not exist.)

## Architecture

Pipeline: `ff-reader → topology → interactors → contacts` (chained by `mypsns run`), then `network`, which is
always run separately.

| Layer | Path | Role |
|---|---|---|
| CLI entry | `mypsns/cli.py` | Typer app; global `-v/-q/--log-file`; registers subcommands |
| Commands | `mypsns/commands/*_cmd.py` | Thin Typer wrappers: collect flags → `build_pipeline_config` → print Rich summary → call service |
| Config | `mypsns/config.py` | Pydantic models per stage (`FFReaderConfig`, `TopologyConfig`, `InteractorsConfig`, `ContactsConfig`, `NetworkConfig`, `PipelineConfig`). Precedence **flags > YAML > defaults**; the resolved config is dumped back to `run_config.yaml` |
| Services | `mypsns/services/*.py` | Stage logic as classes (`Topology_recognition`, `Contacts`, `Network`, …) |
| Utils | `mypsns/utils/*_utilities.py` | Per-stage helper functions (most of the code lives here) |
| Logging | `mypsns/logging_utils.py` | Rich console, banners, `print_summary` tables |
| Resources | `mypsns/resources/sol.csv` | Exact-match water table (workaround for inconsistent ff water parsing) |

Run-folder layout (`<outpath>/`): `forcefield/`, `Topology_recognition/`, `Interactors/` (incl. `found_topologies/*.csv`),
`Adjacency/<inttype>/` (contacts output: `{SB,Hbonds,HC}_adjacency.dat`, energy `.dat`s, `adjacency_labels.dat`),
`Network/` (graphs, node/edge CSVs, `.pml`, `Network/<inttype>/matrices/`).

### Adding or changing an option
1. Add a `Field` to the stage's model in `config.py`. That model is the single source of truth for type, default and constraints.
2. Add a `typer.Option(None, ...)` in the `*_cmd.py`, put it in the overrides dict passed to `build_pipeline_config`,
   add it to the `print_summary` block, and forward it to the service call.
3. Document it in the README flag table for that stage.

## Conventions and gotchas

- **`na_values`**: every `pd.read_csv` uses the shared `DEFAULT_NA_VALUES`, which **excludes `"NA"`** so sodium atoms
  are not read as missing. Keep passing it.
- **Two interactor schemes**: `named_interactors` (default, PyInteraph-compatible, name-based AND/OR rules) and
  `q_interactors` (charge/connectivity based, less validated).
- **`NAMED_INTERACTOR_DEFAULTS` is duplicated** in `config.py` and `utils/interactors_utilities.py`. Edit both.
- User-supplied `named_interactors_rules` **replace** the defaults entirely. They are not merged.
- Interactor roles are written to disk by `interactors`. After changing the rules, re-run `interactors --overwrite`
  and then `contacts`. Stages skip when their output already exists.
- HC centres of mass use **forcefield masses**, not MDAnalysis-guessed ones (guessing turns `CAL` into carbon). The
  salt-bridge charge gate (`--min-residue-charge`, 0.8 e) uses forcefield charges.
- Default cutoffs match PyInteraph 1: SB 4.5 Å, HC CoM 5.0 Å, H-bond D–A 4.5 Å with D–H–A ≥ 120°. Don't change them
  without re-running the validation.
- Topology recognition now **halts** with `TopologyAssociationError` when a residue has no perfect forcefield match
  (`services/topology_recognition.py`). For large custom residues, raise `--max-to-check` toward `-1`.
- Labels: `<source>` names the weighted graph and `<source>_unweighted` the binary one. In `Network.load_static`,
  persistency combines targets with element-wise **max** (union, like a macro-IIN) in both weighting modes. Energy
  combines them with **sum**. This is the default of `load_combined_per_type` (`combine=None` picks the rule from the
  source), so callers should not pass `combine` unless they want to override it.
- The comments and docstrings are long and explain *why*. Match that style.

## Validation against PyInteraph 1

The contacts output (adjacency matrices) has been validated against PyInteraph 1 on a CHARMM test case
(149 residues). Results, re-checked 2026-10-26:

| Type | Edges PyInteraph / myPSNs / shared | Weight difference |
|---|---|---|
| SB | 21 / 21 / 21 | 0.07 % |
| HC | 79 / 79 / 79 | 0.11 % |
| H-bonds | 282 / 306 / 279 | 4.56 % |

The ~4.5 % H-bond gap comes from the detection engine (MDAnalysis `HydrogenBondAnalysis` vs PyInteraph's geometric
test), not from a parameter, so tuning cannot remove it. Treat SB/HC drifting above ~0.1 % as a regression.

- Comparison script: `mypsns/prove.py`. It is an ad-hoc script with hard-coded paths. It compares
  `risultati_paragone/pyinter_test/{hb,hc,sb}-graph.dat` with
  `psns_out/cli_outs/pyinteraph_testcase/Adjacency/named_interactors/*_adjacency.dat`, both under
  `/home/davide/storage/lab_dati/`, and writes `possible_{sele,params}_mismatch.csv` to `lab_dati/`. A second dataset
  (`rec_grk1`, recoverin–GRK1, compared against `risultati_paragone/rec/*_1micros.dat`) is defined above it, but the
  PyInteraph dict overwrites it. Swap the order to compare rec_grk1.
- The test case ran with `named_interactors` and default cutoffs.

## Work in progress (uncommitted as of 2026-10-26)

- `network`: `--weighting {weighted,unweighted,both}` (binary macro-IIN analogue), `--write-matrices` (dumps
  `Network/<inttype>/matrices/*.dat` in Adjacency format via `Network.write_matrices`), threshold semantics
  (weighted `>`, unweighted `≥`). Helpers: `graph_to_matrix` and `binarize_by_threshold`, plus `combine='max'` in
  `load_combined_per_type` (`utils/matrix_utilities.py`).
- `utils/network_utilities.py`: removed dead copies of contacts helpers (`build_graph`, `load_*_contacts`, …). The
  live versions are in `contacts_utilities.py`.
- **Experimental, not wired into the CLI**: `services/compare_variants.py` (`CompareVariants(Network)`) ranks network
  metrics against ClinVar pathogenicity (AUROC/AUPRC) by aligning the structure sequence to a RefSeq protein.
  Helpers: `utils/clinvar_query_utilities.py` (NCBI E-utilities) and `utils/metric_comparison_utilities.py`. It
  needs `requests`, `biopython` and `scipy`: `scipy` is a core dependency, `requests` and `biopython` are the optional
  `variants` extra in `pyproject.toml` (`pip install -e .[variants]`). `clinvar_query_utilities.py`
  defines `suggest_refseq_protein` and `fetch_reference_protein` twice. The planned interface is a
  `mypsns compare-variants` subcommand.
- `note.txt` holds the user's personal notes. Its TODO items have been assessed and moved into the list below.

## TODO

Single development TODO list for the project (`README.md` links here). Update it when items are done or found.
Last reviewed 2026-10-26.

### Housekeeping
- [x] Commit the uncommitted work (network weighting/matrices, topology halt, README/CLAUDE.md). The experimental
  ClinVar files are in their own commit.
- [x] Declare the missing dependencies in `pyproject.toml`: `scipy` is now a core dependency; `requests` and
  `biopython` are the optional `variants` extra. (`environment.yml` does not list `requests`/`biopython` yet.)
- [x] `commands/__init__.py`: imports and `__all__` now include `network_cmd` and `guess_ff_path_cmd`.
- [x] `utils/matrix_utilities.py`: persistency is combined by `max` in both weighting modes, energy by `sum`. The
  code, docstring, README "Edge weights" table and this file agree, and it is now the default of
  `load_combined_per_type`.
- [x] `mypsns/commands/README_CLI.md`: updated (layout, subcommands, `network`, service section). Its sibling
  `commands/pyproject_snippet.toml` is a stale copy of the old drop-in and could be deleted.

### Testing and validation
- [ ] Turn `mypsns/prove.py` into a real regression test or script with CLI arguments instead of hard-coded paths.
  Assert SB/HC ≤ ~0.2 % and H-bonds ≤ ~5 % against the PyInteraph test case. Fix the overwritten `rec_grk1` dict.
- [ ] Add a test that the two `NAMED_INTERACTOR_DEFAULTS` copies (`config.py`, `utils/interactors_utilities.py`) are
  identical, or better, keep a single copy. The README used to claim this test existed.
- [ ] Validate the new `network --write-matrices` output: compare `persistency_unweighted_network_adjacency.dat`
  with PyInteraph's macro-IIN on the test case.
- [ ] Validate `q_interactors` HC labels on a representative system (still listed under open issues).

### Features
- [ ] Wire `CompareVariants` into a `mypsns compare-variants` subcommand (config model + `*_cmd.py`). Remove the
  duplicate `suggest_refseq_protein` / `fetch_reference_protein` definitions in `clinvar_query_utilities.py`.
- [ ] Unmatched-residue workflow: document the manual `possible_topologies → found_topologies` fix with screenshots
  (from `note.txt`) and/or add a helper command that promotes a chosen candidate and writes
  `complete_residue_topologies.csv`.
- [ ] PyMOL: colour intra-community edges by community in `community_metrics/*.pml`, keeping cross-community
  edges grey (from `note.txt`). `_pml_edges` currently draws every edge grey.

### Done (assessed from `note.txt`)
- [x] Node spheres centred on CA (CoM for residues without CA), radius scaled by the node metric (`--sphere-scale`).
- [x] Edges as sticks: solid cylinders with weight-scaled radius (`--stick-scale`).
- [x] One `.pml` file per metric (`node_metrics/`) and per community method (`community_metrics/`).
- [x] Only network nodes are drawn; solitary nodes (community `-1`) are skipped.
- [x] `min_bonds_to_polar` note: documented in the README under `interactors` and *Open issues*.
