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

Tests: `pip install -e .[dev]` (or just `pytest`; it is not in `environment.yml`), then `pytest` from the repo root
(`testpaths = ["tests"]`). They cover the matrix/network helpers, `mypsns compare`, the single
`NAMED_INTERACTOR_DEFAULTS`, and a PyInteraph regression test that reads the lab data (paths below; skipped when
missing; `MYPSNS_VALIDATION_ROOT` points elsewhere). The regression test checks the matrices on disk, so after
changing contact code re-run `contacts` on the test case first. Pipeline stages have no tests: check them by
running the CLI on a real run folder. When doing so, write into a scratch `--outpath` with symlinked
`Adjacency/`/`Interactors/` (and `forcefield/`), because every command re-dumps `<outpath>/run_config.yaml`.

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
| Validation | `mypsns/utils/validation_utilities.py`, `commands/compare_cmd.py`, `tests/` | Matrix comparison (`mypsns compare`) and pytest suite |

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
- **`NAMED_INTERACTOR_DEFAULTS` lives only in `config.py`**. `utils/interactors_utilities.py` imports and re-exports
  it, and `tests/test_named_interactor_defaults.py` fails if a second copy appears. `InteractorsConfig` hands out a
  deep copy, so don't mutate the constant itself.
- User-supplied `named_interactors_rules` **replace** the defaults entirely. They are not merged.
- Interactor roles are written to disk by `interactors`. After changing the rules, re-run `interactors --overwrite`
  and then `contacts`. Stages skip when their output already exists.
- HC centres of mass use **forcefield masses**, not MDAnalysis-guessed ones (guessing turns `CAL` into carbon). The
  salt-bridge charge gate (`--min-residue-charge`, 0.8 e) uses forcefield charges.
- Default cutoffs match PyInteraph 1: SB 4.5 Å, HC CoM 5.0 Å, H-bond D–A 4.5 Å with D–H–A ≥ 120°. Don't change them
  without re-running the validation.
- Topology recognition now **halts** with `TopologyAssociationError` when a residue has no perfect forcefield match
  (`services/topology_recognition.py`). For large custom residues, raise `--max-to-check` toward `-1`. The way past a
  halt without changing the forcefield is `mypsns accept-topology` (`accept_topologies` in the same module). It
  renames the candidate's atoms to the structure's names, because contacts joins found_topologies with the
  structure on `(Resindex, atom_name)`; a plain copy of a candidate passes `interactors` and then fails in
  `contacts`. Recognition removes the other outcome table (`complete_…`/`unmatching_…`) left by an earlier run, so
  `topology --overwrite` discards acceptances. The topology key comes from `topology_key()`, never from
  `.rstrip("_atoms.csv")`.
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

Most of the ~4.5 % H-bond gap comes from the detection engine (MDAnalysis `HydrogenBondAnalysis` vs PyInteraph's
geometric test), and no parameter removes it. The exception is the 3 PyInteraph-only edges: all are Cys86 thiol
donors, which MDAnalysis never pairs because its `d_h_cutoff` default (1.2 Å) is shorter than S–H (1.32 Å). With
1.4 Å they match PyInteraph exactly (51.7 / 22.5 / ~2 %). See the TODO. Treat SB/HC drifting above ~0.1 % as a
regression.

- Definitions: an edge is an off-diagonal pair > 0 (PyInteraph's `hb-graph.dat` has 24 intra-residue diagonal
  entries, ignored). The % is `100·Σ|C−T|/ΣT` over the strict lower triangle.
- Tool: `mypsns compare --preset pyinteraph_testcase -o <dir>` (logic in `utils/validation_utilities.py`). It replaced
  the ad-hoc `mypsns/prove.py`, and its mismatch CSVs are byte-identical to prove.py's for both datasets. Presets
  resolve under `$MYPSNS_VALIDATION_ROOT`, default `/home/davide/storage/lab_dati/`. The test case compares
  `risultati_paragone/pyinter_test/{hb,hc,sb}-graph.dat` with
  `psns_out/cli_outs/pyinteraph_testcase/Adjacency/named_interactors/*_adjacency.dat`. `tests/test_pyinteraph_regression.py`
  asserts SB/HC ≤ 0.2 % and H-bonds ≤ 5 % on it.
- The test case ran with `named_interactors` and default cutoffs.
- `--preset rec_grk1` (recoverin–GRK1 vs `risultati_paragone/rec/*_1micros.dat`) is **not** a validation: H-bonds
  33 %, HC 60 %, SB 44 % (shared-edge values correlate, r 0.82–0.93). The PyInteraph settings of that run are
  unknown.
- Macro-IIN (checked 2026-10-26): `risultati_paragone/rec/rec_grk1_mem_network.dat` is exactly PyInteraph 1.1
  `filter_graph -d hb -d hc -d sb -t 20.5` on the rec `*_1micros.dat` graphs (any cutoff in 20.28–20.71 reproduces
  it). `filter_graph` (installed in the `pyint2` conda env) ORs `matrix > cutoff` and zeroes the diagonal. Feeding
  those PyInteraph graphs through `mypsns network --min-persistency 20.5` reproduces the macro-IIN exactly in
  `persistency_unweighted_network_adjacency.dat` (331/331/331). On the test case, against `filter_graph` output,
  the edge counts PyInteraph / myPSNs / shared are 353/377/350 at cutoff 0 and 192/195/190 at 20.5. Every mismatch is
  an H-bond-only pair, so the network step adds no error of its own. myPSNs uses `≥`, filter_graph uses `>`. No
  value sits exactly on a cutoff in this data, so the difference changes nothing here.

## Recent work (committed 2026-10-26)

- `network`: `--weighting {weighted,unweighted,both}` (binary macro-IIN analogue), `--write-matrices` (dumps
  `Network/<inttype>/matrices/*.dat` in Adjacency format via `Network.write_matrices`), threshold semantics
  (weighted `>`, unweighted `≥`). Helpers: `graph_to_matrix` and `binarize_by_threshold`, plus the per-source `combine`
  default in `load_combined_per_type` (`utils/matrix_utilities.py`).
- `utils/network_utilities.py`: removed dead copies of contacts helpers (`build_graph`, `load_*_contacts`, …). The
  live versions are in `contacts_utilities.py`.
- **Experimental `mypsns compare-variants`** (wired 2026-10-26, not yet committed): `services/compare_variants.py`
  (`CompareVariants(Network)`) ranks network metrics against ClinVar pathogenicity (AUROC/AUPRC/Mann–Whitney). It
  aligns a RefSeq protein to the structure chains to map positions. Helpers: `utils/clinvar_query_utilities.py`
  (NCBI E-utilities) and `utils/metric_comparison_utilities.py`. Things that are easy to get wrong:
  - Config: ClinVar/scoring options are in `CompareVariantsConfig` (`compare_variants:`). The graph-rebuild flags
    override `network:`, so `--config <outpath>/run_config.yaml` ranks the same graphs `network` built.
  - `api_key` is `Field(exclude=True)` and is never dumped. Use `--api-key` or `NCBI_API_KEY`.
  - It does not call `run_analysis`: centralities are computed in memory, so `Network/<inttype>/csv/` is not touched.
  - Only residues aligned to the reference are scored. Partner chains and ions are never negatives.
  - The `[chrpos]` window is GRCh38.
  - Outputs: `Network/<inttype>/variants/<gene>_*.csv`.
  - `requests`/`biopython` are the `variants` extra. The service raises an `ImportError` naming it at import time;
    the command prints it and exits 3, also under `--dry-run`. The `mypsns` conda env has `biopython` but not
    `requests`.
- `note.txt` holds the user's personal notes. Its TODO items have been assessed and moved into the list below.

## TODO

Single development TODO list for the project (`README.md` links here). Update it when items are done or found.
Last reviewed 2026-10-26.

### Housekeeping
- [ ] Delete `mypsns/commands/pyproject_snippet.toml`: a stale copy of the old CLI drop-in, referenced nowhere.
- [ ] `environment.yml`: add `requests` and `biopython` (the `variants` extra in `pyproject.toml` has them, the conda
  env file does not). `scipy` is already there.

### Testing and validation
- [ ] H-bonds: pass `d_h_cutoff` (~1.4 Å) to `HBA(...)` in `run_batch_hba` (`contacts_utilities.py`) so Cys SG–HG1
  is paired. Verified in isolation: Cys86→Asp118/Pro119 go from 0 % to PyInteraph's 51.7/22.5 %. This changes
  validated output, so re-run `contacts` on the test case and `pytest` afterwards. Expect H-bond PyInteraph-only
  edges 3 → 0 and the gap ~4.6 → ~4.0 %. Awaiting the user's go-ahead.
- [ ] `q_interactors` `HC_cluster` is not reproducible (found 2026-10-26, cause confirmed): `add_hc_atoms` builds
  `hcs = list({a[0] for a in true_hcs})`. The set order of the atom names depends on `PYTHONHASHSEED`, and so does
  the greedy split in `get_constrained_conncomps`/`split`. Over 8 seeds, PHE/TRP/LEU each give 2 partitions (TYR/ILE/MET
  are stable). The stored test-case `Interactors/` PHE partition differs from a fresh run. Fix: sort `hcs`, and make the
  split's tie-break deterministic. `named_interactors` is unaffected. Awaiting the user's go-ahead.
- [ ] `q_interactors` HC also yields single-heavy-atom clusters despite `min_num_heavy_ats=2` (LEU `CB` or `CD2`
  alone, and the lone `CB` of ARG/GLN/HSE). Decide whether that is intended (ALA's CB is the motivating case).
- [ ] `get_charge_based_interactors`: its `defaults` dict lacks `min_bonds_to_polar`, so calling it without that
  kwarg raises `KeyError`. The CLI always passes it, so the bug is latent.
- [ ] Unweighted network threshold is `≥`; PyInteraph's `filter_graph` macro-IIN uses `>`. Pick one (switching to `>`
  gives strict parity) and update the README "Edge weights" table and the `binarize_by_threshold` docstring/tests.
- [ ] `network` cannot run from flags alone: `PipelineConfig` requires `forcefield`, which `network` neither uses nor
  exposes (`--forcefield` missing). Make it optional for that command or add the flag.
- [ ] `rec_grk1` vs its PyInteraph graphs is 33/60/44 % (H-bonds/HC/SB). Recover that PyInteraph run's settings
  (frames, selection, cutoffs) before using it as a second validation dataset.
- [ ] Recognition's terminal patching (`get_ter_ffres`) does not rename bond references (found 2026-10-26): the exact
  `NHT3+` Met match in `psns_out/{wt,M654V}_noh2o` has `N` bonded to `H1 H2 H3` while its atoms are `HT1–HT3`, so the
  N–HT bonds are lost. `dangling_bonds()` detects it. Check the effect on N-terminal H-bond donors.

### Features
- [ ] `compare-variants`: run it on a gene that has pathogenic missense sites in ClinVar. On rec_grk1, RCVRN has none
  (18 missense sites: 17 VUS, 1 benign), so every ranking there is empty. The scoring path was only checked with
  VUS relabelled as positives, in a scratch script whose AUROC matched an independent Mann–Whitney U/(n₁n₀).
- [ ] `compare-variants`: positions are parsed from `protein_change` first and from the title's `p.` second. For
  multi-isoform genes, check that `protein_change` refers to the `--reference` isoform (the title's `p.` is tied to
  the NM_ it names). Comparing `structure_aa` with `ref_aa` in `<gene>_clinvar_sites.csv` exposes a mismatch.
- [ ] `clinvar_query_utilities.py`: `clinvar_via_eutilities` + `parse_clinvar_summary` are an unused second ClinVar
  fetcher with a richer, location-aware parser. Fold the parser into `fetch_clinvar_variants` (replacing
  `_parse_summary`), or delete both.
- [ ] README *Unmatched residues*: add the four screenshots at the `[Screenshot n — to add]` placeholders
  (suggested files under `docs/img/`).

### Done
- [x] Commit the network weighting/matrices, topology halt and README/CLAUDE.md work (the experimental ClinVar files
  are in their own commit).
- [x] `mypsns compare-variants` wired: `CompareVariantsConfig`, `commands/compare_variants_cmd.py`, README section.
  The duplicate `suggest_refseq_protein` / `fetch_reference_protein` and the unused `EUTILS` in
  `metric_comparison_utilities.py` were removed. Fixed bugs: the missing `map_reference_positions` method, dict
  iteration, the centrality/ClinVar `label` column clash, NaN edge features for isolated nodes, and `score()`
  crashing when there are no positives. Tested on rec_grk1 (RCVRN, NP_002894.1). Not committed.
- [x] `pyproject.toml`: `scipy` is a core dependency; `requests` and `biopython` are the optional `variants` extra.
- [x] `commands/__init__.py`: imports and `__all__` include `network_cmd` and `guess_ff_path_cmd`.
- [x] Combine rule: persistency by `max` (both weighting modes), energy by `sum`. Code, docstring, README "Edge
  weights" table and this file agree; it is the default of `load_combined_per_type`.
- [x] `mypsns/commands/README_CLI.md` updated (layout, subcommands, `network`, service section).
- [x] `mypsns accept-topology` (2026-10-26, not committed): lists the candidates of unmatched residues (marking
  byte-identical ones) and promotes `TOPOLOGY=<n>` or `TOPOLOGY=<csv>` to `found_topologies/`. It renames
  differently named atoms by mass-checked positional pairing and writes `complete_residue_topologies.csv` once all
  are resolved. The halt message prints the commands. README *Unmatched residues* has the step-by-step workflow.
  Validated end to end on the PyInteraph test case with Gly `HA1`→`HA3` in the structure: after `GLY=0`, the SB/HC/H-bond
  adjacency and energy matrices are byte-identical to the unmodified run. Also fixed on the way: standalone
  `mypsns topology` crashed (`outlog=` kwarg), the halt reported the recognition score as "0/N atom names
  matched", and the halt printed a traceback.
- [x] PyMOL: `community_metrics/*.pml` colour intra-community edges with the community colour (`_community_colour`,
  shared with the spheres); cross-community edges stay `gray60`. `node_metrics/*.pml` are byte-identical to before.
  Checked by regenerating rec_grk1 renders from their own `csv/` files (HEAD code reproduces the existing `.pml`
  exactly).
- [x] Testing and validation (2026-10-26, not committed):
  - `mypsns/prove.py` → `mypsns compare` + `utils/validation_utilities.py`. It takes paths, presets and `TYPE=PATH`
    pairs, writes CSVs to `-o`, and `--max-pct-diff` sets the exit code. Both legacy datasets reproduce prove.py's
    CSVs byte for byte.
  - pytest suite (42 tests): `NAMED_INTERACTOR_DEFAULTS` is now a single copy in `config.py` (also fixed the config
    default sharing inner lists with the constant); `binarize_by_threshold`, `graph_to_matrix`,
    `load_combined_per_type`; the compare logic/CLI; and the PyInteraph regression (skips without data, fails on a
    1 % HC drift). `tests/test_basic.py` (a failing placeholder) was removed.
  - `network --weighting both --write-matrices` validated against PyInteraph macro-IINs (see *Validation*).
  - `q_interactors` HC on the test case: same HC atoms as `named` for ALA/VAL/LEU/ILE/PHE/PRO/MET, but also Tyr ring
    + CB, Thr CG2, Cys CB/SG and the CB of Arg/Gln/HSE; Trp loses CD1/CE2/NE1. Hydrogens differ: `named` also
    keeps backbone HN/HA, as PyInteraph does. HC contacts become 207 edges vs 79 (76 shared). 84 of the 131 extra
    edges involve the added residue types; the other 47 come from per-cluster CoMs between PyInteraph HC
    residues. Follow-ups are listed above.

Assessed from `note.txt`:
- [x] Node spheres centred on CA (CoM for residues without CA), radius scaled by the node metric (`--sphere-scale`).
- [x] Edges as sticks: solid cylinders with weight-scaled radius (`--stick-scale`).
- [x] One `.pml` file per metric (`node_metrics/`) and per community method (`community_metrics/`).
- [x] Only network nodes are drawn; solitary nodes (community `-1`) are skipped.
- [x] `min_bonds_to_polar` note: documented in the README under `interactors` and *Open issues*.
