# Changelog

Changes that matter to users and to run folders made with earlier versions.
Stages skip existing output, so an old run folder keeps its old results until
the affected stages are re-run with `--overwrite`.

## 2026-10-26

### Results change (re-run to get them)

- **Inter-chain H-bonds were missing** from `Hbonds_final.csv`,
  `Hbonds_adjacency.dat` and `combined_adjacency.dat` of every multi-chain run
  (the energy matrices had them). Single-chain runs are unaffected. The
  per-frame `Adjacency/<scheme>/H_bonds/*_merged.csv` are complete: re-run
  `contacts --overwrite`, or rebuild only the matrices from Python with
  `Contacts(topology=…, trj_file=…, outpath=…, interactors_type=[…],
  matrix_type=['persistency']).get_matrices()`.
- **H-bond donor–hydrogen pairing** is now 1.4 Å (`--d-h-cutoff`; MDAnalysis'
  1.2 Å missed Cys S–H, 1.32 Å). Thiol donors now form H-bonds, as in
  PyInteraph: test case 4.56 → 3.93 %, no PyInteraph-only edge left.
- **Hydrophobic centres of mass of terminal residues** now include the
  terminal-patch atoms PyInteraph counts as side chain (e.g. `HA` of an
  N-terminal Met, `HT2` of a COOH C-terminal Leu): the default HC rule selects
  by the new `residue` attribute, and its exclusions add `OT1, OT2`. Runs
  without hydrophobic termini are unchanged. Re-run `interactors --overwrite`
  and `contacts --overwrite`.
- **Unweighted networks** keep values strictly `>` the threshold (was `≥`),
  like PyInteraph's `filter_graph` and like the weighted graphs.
- **Terminal patches** are applied to the bonds of every atom (added atoms are
  bonded to their control atom, deleted ones leave every bonds list), and `NA`
  atom/residue names are no longer read as missing. Re-run from `topology`.
- **`q_interactors` HC clusters** are an exact, deterministic split; earlier
  labels depended on Python's hash seed.

### Added

- `tests/data/`: the PyInteraph 1 test case (structure, trajectory, pruned
  forcefield, PyInteraph graphs and macro-IIN) and the recoverin–GRK1 2018
  PyInteraph graphs, with the configs that reproduce them; presets
  `pyinteraph_testcase`, `rec_grk1`, `rec_grk1_r1` for `mypsns compare`.
- End-to-end regression test on the test case (`pytest -m slow`).
- `mypsns accept-topology` for residues topology recognition cannot match;
  `structure_atom_name` in the recognition tables.
- `mypsns compare` (replaces `mypsns/prove.py`) and the experimental
  `mypsns compare-variants`.
- `network --weighting` and `--write-matrices`; community-coloured edges in
  `community_metrics/*.pml`.
- Flags `--d-h-cutoff`, `--min-residue-charge`, `--da-masses`,
  `--min-num-heavy-ats`.

### Fixed

- `run_config.yaml`: each command now records only its own section(s). Before,
  a stage run from flags rewrote every other section with defaults and dropped
  the forcefield path, so the file could misreport a run (a contacts run at
  SB 3.5 Å read back as 4.5 Å). Check `Adjacency/Contacts.log` (`cutoff=` lines)
  for what an older run really used.
- `contacts` and `mypsns run` never skipped a finished contacts stage (they
  looked for markers that do not exist); they now skip it like the other stages.
- `mypsns report` looked for a `Contacts/` folder that no longer exists; it now
  reports contacts and network per interactor scheme.
- `mypsns run` now passes the configured `named_interactors_rules` and
  `min_residue_charge` to interactors (it used the built-in rules).
- `interactors` and `contacts` no longer require a forcefield path.

### Removed

- `mypsns/prove.py`, `example_network_config.yaml` (merged into
  `example_config.yaml`, which now runs on the test case),
  `mypsns/commands/README_CLI.md` (merged into `README.md`).
