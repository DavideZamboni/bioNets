# CLAUDE.md

myPSNs: Protein Structure Network analysis of MD trajectories, a forcefield-aware re-implementation of
**PyInteraph 1**, validated against it. **`README.md` is the reference**: Part I is the user guide and command
reference, Part II the developer guide (architecture, data flow, extending, interactor rules, testing, validation,
conventions). Read the relevant Part II section before changing behaviour, and update the README with every change
(and `CHANGELOG.md` when results of existing run folders change). This file only adds what the README does not say.

## Commands

```bash
conda activate mypsns                 # or call /home/davide/miniconda3/envs/mypsns/bin/{python,mypsns,pytest}
pytest -m "not slow"                  # ~15 s
pytest                                # + end-to-end PyInteraph regression on tests/data, ~3–5 min
```

## Working rules

- **Never write into `/home/davide/storage/lab_dati/psns_out/`** (the user's runs). Check CLI behaviour on a
  scratch outpath; for a finished run, symlink its `Adjacency/`, `Interactors/`, `forcefield/`,
  `Topology_recognition/` into the scratch folder (every command writes logs and `run_config.yaml` into its
  outpath; `interactors` also rewrites `Topology_recognition/found_topologies/` in place, so copy that one).
- `note.txt` holds the user's personal notes: never commit it.
- Validated defaults (contact cutoffs, `NAMED_INTERACTOR_DEFAULTS`) change only with `pytest` (slow test included)
  passing; dataset-specific settings belong in a preset's `config.yaml`, not in the defaults.
- Comments and docstrings are long and explain *why*; match that.

## Code-level gotchas

- Rules: `resname` is the atom's own forcefield entry (terminal-patch atoms carry the patch name), `residue` the
  residue it belongs to (`residue_name` drops patch rows via `forcefield/final_ff_ters.csv`; set in
  `Interactors.add_inttypes`, not written to disk). HC rules use `residue`, charged groups `resname`.
- `run_config.yaml`: commands call `dump_run_config(cfg, outpath, sections=[...])` with their own section(s);
  `run` passes its four stages. Never dump the whole config from a single stage.
- Stage "done" checks: `contacts_done(cfg)` (every requested scheme has `Adjacency/<scheme>/*_adjacency.dat`).
  Contact files are grouped by folder (`_SUBFOLDER`), never by file-name suffix: inter-chain H-bond files are
  `<c1>_HDon_<c2>_HAcc` / `<c1>_HAcc_<c2>_HDon`.
- `Topology_recognition/found_topologies/*.csv` keep forcefield names in `atom_name`/`bonds` and the structure's
  in `structure_atom_name`; `interactors` renames in `build_wholetop` because contacts joins on
  `(Resindex, atom_name)`. Topology keys come from `topology_key()`, never `.rstrip("_atoms.csv")`.
  `topology --overwrite` removes the other outcome table, i.e. discards acceptances.
- `.tdb` instructions live in the `bonds` column of `final_ff_ters.csv` (`rX`, `-X`, `+X`);
  `is_patch_instruction` separates them from bonds, `has_bond_token` does token (not substring) tests.
- A new `q_interactors` parameter goes in `CHARGE_BASED_PARAMETERS`, `InteractorsConfig` and the
  `interactors`/`run` forwarding (tests check all three).
- `load_combined_per_type(combine=None)` picks the rule from the source (persistency max, energy sum): don't pass
  `combine` unless overriding.
- `api_key` of `CompareVariantsConfig` is `exclude=True`: never dumped. The `variants` extra (`requests`) is not in
  the `mypsns` env; to try `compare-variants` without changing the env, `pip install --target <scratch> requests`
  and set `PYTHONPATH`.

## Lab data (outside the repo)

- `/home/davide/storage/lab_dati/`: `structures/` (`sim.prot.gro`, `traj.xtc`, `rec_grk1_nolip.pdb`,
  `rec_grk1_1micros.xtc`), `forcefields/`, `psns_out/` (user runs, read-only), `risultati_paragone/`.
  `risultati_paragone/pyinter_test/` = the test-case graphs (also in `tests/data/`). `risultati_paragone/rec/*_1micros.dat`
  are **not** graphs of `rec_grk1_1micros.xtc` (no frame window gets below 34 % SB / 53 % HC; provenance unknown);
  only their macro-IIN `rec_grk1_mem_network.dat` was used, to check `network` against `filter_graph`.
- PyInteraph 1 test case source: `~/Desktop/pyinteraph-master/examples/` (`tutorial.sh` has the exact commands).
  PyInteraph 1.1 with `filter_graph` is installed in the `pyint2` conda env.
- rec_grk1 2018 PyInteraph run (source of `tests/data/rec_grk1/`):
  `/home/davide/storage/old/lab_dati/PSNs/complete/LAST/compared_data/ordered/rec/wt_grk1/` (graphs in `network/`,
  `.tpr` and replica xtcs in `top_trj/`). Its settings are in `tests/data/rec_grk1/README.md`.
- `psns_out` runs predate the 2026-10-26 fixes (see `CHANGELOG.md`); `cli_outs/rec_grk1` ran contacts at SB 3.5 Å /
  HC 5.5 Å whatever its `run_config.yaml` says. The user does not need them rebuilt.

## TODO

Open items only (done work is in `CHANGELOG.md` and git history). Last reviewed 2026-10-26.

- [ ] `q_interactors` HC yields single-heavy-atom clusters despite `min_num_heavy_ats=2` (LEU `CB` or `CD2`
  alone, the lone `CB` of ARG/GLN/HSE) where the bonds allow nothing better. Decide whether that is intended (ALA's
  CB is the motivating case).
- [ ] `run_config.yaml` stores paths as given, so a relative path only resolves from the same working directory.
  Decide whether `dump_run_config` should store absolute paths.
- [ ] `compare-variants`: validate on a gene with pathogenic missense variants (RCVRN has none: 18 sites, 17 VUS,
  1 benign). The scoring path was only checked with VUS relabelled as positives, against an independent
  Mann–Whitney U/(n₁n₀).
- [ ] `compare-variants`: positions come from `protein_change` first, the title's `p.` second. For multi-isoform
  genes check that `protein_change` refers to the `--reference` isoform (compare `structure_aa` with `ref_aa` in
  `<gene>_clinvar_sites.csv`).
- [ ] `utils/clinvar_query_utilities.py`: `clinvar_via_eutilities` + `parse_clinvar_summary` are an unused second
  ClinVar fetcher with a richer, location-aware parser. Fold the parser into `fetch_clinvar_variants` (replacing
  `_parse_summary`) or delete both.
- [ ] Optional: screenshots for the README's *When topology recognition stops* (the section shows the terminal
  output as text).
