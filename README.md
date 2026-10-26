# myPSNs

Protein Structure Network analysis from molecular-dynamics trajectories.

myPSNs builds residue–residue interaction networks (PSNs) from a GROMACS
forcefield, a topology, and an MD trajectory. It identifies hydrogen bonds,
salt bridges, and hydrophobic clusters per frame; aggregates them into
persistency- and energy-weighted contact matrices; and constructs graphs with
centrality, community detection, and PyMOL renderings.

The default `named_interactors` role set reproduces **PyInteraph 1**: on a CHARMM
test trajectory, salt bridges and hydrophobic clusters match to ~0.1 %, H-bonds
to ~4.5 % (an engine-limited floor, not a parameter — see *Interactor roles &
PyInteraph reproduction*). The rules remain fully customisable for non-standard
residues, ions, and lipids.

## Pipeline:

ff-reader → topology → interactors → contacts → network

The first four stages produce the contact data; network consumes it. network
is not part of mypsns run — invoke it separately once contacts has
finished.

## Setup

cd <proj_root>/myPSN
conda env create -f environment.yml
conda activate mypsns
pip install -e .

This installs the mypsns console script. Verify with:

mypsns --version
mypsns --help

## Quick start

Run the contact pipeline (one of the two ways below), then build the network.

With a config file (recommended)
## Run

mypsns run     --config example_config.yaml
mypsns network --config example_network_config.yaml

With explicit flags (minimal arguments)

mypsns run \
  --topology   system.gro \
  --forcefield charmm36.ff \
  --trajectory traj.xtc \
  --outpath    run_out

mypsns network \
  --topology system.gro \
  --outpath  run_out \
  --edge-weight-source persistency \
  --edge-weight-source energy

Any config value can be overridden on the command line; flags win over the YAML.
Use --dry-run on any stage to validate inputs and print the plan without
executing.

## CLI commands

| Command              | Purpose
|----------------------|----------------------
| `mypsns run`         | Full pipeline (ff-reader → topology → interactors → contacts).
| `mypsns ff-reader`   | Parse a GROMACS forcefield folder.
| `mypsns topology`    | Match topology residues against the forcefield.
| `mypsns interactors` | Assign interaction roles to topology atoms.
| `mypsns contacts`    | Compute contacts across an MD trajectory.
| `mypsns network`     | Build residue-interaction networks and render PyMOL sessions.
| `mypsns guess_ff`    | Suggest candidate forcefield sub-folders under a path.
| `mypsns validate`    | Validate a config and its referenced files.
| `mypsns report`      | Summarise the outputs in a run folder.

Run `mypsns <command> --help` for the full option list of any stage.

## General notes

## Project structure

```
mypsns/
├── cli.py                  console entry point (Typer)
├── config.py               Pydantic config models
├── commands/               one module per CLI subcommand
├── services/
│   ├── ff_reader.py
│   ├── topology_recognition.py
│   ├── interactors.py
│   ├── contacts.py         contact computation + matrices
│   └── network.py          PSN construction + PyMOL
├── utils/                  per-stage helper functions
└── resources/              bundled lookup tables (e.g. sol.csv)
```

### Modules

- **`cli.py`** — application entry point.
  Configures automatically CLI reading function's parameter names and types (no argparse)
  callback: runs first, reads the global flags (-v, -q, --log-file) and configures logging
  Calls specific commands found in commands/*_cmd.py (through app.command).

- **`config.py`**
  source of truth for what parameters the pipeline accepts, what their types are, what their defaults are, and what constraints apply.
  Any invalid value is caught HERE, before any pipeline code runs.

  Defines global DEFAULT_NA_VALUES (for pandas.read_csv, omitting NA to preserve sodium atom names) and a Pydantic model for each of the main modules.

  Orchestrates whole application setup (flags > .yml > defaults) and writes complete yaml file.
  Also holds `NAMED_INTERACTOR_DEFAULTS` — the PyInteraph-matching `named_interactors` rule
  set used as the default for `InteractorsConfig.named_interactors_rules` (mirrored in
  `utils/interactors_utilities.py` — **keep the two copies identical by hand**; no test checks it yet).

## Services

The `mypsns/services/` folder holds the modules that drive the pipeline.
Expected execution order:

1. `ff_reader` — parse the GROMACS forcefield
2. `topology_recognition` — match topology residues to forcefield entries
3. `interactors` — reconstruct whole-system topology + assign interaction roles
4. `contacts` — compute contacts over a trajectory and build matrices
5. `network` — build residue-interaction networks (standalone; run after `contacts`)

Stages 1–4 are chained by `mypsns run`; `network` is invoked separately.

A note on `na_values`: every stage accepts an optional `na_values` list that is
forwarded to every `pandas.read_csv`. The shared default deliberately **excludes
`"NA"`** so sodium atoms/residues named `NA` are not parsed as missing. There is
no `main.py`; the entry point is the `mypsns` console script (`mypsns/cli.py`),
and per-stage options are set via the config file or CLI flags.

---

### `ff_reader.py`

GROMACS forcefield parser.

**Parameters**
- required:
  - `inpath`: forcefield folder (the `.ff` directory)
  - `outpath`: destination folder (e.g. `charmm36_apr2019_vm` from `charmm36-apr2019_vm.ff`)
- optional:
  - `na_values`: NA-string list (see note above)
  - `nonbonded_name` (default `"ffnonbonded.itp"`): file holding nonbonded parameters
  - `additional_residues` (default `True`): merge extra residues bundled in `resources/`

**What it does**

1. Reads every residue defined in the `.rtp` files. Per atom it records
   `rtp_file`, `resname`, `atom_name`, `atom_type`, `charge`, `chargegroup`
   (integer grouping whose total charge sums to an integer — a GROMACS device
   for Coulomb cutoffs when PME is off), and from `ffnonbonded.itp` the
   `sigma`, `epsilon`, and `mass`.
2. Reads termini definitions from `.n.tdb` / `.c.tdb`. **`.rtp` and `.tdb`
   names must correspond** (e.g. `merged.rtp` ↔ `merged.n.tdb` / `merged.c.tdb`).
   Bond-entry conventions: atoms prefixed `-` (present in every bond entry of a
   terminus) are **found and removed** from the residue; atoms prefixed `r` are
   **found and replaced** with the terminus-defined atom, bonds adjusted.
3. Reads σ and ε from the `.itp` named by `nonbonded_name`.
4. Writes `final_ff_residues.csv` (`residues.csv` ⋈ `atomtypes_se.csv` on
   `atom_type`, adding `mass`, `sigma`, `epsilon`) and `final_ff_ters.csv`
   (the ters table, similarly joined). `atomtypes_se.csv` (atom type → σ, ε) is
   also produced and is later consumed by `contacts` for energy calculations.
   The `mass` column here is the **forcefield mass**, later reused by `contacts`
   for hydrophobic-cluster centres of mass (see *Interactor roles*).

**Expected inputs**: `*.rtp`, `*.tdb` (`.n.tdb` / `.c.tdb`, names matched to the
`.rtp`), and the nonbonded `.itp`.

---

### `topology_recognition.py`

Matches each residue in the input topology to a forcefield entry from
`ff_reader`. The same residue name appearing with different atom sets
(e.g. terminal vs internal) is treated as distinct residue topologies.

**Parameters**
- required: `topology` (PDB/GRO/…), `outpath`, `ff_path` (the `ff_reader` output)
- optional: `na_values`; `outlog` (log folder, defaults to `outpath`)

The constructor calls `parse_initial_topology_info`, writing
`initial_topology_info.csv` — a per-atom table with whichever of `chainID`,
`resindex`, `resid`, `resname`, `index`, `id`, `name`, `type`, `element`,
`mass`, `charge` the topology format supports (columns that raise
`NoDataError` for every atom are omitted). Special case: if the format supplies
`type` but not `element` and the longest type string is ≤2 chars, `type` is
copied into `element` (handles PDBs where atomic symbols land in the type column).

**Methods**

- `get_restypes(verbose=False)` — group atoms into unique residue topologies by
  exact atom-name set. Writes `<outpath>/resid_topologies/<resname>[_<n>]_rixs.txt`
  (Resindex lists) and `<resname>[_<n>]_atoms.csv` (atom names + masses), plus
  `different_residue_topologies.csv` indexing them. The `_<n>` suffix is dropped
  for the first variant of a resname (`ALA_atoms.csv`, then `ALA_1_atoms.csv`, …).
- `get_ff_parameters(verbose=False, maxToCheck=-1)` — associate each topology
  residue with a forcefield entry. Candidates are sorted by name similarity
  (`difflib.SequenceMatcher`); `maxToCheck` caps how many are inspected
  (`-1` = all). Two passes: (1) exact atom-name-set match for candidates with the
  right atom count; (2) on failure, apply each `.tdb` ter substitution and
  re-check. Outcomes split into `found_residues` (exact, with/without ter subs)
  and `best_choices` (ranked by atom overlap, then mass agreement on unmatched
  atoms). Writes `found_topologies/<resname>_ff_atoms.csv` for matches,
  `possible_topologies/<resname>_uncertain_atoms_<n>.csv` for candidates, and
  `complete_residue_topologies.csv` (or `unmatching_residue_topologies.csv` if
  any residue is unmatched).

  If **any** residue topology is left unmatched, the stage logs which residues
  failed (and how many atom names matched the closest candidate) and raises
  `TopologyAssociationError`, halting the pipeline before `interactors` (which
  would otherwise fail later with a confusing missing-file error).

  > Non-standard residues (e.g. N-myristoylated Gly `GLM`, a metal ion `CAL`)
  > must have a forcefield entry to be matched. If recognition fails for one,
  > raise `--max-to-check` toward `-1` (inspect every ff entry) — large custom
  > residues can sort far down the name-similarity list.

---

### `interactors.py`

Bridges `topology_recognition` and `contacts`: reconstructs full topology
(chains, polymer back-bonds, ter chargegroups), deduplicates residues into
unique chemical variants, and assigns interaction roles.

**Parameters**
- required: `topology` (same file as the rest of the pipeline);
  `general_outpath` (must already contain `Topology_recognition/`)
- optional: `na_values`; `outlog` (defaults to `general_outpath`)

**Methods**

- `get_complete_topology()` — fills bidirectional bond info per residue type;
  concatenates per-residue tables in memory; detects chain boundaries from
  `+`-bond connectivity (polymers) and contiguous `rtp_file` runs (non-polymers);
  stamps `Chain` and `Entity` (`P` polymer / `M` monomer); adds polymer
  back-bonds (`-X` references); adjusts ter chargegroups for terminal polymer
  residues; deduplicates residues by atom-level chemistry into unique variants.
  Outputs to `<general_outpath>/Interactors/`:
  - `found_topologies/<variant_id>.csv` — one CSV per unique residue chemistry
    (`<resname>`, then `<resname>_2`, `<resname>_3`, … for subsequent variants).
    Columns: `rtp_file, resname, atom_name, atom_type, charge, chargegroup, bonds, mass`
    (interaction-role columns are appended later by `add_inttypes`).
  - `residue_index.csv` — one row per residue instance:
    `Resindex, resnum, Chain, Entity, variant_id, topology_resname`
    (`topology_resname` is the residue's name in the **input topology**, e.g.
    `HIS`, as opposed to the forcefield/rtp name carried by `variant_id`, e.g.
    `HSE`; it is what selection strings prefer so they match the structure).
  - `atom_index_map.csv` — one row per atom: `Resindex, atom_name, atom_index`.
- `add_inttypes(**kwargs)` — appends interaction-role columns to each variant
  CSV in place. Each variant is analysed in isolation; back-bonds remain in the
  bonds strings but BFS terminates at the variant boundary.
  - **`q_interactors`** (charge-based, PIE-style): H-bond `HAcc`/`HDon`/`HBond`
    from partial charges + connectivity (`max_DA_charge`=−0.3, `min_H_charge`=0.3,
    `DA_masses`=[14,16,19]); salt-bridge `Acidic`/`Basic` from recursive local
    charge (`min_abs_charge`=0.75); `HC` flag + `HC_cluster` label for atoms with
    `|charge| ≤ max_abs_charge` (0.27) and BFS distance `> min_bonds_to_charged`
    (3) from any charged atom (`max_num_heavy_ats`=3, `min_bonds_to_polar`=1).
  - **`named_interactors`** (PyInteraph-compatible, name-based; **default**):
    roles assigned by AND-combined `(attribute, [values])` / `not <attribute>`
    rules — see *Interactor roles & PyInteraph reproduction* for the grammar and
    the default rule set. After role assignment, donor-bound H atoms are found by
    walking from `HDon` atoms via the bonds column and tagged `HBond`; **their
    `HC` tag is preserved** (the amide `HN` and Trp `HE1` are both H-bond
    hydrogens and members of the hydrophobic side-chain CoM, as in PyInteraph).
    The salt-bridge **charge gate** (`min_residue_charge`) is then applied to
    `Basic`/`Acidic` groups.

  Note: at variant level there is no cross-residue context, so the previous
  3-residue BFS widening for `q_interactors` HC classification no longer applies.
  Differences are negligible for the standard amino-acid set, but verify HC
  labels on a representative system if they look off.

---

## Interactor roles & PyInteraph reproduction

The `named_interactors` scheme assigns each atom one or more **roles** that
`contacts` turns into selections: `HC` (hydrophobic-cluster member),
`HDon`/`HAcc` (H-bond donor/acceptor heavy atoms), `HBond` (a donor-bound H),
`Basic`/`Acidic` (charged-group members). The default rule set reproduces
PyInteraph 1's `hydrogen_bonds.ini` / `charged_groups.ini` / `hc_reslist`, so a
run with **no extra config** matches PyInteraph 1 (SB ~0.07 %, HC ~0.11 %,
H-bonds ~4.5 %).

### Rule grammar

`named_interactors_rules` maps each role to a **list of rules**. Within a rule
all conditions hold (AND); across a role's rules any match qualifies (OR).
Conditions are `attribute → value-list`; prefix an attribute with `not ` to
negate.

```yaml
interactors:
  named_interactors_rules:
    HC:
    - resname: [ALA, VAL, LEU, ILE, PHE, PRO, TRP, MET]   # AND
      not atom_name: [N, CA, C, O, H, H1, H2, H3, O1, O2, OXT]
    Basic:
    - resname: [ARG]
      atom_name: [CD, NE, CZ, NH1, HH11, HH12, NH2, HH21, HH22]
    - resname: [LYS, LYSH]                                  # OR (second rule)
      atom_name: [CE, NZ, HZ1, HZ2, HZ3]
```

Supplying `named_interactors_rules` **fully replaces** the defaults (not merged),
so copy the whole block and edit it. Roles are independent: an atom may be both
`HC` and `HAcc` (Met `SD`, Trp `NE1`) or both `HC` and `HBond` (the backbone
amide `HN`); each contact type filters on its own role, so there is no
double-counting.

### The salt-bridge charge gate

`Basic`/`Acidic` rules tag *candidate* charged groups; the `--min-residue-charge`
gate (default **0.8 e**, `0` disables) then keeps only groups whose net
forcefield charge over the `(residue, forcefield-resname)` group clears the
threshold. This makes salt bridges protonation-state–aware without per-residue
hacks: neutral His (HSE/HSD, net ≈ 0) is filtered out while protonated His
(HSP, +1) and metal cations (e.g. Ca²⁺, +2) pass.

### Hydrophobic-cluster centres of mass

HC contacts use the **centre of mass** (not geometry — geometry was tested and
is worse) of the side-chain atoms, weighted by **forcefield masses** loaded from
`Interactors/found_topologies/` (matching PyInteraph, and immune to MDAnalysis's
element-from-name mis-guesses). The amide `HN` and polar side-chain atoms that
are also donors/acceptors (`SD`, `NE1`, `HE1`) are **kept** in the COM, as
PyInteraph does.

### H-bonds are engine-limited

myPSNs uses MDAnalysis `HydrogenBondAnalysis` (keeps a bond when the D–H–A angle
is **≥** the cutoff; 180° = linear). PyInteraph 1 uses its own geometric test
with a different angle convention and measures the H→acceptor distance (its
3.5 Å ≈ 4.5 Å donor–acceptor here). With `angle 120` / `d_a 4.5` the **edge set**
matches (≈279/282 on the test set); the residual ~4.5 % is small per-frame
persistence scatter (~1.8 pts/edge) from the different detectors and **cannot be
tuned away**. SB and HC, by contrast, reach ~0.1 %.

### Custom / non-standard residues

Because rules are name-based, a non-standard residue just needs its own rule.
Two worked examples (recoverin–GRK1):

```yaml
HC:
- resname: [ALA, VAL, LEU, ILE, PHE, PRO, TRP, MET]
  not atom_name: [N, CA, C, O, H, H1, H2, H3, O1, O2, OXT]
# N-myristoylated Gly: keep the myristoyl tail C2–C14, drop the glycine
# backbone and the amide-linker carbonyl (C1, O1).
- resname: [GLM]
  not atom_name: [N, CA, C, O, HN, HA1, HA2, C1, O1]
Basic:
# … standard Basic rules …
# Ca²⁺ ion: single-atom cation; its +2 forcefield charge clears the gate,
# so it pairs with Asp/Glu carboxylates as a salt bridge.
- resname: [CAL]
  atom_name: [CAL]
```

> **Changing rules requires re-running `interactors` with `--overwrite`.**
> Roles are written to disk by the `interactors` stage (into
> `found_topologies/*.csv`); `contacts` only reads them back. The stage **skips**
> if its output already exists and `--overwrite` is not given, so after editing
> `named_interactors_rules` (or any interactor option) you must re-run
> `mypsns interactors … --overwrite`, then `contacts`. Re-running `contacts`
> alone keeps the old roles. Quick check: the new `Contacts.log` selections
> should reflect the change (e.g. a Met `HC` selection containing `SD` and `HN`).

## `contacts`

Computes per-frame contacts and aggregates them into matrices. Reads
`<outpath>/Interactors/`, writes `<outpath>/Adjacency/`.

### Key options

| Flag                      | Default                               | Meaning
|---                        |---                                    |---
| `--nc`                    | 2                                     | Parallel worker processes.
| `--ram-thresh`            | 80                                    | RAM-usage %% that triggers a flush
| `--nframes`               | 500                                   | Initial H-bond batch size (grows adaptively).
| `--beg-frm` / `--end-frm` | 0 / -1                                | Frame range (`-1` = last frame).
| `--interactors-type`      | `named_interactors`                   | Which role columns to consult. Repeatable.
| `--target`                | `HBonds`, `Salt_bridges`, `HCs`       | Contact types to compute. Repeatable.
| `--subject-chain`         | detected polymers                     | Chains that are the focus of contacts. Repeatable.
| `--surrounding-chain`     | all chains                            | Candidate partner chains. Repeatable.
| `--d-a-dist`              | 4.5 Å                                 | H-bond donor–acceptor distance cutoff.
| `--angle`                 | 120°                                  | H-bond D–H–A angle cutoff (≥ cutoff kept; 180° = linear).
| `--sb-cutoff`             | 4.5 Å                                 | Salt-bridge atom-pair distance cutoff.
| `--hc-cutoff`             | 5.0 Å                                 | Hydrophobic-cluster centre-of-mass distance cutoff.
| `--min-residue-charge`    | 0.8 e                                 | Net-charge gate for `Basic`/`Acidic` groups (named scheme). 0 disables.
| `--matrix-type`           | `persistency`, `energy`               | Which matrices to build after contacts. Repeatable.
| `--energy-use-abs / --no-energy-use-abs` | — | Energy → strength %% mapping: `\|E\|` vs attractive-only `max(0, -E)`.

Defaults for `--d-a-dist` / `--angle` / `--sb-cutoff` / `--hc-cutoff` are the
**PyInteraph 1–matching** values. (Earlier releases defaulted to `150°` / `3.5 Å`
SB / `5.5 Å` HC — override the flags for the old behaviour.) These are global
geometric cutoffs, independent of the interactor scheme.

Two interactor schemas can be requested independently:

- **`named_interactors`** (default) — PyInteraph 1–compatible, atom-name-based
  role rules; **validated** to ~0.1 % (SB, HC) against PyInteraph 1.
- **`q_interactors`** — charge- and connectivity-based roles (PIE-style).
  Currently the less-validated of the two.

### Selection strings

Atom selections are emitted in a compact, human-readable form grouped by
`(resname, atom-name set)`, e.g.

(chainID A) and (
    (resname ASP and resid 81 109 116 and name OD1 OD2)
    or (resname GLU and resid 11 14 16 and name OE1 OE2)
)

This assumes **atom names are unique within a residue** (true for standard
protein/lipid topologies). If a residue has duplicate atom names, rename them in
the input PDB / forcefield first.

### Workflow

1. Reassemble a flat per-atom interactors table in memory from the three files
   written by `interactors` (`found_topologies/<variant>.csv`,
   `residue_index.csv`, `atom_index_map.csv`).
2. Write a corrected-topology PDB (`<name>_correct.pdb`) with `Resindex` stamped
   as `resid` and atom IDs aligned to `atom_index`. All subsequent MDAnalysis
   operations use this topology.
3. Split frames across `--nc` worker processes (PBC-aware distance kernels via
   `capped_distance` / `self_capped_distance`).
4. Per chain combination and contact type, dispatch workers; each streams its
   frame chunk to partial CSVs `<name>_upto_<lastframe>.csv`, flushing whenever
   RAM hits the ceiling. HC workers overwrite the Universe's masses with the
   forcefield masses (from the interactors table) before computing each
   side-chain centre of mass.
5. Merge partials into `<name>_merged.csv` per chain combination.
6. Aggregate merged CSVs into per-residue-pair tables `SB_final.csv`,
   `Hbonds_final.csv`, `HC_final.csv`
   (`res_1, res_2, n_frames_contact, mean_dist, chain_1, chain_2, label_1, label_2`).
7. **Persistency matrices** — `{SB,Hbonds,HC}_adjacency.dat` (entry `[i,j]` is
   the contact persistency in %) plus `combined_adjacency.dat` (frame-union).
8. **Energy matrices** — `{SB,Hbonds,HC}_{LJ,Coulomb,total}_energy.dat`
   (mean-per-frame LJ + Coulomb, kJ/mol). Energy is computed with constant
   memory: per-frame rows stream to `<name>_energy.csv` companions, never
   accumulating the whole trajectory in RAM.
9. **Threshold analysis** — for every non-empty subset of the requested types,
   scan the edge threshold and record the largest connected component +
   degree distribution (elbow plots), to help pick a downstream cutoff.

### Output layout

```
<outpath>/Adjacency/
├── <name>_correct.pdb
├── named/                         (and/or q/, per --interactors-type)
│   ├── Salt_bridges/
│   │   ├── *_upto_*.csv           (partial, per-worker)
│   │   ├── *_merged.csv           (merged per chain-pair)
│   │   └── *_sb_energy.csv        (per-frame energy stream)
│   ├── H_bonds/   ( … *_Hbonds_energy.csv )
│   ├── HCs/       ( … *_HC_energy.csv )
│   ├── adjacency_labels.dat
│   ├── SB_final.csv, Hbonds_final.csv, HC_final.csv
│   ├── {SB,Hbonds,HC}_adjacency.dat        (persistency %)
│   ├── combined_adjacency.dat
│   ├── {SB,Hbonds,HC}_{LJ,Coulomb,total}_energy.dat
│   ├── all_total_energy.dat                (element-wise sum across types)
│   ├── graph_report/        (persistency elbow + degree CSVs/plots)
│   └── energy_graph_report/ (energy strength elbow + distributions)
```


## `network`

Builds residue-interaction graphs from the contacts output and renders PyMOL
sessions. Reads `<outpath>/Adjacency/` and `<outpath>/Interactors/`, writes
`<outpath>/Network/`. Standalone — run after `contacts`.

### Key options

| Flag | Default | Meaning
|---|---|---
| `--interactors-type` | `q_interactors`, `named_interactors` | Adjacency sub-folders to consume. Repeatable.
| `--target` | `HBonds`, `Salt_bridges`, `HCs` | Contact types contributing to the graph. Repeatable.
| `--edge-weight-source` | `persistency`, `energy` | Build one network per source. Repeatable. |
| `--energy-use-abs / --no-energy-use-abs` | off | Energy → weight: `\|E\|` vs attractive-only `max(0, -E)`. |
| `--weighting` | `both` | `weighted` (edge weight = value), `unweighted` (edge weight = 1 at/above threshold — binary macro-IIN analogue), or `both`. |
| `--node-metric` | `betweenness_centrality` | Centrality driving CA-sphere scaling. Repeatable. |
| `--community` | `louvain` | Community algorithm. `louvain`, `greedy_modularity`, `label_propagation`, `connected_components`. Repeatable. |
| `--community-seed` | 42 | Seed for stochastic community algorithms. |
| `--community-min-members` | 2 | Dissolve communities smaller than this (members → solitary, excluded from render/count). 1 keeps singletons. |
| `--min-persistency` | 0.0 | Persistency threshold (%). Weighted graphs keep edges `>` it; unweighted keep edges `≥` it. |
| `--min-energy` | 0.0 | Energy threshold (% of strongest edge). Weighted `>`, unweighted `≥`. |
| `--min-degree` | 0 | Drop nodes below this degree. |
| `--sphere-scale` | 2.0 | Max CA-sphere vdw radius (Å). |
| `--stick-scale` | 0.5 | Max edge stick radius (Å). |
| `--write-matrices / --no-write-matrices` | on | Dump each network's adjacency matrices (`.dat`, Adjacency format) under `Network/<inttype>/matrices/`. |
| `--write-graphml / --no-write-graphml` | off | Also write one GraphML per static graph. |
| `--shortest-paths / --no-shortest-paths` | off | All-pairs weighted shortest paths per graph. |
| `--window-size` | — | If set, also build per-window energy networks and track community evolution. |
| `--window-stride` | = window-size | Step between windows (non-overlapping by default). |
| `--frame-start` / `--frame-stop` | — | Restrict to `[start, stop)`. **When set, static networks are re-aggregated from the per-frame energy CSVs over that range** instead of the pre-aggregated matrices.

### Edge weights

- **persistency** — from `{SB,Hbonds,HC}_adjacency.dat` (contact frequency, %).
- **energy** — from `{SB,Hbonds,HC}_total_energy.dat`. Signed energy is mapped
  to a non-negative weight: attractive-only `max(0, -E)` by default, or `|E|`
  with `--energy-use-abs`, then normalised to "% of the strongest edge".

How the requested targets are combined into one matrix per source:

| | combined value | weighted edge | unweighted edge |
|---|---|---|---|
| persistency | element-wise **max** across targets (strongest single type, ≤ 100 %; union of per-type edges) | value, kept if `> --min-persistency` | 1 where value `≥ --min-persistency` |
| energy | **sum** of per-type attractive (or `\|E\|`) weights, normalised to % of strongest | value, kept if `> --min-energy` | 1 where value `≥ --min-energy` |

One graph is built per `(interactors_type × edge_weight_source × weighting)`.
Weighted graphs are labelled `<source>`, unweighted ones `<source>_unweighted`,
so with `--weighting both` the two coexist in every output file name.

### Static vs time-resolved

- **Static (default)** — uses the pre-aggregated contact matrices (fast).
- **Time-resolved (`--window-size N`)** — reads the per-frame energy CSVs,
  aggregates into windows of `N` frames, builds one graph per window, runs
  community detection on each, and writes
  `time_resolved/community_evolution_<source>_<method>.csv`
  (`window_start, resindex, label, community` — pivot externally for a
  community-vs-time heatmap).
- **Frame range** — `--frame-start` / `--frame-stop` apply to both paths. When
  set, the static graphs are rebuilt from the per-frame CSVs restricted to the
  range (energy = mean per frame over the window; persistency = unique frames
  per pair / frames in range), so they honour the sub-trajectory.

### PyMOL rendering

One `.pml` per node metric (`node_metrics/<source>_<metric>.pml`, CA spheres
scaled by that metric) and one per community method
(`community_metrics/<source>_<method>.pml`, uniform CA spheres coloured by
community). Each script loads the reference structure as a translucent grey cartoon, places
one pseudoatom per residue (on CA, or the centre of mass for residues without
one, e.g. ions) and draws network edges as **grey** solid sticks scaled by edge
weight. Residue selectors use the PDB residue numbering from the topology, so
multi-chain offsets and non-contiguous numbering are handled correctly.
Solitary nodes (community `-1`, from `--community-min-members`) are not drawn.

### Network matrices

With `--write-matrices` (default on), every static graph is also saved as
adjacency matrices in the same format as `Adjacency/` (`numpy.savetxt` square
matrix + shared `adjacency_labels.dat`), so they can be compared directly with
PyInteraph graph files:

- `<label>_network_adjacency.dat` — the final graph (after threshold and
  `--min-degree`), rebuilt from the graph itself. Binary for unweighted labels.
- `<label>_combined_adjacency.dat` — combined across targets, before thresholding.
- `<label>_<target>_adjacency.dat` — per-target contribution (matrix path only;
  persistency in %, energy as raw kJ/mol-derived weight).

### Output layout

```
<outpath>/Network/
├── network.log
└── <inttype>/
    ├── csv/
    │   ├── nodes_<label>_<method>.csv     (centrality + community per node)
    │   ├── edges_<label>.csv
    │   └── shortest_paths_<label>.csv     (with --shortest-paths)
    ├── node_metrics/<label>_<metric>.pml
    ├── community_metrics/<label>_<method>.pml
    ├── matrices/                          (with --write-matrices)
    ├── network_<label>.graphml            (with --write-graphml)
    └── time_resolved/community_evolution_<label>_<method>.csv   (with --window-size)
```

## Notes
### `na_values`

Every `pandas.read_csv` in the package shares a default NA-string list. `"NA"`
is deliberately **excluded** so sodium atoms/residues named `NA` are not parsed
as missing values.

### Cutoffs

Defaults now match **PyInteraph 1** out of the box: salt bridges `4.5 Å`,
hydrophobic-cluster CoM `5.0 Å`, H-bonds `4.5 Å` donor–acceptor + `120°` D–H–A.
(Earlier releases defaulted to a tighter `3.5 Å` salt bridge, `5.5 Å` HC, and
`150°` H-bond angle — override the flags if you need the old behaviour.) These
are global geometric cutoffs, independent of the interactor scheme.

### Ions and elements

MDAnalysis guesses elements from atom names, which can mis-assign some ions
(e.g. calcium `CAL` → carbon). Hydrophobic-cluster centres of mass are immune —
they use **forcefield** masses loaded from `Interactors/found_topologies/`, not
the guessed ones. Charge-based logic (the salt-bridge gate) likewise uses
forcefield charges. If you add ions as charged groups, make sure the forcefield
assigns them the expected charge (e.g. CHARMM `CAL` = +2), or the gate drops them.

## Experimental: ClinVar variant comparison

`mypsns/services/compare_variants.py` — `CompareVariants(Network)` overlays
ClinVar pathogenicity onto a built network and ranks node/edge metrics by how
well they separate pathogenic from benign missense sites (AUROC / AUPRC,
Mann–Whitney). Structure residues (`Entity == 'P'`) are mapped to the RefSeq
protein by sequence alignment (Biopython), and ClinVar is queried through NCBI
E-utilities (`utils/clinvar_query_utilities.py`, scoring in
`utils/metric_comparison_utilities.py`).

**Not wired into the CLI yet** — use it from Python after `network`:

```python
cv = CompareVariants(topology_pdb=..., general_outpath=..., ...)   # same constructor as Network
cv.load_static(...)
cv.run_analysis(...)
ranked = cv.rank_against_clinvar(gene="RCVRN", reference="NP_002894.1",
                                 chrom="chr17", start=9898189, stop=9905271)
```

Needs internet plus `requests` and `biopython`, which are not declared
dependencies yet.

## Open issues

- **Water parsing** from the forcefield can return inconsistent duplicate rows.
  Workaround: water is matched against a bundled exact-match table
  (`mypsns/resources/sol.csv`) at startup.
- **H-bonds vs PyInteraph** plateau at ~4.5 % (engine difference, *not* a bug —
  see *Interactor roles & PyInteraph reproduction*). SB/HC reach ~0.1 %.
- **`q_interactors`** is less validated than `named_interactors`. For HC
  classification it uses a `min_bonds_to_polar` / `min_bonds_to_charged` BFS
  rule; `min_bonds_to_polar=1` is needed mainly so ALA's Cβ is classified as HC.
  Verify HC labels on a representative system if they look off.
- **Repulsive interactions** are excluded from energy-weighted networks by
  default (attractive-only). Use `--energy-use-abs` to include them.
- **Unmatched residues need a manual fix.** When topology recognition fails,
  the only way forward is to correct the input or forcefield, or by hand: copy
  the chosen `possible_topologies/<res>_uncertain_atoms_<n>.csv` into
  `found_topologies/` as `<res>_ff_atoms.csv`, fill its `ff_atoms` entry in
  `unmatching_residue_topologies.csv`, rename that file to
  `complete_residue_topologies.csv`, then run `interactors`. This flow is not
  automated or fully documented yet.
- **No automated tests.** `tests/` is a placeholder; validation against
  PyInteraph is done with an ad-hoc script (see *Interactor roles & PyInteraph
  reproduction*).

## TODO

The maintained development TODO list lives in [`CLAUDE.md`](CLAUDE.md#todo).
