screenshots pipeline possible_topologies -> found manuale

pymol
- rappresentazione prop nodi come spehers con vdw centrati su CA modificabile (proporzionale a metrica)
- edges con sticks 
NW COMUNITÀ STICKS CON COLORI DIVERSI (solo edges nel network)
1 FILE = 1 METICA

globalprotect: von.univr.it

  - min_bonds_to_polar: minimum number of bonds separating the possibly HC forming atom from a charged (with partial charge > max_abs_charge) atom
                          to classify it as HC

  SET TO 1 mainly for ALA (otherwise no HC)



## myPSNs
Protein Structure Network analysis from molecular-dynamics trajectories.

myPSNs builds residue–residue interaction networks (PSNs) from a GROMACS
forcefield, a topology, and an MD trajectory. It identifies hydrogen bonds,
salt bridges, and hydrophobic clusters per frame; aggregates them into
persistency- and energy-weighted contact matrices; and constructs graphs with
centrality, community detection, and PyMOL renderings.

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
  Calls specific commands foud in contacts/*_cmd.py (through app.command).

- **`config.py`** 
  source of truth for what parameters the pipeline accepts, what their types are, what their defaults are, and what constraints apply.
  Any invalid value is caught HERE, before any pipeline code runs.

  Defines global DEFAULT_NA_VALUES (for pandas.readcsv, omitting NA to preserve sodium atom names) and a Pydantic model for each of the main modules.

  Orchestrates whole application setup (flags > .yml > defaults) and writes complete yaml file
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
    (`<resname>`, then `<resname>_v2`, `_v3`, … for subsequent variants).
    Columns: `rtp_file, resname, atom_name, atom_type, charge, chargegroup, bonds, mass`.
  - `residue_index.csv` — one row per residue instance:
    `Resindex, resnum, Chain, Entity, variant_id`.
  - `atom_index_map.csv` — one row per atom: `Resindex, atom_name, atom_index`.
- `add_inttypes(**kwargs)` — appends interaction-role columns to each variant
  CSV in place. Each variant is analysed in isolation; back-bonds remain in the
  bonds strings but BFS terminates at the variant boundary.
  - **`q_interactors`** (charge-based, PIE-style): H-bond `HAcc`/`HDon`/`HBond`
    from partial charges + connectivity (`max_DA_charge`=−0.3, `min_H_charge`=0.3,
    `DA_masses`=[14,16,19]); salt-bridge `Acidic`/`Basic` from recursive local
    charge (`min_abs_charge`=0.75); `HC` flag + `HC_cluster` label for atoms with
    `|charge| ≤ max_abs_charge` (0.27) and BFS distance `> min_bonds_to_charged`
    (3) from any charged atom (`max_num_heavy_ats`=3).
  - **`named_interactors`** (PyInteraph-compatible, name-based): AND-combined
    `(attribute, [values])` / `not <attribute>` rules per role, with default rule
    sets matching the original PyInteraph atom definitions; `HBond` H atoms are
    filled afterwards by walking from `HDon` atoms via the bonds column.

  Note: at variant level there is no cross-residue context, so the previous
  3-residue BFS widening for HC classification no longer applies. Differences are
  negligible for the standard amino-acid set, but verify HC labels on a
  representative system if they look off.

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
| `--interactors-type`      | `q_interactors`, `named_interactors`  | Which role columns to consult. Repeatable.
| `--target`                | `HBonds`, `Salt_bridges`, `HCs`       | Contact types to compute. Repeatable.
| `--subject-chain`         | detected polymers                     | Chains that are the focus of contacts. Repeatable.
| `--surrounding-chain`     | all chains                            | Candidate partner chains. Repeatable. 
| `--d-a-dist`              | 4.5 Å                                 | H-bond donor–acceptor distance cutoff.
| `--angle`                 | 150°                                  | H-bond donor–H–acceptor angle cutoff.
| `--sb-cutoff`             | 3.5 Å                                 | Salt-bridge atom-pair distance cutoff.
| `--hc-cutoff`             | 5.5 Å                                 | Hydrophobic-cluster centre-of-mass distance cutoff.
| `--matrix-type`           | `persistency`, `energy`               | Which matrices to build after contacts. Repeatable.
| `--energy-use-abs / --no-energy-use-abs` | Energy → strength %% mapping: `\|E\|` vs attractive-only `max(0, -E)`.

Two interactor schemas can be requested independently:

- **`named_interactors`** — PyInteraph-compatible, atom-name-based role rules.
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
   RAM hits the ceiling.
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
| `--interactors-type` | `q`, `named` | Adjacency sub-folders to consume. Repeatable.
| `--target` | `HBonds`, `Salt_bridges`, `HCs` | Contact types contributing to the graph. Repeatable.
| `--edge-weight-source` | `persistency`, `energy` | Build one network per source. Repeatable. |
| `--energy-use-abs / --no-energy-use-abs` | off | Energy → weight: `\|E\|` vs attractive-only `max(0, -E)`. |
| `--node-metric` | `betweenness_centrality` | Centrality driving CA-sphere scaling. Repeatable. |
| `--community` | `louvain` | Community algorithm. `louvain`, `greedy_modularity`, `label_propagation`, `connected_components`. Repeatable. |
| `--min-frequency` | 0.0 | Drop edges with summed weight at or below this. |
| `--min-degree` | 0 | Drop nodes below this degree. |
| `--sphere-scale` | 2.0 | Max CA-sphere vdw radius (Å). |
| `--stick-scale` | 0.5 | Max edge dash radius (Å). |
| `--write-graphml / --no-write-graphml` | off | Also write one GraphML per static graph. |
| `--shortest-paths / --no-shortest-paths` | off | All-pairs weighted shortest paths per graph. |
| `--window-size` | — | If set, also build per-window energy networks and track community evolution. |
| `--window-stride` | = window-size | Step between windows (non-overlapping by default). |
| `--frame-start` / `--frame-stop` | — | Restrict to `[start, stop)`. **When set, static networks are re-aggregated from the per-frame energy CSVs over that range** instead of the pre-aggregated matrices.

### Edge weights

- **persistency** — from `{SB,Hbonds,HC}_adjacency.dat`, summed across the
  requested targets (contact frequency).
- **energy** — from `{SB,Hbonds,HC}_total_energy.dat`. Signed energy is mapped
  to a non-negative weight: attractive-only `max(0, -E)` by default, or `|E|`
  with `--energy-use-abs`.

One graph is built per `(interactors_type × edge_weight_source)`.

### Static vs time-resolved

- **Static (default)** — uses the pre-aggregated contact matrices (fast).
- **Time-resolved (`--window-size N`)** — reads the per-frame energy CSVs,
  aggregates into windows of `N` frames, builds one graph per window, runs
  community detection on each, and writes
  `community_evolution_<inttype>_<source>_<method>.csv`
  (`window_start, resindex, label, community` — pivot externally for a
  community-vs-time heatmap).
- **Frame range** — `--frame-start` / `--frame-stop` apply to both paths. When
  set, the static graphs are rebuilt from the per-frame CSVs restricted to the
  range (energy = mean per frame over the window; persistency = unique frames
  per pair / frames in range), so they honour the sub-trajectory.

### PyMOL rendering

One `.pml` per `(interactors_type × source × node_metric × community)`. Each
script loads the reference structure as a translucent cartoon, draws one CA
sphere per residue scaled by the node metric and coloured by community, and one
dashed edge per interaction scaled by weight (intra-community edges take the
community colour, cross-community edges are grey). Residue selectors use the PDB
residue numbering from the topology, so multi-chain offsets and non-contiguous
numbering are handled correctly.

### Output layout

```
<outpath>/Network/
├── nodes_<inttype>_<source>_<method>.csv     (centrality + community per node)
├── edges_<inttype>_<source>.csv
├── network_<inttype>_<source>_<metric>_<method>.pml
├── network_<inttype>_<source>.graphml         (with --write-graphml)
├── shortest_paths_<inttype>_<source>.csv      (with --shortest-paths)
└── community_evolution_<inttype>_<source>_<method>.csv   (with --window-size)
```

## Notes
### `na_values`

Every `pandas.read_csv` in the package shares a default NA-string list. `"NA"`
is deliberately **excluded** so sodium atoms/residues named `NA` are not parsed
as missing values.

### Cutoffs

Defaults are tuned for the systems myPSNs was developed against (CHARMM36,
protein + lipids + ions). The salt-bridge default (3.5 Å) is tighter than the
common 4.0–4.5 Å heavy-atom criterion; raise `--sb-cutoff` to 4.5 to match
PyInteraph-style definitions. Hydrophobic clusters use a 5.5 Å CoM cutoff.

### Ions and elements

MDAnalysis guesses elements from atom names, which mis-assigns some ions
(e.g. calcium "CAL" → carbon). The corrected-topology step handles chain/segid
alignment; if you add ions, double-check their element assignment .

## Open issues

- **Water parsing** from the forcefield can return inconsistent duplicate rows.
  Workaround: water is matched against a bundled exact-match table
  (`mypsns/resources/sol.csv`) at startup.
- **`q_interactors`** is less validated than `named_interactors`. For HC
  classification it uses a `min_bonds_to_polar` / `min_bonds_to_charged` BFS
  rule; `min_bonds_to_polar=1` is needed mainly so ALA's Cβ is classified as HC.
  Verify HC labels on a representative system if they look off.
- **Repulsive interactions** are excluded from energy-weighted networks by
  default (attractive-only). Use `--energy-use-abs` to include them.

