# myPSNs

Protein Structure Network (PSN) analysis of molecular-dynamics trajectories.
From a GROMACS forcefield, a structure and a trajectory, myPSNs finds hydrogen
bonds, salt bridges and hydrophobic contacts in every frame, turns them into
residue × residue persistency and interaction-energy matrices, and builds
residue networks with centralities, communities and PyMOL scripts.

myPSNs re-implements **PyInteraph 1** and reads charges, masses, Lennard-Jones
parameters and terminal patches from the forcefield itself. With its default
settings it reproduces PyInteraph 1 on PyInteraph's own test case, which ships
in `tests/data/`: salt bridges differ by 0.07 %, hydrophobic contacts by
0.11 %, H-bonds by 3.9 %, and every PyInteraph edge is found
([Validation](#validation-against-pyinteraph-1)).

```
ff-reader → topology → interactors → contacts      mypsns run
                                         │
                                      network       mypsns network
                                         │
                        compare  ·  compare-variants (experimental)
```

**Part I — Using myPSNs:** [Install](#install) · [Quick start](#quick-start) ·
[Your own system](#your-own-system) · [Command reference](#command-reference)

**Part II — Developer guide:** [Repository layout](#repository-layout) ·
[Architecture](#architecture) · [Data flow and run folder](#data-flow-and-run-folder) ·
[Extending myPSNs](#extending-mypsns) · [Interactor rules](#interactor-rules) ·
[Testing](#testing) · [Validation](#validation-against-pyinteraph-1) ·
[Conventions](#conventions) · [Known issues](#known-issues-and-limitations)

The development TODO list and notes for coding agents are in
[`CLAUDE.md`](CLAUDE.md); [`CHANGELOG.md`](CHANGELOG.md) lists changes that
affect existing run folders.

---

# Part I — Using myPSNs

## Install

From the repository root:

```bash
conda env create -f environment.yml   # Python 3.11, all dependencies, `pip install -e .`
conda activate mypsns
mypsns --version
```

Without conda, `pip install -e .` installs the package and the `mypsns`
command. Two extras exist: `pip install -e '.[dev]'` (pytest, black, ruff) and
`pip install -e '.[variants]'` (`requests`, `biopython`, needed only by the
experimental `compare-variants`). `environment.yml` already includes both.

## Quick start

The PyInteraph 1 test case (one 149-residue chain, 151 frames, CHARMM36) ships
with the repository. From the repository root:

```bash
D=tests/data/pyinteraph_testcase
OUT=/tmp/mypsns_quickstart

# 1. contacts: ff-reader → topology → interactors → contacts   (~3.5 min, most of it parsing CHARMM36)
mypsns run -t $D/structures/sim.prot.gro -x $D/structures/traj.xtc \
           -ff $D/forcefields/charmm36-apr2019_vm.ff -o $OUT --nc 4

# 2. networks, keeping residue pairs in contact in more than 20 % of the frames   (seconds)
mypsns network -t $D/structures/sim.prot.gro -o $OUT --min-persistency 20

# 3. agreement with PyInteraph 1 on the same trajectory
mypsns compare --preset pyinteraph_testcase --computed-dir $OUT/Adjacency/named_interactors
mypsns compare --truth macroIIN=$D/pyinteraph/macro-IIN.dat \
  --computed macroIIN=$OUT/Network/named_interactors/matrices/persistency_unweighted_network_adjacency.dat

# 4. what each stage produced
mypsns report $OUT
```

The first `compare` checks the three contact matrices against the preset's
tolerances (exit code 1 above them):

```text
┃ Type   ┃ Edges T ┃ Edges C ┃ Shared ┃ Only T ┃ Only C ┃ ≠ value ┃ % diff ┃ Status     ┃
│ HBonds │     282 │     309 │    282 │      0 │     27 │     213 │  3.929 │ ok ≤ 4.2 % │
│ HC     │      79 │      79 │     79 │      0 │      0 │      56 │ 0.1118 │ ok ≤ 0.2 % │
│ SB     │      21 │      21 │     21 │      0 │      0 │      17 │ 0.0735 │ ok ≤ 0.2 % │
```

T is PyInteraph (truth), C is myPSNs (computed); `% diff` is
`100·Σ|C − T| / ΣT`. The second compares the binary network with PyInteraph's
macro-IIN at the same 20 % cutoff: 192 PyInteraph edges, all found, plus 5
H-bond pairs that myPSNs puts just above the cutoff.

### What a run folder contains

```
$OUT/
├── run_config.yaml                    every setting each stage used
├── forcefield/                        the parsed forcefield (ff-reader)
├── Topology_recognition/              residue ↔ forcefield matching (topology)
├── Interactors/                       atom roles per residue variant (interactors)
├── Adjacency/                         contacts
│   ├── sim.prot_correct.pdb           the structure renumbered by residue index
│   └── named_interactors/             (and q_interactors/, one folder per role scheme)
│       ├── adjacency_labels.dat       row/column order of every matrix: index, Chain:RES:index
│       ├── {SB,Hbonds,HC}_adjacency.dat          persistency, % of frames
│       ├── {SB,Hbonds,HC}_total_energy.dat       mean interaction energy, kJ/mol
│       ├── {SB,Hbonds,HC}_final.csv              residue pairs, frames in contact, mean distance
│       └── graph_report/plots/*_elbow.png        edge-threshold scans, to choose --min-persistency
└── Network/
    └── named_interactors/
        ├── csv/nodes_<label>_<method>.csv        per residue: centralities, community
        ├── csv/edges_<label>.csv
        ├── node_metrics/<label>_<metric>.pml     PyMOL: spheres sized by a centrality
        ├── community_metrics/<label>_<method>.pml   PyMOL: spheres and edges coloured by community
        └── matrices/<label>_network_adjacency.dat   the final graph, same format as Adjacency/
```

`<label>` is `persistency`, `energy`, `persistency_unweighted` or
`energy_unweighted` (one graph per edge-weight source and weighting). Open a
`.pml` script in PyMOL (`pymol <file>.pml`). Full layouts are in the
[command reference](#command-reference).

The same two steps from a config file (it points at the test case and writes
`/tmp/mypsns_example`):

```bash
mypsns run     --config example_config.yaml
mypsns network --config example_config.yaml
```

### When topology recognition stops

`topology` only accepts a forcefield entry whose atom names match the
structure's exactly, and stops the pipeline otherwise. To see it, rename every
glycine `HA1` to `HA3` in a copy of the test structure:

```bash
sed -E 's/^(.{5}GLY  )  HA1/\1  HA3/' $D/structures/sim.prot.gro > /tmp/gly_ha3.gro
OUT2=/tmp/mypsns_halt
mypsns run -t /tmp/gly_ha3.gro -x $D/structures/traj.xtc \
           -ff $D/forcefields/charmm36-apr2019_vm.ff -o $OUT2 --nc 4      # stops, exit code 1
```

```text
Topology recognition could not associate 1 residue topology/ies to the forcefield:
  - GLY (topology GLY): no exact match; the best of 14 candidate(s) shares 6/7 atom names
...
  1. promote a candidate. List the candidates and how their atoms pair with the structure's:
       mypsns accept-topology -o /tmp/mypsns_halt
     then accept one per topology (by index, or a hand-edited CSV path):
       mypsns accept-topology -o /tmp/mypsns_halt GLY=0
  2. add the missing residue(s) to the forcefield, or correct the structure's atom names,
     and re-run `mypsns topology --overwrite`.
```

List the candidates, check one, accept it, and continue:

```bash
mypsns accept-topology -o $OUT2                    # candidates and their renames
mypsns accept-topology -o $OUT2 GLY=0 --dry-run    # what would change; writes nothing
mypsns accept-topology -o $OUT2 GLY=0
mypsns run -t /tmp/gly_ha3.gro -x $D/structures/traj.xtc \
           -ff $D/forcefields/charmm36-apr2019_vm.ff -o $OUT2 --nc 4      # skips the finished stages
```

```text
              GLY  (GLY, 6 residue(s)) — unresolved
┃ #  ┃ forcefield residue(s) ┃ renames (forcefield → structure) ┃
│ 0  │ GLY                   │ HA1→HA3                          │
│ 1  │ DGLY                  │ HA1→HA3                          │
│ 2  │ GLY                   │ identical to #0                  │
...
```

The contacts of this run equal those of the unmodified one. Choosing a
candidate, hand-editing one and what acceptance writes are described under
[`accept-topology`](#mypsns-accept-topology).

### Ranking network metrics against ClinVar (experimental)

For a human protein, `compare-variants` asks ClinVar for the gene's missense
variants and ranks every network metric by how well it separates pathogenic
sites from the rest (AUROC, AUPRC, Mann–Whitney). It needs internet access and
the `variants` extra. The test case is a yeast protein, so this example uses a
recoverin–GRK1 run:

```bash
mypsns compare-variants --config <outpath>/run_config.yaml --gene RCVRN --reference NP_002894.1
```

See [`compare-variants`](#mypsns-compare-variants-experimental).

## Your own system

**Inputs.**

- **Forcefield:** a GROMACS `.ff` folder. ff-reader reads every `*.rtp`, the
  terminal databases `*.n.tdb` / `*.c.tdb` (named like the `.rtp`, e.g.
  `merged.rtp` ↔ `merged.n.tdb`), and the nonbonded `.itp`
  (`ffnonbonded.itp`). `mypsns guess_ff --parent <dir>` lists the forcefield
  folders under a directory.
- **Structure:** PDB or GRO with the forcefield's residue and atom names. Every
  residue needs a forcefield entry, terminal patches included; non-standard
  residues (a lipidated Gly, an ion) need theirs in the forcefield.
- **Trajectory:** any format MDAnalysis reads, with the structure's atoms in the
  same order.

**Configuration.** Every setting has a default (`mypsns <command> --help`) and
can be set in a YAML file (`--config`, see `example_config.yaml`) or as a flag.
Flags override the file. Each stage records the settings it used in its section
of `<outpath>/run_config.yaml`, leaving the other sections as they are, so the
file always describes the whole run and can be passed back with `--config`.
Paths are stored as given, so a relative path only resolves from the same
working directory. `mypsns validate --config <file> --deep` checks a config,
its files and that MDAnalysis can open them.

**Re-running.** A stage skips when its output exists: `mypsns run` resumes after
the last finished stage. `--overwrite` recomputes. Interaction roles are written
by `interactors`, so after changing the interactor rules or options run
`mypsns interactors --overwrite`, then `mypsns contacts --overwrite`.

**Choosing a cutoff.** `Adjacency/<scheme>/graph_report/` scans the persistency
threshold for every combination of contact types (largest connected component,
degree distribution); its elbow plots help choose `network --min-persistency`.

**Non-standard residues, ions, lipids.** Roles come from name-based rules: add a
rule for the residue (see [Interactor rules](#interactor-rules), with a
myristoylated glycine and a Ca²⁺ example).

**Run time.** ff-reader takes ~3 min for CHARMM36. Contacts scale with frames ×
atoms, and H-bonds dominate: 151 frames of the test case take ~20 s at
`--nc 4`, 20 001 frames of a 219-residue complex ~55 min. Raise `--nc` on a
larger machine; `--beg-frm`/`--end-frm` restrict the frames.

## Command reference

`mypsns --help` lists the commands; `mypsns <command> --help` lists every flag.

| Command | Does |
|---|---|
| `run` | ff-reader → topology → interactors → contacts |
| `ff-reader` | parse a GROMACS forcefield folder |
| `topology` | match every residue of the structure to a forcefield entry |
| `accept-topology` | promote a candidate for residues `topology` could not match |
| `interactors` | give every atom its interaction roles |
| `contacts` | per-frame contacts over the trajectory → matrices |
| `network` | graphs, centralities, communities, PyMOL scripts |
| `compare` | compare residue matrices with reference graphs (PyInteraph) |
| `compare-variants` | rank network metrics against ClinVar (experimental) |
| `validate` | check a config and the files it names |
| `report` | summarise what a run folder contains |
| `guess_ff` | find forcefield folders under a directory |

**Common behaviour.** Global flags go before the command: `-vv` (debug log),
`-q` (warnings only), `--log-file <file>` (debug log to a file), `--version`.
The stage commands (`run`, `ff-reader`, `topology`, `interactors`, `contacts`,
`network`) take `--config`, `--dry-run` (validate and print the plan),
`--overwrite` and, except `ff-reader`, `--resume`; `-t/--topology`,
`-x/--trajectory`, `-ff/--forcefield` and `-o/--outpath` set the paths where
the stage needs them. Exit
codes: 0 success, 1 failure (including a topology halt and a failed `compare`),
2 invalid configuration or missing input, 3 missing optional dependency.

### `mypsns run`

Runs the four contact stages, skipping those already done. Only the flags below
exist; every other setting comes from `--config`.

| Flag | Meaning |
|---|---|
| `--config`, `-t`, `-x`, `-ff`, `-o` | config file and paths |
| `--nc` | worker processes for contacts |
| `--target` | contact types: `HBonds`, `Salt_bridges`, `HCs` (repeatable) |
| `--only` / `--skip` | run only / skip these stages: `ff-reader`, `topology`, `interactors`, `contacts` (repeatable) |
| `--resume` / `--overwrite` / `--dry-run` | skip finished stages (the default) / recompute everything / plan only |

### `mypsns ff-reader`

Parses the forcefield into `<outpath>/forcefield/`: `final_ff_residues.csv`
(one row per atom of every residue: `rtp_file, resname, atom_name, atom_type,
charge, chargegroup, bonds, sigma, epsilon, mass`), `final_ff_ters.csv` (the
terminal patches: `tdb_file` instead of `rtp_file`, and the patch instructions
in `bonds`) and `atomtypes_se.csv` (mass, σ, ε per atom type, for energies).

| Flag | Default | Meaning |
|---|---|---|
| `--nonbonded-name` | `ffnonbonded.itp` | `.itp` with atom types, masses, σ/ε |
| `--no-additional-residues` | off | do not add the residues bundled in `mypsns/resources/` |

### `mypsns topology`

Groups the structure's residues by atom-name set and matches each group to a
forcefield entry, applying terminal patches where needed. `-ff` takes the
original `.ff` folder (parsed first) or an already parsed `forcefield/` folder.
Writes `<outpath>/Topology_recognition/`. If any residue has no exact match it
stops with exit code 1 and prints the [`accept-topology`](#mypsns-accept-topology)
commands that resolve it.

| Flag | Default | Meaning |
|---|---|---|
| `--max-to-check` | 30 | forcefield entries tried per residue, most similar names first; `-1` = all. Raise it for large custom residues. |

### `mypsns accept-topology`

Resolves residues that `topology` could not match. Without choices it lists,
per unmatched residue topology, the candidates recognition wrote to
`Topology_recognition/possible_topologies/` and how their atoms pair with the
structure's; with `TOPOLOGY=CHOICE` arguments it promotes them.

| Argument / flag | Default | Meaning |
|---|---|---|
| `TOPOLOGY=CHOICE` (repeatable) | — | `CHOICE` is a candidate index (`possible_topologies/<TOPOLOGY>_uncertain_atoms_<CHOICE>.csv`) or the path of a CSV |
| `-o` / `--outpath` | — | run folder |
| `--config` | — | only `outpath` and `topology_recognition.na_values` are used; `-o` wins |
| `--dry-run` | off | check the choices and show the renames; write nothing |

**Choosing.** `TOPOLOGY` is the name in the halt message: the residue name, plus
`_<n>` for its second, third… atom set (e.g. a terminal copy). In the candidate
table, *forcefield residue(s)* are the entries the candidate is built from, the
residue and the terminal patch applied to it: choose on that (an N-terminal
L-Met wants `NH3+, MET`, not `GLY-NH3+` or `DMET`). *renames* pairs the atoms
whose names differ, in file order after identical names pair up. A red `!`
marks a pair of different masses: that pairing is wrong and needs a hand-edited
copy. `identical to #n` marks a duplicate file; skip it. A topology with no
candidate (no forcefield entry with its atom count) needs a CSV you write, or a
larger `--max-to-check` and `topology --overwrite`.

**What acceptance writes.** For each choice it copies the candidate to
`found_topologies/<topology>_ff_atoms.csv` and fills its `structure_atom_name`
column with the structure's names (`atom_name` and `bonds` keep the
forcefield's; `interactors` renames them), records the source in
`unmatching_residue_topologies.csv`, and once every row is resolved writes
`complete_residue_topologies.csv`, which `interactors` reads. All choices are
checked before anything is written; each run is logged to
`Topology_recognition/accept_topology.log`.

**Hand-edited candidate.** Copy the closest candidate, edit it, pass its path:

```bash
cp $OUT2/Topology_recognition/possible_topologies/GLY_uncertain_atoms_0.csv gly_fixed.csv
# edit atom_name (the structure's names are simplest) and every bonds entry naming a renamed atom
mypsns accept-topology -o $OUT2 GLY=gly_fixed.csv
```

The file keeps all ten columns (`rtp_file, resname, atom_name, atom_type,
charge, chargegroup, bonds, sigma, epsilon, mass`). A `bonds` entry naming an
atom the residue lacks raises a warning (`-X`/`+X` name the previous/next
residue and are not checked).

**Things to know.** `topology --overwrite` discards acceptances (accept again).
From `interactors` on, the residue's atoms carry the structure's names and the
name-based rules see those: a name a rule treats specially (the HC rule excludes
the AMBER amide `H` but keeps the CHARMM `HN`) changes the residue's roles.
After changing an earlier choice, run `interactors --overwrite` and
`contacts --overwrite`.

### `mypsns interactors`

Rebuilds the whole system from the matched residues (chains, polymer bonds,
terminal charge groups), groups residues into unique chemical variants and
gives every atom its roles under two schemes: `named_interactors`
(PyInteraph-compatible, name-based, the default) and `q_interactors`
(charge/connectivity-based). Both role columns are always written to
`<outpath>/Interactors/`. See [Interactor rules](#interactor-rules).

| Flag | Default | Meaning |
|---|---|---|
| `-m` / `--mode` | `named_interactors` | scheme(s) recorded in the config |
| `--min-residue-charge` | 0.8 e | named: net charge a `Basic`/`Acidic` group needs; 0 = off |
| `--max-da-charge` | −0.3 e | q: max (least negative) charge of an H-bond donor/acceptor |
| `--min-h-charge` | 0.3 e | q: min charge of an H-bond hydrogen |
| `--da-masses` | 14 16 19 | q: rounded masses of donor/acceptor elements (repeatable; replaces the default) |
| `--min-abs-charge` | 0.75 e | q: min \|local charge\| of an `Acidic`/`Basic` group |
| `--max-abs-charge` | 0.27 e | q: max \|charge\| of an HC atom |
| `--max-num-heavy-ats` / `--min-num-heavy-ats` | 3 / 2 | q: heavy atoms per HC cluster (the minimum is soft) |
| `--min-bonds-to-charged` / `--min-bonds-to-polar` | 3 / 1 | q: an HC atom is more than this many bonds from a charged / polar heavy atom |

Custom `named_interactors_rules` are set in the YAML only.

### `mypsns contacts`

Finds the contacts in every frame for each interactor scheme and builds the
matrices. Reads `<outpath>/Interactors/`, writes `<outpath>/Adjacency/`.

| Flag | Default | Meaning |
|---|---|---|
| `--nc` | 2 | worker processes |
| `--ram-thresh` | 80 | % of RAM at which H-bond batches stop growing |
| `--nframes` | 500 | initial H-bond batch, frames |
| `--beg-frm` / `--end-frm` | 0 / −1 | first / last frame, both inclusive (−1 = last) |
| `--interactors-type` | `q_interactors`, `named_interactors` | role schemes, one output folder each |
| `--target` | `HBonds`, `Salt_bridges`, `HCs` | contact types |
| `--subject-chain` / `--surrounding-chain` | polymer chains / all chains | contacts are computed within each subject chain and between it and every surrounding chain |
| `--d-a-dist` | 4.5 Å | H-bond donor–acceptor distance |
| `--angle` | 120° | H-bond D–H–A angle, kept when ≥ |
| `--d-h-cutoff` | 1.4 Å | donor–hydrogen pairing distance (MDAnalysis' 1.2 Å misses Cys S–H, 1.32 Å) |
| `--sb-cutoff` | 4.5 Å | salt bridge: any atom pair of the two charged groups |
| `--hc-cutoff` | 5.0 Å | hydrophobic: side-chain centres of mass |
| `--matrix-type` | `persistency`, `energy` | matrices to build |
| `--energy-use-abs` | off | threshold scans of energies use \|E\| instead of attractive-only max(0, −E) |
| `--thr-step` | 0.1 | step (%) of the threshold scans |

The cutoffs are PyInteraph 1's and apply to both schemes.

```
<outpath>/Adjacency/
├── Contacts.log
├── <name>_correct.pdb                      the structure with residue index as resid
└── <scheme>/                               named_interactors/, q_interactors/
    ├── Salt_bridges/  H_bonds/  HCs/       per-frame contacts of each chain pair:
    │                                       *_upto_<frame>.csv, *_merged.csv, *_energy.csv
    ├── {SB,Hbonds,HC}_final.csv            res_1, res_2, n_frames_contact, mean_dist, chains, labels
    ├── {SB,Hbonds,HC}_final_visited_frames.csv   the frames of every pair
    ├── adjacency_labels.dat                matrix row order
    ├── {SB,Hbonds,HC}_adjacency.dat        persistency, % of frames
    ├── combined_adjacency.dat              % of frames with any contact type
    ├── {SB,Hbonds,HC}_{LJ,Coulomb,total}_energy.dat   mean energy per frame, kJ/mol
    ├── all_total_energy.dat                sum over types
    ├── graph_report/                       persistency threshold scans (csvs/, plots/)
    └── energy_graph_report/                the same for energies
```

### `mypsns network`

Builds graphs from `Adjacency/` and writes `<outpath>/Network/`. Run it after
`contacts`; `mypsns run` does not include it.

| Flag | Default | Meaning |
|---|---|---|
| `--interactors-type` | `q_interactors`, `named_interactors` | schemes to read |
| `--target` | `HBonds`, `Salt_bridges`, `HCs` | contact types combined into each graph |
| `--edge-weight-source` | `persistency`, `energy` | one graph per source |
| `--weighting` | `both` | `weighted` (weight = value), `unweighted` (weight 1), or `both` |
| `--min-persistency` / `--min-energy` | 0 / 0 | keep edges whose value is > this (%, energy as % of the strongest edge) |
| `--min-degree` | 0 | drop nodes below this degree |
| `--energy-use-abs` / `--no-energy-use-abs` | off | energy weight \|E\| instead of attractive-only max(0, −E) |
| `--node-metric` | `betweenness_centrality` | also `degree`, `weighted_degree`, `degree_centrality`, `closeness_centrality`; one `.pml` each (repeatable) |
| `--community` | `louvain` | also `greedy_modularity`, `label_propagation`, `connected_components` (repeatable) |
| `--community-seed` | 42 | seed of the stochastic methods |
| `--community-min-members` | 2 | smaller communities dissolve (members not drawn) |
| `--sphere-scale` / `--stick-scale` | 2.0 / 0.5 Å | PyMOL: largest sphere / edge radius |
| `--reference-pdb` | contacts' corrected PDB | structure the `.pml` scripts load |
| `--write-matrices` / `--write-graphml` / `--shortest-paths` | on / off / off | extra outputs |
| `--window-size` / `--window-stride` | — / = size | also build one energy network per frame window and track communities |
| `--frame-start` / `--frame-stop` | — | use only frames in [start, stop): graphs are rebuilt from the per-frame CSVs |
| `-n` / `--n-frames` | — | informational, logged only |

**Edge weights.** Persistency comes from `*_adjacency.dat`; energy from
`*_total_energy.dat`, mapped to max(0, −E) (or \|E\|) and normalised to % of the
strongest edge. Several targets combine into one matrix per source:

| | combined value | weighted graph | unweighted graph |
|---|---|---|---|
| persistency | element-wise **max** over targets (≤ 100 %) | value, if `> --min-persistency` | 1, if `> --min-persistency` |
| energy | **sum** of the per-target weights, as % of the strongest | value, if `> --min-energy` | 1, if `> --min-energy` |

Both weightings keep the same edges. The unweighted persistency graph is the
analogue of PyInteraph's macro-IIN (`filter_graph -d hb -d hc -d sb -t <cutoff>`,
which also keeps values `>` the cutoff); fed PyInteraph's own graphs it
reproduces a PyInteraph macro-IIN exactly.

**PyMOL scripts.** Each `.pml` loads the reference structure as a grey cartoon,
puts a sphere on every network residue (CA, or the centre of mass when there is
none, e.g. an ion) and draws edges as sticks scaled by weight. In
`node_metrics/` spheres scale with the metric and edges are grey; in
`community_metrics/` spheres and intra-community edges take the community's
colour (16 colours, reused beyond 16) and edges between communities stay grey.

**Time-resolved.** With `--window-size N` the per-frame energies are aggregated
into windows of N frames, one graph is built per window and
`time_resolved/community_evolution_<label>_<method>.csv`
(`window_start, resindex, label, community`) tracks the communities.

```
<outpath>/Network/
├── network.log
└── <scheme>/
    ├── csv/nodes_<label>_<method>.csv      centralities + community per residue
    ├── csv/edges_<label>.csv
    ├── csv/shortest_paths_<label>.csv      with --shortest-paths
    ├── node_metrics/<label>_<metric>.pml
    ├── community_metrics/<label>_<method>.pml
    ├── matrices/                           with --write-matrices: adjacency_labels.dat,
    │                                       <label>_network_adjacency.dat (final graph),
    │                                       <label>_combined_adjacency.dat (before the threshold),
    │                                       <label>_<target>_adjacency.dat (per target)
    ├── network_<label>.graphml             with --write-graphml
    └── time_resolved/                      with --window-size
```

### `mypsns compare`

Compares residue × residue matrices type by type on the strict lower triangle
(the diagonal is ignored: PyInteraph's `hb-graph.dat` keeps intra-residue
H-bonds there). An edge is a value > 0; values closer than `--atol` count as
equal; `% diff` is `100·Σ|C − T| / ΣT`.

| Flag | Meaning |
|---|---|
| `--preset NAME` | a reference shipped in `tests/data/` (below): supplies the truth graphs and per-type tolerances |
| `--truth-dir DIR` + `--truth-suffix S` | reference graphs named like PyInteraph: `{hb,hc,sb}-graph<S>.dat` |
| `--computed-dir DIR` | a myPSNs `Adjacency/<scheme>/` folder |
| `--truth TYPE=PATH` / `--computed TYPE=PATH` | any pair of matrices, matched by `TYPE` (repeatable) |
| `--type` | compare only these types |
| `--max-pct-diff TYPE=PCT` | exit 1 if a type exceeds PCT (overrides the preset's) |
| `--atol` | default 0.01 |
| `-o DIR` / `--labels FILE` | write `possible_sele_mismatch.csv` (edge in one matrix only), `possible_params_mismatch.csv` (both, values differ) and `comparison_summary.csv`; add residue labels from an `adjacency_labels.dat` |

| Preset | Reference | Inputs shipped | Tolerances H-bonds / HC / SB |
|---|---|---|---|
| `pyinteraph_testcase` | PyInteraph 1 example, 149 residues, 151 frames | yes | 4.2 / 0.2 / 0.2 % |
| `rec_grk1` | recoverin–GRK1, PyInteraph 1 run of 2018 on 1 µs (5 × 200 ns) | no (1.3 GB xtc) | 5.0 / 0.25 / 0.2 % |
| `rec_grk1_r1` | the same run, replica r1 (frames 0–20000) | no | 5.0 / 0.2 / 0.2 % |

Each preset has a `config.yaml` with the settings that reproduce it
([Validation](#validation-against-pyinteraph-1)).

### `mypsns compare-variants` (experimental)

Ranks network metrics by how well they separate pathogenic missense sites of a
human gene from the rest of the protein. Needs internet access and the
`variants` extra (`requests`, `biopython`); without them it exits 3 with an
install hint.

1. **Rebuilds the graphs** from `Adjacency/` with the `network:` settings, by
   the same code as `network`. With `--config` pointing at the
   `run_config.yaml` that `network` wrote, it ranks exactly those graphs.
   Centralities are computed in memory; `Network/<scheme>/csv/` is not touched.
2. **Queries ClinVar** (`<gene>[gene]`, optionally narrowed to a GRCh38
   `--chrom/--start/--stop` window) through NCBI E-utilities and keeps missense
   records. A protein position is `non_benign` with at least one pathogenic or
   likely-pathogenic record, `benign` with only benign ones, `vus` otherwise.
3. **Maps positions** by aligning the RefSeq protein (global, BLOSUM62) to every
   protein chain; the best chain above `--min-identity` is the gene's.
4. **Scores each graph:** positives are the `non_benign` sites, negatives every
   other aligned residue (`--negative rest`) or only `benign` sites
   (`--negative benign`). Residues that do not align (partner chains, ions) are
   not scored.

| Flag | Default | Meaning |
|---|---|---|
| `--gene` | required | gene symbol queried in ClinVar |
| `--reference` | required | RefSeq protein the ClinVar `p.` positions refer to (the NP_ of the NM_ in ClinVar titles, e.g. `NM_002903.3` → `NP_002894.1`), or a raw sequence |
| `--chrom` / `--start` / `--stop` | — | optional GRCh38 window (all three or none) |
| `--email` / `--api-key` | env `NCBI_EMAIL` / `NCBI_API_KEY` | NCBI contact and key (10 requests/s instead of 3); the key is never written to `run_config.yaml` |
| `--negative` | `rest` | `rest` or `benign` |
| `--min-identity` | 70 | % identity for a chain to count as the gene's |
| `--interactors-type`, `--target`, `--edge-weight-source`, `--energy-use-abs`, `--weighting`, `--min-persistency`, `--min-energy`, `--min-degree`, `--frame-start`, `--frame-stop` | from `network:` | graph rebuild, as in `network` |
| `--overwrite` / `--dry-run` | off | rewrite this gene's tables / plan only |

```
<outpath>/Network/
├── compare_variants.log
└── <scheme>/variants/
    ├── <gene>_ranking_<label>.csv     per metric: auroc, auroc_oriented, direction, auprc, mwu_p, n_pos, n_neg
    ├── <gene>_nodes_<label>.csv       per aligned residue: metrics, ref_pos, ref_aa, ClinVar label
    ├── <gene>_clinvar_sites.csv       per missense site: counts, label, chain, resindex, resnum, structure vs reference aa
    └── <gene>_clinvar_variants.csv    the raw ClinVar records (a dated snapshot: keep it with the results)
```

Notes: an empty ranking (header only, with a warning) means no pathogenic site
falls on a modelled residue. RCVRN, for example, has 18 missense sites, none
pathogenic. The GRCh38 window can drop records the gene query keeps (for RCVRN,
one missense VUS with no position in the window); omit it when every record
matters. `classify_significance` sends everything not pathogenic or benign to
VUS, including conflicting classifications, `not provided` and empty values.
Mismatches between `structure_aa` and `ref_aa` are expected for orthologues;
many of them suggest the wrong reference isoform. From Python it is a `Network`
subclass: `CompareVariants(...).load_static(...)`, then
`rank_against_clinvar(gene=..., reference=...)` and `write_rankings(...)`.

### `mypsns validate`, `report`, `guess_ff`

- `mypsns validate --config FILE [--deep]`: schema, paths, forcefield content,
  writable outpath; `--deep` also opens the structure and trajectory.
- `mypsns report OUTPATH [--json]`: which stages finished, per scheme.
- `mypsns guess_ff --parent DIR [--recursive-search] [--maximum-depth N]`:
  folders under DIR that look like a GROMACS forcefield.

---

# Part II — Developer guide

## Repository layout

```
mypsns/
├── cli.py                    Typer app: global options, command registration, main()
├── config.py                 Pydantic models per stage, defaults, YAML load/merge/dump
├── logging_utils.py          Rich console, banners, summary tables
├── commands/<name>_cmd.py    one module per command; _helpers.py shared by all
├── services/                 stage logic, one class per stage (FFReader … Network, CompareVariants)
├── utils/<stage>_utilities.py   helper functions (most of the code)
└── resources/sol.csv         exact-match water table
tests/                        pytest suite; tests/data/ = reference datasets (presets)
example_config.yaml           every setting at its default, runnable on the test case
CLAUDE.md                     notes for coding agents and the TODO list
CHANGELOG.md                  changes that affect existing run folders
```

## Architecture

A command is a thin layer: its `*_cmd.py` collects the flags into a nested
overrides dict, `build_pipeline_config` merges them over the YAML and the
defaults (**flags > YAML > defaults**) and validates the result, the command
prints a summary (or stops on `--dry-run`), checks whether the stage is already
done, calls the service and records its section of `run_config.yaml`
(`dump_run_config(cfg, outpath, sections=[...])`).

| Layer | Owns |
|---|---|
| `cli.py` | registration (`app.command(...)`), global flags, logging setup in the callback |
| `config.py` | the parameters: one Pydantic model per stage (`FFReaderConfig`, `TopologyConfig`, `InteractorsConfig`, `ContactsConfig`, `NetworkConfig`, `CompareVariantsConfig`), composed in `PipelineConfig`. Every model is `extra="forbid"`, so a misspelt YAML key fails validation. The single copy of the default rules, `NAMED_INTERACTOR_DEFAULTS`. |
| `commands/` | the only place that knows both the flags and the service signature |
| `services/` | stage logic (`Topology_recognition`, `Interactors`, `Contacts`, `Network`, …) |
| `utils/` | the per-stage functions the services call |

Errors surface as one Rich line and an exit code (see the command reference);
tracebacks appear with `-vv`.

## Data flow and run folder

Every stage reads the previous stage's folder in `<outpath>` and writes its own.

**ff-reader** (`services/ff_reader.py`, `utils/ff_reader_utilities.py`) parses
the `.rtp` residues (atom type, charge, charge group, bonds) and the `.tdb`
terminal patches, joins σ, ε and mass from the nonbonded `.itp`, and writes
`forcefield/final_ff_residues.csv`, `final_ff_ters.csv`, `atomtypes_se.csv`. In
`final_ff_ters.csv` the `bonds` column also carries the patch instructions
(`rX` replace, `-X` delete, `+X` an `[ add ]` control atom outside the patch);
`is_patch_instruction` tells them from bonds. The masses here are the ones used
later for hydrophobic centres of mass.

**topology** (`services/topology_recognition.py`) groups residues by atom-name
set (`resid_topologies/<name>_rixs.txt`, `_atoms.csv`), then tries forcefield
entries in order of name similarity (`--max-to-check`): an exact atom-name match
first, then each terminal patch applied (`get_ter_ffres` replaces, renames,
deletes and adds atoms in every bonds list). Matches go to
`found_topologies/<topology>_ff_atoms.csv`, which keeps the forcefield names in
`atom_name`/`bonds` and the structure's in `structure_atom_name`; the best
non-matches go to `possible_topologies/`. It writes
`complete_residue_topologies.csv`, or `unmatching_residue_topologies.csv` and
raises `TopologyAssociationError`. A patch that leaves one of its atoms with no
bond (a hand-made patch, usually) is reported as a warning. Topology keys come
from `topology_key()`.

**interactors** (`services/interactors.py`) renames atoms to the structure's
names (`build_wholetop`), detects chains (polymers by `+` bonds, `Entity` P;
other molecules by `rtp_file`, M), adds polymer bonds and terminal charge
groups, and deduplicates residues into variants. It writes
`Interactors/found_topologies/<variant>.csv` (`<resname>`, `<resname>_2`, …),
`residue_index.csv` (`Resindex, resnum, Chain, Entity, variant_id,
topology_resname`, the latter being the structure's name) and
`atom_index_map.csv`. `add_inttypes` then appends the `q_interactors`,
`HC_cluster` and `named_interactors` columns to every variant file.

**contacts** (`services/contacts.py`, `utils/contacts_utilities.py`) rebuilds
the per-atom table from those three files, writes the renumbered
`<name>_correct.pdb`, builds compact MDAnalysis selections per chain and role,
and dispatches per chain pair and contact type to `--nc` workers: salt bridges
and HC by PBC-aware capped distances (HC centres of mass with forcefield
masses), H-bonds by MDAnalysis `HydrogenBondAnalysis` in RAM-bounded batches.
Workers stream `*_upto_<frame>.csv`, merged into `*_merged.csv` per chain pair.
`build_results` aggregates them per contact type (by folder) into
`*_final.csv`; then persistency matrices, energies (streamed per frame, constant
memory) and threshold scans. Selections assume atom names are unique within a
residue.

**network** (`services/network.py`, `utils/network_utilities.py`,
`utils/matrix_utilities.py`) combines the per-type matrices per source
(`load_combined_per_type`: persistency by max, energy by sum), thresholds them,
builds `networkx` graphs labelled `<source>` / `<source>_unweighted`, computes
centralities and communities, writes CSVs, `.pml` scripts and matrices. With a
frame range or windows it rebuilds from the per-frame CSVs instead.
`CompareVariants(Network)` reuses `load_static`.

**compare** (`utils/validation_utilities.py`) only reads matrices.

## Extending myPSNs

**An option.** (1) Add a `Field` to the stage's model in `config.py`, with its
default and description: the model is the source of truth for type, default
and constraints. (2) Add a `typer.Option(None, ...)` to the command, put it in
the overrides dict, the `print_summary` block and the service call (and in
`run_cmd.py` if `mypsns run` forwards it). (3) Document it in the command's
table above, and add a test if it changes behaviour. A new `q_interactors`
parameter also goes into `CHARGE_BASED_PARAMETERS`
(`utils/interactors_utilities.py`); the tests check the tuple, the config and
the forwarding.

**A command.** Create `commands/<name>_cmd.py` with one function whose
parameters are `typer.Option`s, import it in `commands/__init__.py` (`__all__`)
and register it in `cli.py` with `app.command(name=..., help=...)`.

**A pipeline stage.** Add its model to `PipelineConfig`, a service class, a
command (which records its own section with `dump_run_config(...,
sections=["<stage>"])` and skips when its output exists), a step in
`run_cmd.py` if `run` chains it, a status check in `report_cmd.py`, and its
reference section and run-folder layout in this README.

## Interactor rules

`named_interactors` gives each atom zero or more roles, which `contacts` turns
into selections: `HC` (hydrophobic side chain), `HDon` / `HAcc` (H-bond donor /
acceptor heavy atom), `HBond` (a donor-bound hydrogen), `Basic` / `Acidic`
(charged-group member). The defaults reproduce PyInteraph 1's
`hydrogen_bonds.ini`, `charged_groups.ini` and hydrophobic residue list:

```yaml
interactors:
  named_interactors_rules:
    HC:
    - residue: [ALA, VAL, LEU, ILE, PHE, PRO, TRP, MET]
      not atom_name: [N, CA, C, O, H, H1, H2, H3, O1, O2, OXT, OT1, OT2]
    HDon:
    - atom_name: [OH, OG, NE2, OG1, NE, N, ND1, NZ, NH1, NH2, OW, ND2, SG, OH2, NE1]
    HAcc:
    - atom_name: [OH, OG, OD1, OD2, OG1, O, ND1, NE2, OE2, OW, SG, OE1, OH2, SD]
    Basic:
    - resname: [NH3+, GLY-NH3+, NHT3+, PRO-NH2+]
      atom_name: [N, H1, H2, H3, HT1, HT2, HT3]
    - resname: [ARG]
      atom_name: [CD, NE, CZ, NH1, HH11, HH12, NH2, HH21, HH22]
    - resname: [HIS, HSD, HSE, HSP, HSA, HSB, HIS1, HISH, HID, HIE, HIP]
      atom_name: [CG, ND1, HD1, CD2, CE1, NE2, HE2]
    - resname: [LYS, LYSH]
      atom_name: [CE, NZ, HZ1, HZ2, HZ3]
    Acidic:
    - resname: [ASP]
      atom_name: [CG, OD1, OD2]
    - resname: [GLU]
      atom_name: [CD, OE1, OE2]
    - resname: [COO-]
      atom_name: [C, O1, O2, O, OXT, OT1, OT2]
    hb_masses: [1]          # rounded mass of the hydrogens tagged HBond
```

**Grammar.** Each role has a list of rules; a rule holds when all its conditions
hold (AND), a role when any rule holds (OR). A condition is
`attribute: [values]`; `not attribute:` negates it. Attributes are columns of
the per-residue tables: `atom_name`, `atom_type`, `resname` and `residue`.
`resname` is the forcefield name of the atom's own entry, so the atoms a
terminal patch adds or replaces carry the patch name (`NH3+`, `NH2`, `COO-`,
`COOH`); `residue` is the residue the atom belongs to, patch atoms included.
The HC rule uses `residue` because PyInteraph's side chain is every atom of the
residue except its backbone list (above): a C-terminal Leu keeps the `HT2` of a
COOH patch. The charged groups use `resname`, since there the patch is the
group. A rule set in the YAML **replaces** the defaults (copy the whole block
from any `run_config.yaml`), and changing it needs `interactors --overwrite`,
then `contacts --overwrite`.

**Charge gate.** `Basic`/`Acidic` rules propose groups; `--min-residue-charge`
(0.8 e, 0 disables) keeps a group only if the forcefield net charge of its
`(residue, resname)` group clears it with the right sign. Neutral His (HSE,
HSD) is therefore dropped, protonated His (HSP) and Ca²⁺ (+2) pass. PyInteraph
has no gate; turn it off to match a PyInteraph run that counted neutral groups
(the `rec_grk1` preset does).

**Hydrophobic centres of mass** are weighted by forcefield masses (MDAnalysis
guesses elements from names and turns `CAL` into carbon). An atom may hold
several roles: the amide `HN` is `HC` and `HBond`, Met `SD` and Trp `NE1` are
`HC` and acceptor/donor, as in PyInteraph.

**Custom residues.** A recoverin–GRK1 example, a myristoylated glycine and Ca²⁺:

```yaml
HC:
- residue: [ALA, VAL, LEU, ILE, PHE, PRO, TRP, MET]
  not atom_name: [N, CA, C, O, H, H1, H2, H3, O1, O2, OXT, OT1, OT2]
- residue: [GLM]                           # myristoyl tail: not the Gly backbone or the C1=O1 linker
  not atom_name: [N, HN, CA, HA1, HA2, C, O, C1, O1]
Basic:
# … the default Basic rules …
- resname: [CAL]                           # Ca²⁺: +2 clears the charge gate, pairs with Asp/Glu
  atom_name: [CAL]
```

**`q_interactors`** derives roles from charges and connectivity instead:
donors/acceptors from partial charges and element masses, charged groups from
the summed charge of an atom's terminal subtree, HC atoms from low charge and
bond distance to charged/polar atoms. A residue's HC atoms are split into
bond-contiguous clusters of `min_num_heavy_ats`–`max_num_heavy_ats` heavy atoms
(`HC1`, `HC2`, … in `HC_cluster`), each with its own centre of mass. The split
is an exact, deterministic search (fewest size violations, then fewest
clusters, then the alphabetically first); past a work budget (~1 s) it falls
back to a greedy split with a warning. It is not PyInteraph-like: use
`named_interactors` for PyInteraph-comparable results.

## Testing

```bash
pytest                    # everything, ~3–5 min
pytest -m "not slow"      # without the end-to-end run, ~15 s
```

| File | Covers |
|---|---|
| `test_pyinteraph_regression.py` (`slow`) | `mypsns run` on the shipped test case into `tmp_path`, compared with PyInteraph through the `pyinteraph_testcase` preset: tolerances, no PyInteraph-only edge, identical SB/HC edge sets |
| `test_validation_utilities.py` | `compare` logic and CLI, presets complete and inside the repo |
| `test_matrix_utilities.py` | thresholding, graph ↔ matrix, combining targets |
| `test_build_results.py` | inter-chain contact files reach the aggregation |
| `test_run_config_dump.py` | each stage updates only its sections of `run_config.yaml` |
| `test_named_interactor_defaults.py` | one copy of the default rules; configs get a deep copy |
| `test_terminal_hc.py` | HC side chain of terminal residues (patch atoms) |
| `test_terminal_patches.py` | `.tdb` patches through ff-reader and recognition, `structure_atom_name` |
| `test_hc_clusters.py` | `q_interactors` HC clusters, identical across hash seeds |
| `test_charge_based_defaults.py` | q-scheme defaults from the config and their CLI forwarding |
| `test_na_values.py` | atoms and residues named `NA` survive every CSV read |
| `test_config_forcefield.py`, `test_stages_without_forcefield.py` | which commands need a forcefield path |

No test reads outside the repository. Reference data live in `tests/data/`
(`pyinteraph_testcase/`, `rec_grk1/`, each with a README on provenance);
`tests/conftest.py` holds synthetic residues. `network` and `compare-variants`
have no end-to-end test: check them by running the CLI on a finished run folder
copied (or symlinked) into a scratch outpath, because every command writes into
its outpath.

## Validation against PyInteraph 1

| Preset | Data | H-bonds | HC | SB |
|---|---|---|---|---|
| `pyinteraph_testcase` | PyInteraph's CHARMM example; config = myPSNs defaults | 3.93 % (282 / 309 / 282) | 0.11 % (79 / 79 / 79) | 0.07 % (21 / 21 / 21) |
| `rec_grk1` | recoverin–GRK1, 1 µs, PyInteraph 1 graphs of 2018 | not run | 0.20 % (322 / 362 / 322) | 0.13 % (126 / 135 / 126) |
| `rec_grk1_r1` | the same, replica r1 | 4.28 % (812 / 1045 / 812) | 0.13 % (223 / 244 / 223) | 0.10 % (99 / 112 / 99) |

Weight difference, then edges PyInteraph / myPSNs / shared; checked 2026-10-26.
Every PyInteraph edge is found in all three. On rec_grk1 no SB/HC pair differs by
more than 0.05 points, PyInteraph's 1-decimal rounding.

**Re-running.** The test case is the slow test (`pytest -m slow`), or by hand as
in the [quick start](#quick-start) with its config:

```bash
D=tests/data/pyinteraph_testcase
mypsns run --config $D/config.yaml -t $D/structures/sim.prot.gro -x $D/structures/traj.xtc \
           -ff $D/forcefields/charmm36-apr2019_vm.ff -o /tmp/tc
mypsns compare --preset pyinteraph_testcase --computed-dir /tmp/tc/Adjacency/named_interactors
```

rec_grk1 needs the lab trajectory (`rec_grk1_1micros.xtc`, `rec_grk1_nolip.pdb`):
`mypsns run --config tests/data/rec_grk1/config.yaml -t … -x … -ff … -o <out>`
(`config_r1.yaml` for replica r1), then `compare --preset rec_grk1` (or
`rec_grk1_r1`). That 2018 run used non-default settings: its own charged groups
(no Arg; Ca²⁺ and neutral COOH termini included; charge gate off), the myristoyl
in HC and an HC cutoff of 5.5 Å. `tests/data/rec_grk1/README.md` explains how
they were read back from its output.

**H-bonds are engine-limited.** myPSNs uses MDAnalysis `HydrogenBondAnalysis`
(D–A ≤ 4.5 Å, D–H–A ≥ 120°); PyInteraph 1 has its own geometric test, which
measures H···A. The edge sets agree, and the remaining ~4 % is small per-frame
scatter between the two detectors that no parameter removes. One parameter did
matter: MDAnalysis pairs donors with hydrogens within `d_h_cutoff`, whose 1.2 Å
default is shorter than the Cys S–H bond (1.32 Å). myPSNs uses 1.4 Å; with 1.2 Å
the test case loses 3 PyInteraph edges (Cys86 thiol donors) and sits at 4.56 %.

**Macro-IIN.** The unweighted persistency network reproduces PyInteraph's
`filter_graph` macro-IINs: exactly when fed PyInteraph's own graphs (a
recoverin–GRK1 macro-IIN at 20.5 %, 331 / 331 edges), and on the test case
192 / 197 / 192 at 20 % from myPSNs' own contacts (the 5 extra pairs are H-bonds
near the cutoff).

## Conventions

- **Defaults are validated.** The contact cutoffs and the default rules
  reproduce PyInteraph 1; change them only with the regression test passing.
- **`na_values`.** Every `pd.read_csv` passes the shared `DEFAULT_NA_VALUES`
  (which excludes `"NA"`) **and** `keep_default_na=False`, otherwise sodium or a
  heme nitrogen named `NA` reads as missing. `test_na_values.py` checks every call.
- **One copy of the default rules**, `NAMED_INTERACTOR_DEFAULTS` in `config.py`;
  `InteractorsConfig` hands out a deep copy, so never mutate the constant.
- **Determinism.** Nothing iterates a set of atom names into a result (HC
  clusters, chain bounds): sort first.
- **`run_config.yaml`.** Each command records only its own section(s); `run`
  records its four stages. A command that skips or fails still leaves the file
  describing the run.
- **Stage outputs are the interface.** Stages communicate only through the run
  folder, so a stage can be re-run alone; keep file names and columns stable,
  or update every reader and this README.
- **Comments explain why.** Docstrings and comments are long and say why a
  choice was made (often a PyInteraph detail); keep that style.

## Known issues and limitations

- **H-bonds** differ from PyInteraph by ~4 % (detector, above); SB and HC by ~0.1 %.
- **`q_interactors`** is less validated than `named_interactors`. On the test
  case its HC set also includes Tyr, Thr CG2, Cys and the CB of Arg/Gln/His,
  drops Trp CD1/CE2/NE1, and uses per-cluster centres of mass, so the HC network
  is much denser (207 edges vs 79). Some residues get single-heavy-atom clusters
  where the bonds allow nothing else (e.g. Leu around CG).
- **Water parsing** from forcefields can produce inconsistent duplicate rows;
  water is matched against the exact table `mypsns/resources/sol.csv` instead.
- **Energy networks** count attractive interactions only unless
  `--energy-use-abs`.
- **Hand-edited terminal patches** must list the bonds of the atoms they add;
  recognition warns when a patch leaves an atom unbonded.
- **Selections** assume atom names are unique within a residue.
- **Run folders made before 2026-10-26** may miss inter-chain H-bonds and
  terminal HC atoms; see [`CHANGELOG.md`](CHANGELOG.md).
