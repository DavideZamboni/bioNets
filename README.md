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
| `mypsns accept-topology` | Promote a proposed forcefield candidate for residues `topology` could not match (see [Unmatched residues](#unmatched-residues-mypsns-accept-topology)).
| `mypsns interactors` | Assign interaction roles to topology atoms.
| `mypsns contacts`    | Compute contacts across an MD trajectory.
| `mypsns network`     | Build residue-interaction networks and render PyMOL sessions.
| `mypsns compare`     | Compare residue matrices with reference ones (e.g. PyInteraph graphs or macro-IIN); see [Checking a run against PyInteraph](#checking-a-run-against-pyinteraph).
| `mypsns compare-variants` | Rank network metrics against ClinVar pathogenicity (experimental; needs the `variants` extra).
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
│   ├── network.py          PSN construction + PyMOL
│   └── compare_variants.py ClinVar ranking of network metrics (experimental)
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
  set used as the default for `InteractorsConfig.named_interactors_rules`. It is the only
  copy: `utils/interactors_utilities.py` imports it, and
  `tests/test_named_interactor_defaults.py` checks that no second copy appears.

## Services

The `mypsns/services/` folder holds the modules that drive the pipeline.
Expected execution order:

1. `ff_reader` — parse the GROMACS forcefield
2. `topology_recognition` — match topology residues to forcefield entries
3. `interactors` — reconstruct whole-system topology + assign interaction roles
4. `contacts` — compute contacts over a trajectory and build matrices
5. `network` — build residue-interaction networks (standalone; run after `contacts`)
6. `compare_variants` — rank network metrics against ClinVar (experimental; `mypsns compare-variants`)

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
- optional: `na_values`

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
  any residue is unmatched). Only the outcome of the current run is kept: a
  `complete_…` or `unmatching_…` file left by an earlier run is removed.

  If **any** residue topology is left unmatched, the stage logs which residues
  failed (and how many atom names the best candidate shares with the structure)
  and raises `TopologyAssociationError`, halting the pipeline before
  `interactors` (which would otherwise fail later with a confusing missing-file
  error). The message prints the `mypsns accept-topology` commands that resolve
  it: see [Unmatched residues](#unmatched-residues-mypsns-accept-topology).

  > Non-standard residues (e.g. N-myristoylated Gly `GLM`, a metal ion `CAL`)
  > must have a forcefield entry to be matched. If recognition fails for one,
  > raise `--max-to-check` toward `-1` (inspect every ff entry) — large custom
  > residues can sort far down the name-similarity list.

#### Unmatched residues: `mypsns accept-topology`

Recognition only accepts a forcefield entry whose atom names match the
structure's **exactly**. A residue that differs by a few names (a hydrogen
naming convention, an N-terminus written as `HT1–HT3` while the forcefield's
`NH3+` patch calls them `H1–H3`) halts the stage, and recognition writes the
closest forcefield entries to `possible_topologies/` as candidates.
`accept-topology` promotes the candidate you choose, so the pipeline can continue.
The fix in the forcefield or the structure (option 2 in the error message) is
still the cleaner one when it is possible.

The topology names used below (`MET`, `GLY`, `MET_1`, …) are the residue
topology names from the error message: the resname, plus `_<n>` for its
second, third, … atom set (e.g. a terminal copy of the residue).

**1. Recognition halts.** `mypsns topology` (or `mypsns run`) stops with a
`TopologyAssociationError` that lists each unmatched topology, how many atom
names its best candidate shares with the structure, and the commands for the
next steps:

```text
Topology recognition could not associate 1 residue topology/ies to the forcefield:
  - MET (topology MET): no exact match; the best of 4 candidate(s) shares 16/19 atom names
...
  1. promote a candidate. List the candidates and how their atoms pair with the structure's:
       mypsns accept-topology -o run/
     then accept one per topology (by index, or a hand-edited CSV path):
       mypsns accept-topology -o run/ MET=0
```

> **[Screenshot 1 — to add]** the halted `mypsns topology` run in the terminal,
> showing the unmatched residue and the suggested `accept-topology` commands.
> Suggested file: `docs/img/accept-topology-1-halt.png`.

**2. List the candidates.** Run the command without choices:

```bash
mypsns accept-topology -o run/
```

```text
             MET  (MET, 2 residue(s)) — unresolved
┏━━━┳━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ # ┃ forcefield residue(s) ┃ renames (forcefield → structure) ┃
┡━━━╇━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ 0 │ GLY-NH3+, MET         │ H1→HT1, H2→HT2, H3→HT3           │
│ 1 │ NH3+, MET             │ H1→HT1, H2→HT2, H3→HT3           │
│ 2 │ GLY-NH3+, DMET        │ H1→HT1, H2→HT2, H3→HT3           │
│ 3 │ NH3+, DMET            │ H1→HT1, H2→HT2, H3→HT3           │
└───┴───────────────────────┴──────────────────────────────────┘
```

Each row is one file, `possible_topologies/<topology>_uncertain_atoms_<#>.csv`:

- **forcefield residue(s)**: the forcefield entries the candidate is built
  from, i.e. the residue and the terminal patch applied to it. **This is what
  you choose on.** Here the structure residue is an N-terminal L-Met, so #1 (the
  generic `NH3+` patch on `MET`) is right; #0 uses the glycine N-terminal patch
  (different charges on CA/HA) and #2/#3 are D-Met.
- **renames**: atoms whose names differ, paired with the structure's atoms the
  way recognition scored the candidate. Atoms with the same name pair with each
  other, and the rest pair in file order. A red `!` marks a pair whose masses
  differ: that pairing is wrong and the candidate needs a hand-edited copy
  (step 4).
- `identical to #n` marks a candidate whose file is identical to an earlier
  one. Recognition often proposes the same topology several times (one unpatched
  copy for each terminal patch that does not apply), so skip those.
- A topology with no candidate at all (none of the forcefield entries inspected
  has its atom count) can only be resolved with a CSV you write (step 4), unless
  raising `--max-to-check` toward `-1` and re-running `topology --overwrite`
  finds one.

To inspect a candidate in full (atom types, charges, bonds), open its CSV.

> **[Screenshot 2 — to add]** the candidate table from `mypsns accept-topology -o run/`,
> and next to it the chosen `possible_topologies/MET_uncertain_atoms_1.csv`
> opened in a spreadsheet. Suggested file: `docs/img/accept-topology-2-candidates.png`.

**3. Accept one candidate per topology.** Check with `--dry-run` first, then apply:

```bash
mypsns accept-topology -o run/ MET=1 --dry-run   # shows the renames, writes nothing
mypsns accept-topology -o run/ MET=1             # several at once: MET=1 GLY_1=0 …
```

For each choice this:

1. copies the candidate to `found_topologies/<topology>_ff_atoms.csv`, renaming
   its differently named atoms to the structure's names (`H1→HT1`, …) in
   `atom_name` **and** in the `bonds` lists. Contacts later joins this table
   with the structure on atom names, so the forcefield names must not survive;
2. fills that row's `ff_atoms` in `unmatching_residue_topologies.csv` and records
   the source in a new `accepted_from` column;
3. once **every** row has `ff_atoms`, writes `complete_residue_topologies.csv`,
   the file `interactors` reads.

All choices are checked before anything is written. If one is invalid (unknown
topology, missing candidate, atoms that cannot be paired), nothing changes and
every problem is listed. Each run is logged to
`Topology_recognition/accept_topology.log`.

> **[Screenshot 3 — to add]** the "Accepted" table and the
> "every residue topology resolved" line. Suggested file:
> `docs/img/accept-topology-3-accepted.png`.

**4. Hand-edited candidate (when no candidate fits as is).** Copy the closest
candidate (or any topology CSV with the same columns), edit it, and pass its
path instead of an index:

```bash
cp run/Topology_recognition/possible_topologies/MET_uncertain_atoms_1.csv met_fixed.csv
# edit met_fixed.csv: atom_name must use the structure's names, and every
# bonds entry that names a renamed atom must be updated too
mypsns accept-topology -o run/ MET=met_fixed.csv
```

The file must keep all ten columns (`rtp_file, resname, atom_name, atom_type,
charge, chargegroup, bonds, sigma, epsilon, mass`). A `bonds` entry that names
an atom the residue does not have raises a warning: usually `atom_name` was
edited but `bonds` was not. (`-X` / `+X` entries name atoms of the previous / next
residue and are not checked.)

> **[Screenshot 4 — to add]** `met_fixed.csv` in a spreadsheet with the edited
> `atom_name` / `bonds` cells highlighted. Suggested file:
> `docs/img/accept-topology-4-hand-edit.png`.

**5. Continue the pipeline.**

```bash
mypsns interactors --config run.yaml   # then: mypsns contacts --config run.yaml
# or: mypsns run --config run.yaml     (skips topology: complete_residue_topologies.csv exists)
```

If `Interactors/` already had outputs (you changed an earlier choice), re-run
`mypsns interactors --overwrite` and then `mypsns contacts --overwrite`.

| Argument / option | Default | Meaning
|-------------------|---------|--------
| `TOPOLOGY=CHOICE` (repeatable) | — | `CHOICE` is a candidate index (`possible_topologies/<TOPOLOGY>_uncertain_atoms_<CHOICE>.csv`) or the path of a CSV. No choices: list the unresolved topologies and their candidates.
| `-o` / `--outpath` | — | Run folder (contains `Topology_recognition/`).
| `--config` | — | YAML config. Only `outpath` and `topology_recognition.na_values` are used; `-o` wins.
| `--dry-run` | off | Check the choices and show the renames; write nothing.

Things to know:

- **Re-running `topology --overwrite` discards acceptances.** Recognition rewrites
  the tables and removes the `complete_residue_topologies.csv` that acceptance
  wrote, so accept again.
- **The interactor rules see the structure's atom names.** After acceptance the
  residue's atoms carry the structure's names, and the `named_interactors` rules
  are name-based. A name the rules treat specially changes the residue's roles:
  for example the HC rule excludes the AMBER amide `H` from the side-chain CoM
  but keeps the CHARMM `HN`. Check the residue's row in
  `Interactors/found_topologies/<variant>.csv` after `interactors`.
- **Validated end to end** on the PyInteraph test case: with every Gly `HA1`
  renamed to `HA3` in the structure, recognition halts on GLY. After
  `accept-topology GLY=0`, `interactors` and `contacts` give SB / HC / H-bond
  adjacency and energy matrices byte-identical to the unmodified run.

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
matches (≈279/282 on the test set); the residual ~4.5 % is mostly small per-frame
persistence scatter (~1.8 pts/edge) from the different detectors, which **cannot be
tuned away**. SB and HC, by contrast, reach ~0.1 %.

One part of the gap *is* fixable: the 3 edges only PyInteraph finds all involve
the Cys thiol donor. MDAnalysis pairs a donor with its hydrogens only within
`d_h_cutoff` (default 1.2 Å), and the S–H bond is ~1.32 Å, so SG–HG1 H-bonds
are never detected. With 1.4 Å they come out at PyInteraph's exact persistence
(Cys86→Asp118 51.7 %, →Pro119 22.5 %). myPSNs does not pass `d_h_cutoff` yet —
see *Open issues*.

### Checking a run against PyInteraph

`mypsns compare` compares residue × residue matrices type by type: edges
(truth / computed / shared), entries whose values differ, and the % weight
difference `100·Σ|computed−truth| / Σtruth` on the strict lower triangle
(diagonal ignored — PyInteraph's `hb-graph.dat` keeps intra-residue H-bonds
there). With `-o DIR` it writes `possible_sele_mismatch.csv` (edge in one matrix
only), `possible_params_mismatch.csv` (both edges, values differ) and
`comparison_summary.csv`; `--labels adjacency_labels.dat` adds residue labels.

```bash
# PyInteraph-named reference folder vs a contacts output folder
mypsns compare --truth-dir ref/ --truth-suffix _1micros \
               --computed-dir run/Adjacency/named_interactors -o cmp/
# any pair of matrices, e.g. a PyInteraph macro-IIN vs the unweighted network
mypsns compare --truth macroIIN=macro.dat \
  --computed macroIIN=run/Network/named_interactors/matrices/persistency_unweighted_network_adjacency.dat
# fail (exit 1) above a tolerance — usable as a CI check
mypsns compare --preset pyinteraph_testcase --max-pct-diff SB=0.2 --max-pct-diff HC=0.2
```

`--preset pyinteraph_testcase` / `rec_grk1` reproduce the two lab datasets
under `$MYPSNS_VALIDATION_ROOT` (default `/home/davide/storage/lab_dati`).
`pytest` runs the unit tests plus a regression test on the PyInteraph test
case (SB/HC ≤ 0.2 %, H-bonds ≤ 5 %), which is skipped when that data is absent.

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
one, e.g. ions) and draws network edges as solid sticks scaled by edge weight.
In `node_metrics/` every edge is **grey**. In `community_metrics/` an edge
between two residues of the same community takes that community's colour (the
colour of its spheres) and only edges that cross communities stay grey, so each
community reads as a coloured sub-network. The palette has 16 colours and wraps
around, so with more than 16 communities two of them share a colour.
Residue selectors use the PDB residue numbering from the topology, so
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

`persistency_unweighted_network_adjacency.dat` is the analogue of PyInteraph's
macro-IIN (`filter_graph -d hb -d hc -d sb -t <cutoff>`). Fed PyInteraph's own
per-type graphs, it reproduces a real PyInteraph macro-IIN exactly
(recoverin–GRK1, cutoff 20.5 %: 331/331 edges). One semantic difference:
`filter_graph` keeps values **>** the cutoff, the unweighted network keeps
values **≥** `--min-persistency`, so a pair sitting exactly on the cutoff is an
edge only here (no such pair on the validation data).

### Output layout

```
<outpath>/Network/
├── network.log
└── <inttype>/
    ├── csv/
    │   ├── nodes_<label>_<method>.csv     (centrality + community per node)
    │   ├── edges_<label>.csv
    │   └── shortest_paths_<label>.csv     (with --shortest-paths)
    ├── node_metrics/<label>_<metric>.pml        (grey edges)
    ├── community_metrics/<label>_<method>.pml   (edges coloured by community, cross-community grey)
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

## Experimental: `compare-variants` (ClinVar variant comparison)

Overlays ClinVar pathogenicity onto the networks and ranks node/edge metrics by
how well they separate pathogenic missense sites from the rest of the protein
(AUROC / AUPRC / Mann–Whitney). Standalone: run it after `network`. It needs
internet access and the optional `variants` extra (`requests` + `biopython`).
Without them the command stops with an install hint (exit code 3).

```bash
pip install -e '.[variants]'
mypsns compare-variants --config <outpath>/run_config.yaml \
    --gene RCVRN --reference NP_002894.1
```

What it does:

1. **Rebuilds the static graphs** from `Adjacency/` with the `network:` settings
   (targets, sources, weighting, thresholds, `min_degree`, frame range), using
   the same code path as `network`. If you point `--config` at the
   `run_config.yaml` that `network` wrote, it ranks exactly those graphs. Network
   flags given on the command line override that section. Centralities are
   computed in memory. Nothing under `Network/<inttype>/csv/` is rewritten.
2. **Queries ClinVar** (`<gene>[gene]`, optionally narrowed to a GRCh38
   `--chrom/--start/--stop` window) through NCBI E-utilities and keeps the
   missense records. Each protein position gets a label: `non_benign` if it has
   at least one pathogenic or likely-pathogenic record, `benign` if it has only
   benign or likely-benign ones, and `vus` otherwise.
3. **Maps positions onto the structure** by aligning the RefSeq protein
   (global, BLOSUM62) to every protein chain (`Entity == 'P'`). The best chain
   above `--min-identity` is taken as the gene's chain, which handles offsets,
   gaps and complexes.
4. **Scores each graph**. Positives are the `non_benign` sites. Negatives are
   every other aligned residue (`--negative rest`) or only the `benign` sites
   (`--negative benign`). Residues that do not align to the reference (partner
   chains, ions, ligands) are not scored, because ClinVar says nothing about
   them for this gene.

### Key options

| Flag | Default | Meaning
|---|---|---
| `--gene` | required | Gene symbol queried in ClinVar.
| `--reference` | required | RefSeq protein the ClinVar `p.` positions refer to: the NP_ of the transcript named in ClinVar titles (e.g. `NM_002903.3` → `NP_002894.1`). A raw protein sequence also works.
| `--chrom` / `--start` / `--stop` | — | Optional GRCh38 window (all three or none). `--gene` alone already returns every record in the gene.
| `--email` | — | Contact email sent to NCBI (env `NCBI_EMAIL`).
| `--api-key` | — | NCBI API key, 10 requests/s instead of 3 (env `NCBI_API_KEY`). **Never written to `run_config.yaml`.**
| `--negative` | `rest` | Negative class: `rest` (every other aligned residue) or `benign`.
| `--min-identity` | 70 | Minimum % identity (over the structure chain) to accept a chain as the gene's.
| `--interactors-type`, `--target`, `--edge-weight-source`, `--energy-use-abs`, `--weighting`, `--min-persistency`, `--min-energy`, `--min-degree`, `--frame-start`, `--frame-stop` | from `network:` | Graph rebuild; same meaning as in `network`.
| `--overwrite` | off | Rewrite this gene's tables if they exist. Without it the stage skips.
| `--dry-run` | off | Validate and print the plan. It also reports a missing `variants` extra.

The ClinVar/scoring options live under `compare_variants:` in the YAML. The
graph-rebuild options live under `network:`.

### Output layout

Every file is prefixed with the gene, so several genes of a complex can share
the folder:

```
<outpath>/Network/
├── compare_variants.log
└── <inttype>/variants/
    ├── <gene>_ranking_<label>.csv     one row per metric: auroc, auroc_oriented,
    │                                  direction, auprc, mwu_p, n_pos, n_neg
    │                                  (best first)
    ├── <gene>_nodes_<label>.csv       per aligned residue: metrics, ref_pos, ref_aa,
    │                                  ClinVar label
    ├── <gene>_clinvar_sites.csv       per missense site: counts, label, chain,
    │                                  resindex, resnum, structure vs reference aa
    └── <gene>_clinvar_variants.csv    raw ClinVar records (dated snapshot)
```

`<label>` follows `network` (`persistency`, `energy`, `*_unweighted`). The two
`clinvar_*` tables are the same in every `<inttype>` folder.

### Notes

- **Empty rankings are a data outcome, not an error.** If no pathogenic missense
  site falls on a modelled residue, every `ranking` table has a header only and
  the command warns. The site tables are still written. Example: RCVRN on the
  rec_grk1 run (October 2026) has 18 missense sites (17 VUS, 1 benign, no
  pathogenic), all mapped to chain A at 87 % identity.
- **ClinVar changes over time.** Keep `<gene>_clinvar_variants.csv` with the
  results.
- Mismatches between `structure_aa` and `ref_aa` in `clinvar_sites` are expected
  for orthologues (e.g. a bovine structure against human RefSeq). Many
  mismatches suggest the wrong reference isoform.
- From Python, it is a `Network` subclass:

```python
from mypsns.services.compare_variants import CompareVariants
cv = CompareVariants(topology_pdb=..., general_outpath=...)   # same constructor as Network
cv.load_static(...)                                          # same arguments as `network`
ranked = cv.rank_against_clinvar(gene="RCVRN", reference="NP_002894.1")
cv.write_rankings(ranked, gene="RCVRN")
```

## Open issues

- **Water parsing** from the forcefield can return inconsistent duplicate rows.
  Workaround: water is matched against a bundled exact-match table
  (`mypsns/resources/sol.csv`) at startup.
- **H-bonds vs PyInteraph** plateau at ~4.5 % (mostly an engine difference —
  see *H-bonds are engine-limited*). SB/HC reach ~0.1 %.
- **Cys thiol H-bonds are missed.** `HydrogenBondAnalysis` is called without
  `d_h_cutoff`, so its 1.2 Å default is shorter than the S–H bond (~1.32 Å) and
  SG–HG1 is never paired as a donor. On the test case this accounts for all 3
  PyInteraph-only H-bond edges. Passing ~1.4 Å fixes it, but changes validated
  output: re-run the PyInteraph comparison when doing so.
- **`network` needs `--config`.** The pipeline config requires `forcefield`, which
  `network` neither uses nor exposes as a flag, so a flags-only
  `mypsns network -t … -o …` fails validation. Pass a YAML with `forcefield:` set.
- **`q_interactors`** is less validated than `named_interactors`. For HC
  classification it uses a `min_bonds_to_polar` / `min_bonds_to_charged` BFS
  rule; `min_bonds_to_polar=1` is needed mainly so ALA's Cβ is classified as HC.
  It is **not** PyInteraph-like: on the PyInteraph test case it agrees with
  `named_interactors` on the HC atoms of ALA/VAL/LEU/ILE/PHE/PRO/MET, but also
  types Tyr (CB + ring), Thr CG2, Cys CB/SG and the CB of Arg/Gln/His as HC,
  and drops Trp CD1/CE2/NE1 (next to the polar NE1). Contacts are measured
  between clusters of ≤ `max_num_heavy_ats` heavy atoms instead of whole side
  chains, so the HC network is far denser (207 vs 79 edges, 76 shared; mean
  persistence of shared edges 55 % vs 25 %). Use `named_interactors` for
  PyInteraph-comparable results. The
  `q_interactors` HC clustering (`HC_cluster`) is **not reproducible between
  runs**: it depends on Python's string-hash seed, so two `interactors` runs on
  the same input can split LEU/PHE/TRP side chains into different clusters (set
  `PYTHONHASHSEED` to compare runs). `named_interactors` is not affected.
- **Repulsive interactions** are excluded from energy-weighted networks by
  default (attractive-only). Use `--energy-use-abs` to include them.
- **Terminal patches can leave stale bond names.** When recognition applies a
  `.tdb` patch that renames atoms, the bonds of the other atoms keep the old
  names. Example: an exactly matched `NHT3+` Met N-terminus whose `N` lists
  `H1 H2 H3` while its atoms are `HT1–HT3`, so the N–HT bonds are missing.
  `accept-topology` warns about such names in the topology it accepts.
- **Tests cover helpers and the PyInteraph regression only.** `pytest` checks
  the matrix/network helpers, `mypsns compare` and the test-case validation
  (skipped without the lab data). Pipeline stages are still checked by running
  the CLI on a real run folder.
- **`rec_grk1` does not match its PyInteraph graphs** (H-bonds 33 %, HC 60 %,
  SB 44 % weight difference; values correlate, r ≈ 0.8–0.9). The PyInteraph
  settings of that run are unknown, so it is not a validation dataset.

## TODO

The maintained development TODO list lives in [`CLAUDE.md`](CLAUDE.md#todo).
