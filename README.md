energy + parser + 

necessario sopprimere exception mdanalysis:
         INFO     The attribute(s) types have already been read from the topology file. The guesser will only   
                  guess empty values for this attribute, if any exists. To overwrite it by completely guessed   
                  values, you can pass the attribute to the force_guess parameter instead of the to_guess one   
         INFO     There is no empty types values. Guesser did not guess any new values for types attribute      
         INFO     attribute masses has been guessed successfully.                                               
         INFO     SB worker: frames 75001–87500                                                                 
         INFO     starting core_sb                                                                              
         INFO     SB worker: frames 87501–100000                                                                
         INFO     The attribute(s) types have already been read from the topology file. The guesser will only   
                  guess empty values for this attribute, if any exists. To overwrite it by completely guessed   
                  values, you can pass the attribute to the force_guess parameter instead of the to_guess one   
         INFO     There is no empty types values. Guesser did not guess any new values for types attribute      
         INFO     attribute masses has been guessed successfully.                                               
         INFO     The attribute(s) types have already been read from the topology file. The guesser will only   
                  guess empty values for this attribute, if any exists. To overwrite it by completely guessed   
                  values, you can pass the attribute to the force_guess parameter instead of the to_guess one   

metti comando x identificazione forcefields in determinato path
screenshots pipeline possible_topologies -> found manuale

doppia threshold x ILE / LYS (Hbonds vs Acidic / Basic) fatto
interaz repulsive? momentaneamente no

rivedi pdb output (catene, ...)

pymol
- rappresentazione prop nodi come spehers con vdw centrati su CA modificabile (proporzionale a metrica)
- edges con sticks 
NW COMUNITÀ STICKS CON COLORI DIVERSI (solo edges nel network)
1 FILE = 1 METICA

globalprotect: von.univr.it

  - min_bonds_to_polar: minimum number of bonds separating the possibly HC forming atom from a charged (with partial charge > max_abs_charge) atom
                          to classify it as HC

  SET TO 1 mainly for ALA (otherwise no HC)



# myPSNs

## Setup
NEW WORKFLOW:
cd <proj_root/myPSN>
conda env create -f environment.yml
conda activate mypsns

```bash
conda create -n mypsns python=3.11 numpy=1.26 pandas=2.2 psutil networkx=3.2 matplotlib=3.8 -c conda-forge
conda activate mypsns
pip install MDAnalysis
pip install --upgrade pip setuptools wheel
cd <proj_root>
pip install -e .
```

## Run

```bash
python -c "from mypsns.main import main; main()"

```

## Open issues

- H2O parsing from the forcefield returns inconsistent rows:
  ```
  merged   HOH   OW    OT   -0.834   0   HW1 HW2   0.315057422683   0.63639     15.9994
  merged   HOH   OW    OT   -0.834   0   HW1 HW2   0.315057422683   0.636386     9.9514
  merged   HOH   HW1   HT    0.417   0   OW HW2    0                0            4.032
  merged   HOH   HW1   HT    0.417   0   OW HW2    0                0            1.008
  merged   HOH   HW1   HT    0.417   0   OW HW2    0.0400013524445  0.19246      4.032
  merged   HOH   HW1   HT    0.417   0   OW HW2    0.0400013524445  0.19246      1.008
  merged   HOH   HW2   HT    0.417   0   OW HW1    0                0            4.032
  merged   HOH   HW2   HT    0.417   0   OW HW1    0                0            1.008
  merged   HOH   HW2   HT    0.417   0   OW HW1    0.0400013524445  0.19246      4.032
  merged   HOH   HW2   HT    0.417   0   OW HW1    0.0400013524445  0.19246      1.008
  ```
  As a temporary workaround, water is added to the resources DB for exact-match lookup at app start: `./resources/sol.csv`.

## General notes

### Project structure

```
mypsns/
├── main.py
├── parser_dispatcher.py
├── services/
│   ├── ff_reader.py
│   ├── topology_recognition.py
│   ├── interactors.py
│   └── contacts.py           
├── utils/
└── resources/
```

### Modules

- **`main.py`** — application entry point. Calls `parser_dispatcher.py`.

  TODO: handle optional parameters. Most important:
  ```python
  na_values: Optional[List[str]] = [
      "", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN",
      "-NaN", "-nan", "1.#IND", "1.#QNAN", "<NA>", "N/A",
      "NULL", "NaN", "n/a", "nan", "null"
  ]
  ```
  Passed to every `pandas.read_csv()` across the app to identify NA/NaN entries.
  Note: `"NA"` has been removed from the default list to preserve sodium atom names.

- **`parser_dispatcher.py`**

  Parameters:
  - **required**:
    - `topology` (str): topology file in any MDAnalysis-compatible format
    - `general_output_folder`: created if it does not exist
    - `ff_folder`: forcefield folder. Output goes to a folder named after `ff_folder` with `ff.` or `.ff` removed and `-` replaced by `_`. If `final_ff_residues.csv` is already present in `ff_folder`, it is treated as an `ff_reader` output folder and `ff_reader` is skipped.

  Calls service-specific processes (the `.py` files in `mypsns/services/`).
  TODO: terminal command parsing.
  TODO: handle optional parameters.

### Services

Folder containing the modules that drive the application pipeline. Expected execution order:

1. `ff_reader` — parse the forcefield
2. `topology_recognition` — match residues in topology to forcefield entries
3. `interactors` — assemble whole-system topology + assign interaction roles
4. `contacts` — compute contact statistics over a trajectory

#### `ff_reader.py`

Parameters:
- **required**:
  - `inpath`: ff_folder from the parser
  - `outpath`: renamed folder from the parser (e.g. `charmm36_apr2019_vm` from `charmm36-apr2019_vm.ff`)
- **optional**:
  - `na_values`: see `main.py`
  - `nonbonded_name` (default `"ffnonbonded.itp"`): file in the input folder containing nonbonded interaction parameters
  - `additional_residues` (default `True`): whether to merge in extra residues from `resources/`

GROMACS forcefield parser.

1. Reads all residues defined in `.rtp` files. For each atom records:
   - `rtp_file`: name of the source `.rtp` file
   - `resname`: residue name (from rtp)
   - `atom_name`: atom name (from rtp)
   - `atom_type`: atom type (from rtp)
   - `charge`: partial atomic charge (from rtp)
   - `chargegroup`: integer grouping atoms whose total charge sums to an integer. Used by GROMACS to speed up Coulombic interactions when PME is off. (from rtp)
   - `sigma`: σ in `C12 = 4εσ¹²` (from `ffnonbonded.itp`)
   - `epsilon`: ε in `C6 = 4εσ⁶` (from `ffnonbonded.itp`)
   - `mass`: atomic mass (from `ffnonbonded.itp`)

2. Reads termini definitions from `.n.tdb` and `.c.tdb` files.
   **The naming of `.rtp` and `.tdb` files must correspond.** E.g. if `merged.rtp` contains residues that can have termini, then `merged.n.tdb` / `merged.c.tdb` must also exist.

   The atom-definition format is the same as in `.rtp`, with the following bond-entry conventions:
   - atoms starting with `-`, present in every bond entry of a ter, are the ones to be **found and removed** from the original residue definition
   - atoms starting with `r`, present in only one specific atom of the ter entry, are the ones to be **found and replaced** with the atom defined in the ter entry; bonds are adjusted accordingly

3. Reads σ, C6 and epsilon C12 from the `.itp` file specified by `nonbonded_name`.
   Reminder: GROMACS uses `U(r) = C12/r¹² − C6/r⁶`.

4. Produces output files:
   - `final_ff_residues.csv`: `residues.csv` LEFT JOIN `atomtypes_se.csv` on `atom_type` (adds `mass`, `sigma`, `epsilon`)
   - `final_ff_ters.csv`: `ters.csv` LEFT JOIN `atomtypes_se.csv` on `atom_type` (adds `sigma`, `epsilon` — `mass` already present in `ters.csv`)

Expected input files:
- `*.rtp`: GROMACS residue topology files
- `*.tdb`: GROMACS termini database files (`.n.tdb` and `.c.tdb`), naming matched to the corresponding `.rtp`
- `.itp` (default `ffnonbonded.itp`): nonbonded parameters



topology_recognition.py
Reads the input topology and matches each residue in it to a forcefield entry produced by ff_reader. Handles cases where the same residue name appears with different atom sets (e.g. terminal vs. internal residues) by treating each unique atom set as a distinct residue topology.
Parameters:

required:

topology (str): topology file (PDB, GRO, …)
outpath (Path): output folder
ff_path (Path): forcefield folder produced by ff_reader


optional:

na_values (List[str]): see main.py
outlog (str | Path): folder for log files. Defaults to outpath.



The constructor immediately calls parse_initial_topology_info, which extracts everything MDAnalysis exposes per atom and writes:

initial_topology_info.csv: per-atom table with whichever of these columns the topology supports — chainID, resindex, resid, resname, index, id, name, type, element, mass, charge. Columns for which every atom raises NoDataError are omitted.

Note: available fields depend on the topology format — see https://userguide.mdanalysis.org/stable/topology_system.html#format-specific-attributes. Special case: if the topology supplies type but not element and the longest type string is ≤2 characters, the type values are copied into element (handles PDB files where atomic symbols land in the type column).
Three public methods:

get_restypes(verbose=False) — group atoms in the topology into unique residue topologies.
For each residue name, iterates over all residue instances and groups them by their exact atom-name set. Multiple distinct atom-name sets → multiple "topologies" for the same residue name (e.g. N-terminal ALA differs from internal ALA in its atom list).
Writes to <outpath>/resid_topologies/:

<resname>[_<n>]_rixs.txt: space-separated list of Resindex values that share this residue topology
<resname>[_<n>]_atoms.csv: atom names and masses defining this topology

The _<n> suffix is omitted for the first topology of a given resname; subsequent variants are _1, _2, … (e.g. ALA_atoms.csv for the first ALA topology found, ALA_1_atoms.csv for the second).
Also writes <outpath>/different_residue_topologies.csv summarising the index. Columns: resname, atoms (atoms CSV filename), rixs (rixs txt filename).
get_ff_parameters(verbose=False, maxToCheck=-1) — associate each unique residue topology with a specific forcefield entry.
Parameters:

verbose: write a separate residues_to_ff_association.log reporting the matching process
maxToCheck (default -1): max number of forcefield entries to inspect per topology residue. -1 = inspect all. Forcefield residue names are sorted by similarity (via difflib.SequenceMatcher) to the topology residue name, so smaller values still cover the closest candidates first.

Two-pass matching:

Pass 1 (exact match): for each candidate ff residue with the right atom count, check whether its atom-name set matches the topology's exactly. First exact match wins.
Pass 2 (ter substitutions): if Pass 1 fails, apply each .tdb ter substitution to every candidate ff residue and re-check. Atoms prefixed with - in the ter entry are removed; atoms prefixed with r are replaced; bonds are adjusted accordingly. Used to recognise N-/C-terminal variants and similar.

Two outcome buckets per topology residue:

found_residues: exact atom-name match (with or without ter substitutions)
best_choices: no exact match; ranked by atom-name overlap and, where overlap is incomplete, by mass agreement on the unmatched atoms

Outputs to <outpath>/:

found_topologies/<resname>_ff_atoms.csv: forcefield topology (atom names, types, charges, bonds, masses, σ, ε) for each successfully matched residue
possible_topologies/<resname>_uncertain_atoms_<n>.csv: top-scoring candidates for each unmatched residue (one file per candidate)
complete_residue_topologies.csv if every residue matched exactly, else unmatching_residue_topologies.csv. Columns: resname, atoms, rixs, ff_atoms, unmatching_masses (+ possible_atoms in the unmatching case).



interactors.py
Bridges topology_recognition and contacts. Reads the per-residue forcefield matches produced by topology_recognition, reconstructs full topology info (chains, polymer back-bonds, ter chargegroups), deduplicates residues into unique chemical variants, and assigns interaction roles to each variant.
Parameters:

required:

topology (str): same topology file as the rest of the pipeline
general_outpath (str | Path): root output folder. Must already contain a Topology_recognition/ subfolder produced by topology_recognition.


optional:

na_values (List[str]): see main.py
outlog (str | Path): folder for log files. Defaults to general_outpath.



Two public methods:

get_complete_topology()

fills bidirectional bond info per residue type
concatenates per-residue dataframes into one whole-topology table (in memory only)
detects chain boundaries from +-bond connectivity (polymer detection) and from contiguous rtp_file runs (non-polymers)
stamps Chain and Entity (P polymer / M monomer)
adds polymer back-bonds (-X references in the bonds column)
adjusts ter chargegroups for terminal polymer residues
deduplicates residues by their atom-level chemistry (after back-bonds and ter adjustments) into unique variants

Outputs to <general_outpath>/Interactors/:

found_topologies/<variant_id>.csv: one CSV per unique residue chemistry. <variant_id> is <resname> for the first variant of a residue type, then <resname>_v2, _v3, … for subsequent ones (e.g. ALA.csv for internal alanine, ALA_v2.csv for N-terminal, ALA_v3.csv for C-terminal). Columns: rtp_file, resname, atom_name, atom_type, charge, chargegroup, bonds, mass.
residue_index.csv: one row per residue instance. Columns: Resindex, resnum, Chain, Entity, variant_id.
atom_index_map.csv: one row per atom in the trajectory. Columns: Resindex, atom_name, atom_index.


add_inttypes(**kwargs)
Iterates over the per-variant CSVs and appends interaction-role columns in place. Each variant is analysed in isolation (no cross-residue context); polymer back-bonds are still present in the bonds strings but BFS terminates at the variant boundary.

q_interactors (charge-based, PIE-style):

HBonds — HAcc / HDon / HBond based on partial charges and connectivity. Tunable kwargs: max_DA_charge (default −0.3), min_H_charge (0.3), DA_masses ([14, 16, 19] for N/O/F).
Salt bridges — Acidic / Basic from the recursive local charge of an atom + its terminal subtree. Tunable: min_abs_charge (0.75).
Hydrophobic clusters — HC flag plus an HC_cluster group label, restricted to atoms with |charge| ≤ max_abs_charge and BFS distance > min_bonds_to_charged from any charged atom. Tunable: max_abs_charge (0.27), max_num_heavy_ats (3), min_bonds_to_charged (3).


named_interactors (PyInteraph-compatible, name-based):

rules of the form (attribute, [values]) or not <attribute>-prefixed, AND-combined per role
default rule sets cover HDon, HAcc, HC, Basic, Acidic matching the original PyInteraph atom definitions
HBond H atoms filled afterwards by walking from HDon atoms via the bonds column

Each variant CSV is read, analysed, and overwritten with two added columns (q_interactors, named_interactors) plus an optional HC_cluster.
Note on cross-residue context: at variant-level analysis, the polymer 3-residue BFS widening that the previous version applied for HC classification is no longer in effect. Differences are negligible for the standard amino-acid set (charged neighbour atoms rarely sit within min_bonds_to_charged=3 bond hops of a candidate HC atom via the backbone), but worth verifying on a representative trajectory if HC labels look off.
#### `contacts.py` (formerly `adjacency.py`)

Computes contact statistics across an MD trajectory using the per-variant tables produced by `interactors`.

Parameters:
- **required**:
  - `topology` (str): same topology file as the rest of the pipeline
  - `trajectory` (str): MD trajectory file (any format MDAnalysis can read)
  - `general_outpath`: same root folder. Reads from `<general_outpath>/Interactors/`, writes to `<general_outpath>/Contacts/`.
- **optional** (with defaults):
  - `nc` (`2`): number of parallel workers
  - `ram_thresh` (`80`): RAM-usage percentage that triggers batch-size freezing for H-bond computation
  - `nframes` (`500`): initial frame batch size for adaptive H-bond batching
  - `beg_frm` (`0`), `end_frm` (`-1`): frame range to analyse (`-1` = last frame)
  - `interactors_type` (`['named_interactors']`): which interaction-role columns to consult on the variant tables. Allowed: `'q_interactors'`, `'named_interactors'`.
  - `targets` (`['HBonds', 'Salt_bridges', 'HCs']`): contact types to compute. Allowed: `'HBonds'`, `'Salt_bridges'`, `'HCs'`.
  - `subject_chains` (defaults to detected polymer chains): chains analysed as the focus of contacts
  - `surrounding_chains` (defaults to all detected chains): chains considered as candidate partners
  - `d_a_dist` (`4.5` Å): donor–acceptor distance cutoff for H-bonds
  - `angle` (`150` °): donor–H–acceptor angle cutoff for H-bonds
  - `sb_cutoff` (`3.5` Å): salt-bridge atom-pair distance cutoff
  - `hc_cutoff` (`5.5` Å): hydrophobic-cluster centre-of-mass distance cutoff

Workflow:

1. Reassembles a flat per-atom interactors table in memory by joining the three files written by `interactors` (`found_topologies/<variant>.csv`, `residue_index.csv`, `atom_index_map.csv`). The legacy `whole_topology_interactors.csv` is no longer used.
2. Writes a corrected-topology PDB (`<name>_correct.pdb`) with `Resindex` stamped as `resid` and atom IDs aligned with the `atom_index` map. This is the topology used by all subsequent MDAnalysis operations.
3. Splits frames across `nc` worker processes.
4. For each chain combination and each target contact type, dispatches workers. Each worker streams its frame chunk and writes partial CSVs `<name>_upto_<lastframe>.csv`. H-bond workers use an adaptive batch size that grows until RAM hits `ram_thresh`, then freezes at the largest fitting size.
5. Merges partial CSVs per chain combination into `<name>_merged.csv`.
6. Aggregates merged CSVs into per-residue-pair contact-frequency tables: `Hbonds_final.csv`, `SB_final.csv`, `HC_final.csv` (columns: `res_1, res_2, n_frames_contact, mean_dist, chain_1, chain_2, label_1, label_2`).
7. Builds adjacency matrices (`<inttype>_adjacency.dat`, `combined_adjacency.dat`) where entry `[i, j]` is the persistency in % of the contact between residues `i` and `j`.
8. Runs persistency-threshold analysis on every non-empty subset of the requested interaction types: scans `p_min` from 0 to 100 in `step` increments, builds the largest connected component at each cutoff, and writes elbow + degree-distribution CSVs and plots. Used to choose the contact-frequency threshold for downstream PSN construction.

Output layout under `<general_outpath>/Contacts/`:

```
Contacts/
├── <name>_correct.pdb
├── named/                       (or q/, depending on interactors_type)
│   ├── H_bonds/
│   │   ├── *_upto_*.csv         (partial, per-worker)
│   │   └── *_merged.csv
│   ├── Salt_bridges/
│   └── HCs/
├── Hbonds_final.csv
├── SB_final.csv
├── HC_final.csv
├── *_adjacency.dat
├── combined_adjacency.dat
├── adjacency_labels.dat
├── csvs/
│   ├── <subset>_persistency_elbow.csv
│   └── <subset>_degree_distribution.csv
└── plots/
    ├── <subset>_persistency_elbow.png
    └── <subset>_degree_distribution.png
```