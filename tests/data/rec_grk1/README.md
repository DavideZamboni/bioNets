# recoverin–GRK1: PyInteraph 1 reference (2018)

Recoverin (myristoylated, 2 Ca²⁺) with the GRK1 N-terminal peptide, 219
residues, simulated in water as five 200 ns replicas (r1–r5, 20 001 frames each,
10 ps). `rec_grk1_1micros.xtc` (lab copy, 1.3 GB, not shipped) is r1 followed by
r2–r5 without their first frame: 100 001 frames. Its topology is
`rec_grk1_nolip.pdb` (3538 atoms, the same as the 2018 `rec_grk1_start.pdb`).

PyInteraph 1 analysed this trajectory in July–August 2018. Its graphs, copied
here into `pyinteraph_2018/`, are the reference:

| File | Content |
|---|---|
| `hc-graph_rec_grk1_1micros.dat`, `sb-graph_rec_grk1_1micros.dat` | whole microsecond (Aug 2018) |
| `{hb,hc,sb}-graph_rec_grk1_r1.dat`, `hb-graph_rec_grk1_r{2..5}.dat` | per replica (Jul 2018) |

Source: `/home/davide/storage/old/lab_dati/PSNs/complete/LAST/compared_data/ordered/rec/wt_grk1/network/`
(next to the 2018 `.tpr`, index file and replica `.xtc`s in `../top_trj/`).

`make_reference.py` turns them into the files the presets compare against:
`*-graph_1micros.dat` (preset `rec_grk1`, config `config.yaml`) and
`*-graph_r1.dat` (preset `rec_grk1_r1`, config `config_r1.yaml`). It merges the
2018 split GLM back into one residue and, for the microsecond H-bonds (no 2018
file), averages the replicas. Its docstring has the details.

## What the 2018 run used

Nothing about the run survives except its outputs, so the settings were read
back from the per-pair tables written next to the graphs (`*_rec_grk1_*.dat`
without `-graph`, which name every group), then confirmed by reproduction:

| | 2018 PyInteraph | myPSNs default | in `config.yaml` |
|---|---|---|---|
| Salt-bridge groups | Asp, Glu, Lys, Ca²⁺ (`sc.cal.CALp`), C-termini (`ctertn`) | + Arg, His, N-termini; 0.8 e charge gate | 2018 groups, gate off |
| C-termini | `ctertn = C, OT1, OT2` matched the neutral COOH termini | dropped by the gate (net 0 e) | `COOH` rule, gate off |
| HC residues | 8 standard + myristoyl (2018 residue A-2 GLM) | 8 standard | + GLM tail (C1–C14, no O1) |
| HC atoms of terminal residues | all but `CA C O N H H1-3 O1 O2 OXT OT1 OT2`: keeps HA of the NH2 N-terminus, HT2 of COOH C-termini | same (default rule selects by `residue`) | default rule |
| HC cutoff | 5.5 Å (`--hc-co 5.5`) | 5.0 Å | 5.5 Å |
| SB cutoff | 4.5 Å | 4.5 Å | 4.5 Å |
| H-bond atoms | default `hydrogen_bonds.ini` | same | same |

The salt-bridge and hydrophobic groups show up in the 2018 tables as
`ASP_sc.asp.COOn`, `CAL_sc.cal.CALp`, `LEU_ctertn`, `GLM_sidechain`, …; no
`ARG_sc.arg.NCNNp` appears anywhere, though Arg is in PyInteraph's default
`charged_groups.ini`.

## Agreement

Weight difference `100·Σ|myPSNs − PyInteraph| / Σ PyInteraph` (strict lower
triangle), and edges PyInteraph / myPSNs / shared, checked 2026-10-26:

| Preset (frames) | SB | HC | H-bonds |
|---|---|---|---|
| `rec_grk1` (100 001) | 0.13 % (126/135/126) | 0.20 % (322/362/322) | not run (~4–5 h) |
| `rec_grk1_r1` (20 001) | 0.10 % (99/112/99) | 0.13 % (223/244/223) | 4.28 % (812/1045/812) |

Every PyInteraph edge is found. For SB and HC no pair differs by more than
0.05 points and every myPSNs-only edge is ≤ 0.05 %, i.e. the residue is
PyInteraph's 1-decimal rounding. H-bonds are at the detector floor seen on the
PyInteraph test case (3.9 %).

Each setting above is needed. Against these references (1 µs), HC is 53.5 % at
5.0 Å, 0.25 % at 5.5 Å without the terminal-atom rules, and 1.15 % with OT1/OT2
also in the centre of mass. The run in `psns_out/cli_outs/rec_grk1` (default
rules, SB 3.5 Å, HC 5.5 Å) gives SB 26 %.
The r1 H-bond check also exposed a myPSNs bug, fixed on 2026-10-26:
inter-chain H-bonds were dropped from `Hbonds_final.csv` and the matrices
(5.58 % with 50 missing GRK1–recoverin edges before the fix).

`risultati_paragone/rec/*_1micros.dat` (in the lab folder, 219 residues) are
**not** graphs of this trajectory: against them no frame window of
`rec_grk1_1micros.xtc` gets below 34 % (SB) / 53 % (HC), while the shared
edges still correlate (r ≈ 0.9). They are kept out of the presets.
