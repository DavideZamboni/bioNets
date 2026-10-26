# PyInteraph 1 test case

The CHARMM example that ships with PyInteraph 1 (`examples/` of the
pyinteraph repository): one 149-residue chain, 151 frames. myPSNs is validated
against PyInteraph's own graphs of this run; `tests/test_pyinteraph_regression.py`
re-runs the whole pipeline on these files and compares.

| Path | Content | Source |
|---|---|---|
| `structures/sim.prot.gro` | topology | `pyinteraph-master/examples/sim.prot.gro` |
| `structures/traj.xtc` | trajectory, 151 frames | `pyinteraph-master/examples/traj.xtc` |
| `structures/sim.prot.A.pdb` | PyInteraph's `-r` reference (same 2354 atoms, chain A); not used by myPSNs | `pyinteraph-master/examples/sim.prot.A.pdb` |
| `pyinteraph/{hb,hc,sb}-graph.dat` | PyInteraph 1 persistency graphs (%) | `pyinteraph-master/examples/` |
| `pyinteraph/macro-IIN.dat` | binary macro-IIN: `filter_graph` union of the three graphs at 20 % | `pyinteraph-master/examples/` |
| `forcefields/charmm36-apr2019_vm.ff/` | CHARMM36 (Apr 2019) for ff-reader | lab copy, pruned |
| `config.yaml` | myPSNs settings of the validation | — |

The structure, trajectory and graphs are byte-identical to the PyInteraph
distribution (md5 `5098bda0…`, `7fcf7ff1…`, `395fd483…`/`83946b62…`/`b16a37b6…`).
PyInteraph's `tutorial.sh` made the graphs with
`pyinteraph -s sim.prot.gro -t traj.xtc -r sim.prot.A.pdb` (`-b --sb-graph`,
`-y --hb-graph`, `-f --hc-graph`, default cutoffs); myPSNs runs on the same
`sim.prot.gro` + `traj.xtc` pair.

The forcefield keeps only what ff-reader reads: every `*.rtp` and `*.tdb` file
(`merged.rtp`, `merged.n.tdb`, `merged.c.tdb`) and `ffnonbonded.itp`. The
other files of the original folder (`ffbonded.itp`, `cmap.itp`, water and ion
`.itp`s, `.hdb`, `.atp`, …) are never opened by the pipeline.

Note PyInteraph's `hb-graph.dat` has 24 non-zero diagonal entries
(intra-residue H-bonds); `mypsns compare` ignores the diagonal.
