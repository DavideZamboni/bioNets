"""Build the rec_grk1 reference graphs from PyInteraph 1's own 2018 output in
pyinteraph_2018/: {hb,hc,sb}-graph_1micros.dat (whole microsecond, preset
rec_grk1) and {hb,hc,sb}-graph_r1.dat (replica r1, preset rec_grk1_r1).

    python tests/data/rec_grk1/make_reference.py

Two things differ from today's myPSNs run of the same trajectory:

* Residue order. The 2018 topology split the N-myristoylated glycine in two
  residues, A-1 GLY (N HN CA HA1 HA2 C O) and A-2 GLM (the myristoyl), so the
  2018 matrices are 220 × 220. Today's structure has one GLM (219 residues).
  Rows/columns 0 and 1 are merged into row 0 by element-wise max; every other
  index k becomes k − 1. The max is exact here: in every graph at most one of
  the two halves has contacts (GLY only H-bonds, the myristoyl only HC).
* H-bonds. There is no 2018 hb-graph of the whole microsecond, only one per
  200 ns replica. rec_grk1_1micros.xtc is r1 followed by r2–r5 without their
  first frame, so the 1 µs persistency is the replica mean up to 4 of 100005
  frames. (For HC and SB, where both exist, 1 µs graph vs replica mean differ
  by 0.20 % / 0.13 %, the 1-decimal rounding of the replica files.)
"""
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SRC  = HERE / 'pyinteraph_2018'


def merge_split_glm(m: np.ndarray) -> np.ndarray:
    """220-residue 2018 order (GLY, myristoyl, ...) → today's 219 (GLM, ...)."""
    out = np.delete(np.delete(m, 0, 0), 0, 1)       # drop GLY: myristoyl → row 0
    glm = np.maximum(m[0, 1:], m[1, 1:])
    out[0, :] = out[:, 0] = glm
    out[0, 0] = 0.0
    return out


def main() -> None:
    sources = {
        'hc-graph_1micros.dat': [SRC / 'hc-graph_rec_grk1_1micros.dat'],
        'sb-graph_1micros.dat': [SRC / 'sb-graph_rec_grk1_1micros.dat'],
        'hb-graph_1micros.dat': sorted(SRC.glob('hb-graph_rec_grk1_r[1-5].dat')),
        **{f'{t}-graph_r1.dat': [SRC / f'{t}-graph_rec_grk1_r1.dat'] for t in ('hb', 'hc', 'sb')},
    }
    for name, files in sources.items():
        assert files and all(f.is_file() for f in files), files
        mean = np.mean([np.loadtxt(f) for f in files], axis=0)
        np.savetxt(HERE / name, merge_split_glm(mean), fmt='%.2f')


if __name__ == '__main__':
    main()
