"""Experiment 7: export i w s Z_plus on the free dofs (symmetrised) for the FSI-in-the-loop re-optimisation in MATLAB.

Z_plus = conj(Z_minus), Z_minus = r0_fsi_zop.assemble(d, 'x', 'Zc') (Richardson L = 4000/8000), e^{+iwt} convention,
symmetrised (Z + Z^T)/2 (reciprocity holds to ~1e-9 relative, see zop_<cell>.log) so that the coupled matrix is complex
symmetric and the pilot's adjoint (same matrix for state and adjoint) stays valid.
s = 1: the pilot's device (rho' c' = rho0 c0);  s = 'phys': rho = 1.21 with c' (s = c'/c0).

Run: python -u tests/r0_fsi_export_z.py <cell> <s|phys> <outfile.mat>
"""
import os, sys
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
import numpy as np
import scipy.io as sio
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from r0_fsi_zop import load_zop, assemble, OUT
from r0_pilot_factor import CELLS, R
from r0_nov_capspec import C0

cell, sarg, fn = sys.argv[1], sys.argv[2], sys.argv[3]
kR = CELLS[cell][1]
s = (2 * np.pi * 100 * R / kR) / C0 if sarg == "phys" else float(sarg)
c = sio.loadmat(os.path.join(OUT, f"ctx_th{CELLS[cell][0]}.mat"), variable_names=["free", "w"])
fr = c["free"].ravel().astype(np.int64) - 1; w = float(c["w"].ravel()[0])
Z = np.conj(assemble(load_zop(cell), "x", "Zc"))[np.ix_(fr, fr)]
Z = 0.5 * (Z + Z.T)
iwZ = 1j * w * s * Z
sio.savemat(fn, dict(iwZ=iwZ, s=s, cell=cell), do_compression=False)
print(f"{cell}: s = {s:.6f}, ||iwZ||_F = {np.linalg.norm(iwZ):.4e} -> {fn}")
