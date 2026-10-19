"""Experiment 10 (robust objective): export the fluid-loading operator i w s Z_plus on the free dofs for MATLAB.

Same operator as tests/r0_fsi_export_z.py (Experiment 7): Z_minus = r0_fsi_zop.assemble(build(cell), 'x', 'Zc')
(Richardson L = 4000/8000, reactive part from rho c Aw^-1, resistive part from the PSD form), Z_plus = conj(Z_minus),
restricted to the free dofs and symmetrised (Z + Z^T)/2 so that the coupled matrix is complex symmetric.
Differences from Experiment 7 (computational only): the zop blocks are built in memory with r0_fsi_zop.build and not
written to disk, and the log goes to results/r0_revision/robust/ (no Experiment-7 file is touched).

Variants written: s = 1 (the pilot's device, rho' c' = rho0 c0) and, for the scaled cells, s = c'/c0 (physical air
density with the effective sound speed c').  At the anchor c' = c0, so the two coincide and only s = 1 is written.

Run: python -u tests/r0_robust_zexport.py <cell> <outdir>   -> <outdir>/iwz_<cell>_dev.mat [, iwz_<cell>_phys.mat]
"""
import os, sys, time
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
import numpy as np
import scipy.io as sio

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from r0_fsi_zop import build, assemble, OUT as FSIOUT
from r0_pilot_factor import CELLS, R
from r0_nov_capspec import C0

ROOT = os.path.dirname(HERE)
ROB = os.path.join(ROOT, "results", "r0_revision", "robust")


def main(cell, outdir):
    os.makedirs(ROB, exist_ok=True); os.makedirs(outdir, exist_ok=True)
    logf = open(os.path.join(ROB, f"zexport_{cell}.log"), "a")

    def log(s):
        print(s, flush=True); logf.write(s + "\n"); logf.flush()
    t0 = time.time()
    th_deg, kR = CELLS[cell]
    log(f"== {time.strftime('%Y-%m-%dT%H:%M:%S')} cell {cell}: theta0 {th_deg}, kR {kR:.6f}")
    d = build(cell, lambda s: None)
    c = sio.loadmat(os.path.join(FSIOUT, f"ctx_th{th_deg}.mat"), variable_names=["free", "w"])
    fr = c["free"].ravel().astype(np.int64) - 1; w = float(c["w"].ravel()[0])
    Z = np.conj(assemble(d, "x", "Zc"))[np.ix_(fr, fr)]
    rec = np.linalg.norm(Z - Z.T) / np.linalg.norm(Z)
    Z = 0.5 * (Z + Z.T)
    cp = 2 * np.pi * 100 * R / kR
    variants = [("dev", 1.0)]
    if abs(cp / C0 - 1) > 1e-6:
        variants.append(("phys", cp / C0))
    for name, s in variants:
        fn = os.path.join(outdir, f"iwz_{cell}_{name}.mat")
        sio.savemat(fn, dict(iwZ=1j * w * s * Z, s=s, cell=cell), do_compression=False)
        log(f"  {name}: s = {s:.6f}, ||iwZ||_F = {np.linalg.norm(1j*w*s*Z):.4e}, reciprocity before symmetrising {rec:.2e} -> {fn}")
    log(f"  runtime {time.time()-t0:.0f} s")
    log("R0_ROBUST_ZEXPORT_DONE")
    logf.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
