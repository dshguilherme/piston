"""Experiment 7: evaluate the FSI-in-the-loop re-optimised designs (reopt_<cell>_<tag>.mat) with the same machinery as
r0_fsi_check.py: in vacuo, device (s=1) and physical density (s=c'/c0), powers with Fx (pilot), unclipped Richardson H
and H at L=8000 (reference spread).  Writes results/r0_revision/fsi/reopt_eval_<cell>.csv and prints a table.

Run: python -u tests/r0_fsi_reopt_eval.py <cell> <tag> [<tag> ...]
"""
import os, sys, csv
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
import numpy as np
import scipy.io as sio
import scipy.linalg as sla
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from r0_fsi_zop import load_zop, assemble, OUT
from r0_fsi_check import Struct, runfile
from r0_pilot_factor import CELLS, R
from r0_nov_capspec import C0

cell, tags = sys.argv[1], sys.argv[2:]
th, kR = CELLS[cell]
sphys = (2 * np.pi * 100 * R / kR) / C0
st = Struct(th); w = st.w; fr = st.free
Fx = sio.loadmat(os.path.join(OUT, "..", "pilot", f"factor_{cell}.mat"), variable_names=["Fx"])["Fx"]
d = load_zop(cell)
Zp = np.conj(assemble(d, "x", "Zc"))[np.ix_(fr, fr)]
Q, offs = d["Q"], d["offs"]
Hx = [2 * d["Hb8000"][i] - d["Hb4000"][i] for i in range(len(offs) - 1)]
H8 = list(d["Hb8000"])


def pq(bl, U):
    QU = Q @ U
    return w ** 2 * sum(np.real(np.vdot(QU[offs[i]:offs[i+1]], bl[i] @ QU[offs[i]:offs[i+1]])) for i in range(len(bl)))


designs = [("uniform", 0.05 * np.ones(960)), ("XU", sio.loadmat(runfile(cell, "XU"), variable_names=["tbest"])["tbest"].ravel()),
           ("XS", sio.loadmat(runfile(cell, "XS"), variable_names=["tbest"])["tbest"].ravel())]
for tg in tags:
    designs.append((f"reopt_{tg}", sio.loadmat(os.path.join(OUT, f"reopt_{cell}_{tg}.mat"), variable_names=["tbest"])["tbest"].ravel()))
rows = []
ref = {}
for name, t in designs:
    K, Mm = st.KM(t)
    Kf = K[fr][:, fr].toarray(); Mf = Mm[fr][:, fr].toarray()
    KD = (1 + 1j * w * st.beta) * Kf - w ** 2 * Mf
    for var, s in (("vac", 0.0), ("dev", 1.0), ("phys", sphys)):
        Uf = sla.solve(KD + 1j * w * s * Zp, st.F[fr])
        U = np.zeros(st.ndof, complex); U[fr] = Uf
        sc = s if var == "phys" else 1.0
        P = (sc * w ** 2 * np.sum(np.abs(Fx @ U) ** 2), sc * pq(Hx, U), sc * pq(H8, U))
        if name == "uniform":
            ref[var] = P[0]
        rows.append(dict(cell=cell, design=name, variant=var, s=s, P_Fx=P[0], P_Hx_unclipped=P[1], P_H8000=P[2],
                         dB_vs_uniform_same_variant=10 * np.log10(P[0] / ref[var]),
                         dB_vs_uniform_H8000=10 * np.log10(P[2] / ref[var]), vol=np.mean(t) / 0.1))
        print(f"{name:14s} {var:5s} P_Fx {P[0]:.4e}  P_Hx {P[1]:+.4e}  P_H8000 {P[2]:.4e}  dB vs U {rows[-1]['dB_vs_uniform_same_variant']:+.2f}"
              f"  (H8000 {rows[-1]['dB_vs_uniform_H8000']:+.2f})", flush=True)
with open(os.path.join(OUT, f"reopt_eval_{cell}.csv"), "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
print("R0_FSI_REOPT_EVAL_DONE")
