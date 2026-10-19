"""Experiment 7 (FSI check), step 1: the full COMPLEX fluid-loading operator of the two-sided spherical cap in the
structural coefficient space of the regime pilot.

Same exact spectral two-sided solver, basis, whitening, L list and acoustic samples as tests/r0_pilot_factor.py
(which produced the pilot's power factor Fx).  Per azimuthal order m (|m| <= 20), with G_m = Psi_m diag(w e^{-i m phi}) B
(N x ndof) and the whitening T_m (r x N):

    Galerkin:  Aw_m c = rho c T_m b,   b = G_m v,  v = normal velocity at the q8 samples,  mu = sum (T_m^T c)_n psi_n e^{i m phi}
    generalised fluid force on the structural dofs (virtual work of -mu n on the shell, both faces; the sign of the
    normal convention cancels):   f = -sum_m G_m^H T_m^T c = -Z v,
    Z_minus = sum_m Q_m^H Zb_m Q_m,   Q_m = T_m G_m,   Zb_m = rho c Aw_m^{-1}            (time convention e^{-i w t})
    power  P = 1/2 Re(v^H Z v) = w^2 U^H Q^H Hb Q U,  Hb_m = 1/2 Herm(Zb_m) = (rho c/4)(Ai + Ai^H) = (rho c/2) Ai Re(Aw) Ai^H
    so the Hermitian part of Z/2 is exactly the pilot's power operator (Fx^H Fx before eigen-clipping).
The structural code uses e^{+i w t} (K_D = (1 + i w beta) K - w^2 M), where the operator is Z_plus = conj(Z_minus) and
the coupled system is (K_D + i w Z_plus) U = F (done in r0_fsi_check.py).

Stored per cell (for L = 4000 and 8000; Richardson x = 2*8000 - 4000 is formed by the user):
    Q (sum r_m x ndof, complex), offsets, m list, Zb_L blocks (rho c Aw^-1), Hb_L blocks (PSD form).
rho c = 1.21 * 343 (the pilot's device: air impedance kept, k = kR/R).

Checks logged: reciprocity (Z = Z^T in the real structural basis), Hermitian part of Zb (direct) vs the PSD form per block,
and the full-matrix 1/2 Herm(Z_x) vs Fx^H Fx (pilot factor, clipped) on ndof x ndof.

Run:  python -u tests/r0_fsi_zop.py <cell> [...]    -> <scratch or out dir>/zop_<cell>.npz, results/r0_revision/fsi/zop_<cell>.log
"""
import os, sys, time
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
import numpy as np
import scipy.io as sio

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from r0_nov_capspec import bessel_terms, whitening, RHO0, C0
from r0_pilot_factor import get_basis, load_samples, CELLS, N, M, R

ROOT = os.path.dirname(HERE)
PIL = os.path.join(ROOT, "results", "r0_revision", "pilot")
OUT = os.path.join(ROOT, "results", "r0_revision", "fsi")
ZDIR = os.environ.get("FSI_ZDIR", OUT)     # big npz files (regenerable); default next to the results
LLIST = (4000, 8000)


def build(cell, log):
    th_deg, kR = CELLS[cell]
    th0 = np.radians(th_deg)
    S = load_samples(th_deg)
    th, ph, w, B = S["th"], S["ph"], S["w"], S["B"]
    basis = get_basis(th_deg)
    _, lam = bessel_terms(kR, max(LLIST))
    npts, ndof = B.shape
    chunks = [slice(s, min(s + 20000, npts)) for s in range(0, npts, 20000)]
    Bc = [B[c].T.tocsr() for c in chunks]
    Qs, offs, ms = [], [0], []
    Zb = {L: [] for L in LLIST}
    Hb = {L: [] for L in LLIST}
    herm_err = []
    rhoc = RHO0 * C0
    for m in range(-M, M + 1):
        E = np.exp(-1j * m * ph) * w
        G = np.zeros((N, ndof), complex)
        for c, Bt in zip(chunks, Bc):
            psi = basis.profile_theta(th[c], m)
            G += (Bt @ (psi * E[None, c]).T).T
        T = whitening(th0, N, m)
        Qs.append(T @ G); offs.append(offs[-1] + T.shape[0]); ms.append(m)
        for L in LLIST:
            F = basis.proj[abs(m)][:N, : L + 1]
            Aw = T @ (R * R * (F * lam[None, : L + 1]) @ F.T) @ T.T
            ReAw = T @ (R * R * (F * lam[: L + 1].real[None, :]) @ F.T) @ T.T
            Ai = np.linalg.inv(Aw)
            Z = rhoc * Ai
            H = 0.5 * rhoc * (Ai @ ReAw @ Ai.conj().T); H = 0.5 * (H + H.conj().T)
            Hd = 0.25 * (Z + Z.conj().T)                       # 1/2 Herm(Z) computed directly
            herm_err.append((m, L, np.linalg.norm(Hd - H) / np.linalg.norm(H), np.linalg.norm(Hd - H) / np.linalg.norm(Z),
                             np.linalg.norm(Z - Z.T) / np.linalg.norm(Z)))
            Zb[L].append(Z); Hb[L].append(H)
    out = dict(Q=np.vstack(Qs), offs=np.array(offs), ms=np.array(ms), theta0_deg=th_deg, kR=kR, rhoc=rhoc, R=R,
               w0=S["w0"], ndof=ndof)
    for L in LLIST:
        out[f"Zb{L}"] = np.array(Zb[L], dtype=object)
        out[f"Hb{L}"] = np.array(Hb[L], dtype=object)
    worst = max(herm_err, key=lambda e: e[2])
    log(f"  per-block 1/2 Herm(rho c Aw^-1) vs PSD form: max rel (to ||H||) {worst[2]:.2e} at m={worst[0]:+d} L={worst[1]}; "
        f"max rel (to ||Z||) {max(e[3] for e in herm_err):.2e}; block symmetry ||Z-Z^T||/||Z|| max {max(e[4] for e in herm_err):.2e}")
    for m, L, e1, e2, e3 in herm_err:
        if L == 8000 and abs(m) in (0, 1, 2, 5, 10, 20) and m >= 0:
            log(f"    m={m:+d} L={L}: rel to ||H|| {e1:.2e}, rel to ||Z|| {e2:.2e}")
    return out


def assemble(d, L="x", part="Z"):
    """Dense ndof x ndof operator Q^H blkdiag(X) Q at L in {4000, 8000, 'x' (Richardson)} (e^{-iwt}).
    part 'Z': X = rho c Aw^-1 (raw); 'H': PSD power blocks; 'Zc': the operator used for FSI,
    X = 2 Hb + (Zb - Zb^H)/2, i.e. reactive (anti-Hermitian) part from Aw^-1 and resistive (Hermitian) part from the
    PSD form, because 1/2 Herm(Aw^-1) computed directly is round-off for |m| >= 5 at low kR (log of r0_fsi_zop)."""
    Q, offs = d["Q"], d["offs"]
    nd = Q.shape[1]
    A = np.zeros((nd, nd), complex)

    def blk(p, i):
        if L == "x":
            return 2 * d[f"{p}b8000"][i] - d[f"{p}b4000"][i]
        return d[f"{p}b{L}"][i]
    for i in range(len(offs) - 1):
        if part == "Zc":
            Zi = blk("Z", i)
            X = 2 * blk("H", i) + 0.5 * (Zi - Zi.conj().T)
        else:
            X = blk(part, i)
        Qi = Q[offs[i]:offs[i + 1]]
        A += Qi.conj().T @ (X @ Qi)
    return A


def load_zop(cell):
    z = np.load(os.path.join(ZDIR, f"zop_{cell}.npz"), allow_pickle=True)
    return {k: z[k] for k in z.files}


def main(cells):
    os.makedirs(OUT, exist_ok=True); os.makedirs(ZDIR, exist_ok=True)
    for cell in cells:
        t0 = time.time()
        logf = open(os.path.join(OUT, f"zop_{cell}.log"), "w")

        def log(s):
            print(s, flush=True); logf.write(s + "\n"); logf.flush()
        th_deg, kR = CELLS[cell]
        log(f"cell {cell}: theta0 = {th_deg} deg, kR = {kR:.6f}, rho c = {RHO0*C0:.3f}, c' = {2*np.pi*100*R/kR:.3f} m/s, "
            f"rho' = {RHO0*C0/(2*np.pi*100*R/kR):.5f} kg/m^3")
        d = build(cell, log)
        np.savez(os.path.join(ZDIR, f"zop_{cell}.npz"), **d)
        # full-matrix checks
        Zx = assemble(d, "x", "Z")
        log(f"  reciprocity ||Zx - Zx^T||_F / ||Zx||_F = {np.linalg.norm(Zx - Zx.T)/np.linalg.norm(Zx):.2e}")
        Hx = assemble(d, "x", "H")
        HdZ = 0.25 * (Zx + Zx.conj().T)
        log(f"  ||1/2 Herm(Zx) - Hx(PSD form)||_F / ||Hx||_F = {np.linalg.norm(HdZ - Hx)/np.linalg.norm(Hx):.2e}; "
            f"||Hx||/||Zx|| = {np.linalg.norm(Hx)/np.linalg.norm(Zx):.2e}")
        del HdZ, Zx
        Fx = sio.loadmat(os.path.join(PIL, f"factor_{cell}.mat"), variable_names=["Fx"])["Fx"]
        FF = Fx.conj().T @ Fx
        del Fx
        log(f"  ||Hx - Fx^H Fx||_F / ||Fx^H Fx||_F = {np.linalg.norm(Hx - FF)/np.linalg.norm(FF):.2e} "
            f"(difference = eigen-clipping of the Richardson blocks)")
        del FF, Hx
        log(f"  runtime {time.time()-t0:.0f} s")
        log("R0_FSI_ZOP_DONE")
        logf.close()


if __name__ == "__main__":
    main(sys.argv[1:])
