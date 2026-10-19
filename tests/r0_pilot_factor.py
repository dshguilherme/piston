"""R0 regime pilot: exact two-sided (rung R4) and spectral-sinc (cross-check of R1) factors per regime cell.

Generalises tests/r0_opt_exact_factor.py (Experiment 5) to any cap half-angle theta0 and any kR.
For structural coefficients U (spline cap, samples exported by tests/r0_pilot_setup.m 'samples'):
    P_exact = w0^2 ||Fx U||^2,     P_sinc_spectral = w0^2 ||Fs U||^2.
Differences from Experiment 5 (mathematically identical, numerically safer at low kR):
    H_L = (rho c / 4)(Aw^-1 + Aw^-H) is evaluated as (rho c / 2) Aw^-1 Re(Aw) Aw^-H, which holds because Aw is
    complex symmetric, and is PSD by construction (Re lambda_l = x^2 j_l'^2 >= 0). Richardson 2 H_8000 - H_4000 and
    the eigen-clipping are as in Experiment 5.
Acoustic medium: rho0 c0 is kept at the air value (1.21 * 343); only k = kR / R is set per cell, i.e. an effective
sound speed c' = w0 R / kR with rho' c' = rho0 c0. All rungs use the same rho0 c0, so ratios, dB and efficiencies are
unaffected by this choice.
Also stores, for the N-convergence check, the exact factor with the first 32 of the 48 edge-basis functions.

Run:  python -u tests/r0_pilot_factor.py <cell> [<cell> ...]     cells: see CELLS
Writes results/r0_revision/pilot/factor_<cell>.mat and factor_<cell>.log
"""
import os, sys, time
import numpy as np
import scipy.io as sio
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from r0_nov_capspec import bessel_terms, legendre_block, legendre_rows, whitening, CapBasis, RHO0, C0

ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "r0_revision", "pilot")
NOV = os.path.join(ROOT, "results", "r0_revision", "novelty")
R = 10.0
N, M, LMAX = 48, 20, 8000
K100R = 2 * np.pi * 100.0 / C0 * R            # 18.3183 (paper anchor)
CELLS = {                                     # name: (theta0 deg, kR)  -- fixed before any result was computed
    "c30_k1": (30, 1.0), "c30_k3": (30, 3.0),
    "c90_k1p5": (90, 1.5), "c90_k7": (90, 7.0),
    "c120_k1p5": (120, 1.5), "c120_k7": (120, 7.0),
    "c30_k18": (30, K100R),
}


def get_basis(th_deg):
    """Legendre projections of the edge basis, |m| <= M, l <= LMAX (frequency independent)."""
    th0 = np.radians(th_deg)
    if th_deg == 30:
        fn = os.path.join(NOV, f"cap_basis_N{N}_M{M}_L{LMAX}.npz")          # Experiment-2 cache (read only)
    else:
        fn = os.path.join(OUT, f"cap_basis_th{th_deg}_N{N}_M{M}_L{LMAX}.npz")
    b = CapBasis(th0, N, 4, [0]); b.L = LMAX
    if os.path.exists(fn):
        d = np.load(fn); b.proj = {int(k[1:]): d[k] for k in d.files}
    else:
        full = CapBasis(th0, N, LMAX, range(0, M + 1)); b.proj = full.proj
        np.savez(fn, **{f"m{m}": F for m, F in b.proj.items()})
    return b


def load_samples(th_deg):
    d = sio.loadmat(os.path.join(OUT, f"samples_th{th_deg}.mat"))
    x = d["x"]; w = d["w"].ravel(); B = sp.csr_matrix(d["B"])
    r = np.linalg.norm(x, axis=0)
    th = np.arccos(np.clip(x[2] / r, -1, 1)); ph = np.arctan2(x[1], x[0])
    return dict(th=th, ph=ph, w=w, B=B, Uu=d["Uu"].ravel(), w0=float(np.ravel(d["w0"])[0]))


def h_psd(F, lam, T):
    A = R * R * (F * lam[None, :]) @ F.T
    Aw = T @ A @ T.T
    ReAw = T @ (R * R * (F * lam.real[None, :]) @ F.T) @ T.T
    Ai = np.linalg.inv(Aw)
    H = 0.5 * RHO0 * C0 * (Ai @ ReAw @ Ai.conj().T)
    return 0.5 * (H + H.conj().T), Aw


def build(cell, log):
    th_deg, kR = CELLS[cell]
    th0 = np.radians(th_deg)
    S = load_samples(th_deg)
    th, ph, w, B = S["th"], S["ph"], S["w"], S["B"]
    assert th.max() <= th0 + 1e-9
    basis = get_basis(th_deg)
    x = kR
    _, lam = bessel_terms(x, LMAX)
    LS = 110 if kR > 10 else int(np.ceil(x)) + 40
    wsinc, _ = bessel_terms(x, LS)
    u = np.cos(th)
    npts, ndof = B.shape
    chunks = [slice(s, min(s + 20000, npts)) for s in range(0, npts, 20000)]
    Bc = [B[c].T.tocsr() for c in chunks]
    Fx, Fx8, Fx32, Fs = [], [], [], []
    negs, conds = [], []
    Uu = S["Uu"]; w0 = S["w0"]
    Pdirect = 0.0
    for m in range(-M, M + 1):
        E = np.exp(-1j * m * ph) * w
        G = np.zeros((N, ndof), complex)
        for c, Bt in zip(chunks, Bc):
            psi = basis.profile_theta(th[c], m)
            G += (Bt @ (psi * E[None, c]).T).T
        for NN, dest in ((N, None), (32, Fx32)):
            T = whitening(th0, NN, m)
            H = {}
            for L in (4000, 8000):
                F = basis.proj[abs(m)][:NN, : L + 1]
                H[L], Aw = h_psd(F, lam[: L + 1], T)
                if NN == N and L == 8000:
                    conds.append((m, T.shape[0], np.linalg.cond(Aw)))
                    # direct dual-series solve for the uniform design (independent of the eigen-factorisation)
                    b = 1j * w0 * (G @ Uu)
                    c_ = np.linalg.solve(Aw, RHO0 * C0 * (T @ b))
                    Pdirect_m = 0.5 * np.real(np.vdot(T @ b, c_))
            Hx = 2 * H[8000] - H[4000]
            targets = ((Hx, Fx, "x"), (H[8000], Fx8, "8000")) if NN == N else ((Hx, Fx32, "x32"),)
            for HH, lst, tag in targets:
                ev, V = np.linalg.eigh(0.5 * (HH + HH.conj().T))
                if tag == "x":
                    neg = ev[ev < 0]
                    negs.append((m, ev.max(), neg.min() if neg.size else 0.0, neg.size))
                ev = np.clip(ev, 0, None)
                lst.append((np.sqrt(ev)[:, None] * V.conj().T) @ (T @ G[:NN]))
            if NN == N:
                Pdirect += Pdirect_m
        if abs(m) <= LS:
            Vrow = np.zeros((LS - abs(m) + 1, ndof), complex)
            for c, Bt in zip(chunks, Bc):
                Pb = legendre_block(u[c], m, LS)
                Vrow += (Bt @ (Pb * (E[c] / (R * R * np.sqrt(2 * np.pi)))[None, :]).T).T
            Fs.append(np.sqrt(RHO0 * C0 * R * R * wsinc[abs(m):])[:, None] * Vrow)
    out = dict(Fx=np.vstack(Fx), Fx8000=np.vstack(Fx8).astype(np.complex64), Fx32=np.vstack(Fx32).astype(np.complex64),
               Fs=np.vstack(Fs).astype(np.complex64))
    for m, emax, nmin, nn in negs:
        if nn:
            log(f"  m={m:+d}: H eig max {emax:.3e}, {nn} negative, most negative {nmin:.2e} ({nmin/emax:.1e} of max) -> clipped")
    log(f"  max cond(Aw) over m: {max(c for _, _, c in conds):.3e}; kept ranks {[r for _, r, _ in conds]}")
    Pu = {k: w0 ** 2 * np.sum(np.abs(v @ Uu) ** 2) for k, v in out.items()}
    Perp = 0.5 * RHO0 * C0 * w0 ** 2 * np.sum(w * np.abs(B @ Uu) ** 2)
    # Pdirect uses Re<mu,v> with mu from the L=8000 operator: compare with the L=8000 factor
    log(f"  uniform design: P_exact {Pu['Fx']:.6e}  P_exact(L8000) {Pu['Fx8000']:.6e}  direct L8000 solve {Pdirect:.6e} "
        f"(factor/direct {Pu['Fx8000']/Pdirect:.8f})  N32/N48 {Pu['Fx32']/Pu['Fx']:.6f}  "
        f"P_sinc_spec {Pu['Fs']:.6e}  sinc/exact {Pu['Fs']/Pu['Fx']:.4f}  P_ERP {Perp:.6e}  sigma_exact {Pu['Fx']/Perp:.4e}")
    out.update(dict(theta0_deg=th_deg, kR=kR, N=N, M=M, LS=LS, Llist=np.array([4000, 8000]), w0=w0,
                    rho0c0=RHO0 * C0, Pu_exact=Pu["Fx"], Pu_sinc_spec=Pu["Fs"], Pu_erp=Perp, Pu_direct8000=Pdirect))
    return out


def anchor_gate(log):
    """Experiment-5 exact factor (Ai + Ai^H form) vs this script's PSD form on the six frozen designs, theta0 = 30."""
    d5 = sio.loadmat(os.path.join(ROOT, "results", "r0_revision", "exact_opt", "cap_samples_q8.mat"), variable_names=["U"])
    e5 = sio.loadmat(os.path.join(ROOT, "results", "r0_revision", "exact_opt", "cap_exact_factor_q8.mat"),
                     variable_names=["Pex_frozen"])
    f = sio.loadmat(os.path.join(OUT, "factor_c30_k18.mat"), variable_names=["Fx", "w0"])
    U = d5["U"]; w0 = float(np.ravel(f["w0"])[0])
    P = w0 ** 2 * np.sum(np.abs(f["Fx"] @ U) ** 2, 0)
    log("anchor gate (c30_k18, six frozen designs): P_exact(PSD form) / P_exact(Experiment 5 factor) = "
        + ", ".join(f"{a:.7f}" for a in P / e5["Pex_frozen"].ravel()))


def main(cells):
    for cell in cells:
        if cell == "gate":
            logf = open(os.path.join(OUT, "factor_anchor_gate.log"), "w")
            anchor_gate(lambda s: (print(s, flush=True), logf.write(s + "\n")))
            logf.close(); continue
        t0 = time.time()
        logf = open(os.path.join(OUT, f"factor_{cell}.log"), "w")

        def log(s):
            print(s, flush=True); logf.write(s + "\n"); logf.flush()
        log(f"cell {cell}: theta0 = {CELLS[cell][0]} deg, kR = {CELLS[cell][1]:.6f}")
        out = build(cell, log)
        sio.savemat(os.path.join(OUT, f"factor_{cell}.mat"), out, do_compression=False)
        log(f"Fx {out['Fx'].shape}, Fs {out['Fs'].shape}; runtime {time.time()-t0:.0f} s")
        log("R0_PILOT_FACTOR_DONE")
        logf.close()
        del out


if __name__ == "__main__":
    main(sys.argv[1:])
