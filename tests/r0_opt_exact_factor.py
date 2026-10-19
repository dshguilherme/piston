"""Reduce the exact spectral two-sided cap operator to a Hermitian PSD form on the STRUCTURAL coefficients.

For a structural displacement coefficient vector U (spline space of the paper's 40x24 cap, 3024 dofs) the
normal velocity is v = i w0 (B U)(x_q) at acoustic quadrature points x_q (exported by r0_opt_export_samples.m).
The Galerkin right-hand sides of the dual-series solver (r0_nov_capspec.py) are linear in v:
    b_m = int psi_n e^{-i m phi} v dS  = i w0 G_m U,   G_m = Psi_m diag(w) B      (N x ndof)
and the two-sided power is P = 1/2 Re rho c b^H A_m^{-1} b (A_m complex symmetric), i.e.
    P = w0^2 sum_m U^H G_m^T T^T H_m T G_m U,   H_m = (rho c / 4) (Aw^{-1} + Aw^{-H}),  Aw = T A_m T^T
(T = whitening of the edge basis, as in the solver).  H_m is Hermitian PSD (power >= 0); we Richardson-extrapolate
H in 1/L (L = 4000, 8000, the same rule as Experiment 2), eigendecompose, clip the (tiny) negative eigenvalues and
write the factor Fx with  P_exact = w0^2 ||Fx U||^2  -- exactly the 'consistent' R0 interface of articleB_r0_eval.m,
so the paper's GCMMA optimiser and its adjoint run unchanged with the exact operator.

Also builds the spectral (Funk-Hecke) sinc factor Fs (P_sinc = w0^2 ||Fs U||^2) as an independent cross-check of the
paper's consistent sinc operator.

Run:  python -u tests/r0_opt_exact_factor.py [nq]      (default nq = 8; nq = 12 for the quadrature check)
Writes results/r0_revision/exact_opt/cap_exact_factor_q<nq>.mat and cap_exact_factor_q<nq>.log
"""
import os, sys, time
import numpy as np
import scipy.io as sio
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from r0_nov_capspec import bessel_terms, legendre_block, whitening, RHO0, C0
from r0_nov_cap_frozen import get_basis, R, TH0, K100

ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "r0_revision", "exact_opt")
N, M, LMAX, LS = 48, 20, 8000, 110


def load_samples(nq):
    d = sio.loadmat(os.path.join(OUT, f"cap_samples_q{nq}.mat"))
    x = d["x"]; w = d["w"].ravel(); B = sp.csr_matrix(d["B"])
    r = np.linalg.norm(x, axis=0)
    th = np.arccos(np.clip(x[2] / r, -1, 1)); ph = np.arctan2(x[1], x[0])
    labels = [str(np.ravel(v)[0]) for v in d["labels"].ravel()]
    return dict(th=th, ph=ph, w=w, B=B, U=d["U"], designs=d["designs"], labels=labels, w0=float(np.ravel(d["w0"])[0]))


def build_factors(S, k=K100, Llist=(4000, 8000), log=print, with_sinc=True):
    th, ph, w, B = S["th"], S["ph"], S["w"], S["B"]
    assert th.max() <= TH0 + 1e-9, "samples outside the cap"
    basis = get_basis(N, M, LMAX)
    x = k * R
    _, lam = bessel_terms(x, max(Llist))
    Fx, Fx8, Fs, info = [], [], [], []
    wsinc, _ = bessel_terms(x, LS)
    u = np.cos(th)
    npts, ndof = B.shape
    chunks = [slice(s, min(s + 20000, npts)) for s in range(0, npts, 20000)]
    Bc = [B[c].T.tocsr() for c in chunks]                              # (ndof, nc) per chunk
    for m in range(-M, M + 1):
        E = np.exp(-1j * m * ph) * w                                   # (np,)
        G = np.zeros((N, ndof), complex)
        for c, Bt in zip(chunks, Bc):
            psi = basis.profile_theta(th[c], m)                        # (N, nc)
            G += (Bt @ (psi * E[None, c]).T).T                         # (N, ndof)
        T = whitening(TH0, N, m)                                       # (r, N)
        H = {}
        for L in Llist:
            F = basis.proj[abs(m)][:, : L + 1]
            A = R * R * (F * lam[None, : L + 1]) @ F.T
            Aw = T @ A @ T.T
            Ai = np.linalg.inv(Aw)
            H[L] = 0.25 * RHO0 * C0 * (Ai + Ai.conj().T)
        Hx = 2 * H[8000] - H[4000]
        for HH, dest, tag in ((Hx, Fx, "x"), (H[8000], Fx8, "8000")):
            ev, V = np.linalg.eigh(0.5 * (HH + HH.conj().T))
            neg = ev[ev < 0]
            if tag == "x":
                info.append((m, T.shape[0], ev.max(), neg.min() if neg.size else 0.0, neg.size))
            ev = np.clip(ev, 0, None)
            dest.append((np.sqrt(ev)[:, None] * V.conj().T) @ (T @ G))
        if with_sinc and abs(m) <= LS:
            Vrow = np.zeros((LS - abs(m) + 1, ndof), complex)
            for c, Bt in zip(chunks, Bc):
                Pb = legendre_block(u[c], m, LS)                       # rows l = |m|..LS
                Vrow += (Bt @ (Pb * (E[c] / (R * R * np.sqrt(2 * np.pi)))[None, :]).T).T
            Fs.append(np.sqrt(RHO0 * C0 * R * R * wsinc[abs(m):])[:, None] * Vrow)
    for m, r, emax, nmin, nn in info:
        if nn:
            log(f"  m={m:+d}: kept {r} of {N}; H eig max {emax:.3e}, {nn} negative, most negative {nmin:.2e} "
                f"({nmin/emax:.1e} of max) -> clipped")
    out = dict(Fx=np.vstack(Fx), Fx8000=np.vstack(Fx8))
    if with_sinc:
        out["Fs"] = np.vstack(Fs)
    return out


def main(nq=8):
    t0 = time.time()
    logf = open(os.path.join(OUT, f"cap_exact_factor_q{nq}.log"), "w")

    def log(s):
        print(s, flush=True); logf.write(s + "\n"); logf.flush()
    S = load_samples(nq)
    log(f"nq={nq}: {S['th'].size} points, ndof {S['B'].shape[1]}, theta max {np.degrees(S['th'].max()):.6f} deg, "
        f"w0 = {S['w0']:.9g}, kR = {K100*R:.6f}")
    Fa = build_factors(S, log=log)
    w0 = S["w0"]; U = S["U"]
    Px = w0 ** 2 * np.sum(np.abs(Fa["Fx"] @ U) ** 2, 0)
    P8 = w0 ** 2 * np.sum(np.abs(Fa["Fx8000"] @ U) ** 2, 0)
    Ps = w0 ** 2 * np.sum(np.abs(Fa["Fs"] @ U) ** 2, 0)
    import csv
    ref = {r["design"]: r for r in csv.DictReader(open(os.path.join(ROOT, "results", "r0_revision", "novelty",
                                                                      "cap_frozen_validation_M20.csv"))) if r["level"] == "8"}
    log("frozen designs: P_exact(factor) / P2s_extrap(m8 P1 trace, Exp.2) | Richardson rel. corr | "
        "P_sinc_spec(factor) / P_sinc_consistent(paper)")
    for i, lab in enumerate(S["labels"]):
        r = ref[lab]
        log(f"  {lab:13s} Pex {Px[i]:.6e}  /trace {Px[i]/float(r['P2s_extrap']):.6f}  (L8000-only rel {P8[i]/Px[i]-1:+.2e})"
            f"   Psinc {Ps[i]:.6e} /consistent {Ps[i]/float(r['P_sinc_consistent']):.6f}   ratio sinc/exact {Ps[i]/Px[i]:.4f}")
    # azimuthal content of the exact power for the frozen designs (truncation check)
    rows_per_m = Fa["Fx"].shape[0] // (2 * M + 1) if Fa["Fx"].shape[0] % (2 * M + 1) == 0 else None
    sio.savemat(os.path.join(OUT, f"cap_exact_factor_q{nq}.mat"),
                dict(Fx=Fa["Fx"], Fx8000=Fa["Fx8000"], nq=nq, N=N, M=M, Llist=np.array([4000, 8000]), kR=K100 * R,
                     Pex_frozen=Px, Psinc_spec_frozen=Ps, labels=np.array(S["labels"], dtype=object)),
                do_compression=False)
    log(f"Fx {Fa['Fx'].shape}, Fs {Fa['Fs'].shape}; runtime {time.time()-t0:.0f} s")
    log("R0_OPT_EXACT_FACTOR_DONE")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 8)
