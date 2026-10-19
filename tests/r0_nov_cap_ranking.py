"""Ranking fidelity of the sinc model for the six frozen cap designs, versus kR (fields frozen at 100 Hz).

Needs (from r0_nov_cap_frozen.py):  results/r0_revision/novelty/cap_m8_projections_M20.npz,
                                     results/r0_revision/novelty/cap_basis_N48_M20_L8000.npz
For each kR:  6x6 power matrices S (sinc) and E (exact two-sided spectral) on the design fields,
  - Kendall tau between diag(S) and diag(E) (design ranking),
  - ratio pencil S x = rho E x on the 6-D span: every combination of the designs has rho_min <= P_s/P_e <= rho_max,
  - certificate: a pair (a,b) is guaranteed to keep its order if P_s(a)/P_s(b) > kappa = rho_max/rho_min;
    we count certified pairs and actual order violations,
  - cross-check of the spectral sinc against the independent spline-consistent sinc sweep
    (results/r0_revision/viz/cap_sinc_freq_sweep.csv, produced by another agent; read only).
Outputs: results/r0_revision/novelty/cap_ranking.csv, cap_ranking_sinc_crosscheck.csv
"""
import os, sys, csv
import numpy as np
from scipy import linalg as sla
from scipy.stats import kendalltau
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from r0_nov_capspec import bessel_terms, whitening, RHO0, C0

ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "r0_revision", "novelty")
RC = RHO0 * C0
R = 10.0
M = 20
L2 = 8000


def load():
    d = np.load(os.path.join(OUT, f"cap_m8_projections_M{M}.npz"))
    labels = [str(x) for x in d["labels"]]
    b = {m: d[f"b{m}"] for m in range(-M, M + 1)}
    V = {m: d[f"V{m}"] for m in range(-M, M + 1)}
    p = np.load(os.path.join(OUT, f"cap_basis_N48_M{M}_L{L2}.npz"))
    proj = {int(k[1:]): p[k] for k in p.files}
    return labels, b, V, proj


_TC = {}


def _T(m):
    if m not in _TC:
        _TC[m] = whitening(np.radians(30.0), 48, m)
    return _TC[m]


def matrices(kR, b, V, proj):
    x = kR
    Ls = min(V[0].shape[0] - 1, int(kR + 60))
    w, _ = bessel_terms(x, Ls)
    S = sum(RC * R * R * V[m][: Ls + 1].conj().T @ (w[: Ls + 1, None] * V[m][: Ls + 1]) for m in V)
    _, lam = bessel_terms(x, L2)
    Es = {}
    for L in (L2 // 2, L2):
        E = 0
        for m in b:
            T = _T(abs(m))
            F = T @ proj[abs(m)][:, : L + 1]
            bm = T @ b[m]
            A = R * R * (F * lam[None, : L + 1]) @ F.T
            Ai = np.linalg.inv(A)
            E = E + RC / 4 * bm.conj().T @ (Ai + Ai.conj().T) @ bm
        Es[L] = E
    E = 2 * Es[L2] - Es[L2 // 2]
    D = Es[L2] - Es[L2 // 2]
    S = 0.5 * (S + S.conj().T); E = 0.5 * (E + E.conj().T)
    return S, E, 0.5 * (D + D.conj().T)


def main():
    labels, b, V, proj = load()
    kRs = np.unique(np.round(np.r_[0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, np.arange(4.0, 40.01, 1.0),
                                   18.318324510727656], 9))
    rows = []
    for kR in kRs:
        if V[0].shape[0] - 1 < kR + 40:
            continue
        S, E, D = matrices(kR, b, V, proj)
        ps, pe = np.real(np.diag(S)), np.real(np.diag(E))
        tau = kendalltau(ps, pe).statistic
        # pencil restricted to the resolved E-subspace: eigenvalues > 100 x ||E_L8000 - E_L4000||_2
        sc = 1 / np.sqrt(np.abs(np.real(np.diag(E))))          # unit exact power per design (span unchanged)
        Pn = np.diag(sc)
        S, E, D = Pn @ S @ Pn, Pn @ E @ Pn, Pn @ D @ Pn
        ev, U = np.linalg.eigh(E)
        keep = ev > 100 * np.linalg.norm(D, 2)
        kept = int(keep.sum())
        if kept == len(labels):                                 # certificate only valid on the FULL span
            Uk = U[:, keep]
            rho = sla.eigh(Uk.conj().T @ S @ Uk, Uk.conj().T @ E @ Uk, eigvals_only=True)
            kappa = rho[-1] / rho[0]
        else:
            rho, kappa = [np.nan, np.nan], np.nan
        cert = viol = 0
        npairs = 0
        min_gap_uncert = np.inf
        for i, j in combinations(range(len(labels)), 2):
            npairs += 1
            hi, lo = (i, j) if ps[i] > ps[j] else (j, i)
            gap = ps[hi] / ps[lo]
            if np.isfinite(kappa) and gap > kappa:
                cert += 1
            if pe[hi] < pe[lo]:
                viol += 1
        order_s = [labels[i] for i in np.argsort(-ps)]
        order_e = [labels[i] for i in np.argsort(-pe)]
        rows.append(dict(kR=kR, freq_Hz=kR / R * C0 / (2 * np.pi), kendall_tau=tau, pencil_dim=kept, E_eig_min_unitpower=ev[0], trunc_err_unitpower=np.linalg.norm(D, 2), rho_min=rho[0], rho_max=rho[-1],
                         kappa=kappa, pairs=npairs, certified_pairs=cert, order_violations=viol,
                         ratio_min=(ps / pe).min(), ratio_max=(ps / pe).max(),
                         order_sinc="|".join(order_s), order_exact="|".join(order_e),
                         **{f"ratio_{l}": ps[i] / pe[i] for i, l in enumerate(labels)},
                         **{f"Psinc_{l}": ps[i] for i, l in enumerate(labels)},
                         **{f"P2s_{l}": pe[i] for i, l in enumerate(labels)}))
        print(f"kR {kR:6.2f} tau {tau:+.3f} viol {viol:2d}/{npairs} rho [{rho[0]:.3f},{rho[-1]:.3f}] "
              f"kappa {kappa:.2f} dim {kept} cert {cert:2d}  ratios {(ps/pe).min():.3f}-{(ps/pe).max():.3f}", flush=True)
    with open(os.path.join(OUT, "cap_ranking.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    # cross-check vs independent spline-consistent sinc sweep (other agent's file, read only)
    fn = os.path.join(ROOT, "results", "r0_revision", "viz", "cap_sinc_freq_sweep.csv")
    if os.path.exists(fn):
        ref = list(csv.DictReader(open(fn)))
        out = []
        for kR in sorted({float(r["kR"]) for r in ref}):
            if V[0].shape[0] - 1 < kR + 40:
                continue
            Ls = int(kR + 60)
            w, _ = bessel_terms(kR, Ls)
            ps = np.real(sum(RC * R * R * np.sum(w[:, None] * np.abs(V[m][: Ls + 1]) ** 2, axis=0) for m in V))
            for r in ref:
                if float(r["kR"]) == kR:
                    i = labels.index(r["design"])
                    out.append(dict(kR=kR, design=r["design"], P_consistent_ref=float(r["P_consistent"]),
                                    P_sinc_spectral_m8=ps[i], rel_diff=ps[i] / float(r["P_consistent"]) - 1))
        with open(os.path.join(OUT, "cap_ranking_sinc_crosscheck.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
        rd = np.array([o["rel_diff"] for o in out])
        print(f"sinc cross-check vs consistent sweep: {len(out)} points, |rel diff| median {np.median(np.abs(rd)):.2e},"
              f" max {np.abs(rd).max():.2e}", flush=True)


if __name__ == "__main__":
    main()



