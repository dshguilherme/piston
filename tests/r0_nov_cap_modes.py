"""Analytic open-spherical-cap study: sinc model vs exact two-sided power for modal velocity fields.

Part A  ratio P_sinc / P_2s versus kR for axisymmetric fields on caps of half-angle theta0:
        v = 1 (uniform, non-zero at the edge) and v_j = cos((j + 1/2) pi theta / theta0), j = 0..7 (zero at edge).
Part B  ratio pencil on the subspace span{v_0..v_{J-1}}:  S c = rho E c.  For every field in the span,
        rho_min <= P_sinc / P_2s <= rho_max; any two fields whose sinc powers differ by more than
        kappa = rho_max / rho_min keep their order under the exact model (ranking certificate).
        Also: the field minimising sinc radiation efficiency vs the field minimising exact efficiency
        (what an optimiser driven by the sinc model would exploit).
Part C  Kendall tau between sinc and exact power over random fields drawn from the span, versus kR.
Outputs: results/r0_revision/novelty/cap_modes_ratio.csv, cap_modes_pencil.csv, cap_modes_tau.csv
"""
import os, sys, csv, time
import numpy as np
from scipy import linalg as sla
from scipy.stats import kendalltau

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from r0_nov_capspec import CapBasis, bessel_terms, legendre_block, RHO0, C0

OUT = os.path.join(os.path.dirname(HERE), "results", "r0_revision", "novelty")
os.makedirs(OUT, exist_ok=True)
RC = RHO0 * C0
R = 1.0
J = 8
LMAX = 8000
KR = np.unique(np.r_[0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.5, np.arange(4, 40.01, 1.0)])


def fields(basis, theta0):
    th = np.arccos(basis.u)
    V = [np.ones_like(th)] + [np.cos((j + 0.5) * np.pi * th / theta0) for j in range(J)]
    return np.array(V)                                         # (J+1, nq)


def operators(basis, Vf, kR, Ls):
    """S (sinc) and E (exact two-sided) power matrices on the field coefficients (fields = rows of Vf)."""
    k = kR / R
    x = kR
    B = R * R * 2 * np.pi * (basis.profile(basis.s, 0) * basis.wu) @ Vf.T        # (N, nf)
    Pb = legendre_block(basis.u, 0, Ls)
    Vl = np.sqrt(2 * np.pi) * (Pb * basis.wu) @ Vf.T                             # (Ls+1, nf)
    w, _ = bessel_terms(x, Ls)
    S = RC * R * R * Vl.T @ (w[:, None] * Vl)
    out = {}
    _, lam = bessel_terms(x, LMAX)
    for L in (LMAX // 2, LMAX):
        F = basis.proj[0][:, : L + 1]
        A = R * R * (F * lam[None, : L + 1]) @ F.T
        Ai = np.linalg.inv(A)
        out[L] = np.real(RC / 4 * B.T @ (Ai + Ai.conj().T) @ B)
    E = 2 * out[LMAX] - out[LMAX // 2]                                           # Richardson in 1/L
    G = R * R * 2 * np.pi * (Vf * basis.wu) @ Vf.T                               # int v_i v_j dS
    return S, E, G, out


def main():
    rng = np.random.default_rng(20261016)
    rowsA, rowsB, rowsC = [], [], []
    for th_deg in (15, 30, 60, 90, 120, 150):
        t = time.time()
        theta0 = np.radians(th_deg)
        basis = CapBasis(theta0, 40, LMAX, [0])
        Vf = fields(basis, theta0)
        names = ["uniform"] + [f"cos{j}" for j in range(J)]
        print(f"theta0={th_deg}: basis {time.time()-t:.0f}s", flush=True)
        # random coefficients for Part C (fixed draw per cap), amplitude decaying with order
        Cr = (rng.standard_normal((400, J)) + 1j * rng.standard_normal((400, J))) / (1 + np.arange(J))[None]
        for kR in KR:
            Ls = int(kR + 60 + 12 / theta0)
            S, E, G, out = operators(basis, Vf, kR, Ls)
            Lchg = np.abs(np.diag(out[LMAX]) - np.diag(out[LMAX // 2])) / np.diag(E)
            for i, nm in enumerate(names):
                rowsA.append(dict(theta0_deg=th_deg, kR=kR, field=nm, P_sinc=S[i, i], P_2s=E[i, i],
                                  ratio=S[i, i] / E[i, i], sigma_sinc=S[i, i] / (0.5 * RC * G[i, i]),
                                  sigma_2s=E[i, i] / (0.5 * RC * G[i, i]), L_rel_change=Lchg[i]))
            # Part B on the edge-vanishing span (cos0..cos{J-1}), restricted to the resolved subspace:
            # keep E-eigendirections whose eigenvalue exceeds 100 x ||E_L8000 - E_L4000||_2 (truncation-error proxy)
            s, e, g = S[1:, 1:], E[1:, 1:], G[1:, 1:]
            e = 0.5 * (e + e.T); s = 0.5 * (s + s.T)
            D = (out[LMAX] - out[LMAX // 2])[1:, 1:]
            errn = np.linalg.norm(0.5 * (D + D.T), 2)
            # work in G-orthonormal coordinates so that eigenvalues are efficiencies
            gl = np.linalg.cholesky(g); gi = np.linalg.inv(gl)
            et = gi @ e @ gi.T; st = gi @ s @ gi.T
            ev, U = np.linalg.eigh(et)
            keep = ev > 100 * np.linalg.norm(gi @ (0.5 * (D + D.T)) @ gi.T, 2)
            ev_rel = ev[0] / ev[-1]
            if keep.sum() >= 2:
                Uk = U[:, keep]
                ek, sk = Uk.T @ et @ Uk, Uk.T @ st @ Uk
                rho = sla.eigh(sk, ek, eigvals_only=True)
                eff_s, Xs = np.linalg.eigh(sk); eff_e, Xe = np.linalg.eigh(ek)
                xs, xe = Xs[:, 0], Xe[:, 0]
                rB = dict(rho_min=rho[0], rho_max=rho[-1], kappa=rho[-1] / rho[0],
                          ratio_sinc_optimal=(xs @ sk @ xs) / (xs @ ek @ xs),
                          ratio_exact_optimal=(xe @ sk @ xe) / (xe @ ek @ xe),
                          exact_regret=(xs @ ek @ xs) / eff_e[0],
                          sinc_min_eff=eff_s[0] / (0.5 * RC), exact_min_eff=eff_e[0] / (0.5 * RC))
            else:
                rB = dict(rho_min=np.nan, rho_max=np.nan, kappa=np.nan, ratio_sinc_optimal=np.nan,
                          ratio_exact_optimal=np.nan, exact_regret=np.nan, sinc_min_eff=np.nan, exact_min_eff=np.nan)
            rho = [rB["rho_min"], rB["rho_max"]]
            rowsB.append(dict(theta0_deg=th_deg, kR=kR, kept_dim=int(keep.sum()), E_min_eig_rel=ev_rel,
                              E_trunc_err_norm=errn, **rB))            # Part C
            ps = np.real(np.einsum("ri,ij,rj->r", Cr.conj(), s, Cr))
            pe = np.real(np.einsum("ri,ij,rj->r", Cr.conj(), e, Cr))
            tau = kendalltau(ps, pe).statistic
            # pairs whose sinc ratio exceeds kappa are certified; count violations among uncertified pairs
            rowsC.append(dict(theta0_deg=th_deg, kR=kR, kendall_tau=tau, ratio_min=(ps / pe).min(),
                              ratio_max=(ps / pe).max(), rho_min=rho[0], rho_max=rho[-1]))
        print(f"   done {time.time()-t:.0f}s; kR=18: " +
              " ".join(f"{r['field']}:{r['ratio']:.3f}" for r in rowsA if r['theta0_deg'] == th_deg and r['kR'] == 18.0),
              flush=True)
    for name, rows in (("cap_modes_ratio.csv", rowsA), ("cap_modes_pencil.csv", rowsB), ("cap_modes_tau.csv", rowsC)):
        with open(os.path.join(OUT, name), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()


