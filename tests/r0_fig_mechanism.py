"""Mechanism data for the JSV paper figure F22 (`fig_mechanism`): normal velocity, surface power density and
far-field directivity of the uniform design (U) and of the exact-driven design X0 at the paper case
(30 deg cap, R = 10 m, 100 Hz, kR = 18.32, steel in air).

Two models, kept apart:
  * VELOCITY  -- in-vacuo structural response stored with the optimization paths: U (displacement coefficients)
    of the uniform design (run S0, iteration 1) and of the end design of run X0 (iteration 600), from
    results/r0_revision/exact_opt/cap_path_states.mat; v_n = i*omega * B U at the 8 x 8 Gauss points of every
    element, B and the weights from results/r0_revision/exact_opt/cap_samples_q8.mat (the same samples the
    exact operator F_x of the runs was built on). Structural convention e^{+i omega t}.
  * POWER     -- exact two-sided operator: the dual-series Galerkin solution of tests/r0_nov_capspec.py (48 functions
    per azimuthal order |m| <= 20, harmonic degree L = 4000 and 8000, basis projections loaded from
    results/r0_revision/novelty/cap_basis_N48_M20_L8000.npz). Solver convention e^{-i omega t}, so the solver
    receives conj(v_n). It returns the pressure jump mu = p+ - p-. From mu:
      surface power density  I(x) = 1/2 Re( mu conj(v) )   [W/m^2], whose integral is the power;
      far field              p ~ e^{ikr}/(kr) sum (-i)^{l+1} a_lm Y_lm,  a_lm = i x^2 j_l'(x) mu_lm,  x = kR,
                             dP/dOmega = |sum (-i)^{l+1} a_lm Y_lm|^2 / (2 rho0 c0 k^2).
    Checks printed (and asserted): the Richardson power 2 P(8000) - P(4000) reproduces P_exact of
    cap_path_powers.csv; the integral of I and the integral of dP/dOmega both equal P(8000).

Outputs (paper_jsv/prototype/figs/data/):
  f22_vmap_<U|X0>.dat   x y vdb dph : element polygons on the developed polar map (as tests/r0_fig_designs.py),
                        vdb = 20 log10(|v_e| / largest |v_e| of both designs), v_e = area-weighted element mean of v_n,
                        clipped at -40 dB; dph = |arg(v_e / v_apex)| in degrees (0 = in phase with the apex ring,
                        180 = antiphase), v_apex = area-weighted mean of v_n over the apex ring (i = 0)
  f22_dpds.dat          s (m), (1/P) dP/ds (1/m) and the running power (1/P) * (power of the cap up to s,
                        summed over the Gauss rows up to and including that row) for U and X0 at the 320
                        meridian Gauss rows; dP/ds = 2 pi R sin(theta) <I>_phi, from the L = 8000 solution.
                        Columns s qU qX cU cX (each design / its own net power) and cUa cXa (both / the net power
                        of U: cUa = cU, cXa = cX * P_exact(X0) / P_exact(U), so the two share one scale)
  f22_dir.dat           polar angle (deg, 0 = +z = apex direction) and azimuthal mean of dP/dOmega in dB re the
                        maximum of U, for U and X0 (L = 8000)
  f22_numbers.tex       TeX macros: powers, radiating/absorbing parts of the surface power, L-sensitivity,
                        P_X0 / P_U and its inverse (magnification of the X0-only view), peak of the X0 running power
Run:  python tests/r0_fig_mechanism.py      (about 2 min, < 1 GB)
"""
import csv
import os
import sys

import numpy as np
import scipy.io as sio
from scipy import special as sps

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from r0_nov_capspec import CapBasis, bessel_terms, legendre_block, RHO0, C0  # noqa: E402
from r0_fig_designs import element_polygons, NMER, NCIRC  # noqa: E402

ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results", "r0_revision")
OUT = os.path.join(ROOT, "paper_jsv", "prototype", "figs", "data")
R = 10.0
TH0 = np.radians(30.0)
N, M, LS = 48, 20, (4000, 8000)
LFF = 80                      # harmonic degrees kept in the far field (j_l'(kR) is negligible beyond)


def load_fields():
    s = sio.loadmat(os.path.join(RES, "exact_opt", "cap_samples_q8.mat"), variable_names=["x", "w", "B", "w0"])
    x, w, B = s["x"], s["w"].ravel(), s["B"].tocsr()
    w0 = float(np.ravel(s["w0"])[0])
    c = sio.loadmat(os.path.join(RES, "exact_opt", "cap_path_states.mat"), squeeze_me=True,
                    variable_names=["runs", "run", "iter", "U"])
    runs, run, it = list(c["runs"]), c["run"], c["iter"]

    def idx(name, i):
        j = np.where((run == runs.index(name) + 1) & (it == i))[0]
        assert j.size == 1
        return int(j[0])

    cols = {"U": idx("S0", 1), "X0": idx("X0", 600)}
    v = {k: 1j * w0 * (B @ c["U"][:, j]) for k, j in cols.items()}
    return x, w, w0, v, cols


def solve_mu(basis, lam, th, ph, w, v_solver, L):
    """mu at the sample points, mu_lm (l <= LFF) and the power, for one harmonic degree L."""
    mu = np.zeros_like(v_solver)
    mulm, P = {}, 0.0
    for m in range(-M, M + 1):
        psi = basis.profile_theta(th, m)
        bm = psi @ (np.exp(-1j * m * ph) * v_solver * w)
        F = basis.proj[abs(m)][:, : L + 1]
        A = R * R * (F * lam[None, : L + 1]) @ F.T
        cc = np.linalg.solve(A, RHO0 * C0 * bm)
        P += 0.5 * np.real(np.vdot(bm, cc))
        mu += (cc @ psi) * np.exp(1j * m * ph)
        mulm[m] = F[:, : LFF + 1].T @ cc
    return mu, mulm, P


def directivity(mulm, xk, k, thg, phg):
    """dP/dOmega on a (theta, phi) grid from mu_lm."""
    l = np.arange(LFF + 1)
    a_fac = (-1j) ** (l + 1) * 1j * xk ** 2 * sps.spherical_jn(l, xk, derivative=True)
    pf = np.zeros((thg.size, phg.size), complex)
    u = np.cos(thg)
    for m in range(-M, M + 1):
        Pb = legendre_block(u, m, LFF)                       # rows l = |m|..LFF
        coef = a_fac[abs(m):] * mulm[m][abs(m):]
        pf += np.outer(coef @ Pb, np.exp(1j * m * phg)) / np.sqrt(2 * np.pi)
    return np.abs(pf) ** 2 / (2 * RHO0 * C0 * k * k)


def main():
    os.makedirs(OUT, exist_ok=True)
    x, w, w0, v, cols = load_fields()
    k = w0 / C0
    xk = k * R
    r = np.linalg.norm(x, axis=0)
    th = np.arccos(np.clip(x[2] / r, -1, 1))
    ph = np.arctan2(x[1], x[0])
    nq = th.size // (NMER * NCIRC)
    elem = np.repeat(np.arange(NMER * NCIRC), nq)           # 64 consecutive points per element
    ring = np.floor(np.degrees(th) / (np.degrees(TH0) / NMER)).astype(int)
    assert np.all(ring == elem % NMER), "quadrature points are not grouped by element in design order"
    d = np.load(os.path.join(RES, "novelty", f"cap_basis_N{N}_M{M}_L8000.npz"))
    basis = CapBasis(TH0, N, 4, [0])
    basis.proj = {int(key[1:]): d[key] for key in d.files}
    _, lam = bessel_terms(xk, max(LS))
    pw = {}
    with open(os.path.join(RES, "exact_opt", "cap_path_powers.csv"), newline="") as f:
        for row in csv.DictReader(f):
            pw[int(row["idx"])] = float(row["P_exact"])
    polys = element_polygons(30)
    thg = np.radians(np.arange(0, 181, 1.0))
    phg = np.radians(np.arange(0, 360, 2.0))
    res = {}
    for key in ("U", "X0"):
        vs = np.conj(v[key])                                  # e^{+iwt} structure -> e^{-iwt} solver
        sol = {L: solve_mu(basis, lam, th, ph, w, vs, L) for L in LS}
        P = {L: sol[L][2] for L in LS}
        Pext = 2 * P[8000] - P[4000]
        assert abs(Pext / pw[cols[key]] - 1) < 1e-5, "power does not reproduce cap_path_powers.csv"
        mu = sol[8000][0]
        I = 0.5 * np.real(mu * np.conj(vs))
        PI = float(np.sum(I * w))
        assert abs(PI / P[8000] - 1) < 1e-8
        D = {L: directivity(sol[L][1], xk, k, thg, phg) for L in LS}
        # integral of dP/dOmega (trapezoid in theta, rectangle in the periodic phi)
        Dm = {L: D[L].mean(axis=1) for L in LS}
        Pff = 2 * np.pi * np.trapezoid(Dm[8000] * np.sin(thg), thg)
        assert abs(Pff / P[8000] - 1) < 2e-3, f"far-field power {Pff} vs {P[8000]}"
        # element means of v and the apex reference
        ve = np.bincount(elem, w * v[key].real) + 1j * np.bincount(elem, w * v[key].imag)
        ve /= np.bincount(elem, w)
        apex = elem % NMER == 0
        vref = np.sum((w * v[key])[apex]) / np.sum(w[apex])
        # surface power per meridian row
        rows, inv = np.unique(np.round(th, 10), return_inverse=True)
        Imean = np.bincount(inv, w * I) / np.bincount(inv, w)
        dpds = 2 * np.pi * R * np.sin(rows) * Imean / P[8000]
        run_p = np.cumsum(np.bincount(inv, w * I)) / P[8000]        # running power up to each row, / P
        res[key] = dict(P=P, Pext=Pext, I=I, D=D, Dm=Dm, Pff=Pff, ve=ve, vref=vref, rows=rows, dpds=dpds, run_p=run_p)
        pos, neg = np.sum(np.clip(I, 0, None) * w), np.sum(np.clip(I, None, 0) * w)
        res[key].update(pos=pos / P[8000], neg=neg / P[8000])
        print(f"{key:3s}: P(4000) {P[4000]:.6e}  P(8000) {P[8000]:.6e}  Richardson {Pext:.6e} "
              f"(csv {pw[cols[key]]:.6e})  far field {Pff:.6e}  radiating {pos / P[8000]:+.2f} P, "
              f"absorbing {neg / P[8000]:+.2f} P")
    vmax = max(np.max(np.abs(res[key]["ve"])) for key in res)   # common reference: largest |v_e| of both
    for key in ("U", "X0"):
        ve, vref = res[key]["ve"], res[key]["vref"]
        vdb = np.maximum(20 * np.log10(np.abs(ve) / vmax), -40.0)
        dph = np.degrees(np.abs(np.angle(ve / vref)))
        with open(os.path.join(OUT, f"f22_vmap_{key}.dat"), "w") as f:
            f.write("x y vdb dph\n")
            for e in range(polys.shape[0]):
                for px, py in polys[e]:
                    f.write(f"{px:.5f} {py:.5f} {vdb[e]:.3f} {dph[e]:.2f}\n")
        print(f"{key:3s}: |v_e| {20 * np.log10(np.abs(ve).min() / vmax):.1f} to "
              f"{20 * np.log10(np.abs(ve).max() / vmax):.1f} dB re max of both; elements in antiphase "
              f"(>90 deg) {np.mean(dph > 90):.2f}")
    # Common scale: both running powers divided by the net power of U. Each shape comes from L = 8000 (cU, cX end
    # at 1); the factor is the ratio of the Richardson powers P_exact, so cXa ends at P_X0 / P_U = 10^(dL/10).
    ratio = pw[cols["X0"]] / pw[cols["U"]]
    cUa = res["U"]["run_p"]
    cXa = res["X0"]["run_p"] * ratio
    with open(os.path.join(OUT, "f22_dpds.dat"), "w") as f:
        f.write("s qU qX cU cX cUa cXa\n")
        for a, b, c, e, g, h, q in zip(R * res["U"]["rows"], res["U"]["dpds"], res["X0"]["dpds"],
                                       res["U"]["run_p"], res["X0"]["run_p"], cUa, cXa):
            f.write(f"{a:.5f} {b:.5f} {c:.5f} {e:.5f} {g:.5f} {h:.5f} {q:.7f}\n")
    ipk = int(np.argmax(cXa))
    print(f"running power on the scale of U: P_X0/P_U {ratio:.6e} ({10 * np.log10(ratio):.4f} dB), magnification "
          f"1/ratio {1 / ratio:.1f}; X0 peak {cXa[ipk]:.5f} (= {res['X0']['run_p'][ipk]:.3f} P_X0) at s "
          f"{R * res['X0']['rows'][ipk]:.3f} m; min of cXa {cXa.min():.2e}; X0 below U at every row: "
          f"{bool(np.all(cXa < cUa))} (smallest margin cUa - cXa {np.min(cUa - cXa):.2e}); "
          f"P8000 ratio {res['X0']['P'][8000] / res['U']['P'][8000]:.6e}")
    Dref = np.max(res["U"]["Dm"][8000])
    dL = {key: np.max(np.abs(10 * np.log10(res[key]["Dm"][8000] / res[key]["Dm"][4000]))) for key in res}
    with open(os.path.join(OUT, "f22_dir.dat"), "w") as f:
        f.write("theta dU dX\n")
        for t, a, b in zip(np.degrees(thg), res["U"]["Dm"][8000], res["X0"]["Dm"][8000]):
            f.write(f"{t:.1f} {10 * np.log10(a / Dref):.3f} {10 * np.log10(b / Dref):.3f}\n")
    with open(os.path.join(OUT, "f22_numbers.tex"), "w", encoding="ascii") as f:
        f.write("% Generated by tests/r0_fig_mechanism.py. Do not edit by hand.\n")
        f.write("\\def\\fmPU{%.3g}\n\\def\\fmPX{%.3g}\n" % (res["U"]["Pext"], res["X0"]["Pext"]))
        f.write("\\def\\fmdL{%.2f}\n" % (10 * np.log10(res["X0"]["Pext"] / res["U"]["Pext"])))
        for key, mac in (("U", "U"), ("X0", "X")):
            f.write("\\def\\fmpos%s{%.1f}\n\\def\\fmneg%s{%.1f}\n" % (mac, res[key]["pos"], mac, res[key]["neg"]))
        f.write("\\def\\fmLdir{%.2f}\n" % max(dL.values()))
        f.write("\\def\\fmkR{%.2f}\n" % xk)
        # net power of X0 on the scale of U, magnification of the X0-only view, peak of the X0 running power
        f.write("\\def\\fmratio{%.4f}\n\\def\\fmratiopct{%.2f}\n\\def\\fmmag{%.0f}\n" % (ratio, 100 * ratio, 1 / ratio))
        f.write("\\def\\fmpeakX{%.1f}\n\\def\\fmpeakXa{%.4f}\n\\def\\fmspeakX{%.1f}\n"
                % (res["X0"]["run_p"][ipk], cXa[ipk], R * res["X0"]["rows"][ipk]))
    print(f"directivity: largest change L 4000 -> 8000 {max(dL.values()):.3f} dB; "
          f"X0 - U at the maximum of U: {10 * np.log10(res['X0']['Dm'][8000][np.argmax(res['U']['Dm'][8000])] / Dref):.1f} dB")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
