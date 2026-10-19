"""Factor-2 (baffle-image) analysis of the sinc model, closed sphere.

Claims tested
 (1) Kernel identity: for ANY surface and ANY source density q, the free-space power of a monopole sheet of strength
     q (volume velocity per unit area) is exactly half the sinc-model power alpha * <q, j0 q>.  Checked independently
     through the far field (angular quadrature of |p_inf|^2), on random points of a non-spherical surface.
 (2) On a sphere, for degree l:  r_l = P_sinc/P_ext = 2 * M_l,  M_l = x^4 j_l(x)^2 |h_l'(x)|^2 = P_monopole/P_ext.
     M_l is checked independently by building the single-layer field (interior j_l, exterior h_l) and computing the
     exterior radial velocity it produces and the exterior power of that velocity.
 (3) Window averages of r_l over x in [x0, x0 + pi] tend to 1 at large x (closed surface: one radiating face),
     while r_l oscillates between ~0 (interior Dirichlet resonances j_l(x)=0) and ~2.
Outputs: results/r0_revision/novelty/sphere_factor2.csv, sphere_factor2.log (stdout), sphere_factor2.png/.pdf
"""
import os, sys, csv
import numpy as np
from scipy import special as sps

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "results", "r0_revision", "novelty")
os.makedirs(OUT, exist_ok=True)
RHO0, C0 = 1.21, 343.0
RC = RHO0 * C0


def h1(l, x, d=False):
    return sps.spherical_jn(l, x, d) + 1j * sps.spherical_yn(l, x, d)


def check_identity(seed=1):
    rng = np.random.default_rng(seed)
    n = 300
    # random points on a wavy ellipsoid-like surface with random weights and complex q
    u = rng.uniform(-1, 1, n); ph = rng.uniform(0, 2 * np.pi, n)
    r = 1 + 0.3 * np.cos(3 * ph) * (1 - u * u)
    X = np.c_[r * np.sqrt(1 - u * u) * np.cos(ph), 0.7 * r * np.sqrt(1 - u * u) * np.sin(ph), 1.3 * r * u]
    w = rng.uniform(0.5, 1.5, n) * 4 * np.pi / n
    q = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    res = []
    for k in (0.3, 1.0, 3.0, 8.0):
        w_ = k * C0
        alpha = RHO0 * w_ ** 2 / (4 * np.pi * C0)
        D = np.linalg.norm(X[:, None] - X[None], axis=2)
        Psinc = alpha * np.real(np.conj(w * q) @ np.sinc(k * D / np.pi) @ (w * q))
        # far field of p = -i w rho0 sum q_i w_i e^{ik|x-y|}/(4 pi |x-y|):  p_inf(s) = -i w rho0/(4pi) sum q w e^{-ik s.y}
        nth = 120; xg, wg = np.polynomial.legendre.leggauss(nth)
        nph = 240; pg = 2 * np.pi * np.arange(nph) / nph
        ct = np.repeat(xg, nph); st = np.sqrt(1 - ct * ct); pp = np.tile(pg, nth)
        S = np.c_[st * np.cos(pp), st * np.sin(pp), ct]
        W = np.repeat(wg, nph) * 2 * np.pi / nph
        pinf = -1j * w_ * RHO0 / (4 * np.pi) * np.exp(-1j * k * S @ X.T) @ (w * q)
        Pmono = np.sum(W * np.abs(pinf) ** 2) / (2 * RC)
        res.append((k, Psinc, Pmono, Psinc / Pmono))
    return res


def M_direct(l, x):
    """Single layer of unit density, degree l, on the unit sphere: p = C j_l(kr) h_l(kR) (r<R), C j_l(kR) h_l(kr) (r>R).
    Jump of dp/dr across r=R equals  +i w rho0 * sigma  (e^{-i w t}: grad p = i w rho0 v; sheet of volume-velocity density sigma = 1).
    Returns P_monopole / P_ext(v = u_out), where u_out = exterior radial velocity produced on r = R+."""
    # dp/dr jump = C k [j_l h_l' - j_l' h_l] = C k * i / x^2     (Wronskian, h = h^(1))
    w_ = 1.0  # cancels
    k = x  # R = 1
    C = (1j * w_ * RHO0 * 1.0) / (k * 1j / x ** 2)         # so that jump = +i w rho0 sigma, sigma = 1
    u_out = C * k * sps.spherical_jn(l, x) * h1(l, x, True) / (1j * w_ * RHO0)
    p_out = C * sps.spherical_jn(l, x) * h1(l, x)
    P_mono = 0.5 * np.real(p_out * np.conj(1.0))            # 1/2 Re p conj(sigma) (interior field carries no power)
    Zl = p_out / u_out                                      # exterior impedance of degree l
    P_ext = 0.5 * np.real(Zl) * np.abs(u_out) ** 2          # same exterior velocity, exact exterior problem
    # exterior Neumann problem with v = 1 (instead of u_out):
    P_ext_unit = 0.5 * np.real(Zl)
    return P_mono / P_ext_unit, P_mono / P_ext


def main():
    log = []
    pr = lambda s: (print(s, flush=True), log.append(s))
    pr("(1) monopole-sheet far-field power vs sinc form on 300 random points of a non-spherical surface")
    for k, Ps, Pm, r in check_identity():
        pr(f"   k={k:4.1f}  P_sinc={Ps:.10e}  P_monopole_farfield={Pm:.10e}  ratio={r:.12f}")
    pr("(2) sphere: r_l = 2 M_l, M_l = x^4 j_l^2 |h_l'|^2; independent single-layer field construction")
    worst = 0; worst_self = 0
    for l in range(0, 9):
        for x in (0.3, 1.0, 3.0, 10.0, 18.318, 30.0):
            Ml = x ** 4 * sps.spherical_jn(l, x) ** 2 * np.abs(h1(l, x, True)) ** 2
            Md, Mself = M_direct(l, x)
            worst = max(worst, abs(Md / Ml - 1)); worst_self = max(worst_self, abs(Mself - 1))
    pr(f"   max |M_direct/M_formula - 1| over l=0..8, x in (0.3..30): {worst:.2e}")
    pr(f"   max |P_mono / P_ext(u_out) - 1| (sheet power = exterior power of the velocity it induces): {worst_self:.2e}")
    pr("   paper Table (kR=1): " + " ".join(f"l={l}:{2*1.0**4*sps.spherical_jn(l,1.0)**2*np.abs(h1(l,1.0,True))**2:.4f}" for l in range(6)))
    # (3) window averages
    xs = np.linspace(0.05, 120, 240001)
    rows = []
    pr("(3) window average of r_l over [x0, x0+pi]")
    for l in range(0, 9):
        r = 2 * xs ** 4 * sps.spherical_jn(l, xs) ** 2 * np.abs(h1(l, xs, True)) ** 2
        for x0 in (1, 3, 10, 18.318 - np.pi / 2, 30, 60, 100):
            m = (xs >= x0) & (xs <= x0 + np.pi)
            rows.append(dict(l=l, x0=x0, mean_r=np.trapezoid(r[m], xs[m]) / np.pi, min_r=r[m].min(), max_r=r[m].max()))
        pr(f"   l={l}: " + " ".join(f"[{d['x0']:.0f}]{d['mean_r']:.3f}({d['min_r']:.2f}-{d['max_r']:.2f})" for d in rows[-7:]))
    with open(os.path.join(OUT, "sphere_factor2.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    with open(os.path.join(OUT, "sphere_factor2.log"), "w") as f:
        f.write("\n".join(log) + "\n")
    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": "#222", "text.color": "#111", "axes.labelcolor": "#111",
                         "xtick.color": "#111", "ytick.color": "#111"})
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.6), facecolor="white")
    xp = np.linspace(0.05, 40, 8000)
    cols = ["#1f4e79", "#c0504d", "#4f8f3a", "#7f5fa0"]
    for c, l in zip(cols, (0, 1, 3, 6)):
        r = 2 * xp ** 4 * sps.spherical_jn(l, xp) ** 2 * np.abs(h1(l, xp, True)) ** 2
        ax[0].plot(xp, r, color=c, lw=1.0, label=f"$\\ell={l}$")
    ax[0].axhline(1, color="#444", ls="--", lw=0.8); ax[0].set_ylim(0, 3.2)
    ax[0].set_xlabel("$kR$"); ax[0].set_ylabel("$P_{\\rm sinc}/P_{\\rm exact}$ (closed sphere)")
    ax[0].legend(frameon=False, ncol=2); ax[0].set_facecolor("white")
    for c, l in zip(cols, (0, 1, 3, 6)):
        r = 2 * xs ** 4 * sps.spherical_jn(l, xs) ** 2 * np.abs(h1(l, xs, True)) ** 2
        nwin = int(np.pi / (xs[1] - xs[0]))
        rm = np.convolve(r, np.ones(nwin) / nwin, mode="valid")
        ax[1].plot(xs[nwin // 2: nwin // 2 + len(rm)], rm, color=c, lw=1.2, label=f"$\\ell={l}$")
    ax[1].axhline(1, color="#444", ls="--", lw=0.8); ax[1].set_xlim(0, 120); ax[1].set_ylim(0, 3.2)
    ax[1].set_xlabel("$kR$"); ax[1].set_ylabel("running mean over $\\Delta kR=\\pi$"); ax[1].set_facecolor("white")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"sphere_factor2.{ext}"), facecolor="white", dpi=200)


if __name__ == "__main__":
    main()

