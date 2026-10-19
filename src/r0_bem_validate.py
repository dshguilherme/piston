"""Validate the two BEM formulations in r0_bem_core (bempp-cl, numba backend).
 (a) closed breathing sphere, Burton-Miller exterior Neumann; sinc-model ratio vs analytic.
 (b) open two-sided pressure-jump screen on a circular disc vs the low-frequency dipole asymptote.
Run:  C:\\Users\\dshgu\\.envs\\bempp\\python.exe src\\r0_bem_validate.py
"""
import os, sys, time, csv
import numpy as np
from scipy.special import spherical_jn

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import bempp_cl.api as bem
bem.enable_console_logging = lambda *a, **k: None
import r0_bem_core as core

RHO0, C0 = core.RHO0, core.C0
OUT_CSV = os.path.join(ROOT, "results", "r0_revision", "bem-validation.csv")
OUT_LOG = os.path.join(ROOT, "paper_jsv", "supplement", "bem-validation.log")
os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
os.makedirs(os.path.dirname(OUT_LOG), exist_ok=True)


class Tee:
    def __init__(self, path):
        self.f = open(path, "w", encoding="utf-8")
        self.o = sys.stdout

    def write(self, s):
        self.o.write(s); self.f.write(s); self.f.flush()

    def flush(self):
        self.o.flush(); self.f.flush()


sys.stdout = Tee(OUT_LOG)
rows = []


def mean_edge(grid):
    V = grid.vertices.T; E = grid.elements.T
    e = np.concatenate([np.linalg.norm(V[E[:, i]] - V[E[:, (i + 1) % 3]], axis=1) for i in range(3)])
    return e.mean()


def p_ext_analytic(x, R=1.0, v0=1.0):
    S = 4 * np.pi * R**2
    return 0.5 * S * RHO0 * C0 * ((x**2) / (1 + x**2)) * abs(v0) ** 2   # Re((x^2+ix)/(1+x^2))


def sinc_analytic(x, R=1.0, v0=1.0):
    k = x / R; om = k * C0
    alpha = RHO0 * om**2 / (4 * np.pi * C0)
    S = 4 * np.pi * R**2
    return alpha * S**2 * spherical_jn(0, x) ** 2 * abs(v0) ** 2


print("bempp-cl", bem.__version__, "| backend: numba (cpu) | precision:", bem.DEFAULT_PRECISION)
print("rho0=%g c0=%g ; convention exp(+ikr)/(4 pi r), e^{-iwt}; dp/dn = i w rho0 v_n" % (RHO0, C0))
print("calibrated signs: dD/dn[p] = %+g*H p ; screen: %+g*H mu = i w rho0 vn ; SIGN_P=%+g"
      % (core.SIGN_H, core.SIGN_SCREEN, core.SIGN_P))

# ----------------------------------------------------------------------------- (a) sphere
# set screen signs here from the derivation: p = D[mu], dD/dn[mu] = SIGN_H*H mu, so SIGN_H*H mu = i w rho vn
core.SIGN_SCREEN = core.SIGN_H
print("\n=== (a) closed breathing sphere, R=1, v0=1 ===")
analytic_ratio = {1.0: 2.83229367309428, 2.0: 2.06705452607952, 4.0: 1.21709}
hdr = "%5s %3s %6s %8s %12s %12s %12s %10s %10s %10s %9s %7s" % (
    "kR", "lev", "n_el", "h", "P_bem", "P_far", "P_sinc_mesh", "ratio", "ratio_ana", "rel_err", "far/bem", "t[s]")
print(hdr)
sphere_pass = {}
for x in (1.0, 2.0, 4.0):
    Pa = p_ext_analytic(x); Ps = sinc_analytic(x)
    ra = Ps / Pa
    for lev in (2, 3, 4, 5):
        g = bem.shapes.regular_sphere(lev)
        k = x
        vn = np.ones(g.number_of_vertices)
        t0 = time.time()
        P, p = core.closed_neumann_bm(g, vn, k)
        t1 = time.time() - t0
        q = 1j * k * C0 * RHO0 * vn
        Pf = core.far_field_power_closed(g, p, q, k)
        Psm = core.sinc_power_mesh(g, vn, k)
        ratio = Psm / P
        err = ratio / ra - 1
        print("%5.1f %3d %6d %8.4f %12.5e %12.5e %12.5e %10.5f %10.5f %+10.2e %9.5f %7.1f" % (
            x, lev, g.number_of_elements, mean_edge(g), P, Pf, Psm, ratio, ra, err, Pf / P, t1))
        rows.append(dict(case="sphere_BM", kR_or_ka=x, level=lev, n_elements=g.number_of_elements,
                         h=mean_edge(g), P_bem=P, P_far=Pf, P_reference=Pa, P_sinc_mesh=Psm,
                         ratio=ratio, rel_error=err, far_over_bem=Pf / P, t_solve_s=t1))
        sphere_pass[(x, lev)] = abs(err)
    print("   analytic: P_ext=%.6e  P_sinc=%.6e  ratio=%.6f (table value %.5f)" % (Pa, Ps, ra, analytic_ratio[x]))

gate = sphere_pass[(1.0, 5)] < 0.01 or sphere_pass[(1.0, 4)] < 0.01
print("\nGATE kR=1 (|ratio/2.8323-1|<1%%): finest level err = %.3e -> %s" % (
    sphere_pass[(1.0, 5)], "PASS" if sphere_pass[(1.0, 5)] < 0.01 else "FAIL"))

# ----------------------------------------------------------------------------- (b) disc
print("\n=== (b) rigid disc a=1, two-sided pressure-jump screen, broadside vn=1 ===")
g0 = core.disc_grid(1.0, 8)
sp_default = bem.function_space(g0, "P", 1)
bnd = core.boundary_vertices(g0)
print("dof check (nr=8): vertices=%d, boundary vertices=%d, interior vertices=%d | "
      "bempp P1 default dofs=%d (include_boundary_dofs=False) ; full P1 dofs=%d" % (
          g0.number_of_vertices, bnd.sum(), (~bnd).sum(), sp_default.global_dof_count,
          core.full_p1(g0)[0].global_dof_count))
print("(bempp default P1 space DOES drop the free-edge vertex dofs (169 = interior count); the solver instead uses the "
      "full space and restricts the matrices itself to interior vertices (same thing, keeps vertex ordering): mu=0 at the free edge.)")
print("disc area check: %.8f vs pi=%.8f" % (np.sum(np.linalg.norm(np.cross(
    g0.vertices.T[g0.elements.T[:, 1]] - g0.vertices.T[g0.elements.T[:, 0]],
    g0.vertices.T[g0.elements.T[:, 2]] - g0.vertices.T[g0.elements.T[:, 0]]), axis=1)) / 2, np.pi))
hdr = "%6s %3s %6s %7s %12s %12s %12s %9s %+9s %9s %7s" % (
    "ka", "nr", "n_el", "h", "P_bem", "P_far", "P_asym", "ratio", "ratio-1", "far/bem", "t[s]")
print(hdr)
t_last = None
for ka in (0.05, 0.1, 0.2):
    for nr in (6, 12, 24, 36, 48):
        g = core.disc_grid(1.0, nr)
        vn = np.ones(g.number_of_vertices)
        t0 = time.time()
        P, mu = core.open_screen_power(g, vn, ka)
        t1 = time.time() - t0
        Pf = core.far_field_power_dl(g, mu, ka)
        Pasym = 8 / (27 * np.pi) * RHO0 * C0 * ka**4 * 1.0**6
        ratio = P / Pasym
        print("%6.2f %3d %6d %7.4f %12.5e %12.5e %12.5e %9.5f %+9.2e %9.5f %7.1f" % (
            ka, nr, g.number_of_elements, mean_edge(g), P, Pf, Pasym, ratio, ratio - 1, Pf / P, t1))
        rows.append(dict(case="disc_screen", kR_or_ka=ka, level=nr, n_elements=g.number_of_elements,
                         h=mean_edge(g), P_bem=P, P_far=Pf, P_reference=Pasym, P_sinc_mesh="",
                         ratio=ratio, rel_error=ratio - 1, far_over_bem=Pf / P, t_solve_s=t1))
        t_last = (g.number_of_elements, t1)

cols = ["case", "kR_or_ka", "level", "n_elements", "h", "P_bem", "P_far", "P_reference", "P_sinc_mesh",
        "ratio", "rel_error", "far_over_bem", "t_solve_s"]
with open(OUT_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader(); w.writerows(rows)
print("\nwrote", OUT_CSV)
print("wrote", OUT_LOG)
