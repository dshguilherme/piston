r"""
Two-sided (unbaffled) rigid circular disc, broadside oscillation: spectral Galerkin in the Hankel domain.

DERIVATION (time convention e^{-i w t}, a = disc radius, z normal, disc in z=0)
-------------------------------------------------------------------------------
Unknown: pressure jump mu(r) = p(0+) - p(0-).  A zero-thickness screen carries a double layer, so
p is odd in z.  Order-0 Hankel transform in r: F(kappa) = int f(r) J0(kappa r) r dr,
f(r) = int F(kappa) J0(kappa r) kappa dkappa.  For z>0 each spectral component is
p_H(kappa,z) = A(kappa) exp(-gamma z), gamma = sqrt(kappa^2-k^2) (kappa>k),  gamma = -i sqrt(k^2-kappa^2) (kappa<k)
(outgoing for e^{-i w t}: exp(i kz z), kz = i*gamma).  Oddness: p(z<0) = -A exp(+gamma z), so mu_H = 2A.
Euler: -i w rho v_z = -dp/dz  ->  v_z(0+) = (dp/dz)/(i w rho) = -gamma A/(i w rho), i.e.
        v_H(kappa) = -(gamma / (2 i w rho)) mu_H(kappa)         [sign fixed by the convention above].
Work done by the disc on the fluid: the net fluid force on the disc per area is -mu, so
P = (1/2) Re int_disc mu conj(v0) dS (as stated in the task); positivity of P is checked numerically.
Basis: mu = sum c_n phi_n, phi_n = (1-r^2/a^2)^{1/2} P_n^{(0,1/2)}(1-2r^2/a^2)  [scipy eval_jacobi(n,0,1/2,.) ; the weight (1-r^2)^{1/2} r dr = (1+x)^{1/2}(1-x)^0 dx at x=1-2r^2, so beta=1/2 is the SECOND index; the task text "(1/2,0)" with argument 1-2r^2 is the index-swapped convention and fails the identity for n>0],
 Phi_n(kappa) = int_0^a phi_n J0(kappa r) r dr = a^2 * 2^{1/2} Gamma(n+3/2)/n! * (kappa a)^{-3/2} J_{2n+3/2}(kappa a)
(identity verified numerically below, also for a != 1).
Galerkin (test with phi_m on the disc, Parseval for Hankel pairs): int_disc phi_m v dS = 2 pi int Phi_m v_H kappa dkappa, so
   sum_n A_mn c_n = -2 i w rho v0 b_m,  A_mn = int_0^inf gamma Phi_m Phi_n kappa dkappa,  b_m = Phi_m(0) = int phi_m r dr.
Only b_0 = a^2/3 is non-zero (the Jacobi weight (1-x)^{1/2} makes phi_n, n>0, orthogonal to phi_0 ... checked from the formula).
A = A_ev - i A_rad, with
  A_rad = int_0^k sqrt(k^2-kappa^2) Phi_m Phi_n kappa dkappa  (kappa = k sin(th): k^3 cos^2 th sin th dth, smooth, Gauss-Legendre)
  A_ev  = int_k^inf sqrt(kappa^2-k^2) Phi_m Phi_n kappa dkappa = int_0^inf s^2 Phi_m Phi_n ds  (kappa^2 = k^2+s^2: removes the
          endpoint sqrt, smooth), Gauss-Legendre panels on [0,S] + analytic tail for s>S from the large-kappa
          mean of Phi_m Phi_n = K_m K_n (-1)^{m+n} kappa^-4 / pi  (oscillating remainder neglected, O(S^-2) relative).
Power: P = (1/2) v0 Re(2 pi sum c_n Phi_n(0)) = 2 pi w rho v0^2 b_0^2 Im[(A^-1)_00].
 Positivity check: low freq A^-1 ~ 1/A_ev + i A_rad/A_ev^2, Im>0.
Far-field (spectral/Parseval) power through the two half-spaces: P_far = pi/(2 w rho) int_0^k |mu_H|^2 sqrt(k^2-kappa^2) kappa dkappa
 with mu_H = sum c_n Phi_n evaluated on an independent (denser) quadrature grid.  NOTE: this is the energy balance of the
 same discrete solution; it tests the quadrature/consistency of the surface-power formula, not the model itself.
Nondimensionalisation in code: kappa-hat = kappa a, k-hat = ka:  P = 2 pi rho c v0^2 a^2 (ka)/9 * Im[(Ahat^-1)_00].
Low-frequency check: Phi_0(0)=a^2/3, A_ev -> pi/6 (a=1), A_rad -> k^3/27 => P -> 8/(27 pi) rho c k^4 a^6 v0^2 (analytic).
"""
import os, sys, time, math
import numpy as np
from scipy import special as sp, integrate
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = r"C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course\results\r0_revision\novelty"
os.makedirs(OUT, exist_ok=True)
logf = open(os.path.join(OUT, "disc_two_sided.log"), "w", encoding="utf-8")
def say(*a):
    s = " ".join(str(x) for x in a)
    print(s); sys.stdout.flush(); logf.write(s + "\n"); logf.flush()

RHO, C0, AR, V0 = 1.21, 343.0, 1.0, 1.0

def Phi_hat(kap, N):
    """Nondimensional (a=1) Hankel transform of phi_n, n=0..N-1, at kappa array (len,) -> (len,N)."""
    n = np.arange(N)
    K = np.sqrt(2.0) * np.exp(sp.gammaln(n + 1.5) - sp.gammaln(n + 1.0))
    kap = np.asarray(kap, float)[:, None]
    return K[None, :] * kap ** -1.5 * sp.jv((2 * n + 1.5)[None, :], kap)

def A_rad(kh, N, nth=300):
    x, w = np.polynomial.legendre.leggauss(nth)
    th = (x + 1) * np.pi / 4; wt = w * np.pi / 4
    kap = kh * np.sin(th)
    wq = wt * kh ** 3 * np.cos(th) ** 2 * np.sin(th)
    Ph = Phi_hat(kap, N)
    return (Ph * wq[:, None]).T @ Ph

_X20, _W20 = np.polynomial.legendre.leggauss(20)
def A_ev(kh, N, S, tail=True, width=4.0, group=200):
    npan = int(np.ceil(S / width)); Sx = npan * width
    A = np.zeros((N, N))
    for p0 in range(0, npan, group):
        p = np.arange(p0, min(npan, p0 + group))
        s = ((p[:, None] * width) + (_X20[None, :] + 1) * width / 2).ravel()
        w = np.tile(_W20 * width / 2, len(p)) * s ** 2
        Ph = Phi_hat(np.sqrt(kh ** 2 + s ** 2), N)
        A += (Ph * w[:, None]).T @ Ph
    if tail:
        n = np.arange(N)
        K = np.sqrt(2.0) * np.exp(sp.gammaln(n + 1.5) - sp.gammaln(n + 1.0))
        T = math.atan(kh / Sx) / (2 * kh) + Sx / (2 * (kh ** 2 + Sx ** 2))
        sg = (-1.0) ** n
        A += np.outer(K * sg, K * sg) * T / np.pi
    return A

def solve_all(Aev, Arad, kh):
    """P(N) for N=1..Nmax, plus function to get x vector."""
    Nmax = Aev.shape[0]
    pref = 2 * np.pi * RHO * C0 * V0 ** 2 * AR ** 2 * kh / 9.0
    P = np.zeros(Nmax)
    for N in range(1, Nmax + 1):
        Ab = Aev[:N, :N] - 1j * Arad[:N, :N]
        e = np.zeros(N, complex); e[0] = 1
        P[N - 1] = pref * np.linalg.solve(Ab, e)[0].imag
    return P

def pick_N(P, tol=1e-6):
    rc = np.abs(np.diff(P)) / np.abs(P[1:])       # rc[j] = change from N=j+1 to N=j+2
    for N in range(3, len(P) + 1):
        if rc[N - 2] < tol and rc[N - 3] < tol:
            return N, rc[N - 2]
    return None, rc[-1]

def Sfor(Nmax): return int(min(30000, max(3000, 15 * Nmax ** 2)))

def power_point(kh, Nstart=None, S=None, tail=True):
    Nmax = Nstart or max(8, int(1.2 * kh) + 8)
    while True:
        Sx = S or Sfor(Nmax)
        Aev = A_ev(kh, Nmax, Sx, tail); Arad = A_rad(kh, Nmax)
        P = solve_all(Aev, Arad, kh)
        N, rc = pick_N(P)
        if N is not None or Nmax >= 120:
            return dict(P=P, N=N or Nmax, rc=rc, Aev=Aev, Arad=Arad, S=Sx, ok=N is not None)
        Nmax += 8

def far_power(kh, Aev, Arad, N):
    Ab = Aev[:N, :N] - 1j * Arad[:N, :N]
    e = np.zeros(N, complex); e[0] = 1
    xv = np.linalg.solve(Ab, e)
    Ar2 = A_rad(kh, N, nth=800)                  # independent, denser grid
    pref = 2 * np.pi * RHO * C0 * V0 ** 2 * AR ** 2 * kh / 9.0
    return pref * float(np.real(np.conj(xv) @ Ar2 @ xv))

def p_piston(kh):
    return 0.5 * RHO * C0 * np.pi * AR ** 2 * V0 ** 2 * (1 - sp.j1(2 * kh) / kh)
def p_dip(kh):
    k = kh / AR
    return 8 / (27 * np.pi) * RHO * C0 * k ** 4 * AR ** 6 * V0 ** 2

say("=== two-sided rigid disc, spectral Galerkin (Hankel) ===")
say("rho0=%.3f c0=%.1f a=%.2f v0=%.1f" % (RHO, C0, AR, V0))

# ---------- CHECK 1: Hankel identity ----------
say("\n[CHECK 1] Hankel identity  int_0^a r J0(kap r)(1-r^2/a^2)^(1/2) P_n^(0,1/2)(1-2r^2/a^2) dr "
    "= a^2 2^(1/2) G(n+3/2)/n! (kap a)^(-3/2) J_(2n+3/2)(kap a)")
worst = 0.0
for a in (1.0, 2.0):
    for n in (0, 1, 2, 5):
        for kap in (0.3, 2.0, 7.5):
            f = lambda r: r * sp.j0(kap * r) * np.sqrt(1 - r ** 2 / a ** 2) * sp.eval_jacobi(n, 0.0, 0.5, 1 - 2 * r ** 2 / a ** 2)
            num = integrate.quad(f, 0, a, epsabs=1e-13, epsrel=1e-13, limit=400)[0]
            ex = a ** 2 * Phi_hat(np.array([kap * a]), n + 1)[0, n]
            err = abs(num - ex) / max(abs(ex), 1e-300)
            worst = max(worst, abs(num - ex) / (1e-10 + 1e-9 * abs(ex)))  # pass if <1 : abs err < 1e-10 + 1e-9|ex| (quad abs noise)
            say("  a=%.0f n=%d kappa=%.1f quad=% .12e formula=% .12e relerr=%.1e" % (a, n, kap, num, ex, err))
say("  worst err / (1e-10+1e-9|formula|): %.2e  -> %s" % (worst, "PASS" if worst < 1 else "FAIL"))
b = [integrate.quad(lambda r: r * np.sqrt(1 - r * r) * sp.eval_jacobi(n, 0., .5, 1 - 2 * r * r), 0, 1)[0] for n in range(4)]
say("  b_n = int phi_n r dr (n=0..3) = ", ["%.3e" % x for x in b], " (expect 1/3,0,0,0)")

# ---------- CHECK 5a: cutoff convergence (and N) ----------
say("\n[CHECK 5] convergence in quadrature cutoff S (fixed N = converged N_used), tail on/off")
for kh in (0.1, 1.0, 5.0, 30.0):
    r0 = power_point(kh); N = r0["N"]
    say("  ka=%g  N_used=%d  default S=%d" % (kh, N, r0["S"]))
    ref = None
    for fac in (0.25, 0.5, 1, 2, 4):
        S = int(r0["S"] * fac)
        Aev = A_ev(kh, N, S, True); Arad = A_rad(kh, N)
        Pn = solve_all(Aev, Arad, kh)[-1]
        Aev0 = A_ev(kh, N, S, False)
        P0 = solve_all(Aev0, Arad, kh)[-1]
        if fac == 4: ref = Pn
        say("    S=%7d  P(tail)=%.9e  P(no tail)=%.9e" % (S, Pn, P0))
    for fac in (0.25, 0.5, 1, 2):
        S = int(r0["S"] * fac)
        Pn = solve_all(A_ev(kh, N, S, True), A_rad(kh, N), kh)[-1]
        say("    S=%7d  rel diff vs S*4 (tail on) = %.2e" % (S, abs(Pn - ref) / ref))

say("\n[CHECK 5] convergence in N")
for kh in (0.1, 1.0, 5.0, 10.0, 30.0):
    r0 = power_point(kh); P = r0["P"]
    sel = [n for n in (2, 4, 6, 8, 12, 16, 24, 32, 40, 48) if n <= len(P)]
    say("  ka=%g  N_used=%d  " % (kh, r0["N"]) + "  ".join("N=%d:%.8e" % (n, P[n - 1]) for n in sel))

# ---------- sweep ----------
must = [0.01, 0.05, 0.1, 0.2, 0.3, 1, 2, 3, 5, 10, 20, 30]
kas = np.unique(np.round(np.concatenate([np.logspace(-2, np.log10(30), 80), must]), 12))
rows = []
say("\n[SWEEP] %d points" % len(kas))
say("   ka        P2            P_far         P_piston      Ppist/P2   P2/dip    N   rel_chg  P2/(rho c pi a2 v0^2)  surf/far-1   t[s]")
for kh in kas:
    t0 = time.time()
    r = power_point(kh); N = r["N"]
    P2 = r["P"][N - 1]
    Pf = far_power(kh, r["Aev"], r["Arad"], N)
    Pp = p_piston(kh); Pd = p_dip(kh)
    rows.append((kh, P2, Pf, Pp, Pp / P2, P2 / Pd, N, r["rc"]))
    say("%8.4f  %.6e  %.6e  %.6e  %.6f  %.6f  %3d  %.1e  %.6f  %.2e  %5.1f%s" % (
        kh, P2, Pf, Pp, Pp / P2, P2 / Pd, N, r["rc"], P2 / (RHO * C0 * np.pi * AR ** 2 * V0 ** 2),
        abs(P2 - Pf) / P2, time.time() - t0, "" if r["ok"] else "  NOT CONVERGED"))
R = np.array(rows)
np.savetxt(os.path.join(OUT, "disc_two_sided.csv"), R, delimiter=",", comments="",
           header="ka,P_two_sided,P_far,P_piston,ratio_sinc_over_two_sided,ratio_two_sided_over_dipole_asymptote,N_used,rel_change_last_N",
           fmt=["%.10g", "%.10e", "%.10e", "%.10e", "%.8f", "%.8f", "%d", "%.3e"])

# ---------- summaries ----------
def at(kh): return R[np.argmin(abs(R[:, 0] - kh))]
say("\n[CHECK 2] low-frequency limit P2 / [8/(27 pi) rho c k^4 a^6 v0^2]")
for kh in (0.01, 0.05, 0.1, 0.2):
    say("  ka=%.2f  ratio=%.6f" % (kh, at(kh)[5]))
say("  (earlier bempp BEM, finest mesh: 0.986, 0.989, 0.998 at ka=0.05, 0.1, 0.2)")
say("\n[CHECK 3] surface vs far-field power: max |P2/P_far - 1| over sweep = %.2e (min P2>0: %s)" % (
    np.max(np.abs(R[:, 1] / R[:, 2] - 1)), bool(np.all(R[:, 1] > 0))))
say("\n[CHECK 4] high-frequency behaviour")
say("   ka    P_piston/P2   P2/P_piston   P2/(rho c pi a^2 v0^2)")
for kh in (10, 20, 30):
    q = at(kh); say("  %5.1f   %.6f     %.6f      %.6f" % (q[0], q[4], 1 / q[4], q[1] / (RHO * C0 * np.pi * AR ** 2 * V0 ** 2)))
say("\nTABLE ratio P_piston(sinc)/P_two_sided:")
say("   ka     ratio")
for kh in (0.1, 0.3, 1, 2, 3, 5, 10, 20, 30):
    q = at(kh); say("  %5.1f   %.6f   (P2/dip=%.6f, P2/P_piston=%.6f)" % (q[0], q[4], q[5], 1 / q[4]))
say("\nlargest N_used: %d ; max rel_change_last_N: %.2e" % (R[:, 6].max(), R[:, 7].max()))

# ---------- plot ----------
plt.rcParams.update({"font.size": 11, "text.color": "#222", "axes.labelcolor": "#222",
                     "xtick.color": "#222", "ytick.color": "#222", "axes.edgecolor": "#222"})
fig, ax = plt.subplots(figsize=(7, 4.6), facecolor="white")
ax.set_facecolor("white")
ax.loglog(R[:, 0], R[:, 4], "-", color="#1f4e9c", lw=2, label=r"$P_{\mathrm{piston}}/P_{2}$ (spectral Galerkin)")
ax.axhline(0.5, color="#a02020", ls="--", lw=1.4, label="0.5")
kk = np.logspace(-2, np.log10(30), 200)
ax.loglog(kk, 27 * np.pi ** 2 / (32 * kk ** 2), ":", color="#444", lw=1.8,
          label=r"low-$ka$ asymptote $27\pi^2/(32\,(ka)^2)$")
ax.set_xlabel("ka"); ax.set_ylabel(r"$P_{\mathrm{piston}}/P_{2\text{-sided}}$")
ax.set_ylim(0.2, 1e4); ax.grid(True, which="both", alpha=0.25)
ax.legend(loc="upper right", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "disc_two_sided.png"), dpi=200, facecolor="white")
fig.savefig(os.path.join(OUT, "disc_two_sided.pdf"), facecolor="white")
say("\nwrote csv/log/png/pdf to " + OUT)
logf.close()


