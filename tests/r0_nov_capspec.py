"""Spectral (semi-analytic) solver for a two-sided, zero-thickness open spherical cap.

Problem (same as the bempp-cl open-screen BEM of the paper): a spherical cap S_c = {r = R, theta < theta0}
in free space, air on both faces, both faces move with the same normal velocity v_n (n = outward radial).
Unknown: pressure jump mu = p(r=R+) - p(r=R-), supported on the cap, mu = 0 at the free edge.

On the full sphere the double-layer potential of mu is diagonal in spherical harmonics:
    p = sum a_lm h_l(kr) Y_lm (r>R),  sum b_lm j_l(kr) Y_lm (r<R)     (h_l = h_l^(1), e^{-i w t})
    continuity of v_r:   a h_l'(x) = b j_l'(x),   dp/dr = i w rho0 v_r
    =>  rho0 c0 v_lm = lambda_l mu_lm,     lambda_l = x^2 j_l'(x) h_l'(x),  x = kR.
The cap problem is the dual-series problem  (Lambda mu)(theta,phi) = rho0 c0 v on the cap, mu = 0 off the cap.
It is solved by Galerkin with basis  psi_{n,m} = sin^{|m|}theta * s * P_n(2 s^2 - 1) * e^{i m phi},
s^2 = (cos theta - cos theta0)/(1 - cos theta0), which has the sqrt(distance) edge behaviour and pole regularity.
All Legendre projections are exact (Gauss rule in s, polynomial integrands).  The only truncation is the
harmonic degree L of the operator sum, whose tail behaves like 1/L; we report L-convergence and a Richardson
extrapolation in 1/L.

Power:  P = 1/2 Re int mu conj(v) dS  (identical definition to the BEM code).
Sinc model on the same cap:  P_sinc = alpha int int conj(v) j0(k|x-y|) v dS dS = rho0 c0 R^2 sum x^2 j_l(x)^2 |V_lm|^2,
V_lm = int_cap v conj(Y_lm) dOmega   (Funk-Hecke; alpha = rho0 w^2 / (4 pi c0)).
Free-space monopole (single-layer, strength v) power = P_sinc / 2  (identity of kernels).
"""
import numpy as np
from scipy import special as sps

RHO0 = 1.21
C0 = 343.0


# ----------------------------------------------------------------------------- spherical Bessel pieces
def bessel_terms(x, L):
    """Return (sinc_l, lam_l) for l = 0..L:
    sinc_l = x^2 j_l(x)^2 (sinc-model modal weight),   lam_l = x^2 j_l'(x) h_l^(1)'(x)  (screen operator eigenvalue).
    Stable for l >> x via log-derivative recurrences and the Wronskian j y' - j' y = 1/x^2."""
    l = np.arange(L + 1)
    j = sps.spherical_jn(l, x)
    jp = sps.spherical_jn(l, x, derivative=True)
    lsw = int(min(L, np.floor(x + 25 + 3 * x ** (1 / 3))))
    lam = np.empty(L + 1, dtype=complex)
    ll = l[: lsw + 1]
    y = sps.spherical_yn(ll, x)
    yp = sps.spherical_yn(ll, x, derivative=True)
    lam[: lsw + 1] = x * x * jp[: lsw + 1] * (jp[: lsw + 1] + 1j * yp)
    if lsw < L:
        # downward ratio r_l = j_l / j_{l-1}
        top = L + 60
        r = x / (2 * top + 3)
        rr = np.empty(L + 2)
        for n in range(top, lsw, -1):
            r = 1.0 / ((2 * n + 1) / x - r)
            if n <= L + 1:
                rr[n] = r
        f = 1.0 / rr[lsw + 1: L + 1] - (l[lsw + 1:] + 1) / x          # j_l'/j_l
        # upward ratio s_l = y_l / y_{l-1}
        s = y[lsw] / y[lsw - 1]
        g = np.empty(L - lsw)
        for i, n in enumerate(range(lsw + 1, L + 1)):
            s = (2 * n - 1) / x - 1.0 / s                                # s_n = y_n/y_{n-1}
            g[i] = 1.0 / s - (n + 1) / x                                 # y_n'/y_n
        jl = j[lsw + 1:]
        lam[lsw + 1:] = x * x * (f * jl) ** 2 + 1j * (f * g / (g - f))
    return x * x * j * j, lam


# ----------------------------------------------------------------------------- normalised assoc. Legendre
def legendre_block(u, m, L):
    """Pbar_l^m(u) for l = m..L (rows), int_{-1}^{1} Pbar^2 du = 1 (no Condon-Shortley sign)."""
    m = abs(m)
    u = np.asarray(u, float)
    st = np.sqrt(np.maximum(0.0, 1 - u * u))
    out = np.zeros((L - m + 1, u.size))
    p = np.full(u.size, 1 / np.sqrt(2.0))
    for i in range(1, m + 1):
        p = np.sqrt((2 * i + 1) / (2.0 * i)) * st * p
    out[0] = p
    if L == m:
        return out
    p1 = np.sqrt(2 * m + 3.0) * u * p
    out[1] = p1
    pm2, pm1 = p, p1
    for l in range(m + 2, L + 1):
        a = np.sqrt((4.0 * l * l - 1) / (l * l - m * m))
        b = np.sqrt(((l - 1.0) ** 2 - m * m) / (4.0 * (l - 1) ** 2 - 1))
        pl = a * (u * pm1 - b * pm2)
        out[l - m] = pl
        pm2, pm1 = pm1, pl
    return out


def legendre_rows(u, m, L):
    """Generator of (l, Pbar_l^m(u)) for l = |m|..L, same normalisation as legendre_block (streaming)."""
    m = abs(m)
    u = np.asarray(u, float)
    st = np.sqrt(np.maximum(0.0, 1 - u * u))
    p = np.full(u.size, 1 / np.sqrt(2.0))
    for i in range(1, m + 1):
        p = np.sqrt((2 * i + 1) / (2.0 * i)) * st * p
    yield m, p
    if L == m:
        return
    p1 = np.sqrt(2 * m + 3.0) * u * p
    yield m + 1, p1
    pm2, pm1 = p, p1
    for l in range(m + 2, L + 1):
        a = np.sqrt((4.0 * l * l - 1) / (l * l - m * m))
        b = np.sqrt(((l - 1.0) ** 2 - m * m) / (4.0 * (l - 1) ** 2 - 1))
        pl = a * (u * pm1 - b * pm2)
        yield l, pl
        pm2, pm1 = pm1, pl


# ----------------------------------------------------------------------------- cap basis
class CapBasis:
    """Galerkin basis on a cap, per azimuthal order m; stores Legendre projections (frequency independent)."""

    def __init__(self, theta0, N, L, mlist):
        self.theta0, self.N, self.L, self.mlist = theta0, N, L, list(mlist)
        self.c0 = np.cos(theta0)
        nq = L + max(abs(m) for m in self.mlist) + N + 8
        xs, ws = np.polynomial.legendre.leggauss(nq)
        self.s = 0.5 * (xs + 1)
        self.ws = 0.5 * ws
        self.u = self.c0 + (1 - self.c0) * self.s ** 2
        self.wu = self.ws * (1 - self.c0) * 2 * self.s           # d(cos theta) weights
        self.proj = {}
        for m in self.mlist:
            if abs(m) in self.proj:
                continue
            pw = np.sqrt(2 * np.pi) * self.profile(self.s, m) * self.wu   # (N, nq)
            F = np.zeros((N, L + 1))
            for l, Pl in legendre_rows(self.u, m, L):
                F[:, l] = pw @ Pl
            self.proj[abs(m)] = F

    def profile(self, s, m):
        u = self.c0 + (1 - self.c0) * s ** 2
        sm = np.sqrt(np.maximum(0, 1 - u * u)) ** abs(m)
        t = 2 * s ** 2 - 1
        return np.array([sm * s * sps.eval_legendre(n, t) for n in range(self.N)])

    def profile_theta(self, theta, m):
        d = (np.cos(theta) - self.c0) / (1 - self.c0)
        s = np.sqrt(np.clip(d, 0, None))
        out = self.profile(s, m)
        out[:, d < 0] = 0
        return out


def whitening(theta0, N, m, nq=600, tol=1e-12):
    """Matrix T (r x N) mapping the raw profiles to an orthonormal, truncated set (Gram eigenvalues > tol*max).
    Galerkin on span(T psi) is well conditioned for every |m|; near-null combinations are dropped."""
    shell = CapBasis.__new__(CapBasis)
    shell.c0, shell.N = np.cos(theta0), N
    xs, ws = np.polynomial.legendre.leggauss(nq)
    s = 0.5 * (xs + 1); wu = 0.5 * ws * (1 - shell.c0) * 2 * s
    P = shell.profile(s, m)
    G = (P * wu) @ P.T
    ev, U = np.linalg.eigh(G)
    keep = ev > tol * ev.max()
    return (U[:, keep] / np.sqrt(ev[keep])).T


def solve_cap(basis, R, k, rhs_by_m, L=None, rho0=RHO0, c0=C0, return_mu=False):
    """rhs_by_m[m] = b_n = int psi_n e^{-i m phi} v dS  (WITHOUT rho0 c0).  Returns P_two_sided."""
    L = basis.L if L is None else L
    x = k * R
    _, lam = bessel_terms(x, L)
    P = 0.0
    coefs = {}
    for m, b in rhs_by_m.items():
        F = basis.proj[abs(m)][:, : L + 1]
        A = R * R * (F * lam[None, :]) @ F.T
        c = np.linalg.solve(A, rho0 * c0 * np.asarray(b, complex))
        coefs[m] = c
        P += 0.5 * np.real(np.vdot(b, c))          # 1/2 Re sum conj(b) c ;  mu = sum c psi
    return (P, coefs) if return_mu else P


def sinc_power_from_V(Vlm_by_m, R, k, rho0=RHO0, c0=C0):
    """Vlm_by_m[m] = array over l = 0..Ls of V_lm = int_cap v conj(Y_lm) dOmega."""
    Ls = max(len(v) for v in Vlm_by_m.values()) - 1
    w, _ = bessel_terms(k * R, Ls)
    return rho0 * c0 * R * R * sum(np.sum(w[: len(v)] * np.abs(v) ** 2) for v in Vlm_by_m.values())


# ----------------------------------------------------------------------------- axisymmetric uniform cap helpers
def uniform_cap(theta0, R, k, N=24, L=3000, Ls=None, basis=None):
    """Uniform v_n = 1 on the cap.  Returns (P_two_sided, P_sinc, basis)."""
    if basis is None:
        basis = CapBasis(theta0, N, L, [0])
    # b_n = int psi_n v dS = R^2 * 2 pi * int psi_n d(cos theta)   (v = 1, m = 0)
    b = R * R * 2 * np.pi * (basis.profile(basis.s, 0) * basis.wu).sum(axis=1)
    P2 = solve_cap(basis, R, k, {0: b})
    Ls = Ls or int(k * R + 60 + 10 / max(theta0, 1e-3))
    Ls = min(Ls, basis.L)
    Pb = legendre_block(basis.u, 0, Ls)
    V = np.sqrt(2 * np.pi) * (Pb * basis.wu).sum(axis=1)       # int_cap conj(Y_l0) dOmega
    Ps = sinc_power_from_V({0: V}, R, k)
    return P2, Ps, basis
