"""Verification of the complex (reactive + resistive) operator of the exact spherical-cap solver r0_nov_capspec.

TEST A  near-closed sphere modal impedance vs closed-sphere rho0 c0 / lam_l  (m = 0, v = P_l(cos theta))
TEST B  flat-disc limit: rigid disc broadside added mass (8/3) rho0 a^3 (small cap, uniform v = 1)
TEST C  sign / convention sanity (Im Zc < 0 mass-like, e^{-i w t})
Outputs: results/r0_revision/fsi/verify.log, verify.csv
"""
import os
os.environ['OMP_NUM_THREADS'] = '4'
os.environ['MKL_NUM_THREADS'] = '4'
os.environ['OPENBLAS_NUM_THREADS'] = '4'
import sys, time, csv
import numpy as np
from scipy import special as sps
from scipy.special import roots_legendre

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
# np.polynomial.legendre.leggauss builds a dense n x n matrix (GBs for n ~ 16000): use scipy's O(n) rule instead
np.polynomial.legendre.leggauss = lambda n: roots_legendre(int(n))
import r0_nov_capspec as cs

OUT = os.path.normpath(os.path.join(HERE, '..', 'results', 'r0_revision', 'fsi'))
os.makedirs(OUT, exist_ok=True)
logf = open(os.path.join(OUT, 'verify.log'), 'w', encoding='utf-8')
rows = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    logf.write(s + '\n')
    logf.flush()


RHO0, C0 = cs.RHO0, cs.C0
T0 = time.time()
checks = []   # (name, value, tol, passed)


def check(name, val, tol):
    ok = bool(abs(val) < tol)
    checks.append((name, val, tol, ok))
    return ok


def zc_whitened(basis, T, R, k, b, L, lam):
    """Zc = int mu conj(v) dS at truncation L using whitened basis."""
    F = basis.proj[0][:, :L + 1]
    A = R * R * (F * lam[None, :L + 1]) @ F.T
    Aw = T @ A @ T.T
    Tb = T @ np.asarray(b, float)
    cw = np.linalg.solve(Aw, RHO0 * C0 * Tb)
    return np.vdot(Tb, cw)


def richardson(z1, z2):
    return 2 * z2 - z1


# ============================================================================ TEST A
log('=' * 100)
log('TEST A: near-closed cap vs closed sphere modal impedance (m=0)')
log('  zcap = Zc / int|v|^2 dS, zsph = rho0 c0 / lam_l, lam_l = x^2 j_l\'(x) h_l^(1)\'(x); Zc via Richardson in 1/L')
log('=' * 100)
R = 1.0
NA = 32
L1A, L2A = 2000, 4000
holes = [10, 5, 2]
kRs = [1.5, 1.0]
ells = [0, 1, 2, 3]

# interior-resonance check
for x in kRs:
    jp = [sps.spherical_jn(l, x, derivative=True) for l in ells]
    log(f'  kR={x}: j_l\'(kR) for l=0..3 = ' + ', '.join(f'{v:+.4f}' for v in jp) +
        '  (interior Neumann resonances at j_l\'=0: l=1 x=2.08, l=2 x=3.34, l=3 x=4.51, l=0 x=4.49 -> none near)')

resA = {}   # (kR, l, hole) -> (zcap_x, ratio_x)
for hole in holes:
    th0 = np.deg2rad(180 - hole)
    t = time.time()
    basis = cs.CapBasis(th0, NA, L2A, [0])
    Tm = cs.whitening(th0, NA, 0)
    log(f'\n  hole half-angle {hole} deg (theta0={180-hole}): basis built in {time.time()-t:.1f}s, whitened rank {Tm.shape[0]}/{NA}')
    for kR in kRs:
        k = kR / R
        _, lam = cs.bessel_terms(kR, L2A)
        _, lam_s = cs.bessel_terms(kR, 3)
        for l in ells:
            Pl = sps.eval_legendre(l, basis.u)
            b = R * R * 2 * np.pi * (basis.profile(basis.s, 0) * Pl * basis.wu).sum(axis=1)
            nrm = R * R * 2 * np.pi * np.sum(Pl ** 2 * basis.wu)
            z1 = zc_whitened(basis, Tm, R, k, b, L1A, lam)
            z2 = zc_whitened(basis, Tm, R, k, b, L2A, lam)
            zx = richardson(z1, z2)
            zsph = RHO0 * C0 / lam_s[l]
            zcap = zx / nrm
            ratio = zcap / zsph
            resA[(kR, l, hole)] = (zcap, zsph, ratio, z1 / nrm, z2 / nrm)
            rows.append(dict(test='A', kR=kR, l=l, hole_deg=hole, L1=L1A, L2=L2A,
                             Zc_L1_re=(z1 / nrm).real, Zc_L1_im=(z1 / nrm).imag,
                             Zc_L2_re=(z2 / nrm).real, Zc_L2_im=(z2 / nrm).imag,
                             zcap_re=zcap.real, zcap_im=zcap.imag, zsph_re=zsph.real, zsph_im=zsph.imag,
                             ratio_re=ratio.real, ratio_im=ratio.imag, ratio_abs=abs(ratio),
                             ratio_phase_deg=np.degrees(np.angle(ratio))))
    del basis

# extra check with 4000/8000 (hole 2 deg, kR=1.5, l=0,1)
log('\n  L-convergence check with L=4000/8000 (hole 2 deg, N=32):')
hole = 2
th0 = np.deg2rad(180 - hole)
basis = cs.CapBasis(th0, NA, 8000, [0])
Tm = cs.whitening(th0, NA, 0)
for kR in (1.5,):
    _, lam = cs.bessel_terms(kR, 8000)
    _, lam_s = cs.bessel_terms(kR, 3)
    for l in (0, 1):
        Pl = sps.eval_legendre(l, basis.u)
        b = R * R * 2 * np.pi * (basis.profile(basis.s, 0) * Pl * basis.wu).sum(axis=1)
        nrm = R * R * 2 * np.pi * np.sum(Pl ** 2 * basis.wu)
        za = zc_whitened(basis, Tm, R, kR, b, 4000, lam)
        zb = zc_whitened(basis, Tm, R, kR, b, 8000, lam)
        zx = richardson(za, zb)
        zxold = resA[(kR, l, hole)][0] * nrm
        log(f'    kR={kR} l={l}: zcap(Rich 2000/4000)={zxold/nrm:.6f}  zcap(Rich 4000/8000)={zx/nrm:.6f}  '
            f'rel diff={abs(zx-zxold)/abs(zx):.2e}')
        rows.append(dict(test='A_Lcheck', kR=kR, l=l, hole_deg=hole, L1=4000, L2=8000,
                         zcap_re=(zx / nrm).real, zcap_im=(zx / nrm).imag,
                         Zc_L1_re=(za / nrm).real, Zc_L1_im=(za / nrm).imag,
                         Zc_L2_re=(zb / nrm).real, Zc_L2_im=(zb / nrm).imag))
del basis

log('\n  Table A1: zcap / zsph (Richardson L=2000/4000)  [|ratio| , phase deg ; re , im]')
for kR in kRs:
    for l in ells:
        zs = resA[(kR, l, 10)][1]
        log(f'  kR={kR}  l={l}   zsph = {zs.real:+.5f} {zs.imag:+.5f}i   (sign Im zsph = {np.sign(zs.imag):+.0f})')
        for hole in holes:
            zc, zsph, ratio, z1n, z2n = resA[(kR, l, hole)]
            log(f'      hole {hole:2d} deg: zcap = {zc.real:+.5f} {zc.imag:+.5f}i | ratio = {ratio.real:+.5f} {ratio.imag:+.5f}i '
                f'| |r|={abs(ratio):.5f} ph={np.degrees(np.angle(ratio)):+.3f} deg '
                f'| raw L1: {z1n.real:+.5f}{z1n.imag:+.5f}i  raw L2: {z2n.real:+.5f}{z2n.imag:+.5f}i')

log('\n  Table A2: extrapolation of ratio to zero hole')
log('   linear from holes (5,2) deg: r0 = r2 - 2 (r5-r2)/3 ;  quadratic from (10,5,2) also shown')
for kR in kRs:
    for l in ells:
        r10, r5, r2 = (resA[(kR, l, h)][2] for h in (10, 5, 2))
        r0 = r2 - 2.0 * (r5 - r2) / 3.0
        # quadratic through (10,r10),(5,r5),(2,r2) at 0 (Lagrange)
        hs = np.array([10., 5., 2.])
        rs = np.array([r10, r5, r2])
        r0q = sum(rs[i] * np.prod([(0 - hs[j]) / (hs[i] - hs[j]) for j in range(3) if j != i]) for i in range(3))
        dev = abs(r0 - 1)
        ok = check(f'A kR={kR} l={l} |r0-1|', dev, 1e-2)
        log(f'   kR={kR} l={l}: ratios 10/5/2 deg = {r10.real:+.4f}{r10.imag:+.4f}i / {r5.real:+.4f}{r5.imag:+.4f}i / '
            f'{r2.real:+.4f}{r2.imag:+.4f}i ; linear r0 = {r0.real:+.5f}{r0.imag:+.5f}i (|r0-1|={dev:.2e}) '
            f'; quad r0 = {r0q.real:+.5f}{r0q.imag:+.5f}i (|r0q-1|={abs(r0q-1):.2e}) -> {"PASS" if ok else "FAIL"} (tol 1e-2)')
        rows.append(dict(test='A_extrap', kR=kR, l=l, ratio_re=r0.real, ratio_im=r0.imag, ratio_abs=abs(r0),
                         dev_from_1=dev, quad_re=r0q.real, quad_im=r0q.imag))
    # power ratio (real part), the Exp 2 style measure
    for l in ells:
        pr = [(resA[(kR, l, h)][0].real / resA[(kR, l, h)][1].real) for h in (10, 5, 2)]
        log(f'   kR={kR} l={l}: POWER ratio Re(zcap)/Re(zsph) hole 10/5/2 = {pr[0]:.4f} / {pr[1]:.4f} / {pr[2]:.4f}')

# ============================================================================ TEST B
log('\n' + '=' * 100)
log('TEST B: flat-disc limit, rigid disc broadside added mass m_a = (8/3) rho0 a^3  (uniform v=1, m=0)')
log('  m_a = -Im(Zc)/w ;  Zc = int mu dS = force.  Resistance vs low-ka dipole formula (see below).')
log('=' * 100)
# Independent derivation of the dipole radiation formula: point force F (peak) radiates P = k^2 |F|^2 /(24 pi rho0 c0).
# Verify it against the closed rigid-sphere exact lam_1 (v = cos theta, F = Zc):
for x in (0.05, 0.1):
    _, ls = cs.bessel_terms(x, 3)
    Zs = RHO0 * C0 / ls[1] * (4 * np.pi / 3)       # R=1, int cos^2 dS = 4 pi/3
    Pexact = 0.5 * Zs.real
    Pdip = x ** 2 * abs(Zs) ** 2 / (24 * np.pi * RHO0 * C0)
    log(f'  dipole-formula sanity on closed sphere l=1, kR={x}: P_exact={Pexact:.6e}, k^2|F|^2/(24 pi rho c)={Pdip:.6e}, ratio={Pexact/Pdip:.5f}'
        f' ; m_a/(0.5 rho V)= {-Zs.imag/(x*C0)/(0.5*RHO0*4*np.pi/3):.5f}')
    check(f'B dipole formula vs sphere kR={x} |ratio-1|', Pexact / Pdip - 1, 5e-2)

NB = 24
Rb = 10.0
resB = {}
for ar, (L1, L2) in [(0.2, (2000, 4000)), (0.1, (4000, 8000)), (0.05, (8000, 16000))]:
    th0 = np.arcsin(ar)
    t = time.time()
    basis = cs.CapBasis(th0, NB, L2, [0])
    Tm = cs.whitening(th0, NB, 0)
    a = Rb * np.sin(th0)
    b = Rb * Rb * 2 * np.pi * (basis.profile(basis.s, 0) * basis.wu).sum(axis=1)
    area = Rb * Rb * 2 * np.pi * np.sum(basis.wu)
    log(f'\n  a/R={ar}: theta0={np.degrees(th0):.4f} deg, a={a:.4f}, area={area:.5f} (pi a^2={np.pi*a*a:.5f}), '
        f'L={L1}/{L2}, basis {time.time()-t:.1f}s, rank {Tm.shape[0]}/{NB}')
    for ka in (0.02, 0.05):
        k = ka / a
        kR = k * Rb
        w = k * C0
        _, lam = cs.bessel_terms(kR, L2)
        z1 = zc_whitened(basis, Tm, Rb, k, b, L1, lam)
        z2 = zc_whitened(basis, Tm, Rb, k, b, L2, lam)
        zx = richardson(z1, z2)
        ma = lambda z: -z.imag / w
        ma_ref = 8.0 / 3.0 * RHO0 * a ** 3
        Rdip = RHO0 * C0 * k ** 4 * a ** 6 * 16 / (27 * np.pi)
        Pdip = k ** 2 * abs(zx) ** 2 / (24 * np.pi * RHO0 * C0)
        resB[(ar, ka)] = dict(ma_ratio=ma(zx) / ma_ref, zx=zx, rdip=Rdip)
        log(f'    ka={ka} (kR={kR:.4f}): Zc raw L1 = {z1.real:.5e}{z1.imag:+.5e}i ; raw L2 = {z2.real:.5e}{z2.imag:+.5e}i ; Rich = {zx.real:.5e}{zx.imag:+.5e}i')
        log(f'       m_a/((8/3) rho a^3): raw L1 {ma(z1)/ma_ref:.6f}  raw L2 {ma(z2)/ma_ref:.6f}  Richardson {ma(zx)/ma_ref:.6f}   Im(Zc)<0: {zx.imag < 0}')
        log(f'       Re(Zc)/[rho c (16/27pi) k^4 a^6] = {zx.real/Rdip:.5f} ; Re(Zc)/2 vs k^2|Zc|^2/(24 pi rho c) ratio = {0.5*zx.real/Pdip:.5f}')
        rows.append(dict(test='B', a_over_R=ar, ka=ka, kR=kR, L1=L1, L2=L2,
                         Zc_L1_re=z1.real, Zc_L1_im=z1.imag, Zc_L2_re=z2.real, Zc_L2_im=z2.imag,
                         Zc_x_re=zx.real, Zc_x_im=zx.imag, ma_ratio_L1=ma(z1) / ma_ref, ma_ratio_L2=ma(z2) / ma_ref,
                         ma_ratio_x=ma(zx) / ma_ref, Re_over_dipole=zx.real / Rdip, P_over_dipP=0.5 * zx.real / Pdip,
                         ImZc_negative=bool(zx.imag < 0)))
    del basis

log('\n  Table B1: m_a / ((8/3) rho0 a^3)  and (a/R)^2 extrapolation (from a/R = 0.1, 0.05)')
for ka in (0.02, 0.05):
    r = {ar: resB[(ar, ka)]['ma_ratio'] for ar in (0.2, 0.1, 0.05)}
    r0 = r[0.05] - (r[0.1] - r[0.05]) / 3.0                                  # (a/R)^2 extrapolation
    r0lin = r[0.05] - (r[0.1] - r[0.05]) * (0.05 / (0.1 - 0.05))              # linear in a/R, for comparison
    log(f'   ka={ka}: ratios a/R=0.2/0.1/0.05: {r[0.2]:.5f} / {r[0.1]:.5f} / {r[0.05]:.5f} ; (a/R)^2 extrap = {r0:.5f} (dev {r0-1:+.2e}) ; linear extrap = {r0lin:.5f}')
    check(f'B ka={ka} |ma ratio -1| (a/R)^2 extrap', r0 - 1, 1e-2)
    rows.append(dict(test='B_extrap', ka=ka, ma_ratio_x=r0, dev_from_1=r0 - 1, ma_ratio_linear_extrap=r0lin))
# ka dependence at fixed a/R (finite-frequency correction)
for ar in (0.2, 0.1, 0.05):
    log(f'   a/R={ar}: ma ratio ka=0.02 -> {resB[(ar,0.02)]["ma_ratio"]:.5f}, ka=0.05 -> {resB[(ar,0.05)]["ma_ratio"]:.5f}')

# ============================================================================ TEST C
log('\n' + '=' * 100)
log('TEST C: sign / convention sanity (e^{-i w t}: mass-like <=> Im Zc < 0)')
log('=' * 100)
allneg = all(resB[key]['zx'].imag < 0 for key in resB)
log(f'  Test B: Im(Zc) < 0 for all {len(resB)} (a/R, ka) cases: {allneg}')
checks.append(('C Im Zc<0 all disc cases', float(allneg), 0.5, allneg))
ok0 = all(resB[key]['zx'].real > 0 for key in resB)
log(f'  Test B: Re(Zc) > 0 (radiated power positive): {ok0}')
checks.append(('C Re Zc>0 all disc cases', float(ok0), 0.5, ok0))
for kR in kRs:
    s = ', '.join(f'l={l}: Im(zsph) {"<" if resA[(kR,l,10)][1].imag<0 else ">"} 0' for l in ells)
    c = ', '.join(f'l={l}: Im(zcap@2deg) {"<" if resA[(kR,l,2)][0].imag<0 else ">"} 0' for l in ells)
    log(f'  closed sphere kR={kR}: {s}')
    log(f'  cap (hole 2 deg) kR={kR}: {c}')
log('  Note: for the closed sphere, l=0 (breathing) and l>=1 modes at low kR have different reactive character '
    '(compressible interior stiffness vs mass-like exterior+interior inertia); signs listed above are those of the exact lam_l.')

# ============================================================================ summary
log('\n' + '=' * 100)
log('SUMMARY')
log('=' * 100)
for name, val, tol, ok in checks:
    log(f'  [{"PASS" if ok else "FAIL"}] {name}: value {val:+.3e}, tol {tol:g}')
npass = sum(1 for c in checks if c[3])
log(f'  {npass}/{len(checks)} checks passed. Runtime {time.time()-T0:.0f} s.')

keys = []
for r in rows:
    for kk in r:
        if kk not in keys:
            keys.append(kk)
with open(os.path.join(OUT, 'verify.csv'), 'w', newline='') as f:
    wr = csv.DictWriter(f, fieldnames=keys)
    wr.writeheader()
    wr.writerows(rows)
log('wrote', os.path.join(OUT, 'verify.csv'))
logf.close()
