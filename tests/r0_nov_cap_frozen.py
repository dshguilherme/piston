"""Frozen spherical-cap fields: spectral two-sided solution vs sinc model, validation vs bempp, and kR sweep.

Uses the exported P1 traces of the frozen cap fields (results/r0_revision/bem_mesh/spherecap_m*.mat, the SAME
files used for the bempp-cl open-screen BEM) and the spectral cap solver in r0_nov_capspec.py.
Outputs (results/r0_revision/novelty/):
  cap_frozen_validation.csv   spectral vs bempp at 100 Hz on the m2/m3 meshes, and m8 (finest trace)
  cap_frozen_sweep.csv        P_sinc, P_two_sided, ratio for the six frozen fields over kR (fields frozen at 100 Hz)
Run:  python -u tests/r0_nov_cap_frozen.py [validate|sweep|all]
"""
import os, sys, csv, time
import numpy as np
import scipy.io as sio

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from r0_nov_capspec import CapBasis, bessel_terms, legendre_rows, RHO0, C0

ROOT = os.path.dirname(HERE)
MESH = os.path.join(ROOT, "results", "r0_revision", "bem_mesh")
OUT = os.path.join(ROOT, "results", "r0_revision", "novelty")
os.makedirs(OUT, exist_ok=True)
R = 10.0
TH0 = np.radians(30.0)
K100 = 2 * np.pi * 100.0 / C0
# Dunavant degree-4 six-point rule (barycentric, weights sum to 1)
_a, _b = 0.445948490915965, 0.091576213509771
_wa, _wb = 0.223381589678011, 0.109951743655322
BARY = np.array([[1 - 2 * _a, _a, _a], [_a, 1 - 2 * _a, _a], [_a, _a, 1 - 2 * _a],
                 [1 - 2 * _b, _b, _b], [_b, 1 - 2 * _b, _b], [_b, _b, 1 - 2 * _b]])
WB = np.array([_wa] * 3 + [_wb] * 3)


def load_mesh(level):
    d = sio.loadmat(os.path.join(MESH, f"spherecap_m{level}.mat"))
    V = d["vertices"].astype(float)
    T = d["triangles"].astype(np.int64)
    if T.min() == 1:
        T = T - 1
    vn = d["vn"]
    labels = [str(x[0]) if hasattr(x, "__len__") else str(x) for x in d["labels"].ravel()]
    return V, T, vn, labels


def quad_points(V, T, vn):
    A = V[T]                                           # (nt,3,3)
    area = 0.5 * np.linalg.norm(np.cross(A[:, 1] - A[:, 0], A[:, 2] - A[:, 0]), axis=1)
    X = np.einsum("qb,tbd->tqd", BARY, A).reshape(-1, 3)
    vq = np.einsum("qb,tbf->tqf", BARY, vn[T]).reshape(-1, vn.shape[1])
    r = np.linalg.norm(X, axis=1)
    w = (area[:, None] * WB[None]).ravel() * (R / r) ** 2      # radial projection onto the sphere
    th = np.arccos(np.clip(X[:, 2] / r, -1, 1))
    ph = np.arctan2(X[:, 1], X[:, 0])
    return th, ph, w, vq


def project(th, ph, w, vq, basis, M, Ls, chunk=40000):
    """b[m] (N x nf) = int psi_n e^{-i m phi} v dS ;  Vlm[m] (Ls+1 x nf) = int v conj(Y_lm) dOmega."""
    nf = vq.shape[1]
    b, Vl = {}, {}
    u = np.cos(th)
    for m in range(-M, M + 1):
        bm = np.zeros((basis.N, nf), complex)
        vm = np.zeros((Ls + 1, nf), complex)
        for s in range(0, len(th), chunk):
            sl = slice(s, s + chunk)
            e = np.exp(-1j * m * ph[sl])[:, None] * vq[sl] * w[sl][:, None]      # (nc, nf)
            psi = basis.profile_theta(th[sl], m)                                  # (N, nc)
            bm += psi @ e
            ee = e / (R * R) / np.sqrt(2 * np.pi)
            for l, Pl in legendre_rows(u[sl], m, Ls):
                vm[l] += Pl @ ee
        b[m], Vl[m] = bm, vm
    return b, Vl


def powers(basis, b, Vl, k, Ls_use=None, Llist=(4000, 8000)):
    """Return P_sinc (nf,), dict L -> P_2s (nf,), and per-|m| sinc/2s contributions."""
    x = k * R
    Ls = max(v.shape[0] for v in Vl.values()) - 1
    wsinc, _ = bessel_terms(x, Ls)
    Psinc = RHO0 * C0 * R * R * sum((wsinc[:, None] * np.abs(v) ** 2).sum(0) for v in Vl.values())
    out = {}
    Lmax = max(Llist)
    _, lam = bessel_terms(x, Lmax)
    for L in Llist:
        P = 0.0
        for m, bm in b.items():
            F = basis.proj[abs(m)][:, : L + 1]
            A = R * R * (F * lam[None, : L + 1]) @ F.T
            c = np.linalg.solve(A, RHO0 * C0 * bm)
            P = P + 0.5 * np.real(np.sum(np.conj(bm) * c, axis=0))
        out[L] = P
    return Psinc, out


_BASIS = {}


def get_basis(N, M, Lmax):
    """CapBasis (level independent), cached in memory and on disk."""
    key = (N, M, Lmax)
    if key in _BASIS:
        return _BASIS[key]
    fn = os.path.join(OUT, f"cap_basis_N{N}_M{M}_L{Lmax}.npz")
    basis = CapBasis(TH0, N, 4, [0])                      # cheap shell, projections replaced below
    basis.L = Lmax
    if os.path.exists(fn):
        d = np.load(fn)
        basis.proj = {int(k[1:]): d[k] for k in d.files}
    else:
        full = CapBasis(TH0, N, Lmax, range(0, M + 1))
        basis.proj = full.proj
        np.savez(fn, **{f"m{m}": F for m, F in basis.proj.items()})
    _BASIS[key] = basis
    return basis


def build(level, N, M, Lmax, Ls):
    V, T, vn, labels = load_mesh(level)
    th, ph, w, vq = quad_points(V, T, vn)
    basis = get_basis(N, M, Lmax)
    b, Vl = project(th, ph, w, vq, basis, M, Ls)
    return basis, b, Vl, labels


BEM = {  # bempp-cl open-screen powers (results/r0_revision/bem_shell_levels.csv / bem_comparison.csv), base quadrature
    2: dict(uniform=5.88699e-7, dyncomp=2.61002e-8, inputpower=2.76502e-8, **{"quadvel-rad": 1.16063e-8,
            "quadvel-I": 1.70919e-8, "quadvel-ones": 5.84103e-8}),
    3: dict(uniform=5.92110e-7, dyncomp=2.61506e-8, inputpower=2.77126e-8, **{"quadvel-rad": 1.16806e-8,
            "quadvel-I": 1.71501e-8, "quadvel-ones": 5.84640e-8}),
    "extrap": dict(uniform=5.98965e-7, dyncomp=2.62519e-8, inputpower=2.78379e-8, **{"quadvel-rad": 1.18298e-8,
                   "quadvel-I": 1.72670e-8, "quadvel-ones": 5.85720e-8}),
}
SINC_CONSISTENT = dict(uniform=3.77918e-7, dyncomp=1.37686e-8, inputpower=1.49608e-8, **{"quadvel-rad": 7.18504e-9,
                       "quadvel-I": 9.64532e-9, "quadvel-ones": 2.87996e-8})


def validate(N=48, M=20, Lmax=8000, Ls=110):
    rows = []
    for level in (2, 3, 8):
        t = time.time()
        basis, b, Vl, labels = build(level, N, M, Lmax, Ls)
        # azimuthal content
        em = {m: float(np.sum(np.abs(b[m]) ** 2)) for m in b}
        etot = sum(em.values())
        tail = sum(v for m, v in em.items() if m != 0) / etot
        Psinc, P2 = powers(basis, b, Vl, K100)
        P2x = 2 * P2[8000] - P2[4000]
        print(f"level m{level}: build {time.time()-t:.0f}s, |b|^2 fraction in m!=0: {tail:.2e}", flush=True)
        for i, lab in enumerate(labels):
            row = dict(level=level, design=lab, N=N, M=M, P_sinc_spec=Psinc[i], P2s_L4000=P2[4000][i],
                       P2s_L8000=P2[8000][i], P2s_extrap=P2x[i], ratio_sinc_over_2s=Psinc[i] / P2x[i],
                       L_rel_change=(P2[8000][i] - P2[4000][i]) / P2x[i],
                       P_bem_same_level=BEM.get(level, {}).get(lab, np.nan), P_bem_extrap=BEM["extrap"][lab],
                       P_sinc_consistent=SINC_CONSISTENT[lab], az_tail_fraction=tail)
            row["spec_over_bem_same_level"] = row["P2s_extrap"] / row["P_bem_same_level"]
            row["spec_over_bem_extrap"] = row["P2s_extrap"] / row["P_bem_extrap"]
            row["sinc_spec_over_consistent"] = row["P_sinc_spec"] / row["P_sinc_consistent"]
            rows.append(row)
            print("  {design:13s} sinc {P_sinc_spec:.5e} (/cons {sinc_spec_over_consistent:.5f})  2s {P2s_extrap:.5e} "
                  "(L-chg {L_rel_change:+.1e})  /bem_lvl {spec_over_bem_same_level:.5f}  /bem_ext {spec_over_bem_extrap:.5f}"
                  "  ratio {ratio_sinc_over_2s:.4f}".format(**row), flush=True)
        if level == 8:
            np.savez(os.path.join(OUT, f"cap_m8_projections_M{M}.npz"), labels=labels,
                     **{f"b{m}": b[m] for m in b}, **{f"V{m}": Vl[m] for m in Vl})
    with open(os.path.join(OUT, f"cap_frozen_validation_M{M}.csv"), "w", newline="") as f:
        wtr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wtr.writeheader(); wtr.writerows(rows)
    # per-design azimuthal energy beyond |m| = M-4 (truncation check), m8
    em = np.array([np.sum(np.abs(b[m]) ** 2, axis=0) for m in range(-M, M + 1)])
    hi = em[np.abs(np.arange(-M, M + 1)) >= M - 4].sum(0) / em.sum(0)
    nz = 1 - em[M] / em.sum(0)
    print("m8 per-design |b|^2 fraction m!=0:", " ".join(f"{l}:{f:.2e}" for l, f in zip(labels, nz)))
    print(f"m8 per-design |b|^2 fraction |m|>={M-4}:", " ".join(f"{l}:{f:.2e}" for l, f in zip(labels, hi)))


def sweep(N=48, M=20, Lmax=8000, Ls=110, kRs=None):
    kRs = np.array(kRs if kRs is not None else
                   np.r_[0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, np.arange(4.0, 40.01, 1.0), 18.318324510727656])
    kRs = np.unique(np.round(kRs, 9))
    assert Ls >= kRs.max() + 60
    fn = os.path.join(OUT, f"cap_m8_projections_M{M}.npz")
    basis = get_basis(N, M, Lmax)
    if os.path.exists(fn):
        d = np.load(fn)
        labels = [str(x) for x in d["labels"]]
        b = {m: d[f"b{m}"] for m in range(-M, M + 1)}
        Vl = {m: d[f"V{m}"][: Ls + 1] for m in range(-M, M + 1)}
    else:
        basis, b, Vl, labels = build(8, N, M, Lmax, Ls)
    rows = []
    for kR in kRs:
        k = kR / R
        Psinc, P2 = powers(basis, b, Vl, k)
        P2x = 2 * P2[8000] - P2[4000]
        for i, lab in enumerate(labels):
            rows.append(dict(kR=kR, freq_Hz=k * C0 / (2 * np.pi), design=lab, P_sinc=Psinc[i], P_2s=P2x[i],
                             P_2s_L8000=P2[8000][i], ratio=Psinc[i] / P2x[i],
                             L_rel_change=(P2[8000][i] - P2[4000][i]) / P2x[i]))
        print(f"kR {kR:7.3f}: " + " ".join(f"{r['ratio']:.3f}" for r in rows[-len(labels):]), flush=True)
    with open(os.path.join(OUT, "cap_frozen_sweep.csv"), "w", newline="") as f:
        wtr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wtr.writeheader(); wtr.writerows(rows)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode in ("validate", "all"):
        validate()
    if mode in ("sweep", "all"):
        sweep()


