"""Parameters of the cap blueprint (JSV paper, figure F19 `fig_cap_geometry`), read from the stored results.

Every number drawn in the blueprint comes from a stored result file, not from the manuscript text:

  results/r0_revision/spherecap_frozen_states.mat   'spec' struct of the structural model that produced the
                                                    frozen fields and the optimization runs: R, alpha (cap
                                                    half-angle), nsub (elements: meridian x circumference),
                                                    deg (spline degree), bc (clamped side), load type,
                                                    direction and magnitude, E, NU, RHO, freq
  results/r0_revision/exact_opt/run_X0.mat          'o' struct of the optimizer: tmin, tmax (thickness bounds),
                                                    volfrac (mean-thickness constraint as a fraction of tmax)
  results/r0_revision/exact_opt/cap_samples_q8.mat  ndof (displacement coefficients), ne (design variables),
                                                    total area of the acoustic quadrature

Conventions taken from the code (tests/obstacle_build_ctx.m, spherecap map and 'pole' load; tests/r0_pilot_setup.m):
  x = (R sin(th) cos(ph), R sin(th) sin(ph), R cos(th)), th = alpha*u (polar angle, apex at u = 0), ph = 2*pi*v;
  O is the origin, the apex is (0, 0, R); the GeoPDEs normal cross(dx/du, dx/dv) is the outward radial unit
  vector e_r (checked here numerically on the quadrature points); bc = 2 is the side u = 1 (the rim) with all
  three displacement components fixed (rotation-free Kirchhoff-Love: rotations are not degrees of freedom);
  the 'pole' load is a force of loadmag along loaddir = +z, shared equally by the control points of the u = 0 row.
  The damping of the runs is beta = 0.1/omega, set in tests/r0_opt_export_samples.m (cap_ctx), not by spec.zeta.

Output: paper_jsv/prototype/figs/data/f19_params.tex (TeX macros read by fig_cap_geometry.tex) and a printout.
Run:   python tests/r0_fig_geometry.py
"""
import os
import numpy as np
import scipy.io as sio

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results", "r0_revision")
OUT = os.path.join(ROOT, "paper_jsv", "prototype", "figs", "data")
C0 = 343.0  # sound speed of the paper (r0_nov_capspec.C0)


def scalar(a):
    return float(np.asarray(a).ravel()[0])


def main():
    spec = sio.loadmat(os.path.join(RES, "spherecap_frozen_states.mat"), squeeze_me=True,
                       struct_as_record=False)["spec"]
    o = sio.loadmat(os.path.join(RES, "exact_opt", "run_X0.mat"), squeeze_me=True, struct_as_record=False)["o"]
    smp = sio.loadmat(os.path.join(RES, "exact_opt", "cap_samples_q8.mat"),
                      variable_names=["ndof", "ne", "w", "x", "B"])
    R, al = float(spec.R), float(spec.alpha)
    n_mer, n_circ = [int(v) for v in np.ravel(spec.nsub)]
    deg = int(np.ravel(spec.deg)[0])
    ndof, ne = int(scalar(smp["ndof"])), int(scalar(smp["ne"]))
    area_q = float(np.sum(smp["w"]))
    # checks against the code conventions
    assert spec.loadtype == "pole" and np.allclose(np.ravel(spec.loaddir), [0, 0, 1]), "load is not +z at the apex"
    assert int(spec.bc) == 2, "support is not on side 2 (u = 1, the rim)"
    assert ne == n_mer * n_circ
    assert ndof == 3 * (n_mer + deg) * n_circ, "coefficient count does not match a periodic degree-2 space"
    x = smp["x"]; B = smp["B"].tocsr(); Ns = ndof // 3
    uz = np.zeros(ndof); uz[2 * Ns:] = 1.0                      # rigid unit translation along z
    nz = B @ uz                                                 # = n . e_z at every quadrature point
    cos_th = x[2] / np.linalg.norm(x, axis=0)
    assert np.max(np.abs(nz - cos_th)) < 1e-6, "normal is not the outward radial vector"
    th0_deg = np.degrees(al)
    a = R * np.sin(al)
    H = R * (1 - np.cos(al))
    arc = R * al
    area = 2 * np.pi * R * R * (1 - np.cos(al))
    assert abs(area - area_q) / area < 1e-9
    freq = float(spec.freq)
    k = 2 * np.pi * freq / C0
    vals = dict(
        fR=R, fThz=th0_deg, fA=a, fH=H, fArc=arc, fArea=area,
        fNmer=n_mer, fNcirc=n_circ, fDmer=th0_deg / n_mer, fDcirc=360.0 / n_circ,
        fLmer=arc / n_mer, fLcircRim=2 * np.pi * a / n_circ,
        fNdof=ndof, fNe=ne, fDeg=deg,
        fTmin=1e3 * float(o.tmin), fTmax=1e3 * float(o.tmax), fTmean=1e3 * float(o.volfrac) * float(o.tmax),
        fLoad=float(spec.loadmag), fFreq=freq, fkR=k * R, fLam=C0 / freq,
        fE=float(spec.E) / 1e9, fNu=float(spec.NU), fRho=float(spec.RHO),
    )
    fmt = dict(fR="{:.0f}", fThz="{:.0f}", fA="{:.2f}", fH="{:.2f}", fArc="{:.2f}", fArea="{:.1f}",
               fNmer="{:d}", fNcirc="{:d}", fDmer="{:.2f}", fDcirc="{:.0f}", fLmer="{:.3f}", fLcircRim="{:.2f}",
               fNdof="{:d}", fNe="{:d}", fDeg="{:d}", fTmin="{:.0f}", fTmax="{:.0f}", fTmean="{:.0f}",
               fLoad="{:.0f}", fFreq="{:.0f}", fkR="{:.2f}", fLam="{:.2f}", fE="{:.0f}", fNu="{:.1f}", fRho="{:.0f}")
    os.makedirs(OUT, exist_ok=True)
    fn = os.path.join(OUT, "f19_params.tex")
    with open(fn, "w", encoding="ascii") as f:
        f.write("% Generated by tests/r0_fig_geometry.py from spherecap_frozen_states.mat (spec),\n"
                "% exact_opt/run_X0.mat (o) and exact_opt/cap_samples_q8.mat. Do not edit by hand.\n")
        for key, v in vals.items():
            f.write("\\def\\%s{%s}\n" % (key, fmt[key].format(int(v) if fmt[key] == "{:d}" else v)))
    for key, v in vals.items():
        print(f"{key:10s} {fmt[key].format(int(v) if fmt[key] == '{:d}' else v)}")
    print("normal = outward radial e_r (max |n.e_z - cos(theta)| = %.1e); wrote %s"
          % (np.max(np.abs(nz - cos_th)), fn))


if __name__ == "__main__":
    main()
