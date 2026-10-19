"""Smoke test: closed-sphere modal formulas, recomputed in Python and compared with the stored MATLAB output.

Stored files (written by tests/r0_sphere_modal_spectrum.m and tests/r0_spherical_exterior_benchmark.m):
  results/r0_revision/sphere-modal-spectrum.csv       ka, l, 4 pi R^2 j_l^2, F_l = Re(-i h_l / h_l'), ratio = 2 x^2 j_l^2 / F_l
  results/r0_revision/sphere-modal-direct-check.csv   analytic quadratic 4 pi R^4 j_l(1)^2 for l = 0..5
  results/r0_revision/spherical-exterior-benchmark.csv analytic sinc integral S^2 j_0^2 and ratio for the breathing sphere
Checks:
  (a) every stored formula value equals the Python value to 1e-12 relative;
  (b) the two forms of the ratio agree: 2 x^2 j_l^2 / F_l = 2 x^4 j_l^2 |h_l'|^2 (Wronskian), as used in
      tests/r0_nov_sphere_factor2.py;
  (c) breathing sphere: F_0 = x^2/(1+x^2) (radiation efficiency), and the paper's kR = 1 values
      l = 0..5: 2.8323 0.9070 0.6850 0.6433 0.6136 0.5932 (sphere_factor2.log, "paper Table").
Prints one line per check and exits 1 on any failure. Runs in about a second.
"""
import csv
import os
import sys

import numpy as np
from scipy import special as sps

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
import config as C  # noqa: E402


def rd(name):
    with open(os.path.join(C.RES, name), newline="") as f:
        return list(csv.DictReader(f))


def modal(ell, x):
    jl = sps.spherical_jn(ell, x)
    hl = jl - 1j * sps.spherical_yn(ell, x)
    hp = sps.spherical_jn(ell, x, True) - 1j * sps.spherical_yn(ell, x, True)
    F = (-1j * hl / hp).real
    return jl, F, hp


fails = []


def check(name, err, tol):
    ok = err <= tol
    print("%-4s %-70s max rel err %.2e (tol %.0e)" % ("ok" if ok else "FAIL", name, err, tol))
    if not ok:
        fails.append(name)


rows = rd("sphere-modal-spectrum.csv")
e_eig = e_F = e_rat = e_wr = 0.0
for r in rows:
    x, ell = float(r["ka"]), int(r["degree_l"])
    jl, F, hp = modal(ell, x)
    e_eig = max(e_eig, abs(4 * np.pi * jl ** 2 / float(r["sinc_integral_operator_eigenvalue_R1"]) - 1))
    e_F = max(e_F, abs(F / float(r["normalized_exterior_radiation_resistance"]) - 1))
    rat = 2 * x ** 2 * jl ** 2 / F
    e_rat = max(e_rat, abs(rat / float(r["sinc_to_exterior_power_ratio"]) - 1))
    e_wr = max(e_wr, abs(2 * x ** 4 * jl ** 2 * abs(hp) ** 2 / rat - 1))
check("sphere-modal-spectrum.csv: 4 pi R^2 j_l^2 (%d rows)" % len(rows), e_eig, 1e-12)
check("sphere-modal-spectrum.csv: F_l = Re(-i h_l/h_l')", e_F, 1e-12)
check("sphere-modal-spectrum.csv: ratio 2 x^2 j_l^2 / F_l", e_rat, 1e-12)
check("Wronskian form 2 x^4 j_l^2 |h_l'|^2 equals 2 x^2 j_l^2 / F_l", e_wr, 1e-12)

e = 0.0
for r in rd("sphere-modal-direct-check.csv"):
    ell = int(r["degree_l"])
    e = max(e, abs(4 * np.pi * sps.spherical_jn(ell, 1.0) ** 2 / float(r["analytic_quadratic"]) - 1))
check("sphere-modal-direct-check.csv: 4 pi R^4 j_l(1)^2, l = 0..5", e, 1e-12)
dq = max(float(r["relative_quadrature_error"]) for r in rd("sphere-modal-direct-check.csv"))
check("stored direct-quadrature error (MATLAB, 20 x 40 rule), l = 0..5", dq, 1e-9)

e1 = e2 = 0.0
for r in rd("spherical-exterior-benchmark.csv"):
    x = float(r["ka"])
    S = 4 * np.pi
    e1 = max(e1, abs(S ** 2 * (np.sin(x) / x) ** 2 / float(r["analytic_sinc_integral"]) - 1))
    jl, F, hp = modal(0, x)
    e2 = max(e2, abs(2 * x ** 2 * jl ** 2 / F / float(r["analytic_ratio"]) - 1))
check("spherical-exterior-benchmark.csv: S^2 j_0^2", e1, 1e-12)
check("spherical-exterior-benchmark.csv: breathing-sphere ratio", e2, 1e-12)

xs = np.array([0.1, 0.5, 1, 2, 4, 8.0])
check("breathing sphere: F_0 = x^2/(1+x^2)", max(abs(modal(0, x)[1] / (x * x / (1 + x * x)) - 1) for x in xs), 1e-12)

paper = [2.8323, 0.9070, 0.6850, 0.6433, 0.6136, 0.5932]
got = [2 * modal(ell, 1.0)[0] ** 2 / modal(ell, 1.0)[1] for ell in range(6)]
err = max(abs(g - p) for g, p in zip(got, paper))
check("kR = 1 ratios l = 0..5 equal the printed 4-decimal values", err, 0.5e-4 + 1e-12)

print("SPHERE_MODAL %s (%d failures)" % ("PASS" if not fails else "FAIL", len(fails)))
sys.exit(1 if fails else 0)
