"""Smoke test: boundary-element validation checked from STORED data only (no BEM solve).

Stored files:
  results/r0_revision/bem-validation.csv   written by src/r0_bem_validate.py (bempp-cl 0.4.2)
  results/r0_revision/novelty/disc_two_sided.csv  Hankel-domain flat-disc solution (tests/r0_nov_disc_hankel.py)
Checks:
  (a) closed breathing sphere, R = 1, v0 = 1: the stored reference power equals the closed form
      P = (1/2) rho0 c0 4 pi R^2 x^2/(1+x^2); the finest-mesh BEM error at kR = 1 is below 1% (the gate of the
      validation script) and the error falls monotonically with the mesh level at each kR; far/bem = 1 to 1e-3;
  (b) rigid disc, a = 1, broadside v = 1: the stored reference is the low-ka dipole asymptote
      P = (1/2) rho0 c0 (16/(27 pi)) k^4 a^6; the Hankel-domain solution at the smallest ka approaches it;
      the BEM screen converges monotonically towards it, and far/bem = 1 to 1e-8 (energy check of the screen).
The re-solve of the screen (nr = 6, 12) needs bempp-cl and is NOT part of the smoke test (brief: no BEM solves);
the audit re-ran it and matched the stored values to 6 digits.
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
import config as C  # noqa: E402

RHO0, C0 = 1.21, 343.0
fails = []


def check(name, ok, detail):
    print("%-4s %-66s %s" % ("ok" if ok else "FAIL", name, detail))
    if not ok:
        fails.append(name)


with open(os.path.join(C.RES, "bem-validation.csv"), newline="") as f:
    rows = list(csv.DictReader(f))
sph = [r for r in rows if r["case"] == "sphere_BM"]
disc = [r for r in rows if r["case"] == "disc_screen"]

e = 0.0
for r in sph:
    x = float(r["kR_or_ka"])
    pref = 0.5 * RHO0 * C0 * 4 * np.pi * x * x / (1 + x * x)
    e = max(e, abs(pref / float(r["P_reference"]) - 1))
check("sphere: stored P_reference = closed form (kR = 1, 2, 4)", e < 1e-12, "max rel err %.1e" % e)
for x in sorted({r["kR_or_ka"] for r in sph}, key=float):
    lv = sorted((int(r["level"]), abs(float(r["rel_error"])), float(r["far_over_bem"])) for r in sph if r["kR_or_ka"] == x)
    errs = [v[1] for v in lv]
    mono = all(a > b for a, b in zip(errs, errs[1:]))
    farok = max(abs(v[2] - 1) for v in lv) < 1e-3
    check("sphere kR=%s: BEM error falls with level; far/bem = 1" % x, mono and farok,
          "errors %s" % " ".join("%.1e" % v for v in errs))
fin = max((r for r in sph if float(r["kR_or_ka"]) == 1.0), key=lambda r: int(r["level"]))
check("sphere kR=1 finest level error < 1% (paper: 1.4e-4)", abs(float(fin["rel_error"])) < 0.01,
      "%.3e at level %s" % (abs(float(fin["rel_error"])), fin["level"]))

with open(os.path.join(C.NOV, "disc_two_sided.csv"), newline="") as f:
    hk = list(csv.DictReader(f))
e = 0.0
for r in disc:
    ka = float(r["kR_or_ka"])
    pas = 0.5 * RHO0 * C0 * 16 / (27 * np.pi) * ka ** 4
    e = max(e, abs(pas / float(r["P_reference"]) - 1))
check("disc: stored P_reference = dipole asymptote (ka = 0.05, 0.1, 0.2)", e < 1e-12, "max rel err %.1e" % e)
k0 = min(hk, key=lambda r: float(r["ka"]))
dev = abs(float(k0["ratio_two_sided_over_dipole_asymptote"]) - 1)
check("disc: Hankel solution -> dipole asymptote at ka = %s" % k0["ka"], dev < 1e-4, "|ratio-1| = %.1e" % dev)
e2 = 0.0
for r in hk:
    ka = float(r["ka"])
    pas = 0.5 * RHO0 * C0 * 16 / (27 * np.pi) * ka ** 4
    e2 = max(e2, abs(float(r["P_two_sided"]) / pas - float(r["ratio_two_sided_over_dipole_asymptote"])))
check("disc: Hankel CSV ratio column = P_two_sided / asymptote (8 printed decimals)", e2 < 5.1e-9, "max abs err %.1e" % e2)
for ka in sorted({r["kR_or_ka"] for r in disc}, key=float):
    lv = sorted((int(r["level"]), float(r["ratio"]), float(r["far_over_bem"])) for r in disc if r["kR_or_ka"] == ka)
    rat = [v[1] for v in lv]
    mono = all(a < b for a, b in zip(rat, rat[1:])) and rat[-1] < 1.0
    farok = max(abs(v[2] - 1) for v in lv) < 1e-8
    check("disc ka=%s: BEM/asymptote rises monotonically to below 1; far/bem = 1" % ka, mono and farok,
          "ratios %s" % " ".join("%.4f" % v for v in rat))
nr6 = next(r for r in disc if r["kR_or_ka"] == "0.05" and r["level"] == "6")
nr12 = next(r for r in disc if r["kR_or_ka"] == "0.05" and r["level"] == "12")
check("disc ka=0.05 nr=6, 12 equal the values the audit re-solved (2.14905e-4, 2.30291e-4)",
      abs(float(nr6["P_bem"]) / 2.14905e-4 - 1) < 5e-6 and abs(float(nr12["P_bem"]) / 2.30291e-4 - 1) < 5e-6,
      "%.6e %.6e" % (float(nr6["P_bem"]), float(nr12["P_bem"])))

print("BEM_STORED %s (%d failures)" % ("PASS" if not fails else "FAIL", len(fails)))
sys.exit(1 if fails else 0)
