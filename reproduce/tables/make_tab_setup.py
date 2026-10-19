"""Table 1 (tab:setup): check each typed setting against the code and logs that define it.

Table 1 is typed in the manuscript; nothing computes it. This driver finds, for every setting, the code
line that sets it (regex search in the original scripts, read only) and checks the counts that follow by
arithmetic. A FAIL means the code no longer says what the table says (or the line moved beyond recognition).
Writes $REPRO_OUT/tables/tab_setup.csv with item, table value, file, line number and the matched code text.
"""
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tabutil as T  # noqa: E402
from tabutil import C  # noqa: E402

PR = "tests/r0_pilot_runs.m"
TR = "tests/r0_pilot_track.m"
XF = "tests/r0_opt_exact_factor.py"
XS = "tests/r0_opt_export_samples.m"
XLOG = "results/r0_revision/exact_opt/cap_exact_factor_q8.log"

CHECKS = [
    ("Geometry: R = 10 m", "R=10 m", PR, r"spec\.R\s*=\s*10\s*;"),
    ("Geometry: theta0 = 30 deg (design point)", "30 deg", XS, r"'alpha',deg2rad\(30\)"),
    ("Geometry: design-point cell kR = 2 pi 100/343 x 10", "kR = 18.32", PR, r"case 'c30_k18',\s*th=30;\s*kR=2\*pi\*100/343\*10"),
    ("Shell model: quadratic splines", "degree 2", PR, r"'deg',\[2 2\]"),
    ("Shell model: 40 x 24 elements", "40x24", PR, r"'nsub',\[40 24\]"),
    ("Material: E, nu, density", "200 GPa, 0.3, 7850", PR, r"'E',2e11,'NU',0\.3,'RHO',7850"),
    ("Damping: beta = 0.1/omega", "beta = 0.1/omega", PR, r"ctx\.beta\s*=\s*0\.1/ctx\.w"),
    ("Load: apex ('pole'), axial, 1 N", "1 N at apex, axial", PR, r"spec\.loadtype='pole';\s*spec\.loaddir=\[0;0;1\]"),
    ("Load magnitude 1 N", "1 N", PR, r"'loadmag',1,"),
    ("Supports: rim (bc = 2)", "rim fixed", PR, r"spec\.bc\s*=\s*2\s*;"),
    ("Frequency 100 Hz", "100 Hz", PR, r"spec\.freq\s*=\s*100\s*;"),
    ("Bounds 25-100 mm, mean 50 mm", "25-100 mm, volfrac 0.5", PR, r"'tmin',0\.025,'tmax',0\.10,'volfrac',0\.5"),
    ("Cold start 50 mm (= volfrac x tmax)", "50 mm", TR, r"t\s*=\s*opt\.volfrac\*opt\.tmax\*ones\(ne,1\)"),
    ("Constraint not area-weighted", "vweights = ones", TR, r"opt\.vweights\s*=\s*ones\(ctx\.ne,1\)"),
    ("Objective: dB, scaled by first iteration", "100 (A + 10 log10 J)/p0", TR, r"f0s\s*=\s*100\*Jdb/p0"),
    ("Objective: A = 100, 10 log10", "Abar 100, dbfac 10", PR, r"'Abar',100,'dbfac',10"),
    ("GCMMA: 20 inner iterations (+ growth 1.3, raa 1e-3, not in Table 1)", "20 inner", TR,
     r"gopt\.kmax\s*=\s*20;\s*gopt\.growth\s*=\s*1\.30;\s*gopt\.raa0eps\s*=\s*1e-3;\s*gopt\.raaeps\s*=\s*1e-3"),
    ("GCMMA: 600 outer iterations", "600", PR, r"nmax\s*=\s*600"),
    ("Stopping rule tolerances", "1e-4, 1e-4, 3, 10, 1e-6", PR,
     r"'xchgtol',1e-4,'nhold',3,'window',10,'feastol',1e-6"),
    ("Stopping rule: objective tolerance 1e-4", "1e-4", PR, r"'stop',true,'ftol',1e-4"),
    ("Sinc operator: spline-consistent, 4x4 Gauss (angular 32x64, not in Table 1)", "4x4", PR,
     r"r0_build_consistent\(ctx,rhoc/cp,cp,\[4 4\],32,64\)"),
    ("Exact operator: N = 48, |m| <= 20, L = 8000", "48, 20", XF, r"N, M, LMAX, LS = 48, 20, 8000"),
    ("Exact operator: Richardson from L = 4000, 8000", "Hx = 2 H8000 - H4000", XF, r"Hx = 2 \* H\[8000\] - H\[4000\]"),
    ("Exact operator: negative eigenvalues set to zero", "clip", XF, r"ev = np\.clip\(ev, 0, None\)"),
    ("Exact operator: 8x8 Gauss = 61 440 points, 3024 dof", "61440, 3024", XLOG, r"nq=8: 61440 points, ndof 3024"),
    ("Exact operator: F_x is 1730 x 3024", "1730x3024", XLOG, r"Fx \(1730, 3024\)"),
]

rows, nfail = [], 0
for item, val, rel, rx in CHECKS:
    path = os.path.join(C.ROOT, *rel.split("/"))
    hit = None
    if os.path.exists(path):
        with open(path, encoding="utf-8", errors="replace") as fh:
            for i, ln in enumerate(fh, 1):
                if re.search(rx, ln):
                    hit = (i, ln.strip())
                    break
    ok = hit is not None
    nfail += not ok
    rows.append([item, val, rel, hit[0] if ok else "", hit[1] if ok else "NOT FOUND", "PASS" if ok else "FAIL"])
    print("%-4s %-72s %s:%s" % ("ok" if ok else "FAIL", item, rel, hit[0] if ok else "-"))

# counts that follow by arithmetic
arith = [("3024 coefficients = 3 x (40 + 2) x 24 (periodic quadratic space)", 3 * (40 + 2) * 24, 3024),
         ("960 thicknesses = 40 x 24", 40 * 24, 960),
         ("61 440 points = 40 x 24 x 8 x 8", 40 * 24 * 64, 61440)]
for item, got, want in arith:
    ok = got == want
    nfail += not ok
    rows.append([item, str(want), "arithmetic", "", str(got), "PASS" if ok else "FAIL"])
    print("%-4s %s" % ("ok" if ok else "FAIL", item))

d = C.out("tables")
with open(os.path.join(d, "tab_setup.csv"), "w", newline="\n") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["item", "table_value", "source", "line", "code_text", "status"])
    w.writerows(rows)
print("tab:setup: %d settings checked, %d FAIL; wrote %s" % (len(rows), nfail, os.path.join(d, "tab_setup.csv")))
sys.exit(1 if nfail else 0)
