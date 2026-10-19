# REPRODUCE-PACKAGE COPY. Do not edit numerics here; edit the original and re-copy.
# Original: paper_jsv/reproduce_plan/rescued_scratch/make_f17.py  (sha256 fbb62598633108d1)
# Changed: original lines 2-4 (base/pil/out) -> config.py (pil = C.PIL, out = C.out("figs","data") = $REPRO_OUT/figs/data);
#          backslash string concatenations (pil + r"\\pilot_runs.csv", out + r"\\f17_regret.dat", ...) replaced by os.path.join
#          (needed for portability; same file names). Added "import os", sys.path insert and "import config as C".
#          open(..., "w") is KEPT in text mode exactly as in the original: on Windows the stored f17_*.dat files are CRLF
#          because of it, and changing it would break byte identity. Numerics, formats and column order unchanged.
# Reads: results/r0_revision/pilot/pilot_runs.csv, results/r0_revision/pilot/pilot_cells_summary.csv
# Writes: $REPRO_OUT/figs/data/f17_regret.dat, f17_pt_erp.dat, f17_pt_sinc.dat, f17_pt_exact.dat
import csv, math, sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
import config as C
pil = C.PIL
out = C.out("figs", "data")
order = ["c30_k1", "c90_k1p5", "c120_k1p5", "c30_k3", "c90_k7", "c120_k7", "c30_k18"]
runs = list(csv.DictReader(open(os.path.join(pil, "pilot_runs.csv"))))
cells = list(csv.DictReader(open(os.path.join(pil, "pilot_cells_summary.csv"))))
R = {(r["cell"], r["run"]): r for r in runs}
# panel (a): regret LB of the from-uniform run per rung; stop status (1 = hit cap)
with open(os.path.join(out, "f17_regret.dat"), "w") as f:
    f.write("idx theta kR erp sinc exact erpcap sinccap exactcap sincpick exactxs\n")
    for i, c in enumerate(order, 1):
        e, s, x = R[(c, "EU")], R[(c, "SU")], R[(c, "XU")]
        pick = [r for r in cells if r["cell"] == c and r["rung"] == "R1_sinc"][0]
        f.write("%d %s %s %.2f %.2f %.2f %d %d %d %.2f %.2f\n" % (
            i, e["theta0_deg"], e["kR"], float(e["regret_bound_dB"]), float(s["regret_bound_dB"]),
            float(x["regret_bound_dB"]), e["stop"] == "nmax", s["stop"] == "nmax", x["stop"] == "nmax",
            float(pick["pick_regret_lowerbound_dB"]), float(R[(c, "XS")]["regret_bound_dB"])))
        print(c, e["regret_bound_dB"], s["regret_bound_dB"], x["regret_bound_dB"], e["stop"], s["stop"], x["stop"], pick["pick_regret_lowerbound_dB"])
# panel (b): every run, pred vs true; group = low kR (<=3) or high (>=7)
for rung, tag in [("R0_ERP", "erp"), ("R1_sinc", "sinc"), ("R4_exact", "exact")]:
    with open(os.path.join(out, "f17_pt_%s.dat" % tag), "w") as f:
        f.write("pred true cap lowkR cell run\n")
        for r in runs:
            if r["rung"] != rung:
                continue
            f.write("%.3f %.3f %d %d %s %s\n" % (float(r["pred_dB_vs_uniform"]), float(r["true_dB_vs_uniform"]),
                    r["stop"] == "nmax", float(r["kR"]) <= 3, r["cell"], r["run"]))
# checks
ncap = sum(r["stop"] == "nmax" for r in runs); print("runs", len(runs), "at cap", ncap,
      "converged rungs", sorted({r["rung"] for r in runs if r["stop"] != "nmax"}))
for lo in (True, False):
    v = [float(r["pred_minus_true_dB"]) for r in runs if r["rung"] == "R1_sinc" and (float(r["kR"]) <= 3) == lo]
    print("sinc pred-true lowkR" if lo else "sinc pred-true highkR", min(v), max(v), "abs", min(map(abs, v)), max(map(abs, v)))
v = [float(r["pred_minus_true_dB"]) for r in runs if r["rung"] == "R0_ERP"]
print("ERP |pred-true|", min(map(abs, v)), max(map(abs, v)))
for c in order:
    xu, sx = R[(c, "XU")], R[(c, "SX")]
    print(c, "SX-XU true loss", float(sx["true_dB_vs_uniform"]) - float(xu["true_dB_vs_uniform"]),
          "sigma XU", xu["sigma_exact"], "XS", R[(c, "XS")]["sigma_exact"], "eL", xu["ref_eL"], R[(c, "XS")]["ref_eL"])
