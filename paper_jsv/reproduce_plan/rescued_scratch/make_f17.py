import csv, math, sys
base = r"C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course"
pil = base + r"\results\r0_revision\pilot"
out = base + r"\paper_jsv\prototype\figs\data"
order = ["c30_k1", "c90_k1p5", "c120_k1p5", "c30_k3", "c90_k7", "c120_k7", "c30_k18"]
runs = list(csv.DictReader(open(pil + r"\pilot_runs.csv")))
cells = list(csv.DictReader(open(pil + r"\pilot_cells_summary.csv")))
R = {(r["cell"], r["run"]): r for r in runs}
# panel (a): regret LB of the from-uniform run per rung; stop status (1 = hit cap)
with open(out + r"\f17_regret.dat", "w") as f:
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
    with open(out + r"\f17_pt_%s.dat" % tag, "w") as f:
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

