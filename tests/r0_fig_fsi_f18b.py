"""Plot data for the JSV paper figure F18b (`fig_fsi`, label fig:fsi): fluid loading, small multiples.

Redesign of F18 (FIG9-FIX-NOTES.md). Reads only stored results; nothing is re-evaluated.
One small panel per cell: 30 deg kR 1, 90 deg kR 1.5, 120 deg kR 1.5, 30 deg kR 3 and the design point
(30 deg kR 18.3). In each panel the x axis is the medium (1 in vacuo, 2 scaled medium, 3 air density) and the
y axis is Delta L, the exact radiated power relative to the uniform design in the same medium (dB).

Column rules (the same as make_f18.py / reproduce/figures/make_fig_fsi_f18.py, so no number changes):
  exact   min over the two exact-driven designs XU, XS of dB_vs_uniform_<variant>   (fsi_summary.csv)
  sinc    dB_vs_uniform_<variant> of SU, the sinc design reached from the uniform start (fsi_summary.csv)
  reopt   30 deg kR 3 only: the exact design re-optimized with fluid loading in that medium, evaluated in the
          same medium, dB_vs_uniform_same_variant of reopt_dev_XU (dev) and reopt_phys_XU (phys)
          (reopt_eval_c30_k3.csv); nan elsewhere

Inputs (results/r0_revision/fsi/): fsi_summary.csv, reopt_eval_c30_k3.csv
Outputs (paper_jsv/prototype/figs/data/): f18b_<cell>.dat with columns  med exact sinc reopt
Run:  python tests/r0_fig_fsi_f18b.py      (standard library only)
"""
import csv
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results", "r0_revision", "fsi")
O = os.path.join(ROOT, "paper_jsv", "prototype", "figs", "data")
CELLS = ["c30_k1", "c90_k1p5", "c120_k1p5", "c30_k3", "c30_k18"]
V = ["vac", "dev", "phys"]

rows = list(csv.DictReader(open(os.path.join(R, "fsi_summary.csv"), newline="")))
reo = {(x["design"], x["variant"]): float(x["dB_vs_uniform_same_variant"])
       for x in csv.DictReader(open(os.path.join(R, "reopt_eval_c30_k3.csv"), newline=""))}

for c in CELLS:
    d = {r["design"]: r for r in rows if r["cell"] == c}
    out = ["med exact sinc reopt"]
    for i, v in enumerate(V, 1):
        ex = min(float(d["XU"]["dB_vs_uniform_" + v]), float(d["XS"]["dB_vs_uniform_" + v]))
        sn = float(d["SU"]["dB_vs_uniform_" + v])
        rp = reo.get(("reopt_%s_XU" % v, v)) if c == "c30_k3" else None
        out.append("%d %.2f %.2f %s" % (i, ex, sn, "nan" if rp is None else "%.2f" % rp))
    with open(os.path.join(O, "f18b_%s.dat" % c), "w", newline="\n") as f:
        f.write("\n".join(out) + "\n")
    print(c, out[1:])
