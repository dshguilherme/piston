"""Table F.1 (tab:fsi): fluid loading of the regime-study designs, and its footnote.

Reads  results/r0_revision/fsi/fsi_summary.csv       (written by tests/r0_fsi_summary.py)
       results/r0_revision/fsi/reopt_eval_c30_k3.csv  (written by tests/r0_fsi_reopt_eval.py)
Writes $REPRO_OUT/tables/tab_fsi.csv, tab_fsi.tex; compares with the typed table and footnote in the manuscript.

Column rules: the same as reproduce/figures/make_fig_fsi_f18.py (copy of the rescued make_f18.py), i.e.
  best exact Delta L = min over XU, XS of dB_vs_uniform_<v>; sinc s = regret_<v>_dB of SU;
  best sinc s = min over SU, SX of regret_<v>_dB; v in (vac, dev = scaled medium, phys = air density).
Footnote: Delta L of the re-optimised designs (reopt_dev_XU under dev, reopt_phys_XU under phys) and the
sinc gap 10 log10(P_SU / P_reopt) in the same medium.
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tabutil as T  # noqa: E402
from tabutil import C  # noqa: E402

ORDER = ["c30_k1", "c90_k1p5", "c120_k1p5", "c30_k3", "c90_k7", "c120_k7", "c30_k18"]
V = ["vac", "dev", "phys"]
summ = T.rd(os.path.join(C.FSI, "fsi_summary.csv"))
rows = []
for c in ORDER:
    d = {r["design"]: r for r in summ if r["cell"] == c}
    ex = [min(float(d["XU"]["dB_vs_uniform_" + v]), float(d["XS"]["dB_vs_uniform_" + v])) for v in V]
    s = [float(d["SU"]["regret_%s_dB" % v]) for v in V]
    b = [min(float(d["SU"]["regret_%s_dB" % v]), float(d["SX"]["regret_%s_dB" % v])) for v in V]
    rows.append([float(d["uniform"]["theta0_deg"]), float(d["uniform"]["kR"])] + ex + s + b)

fmt = [T.f(0), lambda v: "%g" % round(v, 1)] + [T.f(1)] * 9
T.write_table("tab_fsi", ["theta0_deg", "kR"] + ["%s_%s" % (q, v) for q in ("exact_dL", "sinc_s", "bestsinc_s") for v in V],
              rows, fmt)
typed = T.manuscript_table("tab:fsi")
n, bad = T.report("tab:fsi", rows, typed, list(range(11)))

# footnote
ro = {(r["design"], r["variant"]): r for r in T.rd(os.path.join(C.FSI, "reopt_eval_c30_k3.csv"))}
su = next(r for r in summ if r["cell"] == "c30_k3" and r["design"] == "SU")
foot = [float(ro[("reopt_dev_XU", "dev")]["dB_vs_uniform_same_variant"]),
        float(ro[("reopt_phys_XU", "phys")]["dB_vs_uniform_same_variant"]),
        10 * math.log10(float(su["P_dev"]) / float(ro[("reopt_dev_XU", "dev")]["P_Fx"])),
        10 * math.log10(float(su["P_phys"]) / float(ro[("reopt_phys_XU", "phys")]["P_Fx"]))]
print("footnote computed: dL scaled %.2f, air %.2f; s scaled %.2f, air %.2f" % tuple(foot))
line = T.manuscript_line(r"Re-optimized with fluid loading", after_label="tab:fsi")
if line is None:
    print("tab:fsi footnote: comparison SKIP")
else:
    tail = line.split("Re-optimized with fluid loading", 1)[1]
    printed = [m for m in T.NUM.findall(tail) if "." in m][:4]
    for v, p in zip(foot, printed):
        ok = T.agrees(v, p)
        print("  footnote %s vs computed %.3f: %s" % (p, v, "ok" if ok else "MISMATCH"))
        if not ok:
            bad.append(("footnote", p, v))
sys.exit(1 if bad else 0)
