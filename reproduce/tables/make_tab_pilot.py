"""Table 3 (tab:pilot): regime study on seven cells.

Reads  results/r0_revision/pilot/pilot_runs.csv (written by tests/r0_pilot_stats.py)
Writes $REPRO_OUT/tables/tab_pilot.csv, tab_pilot.tex; compares with the typed table in the manuscript.

Column rules (trace comment at ms l.391, P9):
  ERP s = regret_bound_dB of run EU; sinc s = regret_bound_dB of SU;
  exact spread = |regret_bound_dB(XU) - regret_bound_dB(XS)|; sinc e_red = pred_minus_true_dB of SU.
Cell order as in Fig. 7 (reproduce/figures/make_fig_pilot_f17.py).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tabutil as T  # noqa: E402
from tabutil import C  # noqa: E402

ORDER = ["c30_k1", "c90_k1p5", "c120_k1p5", "c30_k3", "c90_k7", "c120_k7", "c30_k18"]
R = {(r["cell"], r["run"]): r for r in T.rd(os.path.join(C.PIL, "pilot_runs.csv"))}
g = lambda c, run, col: float(R[(c, run)][col])

rows = []
for c in ORDER:
    rows.append([float(R[(c, "EU")]["theta0_deg"]), float(R[(c, "EU")]["kR"]),
                 g(c, "EU", "regret_bound_dB"), g(c, "SU", "regret_bound_dB"),
                 abs(g(c, "XU", "regret_bound_dB") - g(c, "XS", "regret_bound_dB")),
                 g(c, "SU", "pred_minus_true_dB")])

stops = {(c, run): R[(c, run)]["stop"] for c in ORDER for run in ("EU", "SU", "XU")}
print("runs from uniform that met the stopping rule:", sorted(k for k, v in stops.items() if v != "nmax"))
fmt = [T.f(0), lambda v: "%g" % round(v, 1), T.f(1), T.f(1), T.f(1), T.f(1)]
T.write_table("tab_pilot", ["theta0_deg", "kR", "erp_s", "sinc_s", "exact_spread", "sinc_e_red"], rows, fmt)

typed = T.manuscript_table("tab:pilot")
n, bad = T.report("tab:pilot", rows, typed, [0, 1, 2, 3, 4, 5])
sys.exit(1 if bad else 0)
