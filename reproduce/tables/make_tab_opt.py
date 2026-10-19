"""Table 2 (tab:opt): cap optimisation at 100 Hz under the sinc and exact two-sided operators.

Reads  results/r0_revision/exact_opt/cap_reopt_summary.csv  (runs S0 X0 XS SX SR XR; written by tests/r0_opt_cap_eval.py)
       results/r0_revision/exact_opt/cap_path_powers.csv    (row label paper_quadvel-rad = frozen design G, and paper_uniform)
Writes $REPRO_OUT/tables/tab_opt.csv, tab_opt.tex; compares with the typed table in the manuscript.

Column rules (the manuscript's caption and its trace comment at ms l.320):
  Delta L_sinc = dB_sinc_vs_uniform, Delta L_2s = dB_exact_vs_uniform, rho = ratio (sinc/exact power),
  e_red = pred_minus_true_dB = Delta L_sinc - Delta L_2s, s = dB_exact(run) - dB_exact(X0).
Data names SR/XR are the paper's SG/XG (design "R" was renamed G; see docs/design_labels.md).
The G row is computed here from cap_path_powers.csv; the audit had it only from the exp5 note.
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tabutil as T  # noqa: E402
from tabutil import C  # noqa: E402

summ = {r["run"]: r for r in T.rd(os.path.join(C.EXO, "cap_reopt_summary.csv"))}
paths = {r["label"]: r for r in T.rd(os.path.join(C.EXO, "cap_path_powers.csv"))}
ref = float(summ["X0"]["dB_exact_vs_uniform"])

ORDER = [("S0", "S0", "sinc", "uniform"), ("SG", "SR", "sinc", "design G"), ("SX", "SX", "sinc", "end of X0"),
         ("X0", "X0", "exact", "uniform"), ("XG", "XR", "exact", "design G"), ("XS", "XS", "exact", "end of S0")]
rows = []
for paper, data, op, start in ORDER:
    r = summ[data]
    dls, dlx = float(r["dB_sinc_vs_uniform"]), float(r["dB_exact_vs_uniform"])
    rows.append([paper, op, start, int(r["niter"]), dls, dlx, float(r["ratio"]), float(r["pred_minus_true_dB"]), dlx - ref])

u, g = paths["paper_uniform"], paths["paper_quadvel-rad"]
dls = 10 * math.log10(float(g["P_sinc_consistent"]) / float(u["P_sinc_consistent"]))
dlx = 10 * math.log10(float(g["P_exact"]) / float(u["P_exact"]))
rows.append(["G", "---", "(earlier)", None, dls, dlx, float(g["ratio_sinc_over_exact"]), dls - dlx, dlx - ref])

ident = lambda v: str(v)
fmt = [ident, ident, ident, lambda v: "---" if v is None else str(v), T.f(2), T.f(2), T.f(3), T.f(2, True), T.f(2)]
T.write_table("tab_opt", ["run", "operator", "start", "iterations", "dL_sinc", "dL_2s", "rho", "e_red", "s"], rows, fmt)

typed = T.manuscript_table("tab:opt")
computed = [[r[3] if r[3] is not None else None, r[4], r[5], r[6], r[7], r[8]] for r in rows]
n, bad = T.report("tab:opt", computed, typed, [3, 4, 5, 6, 7, 8])
print("uniform design: rho = %.3f (caption says 0.636)" % float(u["ratio_sinc_over_exact"]))
sys.exit(1 if bad else 0)
