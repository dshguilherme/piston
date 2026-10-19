"""Compact per-cell x rung summary of the R0 regime pilot -> results/r0_revision/pilot/pilot_cells_summary.csv.
Joins pilot_rung_summary.csv, pilot_fidelity.csv and pilot_reference.csv (written by tests/r0_pilot_stats.py).
Regrets are LOWER bounds on suboptimality (relative to the best exact design found in the cell, not a proven optimum).
"""
import os, csv

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "r0_revision", "pilot")
rd = lambda n: list(csv.DictReader(open(os.path.join(OUT, n))))
summ, fid, ref = rd("pilot_rung_summary.csv"), rd("pilot_fidelity.csv"), {r["cell"]: r for r in rd("pilot_reference.csv")}
F = {(r["cell"], r["rung"], r["set"]): r for r in fid}
rows = []
for s in summ:
    c, g = s["cell"], s["rung"]
    f = lambda st, k: F.get((c, g, st), {}).get(k, "")
    rows.append(dict(cell=c, theta0_deg=s["theta0_deg"], kR=s["kR"], rung=g,
                     optimisation_runs="shared with R1 (R2 = 2 x R1, same argmin)" if g == "R2_sinc_x2" else s["n_runs"],
                     pick_run=s["pick_run"], pick_stop=s["pick_stop"], any_run_at_iteration_cap=s["any_unconverged"],
                     pick_pred_dB_vs_uniform=s["pick_pred_dB"], pick_true_dB_vs_uniform=s["pick_true_dB"],
                     pick_pred_minus_true_dB=s["pick_pred_minus_true_dB"],
                     pick_regret_lowerbound_dB=s["pick_regret_lowerbound_dB"],
                     fromU_regret_lowerbound_dB=s["fromU_regret_lowerbound_dB"],
                     best_exact_dB_vs_uniform=s["best_exact_dB_vs_uniform"], pick_sigma_exact=s["pick_sigma"],
                     ratio_uniform=f("uniform+random", "ratio_uniform"),
                     ratio_median_random=f("uniform+random", "ratio_median"), tau_random=f("uniform+random", "kendall_tau"),
                     ratio_median_all=f("all", "ratio_median"), ratio_p05_all=f("all", "ratio_p05"),
                     ratio_p95_all=f("all", "ratio_p95"), tau_all=f("all", "kendall_tau"),
                     tau_within3dB_best_exact=f("within3dB_best_exact", "kendall_tau"),
                     ref_eL_max_cell=ref[c]["eL_max"], ref_eN_max_cell=ref[c]["eN_max"]))
with open(os.path.join(OUT, "pilot_cells_summary.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(f"wrote pilot_cells_summary.csv: {len(rows)} rows")
