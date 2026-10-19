"""Experiment 8 summary tables from results/r0_revision/longruns (hist_*.csv) and the pilot.

Writes
  longruns/long_runs_A.csv      one row per part-A run and checkpoint iteration N (600, 1200, ...):
                                pick = lowest own-objective iterate up to N; its exact and sinc dB
  longruns/long_cells_A.csv     one row per cell and N: XU/XS/SU/SX exact dB of the pick, best exact
                                found in the cell up to N (any run of pilot + part A), exact spread
                                between starts, sinc regret LB and sinc prediction error
  longruns/long_saturation_A.csv per part-A run: own gain per 100 iterations in each 600 block, and
                                the iteration at which the best-so-far came within 1 and 0.1 dB of
                                its final value
  longruns/long_strategies_B.csv per part-B job: best exact dB at 150/300/.../1200 iterations, wall
                                time, and the iterations / wall time needed to reach reference levels
All dB values are 10 log10(P / P(uniform)) of the same rung (exact unless stated).
"""
import glob
import os
import re

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LR = os.path.join(ROOT, "results", "r0_revision", "longruns")
PILOT = os.path.join(ROOT, "results", "r0_revision", "pilot")
CHECKS = [600, 1200, 1800, 2400, 3000, 3600, 4200, 4800, 5000]
# Feasibility: volume fraction <= 0.5 with a 1e-5 relative tolerance (the pilot accepted 0.500001).
# GCMMA treats the volume constraint softly (penalty c = 1000), and some runs drift above it; every
# "pick" and every "best" below is taken over feasible iterates only.
VOL_MAX = 0.5 * (1 + 1e-5)


def load_hist(kind):
    out = {}
    for fn in sorted(glob.glob(os.path.join(LR, f"hist_{kind}_*.csv"))):
        m = re.match(rf"hist_{kind}_(c\d+_k[0-9p]+)_(\w+)\.csv", os.path.basename(fn))
        if not m:
            continue
        df = pd.read_csv(fn)
        if len(df):
            feas = df["vol"] <= VOL_MAX
            df["feasible"] = feas
            df["own_f"] = df["own_dB"].where(feas, np.inf)
            df["exact_f"] = df["exact_dB"].where(feas, np.inf)
            out[(m.group(1), m.group(2))] = df
    return out


def pilot_pick(cell, run):
    ev = pd.read_csv(os.path.join(PILOT, f"eval_{cell}.csv"))
    u = ev[ev["set"] == "uniform"].iloc[0]
    o = ev[(ev["set"] == "optimum") & (ev["run"] == run)]
    if not len(o):
        return None
    o = o.iloc[0]
    db = lambda col: 10 * np.log10(o[col] / u[col])
    return dict(exact_dB=db("P_exact"), sinc_dB=db("P_sinc"), erp_dB=db("P_erp"), vol=o["vol"],
                exact_L8000_dB=db("P_exact_L8000"), exact_N32_dB=db("P_exact_N32"))


def pilot_best_exact(cell):
    ev = pd.read_csv(os.path.join(PILOT, f"eval_{cell}.csv"))
    u = ev[ev["set"] == "uniform"].iloc[0]
    return 10 * np.log10(ev.loc[ev["vol"] <= VOL_MAX, "P_exact"].min() / u["P_exact"])


def run_status(kind):
    fn = os.path.join(LR, f"chunks_{kind}.csv")
    if not os.path.exists(fn):
        return {}
    ch = pd.read_csv(fn)
    return {(r["cell"], r["name"]): r["status"] for _, r in ch.iterrows()}


def summary_A():
    H = load_hist("A")
    status = run_status("A")
    rows, sat = [], []
    for (cell, run), df in H.items():
        own = df["own_f"].to_numpy()
        it = df["it"].to_numpy()
        p = pilot_pick(cell, run)
        if p is not None:
            rows.append(dict(cell=cell, run=run, rung=run[0], N=int(it[0]) - 1, source="pilot",
                             pick_iter=int(it[0]) - 1, own_dB=np.nan, exact_dB=p["exact_dB"],
                             sinc_dB=p["sinc_dB"], vol=p["vol"]))
        conv = status.get((cell, run)) == "converged"
        horizon = max(int(d["it"].iloc[-1]) for (c2, _), d in H.items() if c2 == cell)
        for N in CHECKS:
            if N <= it[0] - 1 or (N > it[-1] and not conv) or N > max(horizon, it[-1]) + 599:
                continue
            k = np.searchsorted(it, N, side="right")
            j = int(np.argmin(own[:k]))
            src = "long" if N <= it[-1] else f"converged at {int(it[-1])}"
            rows.append(dict(cell=cell, run=run, rung=run[0], N=N, source=src, pick_iter=int(it[j]),
                             own_dB=own[j], exact_dB=df["exact_dB"].iloc[j], sinc_dB=df["sinc_dB"].iloc[j],
                             vol=df["vol"].iloc[j]))
        best = np.minimum.accumulate(own)
        final = best[-1]
        rec = dict(cell=cell, run=run, first_iter=int(it[0]), last_iter=int(it[-1]),
                   status=status.get((cell, run), ""), final_own_dB=final)
        for b0 in range(int(it[0]) - 1, int(it[-1]), 600):
            sel = (it > b0) & (it <= b0 + 600)
            if sel.sum() >= 100:
                bb = best[sel]
                rec[f"gain_per100_{b0}_{b0 + 600}"] = (bb[0] - bb[-1]) / (len(bb) - 1) * 100
        for tol in (1.0, 0.1):
            k = int(np.argmax(best <= final + tol))
            rec[f"iter_within_{tol}dB_of_final"] = int(it[k])
        sat.append(rec)
    A = pd.DataFrame(rows).sort_values(["cell", "run", "N"])
    A.to_csv(os.path.join(LR, "long_runs_A.csv"), index=False, float_format="%.4f")
    pd.DataFrame(sat).to_csv(os.path.join(LR, "long_saturation_A.csv"), index=False, float_format="%.4f")
    crow = []
    HB = load_hist("B")
    for cell in sorted(A["cell"].unique()):
        sub = A[A["cell"] == cell]
        pb = pilot_best_exact(cell)
        best_all = min([pb] + [df["exact_f"].min() for (c2, _), df in list(H.items()) + list(HB.items()) if c2 == cell])
        for N in sorted(sub["N"].unique()):
            r = dict(cell=cell, N=N)
            bestN = pb
            for (c2, run), df in H.items():
                if c2 == cell:
                    bestN = min(bestN, df.loc[df["it"] <= N, "exact_f"].min() if (df["it"] <= N).any() else bestN)
            r["best_exact_found_dB"] = bestN
            for run in ("XU", "XS", "SU", "SX"):
                q = sub[(sub["run"] == run) & (sub["N"] == N)]
                if len(q):
                    r[f"{run}_exact_dB"] = q["exact_dB"].iloc[0]
                    r[f"{run}_own_dB"] = q["own_dB"].iloc[0] if run[0] == "S" else q["exact_dB"].iloc[0]
            if "XU_exact_dB" in r and "XS_exact_dB" in r:
                r["exact_spread_XU_XS_dB"] = abs(r["XU_exact_dB"] - r["XS_exact_dB"])
            r["best_exact_all_runs_final_dB"] = best_all
            if "SU_exact_dB" in r:
                r["SU_regret_LB_dB"] = r["SU_exact_dB"] - bestN
                r["SU_regret_vs_final_best_dB"] = r["SU_exact_dB"] - best_all
                if not np.isnan(r.get("SU_own_dB", np.nan)):
                    r["SU_pred_minus_true_dB"] = r["SU_own_dB"] - r["SU_exact_dB"]
            crow.append(r)
    C = pd.DataFrame(crow)
    C.to_csv(os.path.join(LR, "long_cells_A.csv"), index=False, float_format="%.4f")
    return A, C


def summary_B():
    H = load_hist("B")
    if not H:
        return None
    rows = []
    for (cell, strat), df in H.items():
        bx = np.minimum.accumulate(df["exact_f"].to_numpy())
        it = df["it"].to_numpy()
        wt = np.cumsum(df["time_s"].to_numpy())
        r = dict(cell=cell, strategy=strat, iters=int(it[-1]), walltime_s=wt[-1], final_best_exact_dB=bx[-1],
                 stages=" > ".join(f"{g['rung'].iloc[0]}{len(g)}" for _, g in df.groupby("stage")),
                 max_vol=df["vol"].max(), n_infeasible=int((~df["feasible"]).sum()),
                 first_infeasible_iter=int(df.loc[~df["feasible"], "it"].min()) if (~df["feasible"]).any() else np.nan,
                 best_exact_any_vol_dB=df["exact_dB"].min())
        for N in range(150, 1201, 150):
            if N <= it[-1]:
                r[f"best_exact_dB_at_{N}"] = bx[np.searchsorted(it, N, side="right") - 1]
        rows.append(r)
    B = pd.DataFrame(rows)
    # reference levels: plain exact at 600 and 1200 iterations of the same cell (strategy P)
    for cell in B["cell"].unique():
        if (cell, "P") not in H:
            continue
        dfp = H[(cell, "P")]
        bp = np.minimum.accumulate(dfp["exact_f"].to_numpy())
        for N in (600, 1200):
            level = bp[min(N, len(bp)) - 1]  # P may stop early on the stopping rule: its final value
            for (c2, strat), df in H.items():
                if c2 != cell:
                    continue
                bx = np.minimum.accumulate(df["exact_f"].to_numpy())
                wt = np.cumsum(df["time_s"].to_numpy())
                k = np.argmax(bx <= level + 1e-9) if (bx <= level + 1e-9).any() else None
                m = (B["cell"] == cell) & (B["strategy"] == strat)
                B.loc[m, f"P{N}_level_dB"] = level
                B.loc[m, f"iters_to_P{N}_level"] = int(df["it"].iloc[k]) if k is not None else np.nan
                B.loc[m, f"walltime_to_P{N}_level_s"] = wt[k] if k is not None else np.nan
    B = B.sort_values(["cell", "strategy"])
    B.to_csv(os.path.join(LR, "long_strategies_B.csv"), index=False, float_format="%.4f")
    return B


if __name__ == "__main__":
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    A, C = summary_A()
    print(A.to_string(index=False))
    print(C.to_string(index=False))
    print(pd.read_csv(os.path.join(LR, "long_saturation_A.csv")).to_string(index=False))
    B = summary_B()
    if B is not None:
        print(B.to_string(index=False))
