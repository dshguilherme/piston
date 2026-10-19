"""R0 regime pilot: raw metrics per cell and rung from results/r0_revision/pilot/eval_<cell>.csv and runs_<cell>.csv.

No thresholds are applied here; the script only computes and writes raw metrics.
Rungs: R0 = ERP (P_erp), R0-I = legacy 'I' arm (w^2 ||U||^2, fidelity only), R1 = spline-consistent sinc (P_sinc),
R2 = 2 x R1, R4 = exact (P_exact, Richardson L = 4000/8000).
Outputs (pilot/):
  pilot_fidelity.csv     per cell x rung x design set: n, median/5%/95% of P_rung/P_exact, std of log10 ratio, Kendall tau
  pilot_runs.csv         per optimisation run: stop, tails, predicted and true dB vs uniform, regret, efficiency, ref. errors
  pilot_rung_summary.csv per cell x rung: the rung's pick (lowest own objective over its runs), regret, prediction error
  pilot_reference.csv    per cell: distribution of the reference-accuracy proxies over all evaluated designs
"""
import os, sys, csv
import numpy as np
from scipy.stats import kendalltau, spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "r0_revision", "pilot")
CELLS = ["c30_k1", "c30_k3", "c90_k1p5", "c90_k7", "c120_k1p5", "c120_k7", "c30_k18"]
RUNG_COL = {"R0_ERP": "P_erp", "R0_I": "P_I", "R1_sinc": "P_sinc", "R2_sinc_x2": "P_sinc", "R4_exact": "P_exact"}
RUNG_FAC = {"R0_ERP": 1.0, "R0_I": 1.0, "R1_sinc": 1.0, "R2_sinc_x2": 2.0, "R4_exact": 1.0}
RUN_RUNG = {"X": "R4_exact", "S": "R1_sinc", "E": "R0_ERP"}
OWN_COL = {"X": "P_exact", "S": "P_sinc", "E": "P_erp"}


def read(fn):
    with open(fn) as f:
        return list(csv.DictReader(f))


def db(x):
    return 10 * np.log10(x)


def main():
    fid, runs_out, summ, refrows = [], [], [], []
    for cell in CELLS:
        fe = os.path.join(OUT, f"eval_{cell}.csv")
        if not os.path.exists(fe):
            print(f"{cell}: no eval file, skipped"); continue
        E = read(fe)
        Rr = read(os.path.join(OUT, f"runs_{cell}.csv"))
        th, kR = int(E[0]["theta0_deg"]), float(E[0]["kR"])
        col = lambda c, rows=E: np.array([float(r[c]) for r in rows])
        P = {c: col(c) for c in ["P_erp", "P_I", "P_sinc", "P_sinc_spec", "P_exact", "P_exact_L8000", "P_exact_N32"]}
        dset = np.array([r["set"] for r in E]); run = np.array([r["run"] for r in E]); it = np.array([int(r["iter"]) for r in E])
        iu = np.where(dset == "uniform")[0][0]
        eL = np.abs(P["P_exact"] / P["P_exact_L8000"] - 1); eN = np.abs(P["P_exact_N32"] / P["P_exact"] - 1)
        sig = P["P_exact"] / P["P_erp"]
        es = np.abs(P["P_sinc_spec"] / P["P_sinc"] - 1)
        refrows.append(dict(cell=cell, theta0_deg=th, kR=kR, n=len(E),
                            eL_median=np.median(eL), eL_max=eL.max(), eN_median=np.median(eN), eN_max=eN.max(),
                            sinc_spec_vs_consistent_median=np.median(es), sinc_spec_vs_consistent_max=es.max(),
                            sigma_min=sig.min(), sigma_median=np.median(sig), sigma_uniform=sig[iu],
                            corr_eL_logsigma=spearmanr(eL, np.log10(sig))[0]))
        Pbest = P["P_exact"].min(); ibest = int(np.argmin(P["P_exact"]))
        near3 = db(P["P_exact"] / Pbest) <= 3.0
        sets = {"uniform+random": (dset == "uniform") | (dset == "random"),
                "path": (dset == "path") | (dset == "optimum"),
                "all": np.ones(len(E), bool),
                "within3dB_best_exact": near3}
        for rung, c in RUNG_COL.items():
            ratio = RUNG_FAC[rung] * P[c] / P["P_exact"]
            for sname, m in sets.items():
                if m.sum() < 3:
                    continue
                tau = kendalltau(P[c][m], P["P_exact"][m])[0]
                fid.append(dict(cell=cell, theta0_deg=th, kR=kR, rung=rung, set=sname, n=int(m.sum()),
                                ratio_median=np.median(ratio[m]), ratio_p05=np.percentile(ratio[m], 5),
                                ratio_p95=np.percentile(ratio[m], 95), log10ratio_std=np.std(np.log10(ratio[m])),
                                ratio_uniform=ratio[iu], kendall_tau=tau))
        # runs
        cellruns = []
        for r in Rr:
            code = r["run"]; k = np.where((run == code) & (dset == "optimum"))[0]
            if not k.size:
                continue
            k = k[0]; own = OWN_COL[code[0]]
            pred = db(P[own][k] / P[own][iu]); true = db(P["P_exact"][k] / P["P_exact"][iu])
            vol = float(E[k]["vol"])
            nwarn = 0
            for lf in [f for f in os.listdir(OUT) if f.startswith(f"log_run_{cell}_{code}") and f.endswith(".log")]:
                txt = open(os.path.join(OUT, lf), errors="ignore").read()
                if f"RESULT {cell} {code}" in txt:
                    nwarn = txt.count("close to singular")
            rr = dict(cell=cell, theta0_deg=th, kR=kR, run=code, rung=RUN_RUNG[code[0]], start=r["start"],
                      source=r["source"], niter=int(r["niter"]), stop=r["stop"],
                      feasible=bool(vol <= 0.5 * (1 + 1e-6)), vol=vol, gcmma_singular_warnings=nwarn,
                      P_exact_uniform=P["P_exact"][iu], P_own_uniform=P[own][iu],
                      tail100_dB=float(r["tail100_dB"]),
                      tail50_dB=float(r["tail50_dB"]), walltime_s=float(r["walltime_s"]), final_vol=float(r["final_vol"]),
                      pred_dB_vs_uniform=pred, true_dB_vs_uniform=true, pred_minus_true_dB=pred - true,
                      regret_bound_dB=db(P["P_exact"][k] / Pbest), ratio_sinc_exact=P["P_sinc"][k] / P["P_exact"][k],
                      sigma_exact=sig[k], ref_eL=eL[k], ref_eN=eN[k], P_exact=P["P_exact"][k], P_own=P[own][k])
            runs_out.append(rr); cellruns.append(rr)
        best_run = run[ibest] if dset[ibest] != "random" else "random"
        for rl, rung in RUN_RUNG.items():
            rs = [x for x in cellruns if x["rung"] == rung and x["run"] != "SUn"]
            if not rs:
                continue
            pick = min(rs, key=lambda x: x["P_own"])
            fromU = [x for x in rs if x["start"] == "U"]
            row = dict(cell=cell, theta0_deg=th, kR=kR, rung=rung, n_runs=len(rs), pick_run=pick["run"],
                       pick_regret_lowerbound_dB=pick["regret_bound_dB"], pick_pred_dB=pick["pred_dB_vs_uniform"],
                       pick_true_dB=pick["true_dB_vs_uniform"], pick_pred_minus_true_dB=pick["pred_minus_true_dB"],
                       pick_sigma=pick["sigma_exact"], pick_stop=pick["stop"], pick_tail100_dB=pick["tail100_dB"],
                       fromU_regret_lowerbound_dB=fromU[0]["regret_bound_dB"] if fromU else np.nan,
                       fromU_true_dB=fromU[0]["true_dB_vs_uniform"] if fromU else np.nan,
                       best_exact_dB_vs_uniform=db(Pbest / P["P_exact"][iu]), best_exact_from=f"{dset[ibest]}:{best_run}:{it[ibest]}",
                       any_unconverged=any(x["stop"] != "converged" for x in rs))
            summ.append(row)
            if rung == "R1_sinc":
                r2 = dict(row); r2["rung"] = "R2_sinc_x2"; summ.append(r2)
    for name, rows in (("pilot_fidelity.csv", fid), ("pilot_runs.csv", runs_out), ("pilot_rung_summary.csv", summ),
                       ("pilot_reference.csv", refrows)):
        if not rows:
            continue
        with open(os.path.join(OUT, name), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader()
            for r in rows:
                w.writerow({k: (f"{v:.6g}" if isinstance(v, (float, np.floating)) else v) for k, v in r.items()})
        print(f"wrote {name}: {len(rows)} rows")


if __name__ == "__main__":
    main()
