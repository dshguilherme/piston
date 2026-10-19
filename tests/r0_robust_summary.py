"""Experiment 10: summary statistics of the fragility evaluation (r0_robust_evalstudy.m output). numpy only.

For every cell, design, sigma and power column (exact, exact L=8000, sinc, exact with fluid loading):
  reference  P_ref = power of the NOMINAL uniform design (sigma = 0) under the same column (operator / medium)
  per sample dB_i = 10 log10(P_i / P_ref)  (negative = quieter than the nominal uniform design)
  reported   nominal dB (sigma = 0), dB of the sample-mean power (the expected-power objective),
             mean, median, 5th and 95th percentile of dB_i, the loss of the median against nominal,
             the relative standard deviation of the power and the relative standard error of the
             sample-mean power for N = 8, 16, 32 (std / mean / sqrt(N)).
Run: python tests/r0_robust_summary.py <evalstudy csv> [<out csv>]
"""
import os, sys, csv
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "r0_revision", "robust")
FIELDS = ["cell", "design", "sigma", "column", "n", "nominal_dB", "meanpower_dB", "mean_dB", "median_dB", "p5_dB",
          "p95_dB", "median_loss_dB", "rel_sd_power", "se_rel_N8", "se_rel_N16", "se_rel_N32"]


def load(fn):
    with open(fn) as f:
        r = csv.reader(f)
        head = next(r)
        rows = [x for x in r]
    cols = {h: [x[i] for x in rows] for i, h in enumerate(head)}
    for h in head:
        if h.startswith("P_") or h in ("sigma", "vol"):
            cols[h] = np.array([float(v) for v in cols[h]])
    cols["design"] = np.array(cols["design"]); cols["cell"] = np.array(cols["cell"])
    return head, cols


def summarise(fn):
    head, d = load(fn)
    pcols = [h for h in head if h.startswith("P_")]
    out = []
    for cell in dict.fromkeys(d["cell"]):
        mc = d["cell"] == cell
        iref = np.where(mc & (d["design"] == "uniform") & (d["sigma"] == 0))[0][0]
        for des in dict.fromkeys(d["design"][mc]):
            md = mc & (d["design"] == des)
            inom = np.where(md & (d["sigma"] == 0))[0][0]
            for sg in dict.fromkeys(d["sigma"][md]):
                m = md & (d["sigma"] == sg)
                for c in pcols:
                    p = d[c][m]
                    if np.all(np.isnan(p)):
                        continue
                    ref = d[c][iref]; nom = d[c][inom]
                    db = 10 * np.log10(p / ref)
                    rsd = p.std(ddof=1) / p.mean() if len(p) > 1 else 0.0
                    nd = 10 * np.log10(nom / ref)
                    out.append(dict(cell=cell, design=des, sigma=sg, column=c, n=len(p), nominal_dB=nd,
                                    meanpower_dB=10 * np.log10(p.mean() / ref), mean_dB=db.mean(),
                                    median_dB=np.median(db), p5_dB=np.percentile(db, 5), p95_dB=np.percentile(db, 95),
                                    median_loss_dB=np.median(db) - nd, rel_sd_power=rsd, se_rel_N8=rsd / np.sqrt(8),
                                    se_rel_N16=rsd / np.sqrt(16), se_rel_N32=rsd / np.sqrt(32)))
    return out


def main(fn, outfn=None):
    S = summarise(fn)
    if outfn is None:
        outfn = os.path.join(OUT, os.path.basename(fn).replace("evalstudy_", "evalstudy_summary_"))
    with open(outfn, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        for r in S:
            w.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()})
    for c in dict.fromkeys(r["column"] for r in S):
        print(f"\n== {c}: design sigma | nominal | mean-power | median | p5 | p95 | median loss | rel sd P")
        for r in S:
            if r["column"] == c and r["sigma"] > 0:
                print(f"  {r['design']:>32s} {r['sigma']:6.3f} | {r['nominal_dB']:8.2f} | {r['meanpower_dB']:8.2f} | "
                      f"{r['median_dB']:8.2f} | {r['p5_dB']:8.2f} | {r['p95_dB']:8.2f} | {r['median_loss_dB']:7.2f} | {r['rel_sd_power']:.3g}")
    print(f"-> {outfn}")


if __name__ == "__main__":
    main(*sys.argv[1:])
