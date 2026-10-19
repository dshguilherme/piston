"""Test of cheap corrections to the sinc model, and of a cancellation index as failure predictor.

Data (all produced by the r0_nov_* scripts):
  cap_frozen_sweep.csv  six frozen cap fields, R = 10 m, theta0 = 30 deg, kR = 0.1..40 (spectral two-sided reference)
  cap_modes_ratio.csv   analytic cap modes (uniform, cos0..cos7) on caps theta0 = 15..150 deg
Corrections (NOT fitted to the frozen fields):
  C1  two-face factor:      P_corr = 2 P_sinc
  C2  canonical-cap factor: P_corr = P_sinc / r_u(kR), r_u = sinc/two-sided ratio of the UNIFORM-velocity cap of the same
                            geometry (one canonical exact solve per frequency; no information from the field itself)
Error measure: e = P_corr / P_2s - 1.   Regimes: kR < 4, 4 <= kR < 10, kR >= 10 (cap: k a = kR sin(theta0)).
Cancellation index: sinc radiation efficiency sigma_s = P_sinc / (0.5 rho c int |v|^2 dS); Spearman correlation between
|log(2 P_sinc / P_2s)| and log sigma_s.
Outputs: results/r0_revision/novelty/correction_summary.csv, correction.log
"""
import os, sys, csv
import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(os.path.dirname(HERE), "results", "r0_revision", "novelty")
RHO0, C0 = 1.21, 343.0


def rd(name):
    return list(csv.DictReader(open(os.path.join(OUT, name))))


def v2_integrals():
    """int |v|^2 dS for the six frozen fields on the m8 trace (6-point rule, P1)."""
    import r0_nov_cap_frozen as cf
    V, T, vn, labels = cf.load_mesh(8)
    th, ph, w, vq = cf.quad_points(V, T, vn)
    return dict(zip(labels, (w[:, None] * np.abs(vq) ** 2).sum(0)))


def stats(e):
    e = np.asarray(e)
    return dict(n=len(e), median=np.median(e), p05=np.percentile(e, 5), p95=np.percentile(e, 95),
                max_abs=np.abs(e).max(), frac_within_20pct=np.mean(np.abs(e) <= 0.2))


def regime(kR):
    return "kR<4" if kR < 4 else ("4<=kR<10" if kR < 10 else "kR>=10")


def main():
    log = []
    pr = lambda s: (print(s, flush=True), log.append(s))
    modes = rd("cap_modes_ratio.csv")
    ru = {(int(m["theta0_deg"]), round(float(m["kR"]), 6)): float(m["ratio"]) for m in modes if m["field"] == "uniform"}
    rows = []
    # frozen fields
    fr = rd("cap_frozen_sweep.csv")
    v2 = v2_integrals()
    data = {}
    for r in fr:
        kR = float(r["kR"]); key = (30, round(kR, 6))
        ratio = float(r["ratio"])
        if key not in ru:
            continue
        sig = float(r["P_sinc"]) / (0.5 * RHO0 * C0 * v2[r["design"]])
        data.setdefault(regime(kR), []).append((ratio, ru[key], sig))
    pr("FROZEN CAP FIELDS (6 fields, theta0=30 deg)")
    for reg in ("kR<4", "4<=kR<10", "kR>=10"):
        d = np.array(data[reg])
        for name, e in (("none", d[:, 0] - 1), ("C1 x2", 2 * d[:, 0] - 1), ("C2 /r_u", d[:, 0] / d[:, 1] - 1)):
            s = stats(e); rows.append(dict(dataset="frozen_cap", regime=reg, correction=name, **s))
            pr(f"  {reg:9s} {name:8s} n={s['n']:3d} median {s['median']:+.3f}  5-95% [{s['p05']:+.3f},{s['p95']:+.3f}]"
               f"  max|e| {s['max_abs']:.3f}  within 20%: {s['frac_within_20pct']:.2f}")
    # modes
    pr("ANALYTIC CAP MODES (uniform excluded from C2 by construction)")
    md = {}
    for m in modes:
        th, kR = int(m["theta0_deg"]), float(m["kR"])
        if m["field"] == "uniform":
            continue
        md.setdefault((th, regime(kR)), []).append((float(m["ratio"]), ru[(th, round(kR, 6))], float(m["sigma_sinc"]),
                                                     m["field"]))
    for th in (15, 30, 60, 90, 120, 150):
        for reg in ("kR<4", "4<=kR<10", "kR>=10"):
            d = np.array([x[:3] for x in md[(th, reg)]], float)
            for name, e in (("C1 x2", 2 * d[:, 0] - 1), ("C2 /r_u", d[:, 0] / d[:, 1] - 1)):
                s = stats(e); rows.append(dict(dataset=f"modes_theta{th}", regime=reg, correction=name, **s))
            s1 = stats(2 * d[:, 0] - 1); s2 = stats(d[:, 0] / d[:, 1] - 1)
            pr(f"  th0={th:3d} {reg:9s} C1: median {s1['median']:+.3f} [{s1['p05']:+.3f},{s1['p95']:+.3f}] w20 {s1['frac_within_20pct']:.2f}"
               f" | C2: median {s2['median']:+.3f} [{s2['p05']:+.3f},{s2['p95']:+.3f}] w20 {s2['frac_within_20pct']:.2f}")
    # cancellation index as predictor (kR >= 4, where the dipole regime is excluded)
    pr("CANCELLATION INDEX: Spearman( |log(2 P_sinc/P_2s)| , log sigma_sinc ), kR >= 4")
    for label, pts in (("frozen", [x for reg in ("4<=kR<10", "kR>=10") for x in data[reg]]),
                       ("modes", [x[:3] for k, v in md.items() if k[1] != "kR<4" for x in v])):
        d = np.array(pts, float)
        rs = spearmanr(np.abs(np.log(2 * d[:, 0])), np.log(d[:, 2]))
        pr(f"  {label:7s} n={len(d)} rho_s={rs.statistic:+.3f} p={rs.pvalue:.2e}")
    with open(os.path.join(OUT, "correction_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    with open(os.path.join(OUT, "correction.log"), "w") as f:
        f.write("\n".join(log) + "\n")


if __name__ == "__main__":
    main()
