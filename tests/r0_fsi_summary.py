"""Combine per-cell FSI results (results/r0_revision/fsi) into summary CSVs, markdown tables and figures.
Run: python -u tests/r0_fsi_summary.py
"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "2"
import csv, math, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(HERE), "results", "r0_revision", "fsi")
CELLS = ["c30_k1", "c30_k3", "c90_k1p5", "c90_k7", "c120_k1p5", "c120_k7", "c30_k18"]
PLOT_ORDER = ["c30_k1", "c90_k1p5", "c120_k1p5", "c30_k3", "c90_k7", "c120_k7", "c30_k18"]
DES = ["uniform", "XU", "XS", "SU", "SX", "EU", "EX"]
VAR = ["vac", "dev", "phys"]
CP = 6283.19 / 343.0  # s_phys = CP / kR


def rd(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def fl(x):
    try:
        return float(x)
    except Exception:
        return float("nan")


def db(x):
    return 10 * math.log10(x) if x > 0 else float("nan")


def label(c):
    th, k = c[1:].split("_")
    k = k[1:].replace("p", ".")
    return f"{th.replace('', '')} deg, kR {k}" if False else f"{c[1:].split('_')[0]} deg, kR {k}"


# ---------------- load
rows, sweeps, modes = {}, {}, {}
for c in CELLS:
    rows[c] = {(r["design"], r["variant"]): r for r in rd(os.path.join(D, f"fsi_{c}.csv"))}
    sweeps[c] = rd(os.path.join(D, f"fsi_sweep_{c}.csv"))
    modes[c] = rd(os.path.join(D, f"fsi_modes_{c}.csv"))

P = {}  # P[c][v][d]
for c in CELLS:
    P[c] = {v: {d: fl(rows[c][(d, v)]["P_rad_Fx"]) for d in DES} for v in VAR}

# ---------------- main summary
out = []
rank, regret = {}, {}
for c in CELLS:
    for v in VAR:
        order = sorted(DES, key=lambda d: P[c][v][d])
        mn = P[c][v][order[0]]
        for d in DES:
            rank[(c, v, d)] = order.index(d) + 1
            regret[(c, v, d)] = db(P[c][v][d] / mn)
hdr = ["cell", "theta0_deg", "kR", "design", "P_vac", "P_dev", "P_phys", "ratio_dev", "ratio_phys",
       "dB_dev_over_vac", "dB_phys_over_vac", "ratio_dev_L8000", "sigma_vac", "sigma_dev", "sigma_phys",
       "dB_vs_uniform_vac", "dB_vs_uniform_dev", "dB_vs_uniform_phys", "rank_vac", "rank_dev", "rank_phys",
       "regret_vac_dB", "regret_dev_dB", "regret_phys_dB", "Prad_over_Pin_dev", "balance_resid_dev",
       "KEa_over_KE_dev", "KEa_over_KE_phys", "Hx_over_Fx_vac", "vel_frac_in_fluid_basis"]
main = {}
for c in CELLS:
    for d in DES:
        r = {v: rows[c][(d, v)] for v in VAR}
        rl = rows[c][(d, "dev_L8000")]
        o = dict(cell=c, theta0_deg=fl(r["vac"]["theta0_deg"]), kR=fl(r["vac"]["kR"]), design=d,
                 P_vac=P[c]["vac"][d], P_dev=P[c]["dev"][d], P_phys=P[c]["phys"][d],
                 ratio_dev=fl(r["dev"]["ratio_fsi_over_vac"]), ratio_phys=fl(r["phys"]["ratio_fsi_over_vac"]),
                 dB_dev_over_vac=fl(r["dev"]["dB_fsi_over_vac"]), dB_phys_over_vac=fl(r["phys"]["dB_fsi_over_vac"]),
                 ratio_dev_L8000=fl(rl["ratio_fsi_over_vac"]),
                 sigma_vac=fl(r["vac"]["sigma"]), sigma_dev=fl(r["dev"]["sigma"]), sigma_phys=fl(r["phys"]["sigma"]),
                 Prad_over_Pin_dev=P[c]["dev"][d] / fl(r["dev"]["P_in"]),
                 balance_resid_dev=fl(r["dev"]["balance_resid_rel_Pin"]),
                 KEa_over_KE_dev=fl(r["dev"]["KE_added_over_KE"]), KEa_over_KE_phys=fl(r["phys"]["KE_added_over_KE"]),
                 Hx_over_Fx_vac=fl(r["vac"]["P_rad_Hx_unclipped"]) / P[c]["vac"][d],
                 vel_frac_in_fluid_basis=fl(r["vac"]["vel_frac_in_fluid_basis"]))
        for v in VAR:
            o[f"dB_vs_uniform_{v}"] = fl(r[v]["dB_vs_uniform_same_variant"])
            o[f"rank_{v}"] = rank[(c, v, d)]
            o[f"regret_{v}_dB"] = regret[(c, v, d)]
        main[(c, d)] = o
        out.append(o)
with open(os.path.join(D, "fsi_summary.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=hdr); w.writeheader(); w.writerows(out)

# ---------------- verdict
def best(c, v, ds):
    d = min(ds, key=lambda x: P[c][v][x]); return d, P[c][v][d]


vh = ["cell", "theta0_deg", "kR"]
for v in VAR:
    vh += [f"best_exact_{v}", f"best_sinc_{v}", f"best_erp_{v}", f"sinc_regret_{v}_dB", f"sinc_regret_fromU_{v}_dB",
           f"exact_depth_vs_uniform_{v}_dB", f"sinc_depth_vs_uniform_{v}_dB", f"erp_depth_vs_uniform_{v}_dB",
           f"quietest_of_7_{v}"]
vh += ["vac_exact_best_still_quietest_dev", "vac_exact_best_still_quietest_phys"]
verd = []
V = {}
for c in CELLS:
    o = dict(cell=c, theta0_deg=main[(c, "uniform")]["theta0_deg"], kR=main[(c, "uniform")]["kR"])
    for v in VAR:
        be, pe = best(c, v, ["XU", "XS"]); bs, ps = best(c, v, ["SU", "SX"]); bp, pp = best(c, v, ["EU", "EX"])
        pu = P[c][v]["uniform"]
        o[f"best_exact_{v}"], o[f"best_sinc_{v}"], o[f"best_erp_{v}"] = be, bs, bp
        o[f"sinc_regret_{v}_dB"] = db(ps / pe)
        o[f"sinc_regret_fromU_{v}_dB"] = db(P[c][v]["SU"] / P[c][v]["XU"])
        o[f"exact_depth_vs_uniform_{v}_dB"] = db(pe / pu)
        o[f"sinc_depth_vs_uniform_{v}_dB"] = db(ps / pu)
        o[f"erp_depth_vs_uniform_{v}_dB"] = db(pp / pu)
        o[f"quietest_of_7_{v}"] = min(DES, key=lambda x: P[c][v][x])
    for v in ("dev", "phys"):
        o[f"vac_exact_best_still_quietest_{v}"] = (o[f"quietest_of_7_{v}"] == o["best_exact_vac"])
    verd.append(o); V[c] = o
with open(os.path.join(D, "fsi_cell_verdict.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=vh); w.writeheader(); w.writerows(verd)

# ---------------- modes
mh = ["cell", "design", "n_modes", "f1_vac_Hz", "shift_dev_pct_median", "shift_dev_pct_min", "shift_dev_pct_max",
      "shift_phys_pct_median", "shift_phys_pct_min", "shift_phys_pct_max"]
mrows = []
for c in CELLS:
    for d in DES:
        m = [r for r in modes[c] if r["design"] == d]
        if not m:
            continue
        sd = np.array([fl(r["shift_dev_pct"]) for r in m]); sp_ = np.array([fl(r["shift_phys_pct"]) for r in m])
        f1 = min(fl(r["f_vac"]) for r in m if int(float(r["mode"])) == 1)
        mrows.append(dict(cell=c, design=d, n_modes=len(m), f1_vac_Hz=f1,
                          shift_dev_pct_median=np.median(sd), shift_dev_pct_min=sd.min(), shift_dev_pct_max=sd.max(),
                          shift_phys_pct_median=np.median(sp_), shift_phys_pct_min=sp_.min(), shift_phys_pct_max=sp_.max()))
with open(os.path.join(D, "fsi_modes_summary.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=mh); w.writeheader(); w.writerows(mrows)

# ---------------- log tables
L = []
def g(x, n=3):
    return f"{x:.{n}g}"
def table(head, body):
    L.append("| " + " | ".join(head) + " |"); L.append("|" + "|".join("---" for _ in head) + "|")
    for b in body:
        L.append("| " + " | ".join(b) + " |")
    L.append("")

L.append("## Main table (P in the pilot's units; dB vs uniform in the same variant; negative = quieter)")
body = []
for c in CELLS:
    for d in DES:
        o = main[(c, d)]
        body.append([c, d, g(o["P_vac"]), g(o["ratio_dev"], 4), g(o["ratio_phys"], 4),
                     f"{o['dB_vs_uniform_vac']:.2f}", f"{o['dB_vs_uniform_dev']:.2f}", f"{o['dB_vs_uniform_phys']:.2f}",
                     g(o["sigma_vac"]), g(o["sigma_dev"]), f"{o['rank_vac']}/{o['rank_dev']}/{o['rank_phys']}"])
table(["cell", "design", "P_vac", "ratio_dev", "ratio_phys", "dB vs U vac", "dB vs U dev", "dB vs U phys",
       "sigma vac", "sigma dev", "rank vac/dev/phys"], body)

L.append("## Verdict per cell")
body = []
for c in CELLS:
    o = V[c]
    body.append([c, f"{o['best_exact_vac']}/{o['best_exact_dev']}/{o['best_exact_phys']}",
                 f"{o['best_sinc_vac']}/{o['best_sinc_dev']}/{o['best_sinc_phys']}",
                 " / ".join(f"{o[f'sinc_regret_{v}_dB']:.2f}" for v in VAR),
                 " / ".join(f"{o[f'sinc_regret_fromU_{v}_dB']:.2f}" for v in VAR),
                 " / ".join(f"{o[f'exact_depth_vs_uniform_{v}_dB']:.2f}" for v in VAR),
                 f"{o['quietest_of_7_vac']}/{o['quietest_of_7_dev']}/{o['quietest_of_7_phys']}",
                 f"{o['vac_exact_best_still_quietest_dev']}/{o['vac_exact_best_still_quietest_phys']}"])
table(["cell", "best exact v/d/p", "best sinc v/d/p", "sinc regret dB v/d/p", "sinc regret from U dB v/d/p",
       "exact depth dB v/d/p", "quietest of 7 v/d/p", "vac-exact still quietest dev/phys"], body)

L.append("## Modes summary (shift in %, 40 modes)")
table(["cell", "design", "f1 vac Hz", "dev med [min,max]", "phys med [min,max]"],
      [[r["cell"], r["design"], f"{r['f1_vac_Hz']:.2f}",
        f"{r['shift_dev_pct_median']:.3f} [{r['shift_dev_pct_min']:.3f}, {r['shift_dev_pct_max']:.3f}]",
        f"{r['shift_phys_pct_median']:.3f} [{r['shift_phys_pct_min']:.3f}, {r['shift_phys_pct_max']:.3f}]"] for r in mrows])

L.append("## Diagnostics")
body = []
for c in CELLS:
    for d in DES:
        o = main[(c, d)]
        body.append([c, d, g(o["Prad_over_Pin_dev"]), g(o["balance_resid_dev"]), g(o["KEa_over_KE_dev"]),
                     g(o["KEa_over_KE_phys"]), g(o["ratio_dev_L8000"], 5), g(o["Hx_over_Fx_vac"], 4), g(o["vel_frac_in_fluid_basis"], 6)])
table(["cell", "design", "Prad/Pin dev", "balance resid dev", "KEa/KE dev", "KEa/KE phys", "ratio dev L8000",
       "Hx/Fx vac", "vel frac"], body)
txt = "\n".join(L)
open(os.path.join(D, "fsi_summary.log"), "w", encoding="utf-8").write(txt)
print(txt)

# ---------------- figures
plt.rcParams.update({"font.size": 12, "axes.labelsize": 13, "axes.titlesize": 12, "legend.fontsize": 11,
                     "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white"})
COL = {"uniform": "#7f7f7f", "XU": "#0072B2", "XS": "#56B4E9", "SU": "#D55E00", "SX": "#E69F00",
       "EU": "#009E73", "EX": "#44BBA4"}
LS = {"uniform": "-", "XU": "-", "XS": "--", "SU": "-", "SX": "--"}


def save(fig, name):
    for ext, kw in (("png", dict(dpi=200)), ("pdf", {})):
        fig.savefig(os.path.join(D, f"{name}.{ext}"), facecolor="white", bbox_inches="tight", **kw)
    plt.close(fig)


# (a) depth
fig, ax = plt.subplots(figsize=(11, 5.2))
x0 = np.arange(len(PLOT_ORDER))
fam = [("exact", ["XU", "XS"], "#0072B2", -0.22, "Exact (best of XU, XS)"),
       ("sinc", ["SU", "SX"], "#D55E00", 0.0, "Sinc (best of SU, SX)"),
       ("ERP", ["EU", "EX"], "#009E73", 0.22, "ERP (best of EU, EX)")]
mk = {"vac": ("o", True), "dev": ("o", False), "phys": ("x", None)}
for nm, ds, col, off, lab in fam:
    for v in VAR:
        y = [-db(min(P[c][v][d] for d in ds) / P[c][v]["uniform"]) * -1 for c in PLOT_ORDER]
        y = [-db(P[c][v]["uniform"] / min(P[c][v][d] for d in ds)) * -1 for c in PLOT_ORDER]
        y = [db(P[c][v]["uniform"] / min(P[c][v][d] for d in ds)) for c in PLOT_ORDER]
        m, filled = mk[v]
        kw = dict(marker=m, ms=9, color=col, ls="none")
        if filled is True:
            kw.update(mfc=col)
        elif filled is False:
            kw.update(mfc="white", mew=1.8)
        else:
            kw.update(mew=2.2, ms=10)
        ax.plot(x0 + off, y, **kw)
    for i, c in enumerate(PLOT_ORDER):  # connect vac->phys per cell for the family
        ys = [db(P[c][v]["uniform"] / min(P[c][v][d] for d in ds)) for v in VAR]
        ax.plot([i + off] * 3, ys, "-", color=col, lw=1, alpha=0.5)
ax.set_xticks(x0); ax.set_xticklabels([label(c) for c in PLOT_ORDER], rotation=20, ha="right")
ax.set_ylabel("Reduction vs uniform (dB)"); ax.grid(axis="y", alpha=0.3); ax.axhline(0, color="k", lw=0.6)
from matplotlib.lines import Line2D
h = [Line2D([], [], marker="s", ls="none", color=c, ms=9, label=l) for _, _, c, _, l in fam]
h += [Line2D([], [], marker="o", ls="none", color="k", ms=9, label="in vacuo"),
      Line2D([], [], marker="o", ls="none", color="k", mfc="white", mew=1.8, ms=9, label="device (fluid)"),
      Line2D([], [], marker="x", ls="none", color="k", mew=2.2, ms=10, label="physical density")]
ax.legend(handles=h, ncol=2, loc="upper left", frameon=False)
ax.set_ylim(top=ax.get_ylim()[1] * 1.45)
save(fig, "fsi_depth")

# (b) ratio sweep
fig, axs = plt.subplots(2, 4, figsize=(17, 8))
axs = axs.ravel()
sd_ = ["uniform", "XU", "XS", "SU", "SX"]
for i, c in enumerate(PLOT_ORDER):
    ax = axs[i]; kR = main[(c, "uniform")]["kR"]
    for d in sd_:
        rr = [r for r in sweeps[c] if r["design"] == d]
        s = np.array([fl(r["s"]) for r in rr]); y = np.array([fl(r["ratio_fsi_over_vac"]) for r in rr])
        o = np.argsort(s)
        ax.plot(s[o], y[o], LS[d], color=COL[d], lw=2, marker="o", ms=3.5, label=d)
    ax.axvline(1.0, color="k", ls="--", lw=1)
    ax.axvline(CP / kR, color="k", ls=":", lw=1.5)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_title(label(c))
    ax.set_xlabel("Fluid-loading scale s"); ax.grid(alpha=0.3, which="both", lw=0.4)
    if i % 4 == 0:
        ax.set_ylabel(r"$P_{FSI}/P_{vac}$")
for ax in axs[7:]:
    ax.axis("off")
h = [Line2D([], [], color=COL[d], ls=LS[d], lw=2, marker="o", ms=4, label=d) for d in sd_]
h += [Line2D([], [], color="k", ls="--", lw=1, label="s = 1 (device)"),
      Line2D([], [], color="k", ls=":", lw=1.5, label=r"s = c'/c$_0$ (physical density)")]
axs[7].legend(handles=h, loc="center", frameon=False, fontsize=13)
fig.tight_layout()
save(fig, "fsi_ratio_sweep")

# (c) regret
fig, ax = plt.subplots(figsize=(10, 4.8))
for v, (m, filled), off, lab in zip(VAR, [mk[v] for v in VAR], (-0.15, 0, 0.15),
                                    ("in vacuo", "device (fluid)", "physical density")):
    y = [V[c][f"sinc_regret_{v}_dB"] for c in PLOT_ORDER]
    kw = dict(marker=m, ms=10, color="#D55E00", ls="none", label=lab)
    if filled is True:
        kw.update(mfc="#D55E00")
    elif filled is False:
        kw.update(mfc="white", mew=1.8)
    else:
        kw.update(mew=2.2)
    ax.plot(x0 + off, y, **kw)
ax.set_xticks(x0); ax.set_xticklabels([label(c) for c in PLOT_ORDER], rotation=20, ha="right")
ax.set_ylabel("Sinc regret: best sinc vs best exact (dB)"); ax.grid(axis="y", alpha=0.3); ax.axhline(0, color="k", lw=0.6)
ax.legend(frameon=False)
save(fig, "fsi_regret")
print("DONE")

