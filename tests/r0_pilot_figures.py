"""Figures for the radiation-model regime pilot (reads CSVs in results/r0_revision/pilot).

Run from harmonic_shell_course:  python -W ignore tests\\r0_pilot_figures.py
Missing cells are drawn as grey 'n/a'; simply re-run when more data exist.
"""
import csv
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
PILOT = os.path.normpath(os.path.join(HERE, "..", "results", "r0_revision", "pilot"))
CELLS = ["c30_k1", "c30_k3", "c30_k18", "c90_k1p5", "c90_k7", "c120_k1p5", "c120_k7"]
# fallback labels only for cells with no data (parsed from the cell name)
RUNGS = ["R0_ERP", "R1_sinc", "R2_sinc_x2", "R4_exact"]
RUNG_LABEL = {"R0_ERP": "R0 ERP", "R1_sinc": "R1 sinc", "R2_sinc_x2": "R2 2xsinc", "R4_exact": "R4 exact"}

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    "text.color": "#111111", "axes.labelcolor": "#111111", "xtick.color": "#111111",
    "ytick.color": "#111111", "axes.edgecolor": "#333333", "font.size": 10,
    "axes.titlesize": 11, "legend.fontsize": 9,
})


def read(name):
    p = os.path.join(PILOT, name)
    if not os.path.exists(p):
        return []
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def fl(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


summ = read("pilot_rung_summary.csv")
fid = read("pilot_fidelity.csv")
runs = read("pilot_runs.csv")

# cell metadata from CSVs
meta = {}
for r in summ + fid + runs:
    meta[r["cell"]] = (fl(r["theta0_deg"]), fl(r["kR"]))


def cell_label(c):
    if c in meta:
        th, kr = meta[c]
        return "%g°\nkR=%s" % (th, ("%.3g" % kr))
    try:  # no data yet: parse from the cell name, e.g. c30_k1p5
        a, b = c.split("_")
        return "%s°\nkR=%s" % (a[1:], b[1:].replace("p", "."))
    except Exception:
        return c.replace("_", "\n")


S = {(r["cell"], r["rung"]): r for r in summ}
F = {(r["cell"], r["rung"], r["set"]): r for r in fid}


def text_color(rgba):
    lum = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
    return "white" if lum < 0.5 else "black"


def fmt_ratio(v):
    if not np.isfinite(v):
        return "n/a"
    if v > 1e3 or v < 1e-2:
        return "%.2e" % v
    return "%.3g" % v


def draw_map(ax, data, annot, cmap, vmin, vmax, title, cbar_label, fig):
    nr, nc = len(RUNGS), len(CELLS)
    M = np.full((nr, nc), np.nan)
    for i, rg in enumerate(RUNGS):
        for j, c in enumerate(CELLS):
            M[i, j] = data(c, rg)
    cm = plt.get_cmap(cmap)
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    for i in range(nr):
        for j in range(nc):
            v = M[i, j]
            if np.isfinite(v):
                col = cm(norm(v))
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, facecolor=col, edgecolor="white", lw=2))
                main, sub = annot(CELLS[j], RUNGS[i])
                tc = text_color(col)
                ax.text(j, i - 0.08, main, ha="center", va="center", color=tc, fontsize=11, fontweight="bold")
                ax.text(j, i + 0.27, sub, ha="center", va="center", color=tc, fontsize=7.5)
            else:
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, facecolor="#d9d9d9", edgecolor="white", lw=2))
                ax.text(j, i, "n/a", ha="center", va="center", color="#555555", fontsize=10)
    ax.set_xlim(-.5, nc - .5)
    ax.set_ylim(nr - .5, -.5)
    ax.set_xticks(range(nc))
    ax.set_xticklabels([cell_label(c) for c in CELLS])
    ax.set_yticks(range(nr))
    ax.set_yticklabels([RUNG_LABEL[r] for r in RUNGS])
    ax.xaxis.tick_top()
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(title, loc="left", pad=36)
    sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cm)
    cb = fig.colorbar(sm, ax=ax, pad=0.015, fraction=0.04)
    cb.set_label(cbar_label)


def fig_map():
    fig, axs = plt.subplots(2, 1, figsize=(12, 10))
    allreg = [fl(r["pick_regret_lowerbound_dB"]) for r in summ]
    allreg = [v for v in allreg if np.isfinite(v)]
    vmax = max(allreg + [1.0])

    def reg(c, rg):
        r = S.get((c, rg))
        return fl(r["pick_regret_lowerbound_dB"]) if r else np.nan

    def reg_annot(c, rg):
        r = S[(c, rg)]
        star = "*" if str(r["any_unconverged"]).strip().lower() == "true" else ""
        main = "%.1f%s" % (fl(r["pick_regret_lowerbound_dB"]), star)
        sub = "pred %.1f / true %.1f dB" % (fl(r["pick_pred_dB"]), fl(r["pick_true_dB"]))
        return main, sub

    draw_map(axs[0], reg, reg_annot, "magma_r", 0, vmax,
             "(a) Quality of each rung's chosen design",
             "regret lower bound of the rung's chosen design\n(dB, true exact power vs best exact design found)", fig)

    def tau(c, rg):
        r = F.get((c, rg, "all"))
        return fl(r["kendall_tau"]) if r else np.nan

    def tau_annot(c, rg):
        r = F[(c, rg, "all")]
        return "%.2f" % fl(r["kendall_tau"]), "median ratio " + fmt_ratio(fl(r["ratio_median"]))

    draw_map(axs[1], tau, tau_annot, "RdBu", -1, 1,
             "(b) Rank agreement with the exact operator (set 'all')",
             "Kendall tau vs exact power (set 'all')", fig)
    fig.text(0.01, 0.005,
             "* = at least one run of that rung did not converge (stopped at the iteration cap).  "
             "R2 (2xsinc) has the same optimisation result as R1 (sinc); only the fidelity ratio differs.\n"
             "Panel (a) sub-text: the rung's own predicted vs exact-evaluated reduction relative to the uniform design.  "
             "Panel (b) sub-text: median of (rung power / exact power).  Grey = no data yet.",
             fontsize=8.5, va="bottom", ha="left")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(PILOT, "pilot_regime_map." + ext), dpi=200)
    plt.close(fig)


def fig_pred_true():
    fig, ax = plt.subplots(figsize=(8, 7))
    cells = [c for c in CELLS if any(r["cell"] == c for r in runs)]
    cmap = plt.get_cmap("tab10")
    mk = {"R0_ERP": "o", "R1_sinc": "s", "R4_exact": "^"}
    lab = {"R0_ERP": "E (R0 ERP)", "R1_sinc": "S (R1 sinc)", "R4_exact": "X (R4 exact)"}
    xs, ys = [], []
    for r in runs:
        x, y = fl(r["pred_dB_vs_uniform"]), fl(r["true_dB_vs_uniform"])
        if not (np.isfinite(x) and np.isfinite(y)):
            continue
        ci = CELLS.index(r["cell"]) if r["cell"] in CELLS else 0
        ax.scatter(x, y, marker=mk.get(r["rung"], "D"), s=70, color=cmap(ci), edgecolor="black", lw=0.7, zorder=3)
        ax.annotate(r["run"], (x, y), textcoords="offset points", xytext=(6, 5), fontsize=8)
        xs.append(x); ys.append(y)
    if xs:
        lo = min(xs + ys); hi = max(xs + ys)
        pad = 0.05 * (hi - lo + 1e-9)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "--", color="#444444", lw=1, zorder=1)
        ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("predicted reduction vs uniform (dB, own model)")
    ax.set_ylabel("true reduction vs uniform (dB, exact)")
    ax.set_title("Predicted vs exact-evaluated reduction (labels = run id; below diagonal = model overpromises)", fontsize=9.5)
    ax.grid(alpha=0.3)
    h1 = [Line2D([], [], marker=mk[k], ls="", color="w", markerfacecolor="#999999", markeredgecolor="black",
                 markersize=8, label=lab[k]) for k in mk]
    h2 = [Line2D([], [], marker="o", ls="", color="w", markerfacecolor=cmap(CELLS.index(c)),
                 markeredgecolor="black", markersize=8, label=cell_label(c).replace("\n", " ")) for c in cells]
    h3 = [Line2D([], [], ls="--", color="#444444", label="y = x")]
    ax.legend(handles=h1 + h2 + h3, loc="best", framealpha=0.9)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(PILOT, "pilot_pred_vs_true." + ext), dpi=200)
    plt.close(fig)


def fig_grad():
    gcells = []
    G = {}
    for c in CELLS:
        rows = read("grad_%s.csv" % c)
        if rows:
            gcells.append(c)
            G[c] = rows
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.6))
    pts = ["uniform", "XU", "SU", "EU"]
    cmap = plt.get_cmap("tab10")
    mk = {"E": "o", "S": "s", "X": "^"}
    n = max(len(gcells), 1)
    offs = np.linspace(-0.25, 0.25, n) if n > 1 else [0.0]
    # panel a: cos angle, rungs E and S
    for k, c in enumerate(gcells):
        seen = set()
        for r in G[c]:
            if r["rung"] not in ("E", "S") or r["point"] not in pts:
                continue
            key = (r["point"], r["rung"])
            if key in seen:
                continue
            seen.add(key)
            x = pts.index(r["point"]) + offs[k] + (-0.07 if r["rung"] == "E" else 0.07)
            axs[0].scatter(x, fl(r["cos_angle_vs_exact"]), marker=mk[r["rung"]], s=80,
                           color=cmap(CELLS.index(c)), edgecolor="black", lw=0.7, zorder=3)
    axs[0].axhline(1, color="#888888", lw=0.8, ls=":")
    axs[0].axhline(0, color="#333333", lw=0.8)
    axs[0].set_ylim(-1.05, 1.1)
    axs[0].set_xticks(range(len(pts)))
    axs[0].set_xticklabels(["uniform", "XU\n(exact from U)", "SU\n(sinc from U)", "EU\n(ERP from U)"])
    axs[0].set_xlabel("design point where gradients are compared")
    axs[0].set_ylabel("cosine of angle between own and exact gradient")
    axs[0].set_title("(a) Gradient direction agreement with exact operator")
    axs[0].grid(alpha=0.3)
    # panel b: exact dB change at step 5e-4
    allv = []
    for k, c in enumerate(gcells):
        seen = set()
        for r in G[c]:
            if r["point"] not in pts or abs(fl(r["step_maxdt"]) - 5e-4) > 1e-9:
                continue
            key = (r["point"], r["rung"])
            if key in seen:
                continue
            seen.add(key)
            v = fl(r["exact_dB_change"])
            allv.append(v)
            x = pts.index(r["point"]) + offs[k] + {"E": -0.07, "S": 0.0, "X": 0.07}[r["rung"]]
            axs[1].scatter(x, v, marker=mk[r["rung"]], s=80, color=cmap(CELLS.index(c)),
                           edgecolor="black", lw=0.7, zorder=3)
    axs[1].set_yscale("symlog", linthresh=0.01)
    axs[1].axhline(0, color="#333333", lw=0.8)
    axs[1].set_xticks(range(len(pts)))
    axs[1].set_xticklabels(["uniform", "XU", "SU", "EU"])
    axs[1].set_xlabel("design point")
    axs[1].set_ylabel("exact power change after a step (dB, step_maxdt = 5e-4)")
    axs[1].set_title("(b) Exact effect of a step along each rung's own descent direction")
    axs[1].grid(alpha=0.3)
    fig.text(0.99, 0.005, "(b): negative = exact power decreases (good); positive = step makes exact power worse; symlog axis",
             fontsize=8.5, va="bottom", ha="right")
    h1 = [Line2D([], [], marker=mk[k], ls="", color="w", markerfacecolor="#999999", markeredgecolor="black",
                 markersize=8, label=l) for k, l in (("E", "E (R0 ERP)"), ("S", "S (R1 sinc)"), ("X", "X (exact)"))]
    h2 = [Line2D([], [], marker="o", ls="", color="w", markerfacecolor=cmap(CELLS.index(c)),
                 markeredgecolor="black", markersize=8, label=cell_label(c).replace("\n", " ")) for c in gcells]
    axs[1].legend(handles=h1 + h2, loc="upper right", framealpha=0.9)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(PILOT, "pilot_gradient_test." + ext), dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    fig_map()
    fig_pred_true()
    fig_grad()
    print("done")
