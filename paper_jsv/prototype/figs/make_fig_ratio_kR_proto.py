"""Prototype restyle of fig_ratio_kR (figure-first rules): greyscale + line style, 137 mm wide,
8 pt sans words / math in CM, named regime bands, top ka axis, reference lines labelled at their end.
Reads only existing CSVs (same sources as tests/r0_d2_paper_figures.py: fig_ratio_kR). No solver runs."""
import os, csv, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple

ROOT = r"C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course"
NOV = os.path.join(ROOT, "results", "r0_revision", "novelty")
VIZ = os.path.join(ROOT, "results", "r0_revision", "viz")
OUT = os.path.join(ROOT, "paper_jsv", "prototype", "figs")

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["cmss10", "DejaVu Sans"], "mathtext.fontset": "cm",
    "font.size": 8, "axes.labelsize": 8, "legend.fontsize": 7, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white", "savefig.transparent": False,
    "axes.linewidth": 0.6, "lines.linewidth": 0.9, "legend.frameon": False,
    "pdf.fonttype": 42, "axes.unicode_minus": False,
})

def rd(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))

def fl(x):
    try: return float(x)
    except Exception: return np.nan

DES = ["uniform", "quadvel-ones", "inputpower", "dyncomp", "quadvel-I", "quadvel-rad"]
KR100 = 2 * np.pi * 100 / 343.0 * 10.0
sw = rd(os.path.join(NOV, "cap_frozen_sweep.csv"))
mo = [r for r in rd(os.path.join(NOV, "cap_modes_ratio.csv")) if r["theta0_deg"] == "30" and r["field"] == "uniform"]
dc = rd(os.path.join(NOV, "disc_two_sided.csv"))
f4 = rd(os.path.join(VIZ, "numbers_fig4.csv"))

mm = 1 / 25.4
fig = plt.figure(figsize=(137 * mm, 84 * mm))
ax = fig.add_axes([0.095, 0.29, 0.815, 0.58])

# named regimes (annotation, not data): light band where rho crosses one and reaches its minimum
ax.axvspan(1.0, 10.0, color="0.93", lw=0, zorder=0)
for x, s in ((0.32, "dipole regime\n$\\rho\\gg1$"), (3.16, "transition and minimum"), (20.0, "two-face plateau\n$\\rho\\to0.5$")):
    ax.text(x, 1.6e3, s, ha="center", va="top", fontsize=7, color="0.25")

others = [d for d in DES if d != "quadvel-rad"]
for d in others:
    rr = sorted((fl(r["kR"]), fl(r["ratio"])) for r in sw if r["design"] == d)
    ax.plot(*zip(*rr), "-", color="0.55", lw=0.7, zorder=2)
rr = sorted((fl(r["kR"]), fl(r["ratio"])) for r in sw if r["design"] == "quadvel-rad")
ax.plot(*zip(*rr), "-", color="black", lw=1.3, zorder=3)
mo = sorted((fl(r["kR"]), fl(r["ratio"])) for r in mo)
ax.plot(*zip(*mo), color="black", ls=(0, (4, 1.6)), lw=0.9, zorder=3)
dd = sorted((2 * fl(r["ka"]), fl(r["ratio_sinc_over_two_sided"])) for r in dc if 0.1 <= 2 * fl(r["ka"]) <= 40)
ax.plot(*zip(*dd), color="black", ls=(0, (1, 1.4)), lw=1.1, zorder=3)

for d, mk, fc in (("quadvel-rad", "o", "black"), ("uniform", "s", "0.6")):
    rows = sorted((fl(r["kR"]), fl(r["value"]), int(fl(r["resolved"]))) for r in f4
                  if r["quantity"] == "ratio_consistent_over_bem_far" and r["design"] == d)
    k = np.array([r[0] for r in rows]); y = np.array([r[1] for r in rows]); res = np.array([r[2] for r in rows], bool)
    ax.plot(k[res], y[res], mk, mfc=fc, mec="black", mew=0.5, ms=4.2, ls="none", zorder=4)
    ax.plot(k[~res], y[~res], mk, mfc="white", mec="black", mew=0.8, ms=4.2, ls="none", zorder=4)

ax.axhline(1, color="0.35", lw=0.6, zorder=1); ax.axhline(0.5, color="0.35", lw=0.6, ls=(0, (6, 2)), zorder=1)
ax.text(41.5, 1.0, "$\\rho=1$", ha="left", va="center", fontsize=7)
ax.text(41.5, 0.5, "$\\rho=0.5$", ha="left", va="center", fontsize=7)
ax.axvline(KR100, ymax=0.76, color="0.35", lw=0.6, ls=(0, (1, 1.5)), zorder=1)
ax.text(KR100 * 0.96, 0.13, "100 Hz", ha="right", va="bottom", fontsize=7)

ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.1, 40); ax.set_ylim(0.1, 3000)
ax.set_xlabel("$kR$ (cap radius of curvature $R$ = 10 m)")
ax.set_ylabel("$\\rho=P_{\\rm sinc}/P_{\\rm 2s}$")
fmt = mticker.FuncFormatter(lambda v, _: f"{v:g}")
ax.xaxis.set_major_formatter(fmt)
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: "$10^{%d}$" % round(np.log10(v))))
top = ax.secondary_xaxis("top", functions=(lambda x: x / 2, lambda x: 2 * x))
top.set_xlabel("$ka=kR/2$ (rim radius $a$ = 5 m)", labelpad=3)
top.xaxis.set_major_formatter(fmt)
for s in ("right",):
    ax.spines[s].set_visible(False)

h = [Line2D([], [], color="black", lw=1.3), Line2D([], [], color="0.55", lw=0.7),
     Line2D([], [], color="black", ls=(0, (4, 1.6)), lw=0.9), Line2D([], [], color="black", ls=(0, (1, 1.4)), lw=1.1),
     Line2D([], [], marker="o", ls="none", mfc="black", mec="black", ms=4.2),
     Line2D([], [], marker="s", ls="none", mfc="0.6", mec="black", mew=0.5, ms=4.2),
     (Line2D([], [], marker="o", ls="none", mfc="white", mec="black", mew=0.8, ms=4.2),
      Line2D([], [], marker="s", ls="none", mfc="white", mec="black", mew=0.8, ms=4.2))]
l = ["design G, frozen", "other five designs, frozen", "uniform-velocity cap", "flat disc, same rim radius",
     "BEM far field, design G", "BEM far field, uniform-thickness design", "BEM, under-resolved mesh"]
fig.legend(h, l, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 0.0), columnspacing=1.4, handlelength=2.4,
           handletextpad=0.5, labelspacing=0.35,
           handler_map={tuple: HandlerTuple(ndivide=None, pad=0.6)})
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, f"fig_ratio_kR_proto.{ext}"), dpi=200)
print("ok", fig.get_size_inches() * 25.4)
