"""Experiment 9 figures (white background, dark text). Reads results/r0_revision/passrule/*.csv written by
tests/r0_passrule_stats.py. Writes passrule_scatter.{png,pdf} and passrule_bars.{png,pdf} in the same folder."""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "r0_revision", "passrule")
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
                     "text.color": "#111111", "axes.labelcolor": "#111111", "xtick.color": "#111111",
                     "ytick.color": "#111111", "axes.edgecolor": "#333333", "font.size": 9})


def read(fn):
    with open(os.path.join(OUT, fn), newline="") as f:
        return list(csv.DictReader(f))


D = read("passrule_designs.csv")
ENV = {r["case"]: (float(r["a"]), float(r["b"])) for r in read("passrule_envelope.csv") if r["kind"] == "hull"}
G = read("passrule_guarantee.csv")
f = lambda r, k: float(r[k])

STYLE = {  # population: (colour, marker, size, label)
    "R": ("#7f7f7f", "o", 14, "random smooth fields"),
    "X": ("#1f77b4", "o", 6, "exact-driven iterates"),
    "E": ("#2ca02c", "o", 6, "ERP-driven iterates"),
    "P": ("#9467bd", "s", 22, "paper designs"),
    "S": ("#d62728", "o", 6, "sinc-driven iterates"),
}

# ------------------------------------------------------------------ figure 1: e against R_pred
fig, axs = plt.subplots(2, 2, figsize=(10, 8.2))
panels = [("paper", "(a) paper case (30°, kR 18.3): 3120 path/paper designs + 40 random"),
          ("c30_k3", "(b) 30°, kR 3 (pilot)"), ("c120_k1p5", "(c) 120°, kR 1.5 (pilot)")]
for ax, (case, title) in zip(axs.flat[:3], panels):
    A = [d for d in D if d["case"] == case]
    for drv in ["R", "E", "X", "P", "S"]:
        pts = [d for d in A if d["driver"] == drv]
        if not pts:
            continue
        c, m, s, lab = STYLE[drv]
        ax.scatter([f(d, "R_sinc") for d in pts], [f(d, "e_sinc") for d in pts], s=s, c=c, marker=m,
                   alpha=0.6 if drv in "XSE" else 0.9, lw=0, label=f"{lab} (n={len(pts)})")
    ends = [d for d in A if d["is_end"] == "True" and d["driver"] == "S"]
    ax.scatter([f(d, "R_sinc") for d in ends], [f(d, "e_sinc") for d in ends], s=70, facecolor="none",
               edgecolor="k", lw=1.3, label="sinc end designs", zorder=5)
    xs = np.array([f(d, "R_sinc") for d in A])
    xg = np.linspace(min(0, xs.min()), xs.max() * 1.05 + 1, 50)
    a, b = ENV[case]
    S = [f(d, "R_sinc") for d in A if d["driver"] == "S"]
    xe = np.linspace(min(S), max(S), 20)
    ax.plot(xe, a + b * xe, color="#d62728", lw=1.6, ls="-", label=f"envelope on sinc iterates: {a:.1f} + {b:.2f} R")
    emax_r = max(f(d, "e_sinc") for d in A if d["driver"] == "R")
    ax.axhline(emax_r, color="#7f7f7f", lw=1.0, ls=":", label=f"e_max random = {emax_r:.2f} dB")
    for T, ls in ((20, "--"), (40, "-.")):
        ax.plot(xg, xg - T, color="#111111", lw=1.0, ls=ls, label=f"R_true = {T} dB (T = {T})")
    ax.axhline(0, color="#bbbbbb", lw=0.6)
    yv = np.array([f(d, "e_sinc") for d in A])
    ax.set_ylim(max(yv.min() - 2, -60), yv.max() + 3 + (0.5 if case == "paper" else 0))
    ax.set_xlim(xg[0], xg[-1])
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("predicted reduction R_sinc (dB vs uniform)")
    ax.set_ylabel("e = R_sinc − R_exact (dB), > 0: sinc overstates")
    ax.legend(fontsize=6.3, loc="lower right" if case != "paper" else "lower left", framealpha=0.95)
    ax.grid(alpha=0.25)
axs[0, 0].set_ylim(-2.6, 1.4)
axs[0, 1].set_ylim(-65, 25)  # all points shown (min e = -60.6 dB)
# panel (d): sinc-driven iterates of all cells, envelope per kR class
ax = axs[1, 1]
cols = {"c30_k1": "#d62728", "c30_k3": "#ff7f0e", "c90_k1p5": "#8c564b", "c120_k1p5": "#e377c2",
        "c90_k7": "#1f77b4", "c120_k7": "#17becf", "c30_k18": "#2ca02c"}
for case, c in cols.items():
    pts = [d for d in D if d["case"] == case and d["driver"] == "S"]
    ax.scatter([f(d, "R_sinc") for d in pts], [f(d, "e_sinc") for d in pts], s=5, c=c, lw=0, alpha=0.7, label=case)
for key, c, lab in (("pooled_lowkR", "#d62728", "kR ≤ 3"), ("pooled_highkR", "#1f77b4", "kR ≥ 7")):
    a, b = ENV[key]
    xg = np.linspace(0, 75, 20)
    ax.plot(xg, a + b * xg, color=c, lw=1.6, label=f"envelope {lab}: {a:.1f} + {b:.3f} R")
for T, ls in ((20, "--"), (40, "-.")):
    xg = np.linspace(0, 75, 20)
    ax.plot(xg, xg - T, color="#111111", lw=1.0, ls=ls, label=f"R_true = {T} dB")
ax.set_xlim(0, 75); ax.set_ylim(-3, 40)
ax.set_title("(d) sinc-driven iterates, all 7 pilot cells", fontsize=9)
ax.set_xlabel("predicted reduction R_sinc (dB vs uniform)")
ax.set_ylabel("e = R_sinc − R_exact (dB)")
ax.legend(fontsize=6.3, loc="upper left", ncol=2, framealpha=0.95)
ax.grid(alpha=0.25)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, f"passrule_scatter.{ext}"), dpi=170)
plt.close(fig)

# ------------------------------------------------------------------ figure 2: R_guaranteed vs true reduction
order = [("c30_k1", "SU"), ("c30_k1", "SX"), ("c30_k3", "SU"), ("c30_k3", "SX"), ("c90_k1p5", "SU"),
         ("c90_k1p5", "SX"), ("c120_k1p5", "SU"), ("c120_k1p5", "SX"), ("c90_k7", "SU"), ("c90_k7", "SX"),
         ("c120_k7", "SU"), ("c120_k7", "SX"), ("paper", "S0"), ("paper", "SX"), ("paper", "SR")]
srcs = [("random_max", "#9e9e9e", "e_max over 40 random fields"),
        ("sinc_paths_LOO_max", "#ff7f0e", "e_max over the other sinc run's iterates"),
        ("envelope_same_class_other_cells", "#1f77b4", "envelope e(R_pred), other cells of the same kR class")]
fig, ax = plt.subplots(figsize=(11, 4.8))
x = np.arange(len(order)); w = 0.26
look = {(r["case"], r["run"], r["emax_source"]): r for r in G}
for i, (src, c, lab) in enumerate(srcs):
    vals = [f(look[(cs, rn, src)], "R_guaranteed") for cs, rn in order]
    ax.bar(x + (i - 1) * w, vals, w, color=c, label=f"R_guaranteed: {lab}", edgecolor="white", lw=0.4)
Rp = [f(look[(cs, rn, "random_max")], "R_pred") for cs, rn in order]
Rt = [f(look[(cs, rn, "random_max")], "R_true") for cs, rn in order]
Rf = [f(look[(cs, rn, "random_max")], "R_true_fsi_phys") for cs, rn in order]
ax.scatter(x, Rp, marker="_", s=300, c="#d62728", lw=2, label="R_sinc predicted", zorder=6)
ax.scatter(x, Rt, marker="D", s=34, c="#111111", label="R_exact true (in vacuo)", zorder=7)
ax.scatter(x, Rf, marker="D", s=34, facecolor="none", edgecolor="#111111", lw=1.0,
           label="R_true with fluid loading, physical density (exp7)", zorder=7)
for T, ls in ((20, "--"), (40, "-.")):
    ax.axhline(T, color="#111111", lw=1.0, ls=ls)
    ax.text(len(order) - 0.45, T + 0.6, f"T = {T} dB", fontsize=8, ha="right")
ax.set_xticks(x)
ax.set_xticklabels([f"{('anchor' if cs == 'paper' else cs)}\n{rn}" for cs, rn in order], fontsize=7.5)
for xv in (7.5, 11.5):
    ax.axvline(xv, color="#cccccc", lw=0.8)
ax.text(3.5, 76, "kR ≤ 3", ha="center", fontsize=9); ax.text(9.5, 76, "kR = 7", ha="center", fontsize=9)
ax.text(13, 76, "paper case (kR 18.3)", ha="center", fontsize=9)
ax.set_ylim(0, 80)
ax.set_ylabel("reduction vs uniform (dB)")
ax.set_title("Pass rule R_guaranteed = R_sinc − e_max for every sinc end design, against the true reduction", fontsize=9)
ax.legend(fontsize=7, loc="upper right", bbox_to_anchor=(1.0, 0.93), framealpha=0.95)
ax.grid(axis="y", alpha=0.25)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, f"passrule_bars.{ext}"), dpi=170)
plt.close(fig)
print("figures written")
