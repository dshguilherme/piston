"""Figure: sinc model versus two-sided BEM for frozen shell fields.
Run: C:\\Users\\dshgu\\.envs\\bempp\\python.exe r0_bem_comparison_figure.py
Reads bem_comparison.csv and bem_shell_levels.csv, writes bem-comparison.pdf/.png.
"""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course"
OUT = os.path.join(ROOT, "results", "r0_revision")


def rd(name):
    with open(os.path.join(OUT, name)) as f:
        return list(csv.DictReader(f))


comp = rd("bem_comparison.csv")
lev = rd("bem_shell_levels.csv")

plt.rcParams.update({"font.size": 9, "axes.labelsize": 9, "legend.fontsize": 8,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "font.family": "serif",
                     "figure.facecolor": "white", "axes.facecolor": "white",
                     "savefig.facecolor": "white", "text.color": "black",
                     "axes.labelcolor": "black", "xtick.color": "black", "ytick.color": "black",
                     "axes.edgecolor": "black"})

fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.4, 4.8), gridspec_kw={"width_ratios": [1.35, 1]})
fig.patch.set_facecolor("white")

# ---- (a) ratios
order = ([c for c in comp if c["geometry"] == "spherecap"] + [c for c in comp if c["geometry"] == "cylpanel"])
x = np.arange(len(order))
w = 0.38
rc = np.array([float(c["ratio_centroid_over_bem"]) for c in order])
rs = np.array([float(c["ratio_consistent_over_bem"]) for c in order])
unc = np.array([float(c["P_bem_uncertainty"]) for c in order])
pb = np.array([float(c["P_bem_extrap"]) for c in order])
rel = unc / pb
ax.bar(x - w / 2, rc, w, color="#9aa7b4", edgecolor="black", linewidth=0.5, label="sinc (centroid)")
ax.bar(x + w / 2, rs, w, color="#2b6a99", edgecolor="black", linewidth=0.5, label="sinc (spline-consistent)")
ax.errorbar(x + w / 2, rs, yerr=[rs * rel / (1 + 0 * rel), rs * rel], fmt="none", ecolor="black", capsize=2, lw=0.8)
ax.errorbar(x - w / 2, rc, yerr=[rc * rel, rc * rel], fmt="none", ecolor="black", capsize=2, lw=0.8)
ax.axhline(1.0, color="black", lw=0.8, ls="--")
ax.set_yscale("log")
ax.set_xticks(x)
ax.set_xticklabels([("cap " if c["geometry"] == "spherecap" else "panel ") + c["design"].replace("_", " ") for c in order],
                   rotation=40, ha="right")
ax.set_ylabel(r"$P_{\mathrm{sinc}}/P_{\mathrm{BEM}}$")
ax.set_title("(a) Power ratio to two-sided BEM", fontsize=9)
ax.legend(frameon=False, loc="upper center", ncol=2, fontsize=7.5)
ax.grid(True, which="major", axis="y", color="0.85", lw=0.5)
ax.set_axisbelow(True)

# ---- (b) convergence
focus = [("spherecap", "quadvel-rad", "cap quadvel-rad", "o"),
         ("cylpanel", "legacy_long_horizon", "panel long-horizon", "s")]
colors = {"spherecap": "#b3452c", "cylpanel": "#2b6a99"}
for g, d, lab, mk in focus:
    # base-quadrature levels only (reg4_sing4)
    L = sorted([l for l in lev if l["geometry"] == g and l["design"] == d and l["quad_setting"] == "reg4_sing4"],
               key=lambda l: int(l["m"]))
    h = np.array([float(l["h_max"]) for l in L])
    pbem = np.array([float(l["P_bem"]) for l in L])
    psm = np.array([float(l["P_sinc_mesh"]) for l in L])
    cc = [c for c in comp if c["geometry"] == g and c["design"] == d][0]
    eb = float(cc["P_bem_extrap"])
    es = float(cc["P_consistent_sinc"])
    col = colors[g]
    bx.plot(h, pbem / eb, ls="-", marker=mk, color=col, mfc=col, ms=4, lw=1.0, label=f"two-sided BEM ({lab})")
    bx.plot(h, psm / es, ls="--", marker=mk, color=col, mfc="white", ms=4, lw=1.0,
            label=f"sinc (spline-consistent), on mesh ({lab})")
    bx.plot([-0.025 if g == "spherecap" else 0.025], [1], marker="*", color=col, ms=9, mec="black", mew=0.5, ls="none")
bx.axhline(1.0, color="black", lw=0.6, ls=":")
bx.set_xlim(left=-0.03)
bx.set_xlabel(r"element size $h_{\max}$ (m)")
bx.set_ylabel("power / extrapolated value")
bx.set_title("(b) Mesh convergence", fontsize=9)
bx.grid(True, color="0.85", lw=0.5)
bx.set_axisbelow(True)
h_, l_ = bx.get_legend_handles_labels()
bx.legend(h_, l_, frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.2))
bx.annotate(r"stars: extrapolated values ($h\to0$)", xy=(0.0, 1.0), xytext=(0.12, 0.2),
            textcoords="axes fraction", fontsize=7, color="black",
            arrowprops=dict(arrowstyle="-", lw=0.5, color="black"))
ax.set_ylim(0.01, 6)

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, f"bem-comparison.{ext}"), dpi=300, facecolor="white")
print("saved")


