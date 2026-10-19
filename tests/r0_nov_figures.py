"""Figures for the novelty study (white background, dark text). Reads only CSVs written by r0_nov_* scripts and
(read-only) the BEM CSVs already in results/r0_revision.  Writes PNG+PDF into results/r0_revision/novelty/."""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "r0_revision", "novelty")
plt.rcParams.update({"font.size": 10, "axes.edgecolor": "#222", "text.color": "#111", "axes.labelcolor": "#111",
                     "xtick.color": "#111", "ytick.color": "#111", "axes.facecolor": "white",
                     "figure.facecolor": "white", "savefig.facecolor": "white", "legend.frameon": False})
DES = ["uniform", "quadvel-ones", "inputpower", "dyncomp", "quadvel-I", "quadvel-rad"]
LET = {"uniform": "uniform", "quadvel-ones": "B", "inputpower": "C", "dyncomp": "D", "quadvel-I": "E",
       "quadvel-rad": "radiation obj."}
COL = ["#1f4e79", "#c0504d", "#4f8f3a", "#7f5fa0", "#d08a1e", "#2a9d9d"]


def rd(path):
    return list(csv.DictReader(open(path)))


def save(fig, name):
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), dpi=200)
    plt.close(fig)


def fig_validation():
    v = [r for r in rd(os.path.join(OUT, "cap_frozen_validation_M20.csv"))]
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    x = np.arange(len(DES))
    for off, lev, mk, lab in ((-0.2, "2", "o", "same trace, BEM m=2"), (0.0, "3", "s", "same trace, BEM m=3")):
        y = [float(next(r for r in v if r["level"] == lev and r["design"] == d)["spec_over_bem_same_level"]) for d in DES]
        ax.plot(x + off, y, mk, color="#1f4e79" if lev == "2" else "#c0504d", label=lab)
    y = [float(next(r for r in v if r["level"] == "8" and r["design"] == d)["spec_over_bem_extrap"]) for d in DES]
    ax.plot(x + 0.2, y, "^", color="#4f8f3a", label="finest trace vs extrapolated BEM")
    ax.axhline(1, color="#444", lw=0.8, ls="--")
    ax.set_xticks(x); ax.set_xticklabels([LET[d] for d in DES], rotation=20)
    ax.set_ylabel("spectral / bempp two-sided power"); ax.set_ylim(0.975, 1.008)
    ax.legend(fontsize=8, loc="lower left")
    save(fig, "nov_cap_validation")


def fig_ratio():
    sw = rd(os.path.join(OUT, "cap_frozen_sweep.csv"))
    mo = [r for r in rd(os.path.join(OUT, "cap_modes_ratio.csv")) if r["theta0_deg"] == "30" and r["field"] == "uniform"]
    dc = rd(os.path.join(OUT, "disc_two_sided.csv"))
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    for c, d in zip(COL, DES):
        rr = sorted([(float(r["kR"]), float(r["ratio"])) for r in sw if r["design"] == d])
        ax.plot(*zip(*rr), "-", color=c, lw=1.1, label=LET[d])
    ax.plot([float(r["kR"]) for r in mo], [float(r["ratio"]) for r in mo], "k--", lw=1.4, label="cap, uniform $v_n$")
    # flat disc of the same rim radius a = R sin(30 deg) = R/2  ->  kR = 2 ka
    ax.plot([2 * float(r["ka"]) for r in dc], [float(r["ratio_sinc_over_two_sided"]) for r in dc], ":", color="#555",
            lw=1.6, label="flat disc, same rim radius")
    bem = rd(os.path.join(ROOT, "results", "r0_revision", "bem_comparison.csv"))
    try:
        pts = [float(r["ratio_consistent_over_bem"]) for r in bem if r["geometry"] == "spherecap"]
        ax.plot([18.318] * len(pts), pts, "kx", ms=6, label="bempp (paper, 100 Hz)")
    except Exception:
        pass
    ax.axhline(0.5, color="#888", lw=0.8); ax.axhline(1.0, color="#bbb", lw=0.8)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.1, 40); ax.set_ylim(0.1, 3000)
    ax.set_xlabel("$kR$ (cap: $R=10$ m, $\\theta_0=30^\\circ$; fields frozen at 100 Hz)")
    ax.set_ylabel("$P_{\\rm sinc}/P_{\\rm two\\text{-}sided}$")
    ax.legend(fontsize=7.5, ncol=2, loc="upper right")
    save(fig, "nov_cap_ratio_vs_kR")


def fig_ranking():
    rk = rd(os.path.join(OUT, "cap_ranking.csv"))
    kR = np.array([float(r["kR"]) for r in rk]); tau = np.array([float(r["kendall_tau"]) for r in rk])
    lo = np.array([float(r["rho_min"]) for r in rk]); hi = np.array([float(r["rho_max"]) for r in rk])
    rmin = np.array([float(r["ratio_min"]) for r in rk]); rmax = np.array([float(r["ratio_max"]) for r in rk])
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.5))
    ax[0].plot(kR, tau, "o-", color="#1f4e79", ms=3.5)
    ax[0].axvline(18.318, color="#c0504d", lw=0.8, ls="--"); ax[0].text(18.9, 0.42, "100 Hz", color="#c0504d", fontsize=8)
    ax[0].set_xscale("log"); ax[0].set_ylim(0.4, 1.05); ax[0].set_xlabel("$kR$")
    ax[0].set_ylabel("Kendall $\\tau$ (six cap designs)")
    ax[1].fill_between(kR, rmin, rmax, color="#1f4e79", alpha=0.6, label="range over the 6 designs")
    ax[1].axhline(0.5, color="#555", lw=0.8)
    ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlim(0.1, 40); ax[1].set_ylim(0.1, 3000)
    ax[1].xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter()); ax[1].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax[1].set_xticks([0.1, 1, 4, 10, 40])
    ax[1].set_xlabel("$kR$"); ax[1].set_ylabel("$P_{\\rm sinc}/P_{\\rm two\\text{-}sided}$")
    ax[1].legend(fontsize=8, loc="upper right")
    save(fig, "nov_cap_ranking")


def fig_modes():
    mo = rd(os.path.join(OUT, "cap_modes_ratio.csv"))
    ths = [15, 30, 60, 90, 120, 150]
    fields = ["cos%d" % j for j in range(8)]
    M = np.zeros((len(ths), len(fields)))
    for i, th in enumerate(ths):
        for j, f in enumerate(fields):
            v = [float(r["ratio"]) for r in mo if int(r["theta0_deg"]) == th and r["field"] == f and float(r["kR"]) >= 10]
            M[i, j] = np.exp(np.mean(np.log(v)))
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    im = ax.imshow(np.log2(M), cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    for i in range(len(ths)):
        for j in range(len(fields)):
            ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=7.5, color="#111")
    ax.set_xticks(range(len(fields))); ax.set_xticklabels([f"$j={j}$" for j in range(8)])
    ax.set_yticks(range(len(ths))); ax.set_yticklabels([f"{t}$^\\circ$" for t in ths])
    ax.set_xlabel("axisymmetric mode $v_j=\\cos((j+1/2)\\pi\\theta/\\theta_0)$"); ax.set_ylabel("cap half-angle $\\theta_0$")
    cb = fig.colorbar(im, ax=ax); cb.set_label("$\\log_2$ of geometric-mean ratio, $10\\leq kR\\leq 40$")
    save(fig, "nov_cap_mode_map")


if __name__ == "__main__":
    fig_validation(); fig_ratio(); fig_ranking(); fig_modes()



