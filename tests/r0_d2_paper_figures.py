"""Publication figures (direction 2) for the radiation-operator manuscript.

Plots only data already stored under results/r0_revision/ (CSV / NPZ / MAT) plus closed-form curves (scipy).
No BEM assembly or solve, no spectral solver run.  Writes PDF + PNG (200 dpi) into results/r0_revision/paper_d2/.
Run:  python tests/r0_d2_paper_figures.py [stem ...]     (no argument = all figures)
"""
import os, sys, csv
import numpy as np
import scipy.io as sio
from scipy import special as sps
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, "results", "r0_revision")
NOV = os.path.join(RES, "novelty")
VIZ = os.path.join(RES, "viz")
EXO = os.path.join(RES, "exact_opt")
OUT = os.path.join(RES, "paper_d2")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["cmr10", "DejaVu Serif"], "mathtext.fontset": "cm",
    "axes.formatter.use_mathtext": True, "axes.unicode_minus": False,
    "font.size": 9.5, "axes.labelsize": 9.5, "legend.fontsize": 8, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white", "savefig.transparent": False,
    "text.color": "black", "axes.labelcolor": "black", "xtick.color": "black", "ytick.color": "black",
    "axes.edgecolor": "black", "legend.frameon": False, "axes.linewidth": 0.7, "lines.linewidth": 1.2,
    "pdf.fonttype": 42, "ps.fonttype": 42, "figure.dpi": 100,
})

# internal design keys -> plain labels (order used everywhere)
DES = ["uniform", "quadvel-ones", "inputpower", "dyncomp", "quadvel-I", "quadvel-rad"]
LAB = {"uniform": "uniform", "quadvel-ones": "B", "inputpower": "C", "dyncomp": "D", "quadvel-I": "E", "quadvel-rad": "R"}
# colourblind-safe (Okabe-Ito)
DCOL = {"uniform": "#0072B2", "quadvel-ones": "#E69F00", "inputpower": "#009E73", "dyncomp": "#56B4E9",
        "quadvel-I": "#CC79A7", "quadvel-rad": "#D55E00"}
RUNS = ["S0", "SR", "SX", "X0", "XR", "XS"]
RLAB = {"S0": "sinc, from uniform", "SR": "sinc, from design R", "SX": "sinc, from exact optimum",
        "X0": "exact, from uniform", "XR": "exact, from design R", "XS": "exact, from sinc optimum"}
RCOL = {"S0": "#0072B2", "SR": "#009E73", "SX": "#56B4E9", "X0": "#D55E00", "XR": "#CC79A7", "XS": "#E69F00"}
KR100 = 2 * np.pi * 100 / 343.0 * 10.0          # 18.32
GREY = "#555555"


def rd(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def fl(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"{stem}.{ext}"), dpi=200, facecolor="white", transparent=False,
                    bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    print("wrote", stem)


def panel_tag(ax, s, x=-0.02, y=1.02, **kw):
    kw.setdefault("va", "bottom")
    ax.text(x, y, s, transform=ax.transAxes, ha="left", fontsize=10, **kw)


# --------------------------------------------------------------------------------------------------------------
def h1p(l, x):
    return sps.spherical_jn(l, x, True) + 1j * sps.spherical_yn(l, x, True)


def r_l(l, x):
    return 2 * x ** 4 * sps.spherical_jn(l, x) ** 2 * np.abs(h1p(l, x)) ** 2


def fig_sphere_modes():
    chk = [r_l(l, 1.0) for l in (0, 1, 2)]
    print("  r_l(kR=1), l=0,1,2:", ["%.4f" % c for c in chk], "(expected 2.8323 0.9070 0.6850)")
    assert np.allclose(chk, [2.8323, 0.9070, 0.6850], atol=6e-5)
    ls = (0, 1, 3, 6)
    cols = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.7), gridspec_kw={"width_ratios": [1, 1]})
    x = np.linspace(0.02, 40, 20000)
    for l, c in zip(ls, cols):
        ax[0].plot(x, r_l(l, x), color=c, lw=1.0, label=f"$\\ell={l}$")
        lim = 2 * (l + 1) ** 2 / (2 * l + 1) ** 2
        ax[0].plot([0, 4.0], [lim, lim], color=c, lw=1.5, ls=(0, (2.5, 1.5)), alpha=0.95)
    ax[0].axhline(1, color=GREY, lw=0.8)
    ax[0].set_xlim(0, 40); ax[0].set_ylim(0, 3.7)
    ax[0].set_xlabel("$kR$"); ax[0].set_ylabel("$P_{\\rm sinc}/P_{\\rm exact}$ (closed sphere)")
    ax[0].legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.0), columnspacing=1.0, handlelength=1.4)
    panel_tag(ax[0], "(a)")
    xs = np.linspace(0.02, 120, 60001)
    dx = xs[1] - xs[0]
    nwin = int(round(np.pi / dx))
    for l, c in zip(ls, cols):
        r = r_l(l, xs)
        cs = np.concatenate([[0], np.cumsum(r)])
        rm = (cs[nwin:] - cs[:-nwin]) / nwin
        xc = xs[: len(rm)] + 0.5 * nwin * dx
        ax[1].plot(xc, rm, color=c, lw=1.1, label=f"$\\ell={l}$")
    ax[1].axhline(1, color=GREY, lw=0.8)
    ax[1].set_xlim(0, 120); ax[1].set_ylim(0, 3.7)
    ax[1].set_xlabel("$kR$"); ax[1].set_ylabel("running mean over $\\Delta kR=\\pi$")
    panel_tag(ax[1], "(b)")
    fig.tight_layout(w_pad=1.5)
    save(fig, "fig_sphere_modes")


# --------------------------------------------------------------------------------------------------------------
def fig_validation():
    v = rd(os.path.join(NOV, "cap_frozen_validation_M20.csv"))
    fe = rd(os.path.join(NOV, "cap_flat_extrap.csv"))
    dsc = {round(fl(r["ka"]), 6): fl(r["ratio_two_sided_over_dipole_asymptote"]) for r in rd(os.path.join(NOV, "disc_two_sided.csv"))}
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.9), gridspec_kw={"width_ratios": [1.15, 1]})
    a = ax[0]
    x = np.arange(len(DES))
    spec = [("2", "spec_over_bem_same_level", -0.22, "o", "#0072B2", "coarse mesh, same trace"),
            ("3", "spec_over_bem_same_level", 0.0, "s", "#D55E00", "fine mesh, same trace"),
            ("8", "spec_over_bem_extrap", 0.22, "^", "#009E73", "finest trace, extrapolated BEM")]
    for lev, col, off, mk, c, lab in spec:
        y = [fl(next(r for r in v if r["level"] == lev and r["design"] == d)[col]) for d in DES]
        a.plot(x + off, y, mk, color=c, ms=5, label=lab)
    a.axhline(1, color=GREY, lw=0.8, ls="--")
    a.set_xticks(x); a.set_xticklabels([LAB[d] for d in DES])
    a.set_ylabel("exact / BEM two-sided power")
    a.set_xlabel("cap design"); a.set_ylim(0.975, 1.012); a.set_xlim(-0.6, 5.6)
    a.legend(loc="upper right", fontsize=7.5, handletextpad=0.3, borderaxespad=0.2)
    panel_tag(a, "(a)")
    b = ax[1]
    cols = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]
    for r, c in zip(fe, cols):
        ka = fl(r["ka"]); d = dsc[round(ka, 6)]
        pts = sorted((float(p.split(":")[0]), float(p.split(":")[1])) for p in r["cap_values"].split(";"))
        xx = [0.0] + [p[0] ** 2 for p in pts]
        yy = [fl(r["cap_extrap"]) / d] + [p[1] / d for p in pts]
        b.plot(xx[1:], yy[1:], "-o", color=c, ms=3.5, lw=1.0, label=f"$ka={ka:g}$")
        b.plot(xx[0], yy[0], "D", color=c, ms=5.5, mec="black", mew=0.6)
    b.axhline(1, color=GREY, lw=0.8, ls="--")
    b.set_xlabel("$(a/R)^2$"); b.set_ylabel("cap / flat disc, two-sided power")
    b.set_xlim(-0.012, 0.175); b.set_ylim(0.98, 1.2)
    h, l = b.get_legend_handles_labels()
    h.append(Line2D([], [], marker="D", ls="none", color="white", mfc="#999999", mec="black", ms=5.5))
    l.append("extrapolated to $(a/R)^2\\to0$")
    b.legend(h, l, loc="upper right", fontsize=7.5, handletextpad=0.4, borderaxespad=0.2)
    panel_tag(b, "(b)")
    fig.tight_layout(w_pad=1.6)
    save(fig, "fig_validation")


# --------------------------------------------------------------------------------------------------------------
def fig_ratio_kR():
    sw = rd(os.path.join(NOV, "cap_frozen_sweep.csv"))
    mo = [r for r in rd(os.path.join(NOV, "cap_modes_ratio.csv")) if r["theta0_deg"] == "30" and r["field"] == "uniform"]
    dc = rd(os.path.join(NOV, "disc_two_sided.csv"))
    f4 = rd(os.path.join(VIZ, "numbers_fig4.csv"))
    fig, ax = plt.subplots(figsize=(6.5, 3.9))
    for d in DES:
        rr = sorted((fl(r["kR"]), fl(r["ratio"])) for r in sw if r["design"] == d)
        ax.plot(*zip(*rr), "-", color=DCOL[d], lw=1.2, label=LAB[d])
    mo = sorted((fl(r["kR"]), fl(r["ratio"])) for r in mo)
    ax.plot(*zip(*mo), color="black", ls="--", lw=1.3, label="uniform-velocity cap")
    dd = sorted((2 * fl(r["ka"]), fl(r["ratio_sinc_over_two_sided"])) for r in dc if 0.1 <= 2 * fl(r["ka"]) <= 40)
    ax.plot(*zip(*dd), color="black", ls=":", lw=1.6, label="flat disc, same rim radius")
    for d, mk in (("quadvel-rad", "o"), ("uniform", "s")):
        rows = sorted((fl(r["kR"]), fl(r["value"]), int(fl(r["resolved"]))) for r in f4
                      if r["quantity"] == "ratio_consistent_over_bem_far" and r["design"] == d)
        k = np.array([r[0] for r in rows]); y = np.array([r[1] for r in rows]); res = np.array([r[2] for r in rows], bool)
        ax.plot(k[res], y[res], mk, color=DCOL[d], mfc=DCOL[d], mec="black", mew=0.5, ms=5, ls="none")
        ax.plot(k[~res], y[~res], mk, color=DCOL[d], mfc="white", mec=DCOL[d], mew=1.3, ms=5.5, ls="none")
    ax.axhline(1, color=GREY, lw=0.8); ax.axhline(0.5, color=GREY, lw=0.8)
    ax.axvline(KR100, color=GREY, lw=0.9, ls=":")
    ax.text(KR100 * 1.04, 2.4e3, "100 Hz", fontsize=8.5, ha="left", va="top")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.1, 40); ax.set_ylim(0.1, 3000)
    ax.set_xlabel("$kR$"); ax.set_ylabel("$P_{\\rm sinc}\\,/\\,P_{\\rm two\\text{-}sided}$")
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.text(0.16, 40, "dipole regime", fontsize=8.5, rotation=-58, ha="center", va="center", color=GREY)
    ax.text(22, 0.31, "two-face plateau", fontsize=8.5, ha="center", va="top", color=GREY)
    h, l = ax.get_legend_handles_labels()
    h += [Line2D([], [], marker="o", ls="none", mfc=DCOL["quadvel-rad"], mec="black", mew=0.5, ms=5),
          Line2D([], [], marker="s", ls="none", mfc=DCOL["uniform"], mec="black", mew=0.5, ms=5),
          Line2D([], [], marker="o", ls="none", mfc="white", mec=GREY, mew=1.2, ms=5.5)]
    l += ["BEM, design R", "BEM, uniform", "BEM, unresolved"]
    ax.legend(h, l, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.16), fontsize=7.8, columnspacing=1.2,
              handlelength=1.8, handletextpad=0.4, borderaxespad=0.2, labelspacing=0.3)
    fig.tight_layout()
    save(fig, "fig_ratio_kR")


# --------------------------------------------------------------------------------------------------------------
def fig_regime_map():
    mo = rd(os.path.join(NOV, "cap_modes_ratio.csv"))
    ths = [15, 30, 60, 90, 120, 150]
    M = np.zeros((len(ths), 8))
    for i, th in enumerate(ths):
        for j in range(8):
            v = [fl(r["ratio"]) for r in mo if int(r["theta0_deg"]) == th and r["field"] == f"cos{j}" and fl(r["kR"]) >= 10]
            M[i, j] = np.exp(np.mean(np.log(v)))
    fig, ax = plt.subplots(figsize=(4.6, 3.3))
    im = ax.imshow(np.log2(M), cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    for i in range(len(ths)):
        for j in range(8):
            ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=7.5, color="black")
    ax.set_xticks(range(8)); ax.set_xticklabels([str(j) for j in range(8)])
    ax.set_yticks(range(len(ths))); ax.set_yticklabels([f"${t}^\\circ$" for t in ths])
    ax.set_xlabel("axisymmetric mode index $j$"); ax.set_ylabel("cap half-angle $\\theta_0$")
    cb = fig.colorbar(im, ax=ax, pad=0.03, fraction=0.06)
    cb.set_ticks([-2, -1, 0, 1, 2]); cb.set_ticklabels(["1/4", "1/2", "1", "2", "4"])
    cb.set_label("geometric-mean $P_{\\rm sinc}/P_{\\rm two\\text{-}sided}$, $10\\leq kR\\leq40$")
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    print("  map range:", M.min(), M.max())
    fig.tight_layout()
    save(fig, "fig_regime_map")


# --------------------------------------------------------------------------------------------------------------
def load_cap_npz():
    Z = np.load(os.path.join(VIZ, "cap_m2_100hz_bem.npz"), allow_pickle=True)
    return Z, [str(s) for s in Z["labels"]]


def fig_face_split():
    Z, labels = load_cap_npz()
    V, T = Z["vertices"], Z["triangles"]
    tri = mtri.Triangulation(V[:, 0], V[:, 1], T)
    FW, FH = 6.5, 8.5
    fig = plt.figure(figsize=(FW, FH))
    pos = lambda l_, b_, w_, h_: [l_ / FW, b_ / FH, w_ / FW, h_ / FH]
    a = fig.add_axes(pos(0.75, 6.85, 5.65, 1.5))
    x = np.arange(len(DES)); w = 0.26
    out_s, in_s, sinc_s = [], [], []
    for d in DES:
        j = labels.index(d)
        Pb = Z["P_bem"][j]
        out_s.append(Z["P_plus"][j] / Pb); in_s.append(Z["P_minus"][j] / Pb); sinc_s.append(Z["P_sinc_mesh"][j] / Pb)
    print("  shares (outer, inner, sinc):", [(LAB[d], round(o, 3), round(i, 3), round(s, 3)) for d, o, i, s in zip(DES, out_s, in_s, sinc_s)])
    bars = [(out_s, "outer-face share $P^{+}/P$", "#0072B2", None), (in_s, "inner-face share $P^{-}/P$", "#E69F00", None),
            (sinc_s, "$P_{\\rm sinc}/P$", "#999999", "//")]
    for k, (vals, lab, c, hat) in enumerate(bars):
        a.bar(x + (k - 1) * w, vals, w, color=c, edgecolor="black", linewidth=0.5, hatch=hat, label=lab)
        for xx, v in zip(x + (k - 1) * w, vals):
            a.text(xx, v + 0.012, f"{v:.2f}", ha="center", va="bottom", fontsize=6.3)
    a.axhline(0.5, color=GREY, lw=0.7, ls=":")
    a.set_xticks(x); a.set_xticklabels([LAB[d] for d in DES]); a.set_ylim(0, 0.9)
    a.set_ylabel("fraction of two-sided BEM power"); a.set_xlabel("cap design")
    a.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.0), fontsize=7.8, columnspacing=1.2, handlelength=1.6)
    panel_tag(a, "(a)", x=-0.13)
    tags = iter(["(b)", "(c)", "(d)", "(e)"])
    for r, d in enumerate(("quadvel-rad", "uniform")):
        j = labels.index(d)
        Ip, Im = Z["I_plus"][:, j], Z["I_minus"][:, j]
        lim = max(np.abs(Ip).max(), np.abs(Im).max())
        e = int(np.floor(np.log10(lim)))
        sc = 10.0 ** e
        for c, (f, face) in enumerate(((Ip, "outer face"), (Im, "inner face"))):
            ax = fig.add_axes(pos(0.75 + c * 2.55, 3.45 - r * 2.95, 2.35, 2.35))
            h = ax.tripcolor(tri, f / sc, shading="gouraud", cmap="RdBu_r", vmin=-lim / sc, vmax=lim / sc)
            ax.set_aspect("equal"); ax.set_xlim(-5.3, 5.3); ax.set_ylim(-5.3, 5.3)
            ax.set_xticks([-4, 0, 4]); ax.set_yticks([-4, 0, 4])
            ax.set_xlabel("$x$ (m)")
            if c == 0:
                ax.set_ylabel("$y$ (m)")
            else:
                ax.set_yticklabels([])
            nm = "design R" if d == "quadvel-rad" else "uniform"
            ax.set_title(f"{nm}, {face}", fontsize=9, pad=3)
            panel_tag(ax, next(tags), x=0.03, y=0.97, va="top")
        cax = fig.add_axes(pos(5.95, 3.45 - r * 2.95 + 0.3, 0.12, 1.75))
        cb = fig.colorbar(h, cax=cax)
        cb.set_label(f"active intensity ($10^{{{e}}}$ W/m$^2$)", fontsize=8)
        cb.ax.tick_params(labelsize=8)
    save(fig, "fig_face_split")


# --------------------------------------------------------------------------------------------------------------
def fig_directivity():
    Z, labels = load_cap_npz()
    DR = 40.0
    fig = plt.figure(figsize=(6.5, 3.5))
    for k, (d, nm, tag) in enumerate((("quadvel-rad", "design R", "(a)"), ("uniform", "uniform", "(b)"))):
        key = d.replace("-", "_")
        B, S = Z[f"dPdO_bem__{key}"], Z[f"dPdO_sinc__{key}"]
        mx = B.max()
        dB = lambda a_: 10 * np.log10(np.maximum(a_ / mx, 1e-30))
        al = np.arange(0, 361, 1.0)

        def cut(arr):
            a_ = dB(arr)
            v = np.where(al <= 180, a_[np.minimum(al, 180).astype(int), 0], a_[np.clip(360 - al, 0, 180).astype(int), 90])
            return np.deg2rad(al), np.maximum(v, -DR)
        ax = fig.add_subplot(1, 2, k + 1, projection="polar")
        a1, v1 = cut(B); a2, v2 = cut(S)
        ax.plot(a1, v1, color="black", lw=1.4, label="BEM")
        ax.plot(a2, v2, color="#D55E00", lw=1.4, ls="--", label="sinc model")
        ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
        ax.set_rlim(-DR, 0); ax.set_rticks([-30, -20, -10]); ax.set_rlabel_position(180)
        ax.set_thetagrids(np.arange(0, 360, 45), labels=[""] * 8)
        ax.tick_params(labelsize=7.5)
        ax.grid(color="#bbbbbb", lw=0.5)
        ax.set_title(f"{tag} {nm}", fontsize=9.5, pad=14)
        ax.text(0, 1.0, "$+z$ (outer)", transform=ax.transAxes, ha="center", va="bottom", fontsize=8,
                clip_on=False) if False else None
        ax.annotate("$+z$ (outer, convex)", xy=(0.5, 1.0), xycoords="axes fraction", xytext=(0.5, 1.015),
                    textcoords="axes fraction", ha="center", va="bottom", fontsize=8, annotation_clip=False)
        ax.annotate("$-z$ (inner)", xy=(0.5, 0.0), xycoords="axes fraction", xytext=(0.5, -0.03),
                    textcoords="axes fraction", ha="center", va="top", fontsize=8, annotation_clip=False)
        if k == 0:
            ax.legend(loc="lower left", bbox_to_anchor=(-0.22, -0.2), fontsize=8, handlelength=1.8)
        ax.text(0.5, -0.14, "dB re BEM maximum", transform=ax.transAxes, ha="center", va="top", fontsize=8)
    fig.subplots_adjust(left=0.05, right=0.95, top=0.82, bottom=0.14, wspace=0.28)
    save(fig, "fig_directivity")


# --------------------------------------------------------------------------------------------------------------
def fig_power_freq():
    f4 = rd(os.path.join(VIZ, "numbers_fig4.csv"))
    sn = rd(os.path.join(VIZ, "cap_sinc_freq_sweep.csv"))
    nv = rd(os.path.join(NOV, "cap_frozen_sweep.csv"))
    KRMAX = [fl(r["value"]) for r in f4 if r["quantity"] == "kR_resolution_limit_6_per_lambda_m2"][0]
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 3.3), gridspec_kw={"width_ratios": [1.55, 1]})
    a = ax[0]
    a.axvspan(KRMAX, 40, color="#e3e3e3", lw=0)
    for d in ("uniform", "quadvel-rad"):
        c = DCOL[d]
        D = sorted((fl(r["kR"]), fl(r["P_consistent"])) for r in sn if r["design"] == d and int(fl(r["dense_flag"])) == 1)
        a.plot(*zip(*D), "-", color=c, lw=1.3)
        N = sorted((fl(r["kR"]), fl(r["P_2s"])) for r in nv if r["design"] == d)
        a.plot(*zip(*N), "--", color=c, lw=1.3)
        rows = sorted((fl(r["kR"]), fl(r["value"]), int(fl(r["resolved"]))) for r in f4
                      if r["quantity"] == "P_bem_far_m2" and r["design"] == d)
        k = np.array([q[0] for q in rows]); y = np.array([q[1] for q in rows]); res = np.array([q[2] for q in rows], bool)
        a.plot(k[res], y[res], "o", mfc=c, mec="black", mew=0.4, ms=4.5, ls="none")
        a.plot(k[~res], y[~res], "o", mfc="white", mec=c, mew=1.2, ms=4.8, ls="none")
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(0.1, 40)
    a.axvline(KR100, color=GREY, lw=0.8, ls=":")
    a.text(KR100 / 1.06, 0.97, "100 Hz", transform=matplotlib.transforms.blended_transform_factory(a.transData, a.transAxes),
           ha="right", va="top", fontsize=8)
    a.set_xlabel("$kR$"); a.set_ylabel("radiated power $P$ (W)")
    a.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    a.text(0.5 * (KRMAX + 40) ** 1.0 * 0 + np.sqrt(KRMAX * 40) * 1.02, 0.5, "fewer than 6 elements per wavelength",
           transform=matplotlib.transforms.blended_transform_factory(a.transData, a.transAxes),
           rotation=90, ha="center", va="center", fontsize=7.5, color="black")
    h = [Line2D([], [], color="black", lw=1.3), Line2D([], [], color="black", lw=1.3, ls="--"),
         Line2D([], [], marker="o", ls="none", mfc="#999999", mec="black", mew=0.4, ms=4.5),
         Line2D([], [], marker="o", ls="none", mfc="white", mec="black", mew=1.1, ms=4.8)]
    leg1 = a.legend(h, ["consistent sinc", "exact (dual series)", "BEM far field", "BEM, unresolved"],
                    loc="upper left", fontsize=7.5, handlelength=1.8, borderaxespad=0.2)
    a.add_artist(leg1)
    a.legend([Patch(color=DCOL["quadvel-rad"]), Patch(color=DCOL["uniform"])], ["design R", "uniform"],
             loc="lower right", bbox_to_anchor=(0.80, 0.0), fontsize=7.5, handlelength=1.2, borderaxespad=0.2)
    panel_tag(a, "(a)", x=-0.03)
    b = ax[1]
    ks = [1, 4, 10, 16]
    P = {q: sorted((fl(r["kR"]), fl(r["value"])) for r in f4 if r["quantity"] == q and r["design"] == "quadvel-rad")
         for q in ("P_plus_over_P_surface", "P_minus_over_P_surface")}
    print("  face shares R:", P)
    b.axhline(0, color=GREY, lw=0.7)
    b.plot(*zip(*P["P_plus_over_P_surface"]), "-^", color="#0072B2", ms=4.5, label="outer face $P^{+}/P$")
    b.plot(*zip(*P["P_minus_over_P_surface"]), "-v", color="#E69F00", ms=4.5, label="inner face $P^{-}/P$")
    b.set_xscale("log"); b.set_yscale("symlog", linthresh=1.0, linscale=1.0)
    b.set_xlim(0.8, 20); b.set_ylim(-8, 8)
    b.set_xticks(ks); b.set_xticklabels([str(k) for k in ks])
    b.xaxis.set_minor_formatter(mticker.NullFormatter())
    b.set_yticks([-5, -1, 0, 1, 5]); b.set_yticklabels(["$-5$", "$-1$", "0", "1", "5"])
    b.set_xlabel("$kR$"); b.set_ylabel("share of BEM power")
    b.legend(loc="lower right", fontsize=7.5, handlelength=1.6, borderaxespad=0.2)
    panel_tag(b, "(b)", x=-0.02)
    fig.tight_layout(w_pad=1.4)
    save(fig, "fig_power_freq")


# --------------------------------------------------------------------------------------------------------------
def fig_bem_comparison():
    comp = rd(os.path.join(RES, "bem_comparison.csv"))
    lev = rd(os.path.join(RES, "bem_shell_levels.csv"))
    ex = rd(os.path.join(NOV, "cap_frozen_validation_M20.csv"))
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.5, 3.6), gridspec_kw={"width_ratios": [1.4, 1]})
    bem_r = [fl(next(c for c in comp if c["geometry"] == "spherecap" and c["design"] == d)["ratio_consistent_over_bem"]) for d in DES]
    ex_r = []
    for d in DES:
        r = next(q for q in ex if q["level"] == "8" and q["design"] == d)
        ex_r.append(fl(r["P_sinc_consistent"]) / fl(r["P2s_extrap"]))
    print("  BEM ratios:", np.round(bem_r, 3), " exact ratios:", np.round(ex_r, 3))
    pan = [fl(next(c for c in comp if c["geometry"] == "cylpanel" and c["design"] == d)["ratio_consistent_over_bem"])
           for d in ("uniform", "legacy_long_horizon")]
    print("  panel BEM ratios:", pan)
    w = 0.36
    x = np.arange(6)
    xp = np.array([7.3, 9.4])
    cb, cx = "#0072B2", "#E69F00"
    a.bar(x - w / 2, bem_r, w, color=cb, edgecolor="black", linewidth=0.5, label="BEM (extrapolated)")
    a.bar(x + w / 2, ex_r, w, color=cx, edgecolor="black", linewidth=0.5, label="exact (dual series)")
    a.bar(xp, pan, 0.55, color=cb, edgecolor="black", linewidth=0.5)
    a.axhline(1, color="black", lw=0.8, ls="--"); a.axhline(0.5, color=GREY, lw=0.8, ls=":")
    a.set_yscale("log"); a.set_ylim(0.01, 2.8)
    a.set_xticks(list(x) + list(xp))
    a.set_xticklabels([LAB[d] for d in DES] + ["panel,\nuniform", "panel,\nlong horizon"], rotation=0, fontsize=7.5)
    a.set_ylabel("$P_{\\rm sinc}/P_{\\rm two\\text{-}sided}$")
    a.set_xlim(-0.7, 10.5)
    a.legend(loc="upper left", ncol=2, fontsize=7.3, borderaxespad=0.2, columnspacing=1.0, handlelength=1.4)
    a.set_yticks([0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]); a.yaxis.set_minor_locator(mticker.NullLocator())
    a.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    panel_tag(a, "(a)", x=-0.12)
    focus = [("spherecap", "quadvel-rad", "design R", "o", "#D55E00"), ("cylpanel", "legacy_long_horizon", "panel, long horizon", "s", "#009E73")]
    handles, labels = [], []
    for g, d, lab, mk, col in focus:
        L = sorted([l for l in lev if l["geometry"] == g and l["design"] == d and l["quad_setting"] == "reg4_sing4"], key=lambda l: int(l["m"]))
        h = np.array([fl(l["h_max"]) for l in L]); pb = np.array([fl(l["P_bem"]) for l in L]); ps = np.array([fl(l["P_sinc_mesh"]) for l in L])
        cc = next(c for c in comp if c["geometry"] == g and c["design"] == d)
        eb, es = fl(cc["P_bem_extrap"]), fl(cc["P_consistent_sinc"])
        l1, = b.plot(h, pb / eb, "-", marker=mk, color=col, ms=4.5, mfc=col)
        l2, = b.plot(h, ps / es, "--", marker=mk, color=col, ms=4.5, mfc="white")
        b.plot([0], [1], marker="*", color=col, ms=9, mec="black", mew=0.5, ls="none")
        handles += [l1, l2]; labels += [f"BEM, {lab}", f"sinc on mesh, {lab}"]
    b.axhline(1, color=GREY, lw=0.7, ls=":")
    b.set_xlim(-0.06, 1.25); b.set_ylim(0.94, 1.6)
    b.set_xlabel("element size $h_{\\max}$ (m)"); b.set_ylabel("power / extrapolated value")
    handles.append(Line2D([], [], marker="*", ls="none", color="#999999", mec="black", mew=0.5, ms=9))
    labels.append("extrapolated ($h\\to0$)")
    b.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1, fontsize=7.2, handlelength=2.0,
             labelspacing=0.25)
    panel_tag(b, "(b)", x=-0.02)
    fig.tight_layout(w_pad=1.5)
    save(fig, "fig_bem_comparison")


# --------------------------------------------------------------------------------------------------------------
def load_paths():
    rows = rd(os.path.join(EXO, "cap_path_powers.csv"))
    lab = np.array([r["label"] for r in rows]); it = np.array([int(r["iter"]) for r in rows])
    Ps = np.array([fl(r["P_sinc_consistent"]) for r in rows]); Px = np.array([fl(r["P_exact"]) for r in rows])
    iu = np.where(lab == "paper_uniform")[0][0]
    dbs = 10 * np.log10(Ps / Ps[iu]); dbx = 10 * np.log10(Px / Px[iu])
    return lab, it, Ps, Px, dbs, dbx, rows


def fig_opt_history():
    lab, it, Ps, Px, dbs, dbx, _ = load_paths()
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 3.0), sharey=True)
    for a, names, tag in ((ax[0], ["S0", "SR", "SX"], "(a) sinc-driven"), (ax[1], ["X0", "XR", "XS"], "(b) exact-driven")):
        for n in names:
            m = (lab == n) & (it > 0)
            k = it[m]; o = np.argsort(k)
            a.plot(k[o], dbs[m][o], "-", color=RCOL[n], lw=1.2)
            a.plot(k[o], dbx[m][o], ":", color=RCOL[n], lw=1.8)
        a.set_xscale("log"); a.set_xlabel("outer iteration"); a.set_xlim(1, 700)
        a.grid(alpha=0.25, lw=0.5)
        h = [Line2D([], [], color=RCOL[n], lw=1.4) for n in names] + [Line2D([], [], color="black", lw=1.2),
                                                                       Line2D([], [], color="black", lw=1.8, ls=":")]
        a.legend(h, [RLAB[n] for n in names] + ["sinc-model power", "exact power"], loc="upper right", fontsize=7.2,
                 handlelength=2.0, labelspacing=0.25, borderaxespad=0.2)
        panel_tag(a, tag, x=0.0)
    ax[0].set_ylabel("power change vs uniform (dB)")
    ax[0].set_ylim(-32, 2)
    fig.tight_layout(w_pad=1.0)
    save(fig, "fig_opt_history")


def best_designs(lab, it, Ps, Px):
    opt = {n: np.where((lab == n) & (it == 0))[0][0] for n in RUNS}
    js = min((opt[n] for n in RUNS if n[0] == "S"), key=lambda j: Ps[j])
    jx = min((opt[n] for n in RUNS if n[0] == "X"), key=lambda j: Px[j])
    return opt, js, jx


def fig_opt_paths():
    lab, it, Ps, Px, dbs, dbx, _ = load_paths()
    opt, js, jx = best_designs(lab, it, Ps, Px)
    print("  best sinc-driven:", lab[js], dbs[js], dbx[js], " best exact-driven:", lab[jx], dbs[jx], dbx[jx])
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 3.2), gridspec_kw={"width_ratios": [1.15, 1]})
    six = {d: np.where(lab == "paper_" + d)[0][0] for d in DES}
    for pi, a in enumerate(ax):
        for n in RUNS:
            m = (lab == n) & (it > 0)
            a.plot(dbs[m], dbx[m], ".", ms=2.2 if pi == 0 else 3.0, color=RCOL[n], alpha=0.55, rasterized=True)
        a.plot(dbs[list(six.values())], dbx[list(six.values())], "s", mfc="none", mec="black", ms=5.5, mew=0.9)
        a.plot(dbs[js], dbx[js], "*", ms=12, mfc=RCOL[lab[js]], mec="black", mew=0.8)
        a.plot(dbs[jx], dbx[jx], "*", ms=12, mfc=RCOL[lab[jx]], mec="black", mew=0.8)
        a.set_xlabel("sinc-model power change vs uniform (dB)")
        a.grid(alpha=0.25, lw=0.5)
    a = ax[0]
    lo, hi = -35, 3
    a.plot([lo, hi], [lo, hi], color="black", lw=0.8, ls="--")
    a.set_xlim(-30, 2.5); a.set_ylim(-30, 2.5)
    a.set_ylabel("exact power change vs uniform (dB)")
    off = {"uniform": (4, -9), "quadvel-ones": (-3, 5), "inputpower": (-4, 7), "dyncomp": (-14, 6), "quadvel-I": (4, -9), "quadvel-rad": (5, -9)}
    for d, j in six.items():
        if d == "uniform":
            a.annotate(LAB[d], (dbs[j], dbx[j]), xytext=(-4, 5), textcoords="offset points", fontsize=7.5, ha="right")
        elif d in ("dyncomp", "inputpower"):
            continue
        else:
            a.annotate(LAB[d], (dbs[j], dbx[j]), xytext=off[d], textcoords="offset points", fontsize=7.5)
    a.annotate("C, D", (0.5 * (dbs[six['inputpower']] + dbs[six['dyncomp']]), 0.5 * (dbx[six['inputpower']] + dbx[six['dyncomp']])),
               xytext=(-14, 6), textcoords="offset points", fontsize=7.5)
    panel_tag(a, "(a)", x=0.0)
    # zoom: region within ~3 dB of the optima
    b = ax[1]
    x0, x1 = min(dbs[js], dbs[jx]) - 0.5, max(dbs[js], dbs[jx]) + 3.0
    y0, y1 = dbx[jx] - 0.5, dbx[jx] + 3.5
    b.plot([-35, 5], [-35, 5], color="black", lw=0.8, ls="--")
    b.set_xlim(x0, x1); b.set_ylim(y0, y1)
    panel_tag(b, "(b)", x=0.0)
    b.annotate("best sinc-driven", (dbs[js], dbx[js]), xytext=(10, -2), textcoords="offset points", fontsize=7.8,
               ha="left", va="top")
    b.annotate("best exact-driven", (dbs[jx], dbx[jx]), xytext=(10, -6), textcoords="offset points", fontsize=7.8,
               ha="left", va="top")
    from matplotlib.patches import Rectangle
    a.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec=GREY, lw=0.7))
    h = [Line2D([], [], marker=".", ls="none", color=RCOL[n], ms=6) for n in RUNS]
    h += [Line2D([], [], marker="s", ls="none", mfc="none", mec="black", ms=5.5),
          Line2D([], [], marker="*", ls="none", mfc="#bbbbbb", mec="black", ms=10),
          Line2D([], [], color="black", ls="--", lw=0.8)]
    l = [RLAB[n] for n in RUNS] + ["earlier designs", "best design of each driver", "predicted = true"]
    fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=4, fontsize=7.3, columnspacing=1.0,
               handletextpad=0.3, labelspacing=0.3)
    fig.tight_layout(w_pad=1.2)
    save(fig, "fig_opt_paths")


# --------------------------------------------------------------------------------------------------------------
def fig_opt_designs():
    st = sio.loadmat(os.path.join(EXO, "cap_path_states.mat"), variable_names=["T", "run", "iter", "runs", "frozen_labels"])
    T = st["T"]; run = st["run"].ravel().astype(int); it = st["iter"].ravel().astype(int)
    runs = [str(np.ravel(r)[0]) for r in st["runs"].ravel()]
    flab = [str(np.ravel(v)[0]) for v in st["frozen_labels"].ravel()]
    xs = sio.loadmat(os.path.join(EXO, "cap_samples_q8.mat"), variable_names=["x"])["x"]
    r = np.linalg.norm(xs, axis=0)
    th = np.arccos(np.clip(xs[2] / r, -1, 1)); ph = np.arctan2(xs[1], xs[0])
    th = th.reshape(-1, 64).mean(1)
    ph = np.arctan2(np.sin(ph).reshape(-1, 64).mean(1), np.cos(ph).reshape(-1, 64).mean(1))
    n1, n2 = 40, 24
    TH = np.degrees(th.reshape(n2, n1).T); PH = ph.reshape(n2, n1).T
    th_c = TH.mean(1); order = np.argsort(PH[0]); ph_c = PH[0][order]
    the = np.r_[0, 0.5 * (th_c[1:] + th_c[:-1]), 30.0]
    phe = np.r_[ph_c[0] - np.pi / n2, 0.5 * (ph_c[1:] + ph_c[:-1]), ph_c[-1] + np.pi / n2]
    ir = np.where((run == 0) & (it == -1 - flab.index("quadvel-rad")))[0][0]
    best = {n: np.where((run == runs.index(n) + 1) & (it == 0))[0][0] for n in ("S0", "X0", "SX")}
    items = [("design R", ir), ("sinc,\nfrom uniform", best["S0"]), ("exact,\nfrom uniform", best["X0"]),
             ("sinc, from\nexact optimum", best["SX"])]
    fig, axs = plt.subplots(1, 4, figsize=(6.5, 1.95), subplot_kw=dict(projection="polar"))
    for ax, (ttl, j) in zip(axs, items):
        t = T[:, j].reshape(n2, n1).T[:, order]
        pc = ax.pcolormesh(phe, the, t * 1000, cmap="viridis", vmin=25, vmax=100, shading="flat")
        ax.set_title(ttl, fontsize=8.5, pad=3)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_ylim(0, 30)
        ax.set_theta_zero_location("N"); ax.grid(False)
        ax.spines["polar"].set_linewidth(0.5)
    fig.subplots_adjust(left=0.005, right=0.90, top=0.80, bottom=0.01, wspace=0.04)
    cax = fig.add_axes([0.925, 0.12, 0.018, 0.66])
    cb = fig.colorbar(pc, cax=cax); cb.set_label("thickness (mm)", fontsize=8.5); cb.ax.tick_params(labelsize=8)
    print("  thickness range (mm):", T[:, [ir, best['S0'], best['X0'], best['SX']]].min() * 1000, T[:, [ir, best['S0'], best['X0'], best['SX']]].max() * 1000)
    save(fig, "fig_opt_designs")


# --------------------------------------------------------------------------------------------------------------
def fig_ranking():
    rk = rd(os.path.join(NOV, "cap_ranking.csv"))
    kR = np.array([fl(r["kR"]) for r in rk]); tau = np.array([fl(r["kendall_tau"]) for r in rk])
    fig, ax = plt.subplots(figsize=(4.0, 2.9))
    ax.plot(kR, tau, "-o", color="#0072B2", ms=3.5, lw=1.1)
    ax.axvline(KR100, color=GREY, lw=0.9, ls=":")
    ax.text(KR100 * 1.05, 0.42, "100 Hz", fontsize=8.5, ha="left", va="bottom")
    ax.set_xscale("log"); ax.set_xlim(0.1, 40); ax.set_ylim(0.4, 1.05)
    ax.set_xticks([0.1, 1, 10]); ax.set_xticklabels(["0.1", "1", "10"])
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_xlabel("$kR$"); ax.set_ylabel("Kendall $\\tau$, sinc vs exact (six designs)")
    ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout()
    save(fig, "fig_ranking")


ALL = ["fig_sphere_modes", "fig_validation", "fig_ratio_kR", "fig_regime_map", "fig_face_split", "fig_directivity",
       "fig_power_freq", "fig_bem_comparison", "fig_opt_history", "fig_opt_paths", "fig_opt_designs", "fig_ranking"]

if __name__ == "__main__":
    todo = sys.argv[1:] or ALL
    for s in todo:
        print(s)
        globals()[s]()





