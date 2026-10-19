"""Robust-optimisation figures (fragility, FSI fragility, histories, N study).
Re-runnable: picks up *_robust.csv and new hist_*.csv files when present."""
import csv, glob, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = r"C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course\results\r0_revision\robust"
written = []
plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "legend.fontsize": 8.5,
                     "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white"})

CELLS = {"c30_k18": "paper case, 30 deg cap, kR 18.3", "c30_k3": "scaled case, 30 deg cap, kR 3"}
DESIGNS = {"c30_k18": ["S0", "SX", "X0", "XS"],
           "c30_k3": ["XU", "XS", "SU", "SX", "Fdev", "LXU_it1568"]}
BLUES = ["#08306b", "#2171b5", "#6baed6", "#41b6c4", "#9ecae1"]
ORANGES = ["#d94801", "#e6550d", "#fdae6b", "#a63603", "#fd8d3c"]
RP = ["#238b45", "#6a3d9a", "#74c476", "#9e9ac8", "#00441b", "#54278f"]


def is_exact(n):
    return n[0] in "XF" or n.startswith("LX")


def load(path):
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            rows.append(r)
    return rows


def col(rows, name):
    return np.array([float(r[name]) if r.get(name, "") not in ("", None) else np.nan for r in rows])


def read_cell(cell):
    """Return dict design -> {col -> {sigma -> array}} merged from main + robust files."""
    out, ref = {}, {}
    for suffix, tag in (("", "main"), ("_robust", "rob")):
        p = os.path.join(D, f"evalstudy_{cell}{suffix}.csv")
        if not os.path.exists(p):
            continue
        rows = load(p)
        des = np.array([r["design"] for r in rows])
        sig = col(rows, "sigma")
        cols = [c for c in rows[0].keys() if c.startswith("P_")]
        data = {c: col(rows, c) for c in cols}
        if tag == "main":
            m = (des == "uniform") & (sig == 0)
            ref = {c: data[c][m][0] for c in cols}
        for d in dict.fromkeys(des):
            if tag == "rob" and (d == "uniform" or not d.startswith("R")):
                continue
            out.setdefault(d, {})
            for c in cols:
                out[d][c] = {}
                for s in np.unique(sig[des == d]):
                    out[d][c][s] = data[c][(des == d) & (sig == s)]
    return out, ref


def save(fig, name):
    for ext in ("png", "pdf"):
        p = os.path.join(D, f"{name}.{ext}")
        fig.savefig(p, dpi=200, facecolor="white", bbox_inches="tight")
        written.append(p)
    plt.close(fig)


XNOM = 0.05
SIGS = [0.001, 0.003, 0.01, 0.03]


def draw_design(ax, vals, ref, color, label, ls="-", marker="o", band=True):
    """vals: sigma -> samples. Plots median, band, mean-dashed, nominal marker."""
    xs, med, lo, hi, mn = [], [], [], [], []
    for s in SIGS:
        v = vals.get(s)
        if v is None:
            continue
        v = v[np.isfinite(v)]
        if v.size == 0:
            continue
        dB = 10 * np.log10(v / ref)
        xs.append(s * 100)
        med.append(np.median(dB)); lo.append(np.percentile(dB, 5)); hi.append(np.percentile(dB, 95))
        mn.append(10 * np.log10(np.mean(v) / ref))
    if xs:
        ax.plot(xs, med, ls=ls, marker=marker, color=color, ms=4.5, lw=1.6, label=label)
        if band:
            ax.fill_between(xs, lo, hi, color=color, alpha=0.13, lw=0)
        ax.plot(xs, mn, ls="--", color=color, lw=1.0, alpha=0.9)
    v0 = vals.get(0.0)
    if v0 is not None and np.isfinite(v0[0]):
        ax.plot([XNOM], [10 * np.log10(v0[0] / ref)], marker=marker, color=color, ms=6,
                mfc="white" if ls != "-" else color, ls="none")


def style(ax, title):
    ax.set_xscale("log")
    ax.set_xticks([XNOM, 0.1, 0.3, 1, 3])
    ax.set_xticklabels(["nominal", "0.1", "0.3", "1", "3"])
    ax.minorticks_off()
    ax.set_xlim(0.035, 4.5)
    ax.axhline(-20, color="k", ls=":", lw=1)
    ax.set_xlabel(r"relative thickness error $\sigma$ (%)")
    ax.set_title(title)
    ax.grid(alpha=0.25)


def legend_below(fig, axes, extra_note=True):
    h, l = [], []
    for ax in axes:
        hh, ll = ax.get_legend_handles_labels()
        h += hh; l += ll
    return h, l


def fragility(column, fname, fsi=False):
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.6), sharey=False)
    for ax, cell in zip(axes, CELLS):
        data, ref = read_cell(cell)
        if not data:
            continue
        names = [d for d in DESIGNS[cell] if d in data]
        names += [d for d in data if d.startswith("R") and d not in names]
        ie = ir = 0
        if not fsi:
            for d in names:
                if d.startswith("R"):
                    c = RP[ir % len(RP)]; ir += 1
                    lab = re.sub(r"_it\d+$", "", d)
                elif is_exact(d):
                    c = BLUES[ie % len(BLUES)]; ie += 1; lab = d
                else:
                    c = ORANGES[(len(names) - 1 - ie) % len(ORANGES)] if False else None
                    lab = d
                    c = ORANGES[sum(1 for x in names[:names.index(d)] if not is_exact(x) and not x.startswith("R")) % len(ORANGES)]
                if column in data[d]:
                    draw_design(ax, data[d][column], ref[column], c, lab)
        else:
            specs = [("P_exact_fsi_dev", "-", "o", "dev")]
            if cell == "c30_k3":
                specs.append(("P_exact_fsi_phys", "-.", "s", "phys"))
            ie = ir = io = 0
            for d in names:
                if d.startswith("R"):
                    c = RP[ir % len(RP)]; ir += 1
                elif is_exact(d):
                    c = BLUES[ie % len(BLUES)]; ie += 1
                else:
                    c = ORANGES[io % len(ORANGES)]; io += 1
                for cn, ls, mk, tag in specs:
                    if cn in data[d] and cn in ref:
                        lab = re.sub(r"_it\d+$", "", d) if d.startswith("R") else d
                        draw_design(ax, data[d][cn], ref[cn], c, f"{lab} ({tag})", ls=ls, marker=mk)
        style(ax, f"{cell}: {CELLS[cell]}")
        ax.set_ylabel("reduction vs nominal uniform (dB)")
        ax.text(3.9, -20, "20 dB target", ha="right", va="bottom", fontsize=9)
        lg = ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, title="median, 5-95% band;\ndashed = mean")
        lg.get_title().set_fontsize(8)
    fig.tight_layout()
    save(fig, fname)


def history():
    files = sorted(glob.glob(os.path.join(D, "hist_*.csv")))
    bycell = {}
    for f in files:
        m = re.match(r"hist_(c\d+_k\d+)_(.*)\.csv", os.path.basename(f))
        if m:
            bycell.setdefault(m.group(1), []).append((m.group(2), f))
    if not bycell:
        return
    cells = sorted(bycell)
    fig, axes = plt.subplots(1, len(cells), figsize=(6.8 * len(cells), 5.4), squeeze=False)
    cmap = plt.get_cmap("tab10")
    for ax, cell in zip(axes[0], cells):
        for i, (lab, f) in enumerate(bycell[cell]):
            r = load(f)
            it = col(r, "it")
            c = cmap(i % 10)
            ax.plot(it, col(r, "rob_dB"), "-", color=c, lw=1.6, label=lab)
            ax.plot(it, col(r, "nom_exact_dB"), "--", color=c, lw=1.0)
        ax.set_title(f"{cell}: {CELLS.get(cell, '')}")
        ax.set_xlabel("iteration"); ax.set_ylabel("dB vs nominal uniform")
        ax.grid(alpha=0.25)
        ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False,
                  title="solid: robust objective\ndashed: nominal exact")
    fig.tight_layout()
    save(fig, "robust_history")


def nstudy():
    files = sorted(glob.glob(os.path.join(D, "hist_c30_k3_*XX_N*.csv")))
    if not files:
        return
    def nval(f):
        return int(re.search(r"XX_N(\d+)", f).group(1))
    files.sort(key=nval)
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    cmap = plt.get_cmap("viridis")
    for i, f in enumerate(files):
        r = load(f); it = col(r, "it"); c = cmap(i / max(1, len(files) - 1) * 0.85)
        ax.plot(it, col(r, "rob_dB"), "-", color=c, lw=1.7, label=f"N={nval(f)} robust objective")
        ax.plot(it, col(r, "nom_exact_dB"), "--", color=c, lw=1.0, label=f"N={nval(f)} nominal exact")
    data, ref = read_cell("c30_k3")
    if "XU" in data and 0.01 in data["XU"].get("P_exact", {}):
        y = 10 * np.log10(np.mean(data["XU"]["P_exact"][0.01]) / ref["P_exact"])
        ax.axhline(y, color="0.45", lw=1.4)
        ax.text(0.99, y, "nominal XU, 200 fresh samples", transform=ax.get_yaxis_transform(),
                ha="right", va="bottom", color="0.3", fontsize=9)
    ax.set_xlabel("iteration"); ax.set_ylabel("dB vs nominal uniform")
    ax.set_title("c30_k3 sample-count study (sigma 1%)")
    ax.grid(alpha=0.25)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
    fig.tight_layout()
    save(fig, "robust_Nstudy")


fragility("P_exact", "robust_fragility")
fragility(None, "robust_fragility_fsi", fsi=True)
history()
nstudy()
print("Files written:")
for p in written:
    print(" ", p)
