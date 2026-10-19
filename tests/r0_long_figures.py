"""Figures for the r0 long-run study (part A: restarts from pilot; part B: strategies).

Run: C:/Users/dshgu/miniconda3/python.exe r0_long_figures.py
Skips missing cells/runs silently (prints a summary).
"""
import os
import re
import glob
import numpy as np
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = r"C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course\results\r0_revision"
LONG = os.path.join(ROOT, "longruns")
PILOT = os.path.join(ROOT, "pilot")

CELLS = ["c30_k3", "c90_k1p5", "c120_k1p5", "c30_k1", "c30_k18", "c90_k7", "c120_k7"]
LABEL = {"c30_k3": "30°, kR 3", "c90_k1p5": "90°, kR 1.5", "c120_k1p5": "120°, kR 1.5",
         "c30_k1": "30°, kR 1", "c30_k18": "30°, kR 18.3", "c90_k7": "90°, kR 7",
         "c120_k7": "120°, kR 7"}
RUNS = ["XU", "XS", "SU", "SX"]
RUN_LABEL = {"XU": "XU (exact from uniform)", "XS": "XS (exact from sinc opt.)",
             "SU": "SU (sinc from uniform)", "SX": "SX (sinc from exact opt.)"}
# Okabe-Ito
RUN_COL = {"XU": "#0072B2", "XS": "#D55E00", "SU": "#009E73", "SX": "#CC79A7"}
STRATS = ["P", "S1", "S2", "S3", "S4", "S5"]
STRAT_LABEL = {"P": "P plain exact", "S1": "S1 ERP>sinc>exact", "S2": "S2 sinc>exact",
               "S3": "S3 sinc>exact+trust region", "S4": "S4 kR continuation",
               "S5": "S5 switching rule"}
STRAT_COL = {"P": "#000000", "S1": "#E69F00", "S2": "#56B4E9", "S3": "#009E73",
             "S4": "#D55E00", "S5": "#CC79A7"}
RESTART = 600
PLAIN_END = 1200   # plain strategy P (stopped by its rule) is continued flat to here
VOL_TOL = 0.5 * (1 + 1e-5)


def feasible(vol):
    v = np.asarray(vol, float)
    return v <= VOL_TOL   # NaN -> False

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    "font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11, "legend.fontsize": 9.5,
    "axes.grid": True, "grid.color": "#d9d9d9", "grid.linewidth": 0.7,
    "axes.spines.top": False, "axes.spines.right": False,
})


def db(p, p0):
    return 10.0 * np.log10(np.asarray(p, float) / float(p0))


def read_csv(f):
    """Minimal CSV reader (no pandas): dict col -> array (float where possible, else str)."""
    with open(f, newline="") as fh:
        rows = list(csv.reader(fh))
    if len(rows) < 2:
        return {}
    hdr, body = rows[0], [r for r in rows[1:] if len(r) == len(rows[0])]
    out = {}
    for j, h in enumerate(hdr):
        col = [r[j] for r in body]
        try:
            out[h] = np.array([float(x) if x.strip() != "" else np.nan for x in col])
        except ValueError:
            out[h] = np.array(col, dtype=object)
    return out


def nrows(d):
    return len(next(iter(d.values()))) if d else 0


def dedup_last(it, *cols):
    """sort by it (stable) and keep the last row of each duplicated it."""
    o = np.argsort(it, kind="stable")
    it = it[o]
    keep = np.ones(len(it), bool)
    keep[:-1] = it[:-1] != it[1:]
    return (it[keep],) + tuple(c[o][keep] for c in cols)


# ---------------------------------------------------------------- loading
def load_pilot(cell):
    f = os.path.join(PILOT, "eval_%s.csv" % cell)
    if not os.path.isfile(f):
        return None
    try:
        d = read_csv(f)
        return d if d else None
    except Exception as e:
        print("  could not read", f, e)
        return None


def _uniform(pdf):
    m = pdf["set"] == "uniform"
    if not m.any():
        return None
    i = np.where(m)[0][0]
    return {k: pdf[k][i] for k in ("P_exact", "P_sinc")}


def pilot_all_exact_dB(pdf):
    if pdf is None:
        return np.array([])
    u = _uniform(pdf)
    if u is None:
        return np.array([])
    v = pdf["P_exact"]
    ok = np.isfinite(v) & feasible(pdf["vol"])
    return db(v[ok], u["P_exact"])


def pilot_run(pdf, run):
    """dict it, own_dB, exact_dB, sinc_dB for pilot iterations (path + optimum as 600)."""
    if pdf is None:
        return None
    u = _uniform(pdf)
    if u is None:
        return None
    runcol = pdf["run"]
    mp = (pdf["set"] == "path") & (runcol == run)
    mo = (pdf["set"] == "optimum") & (runcol == run)
    if not (mp.any() or mo.any()):
        return None
    it = np.concatenate([pdf["iter"][mp].astype(int), np.full(int(mo.sum()), RESTART)])
    pe = np.concatenate([pdf["P_exact"][mp], pdf["P_exact"][mo]])
    ps = np.concatenate([pdf["P_sinc"][mp], pdf["P_sinc"][mo]])
    vo = np.concatenate([pdf["vol"][mp], pdf["vol"][mo]])
    # path rows first, optimum last -> optimum wins at it==600
    it, pe, ps, vo = dedup_last(it, pe, ps, vo)
    own = pe if run[0] == "X" else ps
    p0 = u["P_exact"] if run[0] == "X" else u["P_sinc"]
    return {"it": it, "own_dB": db(own, p0), "exact_dB": db(pe, u["P_exact"]),
            "sinc_dB": db(ps, u["P_sinc"]), "vol": vo}


def _load_hist(prefix, pat):
    out = {}
    for f in sorted(glob.glob(os.path.join(LONG, prefix + "_*.csv"))):
        m = re.match(pat, os.path.basename(f))
        if not m:
            continue
        try:
            d = read_csv(f)
        except Exception as e:
            print("  could not read", f, e)
            continue
        if nrows(d) == 0:
            continue
        out.setdefault(m.group(1), {})[m.group(2)] = d
    return out


def load_hist_A():
    return _load_hist("hist_A", r"hist_A_(.+)_(XU|XS|SU|SX)\.csv$")


def load_hist_B():
    return _load_hist("hist_B", r"hist_B_(.+)_(P|S\d)\.csv$")


def combined(pil, long_d):
    keys = ("own_dB", "exact_dB", "sinc_dB", "vol")
    parts = [p for p in (pil, long_d) if p is not None]
    if not parts:
        return None
    it = np.concatenate([np.asarray(p["it"], float) for p in parts])
    cols = [np.concatenate([np.asarray(p[k], float) for p in parts]) for k in keys]
    res = dedup_last(it, *cols)
    return dict(zip(("it",) + keys, res))


def pick(d):
    """Running argmin of own_dB -> exact_dB and sinc_dB of the pick."""
    own = d["own_dB"]
    ok = feasible(d["vol"])
    best, bi = np.inf, -1
    idx = np.zeros(len(own), int)
    for i, v in enumerate(own):
        if ok[i] and np.isfinite(v) and v < best:
            best, bi = v, i
        idx[i] = bi
    ex = np.where(idx >= 0, d["exact_dB"][np.clip(idx, 0, None)], np.nan)
    sn = np.where(idx >= 0, d["sinc_dB"][np.clip(idx, 0, None)], np.nan)
    return {"it": d["it"], "exact": ex, "sinc": sn}


def feas_cummin(ex, vol):
    """Running minimum of exact dB over feasible iterates only (NaN before the first)."""
    return np.fmin.accumulate(np.where(feasible(vol), ex, np.nan))


def ffill_on(grid, it, val):
    """value of step function (last it' <= N) on grid; NaN before first."""
    o = np.argsort(it, kind="stable")
    it, val = np.asarray(it)[o], np.asarray(val)[o]
    k = np.searchsorted(it, grid, side="right") - 1
    out = np.where(k >= 0, val[np.clip(k, 0, None)], np.nan)
    return out

def grid_axes(n, w=5.2, h=3.9, **kw):
    nc = min(n, 3)
    nr = int(np.ceil(n / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(w * nc, h * nr), squeeze=False, **kw)
    flat = list(axs.ravel())
    for a in flat[n:]:
        a.set_visible(False)
    return fig, flat[:n]


def save(fig, name):
    for ext in ("png", "pdf"):
        p = os.path.join(LONG, "%s.%s" % (name, ext))
        fig.savefig(p, dpi=160, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print("  wrote", os.path.join(LONG, name + ".png/.pdf"))


def restart_line(ax, label=False):
    ax.axvline(RESTART, color="grey", ls=":", lw=1.2)
    if label:
        ax.text(RESTART, 1.0, " restart", transform=ax.get_xaxis_transform(),
                color="grey", fontsize=9, va="top", ha="left")


# ---------------------------------------------------------------- figures
def fig_convergence(HA, pil):
    cells = [c for c in CELLS if c in HA]
    if not cells:
        return
    fig, axs = grid_axes(len(cells))
    for ax, c in zip(axs, cells):
        for r in RUNS:
            d = combined(pilot_run(pil.get(c), r), HA[c].get(r))
            if d is None or (r not in HA[c] and pilot_run(pil.get(c), r) is None):
                continue
            pk = pick(d)
            ax.plot(pk["it"], pk["exact"], color=RUN_COL[r], lw=1.8, label=RUN_LABEL[r] + ", pick exact")
            if r[0] == "S":
                ax.plot(pk["it"], pk["sinc"], color=RUN_COL[r], lw=1.4, ls="--",
                        label=r + " pick, sinc prediction")
        restart_line(ax, True)
        ax.set_title(LABEL[c])
        ax.set_xlabel("iteration")
        ax.set_ylabel("dB vs uniform")
    h, l = [], []
    for a in axs:
        for x, y in zip(*a.get_legend_handles_labels()):
            if y not in l:
                h.append(x); l.append(y)
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.04))
    fig.suptitle("Convergence of the running-best pick (solid: exact model, dashed: sinc prediction)", y=1.0)
    fig.tight_layout(rect=(0, 0.08, 1, 0.98))
    save(fig, "long_convergence_A")


def fig_regret(HA, pil):
    cells = [c for c in CELLS if c in HA and ("SU" in HA[c] or ("XU" in HA[c] and "XS" in HA[c]))]
    if not cells:
        return
    fig, axs = grid_axes(len(cells))
    for ax, c in zip(axs, cells):
        runs = HA[c]
        grid = np.unique(np.concatenate([runs[r]["it"] for r in runs]))
        grid = grid[grid >= RESTART + 1]
        pm = pilot_all_exact_dB(pil.get(c))
        best = np.full(len(grid), np.nanmin(pm) if pm.size else np.inf)
        for r in runs:
            cm = feas_cummin(runs[r]["exact_dB"], runs[r]["vol"])
            v = ffill_on(grid, runs[r]["it"], cm)
            best = np.fmin(best, np.where(np.isnan(v), np.inf, v))
        picks = {}
        for r in ("SU", "XU", "XS"):
            if r in runs:
                pk = pick(combined(pilot_run(pil.get(c), r), runs[r]))
                picks[r] = ffill_on(grid, pk["it"], pk["exact"])
        if "SU" in picks:
            ax.plot(grid, picks["SU"] - best, color=RUN_COL["SU"], lw=1.8,
                    label="sinc regret lower bound (SU pick − best known)")
        if "XU" in picks and "XS" in picks:
            ax.plot(grid, np.abs(picks["XU"] - picks["XS"]), color="#444444", lw=1.6,
                    ls="--", label="|XU pick − XS pick| (start-to-start spread)")
        ax.set_title(LABEL[c])
        ax.set_xlabel("iteration N")
        ax.set_ylabel("dB")
        ax.set_ylim(bottom=min(0, ax.get_ylim()[0]))
        if not ax.get_legend_handles_labels()[0]:
            ax.text(0.5, 0.5, "need SU, or XU and XS", transform=ax.transAxes, ha="center", color="grey")
    h, l = [], []
    for a in axs:
        for x, y in zip(*a.get_legend_handles_labels()):
            if y not in l:
                h.append(x); l.append(y)
    if h:
        fig.legend(h, l, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.03))
    fig.suptitle("Sinc regret lower bound and start-to-start spread (exact model)", y=1.0)
    fig.tight_layout(rect=(0, 0.07, 1, 0.98))
    save(fig, "long_regret_A")


def fig_strategies(HB, HA, pil):
    cells = [c for c in CELLS if c in HB]
    if not cells:
        return
    fig, axs = plt.subplots(len(cells), 2, figsize=(11.5, 4.4 * len(cells)), squeeze=False)
    for row, c in zip(axs, cells):
        axl, axr = row
        allv, any_inf = [], False
        for s in STRATS:
            df = HB[c].get(s)
            if df is None:
                continue
            ex = df["exact_dB"]
            rm = feas_cummin(ex, df["vol"])
            ru = np.fmin.accumulate(ex)
            it = df["it"]
            tm = np.cumsum(np.nan_to_num(df["time_s"])) / 60.0
            stg = df["stage"]
            ch = np.where(stg[1:] != stg[:-1])[0] + 1
            kw = dict(color=STRAT_COL[s], lw=1.7)
            axl.plot(it, rm, label=STRAT_LABEL[s], **kw)
            axr.plot(tm, rm, **kw)
            allv.append(rm)
            if np.any(~feasible(df["vol"])):
                kd = dict(color=STRAT_COL[s], lw=0.9, ls=":")
                axl.plot(it, ru, **kd)
                axr.plot(tm, ru, **kd)
                any_inf = True
            if len(ch):
                mk = dict(ls="none", marker="o", ms=5, mfc="white", mec=STRAT_COL[s], mew=1.4)
                good = ch[np.isfinite(rm[ch])]
                axl.plot(it[good], rm[good], **mk)
                axr.plot(tm[good], rm[good], **mk)
            if s == "P" and it[-1] < PLAIN_END and np.isfinite(rm[-1]):
                axl.plot([it[-1], PLAIN_END], [rm[-1], rm[-1]], color=STRAT_COL[s], lw=1.7, ls="--")
                axl.plot([it[-1]], [rm[-1]], marker="s", ms=6, ls="none", color=STRAT_COL[s],
                         label="P converged (stopped at %d)" % it[-1])
        if "XU" in HA.get(c, {}):
            d = combined(pilot_run(pil.get(c), "XU"), HA[c]["XU"])
            if d is not None:
                g = feas_cummin(d["exact_dB"], d["vol"])
                axl.plot(d["it"], g, color="#999999", lw=1.0, label="XU pilot + restart at 600")
                allv.append(g)
        vals = np.concatenate(allv) if allv else np.array([])
        vals = vals[np.isfinite(vals)]
        if vals.size:
            lo = vals.min() - 3
            hi = max(0.0, vals.max()) + 2
            axl.set_ylim(lo, hi); axr.set_ylim(lo, hi)
        dot = "  (dotted: ignoring the volume constraint)" if any_inf else ""
        axl.set_title("%s: running minimum vs iteration%s" % (LABEL[c], ""), fontsize=11)
        axr.set_title("%s: running minimum vs wall time" % LABEL[c], fontsize=11)
        if any_inf:
            axr.plot([], [], color="#555555", lw=0.9, ls=":", label="dotted: ignoring the volume constraint")
            axr.legend(loc="upper right", frameon=True, framealpha=0.9, facecolor="white",
                       edgecolor="#cccccc", fontsize=8.5)
        axl.set_xlabel("iteration"); axr.set_xlabel("cumulative wall time (min)")
        axl.set_ylabel("best exact dB vs uniform"); axr.set_ylabel("best exact dB vs uniform")
        axl.legend(loc="upper right", frameon=True, framealpha=0.9, facecolor="white", edgecolor="#cccccc")
    fig.text(0.5, 0.0, "solid: best feasible (vol <= 0.5(1+1e-5)); open circles: stage changes; "
             "dotted: ignoring the volume constraint; dashed: P continued flat after stopping",
             ha="center", fontsize=9, color="#444444")
    fig.tight_layout(rect=(0, 0.015, 1, 1))
    save(fig, "long_strategies_B")


def fig_gap(HA, pil):
    cells = [c for c in CELLS if c in HA and any(r in HA[c] for r in ("SU", "SX"))]
    if not cells:
        return
    fig, axs = grid_axes(len(cells))
    for ax, c in zip(axs, cells):
        for r in ("SU", "SX"):
            if r not in HA[c]:
                continue
            d = combined(pilot_run(pil.get(c), r), HA[c][r])
            ax.plot(d["it"], d["sinc_dB"] - d["exact_dB"], color=RUN_COL[r], lw=1.5, label=RUN_LABEL[r])
        ax.axhline(0, color="#888888", lw=0.8)
        restart_line(ax, True)
        ax.set_title(LABEL[c])
        ax.set_xlabel("iteration")
        ax.set_ylabel("sinc dB − exact dB")
        ax.legend(loc="best", frameon=True, framealpha=0.9, facecolor="white", edgecolor="#cccccc")
    fig.suptitle("Gap between sinc and exact prediction along the sinc runs", y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    save(fig, "long_gap_A")


def main():
    HA, HB = load_hist_A(), load_hist_B()
    pil = {c: load_pilot(c) for c in CELLS}
    print("Part A cells/runs found:")
    for c in CELLS:
        if c in HA:
            print("  %-10s %s" % (c, ", ".join("%s(%d rows, it %d..%d)" % (
                r, nrows(d), d["it"].min(), d["it"].max()) for r, d in sorted(HA[c].items()))))
    if not HA:
        print("  none")
    print("Part B cells/strategies found:")
    for c in CELLS:
        if c in HB:
            print("  %-10s %s" % (c, ", ".join("%s(%d rows)" % (s, nrows(d)) for s, d in sorted(HB[c].items()))))
    if not HB:
        print("  none")
    print("Pilot eval files missing:", [c for c in CELLS if pil[c] is None] or "none")
    fig_convergence(HA, pil)
    fig_regret(HA, pil)
    fig_strategies(HB, HA, pil)
    fig_gap(HA, pil)


if __name__ == "__main__":
    main()

