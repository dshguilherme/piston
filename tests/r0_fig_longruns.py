"""Plot data for the JSV paper figure F27 (`fig_longruns`, Section 4.4 of manuscript_eng.tex):
level against iteration for the continued runs of Experiment 8, two panels.

  (a) design point (30 deg cap, kR = 18.3): exact-driven run X0 (pilot run c30_k18 XU), 1 to 4355 iterations
  (b) 90 deg cap, kR = 1.5: exact runs from uniform (XU, stopped at 944) and from the sinc design (XS, 3000),
      sinc run from uniform (SU, 3000): its true (exact) and its predicted (sinc) level

Level convention of the paper (Sec. 3): dL = 10 log10(P(d)/P(d0)) in dB, d0 = uniform design of the same
case and operator; 0 = uniform, negative = quieter. Every value is the running "pick" of the run: the
feasible iterate (vol <= 0.5 (1 + 1e-5)) with the lowest own objective so far (the design the optimizer
returns, Sec. 2.5); the exact and sinc levels are those of the pick. Iterations 1-600 come from the
regime-study evaluation (every fifth iterate plus the returned design at 600), 601 onwards from the long-run
history (every iteration; thinned here to every fifth plus the last). Counts include the first 600.

Inputs (read only): results/r0_revision/pilot/eval_<cell>.csv, results/r0_revision/longruns/hist_A_<cell>_<run>.csv
Outputs (new files, paper_jsv/prototype/figs/data/):
  f27_<panel>_<run>.dat   columns: it dLx dLs  (dLs = sinc level of the pick; dL of the own objective = dLx
                          for exact runs, dLs for the sinc run)
  f27_numbers.tex         TeX macros used by the figure (end iterations and levels), checked against
                          longruns/long_runs_A.csv below
Run:  python tests/r0_fig_longruns.py      (numpy only; no solver)
"""
import csv
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PIL = os.path.join(ROOT, "results", "r0_revision", "pilot")
LONG = os.path.join(ROOT, "results", "r0_revision", "longruns")
OUT = os.path.join(ROOT, "paper_jsv", "prototype", "figs", "data")
HEAD = "% Do not edit by hand: written by tests/r0_fig_longruns.py from results/r0_revision/{pilot,longruns}/\n"
VOLTOL = 0.5 * (1 + 1e-5)


def rcsv(f):
    with open(f, newline="") as fh:
        return list(csv.DictReader(fh))


def num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return np.nan


def pilot_part(cell, run):
    rows = rcsv(os.path.join(PIL, "eval_%s.csv" % cell))
    uni = [r for r in rows if r["set"] == "uniform"][0]
    pe0, ps0 = num(uni["P_exact"]), num(uni["P_sinc"])
    path = [r for r in rows if r["set"] == "path" and r["run"] == run]
    opt = [r for r in rows if r["set"] == "optimum" and r["run"] == run]
    d = {}
    for r in path:
        d[int(r["iter"])] = r
    for r in opt:          # the returned design is stored at iteration 600
        d[600] = r
    it = np.array(sorted(d))
    pe = np.array([num(d[i]["P_exact"]) for i in it])
    ps = np.array([num(d[i]["P_sinc"]) for i in it])
    vo = np.array([num(d[i]["vol"]) for i in it])
    return it, 10 * np.log10(pe / pe0), 10 * np.log10(ps / ps0), vo


def long_part(cell, run):
    f = os.path.join(LONG, "hist_A_%s_%s.csv" % (cell, run))
    rows = rcsv(f)
    it = np.array([int(r["it"]) for r in rows])
    ex = np.array([num(r["exact_dB"]) for r in rows])
    si = np.array([num(r["sinc_dB"]) for r in rows])
    vo = np.array([num(r["vol"]) for r in rows])
    # a restarted chunk may repeat an iteration: keep the last row of each
    keep = np.r_[it[:-1] != it[1:], True]
    return it[keep], ex[keep], si[keep], vo[keep]


def picks(cell, run):
    a = pilot_part(cell, run)
    # iteration 1 is the start design (as in the long-run histories): the uniform design (0 dB) for runs from
    # uniform, the returned design of the sinc run from uniform (its pick at 600) for XS
    if run in ("XU", "SU"):
        s = (1, 0.0, 0.0, 0.5)
    else:
        su = pilot_part(cell, "SU")
        k = int(np.where(su[0] == 600)[0][0])
        s = (1, su[1][k], su[2][k], su[3][k])
    if a[0][0] > 1:
        a = tuple(np.r_[s[j], a[j]] for j in range(4))
    b = long_part(cell, run)
    it = np.r_[a[0], b[0]]
    ex = np.r_[a[1], b[1]]
    si = np.r_[a[2], b[2]]
    vo = np.r_[a[3], b[3]]
    own = ex if run[0] == "X" else si
    best = np.inf
    px, ps = np.nan, np.nan
    out = []
    for k in range(len(it)):
        if vo[k] <= VOLTOL and np.isfinite(own[k]) and own[k] < best:
            best, px, ps = own[k], ex[k], si[k]
        if np.isfinite(best):
            out.append((it[k], px, ps))
    out = np.array(out)
    # thin the long part to every fifth iteration, keep the first 600 as stored and the last row
    m = (out[:, 0] <= 600) | (out[:, 0] % 5 == 0)
    m[-1] = True
    return out[m]


def write(name, arr):
    with open(os.path.join(OUT, name), "w", newline="\n") as fh:
        fh.write(HEAD)
        fh.write("it dLx dLs\n")
        for r in arr:
            fh.write("%d %.4f %.4f\n" % (r[0], r[1], r[2]))


def main():
    runs = {("a", "X0"): ("c30_k18", "XU"),
            ("b", "XU"): ("c90_k1p5", "XU"),
            ("b", "XS"): ("c90_k1p5", "XS"),
            ("b", "SU"): ("c90_k1p5", "SU")}
    ref = {(r["cell"], r["run"], int(r["N"])): r for r in rcsv(os.path.join(LONG, "long_runs_A.csv"))}
    ends = {}
    for (panel, tag), (cell, run) in runs.items():
        p = picks(cell, run)
        write("f27_%s_%s.dat" % (panel, tag), p)
        last = p[-1]
        ends[(panel, tag)] = last
        # check the last pick and the pick at 600 against the Experiment 8 summary table
        cand = [v for (c, rr, _), v in ref.items() if c == cell and rr == run]
        r = max(cand, key=lambda v: int(v["N"])) if cand else None   # last checkpoint of the run
        col = "exact_dB" if run[0] == "X" else "sinc_dB"
        stored = num(r[col]) if r else np.nan
        mine = last[1] if run[0] == "X" else last[2]
        print("%-10s %-3s last it %5d  dLx %.4f  dLs %.4f  summary %s %.4f  diff %.1e" %
              (cell, run, last[0], last[1], last[2], col, stored, abs(mine - stored)))
        assert abs(mine - stored) < 1e-3, "mismatch with long_runs_A.csv"
        at600 = p[p[:, 0] == 600][0]
        r6 = ref.get((cell, run, 600))
        s6 = num(r6["exact_dB"]) if run[0] == "X" else num(r6["sinc_dB"])
        m6 = at600[1] if run[0] == "X" else at600[2]
        assert abs(m6 - s6) < 1e-3, "mismatch at 600"
    with open(os.path.join(OUT, "f27_numbers.tex"), "w", newline="\n") as fh:
        fh.write(HEAD)
        def mac(name, val, fmt):
            fh.write("\\def\\%s{%s}\n" % (name, fmt % val))
        mac("lrXzEndIt", ends[("a", "X0")][0], "%d")
        mac("lrXzEnd", ends[("a", "X0")][1], "%.2f")
        a = ends[("a", "X0")]
        p = picks("c30_k18", "XU")
        mac("lrXzSix", p[p[:, 0] == 600][0][1], "%.2f")
        mac("lrXUEndIt", ends[("b", "XU")][0], "%d")
        mac("lrXUEnd", ends[("b", "XU")][1], "%.1f")
        mac("lrXSEnd", ends[("b", "XS")][1], "%.1f")
        mac("lrSUEndTrue", ends[("b", "SU")][1], "%.1f")
        mac("lrSUEndPred", ends[("b", "SU")][2], "%.1f")
    print("written to", OUT)


if __name__ == "__main__":
    main()
