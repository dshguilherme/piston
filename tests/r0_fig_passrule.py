"""Plot data for the JSV paper figures F24 (`fig_passrule_scatter`) and F25 (`fig_passrule_bars`).

Experiment 9 (the pass rule for the sinc model). Reads only the stored Experiment-9 tables; no design is
re-evaluated and no envelope is re-fitted here, except that the envelope values written for the figures are
recomputed from the stored coefficients a, b.

Conventions (as in tests/r0_passrule_stats.py and the manuscript, Sec. 3.2):
  R      reduction against the uniform design, dB, positive = quieter
  e      = R_sinc - R_exact (dB); e > 0 means the sinc model overstates the reduction
  R_g    = R_pred - e_max, accepted at target T when R_g >= T

Inputs (results/r0_revision/passrule/):
  passrule_designs.csv    one row per evaluated design: case, set, driver (S sinc, X exact, E ERP, R random,
                          P paper design), R_sinc, R_exact, e_sinc
  passrule_envelope.csv   upper linear envelopes e <= a + b R_pred (kind 'hull'); pooled_lowkR, pooled_highkR
  passrule_guarantee.csv  the rule applied to every sinc end design, per e_max source; the figures use the
                          source 'envelope_same_class_other_cells' (envelope fitted on the other cells of the
                          same kR class, evaluated at the design's own R_pred)

Outputs (paper_jsv/prototype/figs/data/):
  f24_<case>_<sinc|exact|random>.dat   columns Rpred e, for case = paper, c30k3, c120k1p5
  f24_numbers.tex                      TeX macros: pooled envelope coefficients, R_pred needed for T = 20 dB,
                                       counts and maxima quoted in the caption
  f25_bars.dat                         one row per sinc end design of the seven regime-study cells (SU, SX):
                                       idx theta kR start Rpred Rx Rg Rxair extrap holds holdsair
  f25_numbers.tex                      TeX macros: counts quoted in the caption

Run:  python tests/r0_fig_passrule.py      (numpy only)
"""
import csv
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PR = os.path.join(ROOT, "results", "r0_revision", "passrule")
OUT = os.path.join(ROOT, "paper_jsv", "prototype", "figs", "data")
HEAD = "% Do not edit by hand: written by tests/r0_fig_passrule.py from results/r0_revision/passrule/\n"

# cell order of Table tab:pilot and Fig. fig:pilot (by kR, then theta0)
CELLS = [("c30_k1", 30, 1.0), ("c90_k1p5", 90, 1.5), ("c120_k1p5", 120, 1.5), ("c30_k3", 30, 3.0),
         ("c90_k7", 90, 7.0), ("c120_k7", 120, 7.0), ("c30_k18", 30, 18.3)]
SOURCE = "envelope_same_class_other_cells"


def read(name):
    with open(os.path.join(PR, name), newline="") as f:
        return list(csv.DictReader(f))


def write_dat(path, header, rows, fmt):
    with open(path, "w", newline="\n") as f:
        f.write(header + "\n")
        for r in rows:
            f.write(fmt % tuple(r) + "\n")


def main():
    os.makedirs(OUT, exist_ok=True)
    designs = read("passrule_designs.csv")
    env = {r["case"]: r for r in read("passrule_envelope.csv") if r["kind"] == "hull"}
    guar = [r for r in read("passrule_guarantee.csv") if r["emax_source"] == SOURCE]

    # ---------------------------------------------------------------- F24 scatter data
    groups = {"sinc": "S", "exact": "X", "random": "R"}
    stats = {}
    for case, tag in (("paper", "paper"), ("c30_k3", "c30k3"), ("c120_k1p5", "c120k1p5")):
        for gname, drv in groups.items():
            rows = [(float(d["R_sinc"]), float(d["e_sinc"])) for d in designs
                    if d["case"] == case and d["driver"] == drv]
            write_dat(os.path.join(OUT, "f24_%s_%s.dat" % (tag, gname)), "Rpred e", rows, "%.4f %.4f")
            arr = np.array(rows)
            stats[(tag, gname)] = (len(rows), arr[:, 1].min(), arr[:, 1].max())

    lo, hi = env["pooled_lowkR"], env["pooled_highkR"]
    alo, blo, ahi, bhi = (float(lo["a"]), float(lo["b"]), float(hi["a"]), float(hi["b"]))
    # sinc-driven paper-case designs above the pooled high-kR envelope (the envelope was fitted on the
    # pilot anchor cell, every 5th iterate, not on the full paper paths)
    paper_s = np.array([(float(d["R_sinc"]), float(d["e_sinc"])) for d in designs
                        if d["case"] == "paper" and d["driver"] == "S"])
    exceed = paper_s[:, 1] - (ahi + bhi * paper_s[:, 0])
    paper_p = np.array([(float(d["R_sinc"]), float(d["e_sinc"])) for d in designs
                        if d["case"] == "paper" and d["driver"] == "P"])
    exceed_p = paper_p[:, 1] - (ahi + bhi * paper_p[:, 0])
    with open(os.path.join(OUT, "f24_numbers.tex"), "w", newline="\n") as f:
        f.write(HEAD)
        f.write("\\def\\envLoA{%.2f}\\def\\envLoB{%.3f}\\def\\envHiA{%.2f}\\def\\envHiB{%.3f}\n"
                % (alo, blo, ahi, bhi))
        f.write("\\def\\needLoTtwenty{%.1f}\\def\\needHiTtwenty{%.1f}\n"
                % (float(lo["Rpred_needed_T20"]), float(hi["Rpred_needed_T20"])))
        for (tag, gname), (n, emin, emax) in stats.items():
            key = (tag + gname).replace("c30k3", "ca").replace("c120k1p5", "cb")
            f.write("\\def\\n%s{%d}\\def\\emin%s{%.2f}\\def\\emax%s{%.2f}\n" % (key, n, key, emin, key, emax))
        f.write("%% paper-case sinc-driven designs above the pooled high-kR envelope: %d of %d, by at most %.3f dB\n"
                % ((exceed > 0).sum(), len(exceed), exceed.max()))
        f.write("%% paper designs (G, B-E) above it: %d of %d, by at most %.3f dB\n"
                % ((exceed_p > 0).sum(), len(exceed_p), exceed_p.max()))
    print("paper sinc-driven above pooled high-kR envelope: %d of %d, max %.4f dB"
          % ((exceed > 0).sum(), len(exceed), exceed.max()))
    print("paper designs above it: %d of %d, max %.4f dB" % ((exceed_p > 0).sum(), len(exceed_p), exceed_p.max()))

    # ---------------------------------------------------------------- F25 bar data
    rows = []
    idx = 0
    nhold = nhold_air = nacc20 = nacc40 = ntrue20 = ntrue40 = 0
    margins = []
    for cell, th, kr in CELLS:
        for run, start in (("SU", "U"), ("SX", "X")):
            g = [r for r in guar if r["case"] == cell and r["run"] == run]
            assert len(g) == 1, (cell, run)
            g = g[0]
            idx += 1
            rpred, rx, rg = float(g["R_pred"]), float(g["R_true"]), float(g["R_guaranteed"])
            rxair = float(g["R_true_fsi_phys"])
            ext = 1 if g["note"].strip() == "extrapolated" else 0
            hold = 1 if g["guarantee_holds"] == "yes" else 0
            hair = 1 if g["guarantee_holds_fsi_phys"] == "yes" else 0
            nhold += hold
            nhold_air += hair
            nacc20 += rg >= 20
            nacc40 += rg >= 40
            ntrue20 += rx >= 20
            ntrue40 += rx >= 40
            margins.append(rx - rg)
            rows.append((idx, th, kr, start, rpred, rx, rg, rxair, ext, hold, hair))
    write_dat(os.path.join(OUT, "f25_bars.dat"),
              "idx theta kR start Rpred Rx Rg Rxair extrap holds holdsair", rows,
              "%d %d %.1f %s %.3f %.3f %.3f %.3f %d %d %d")
    with open(os.path.join(OUT, "f25_numbers.tex"), "w", newline="\n") as f:
        f.write(HEAD)
        f.write("\\def\\barsN{%d}\\def\\barsHold{%d}\\def\\barsHoldAir{%d}\n" % (len(rows), nhold, nhold_air))
        f.write("\\def\\barsAccTwenty{%d}\\def\\barsTrueTwenty{%d}\\def\\barsAccForty{%d}\\def\\barsTrueForty{%d}\n"
                % (nacc20, ntrue20, nacc40, ntrue40))
        f.write("\\def\\barsMarginMin{%.1f}\\def\\barsMarginMax{%.1f}\n" % (min(margins), max(margins)))
    print("bars: n=%d holds=%d holds_air=%d accepted20=%d true20=%d accepted40=%d true40=%d margin %.2f..%.2f"
          % (len(rows), nhold, nhold_air, nacc20, ntrue20, nacc40, ntrue40, min(margins), max(margins)))


if __name__ == "__main__":
    main()
