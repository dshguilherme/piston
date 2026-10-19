"""Plot data for the JSV paper figures F24 v2 (`fig_passrule_scatter_v2`) and F25 v2 (`fig_passbars_v2`),
written in the level convention of the paper (2026-10-18 convention pass).

Copy of tests/r0_fig_passrule.py (left unchanged, it still writes f24_* and f25_*). This version writes NEW
files only (f24v2_*, f25v2_*), with every level in the convention of manuscript_eng.tex, Sec. 3:

  dL     level relative to the uniform design, dB, 10 log10(P(d)/P(d0)); 0 = uniform, negative = quieter
         (dL = -R of the stored tables, which hold reductions R, positive = quieter)
  e      = dL_x - dL_pred (dB); numerically identical to the stored e = R_sinc - R_exact.
         e > 0: the design is louder than the sinc model predicts (the model overstates the quietness)
  dL_g   = dL_pred + e_max(dL_pred), guaranteed level; a design passes at target T when dL_g <= -T
  e_max  envelope e_max(dL_pred) = a - b dL_pred (same a, b as the stored e <= a + b R_pred)

No design is re-evaluated and no envelope is re-fitted; the stored values are only sign-flipped where they
are levels, and e_max of each bar is R_pred - R_guaranteed of the stored row.

Inputs (results/r0_revision/passrule/): passrule_designs.csv, passrule_envelope.csv, passrule_guarantee.csv
(source 'envelope_same_class_other_cells'), as in r0_fig_passrule.py.

Outputs (paper_jsv/prototype/figs/data/), all new:
  f24v2_<case>_<sinc|exact|random>.dat  columns dLpred e, for case = paper, c30k3, c120k1p5
  f24v2_numbers.tex                     TeX macros: envelope coefficients, dL_pred at which the rule certifies
                                        the -20 dB target (negative), counts and extremes
  f25v2_bars.dat                        one row per sinc end design of the seven regime-study cases (SU, SX):
                                        idx theta kR start dLpred dLx dLg emax extrap holds
  f25v2_numbers.tex                     TeX macros: counts quoted in the caption

Run:  python tests/r0_fig_passrule_v2.py      (numpy only)
"""
import csv
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PR = os.path.join(ROOT, "results", "r0_revision", "passrule")
OUT = os.path.join(ROOT, "paper_jsv", "prototype", "figs", "data")
HEAD = "% Do not edit by hand: written by tests/r0_fig_passrule_v2.py from results/r0_revision/passrule/\n"

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

    # ---------------------------------------------------------------- F24 v2 scatter data (levels)
    groups = {"sinc": "S", "exact": "X", "random": "R"}
    stats = {}
    for case, tag in (("paper", "paper"), ("c30_k3", "c30k3"), ("c120_k1p5", "c120k1p5")):
        for gname, drv in groups.items():
            rows = [(-float(d["R_sinc"]), float(d["e_sinc"])) for d in designs
                    if d["case"] == case and d["driver"] == drv]
            write_dat(os.path.join(OUT, "f24v2_%s_%s.dat" % (tag, gname)), "dLpred e", rows, "%.4f %.4f")
            arr = np.array(rows)
            stats[(tag, gname)] = (len(rows), arr[:, 0].min(), arr[:, 0].max(), arr[:, 1].min(), arr[:, 1].max())

    lo, hi = env["pooled_lowkR"], env["pooled_highkR"]
    alo, blo, ahi, bhi = (float(lo["a"]), float(lo["b"]), float(hi["a"]), float(hi["b"]))
    with open(os.path.join(OUT, "f24v2_numbers.tex"), "w", newline="\n") as f:
        f.write(HEAD)
        f.write("% envelope e_max(dL_pred) = a - b dL_pred, pooled over the sinc-driven iterates of each class\n")
        f.write("\\def\\envLoA{%.2f}\\def\\envLoB{%.3f}\\def\\envHiA{%.2f}\\def\\envHiB{%.3f}\n"
                % (alo, blo, ahi, bhi))
        f.write("% dL_pred at or below which the rule certifies the -20 and -40 dB targets (negative levels)\n")
        f.write("\\def\\needLoTtwenty{%.1f}\\def\\needHiTtwenty{%.1f}\\def\\needLoTforty{%.1f}\n"
                % (-float(lo["Rpred_needed_T20"]), -float(hi["Rpred_needed_T20"]), -float(lo["Rpred_needed_T40"])))
        for (tag, gname), (n, xmin, xmax, emin, emax) in stats.items():
            key = (tag + gname).replace("c30k3", "ca").replace("c120k1p5", "cb")
            f.write("\\def\\n%s{%d}\\def\\emin%s{%.2f}\\def\\emax%s{%.2f}\\def\\xmin%s{%.2f}\\def\\xmax%s{%.2f}\n"
                    % (key, n, key, emin, key, emax, key, xmin, key, xmax))
    print("envelopes: low a=%.5f b=%.6f; high a=%.5f b=%.6f" % (alo, blo, ahi, bhi))
    print("certify -20 dB at dL_pred <= %.4f (kR<=3), %.4f (kR>=7); -40 dB at %.4f (kR<=3)"
          % (-float(lo["Rpred_needed_T20"]), -float(hi["Rpred_needed_T20"]), -float(lo["Rpred_needed_T40"])))

    # ---------------------------------------------------------------- F25 v2 bar data (levels)
    rows = []
    idx = 0
    nhold = nacc20 = nacc40 = ntrue20 = ntrue40 = 0
    margins = []
    for cell, th, kr in CELLS:
        for run, start in (("SU", "U"), ("SX", "X")):
            g = [r for r in guar if r["case"] == cell and r["run"] == run]
            assert len(g) == 1, (cell, run)
            g = g[0]
            idx += 1
            dlpred, dlx, dlg = -float(g["R_pred"]), -float(g["R_true"]), -float(g["R_guaranteed"])
            emax = dlg - dlpred
            assert abs(emax - float(g["e_max"])) < 1e-3, (cell, run, emax, g["e_max"])
            ext = 1 if g["note"].strip() == "extrapolated" else 0
            hold = 1 if g["guarantee_holds"] == "yes" else 0
            assert hold == (dlg >= dlx), (cell, run)
            nhold += hold
            nacc20 += dlg <= -20
            nacc40 += dlg <= -40
            ntrue20 += dlx <= -20
            ntrue40 += dlx <= -40
            margins.append(dlg - dlx)
            rows.append((idx, th, kr, start, dlpred, dlx, dlg, emax, ext, hold))
    write_dat(os.path.join(OUT, "f25v2_bars.dat"),
              "idx theta kR start dLpred dLx dLg emax extrap holds", rows,
              "%d %d %.1f %s %.3f %.3f %.3f %.3f %d %d")
    with open(os.path.join(OUT, "f25v2_numbers.tex"), "w", newline="\n") as f:
        f.write(HEAD)
        f.write("\\def\\barsN{%d}\\def\\barsHold{%d}\n" % (len(rows), nhold))
        f.write("\\def\\barsAccTwenty{%d}\\def\\barsTrueTwenty{%d}\\def\\barsAccForty{%d}\\def\\barsTrueForty{%d}\n"
                % (nacc20, ntrue20, nacc40, ntrue40))
        f.write("\\def\\barsMarginMin{%.1f}\\def\\barsMarginMax{%.1f}\n" % (min(margins), max(margins)))
    print("bars: n=%d holds=%d accepted(-20)=%d true(-20)=%d accepted(-40)=%d true(-40)=%d margin dLg-dLx %.2f..%.2f"
          % (len(rows), nhold, nacc20, ntrue20, nacc40, ntrue40, min(margins), max(margins)))
    for r in rows:
        print("  %2d %3d %5.1f %s pred %8.3f true %8.3f guar %8.3f emax %6.3f ext %d" % r[:9])


if __name__ == "__main__":
    main()
