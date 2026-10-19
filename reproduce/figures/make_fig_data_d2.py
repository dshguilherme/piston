# REPRODUCE-PACKAGE COPY. Do not edit numerics here; edit the original and re-copy.
# Original: paper_jsv/prototype/figs/make_fig_data_d2.py  (sha256 484c378ee4aa4cb4)
# Changed: original lines 15-19 (ROOT/RES/NOV,VIZ,EXO/OUT/makedirs) -> config.py (C.RES, C.NOV, C.VIZ, C.EXO);
#          output now goes to $REPRO_OUT/figs/data (the original wrote into paper_jsv/prototype/figs/data).
#          Everything else (tables, formats, column order, file names) is unchanged.
# Reads: results/r0_revision/novelty/{cap_frozen_validation_M20,cap_flat_extrap,cap_modes_ratio,cap_frozen_sweep}.csv,
#        results/r0_revision/viz/numbers_cap_100hz.csv, results/r0_revision/exact_opt/{cap_path_powers,cap_reopt_summary}.csv
"""Extract plotting tables for the pgfplots figures F6, F8, F9 and F13 of the prototype.

Reads only existing result files under results/r0_revision/ (no solver runs) and writes small
whitespace-separated tables into prototype/figs/data/. Every output column is a value copied from,
or a stated arithmetic operation on, the named source file.

  F6  fig_verification : novelty/cap_frozen_validation_M20.csv, novelty/cap_flat_extrap.csv
  F8  fig_regime_map   : novelty/cap_modes_ratio.csv, novelty/cap_frozen_sweep.csv
  F9  fig_face_power   : viz/numbers_cap_100hz.csv
  F13 fig_opt_paths    : exact_opt/cap_path_powers.csv, exact_opt/cap_reopt_summary.csv
Run: python make_fig_data_d2.py
"""
import os, csv, math
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
import config as C

RES = C.RES
NOV, VIZ, EXO = C.NOV, C.VIZ, C.EXO
OUT = C.out("figs", "data")

DES = ["uniform", "quadvel-ones", "inputpower", "dyncomp", "quadvel-I", "quadvel-rad"]
LAB = {"uniform": "U", "quadvel-ones": "B", "inputpower": "C", "dyncomp": "D", "quadvel-I": "E", "quadvel-rad": "G"}  # design R renamed G (symbol clash with radius R)
THS = [15, 30, 60, 90, 120, 150]


def rd(name):
    with open(name, newline="") as f:
        return [r for r in csv.DictReader(f) if not str(list(r.values())[0]).startswith("#")]


def wr(stem, header, rows):
    p = os.path.join(OUT, stem + ".dat")
    with open(p, "w", newline="\n") as f:
        f.write(" ".join(header) + "\n")
        for r in rows:
            f.write(" ".join(str(v) for v in r) + "\n")
    print("wrote", p, len(rows), "rows")


# ------------------------------------------------------------------ F6 verification
v = rd(os.path.join(NOV, "cap_frozen_validation_M20.csv"))
rows = []
for i, d in enumerate(DES):
    g = lambda lev, col: float(next(r for r in v if r["level"] == lev and r["design"] == d)[col])
    rows.append([i, LAB[d], "%.6f" % g("2", "spec_over_bem_same_level"), "%.6f" % g("3", "spec_over_bem_same_level"),
                 "%.6f" % g("8", "spec_over_bem_extrap")])
# coarse = level 2 (5.2 el./wavelength), fine = level 3 (7.9), finest trace = level 8 vs extrapolated BEM
wr("f6_bem", ["idx", "design", "coarse", "fine", "finest_vs_extrap"], rows)

fe = rd(os.path.join(NOV, "cap_flat_extrap.csv"))
for r in fe:
    ka = float(r["ka"]); disc = float(r["disc_hankel"])
    pts = sorted((float(a), float(b)) for a, b in (p.split(":") for p in r["cap_values"].split(";")))
    rows = [["%.4f" % (a * a), "%.6f" % (b / disc)] for a, b in pts]
    rows.insert(0, ["0.0000", "%.6f" % (float(r["cap_extrap"]) / disc)])  # first row: Richardson extrapolation
    wr("f6_flat_ka%g" % ka, ["aR2", "cap_over_disc"], rows)

# ------------------------------------------------------------------ F8 regime map
mo = rd(os.path.join(NOV, "cap_modes_ratio.csv"))
ks = sorted({float(r["kR"]) for r in mo})
edges = [ks[0] / math.sqrt(ks[1] / ks[0])] + [math.sqrt(a * b) for a, b in zip(ks[:-1], ks[1:])] + [ks[-1] * math.sqrt(ks[-1] / ks[-2])]
edges[-1] = 40.0 * math.sqrt(40.0 / 39.0)
R = {(int(r["theta0_deg"]), r["field"], float(r["kR"])): float(r["ratio"]) for r in mo}
CLIP = 3.0
for fld in ("cos0", "cos7"):
    # one rectangle per cell: x0 x1 (kR edges, geometric midpoints), row index, log2 ratio clipped to +-3, raw ratio
    rows = []
    for iy, th in enumerate(THS):
        for ik, k in enumerate(ks):
            q = R[(th, fld, k)]
            rows.append(["%.5g" % edges[ik], "%.5g" % edges[ik + 1], iy, "%.4f" % max(-CLIP, min(CLIP, math.log2(q))), "%.5g" % q])
    wr("f8_cells_%s" % fld, ["x0", "x1", "row", "l2c", "ratio"], rows)
    # staircase boundary between cells with ratio > 1 and ratio < 1 (no interpolation); nan separates segments
    seg = []
    for iy, th in enumerate(THS):
        for ik, k in enumerate(ks):
            s = R[(th, fld, k)] > 1
            if ik + 1 < len(ks) and (R[(th, fld, ks[ik + 1])] > 1) != s:   # vertical edge
                seg += [["%.5g" % edges[ik + 1], iy - 0.5], ["%.5g" % edges[ik + 1], iy + 0.5], ["nan", "nan"]]
            if iy + 1 < len(THS) and (R[(THS[iy + 1], fld, k)] > 1) != s:  # horizontal edge
                seg += [["%.5g" % edges[ik], iy + 0.5], ["%.5g" % edges[ik + 1], iy + 0.5], ["nan", "nan"]]
    wr("f8_rho1_%s" % fld, ["x", "y"], seg)

# geometric mean over 10 <= kR <= 40 (31 values) for each opening and each analytic field
fields = ["uniform"] + ["cos%d" % j for j in range(8)]
rows = []
for iy, th in enumerate(THS):
    for ix, fld in enumerate(fields):
        vals = [R[(th, fld, k)] for k in ks if 10 <= k <= 40]
        gm = math.exp(sum(math.log(x) for x in vals) / len(vals))
        rows.append([ix, iy, "%.4f" % max(-CLIP, min(CLIP, math.log2(gm))), "%.2f" % gm, len(vals)])
wr("f8_geomean_modes", ["col", "row", "l2c", "gm", "n"], rows)

sw = rd(os.path.join(NOV, "cap_frozen_sweep.csv"))
rows = []
for ix, d in enumerate(DES):
    vals = [float(r["ratio"]) for r in sw if r["design"] == d and 10 <= float(r["kR"]) <= 40]
    gm = math.exp(sum(math.log(x) for x in vals) / len(vals))
    rows.append([ix, LAB[d], "%.4f" % math.log2(gm), "%.2f" % gm, len(vals)])
wr("f8_geomean_designs30", ["col", "design", "l2c", "gm", "n"], rows)

# ------------------------------------------------------------------ F9 face power (coarse BEM mesh, 100 Hz)
nc = rd(os.path.join(VIZ, "numbers_cap_100hz.csv"))
q = {(r["quantity"], r["design"]): float(r["value"]) for r in nc}
rows = []
for i, d in enumerate(DES):
    rows.append([i, LAB[d], "%.4f" % q[("P_plus/P_bem", d)], "%.4f" % q[("P_minus/P_bem", d)],
                 "%.4f" % q[("P_sinc_mesh/P_bem", d)],
                 "%.4f" % q[("north_fraction_bem", d)] if ("north_fraction_bem", d) in q else "nan"])
wr("f9_shares", ["idx", "design", "outer", "inner", "sinc", "farfield_outer"], rows)

# ------------------------------------------------------------------ F13 optimisation paths
pp = rd(os.path.join(EXO, "cap_path_powers.csv"))
summ = {r["run"]: r for r in rd(os.path.join(EXO, "cap_reopt_summary.csv"))}
for run in ("S0", "SR", "SX", "X0", "XR", "XS"):
    its = sorted((int(r["iter"]), float(r["dB_sinc_vs_uniform"]), float(r["dB_exact_vs_uniform"])) for r in pp
                 if r["label"] == run and int(r["iter"]) >= 1)
    wr("f13_%s" % run, ["iter", "dB_sinc", "dB_exact"], [[i, "%.4f" % a, "%.4f" % b] for i, a, b in its])
    s = summ[run]
    print("  %s niter=%s end sinc=%.3f exact=%.3f (summary %.3f %.3f)" % (run, s["niter"], its[-1][1], its[-1][2],
          float(s["dB_sinc_vs_uniform"]), float(s["dB_exact_vs_uniform"])))
