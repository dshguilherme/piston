# REPRODUCE-PACKAGE COPY. Do not edit numerics here; edit the original and re-copy.
# Original: tests/r0_fig_designs.py  (sha256 b7232fe472eba11c)
# Changed: original lines 43-45 (ROOT/RES/OUT) -> config.py: RES = C.RES (reads), OUT = C.out("figs","data") (= $REPRO_OUT/figs/data).
#          ROOT is no longer needed. Added "import sys", sys.path insert and "import config as C". Everything else unchanged
#          (the .tex headers written into f20_numbers.tex / f21_numbers.tex still name tests/r0_fig_designs.py, as in the original).
# Scope: in-progress figures F20 (fig_designs_paper) and F21 (fig_designs_pilot); NOT part of the audited manuscript.
# Reads: results/r0_revision/exact_opt/{cap_path_states.mat, cap_path_powers.csv, cap_reopt_summary.csv},
#        results/r0_revision/spherecap_frozen_states.mat,
#        results/r0_revision/pilot/{pilot_runs.csv, samples_th30.mat, samples_th120.mat,
#        run_c30_k3_XU.mat, run_c30_k3_SU.mat, run_c120_k1p5_XU.mat, run_c120_k1p5_SU.mat}
# Writes: $REPRO_OUT/figs/data/f20_map_{U,G,S0,X0}.dat, f20_prof.dat, f20_numbers.tex, f21_map_<cell>_<XU|SU>.dat, f21_numbers.tex
"""Thickness-map data for the JSV paper figures F20 (`fig_designs_paper`) and F21 (`fig_designs_pilot`).

The optimizer's design variable is one thickness per element (960 elements: 40 along the meridian x 24 around
the axis, element index = i_meridian + 40 * j_azimuth, apex at i = 0; see tests/obstacle_build_ctx.m, the
'spherecap' map th = alpha*u, ph = 2*pi*v, and samples_th*.mat 'ecth'/'ecph', checked below). Each element is
written as one flat-coloured polygon, without smoothing or interpolation, on a developed polar map of the cap:

    X = s cos(ph), Y = s sin(ph),   s = R * th   (meridian arc length from the apex, metres)

i.e. an azimuthal-equidistant map centred on the apex. It is one-to-one for every half-angle, including the
120 deg cap that a plan view would fold over at the equator. The polar-angle edges of an element map to circles
of radius s; each 15 deg arc of an element edge is drawn as a 3-segment polyline (sagitta < 0.1 % of s).

Inputs (all stored results, nothing recomputed):
  results/r0_revision/exact_opt/cap_path_states.mat   T (thickness of every iterate), run, iter, runs:
                                                      end designs of S0 (sinc, from uniform) and X0 (exact,
                                                      from uniform); frozen design 'quadvel-rad' = design G
                                                      (stored as run 0, iter -4), uniform = run 0, iter -1
  results/r0_revision/exact_opt/cap_path_powers.csv   dB_exact_vs_uniform / dB_sinc_vs_uniform of those rows
  results/r0_revision/exact_opt/cap_reopt_summary.csv niter of S0, X0
  results/r0_revision/pilot/run_<cell>_<XU|SU>.mat    tbest (= last iterate) of the exact (XU) and sinc (SU)
                                                      runs from the uniform design, cells c30_k3, c120_k1p5
  results/r0_revision/pilot/pilot_runs.csv            true_dB_vs_uniform (exact), pred_dB_vs_uniform (own
                                                      objective), regret_bound_dB (s), niter, stop
  results/r0_revision/pilot/samples_th<deg>.mat       element centroids, used only to check the ordering

Outputs (paper_jsv/prototype/figs/data/):
  f20_map_<U|G|S0|X0>.dat     x y h_mm, 8 rows per element (pgfplots patch, patch type=polygon, vertex count=8)
  f20_prof.dat                s_m and h_mm along meridians A (element column 0, ph = 7.5 deg) and
                              B (column 1, ph = 22.5 deg) for U, G, S0, X0; 41 rows for a const plot
  f20_numbers.tex             TeX macros: exact and sinc dB vs uniform, iterations
  f21_map_<cell>_<XU|SU>.dat  as f20_map for the two pilot cells
  f21_numbers.tex             TeX macros: dB values, s, kR, theta0, arc lengths, equator radius

Run:  python tests/r0_fig_designs.py
"""
import csv
import os
import sys

import numpy as np
import scipy.io as sio

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir))
import config as C

RES = C.RES
OUT = C.out("figs", "data")
R = 10.0
NMER, NCIRC = 40, 24
SUB = 3                         # segments per 15 deg azimuthal edge


def element_polygons(theta0_deg):
    """(960, 6, 2) polygon vertices of every element on the developed polar map, element order of the code."""
    ds = R * np.radians(theta0_deg) / NMER
    dph = 2 * np.pi / NCIRC
    polys = np.zeros((NMER * NCIRC, 2 * (SUB + 1), 2))
    for j in range(NCIRC):
        ph = j * dph + np.linspace(0, dph, SUB + 1)
        for i in range(NMER):
            s0, s1 = i * ds, (i + 1) * ds
            # inner arc forward, outer arc backward, both as SUB-segment polylines, so that neighbouring
            # rings share their vertices exactly (no slivers); at the apex the inner points coincide
            inner = [(s0 * np.cos(p), s0 * np.sin(p)) for p in ph]
            outer = [(s1 * np.cos(p), s1 * np.sin(p)) for p in ph[::-1]]
            polys[i + NMER * j] = np.array(inner + outer)
    return polys


def check_ordering(theta0_deg):
    s = sio.loadmat(os.path.join(RES, "pilot", f"samples_th{theta0_deg}.mat"), variable_names=["ecth", "ecph"])
    th = np.degrees(s["ecth"].ravel())
    ph = np.mod(np.degrees(s["ecph"].ravel()), 360)
    e = np.arange(NMER * NCIRC)
    i, j = e % NMER, e // NMER
    dth = theta0_deg / NMER
    assert np.all((th > i * dth) & (th < (i + 1) * dth)), "meridian index is not i = e mod 40"
    assert np.all((ph > j * 15) & (ph < (j + 1) * 15)), "azimuth index is not j = e div 40"


def write_map(fn, polys, h):
    with open(fn, "w") as f:
        f.write("x y h\n")
        for e in range(polys.shape[0]):
            for x, y in polys[e]:
                f.write(f"{x:.5f} {y:.5f} {1e3 * h[e]:.3f}\n")


def path_designs():
    c = sio.loadmat(os.path.join(RES, "exact_opt", "cap_path_states.mat"), squeeze_me=True,
                    variable_names=["runs", "T", "run", "iter"])
    runs, T, run, it = list(c["runs"]), c["T"], c["run"], c["iter"]

    def pick(r, i):
        idx = np.where((run == r) & (it == i))[0]
        assert idx.size == 1
        return T[:, idx[0]], int(idx[0])

    out = {}
    for key in ("S0", "X0"):
        r = runs.index(key) + 1
        last = int(np.max(it[run == r]))
        out[key] = pick(r, last)
    out["U"] = pick(0, -1)
    out["G"] = pick(0, -4)
    # identity checks against the frozen-state file (G = quadvel-rad, U = uniform)
    fz = sio.loadmat(os.path.join(RES, "spherecap_frozen_states.mat"), squeeze_me=True)
    lab = list(fz["labels"])
    assert np.allclose(out["G"][0], fz["designs"][:, lab.index("quadvel-rad")])
    assert np.allclose(out["U"][0], 0.05)
    return out


def read_csv(fn):
    with open(fn, newline="") as f:
        return list(csv.DictReader(f))


def main():
    os.makedirs(OUT, exist_ok=True)
    # ------------------------------------------------------------------ F20, paper case (30 deg, kR 18.3)
    check_ordering(30)
    polys30 = element_polygons(30)
    des = path_designs()
    pw = read_csv(os.path.join(RES, "exact_opt", "cap_path_powers.csv"))
    byidx = {int(r["idx"]): r for r in pw}
    summ = {r["run"]: r for r in read_csv(os.path.join(RES, "exact_opt", "cap_reopt_summary.csv"))}
    nums = {}
    for key, (h, idx) in des.items():
        write_map(os.path.join(OUT, f"f20_map_{key}.dat"), polys30, h)
        row = byidx[idx]
        nums[key] = (float(row["dB_exact_vs_uniform"]), float(row["dB_sinc_vs_uniform"]))
        print(f"F20 {key:3s}: idx {idx:4d} label {row['label']:20s} exact {nums[key][0]:+8.3f} dB  "
              f"sinc {nums[key][1]:+8.3f} dB  mean h {1e3 * h.mean():.2f} mm  "
              f"at bounds {np.mean((h < 0.0251) | (h > 0.0999)):.2f}")
    for key in ("S0", "X0"):   # the path file and the run summary must agree
        assert abs(nums[key][0] - float(summ[key]["dB_exact_vs_uniform"])) < 1e-3
    ds = R * np.radians(30) / NMER
    with open(os.path.join(OUT, "f20_prof.dat"), "w") as f:
        f.write("s " + " ".join(f"{k}A {k}B" for k in ("U", "G", "S0", "X0")) + "\n")
        for i in range(NMER + 1):
            ii = min(i, NMER - 1)      # last row repeats the rim element for the const plot
            vals = []
            for k in ("U", "G", "S0", "X0"):
                h = des[k][0]
                vals += [1e3 * h[ii + NMER * 0], 1e3 * h[ii + NMER * 1]]
            f.write(f"{i * ds:.5f} " + " ".join(f"{v:.3f}" for v in vals) + "\n")
    with open(os.path.join(OUT, "f20_numbers.tex"), "w", encoding="ascii") as f:
        f.write("% Generated by tests/r0_fig_designs.py from exact_opt/cap_path_powers.csv and\n"
                "% exact_opt/cap_reopt_summary.csv. Do not edit by hand.\n")
        for key, mac in (("G", "G"), ("S0", "S"), ("X0", "X")):
            f.write("\\def\\fdbx%s{%.2f}\n" % (mac, nums[key][0]))
            f.write("\\def\\fdbs%s{%.2f}\n" % (mac, nums[key][1]))
        f.write("\\def\\fitS{%d}\n\\def\\fitX{%d}\n" % (int(summ["S0"]["niter"]), int(summ["X0"]["niter"])))
    # ------------------------------------------------------------------ F21, two pilot cells
    pr = read_csv(os.path.join(RES, "pilot", "pilot_runs.csv"))
    cells = {"c30_k3": (30, "ca"), "c120_k1p5": (120, "cb")}
    with open(os.path.join(OUT, "f21_numbers.tex"), "w", encoding="ascii") as f:
        f.write("% Generated by tests/r0_fig_designs.py from pilot/pilot_runs.csv and pilot/run_*.mat.\n"
                "% Do not edit by hand. Macro prefix ca = 30 deg, kR 3; cb = 120 deg, kR 1.5.\n")
        for cell, (th0, mac) in cells.items():
            check_ordering(th0)
            polys = element_polygons(th0)
            rows = {r["run"]: r for r in pr if r["cell"] == cell}
            kR = float(rows["XU"]["kR"])
            f.write("\\def\\f%sth{%d}\n\\def\\f%skR{%g}\n" % (mac, th0, mac, kR))
            f.write("\\def\\f%sarc{%.2f}\n" % (mac, R * np.radians(th0)))
            if th0 > 90:
                f.write("\\def\\f%seq{%.2f}\n" % (mac, R * np.pi / 2))
            best = min(float(r["true_dB_vs_uniform"]) for r in rows.values())
            for run in ("XU", "SU"):
                m = sio.loadmat(os.path.join(RES, "pilot", f"run_{cell}_{run}.mat"), squeeze_me=True,
                                struct_as_record=False)
                h = m["tbest"]
                assert np.allclose(h, m["hist"].T[:, -1]), "tbest is not the end design"
                write_map(os.path.join(OUT, f"f21_map_{cell}_{run}.dat"), polys, h)
                r = rows[run]
                tag = "x" if run == "XU" else "s"
                f.write("\\def\\f%sdb%s{%.1f}\n" % (mac, tag, float(r["true_dB_vs_uniform"])))
                f.write("\\def\\f%spred%s{%.1f}\n" % (mac, tag, float(r["pred_dB_vs_uniform"])))
                f.write("\\def\\f%sgap%s{%.1f}\n" % (mac, tag, float(r["regret_bound_dB"])))
                f.write("\\def\\f%sit%s{%d}\n" % (mac, tag, int(r["niter"])))
                print(f"F21 {cell:10s} {run}: exact {float(r['true_dB_vs_uniform']):+7.2f} dB, own "
                      f"{float(r['pred_dB_vs_uniform']):+7.2f} dB, s {float(r['regret_bound_dB']):6.2f} dB, "
                      f"{r['niter']} it ({r['stop']}), h {1e3 * h.min():.1f}-{1e3 * h.max():.1f} mm, "
                      f"mean {1e3 * h.mean():.2f} mm")
            f.write("\\def\\f%sbest{%.1f}\n" % (mac, best))
            print(f"F21 {cell}: best exact design in the cell {best:+.2f} dB")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
