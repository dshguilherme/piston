"""Diagnostic 4: radiated power vs kR for the frozen cap field.
Curves: sinc (MATLAB centroid + spline-consistent), sinc on the BEM mesh, BEM m=2 (FAR-FIELD power primary;
surface power shown faint because it is wrong at low kR - see novelty/bem_lowk_check.csv), the independent exact
spectral two-sided solution (novelty/cap_frozen_sweep.csv, read only), and the breathing-sphere sinc/exterior
benchmark ratio. Reads CSVs only.
Outputs: results/r0_revision/viz/fig4_cap_power_vs_frequency.{png,pdf}, numbers_fig4.csv
"""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

ROOT = r"C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course"
RES = os.path.join(ROOT, "results", "r0_revision")
VIZ = os.path.join(RES, "viz")
NOV = os.path.join(RES, "novelty")
R_CAP = 10.0
H_MAX = 0.6540312923078148          # m=2 cap max edge (bem_shell_runs/spherecap_m2_q0.json)
KR_MAX = 2 * np.pi * R_CAP / (6 * H_MAX)   # >= 6 elements per wavelength
KR_DESIGN = 2 * np.pi * 100 / 343 * R_CAP
SRC_B = "results/r0_revision/viz/cap_m2_freq_sweep.csv"
SRC_S = "results/r0_revision/viz/cap_sinc_freq_sweep.csv"
SRC_N = "results/r0_revision/novelty/cap_frozen_sweep.csv"


def rows(p):
    with open(p) as f:
        return list(csv.DictReader(f))


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def match(kq, kref, vref, tol=1e-3):
    out = np.full(len(kq), np.nan)
    for i, x in enumerate(kq):
        if len(kref):
            j = np.argmin(abs(kref - x))
            if abs(kref[j] - x) < tol:
                out[i] = vref[j]
    return out


bem = rows(os.path.join(VIZ, "cap_m2_freq_sweep.csv"))
snc = rows(os.path.join(VIZ, "cap_sinc_freq_sweep.csv"))
sph = rows(os.path.join(RES, "spherical-exterior-benchmark.csv"))
nov = rows(os.path.join(NOV, "cap_frozen_sweep.csv")) if os.path.exists(os.path.join(NOV, "cap_frozen_sweep.csv")) else []

plt.rcParams.update({"font.size": 10, "text.color": "#111111", "axes.labelcolor": "#111111",
                     "xtick.color": "#111111", "ytick.color": "#111111", "axes.edgecolor": "#333333"})
fig, axs = plt.subplots(3, 1, figsize=(13.0, 12.5), sharex=True, gridspec_kw=dict(height_ratios=[1.35, 1.0, 0.75]))
fig.patch.set_facecolor("white")
for a in axs:
    a.set_facecolor("white")
    a.axvspan(KR_MAX, 30, color="#d9d9d9", alpha=0.6, lw=0)
    a.axvline(KR_DESIGN, color="#555555", ls=":", lw=1)
    a.grid(True, which="both", color="#e6e6e6", lw=0.6)

numbers = []
styles = {"quadvel-rad": dict(c="#1f4e9c", lw=2.0, al=1.0), "uniform": dict(c="#b5651d", lw=1.2, al=0.85)}
for design, st in styles.items():
    c, lw, al = st["c"], st["lw"], st["al"]
    B = sorted([r for r in bem if r["design"] == design], key=lambda r: fnum(r["kR"]))
    kb = np.array([fnum(r["kR"]) for r in B])
    Ps = np.array([fnum(r["P_bem"]) for r in B])          # surface pairing
    Pf = np.array([fnum(r["P_far"]) for r in B])          # far field (primary)
    Pm = np.array([fnum(r["P_sinc_mesh"]) for r in B])
    Pp = np.array([fnum(r.get("P_plus")) for r in B]); Pn = np.array([fnum(r.get("P_minus")) for r in B])
    res = kb <= KR_MAX + 1e-9
    bad = np.abs(Pf / Ps - 1) > 0.02
    D = sorted([r for r in snc if r["design"] == design and int(float(r["dense_flag"])) == 1], key=lambda r: fnum(r["kR"]))
    S = sorted([r for r in snc if r["design"] == design and int(float(r["dense_flag"])) == 0
                and int(float(r["angular_order"])) == 64], key=lambda r: fnum(r["kR"]))
    kd = np.array([fnum(r["kR"]) for r in D]); kc = np.array([fnum(r["kR"]) for r in S])
    Pc = np.array([fnum(r["P_consistent"]) for r in S]); Pl = np.array([fnum(r["P_centroid"]) for r in S])
    N = sorted([r for r in nov if r["design"] == design], key=lambda r: fnum(r["kR"]))
    kn = np.array([fnum(r["kR"]) for r in N]); Pnv = np.array([fnum(r["P_2s"]) for r in N])
    Psn = np.array([fnum(r["P_sinc"]) for r in N])
    a = axs[0]
    a.loglog(kd, [fnum(r["P_consistent"]) for r in D], "-", color=c, lw=lw, alpha=al, label=f"sinc, spline-consistent [{design}]")
    a.loglog(kd, [fnum(r["P_centroid"]) for r in D], "--", color=c, lw=lw * 0.7, alpha=al, label=f"sinc, centroid [{design}]")
    a.loglog(kb, Pm, "x", color=c, ms=6, alpha=al, label=f"sinc on BEM mesh m=2 [{design}]")
    if len(kn):
        a.loglog(kn, Pnv, "-.", color="#2a2a2a" if design == "quadvel-rad" else "#8a8a8a", lw=1.3,
                 label=f"exact spectral two-sided (independent, novelty) [{design}]")
    a.loglog(kb[res], Pf[res], "o", color=c, mfc=c, ms=6.5, label=f"BEM two-sided m=2, far-field power [{design}]")
    a.loglog(kb[~res], Pf[~res], "o", color=c, mfc="white", mew=1.5, ms=6.5, label=f"BEM far-field, UNRESOLVED (<6 el./wl) [{design}]")
    a.loglog(kb, Ps, "s", color=c, mfc="none", ms=4, alpha=0.5, label=f"BEM surface power (not used) [{design}]")
    a.loglog(kb[bad], Ps[bad], "s", mfc="none", mec="red", ms=10, mew=1.0,
             label="surface power off far field by >2 % (unreliable)" if design == "quadvel-rad" else None)
    rc = match(kb, kc, Pc) / Pf
    rl = match(kb, kc, Pl) / Pf
    a = axs[1]
    a.loglog(kb, rc, "-o", color=c, lw=lw, ms=5, alpha=al, label=f"consistent sinc / BEM far field [{design}]")
    a.loglog(kb, rl, "--s", color=c, lw=lw * 0.7, ms=4, mfc="white", alpha=al, label=f"centroid sinc / BEM far field [{design}]")
    for x, y in zip(kb[~res], rc[~res]):
        a.plot(x, y, "o", mfc="white", mec="k", ms=9, mew=1.3, zorder=5)
    if len(kn):
        a.loglog(kn, Psn / Pnv, "-.", color="#2a2a2a" if design == "quadvel-rad" else "#8a8a8a", lw=1.2,
                 label=f"spectral sinc / spectral two-sided [{design}]")
    a = axs[2]
    ok = np.isfinite(Pp)
    if ok.any():
        a.semilogx(kb[ok], Pp[ok] / Ps[ok], "-^", color=c, lw=lw, alpha=al, label=f"P+ / P (outer, convex face) [{design}]")
        a.semilogx(kb[ok], Pn[ok] / Ps[ok], "-v", color=c, lw=lw * 0.7, mfc="white", alpha=al, label=f"P- / P (inner, concave face) [{design}]")
        a.semilogx(kb[ok], (Pm / Ps)[ok], ":x", color=c, lw=1, ms=7, alpha=al, label=f"P_sinc,mesh / P (surface) [{design}]")
    for i in range(len(kb)):
        base = dict(figure="fig4", design=design, geometry="spherecap", kR=kb[i], resolved=int(res[i]))
        numbers += [dict(base, quantity="P_bem_far_m2", value=Pf[i], unit="W", source_file=SRC_B),
                    dict(base, quantity="P_bem_surface_m2", value=Ps[i], unit="W", source_file=SRC_B),
                    dict(base, quantity="far_over_surface", value=Pf[i] / Ps[i], unit="-", source_file=SRC_B),
                    dict(base, quantity="P_sinc_mesh_m2", value=Pm[i], unit="W", source_file=SRC_B),
                    dict(base, quantity="ratio_consistent_over_bem_far", value=rc[i], unit="-", source_file=SRC_B + "; " + SRC_S),
                    dict(base, quantity="ratio_centroid_over_bem_far", value=rl[i], unit="-", source_file=SRC_B + "; " + SRC_S),
                    dict(base, quantity="bem_far_over_spectral", value=Pf[i] / match([kb[i]], kn, Pnv)[0], unit="-", source_file=SRC_B + "; " + SRC_N)]
        if np.isfinite(Pp[i]):
            numbers += [dict(base, quantity="P_plus_over_P_surface", value=Pp[i] / Ps[i], unit="-", source_file=SRC_B),
                        dict(base, quantity="P_minus_over_P_surface", value=Pn[i] / Ps[i], unit="-", source_file=SRC_B)]
    for x, p1, p2 in zip(kc, Pc, Pl):
        numbers += [dict(figure="fig4", design=design, geometry="spherecap", kR=x, quantity="P_consistent_sinc", value=p1, unit="W", source_file=SRC_S),
                    dict(figure="fig4", design=design, geometry="spherecap", kR=x, quantity="P_centroid_sinc", value=p2, unit="W", source_file=SRC_S)]

x = np.logspace(-1, np.log10(25), 800)
axs[1].loglog(x, 2 * (1 + x ** 2) * (np.sin(x) / x) ** 2, color="#7a7a7a", lw=0.9,
              label="breathing sphere: sinc / exact exterior (analytic)")
ks = np.array([fnum(r["ka"]) for r in sph]); rs = np.array([fnum(r["prescribed_sinc_to_exterior_power_ratio"]) for r in sph])
axs[1].loglog(ks, rs, "D", color="#7a7a7a", ms=4, label="spherical-exterior-benchmark.csv")
for kk, rr in zip(ks, rs):
    numbers.append(dict(figure="fig4", design="breathing-sphere", kR=kk, quantity="sphere_sinc_over_exterior",
                        value=rr, unit="-", source_file="results/r0_revision/spherical-exterior-benchmark.csv"))
for v in (0.5, 1.0, 2.0):
    axs[1].axhline(v, color="#222222", lw=0.8, ls="-" if v == 1 else "-.")
axs[2].axhline(0.5, color="#222222", lw=0.8, ls="-.")
axs[2].axhline(0.0, color="#222222", lw=0.6)
numbers.append(dict(figure="fig4", geometry="spherecap", kR=KR_MAX, quantity="kR_resolution_limit_6_per_lambda_m2",
                    value=KR_MAX, unit="-", source_file="bem_shell_runs/spherecap_m2_q0.json (h_max)"))

axs[0].set_ylabel("radiated power P (W)\n(normal velocity frozen at the 100 Hz field)")
axs[1].set_ylabel("power ratio sinc / two-sided (-)")
axs[1].set_ylim(5e-2, 2e3)
axs[2].set_ylabel("share of two-sided power (-)\n(symlog axis, linear for |y|<1)"); axs[2].set_yscale("symlog", linthresh=1.0); axs[2].set_ylim(-10, 30)
axs[2].set_xlabel("kR   (R = 10 m radius of curvature of the cap;  f = 343 kR / (20 pi) Hz)")
axs[2].set_xlim(0.08, 25)
axs[0].set_title(f"Frozen spherical-cap field: radiated power vs frequency (BEM mesh m=2)\n"
                 f"grey band: kR > {KR_MAX:.1f}, fewer than 6 elements per wavelength (unresolved)\n"
                 f"dotted line: design point kR = {KR_DESIGN:.2f} (100 Hz)", fontsize=9.5)
axs[0].legend(fontsize=7.5, ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0), framealpha=1, facecolor="white")
axs[1].legend(fontsize=7.5, ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0), framealpha=1, facecolor="white")
if axs[2].get_legend_handles_labels()[0]:
    axs[2].legend(fontsize=7.5, ncol=1, loc="upper left", bbox_to_anchor=(1.01, 1.0), framealpha=1, facecolor="white")
axs[1].text(0.012, 0.04, "black open circles: BEM point with < 6 elements per wavelength", transform=axs[1].transAxes, fontsize=8)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(VIZ, f"fig4_cap_power_vs_frequency.{ext}"), dpi=200, facecolor="white")
im = np.asarray(Image.open(os.path.join(VIZ, "fig4_cap_power_vs_frequency.png")).convert("RGB"))
print("corners", [tuple(int(v) for v in im[i, j]) for i in (0, -1) for j in (0, -1)])
cols = ["figure", "quantity", "design", "geometry", "kR", "resolved", "value", "unit", "source_file"]
with open(os.path.join(VIZ, "numbers_fig4.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    for n in numbers:
        w.writerow({k: n.get(k, "") for k in cols})
print("kR_max", KR_MAX, "rows", len(numbers))
