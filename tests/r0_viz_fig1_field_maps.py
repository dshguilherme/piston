"""Fig 1 / 1b: field maps of the exported normal-velocity trace (|vn|, arg vn, zoom with mesh)."""
import os, csv
import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Rectangle
from PIL import Image

R = r"C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course"
MD = os.path.join(R, "results", "r0_revision", "bem_mesh")
OUT = os.path.join(R, "results", "r0_revision", "viz")
os.makedirs(OUT, exist_ok=True)
C0 = 343.0
plt.rcParams.update({"font.size": 10, "axes.titlesize": 10, "axes.labelsize": 10,
                     "text.color": "#1a1a1a", "axes.labelcolor": "#1a1a1a",
                     "xtick.color": "#1a1a1a", "ytick.color": "#1a1a1a",
                     "axes.edgecolor": "#1a1a1a", "figure.facecolor": "white",
                     "axes.facecolor": "white", "savefig.facecolor": "white"})

CASES = [("cap", "spherecap", "quadvel-rad"), ("cap", "spherecap", "uniform"),
         ("panel", "cylpanel", "uniform"), ("panel", "cylpanel", "legacy_long_horizon")]
GEOM_NAME = {"cap": "Spherical cap (R=10 m, 30 deg)", "panel": "Cylindrical panel (R=25 m, +/-40 deg)"}


def load(stem, m):
    fn = f"{stem}_m{m}.mat"
    d = sio.loadmat(os.path.join(MD, fn))
    V = np.asarray(d["vertices"], float)
    T = np.asarray(d["triangles"]).astype(np.int64) - 1
    labels = [str(x[0]) for x in d["labels"].ravel()]
    return dict(V=V, T=T, vn=np.asarray(d["vn"]), labels=labels,
                freq=float(np.ravel(d["freq"])[0]), file=fn)


def project(geom, V):
    if geom == "cap":
        return V[:, 0], V[:, 1]
    return 25.0 * np.arctan2(V[:, 0], V[:, 1]), V[:, 2]


def stats(V, T, u):
    a, b, c = V[T[:, 0]], V[T[:, 1]], V[T[:, 2]]
    A = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    u0, u1, u2 = u[T[:, 0]], u[T[:, 1]], u[T[:, 2]]
    # exact P1 integral of |u|^2 per triangle
    i2 = A / 12.0 * (np.abs(u0) ** 2 + np.abs(u1) ** 2 + np.abs(u2) ** 2 + np.abs(u0 + u1 + u2) ** 2)
    area = A.sum()
    rms = np.sqrt(i2.sum() / area)
    e = np.stack([np.linalg.norm(b - a, axis=1), np.linalg.norm(c - b, axis=1),
                  np.linalg.norm(a - c, axis=1)])
    return area, rms, e.max()


def zoom_window(x, y, T, amp, frac):
    """Window of 1/6 extent centred on the strongest |vn| gradient (triangle-wise, interior)."""
    xa, ya = x[T], y[T]
    ar = amp[T]
    # gradient of P1 field per triangle in projected plane
    x0, x1, x2 = xa.T; y0, y1, y2 = ya.T
    det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    gx = ((ar[:, 1] - ar[:, 0]) * (y2 - y0) - (ar[:, 2] - ar[:, 0]) * (y1 - y0)) / det
    gy = (-(ar[:, 1] - ar[:, 0]) * (x2 - x0) + (ar[:, 2] - ar[:, 0]) * (x1 - x0)) / det
    g = np.hypot(gx, gy)
    # normalise by local mean to pick strong *relative* variation among strong-field regions
    m_ = ar.mean(axis=1)
    score = np.where(m_ > 0.05 * ar.max(), g / np.maximum(m_, 1e-300), 0.0)
    k = np.argmax(score)
    cx, cy = xa[k].mean(), ya[k].mean()
    w = frac * (x.max() - x.min()); h = frac * (y.max() - y.min())
    cx = np.clip(cx, x.min() + w / 2, x.max() - w / 2)
    cy = np.clip(cy, y.min() + h / 2, y.max() - h / 2)
    return (cx - w / 2, cx + w / 2, cy - h / 2, cy + h / 2)


def make_norm(amp):
    mn, mx = amp[amp > 0].min(), amp.max()
    if mx / mn > 100:
        vmin = max(mn, mx * 1e-4)
        return LogNorm(vmin=vmin, vmax=mx), True, mx / mn
    return Normalize(vmin=0, vmax=mx), False, mx / mn


def build(ms, with_phase, fname, rows_data, zooms):
    ncol = 3 if with_phase else 2
    fig, axes = plt.subplots(4, ncol, figsize=(5.6 * ncol + 1.0, 21.0 if with_phase else 21.0),
                             squeeze=False)
    fig.patch.set_facecolor("white")
    records = []
    for r, (geom, stem, design) in enumerate(CASES):
        D = load(stem, ms)
        idx = D["labels"].index(design)
        u = D["vn"][:, idx]
        amp = np.abs(u)
        x, y = project(geom, D["V"])
        tri = mtri.Triangulation(x, y, D["T"])
        area, rms, hmax = stats(D["V"], D["T"], u)
        lam = C0 / D["freq"]
        epw = lam / hmax
        norm, islog, ratio = make_norm(amp)
        if r not in zooms:
            zooms[r] = zoom_window(x, y, D["T"], amp, 1.0 / 6.0)
        z = zooms[r]
        xl, yl = ("x (m)", "y (m)") if geom == "cap" else ("arc length s (m)", "z (m)")
        ax = axes[r, 0]
        pc = ax.tripcolor(tri, amp, shading="gouraud", cmap="viridis", norm=norm)
        ax.set_aspect("equal"); ax.set_xlabel(xl); ax.set_ylabel(yl)
        ax.set_facecolor("white")
        ax.add_patch(Rectangle((z[0], z[2]), z[1] - z[0], z[3] - z[2], fill=False, ec="red", lw=1.6))
        cb = fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label(r"$|v_n|$ (m/s)" + (", log scale" if islog else ", linear"))
        ax.set_title(f"|v_n|" if r else "|v_n|", fontsize=10)
        ax.text(0.02, 0.02, f"max|v_n| = {amp.max():.3e} m/s\nRMS(v_n) = {rms:.3e} m/s",
                transform=ax.transAxes, fontsize=9, va="bottom", ha="left",
                bbox=dict(fc="white", ec="#888888", alpha=0.9, pad=2))
        col = 1
        if with_phase:
            ax2 = axes[r, 1]
            # flat shading per triangle: phase of the mean phasor of the three vertex values
            ph = np.degrees(np.angle(u[D["T"]].mean(axis=1) / (amp[D["T"]].mean(axis=1) + 1e-300)))
            # use unit phasors so tiny-amplitude vertices do not dominate
            ph = np.degrees(np.angle((u[D["T"]] / np.maximum(amp[D["T"]], 1e-300)).mean(axis=1)))
            pc2 = ax2.tripcolor(tri, facecolors=ph, cmap="twilight", vmin=-180, vmax=180)
            ax2.set_aspect("equal"); ax2.set_xlabel(xl); ax2.set_ylabel(yl)
            ax2.set_facecolor("white")
            cb2 = fig.colorbar(pc2, ax=ax2, fraction=0.046, pad=0.03, ticks=[-180, -90, 0, 90, 180])
            cb2.set_label(r"$\arg v_n$ (deg)")
            ax2.set_title("arg v_n (per-triangle, flat)", fontsize=10)
            col = 2
        ax3 = axes[r, col]
        sel = np.where(((x[D["T"]] >= z[0]) & (x[D["T"]] <= z[1]) &
                        (y[D["T"]] >= z[2]) & (y[D["T"]] <= z[3])).any(axis=1))[0]
        tri_z = mtri.Triangulation(x, y, D["T"][sel])
        pc3 = ax3.tripcolor(tri_z, amp, shading="gouraud", cmap="viridis", norm=norm)
        ax3.triplot(tri_z, color="#222222", lw=0.35)
        ax3.set_xlim(z[0], z[1]); ax3.set_ylim(z[2], z[3])
        ax3.set_aspect("equal"); ax3.set_xlabel(xl); ax3.set_ylabel(yl)
        ax3.set_facecolor("white")
        for s in ax3.spines.values():
            s.set_edgecolor("red"); s.set_linewidth(1.6)
        cb3 = fig.colorbar(pc3, ax=ax3, fraction=0.046, pad=0.03)
        cb3.set_label(r"$|v_n|$ (m/s)" + (", log scale" if islog else ", linear"))
        ax3.set_title(f"zoom (red box), {len(sel)} triangles shown, all edges drawn", fontsize=10)
        # row title
        axes[r, 0].annotate(
            f"{GEOM_NAME[geom]} | design: {design} | f = {D['freq']:.5g} Hz\nm={ms} | "
            f"{len(D['T'])} triangles | h_max = {hmax:.3f} m | lambda/h_max = {epw:.2f} "
            f"(lambda = {lam:.3f} m)",
            xy=(0.0, 1.12), xycoords="axes fraction", fontsize=10.5, fontweight="bold",
            ha="left", va="bottom", annotation_clip=False)
        for q, v, un in [("max_abs_vn", amp.max(), "m/s"), ("rms_vn_area_weighted", rms, "m/s"),
                         ("min_abs_vn", amp.min(), "m/s"),
                         ("surface_area", area, "m^2"), ("n_triangles", len(D["T"]), "1"),
                         ("h_max", hmax, "m"), ("lambda", lam, "m"),
                         ("elements_per_wavelength", epw, "1"),
                         ("frequency", D["freq"], "Hz"),
                         ("zoom_xmin", z[0], "m"), ("zoom_xmax", z[1], "m"),
                         ("zoom_ymin", z[2], "m"), ("zoom_ymax", z[3], "m"),
                         ("color_scale_log(1=yes)", int(islog), "1"),
                         ("max_over_min_abs_vn", ratio, "1")]:
            records.append([("fig1" if with_phase else "fig1b") + f"_m{ms}", q, design, geom, v, un, D["file"]])
        print(f"m={ms} {geom:5s} {design:20s} ntri={len(D['T'])} area={area:.4f} hmax={hmax:.4f} "
              f"lam={lam:.4f} epw={epw:.3f} max={amp.max():.4e} min={amp.min():.3e} rms={rms:.4e} "
              f"log={islog} ratio={ratio:.1f} zoom={np.round(z,3)}")
    fig.subplots_adjust(left=0.07, right=0.93, top=0.95, bottom=0.03, hspace=0.5, wspace=0.32)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, fname + "." + ext), facecolor="white", dpi=200)
    plt.close(fig)
    return records


zooms = {}
rec = build(2, True, "fig1_field_maps", None, zooms)
rec += build(3, False, "fig1b_field_maps_m3", None, zooms)

with open(os.path.join(OUT, "numbers_fig1.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["figure", "quantity", "design", "geometry", "value", "unit", "source_file"])
    for r_ in rec:
        w.writerow(r_)

for fn in ("fig1_field_maps.png", "fig1b_field_maps_m3.png"):
    im = Image.open(os.path.join(OUT, fn)).convert("RGB")
    W, H = im.size
    print(fn, im.size, [im.getpixel(p) for p in [(0, 0), (W - 1, 0), (0, H - 1), (W - 1, H - 1)]])



