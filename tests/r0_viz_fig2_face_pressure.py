"""Fig 2: per-face pressure on the cap (BEM, m=2, 100 Hz). Reads the npz only."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.ticker
from PIL import Image

ROOT = r"C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course"
VIZ = os.path.join(ROOT, "results", "r0_revision", "viz")
Z = np.load(os.path.join(VIZ, "cap_m2_100hz_bem.npz"), allow_pickle=True)
labels = [str(s) for s in Z["labels"]]
V, T = Z["vertices"], Z["triangles"]
tri = mtri.Triangulation(V[:, 0], V[:, 1], T)
plt.rcParams.update({"font.size": 10, "text.color": "#222222", "axes.labelcolor": "#222222",
                     "xtick.color": "#222222", "ytick.color": "#222222", "axes.edgecolor": "#222222"})

designs = ["quadvel-rad", "uniform"]
fig = plt.figure(figsize=(22, 8.4), facecolor="white")
gs = fig.add_gridspec(2, 6, width_ratios=[1, 1, 1, 1, 1, 0.95], left=0.04, right=0.995, top=0.875, bottom=0.075,
                      wspace=0.62, hspace=0.36)
for r, lab in enumerate(designs):
    j = labels.index(lab)
    pp, pm, mu = np.abs(Z["p_plus"][:, j]), np.abs(Z["p_minus"][:, j]), np.abs(Z["mu"][:, j])
    Ip, Im = Z["I_plus"][:, j], Z["I_minus"][:, j]
    vmax = max(pp.max(), pm.max(), mu.max())
    lim = max(np.abs(Ip).max(), np.abs(Im).max())
    panels = [(pp, "|p+| (convex/outer face, +n side) [Pa]", "viridis", 0, vmax),
              (pm, "|p-| (concave/inner face) [Pa]", "viridis", 0, vmax),
              (mu, "|mu| = |p+ - p-| (jump) [Pa]", "viridis", 0, vmax),
              (Ip, "I+ = Re(p+ conj(vn))/2 [W/m$^2$]", "RdBu_r", -lim, lim),
              (Im, "I- = -Re(p- conj(vn))/2 [W/m$^2$]", "RdBu_r", -lim, lim)]
    for c, (f, ttl, cm, lo, hi) in enumerate(panels):
        ax = fig.add_subplot(gs[r, c]); ax.set_facecolor("white")
        h = ax.tripcolor(tri, f, shading="gouraud", cmap=cm, vmin=lo, vmax=hi)
        ax.set_aspect("equal"); ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
        ax.set_title(f"{lab}\n{ttl}", fontsize=10)
        cb = fig.colorbar(h, ax=ax, fraction=0.046, pad=0.04)
        if c < 3:
            cb.set_label("Pa")
        else:
            e = int(np.floor(np.log10(hi))); cb.ax.yaxis.offsetText.set_visible(False)
            cb.formatter = matplotlib.ticker.FuncFormatter(lambda v, _, e=e: f"{v/10**e:.1f}"); cb.update_ticks()
            cb.set_label(f"W/m$^2$ (tick values x 1e{e})")
        cb.ax.tick_params(labelsize=9)
    ax = fig.add_subplot(gs[r, 5]); ax.axis("off"); ax.set_facecolor("white")
    Pb, Pp, Pm, Ps = Z["P_bem"][j], Z["P_plus"][j], Z["P_minus"][j], Z["P_sinc_mesh"][j]
    txt = (f"{lab}\nP_bem = {Pb:.4e} W\nP_plus/P_bem  = {Pp/Pb:.4f}\nP_minus/P_bem = {Pm/Pb:.4f}\n"
           f"P_sinc/P_bem  = {Ps/Pb:.4f}\n2 P_sinc/P_bem = {2*Ps/Pb:.4f}\nP_sinc/P_plus = {Ps/Pp:.4f}\nP_sinc/P_minus= {Ps/Pm:.4f}")
    ax.text(0.0, 0.5, txt, transform=ax.transAxes, va="center", ha="left", fontsize=10, family="monospace",
            bbox=dict(boxstyle="round", facecolor="white", edgecolor="#555555"))
fig.suptitle("Spherical cap (R=10 m, half-angle 30 deg), BEM mesh m=2, f=100 Hz, kR=18.3: per-face pressure, top view (x,y).\n"
             "Caution: m=2 has only 5.24 elements per wavelength (below 6). Jump-sign s = %+d.\nI+/I- panels: linear scale, symmetric limits = max|I| over both faces of the row (not saturated)." % int(Z["jump_sign"]),
             fontsize=12, color="#222222")
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(VIZ, f"fig2_cap_face_pressure.{ext}"), facecolor="white", dpi=200)
im = np.array(Image.open(os.path.join(VIZ, "fig2_cap_face_pressure.png")).convert("RGB"))
print("corners", im[0, 0], im[0, -1], im[-1, 0], im[-1, -1])




