"""Fig 3: far-field directivity of cap, BEM vs sinc. Reads the npz only. usage: [design key list]"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

ROOT = r"C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course"
VIZ = os.path.join(ROOT, "results", "r0_revision", "viz")
Z = np.load(os.path.join(VIZ, "cap_m2_100hz_bem.npz"), allow_pickle=True)
labels = [str(s) for s in Z["labels"]]
th, phd = Z["theta_deg"], Z["phi_deg"]
plt.rcParams.update({"font.size": 10, "text.color": "#222222", "axes.labelcolor": "#222222",
                     "xtick.color": "#222222", "ytick.color": "#222222", "axes.edgecolor": "#222222"})
DR = 40.0


def make(lab, fname):
    key = lab.replace("-", "_"); j = labels.index(lab)
    B, S = Z[f"dPdO_bem__{key}"], Z[f"dPdO_sinc__{key}"]      # (181 theta, 181 phi) W/sr
    mx = B.max()
    dB = lambda a: 10 * np.log10(np.maximum(a / mx, 1e-30))
    Bd, Sd = dB(B), dB(S)
    fig = plt.figure(figsize=(17, 10.5), facecolor="white")
    gs = fig.add_gridspec(2, 3, left=0.05, right=0.97, top=0.86, bottom=0.07, wspace=0.28, hspace=0.38)

    def cut(a, ip0, ip1):
        # full 360 deg cut: alpha from +z; alpha<=180 uses phi index ip0, alpha>180 uses ip1
        al = np.arange(0, 361, 1.0)
        v = np.where(al <= 180, a[np.minimum(al, 180).astype(int), ip0], a[np.clip(360 - al, 0, 180).astype(int), ip1])
        return np.deg2rad(al), v

    def polar(spec, ip0, ip1, title):
        ax = fig.add_subplot(spec, projection="polar"); ax.set_facecolor("white")
        for arr, c, ls, nm in ((Bd, "#1b4f9c", "-", "BEM"), (Sd, "#d1495b", "--", "sinc model")):
            a, v = cut(arr, ip0, ip1)
            ax.plot(a, np.maximum(v, -DR), color=c, ls=ls, lw=2, label=nm)
        ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
        ax.set_rlim(-DR, 0); ax.set_rticks([-40, -30, -20, -10, 0]); ax.set_rlabel_position(157)
        ax.tick_params(labelsize=9)
        ax.set_title(title, fontsize=10, pad=42)
        ax.text(0, 6, "+z (convex/outer)", ha="center", va="bottom", fontsize=9, clip_on=False)
        ax.text(np.pi, 6, "-z (concave)", ha="center", va="top", fontsize=9, clip_on=False)
        ax.legend(loc="lower left", bbox_to_anchor=(-0.18, -0.12), fontsize=9, frameon=False)
        ax.text(np.deg2rad(135), -DR - 7, "dB re max BEM", fontsize=9, ha="center", va="center", clip_on=False)

    polar(gs[0, 0], 0, 90, "(a) cut phi = 0 / 180 (x-z plane; +x on the right)\n10log10(dP/d$\\Omega$ / max BEM) [dB]")
    polar(gs[0, 1], 45, 135, "(b) cut phi = 90 / 270 (y-z plane; +y on the right)\n10log10(dP/d$\\Omega$ / max BEM) [dB]")

    def amap(spec, arr, ttl, cm, lo, hi, cbl, mask=None):
        ax = fig.add_subplot(spec); ax.set_facecolor("white")
        a = np.ma.array(arr, mask=mask)
        h = ax.pcolormesh(phd, th, a, shading="auto", cmap=cm, vmin=lo, vmax=hi)
        ax.set_ylim(180, 0); ax.set_xlabel("phi [deg]"); ax.set_ylabel("theta from +z [deg]")
        ax.set_xticks(np.arange(0, 361, 60)); ax.set_yticks(np.arange(0, 181, 30)); ax.set_title(ttl, fontsize=10)
        cb = fig.colorbar(h, ax=ax, pad=0.03); cb.set_label(cbl); cb.ax.tick_params(labelsize=9)

    amap(gs[1, 0], np.maximum(Bd, -DR), "(c) BEM", "viridis", -DR, 0, "10log10(dP/d$\\Omega$ / max BEM) [dB]")
    amap(gs[1, 1], np.maximum(Sd, -DR), "(d) sinc model (same colour scale)", "viridis", -DR, 0,
         "10log10(dP/d$\\Omega$ / max BEM) [dB]")
    both_low = (Bd < -DR) & (Sd < -DR)
    diff = 10 * np.log10(np.maximum(S, 1e-300) / np.maximum(B, 1e-300))
    lim = float(np.ceil(np.abs(np.ma.array(diff, mask=both_low)).max()))
    lim = min(lim, 30.0)
    amap(gs[1, 2], np.clip(diff, -lim, lim), "(e) 10log10(sinc/BEM); masked where both < -40 dB", "RdBu_r", -lim, lim,
         "sinc minus BEM [dB]", mask=both_low)

    ax = fig.add_subplot(gs[0, 2]); ax.axis("off")
    Pb, Ps = Z["P_bem"][j], Z["P_sinc_mesh"][j]
    Ib, Is = float(Z[f"P_int_bem__{key}"]), float(Z[f"P_int_sinc__{key}"])
    nb, sb = float(Z[f"P_north_bem__{key}"]), float(Z[f"P_south_bem__{key}"])
    ns, ss = float(Z[f"P_north_sinc__{key}"]), float(Z[f"P_south_sinc__{key}"])
    txt = (f"design: {lab}  (m=2, 100 Hz, kR=18.3)\n\n"
           f"BEM:  integral of pattern = {Ib:.5e} W\n      P_bem (surface)     = {Pb:.5e} W\n      rel. err = {Ib/Pb-1:+.2e}\n\n"
           f"sinc: integral of pattern = {Is:.5e} W\n      sinc_power_mesh     = {Ps:.5e} W\n      rel. err = {Is/Ps-1:+.2e}\n\n"
           f"hemisphere split (theta<90 | theta>90):\n"
           f"  BEM : {nb:.4e} | {sb:.4e} W\n        ({nb/Ib:.3f} | {sb/Ib:.3f})\n"
           f"  sinc: {ns:.4e} | {ss:.4e} W\n        ({ns/Is:.3f} | {ss/Is:.3f})\n\n"
           f"P_sinc/P_bem = {Is/Ib:.4f}\n"
           f"max dP/d$\\Omega$ (BEM) = {mx:.3e} W/sr\nmax dP/d$\\Omega$ (sinc) = {S.max():.3e} W/sr")
    ax.text(0.0, 0.5, txt, transform=ax.transAxes, va="center", ha="left", fontsize=10, family="monospace",
            bbox=dict(boxstyle="round", facecolor="white", edgecolor="#555555"))
    fig.suptitle(f"Far-field directivity of the spherical cap, {lab}: BEM (two-sided screen) vs sinc model; BEM mesh m=2 "
                 "(5.24 elements per wavelength, below 6)", fontsize=12, color="#222222", y=0.995)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(VIZ, f"{fname}.{ext}"), facecolor="white", dpi=200)
    im = np.array(Image.open(os.path.join(VIZ, f"{fname}.png")).convert("RGB"))
    print(fname, "corners", im[0, 0], im[0, -1], im[-1, 0], im[-1, -1], "diff-map lim", lim)
    plt.close(fig)


make("quadvel-rad", "fig3_cap_directivity")
make("uniform", "fig3b_cap_directivity_uniform")


