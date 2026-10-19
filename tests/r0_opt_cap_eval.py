"""Evaluate every design on the stored optimisation paths (r0_opt_cap_reopt.m) under the sinc and exact operators.

Inputs (results/r0_revision/exact_opt/): cap_path_states.mat (designs T, structural states U, run/iter tags),
cap_exact_factor_q8.mat (exact factor Fx, Richardson L=4000/8000; Fx8000 = L=8000 only), cap_samples_q8.mat.
Sinc power: Psinc_q4 = paper's spline-consistent operator (surface 4x4, angular 32x64, the operator the sinc runs
optimised) and, independently, the spectral Funk-Hecke sinc factor (r0_opt_exact_factor.build_factors).
Outputs: cap_path_powers.csv (every design), cap_reopt_summary.csv (optima), cap_reopt_eval.log,
figures exactopt_paths.{png,pdf}, exactopt_history.{png,pdf}, exactopt_designs.{png,pdf}.
Run:  python -u tests/r0_opt_cap_eval.py
"""
import os, sys, csv
import numpy as np
import scipy.io as sio
from scipy.stats import kendalltau, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from r0_opt_exact_factor import load_samples, build_factors, OUT

RUNS = ["S0", "X0", "XS", "SX", "SR", "XR"]
DESC = {"S0": "sinc, start uniform", "X0": "exact, start uniform", "XS": "exact, start sinc optimum",
        "SX": "sinc, start exact optimum", "SR": "sinc, start paper quadvel-rad", "XR": "exact, start paper quadvel-rad"}


def db(a, b):
    return 10 * np.log10(a / b)


def main():
    logf = open(os.path.join(OUT, "cap_reopt_eval.log"), "w")

    def log(s=""):
        print(s, flush=True); logf.write(s + "\n"); logf.flush()
    d = sio.loadmat(os.path.join(OUT, "cap_path_states.mat"))
    T, U = d["T"], d["U"]
    run = d["run"].ravel().astype(int); it = d["iter"].ravel().astype(int)
    runs = [str(np.ravel(r)[0]) for r in d["runs"].ravel()]
    w0 = float(np.ravel(d["w0"])[0]); Pq4 = d["Psinc_q4"].ravel()
    flab = [str(np.ravel(v)[0]) for v in d["frozen_labels"].ravel()]
    F = sio.loadmat(os.path.join(OUT, "cap_exact_factor_q8.mat"))
    S = load_samples(8)
    Fs = build_factors(S, log=lambda s: None)["Fs"]
    nd = U.shape[1]
    Pex = np.empty(nd); Pex8 = np.empty(nd); Psp = np.empty(nd)
    for s in range(0, nd, 200):
        sl = slice(s, s + 200)
        Pex[sl] = w0 ** 2 * np.sum(np.abs(F["Fx"] @ U[:, sl]) ** 2, 0)
        Pex8[sl] = w0 ** 2 * np.sum(np.abs(F["Fx8000"] @ U[:, sl]) ** 2, 0)
        Psp[sl] = w0 ** 2 * np.sum(np.abs(Fs @ U[:, sl]) ** 2, 0)
    del Fs
    log(f"{nd} designs; sinc(q4 consistent)/sinc(spectral) over all designs: "
        f"[{(Pq4/Psp).min():.7f}, {(Pq4/Psp).max():.7f}]; exact L8000-only/Richardson: "
        f"[{(Pex8/Pex).min():.6f}, {(Pex8/Pex).max():.6f}]")
    Ps = Pq4                                       # the operator the sinc runs actually optimised
    vol = T.mean(0) / 0.10
    iu = np.where((run == 0) & (it == -1 - flab.index("uniform")))[0][0]
    ir = np.where((run == 0) & (it == -1 - flab.index("quadvel-rad")))[0][0]
    Ps_u, Px_u = Ps[iu], Pex[iu]
    log(f"uniform t=0.05: P_sinc {Ps_u:.6e} W, P_exact {Px_u:.6e} W (ratio {Ps_u/Px_u:.4f})")
    # ---------------------------------------------------------------- per-design CSV
    with open(os.path.join(OUT, "cap_path_powers.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["idx", "run", "iter", "label", "vol_frac", "P_sinc_consistent", "P_sinc_spectral", "P_exact",
                    "P_exact_L8000", "ratio_sinc_over_exact", "dB_sinc_vs_uniform", "dB_exact_vs_uniform"])
        for j in range(nd):
            lab = runs[run[j] - 1] if run[j] > 0 else "paper_" + flab[-it[j] - 1]
            w.writerow([j, run[j], it[j], lab, f"{vol[j]:.6f}", f"{Ps[j]:.9e}", f"{Psp[j]:.9e}", f"{Pex[j]:.9e}",
                        f"{Pex8[j]:.9e}", f"{Ps[j]/Pex[j]:.6f}", f"{db(Ps[j], Ps_u):.4f}", f"{db(Pex[j], Px_u):.4f}"])
    # ---------------------------------------------------------------- optima
    best = {}
    rows = []
    log("\nrun  description                      niter  P_sinc(W)     P_exact(W)    ratio   dB_sinc  dB_exact  pred-true(dB) vol")
    for r, name in enumerate(runs, start=1):
        sel = np.where((run == r) & (it == 0))[0]
        if not sel.size:
            continue
        j = sel[0]; niter = int((run == r).sum() - 1)
        best[name] = j
        row = dict(run=name, description=DESC[name], niter=niter, P_sinc=Ps[j], P_exact=Pex[j],
                   ratio=Ps[j] / Pex[j], dB_sinc_vs_uniform=db(Ps[j], Ps_u), dB_exact_vs_uniform=db(Pex[j], Px_u),
                   vol_frac=vol[j])
        row["pred_minus_true_dB"] = row["dB_sinc_vs_uniform"] - row["dB_exact_vs_uniform"]
        rows.append(row)
        log(f"{name:4s} {DESC[name]:32s} {niter:5d}  {Ps[j]:.5e}  {Pex[j]:.5e}  {Ps[j]/Pex[j]:.4f}  "
            f"{row['dB_sinc_vs_uniform']:+7.3f}  {row['dB_exact_vs_uniform']:+7.3f}  {row['pred_minus_true_dB']:+6.3f}   {vol[j]:.4f}")
    for lab in flab:
        j = np.where((run == 0) & (it == -1 - flab.index(lab)))[0][0]
        log(f"paper {lab:13s}                          {Ps[j]:.5e}  {Pex[j]:.5e}  {Ps[j]/Pex[j]:.4f}  "
            f"{db(Ps[j], Ps_u):+7.3f}  {db(Pex[j], Px_u):+7.3f}  {db(Ps[j], Ps_u)-db(Pex[j], Px_u):+6.3f}   {vol[j]:.4f}")
    # best sinc design = lowest P_sinc among sinc-run optima; best exact = lowest P_exact among exact-run optima
    sruns = [n for n in best if n[0] == "S"]; xruns = [n for n in best if n[0] == "X"]
    js = min((best[n] for n in sruns), key=lambda j: Ps[j]); jx = min((best[n] for n in xruns), key=lambda j: Pex[j])
    allj = np.arange(nd)
    j_ex_min = allj[np.argmin(Pex)]; j_s_min = allj[np.argmin(Ps)]
    log(f"\nsinc-optimal design: {runs[run[js]-1]}  exact-optimal design: {runs[run[jx]-1]}")
    log(f"lowest P_exact over ALL {nd} evaluated designs: idx {j_ex_min} ({'run '+runs[run[j_ex_min]-1] if run[j_ex_min] else 'paper'}, "
        f"iter {it[j_ex_min]}) {Pex[j_ex_min]:.6e};  lowest P_sinc over all: idx {j_s_min} {Ps[j_s_min]:.6e}")
    regret = Pex[js] / Pex[jx]
    log(f"REGRET  P_exact(sinc-opt)/P_exact(exact-opt) = {regret:.5f}  ({db(Pex[js], Pex[jx]):+.3f} dB)")
    log(f"reverse: P_sinc(exact-opt)/P_sinc(sinc-opt) = {Ps[jx]/Ps[js]:.5f} ({db(Ps[jx], Ps[js]):+.3f} dB)")
    log(f"predicted reduction of sinc-opt (sinc model): {db(Ps[js], Ps_u):+.3f} dB; true (exact) reduction of the same "
        f"design: {db(Pex[js], Px_u):+.3f} dB; best achievable (exact-opt): {db(Pex[jx], Px_u):+.3f} dB")
    log(f"paper quadvel-rad design: exact {db(Pex[ir], Px_u):+.3f} dB (regret vs exact-opt {db(Pex[ir], Pex[jx]):+.3f} dB)")
    ts, tx = T[:, js], T[:, jx]
    log(f"design difference ||t_s - t_x||/||t_x|| = {np.linalg.norm(ts-tx)/np.linalg.norm(tx):.4f}; "
        f"corr(t_s, t_x) = {np.corrcoef(ts, tx)[0,1]:.4f}; fraction of elements at bounds: sinc {np.mean((ts<0.0251)|(ts>0.0999)):.3f}, "
        f"exact {np.mean((tx<0.0251)|(tx>0.0999)):.3f}")
    # ---------------------------------------------------------------- rank agreement
    log("\nKendall tau (sinc vs exact power) over sampled designs")
    paths = run > 0
    def tau(mask, tag):
        if mask.sum() < 3:
            return
        tk = kendalltau(Ps[mask], Pex[mask]).statistic; sr = spearmanr(Ps[mask], Pex[mask]).statistic
        log(f"  {tag:58s} n={mask.sum():5d}  tau={tk:.4f}  spearman={sr:.4f}")
        return tk
    tau(paths, "all path designs (all runs)")
    for r, name in enumerate(runs, start=1):
        tau((run == r) & (it > 0), f"run {name} path")
    near = paths & (Pex <= 2 * Pex[jx])
    tau(near, "path designs within 3 dB of the exact optimum (exact)")
    near1 = paths & (Pex <= 10 ** 0.1 * Pex[jx])
    tau(near1, "path designs within 1 dB of the exact optimum (exact)")
    nears = paths & (Ps <= 2 * Ps[js])
    tau(nears, "path designs within 3 dB of the sinc optimum (sinc)")
    tau(run == 0, "the six paper designs")
    finals = np.array(sorted(best.values()) + [j for j in range(nd) if run[j] == 0])
    m = np.zeros(nd, bool); m[finals] = True
    tau(m, "run optima + six paper designs")
    # ---------------------------------------------------------------- summary CSV
    with open(os.path.join(OUT, "cap_reopt_summary.csv"), "w", newline="") as f:
        wtr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wtr.writeheader(); wtr.writerows(rows)
        f.write(f"# uniform P_sinc={Ps_u:.9e} P_exact={Px_u:.9e}\n")
        f.write(f"# sinc_opt={runs[run[js]-1]} exact_opt={runs[run[jx]-1]} regret={regret:.6f} ({db(Pex[js], Pex[jx]):+.4f} dB)\n")
    np.savez(os.path.join(OUT, "cap_path_powers.npz"), Ps=Ps, Psp=Psp, Pex=Pex, Pex8=Pex8, run=run, it=it, vol=vol,
             runs=np.array(runs), js=js, jx=jx, iu=iu, ir=ir)
    figures(T, Ps, Pex, run, it, runs, best, js, jx, iu, ir, Ps_u, Px_u, S)
    log("R0_OPT_CAP_EVAL_DONE")


def figures(T, Ps, Pex, run, it, runs, best, js, jx, iu, ir, Ps_u, Px_u, S):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.facecolor": "white", "figure.facecolor": "white",
                         "savefig.facecolor": "white", "text.color": "black", "axes.labelcolor": "black",
                         "xtick.color": "black", "ytick.color": "black"})
    col = {"S0": "#1f77b4", "SX": "#6baed6", "SR": "#08306b", "X0": "#d62728", "XS": "#fc9272", "XR": "#67000d"}
    # (1) scatter of all path designs: dB vs uniform under each operator
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    for r, name in enumerate(runs, start=1):
        m = (run == r) & (it > 0)
        if m.any():
            ax.plot(10*np.log10(Ps[m]/Ps_u), 10*np.log10(Pex[m]/Px_u), ".", ms=2.5, color=col[name], alpha=0.6,
                    label=f"{name}: {DESC[name]}")
    m0 = run == 0
    ax.plot(10*np.log10(Ps[m0]/Ps_u), 10*np.log10(Pex[m0]/Px_u), "ks", mfc="none", ms=5, label="six paper designs")
    ax.plot(10*np.log10(Ps[js]/Ps_u), 10*np.log10(Pex[js]/Px_u), "b*", ms=11, mec="k", label="sinc optimum")
    ax.plot(10*np.log10(Ps[jx]/Ps_u), 10*np.log10(Pex[jx]/Px_u), "r*", ms=11, mec="k", label="exact optimum")
    lo = min(10*np.log10(Ps[run > 0].min()/Ps_u), 10*np.log10(Pex[run > 0].min()/Px_u)) - 0.5
    ax.plot([lo, 1], [lo, 1], "k--", lw=0.8, label="predicted = true")
    ax.set_xlabel("sinc-model power change vs uniform (dB)")
    ax.set_ylabel("exact two-sided power change vs uniform (dB)")
    ax.set_title("Cap, 100 Hz: designs along all optimisation paths")
    ax.grid(alpha=0.3); ax.legend(fontsize=6, loc="upper left")
    fig.tight_layout()
    for e in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"exactopt_paths.{e}"), dpi=200)
    plt.close(fig)
    # (2) histories
    fig, axs = plt.subplots(1, 2, figsize=(8.0, 3.4), sharey=True)
    for ax, names, ttl in ((axs[0], ["S0", "SR", "SX"], "sinc-driven runs"), (axs[1], ["X0", "XR", "XS"], "exact-driven runs")):
        for name in names:
            if name not in runs:
                continue
            r = runs.index(name) + 1
            m = (run == r) & (it > 0)
            if not m.any():
                continue
            k = it[m]
            ax.plot(k, 10*np.log10(Ps[m]/Ps_u), "-", color=col[name], lw=1, label=f"{name} sinc model")
            ax.plot(k, 10*np.log10(Pex[m]/Px_u), ":", color=col[name], lw=1.4, label=f"{name} exact")
        ax.set_xscale("log"); ax.set_xlabel("outer iteration"); ax.set_title(ttl); ax.grid(alpha=0.3)
        ax.legend(fontsize=6)
    axs[0].set_ylabel("power change vs uniform (dB)")
    fig.tight_layout()
    for e in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"exactopt_history.{e}"), dpi=200)
    plt.close(fig)
    # (3) thickness maps (polar view of the cap, element centroids from the q8 samples)
    th = S["th"].reshape(-1, 64).mean(1); ph = np.arctan2(np.sin(S["ph"]).reshape(-1, 64).mean(1),
                                                            np.cos(S["ph"]).reshape(-1, 64).mean(1))
    n1, n2 = 40, 24
    TH = np.degrees(th.reshape(n2, n1).T); PH = ph.reshape(n2, n1).T
    th_c = TH.mean(1); order = np.argsort(PH[0])
    ph_c = PH[0][order]
    the = np.r_[0, 0.5*(th_c[1:]+th_c[:-1]), 30.0]
    phe = np.r_[ph_c[0]-np.pi/n2, 0.5*(ph_c[1:]+ph_c[:-1]), ph_c[-1]+np.pi/n2]
    fig, axs = plt.subplots(1, 4, figsize=(11.0, 3.4), subplot_kw=dict(projection="polar"))
    ttls = [("paper quadvel-rad (legacy sinc)", ir), ("S0: sinc, from uniform", best["S0"]),
            ("X0: exact, from uniform", jx if runs[run[jx]-1] == "X0" else best["X0"]),
            ("SX: sinc, from X0 optimum", js if runs[run[js]-1] == "SX" else best["SX"])]
    for ax, (ttl, j) in zip(axs, ttls):
        t = T[:, j].reshape(n2, n1).T[:, order]
        pc = ax.pcolormesh(phe, the, t * 1000, cmap="viridis", vmin=25, vmax=100, shading="flat")
        ax.set_title(ttl, fontsize=8, pad=10); ax.set_yticks([10, 20, 30]); ax.set_yticklabels(["10°", "20°", "30°"], fontsize=6)
        ax.set_xticklabels([])
    cb = fig.colorbar(pc, ax=axs, shrink=0.8, pad=0.04); cb.set_label("thickness (mm)")
    for e in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"exactopt_designs.{e}"), dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
