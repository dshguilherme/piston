"""Experiment 9: quantify the 'enough' pass rule for the sinc (piston) model with data already on disk.

Rule (interpretation, sign convention to confirm with the author):
    R = -(change of radiated power level vs the uniform design), dB, positive = reduction.
    e = R_sinc_pred - R_exact  (dB), positive = the sinc model OVERSTATES the reduction.
    R_guaranteed = R_sinc_pred(design) - e_max  >=  T  (T = 10, 20, 30, 40 dB).

Inputs (read only):
    results/r0_revision/pilot/eval_<cell>.csv          every evaluated design, all rungs
    results/r0_revision/exact_opt/cap_path_powers.csv  paper case (anchor), every iterate of 6 runs + 6 paper designs
    results/r0_revision/longruns/hist_*.csv            c30_k3 long runs (every iterate, sinc and exact dB)
    results/r0_revision/fsi/fsi_summary.csv            fluid-loaded true reductions of the run end designs
Outputs (results/r0_revision/passrule/):
    passrule_designs.csv      one row per design: R_sinc, R_erp, R_exact, e_sinc, e_erp
    passrule_population.csv   error statistics per case x population x model
    passrule_envelope.csv     upper linear envelopes e <= a + b R_pred and binned maxima
    passrule_guarantee.csv    the rule applied to every sinc end design, per e_max source, per T
    passrule_expected.csv     expected-value version over the sinc runs of each case
Only numpy, scipy (linprog) and the csv module are used (pandas is not installed in the base env).
"""
import os, csv
import numpy as np
from scipy.optimize import linprog

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results", "r0_revision")
OUT = os.path.join(RES, "passrule")
os.makedirs(OUT, exist_ok=True)
CELLS = ["c30_k1", "c30_k3", "c90_k1p5", "c120_k1p5", "c90_k7", "c120_k7", "c30_k18"]
TS = [10, 20, 30, 40]


def read(fn):
    with open(fn, newline="") as f:
        return list(csv.DictReader(f))


def write(fn, rows, cols=None):
    cols = cols or list(rows[0].keys())
    with open(os.path.join(OUT, fn), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})


def red(P, P0):  # reduction in dB, positive = quieter than uniform
    return -10.0 * np.log10(P / P0)


# ---------------------------------------------------------------- designs
D = []  # dicts: case, src, pop_tag, run, family, iter, R_sinc, R_erp, R_exact, R_exact_L8000, is_end
for cell in CELLS:
    rows = read(os.path.join(RES, "pilot", f"eval_{cell}.csv"))
    u = [r for r in rows if r["set"] == "uniform"][0]
    P0 = {k: float(u[k]) for k in ("P_sinc", "P_erp", "P_exact", "P_exact_L8000")}
    for r in rows:
        if r["set"] == "uniform":
            continue
        run = r["run"]
        D.append(dict(case=cell, src="pilot_eval", set=r["set"], run=run, family=run,
                      driver=run[0] if run else "R", iter=int(r["iter"]),
                      R_sinc=red(float(r["P_sinc"]), P0["P_sinc"]), R_erp=red(float(r["P_erp"]), P0["P_erp"]),
                      R_exact=red(float(r["P_exact"]), P0["P_exact"]),
                      R_exact_L8000=red(float(r["P_exact_L8000"]), P0["P_exact_L8000"]),
                      is_end=r["set"] == "optimum"))

# paper case: every iterate of S0, X0, XS, SX, SR, XR plus the six paper designs; ERP not available here
rows = read(os.path.join(RES, "exact_opt", "cap_path_powers.csv"))
u = [r for r in rows if r["label"] == "paper_uniform"][0]
P0s, P0x, P0x8 = float(u["P_sinc_consistent"]), float(u["P_exact"]), float(u["P_exact_L8000"])
paper = []
for r in rows:
    if r["label"] == "paper_uniform":
        continue
    lab = r["label"]
    is_run = r["run"] != "0"
    paper.append(dict(case="paper", src="cap_path_powers", set="path" if is_run else "paper_design",
                      run=lab if is_run else lab, family=lab, driver=lab[0] if is_run else "P", iter=int(r["iter"]),
                      R_sinc=red(float(r["P_sinc_consistent"]), P0s), R_erp=np.nan,
                      R_exact=red(float(r["P_exact"]), P0x), R_exact_L8000=red(float(r["P_exact_L8000"]), P0x8),
                      is_end=False))
# end design of each run = best own objective along the run (matches cap_reopt_summary.csv)
for lab in ["S0", "X0", "XS", "SX", "SR", "XR"]:
    pr = [d for d in paper if d["run"] == lab]
    key = "R_sinc" if lab[0] == "S" else "R_exact"
    best = max(pr, key=lambda d: d[key])
    best["is_end"] = True
    best["set"] = "optimum"
D += paper
# the paper case has no random fields; the anchor pilot cell is the same cap, so its 40 random fields are borrowed
for d in [d for d in D if d["case"] == "c30_k18" and d["set"] == "random"]:
    e = dict(d); e["case"] = "paper"; e["src"] = "pilot_eval(c30_k18)"; D.append(e)

# c30_k3 long runs: pilot c30_k3 designs + every iterate of the long histories (sinc stages at the target kR only)
long = []
for d in [d for d in D if d["case"] == "c30_k3"]:
    e = dict(d); e["case"] = "c30_k3_long"; long.append(e)
LDIR = os.path.join(RES, "longruns")
fam = {"A_SU": "SU", "A_XU": "XU", "A_XS": "XS", "B_P": "B_P", "B_S1": "B_S1", "B_S2": "B_S2", "B_S3": "B_S3"}
for name in ["A_c30_k3_SU", "A_c30_k3_XU", "A_c30_k3_XS", "B_c30_k3_P", "B_c30_k3_S1", "B_c30_k3_S2", "B_c30_k3_S3"]:
    h = read(os.path.join(LDIR, f"hist_{name}.csv"))
    short = name.replace("_c30_k3", "")
    for stage in sorted(set(r["stage"] for r in h)):
        hs = [r for r in h if r["stage"] == stage]
        rung = hs[0]["rung"]
        if rung == "E":
            continue
        # B_S2 and B_S3 sinc stages (sinc from uniform, 300 its) reproduce the pilot SU path to all printed digits
        # (checked: R_sinc, R_exact at its 100/200/300 identical; the S3 move limit never binds), so they are
        # duplicates of SU and are not added; their exact stages (started at SU@300) are kept.
        if rung == "S" and short in ("B_S2", "B_S3"):
            continue
        run = f"{short}:{rung}{stage}"
        recs = [dict(case="c30_k3_long", src=f"hist_{name}", set="path", run=run, family=fam[short], driver=rung,
                     iter=int(r["it"]), R_sinc=-float(r["sinc_dB"]), R_erp=-float(r["erp_dB"]),
                     R_exact=-float(r["exact_dB"]), R_exact_L8000=np.nan, is_end=False) for r in hs]
        key = "R_sinc" if rung == "S" else "R_exact"
        best = max(recs, key=lambda d: d[key]); best["is_end"] = True; best["set"] = "optimum"
        long += recs
D += long
CASES = CELLS + ["paper", "c30_k3_long"]

for d in D:
    d["e_sinc"] = d["R_sinc"] - d["R_exact"]
    d["e_erp"] = d["R_erp"] - d["R_exact"]
    d["e_sinc_L8000"] = d["R_sinc"] - d["R_exact_L8000"] if np.isfinite(d["R_exact_L8000"]) else np.nan
write("passrule_designs.csv", D)


# ---------------------------------------------------------------- populations
def pops(case):
    A = [d for d in D if d["case"] == case]
    out = {
        "all": A,
        "random": [d for d in A if d["set"] == "random"],
        "sinc_paths": [d for d in A if d["driver"] == "S"],
        "exact_paths": [d for d in A if d["driver"] == "X"],
        "erp_paths": [d for d in A if d["driver"] == "E"],
        "end_designs": [d for d in A if d["is_end"]],
        "sinc_ends": [d for d in A if d["is_end"] and d["driver"] == "S"],
    }
    if case == "paper":
        out["paper_designs"] = [d for d in A if d["set"] == "paper_design"]
    return out


def stats(x):
    x = np.asarray([v for v in x if np.isfinite(v)])
    if len(x) == 0:
        return dict(n=0, mean=np.nan, std=np.nan, p05=np.nan, p50=np.nan, p95=np.nan, max=np.nan, min=np.nan)
    return dict(n=len(x), mean=float(x.mean()), std=float(x.std(ddof=1)) if len(x) > 1 else np.nan,
                p05=float(np.percentile(x, 5)), p50=float(np.median(x)), p95=float(np.percentile(x, 95)),
                max=float(x.max()), min=float(x.min()))


POP = []
for case in CASES:
    for pname, P in pops(case).items():
        for model, col in (("sinc", "e_sinc"), ("erp", "e_erp"), ("sinc_L8000", "e_sinc_L8000")):
            s = stats([d[col] for d in P])
            if s["n"] == 0:
                continue
            rp = stats([d["R_sinc" if model != "erp" else "R_erp"] for d in P])
            POP.append(dict(case=case, population=pname, model=model, **s,
                            Rpred_min=rp["min"], Rpred_max=rp["max"], small=("yes" if s["n"] < 10 else "")))
write("passrule_population.csv", POP)


# ---------------------------------------------------------------- envelopes e <= a + b R_pred
def upper_line(x, y):
    """Lowest line above all points (min sum of gaps), slope constrained >= 0 (the overstatement is not
    allowed to shrink with depth; without this a narrow training range extrapolates to nonsense). Returns a, b."""
    x = np.asarray(x); y = np.asarray(y)
    if len(x) < 2:
        return np.nan, np.nan
    c = np.array([len(x), x.sum()])
    res = linprog(c, A_ub=-np.column_stack([np.ones_like(x), x]), b_ub=-y, bounds=[(None, None), (0, None)],
                  method="highs")
    return (res.x[0], res.x[1]) if res.success else (np.nan, np.nan)


def qline(x, y, q=0.95):
    """Linear quantile regression (pinball loss) via LP."""
    x = np.asarray(x); y = np.asarray(y); n = len(x)
    if n < 5:
        return np.nan, np.nan
    # variables a, b, u+ (n), u- (n): y - a - b x = u+ - u-
    c = np.concatenate([[0, 0], q * np.ones(n), (1 - q) * np.ones(n)])
    Aeq = np.hstack([np.column_stack([np.ones(n), x]), np.eye(n), -np.eye(n)])
    res = linprog(c, A_eq=Aeq, b_eq=y, bounds=[(None, None)] * 2 + [(0, None)] * (2 * n), method="highs")
    return (res.x[0], res.x[1]) if res.success else (np.nan, np.nan)


ENV = []
envfit = {}
BINS = np.arange(-10, 95, 5.0)
LOWKR = {"c30_k1", "c30_k3", "c90_k1p5", "c120_k1p5"}  # kR <= 3 (known before any optimisation)
for case in CASES + ["pooled_pilot", "pooled_lowkR", "pooled_highkR"]:
    if case == "pooled_pilot":
        S = [d for d in D if d["case"] in CELLS and d["driver"] == "S"]
    elif case == "pooled_lowkR":
        S = [d for d in D if d["case"] in LOWKR and d["driver"] == "S"]
    elif case == "pooled_highkR":
        S = [d for d in D if d["case"] in CELLS and d["case"] not in LOWKR and d["driver"] == "S"]
    else:
        S = [d for d in D if d["case"] == case and d["driver"] == "S"]
    x = np.array([d["R_sinc"] for d in S]); y = np.array([d["e_sinc"] for d in S])
    a, b = upper_line(x, y); qa, qb = qline(x, y)
    envfit[case] = (a, b)
    # how tight is a single e_max compared with the envelope on the same points
    emax = y.max()
    ENV.append(dict(case=case, kind="hull", n=len(x), a=a, b=b, qa95=qa, qb95=qb, e_max=emax,
                    mean_slack_emax=float(np.mean(emax - y)), mean_slack_env=float(np.mean(a + b * x - y)),
                    corr_e_Rpred=float(np.corrcoef(x, y)[0, 1]),
                    # predicted reduction needed so that R_pred - (a + b R_pred) >= T
                    **{f"Rpred_needed_T{T}": (float((T + a) / (1 - b)) if b < 1 else np.inf) for T in TS},
                    **{f"Rpred_needed_emax_T{T}": float(T + emax) for T in TS}))
    for lo in BINS:
        m = (x >= lo) & (x < lo + 5)
        if m.sum() == 0:
            continue
        ENV.append(dict(case=case, kind="bin", n=int(m.sum()), a=lo, b=lo + 5, e_max=float(y[m].max()),
                        qa95=float(np.percentile(y[m], 95)), qb95=float(y[m].mean())))
write("passrule_envelope.csv", ENV,
      ["case", "kind", "n", "a", "b", "qa95", "qb95", "e_max", "mean_slack_emax", "mean_slack_env", "corr_e_Rpred"]
      + [f"Rpred_needed_T{T}" for T in TS] + [f"Rpred_needed_emax_T{T}" for T in TS])

# ---------------------------------------------------------------- FSI true reductions of the end designs
fsi = {}
for r in read(os.path.join(RES, "fsi", "fsi_summary.csv")):
    fsi[(r["cell"], r["design"])] = (-float(r["dB_vs_uniform_dev"]), -float(r["dB_vs_uniform_phys"]))
for k in list(fsi):
    if k[0] == "c30_k18":  # pilot anchor SU/SX/XU/XS are the exp5 S0/SX/X0/XS
        m = {"SU": "S0", "SX": "SX", "XU": "X0", "XS": "XS"}.get(k[1])
        if m:
            fsi[("paper", m)] = fsi[k]
    if k[0] == "c30_k3":
        fsi[("c30_k3_long", k[1])] = fsi[k]


# ---------------------------------------------------------------- the rule applied to every sinc end design
def emax_sources(case, end):
    """e_max candidates for a given sinc end design. Returns {name: (e_bound, note)}."""
    P = pops(case)
    out = {}
    out["random_max"] = max(d["e_sinc"] for d in P["random"])
    out["random_p95"] = float(np.percentile([d["e_sinc"] for d in P["random"]], 95))
    out["exact_paths_max"] = max(d["e_sinc"] for d in P["exact_paths"])
    out["all_max(circular)"] = max(d["e_sinc"] for d in P["all"])
    out["sinc_paths_max(circular)"] = max(d["e_sinc"] for d in P["sinc_paths"])
    def env(pts):
        xs = [d["R_sinc"] for d in pts]
        a, b = upper_line(xs, [d["e_sinc"] for d in pts])
        extrap = "extrapolated" if not (min(xs) <= end["R_sinc"] <= max(xs)) else ""
        return a + b * end["R_sinc"], extrap
    loo = [d for d in P["sinc_paths"] if d["family"] != end["family"]]
    if loo:
        out["sinc_paths_LOO_max"] = max(d["e_sinc"] for d in loo)
        out["sinc_paths_LOO_p95"] = float(np.percentile([d["e_sinc"] for d in loo], 95))
        out["envelope_LOO"] = env(loo)
    base = case if case in CELLS else ("c30_k18" if case == "paper" else "c30_k3")
    other = [d for d in D if d["case"] in CELLS and d["case"] != base and d["driver"] == "S"]
    out["other_cells_sinc_max"] = max(d["e_sinc"] for d in other)
    out["envelope_other_cells"] = env(other)
    same = [d for d in other if (d["case"] in LOWKR) == (base in LOWKR)]
    out["same_class_other_cells_max"] = max(d["e_sinc"] for d in same)
    out["envelope_same_class_other_cells"] = env(same)
    return out


G = []
for case in CASES:
    ends = [d for d in D if d["case"] == case and d["is_end"] and d["driver"] == "S"]
    for end in ends:
        fd = fsi.get((case, end["run"]), (np.nan, np.nan))
        for src, em in emax_sources(case, end).items():
            note = ""
            if isinstance(em, tuple):
                em, note = em
            Rg = end["R_sinc"] - em
            row = dict(case=case, run=end["run"], iter=end["iter"], R_pred=end["R_sinc"], R_true=end["R_exact"],
                       e_self=end["e_sinc"], R_true_fsi_dev=fd[0], R_true_fsi_phys=fd[1], emax_source=src,
                       note=note, e_max=em, R_guaranteed=Rg, guarantee_holds=("yes" if end["R_exact"] >= Rg - 1e-9 else "no"),
                       guarantee_holds_fsi_phys=("" if not np.isfinite(fd[1]) else ("yes" if fd[1] >= Rg else "no")))
            for T in TS:
                row[f"rule_pass_T{T}"] = "pass" if Rg >= T else "fail"
                row[f"true_pass_T{T}"] = "pass" if end["R_exact"] >= T else "fail"
            G.append(row)
write("passrule_guarantee.csv", G)

# ---------------------------------------------------------------- expected-value version over the sinc runs of a case
EV = []
for case in CASES:
    ends = [d for d in D if d["case"] == case and d["is_end"] and d["driver"] == "S"]
    byfam = {}
    for d in ends:  # one end design per independent run family (the longest continuation)
        if d["family"] not in byfam or d["iter"] > byfam[d["family"]]["iter"]:
            byfam[d["family"]] = d
    ends = list(byfam.values())
    if len(ends) < 2:
        continue
    Rp = np.array([d["R_sinc"] for d in ends]); Rt = np.array([d["R_exact"] for d in ends])
    P = pops(case)
    for src, em in (("random_max", max(d["e_sinc"] for d in P["random"])),
                    ("exact_paths_max", max(d["e_sinc"] for d in P["exact_paths"])),
                    ("sinc_paths_max(circular)", max(d["e_sinc"] for d in P["sinc_paths"])),
                    ("sinc_ends_max(circular)", float(max(Rp - Rt)))):
        Rg = Rp.mean() - em
        EV.append(dict(case=case, n_sinc_runs=len(ends), runs="|".join(d["run"] for d in ends),
                       E_R_pred=float(Rp.mean()), E_R_true=float(Rt.mean()), min_R_true=float(Rt.min()),
                       emax_source=src, e_max=float(em), R_guaranteed=float(Rg),
                       holds_for_mean=("yes" if Rt.mean() >= Rg else "no"),
                       holds_for_every_run=("yes" if Rt.min() >= Rg else "no"),
                       **{f"rule_pass_T{T}": ("pass" if Rg >= T else "fail") for T in TS}))
write("passrule_expected.csv", EV)

# ---------------------------------------------------------------- console summary
print("case population model n mean std p95 max")
for r in POP:
    if r["model"] in ("sinc", "erp"):
        print(f"{r['case']:12s} {r['population']:14s} {r['model']:5s} {r['n']:5d} {r['mean']:8.2f} {r['std']:7.2f} "
              f"{r['p95']:8.2f} {r['max']:8.2f}  Rpred [{r['Rpred_min']:.1f},{r['Rpred_max']:.1f}]")
print()
for r in ENV:
    if r["kind"] == "hull":
        print(f"env {r['case']:12s} n={r['n']:5d} e<= {r['a']:.2f} + {r['b']:.3f} R  q95: {r['qa95']:.2f}+{r['qb95']:.3f}R "
              f"emax={r['e_max']:.2f} slack emax {r['mean_slack_emax']:.2f} env {r['mean_slack_env']:.2f} "
              f"corr {r['corr_e_Rpred']:.2f} need(env) T20 {r['Rpred_needed_T20']:.1f} T40 {r['Rpred_needed_T40']:.1f}")

# ---------------------------------------------------------------- per e_max source: how often the guarantee holds,
# and the rule's decisions against the truth, over the 14 sinc end designs of the 7 pilot cells (SU, SX)
SS = []
srcs = list(dict.fromkeys(r["emax_source"] for r in G))
for src in srcs:
    for grp, cs in (("all7", CELLS), ("lowkR", sorted(LOWKR)), ("highkR", ["c90_k7", "c120_k7", "c30_k18"])):
        R = [r for r in G if r["emax_source"] == src and r["case"] in cs]
        row = dict(emax_source=src, group=grp, n=len(R), holds=sum(r["guarantee_holds"] == "yes" for r in R),
                   holds_fsi_phys=sum(r["guarantee_holds_fsi_phys"] == "yes" for r in R),
                   mean_margin=float(np.mean([r["R_true"] - r["R_guaranteed"] for r in R])),
                   min_margin=float(np.min([r["R_true"] - r["R_guaranteed"] for r in R])))
        for T in TS:
            row[f"false_accept_T{T}"] = sum(r[f"rule_pass_T{T}"] == "pass" and r[f"true_pass_T{T}"] == "fail" for r in R)
            row[f"false_reject_T{T}"] = sum(r[f"rule_pass_T{T}"] == "fail" and r[f"true_pass_T{T}"] == "pass" for r in R)
        SS.append(row)
write("passrule_source_summary.csv", SS)
print()
for r in SS:
    print(" ".join(f"{k}={v:.2f}" if isinstance(v, float) else f"{k}={v}" for k, v in r.items()))
print()
for r in EV:
    print(f"{r['case']:12s} n={r['n_sinc_runs']} {r['runs']:28s} E[Rp]={r['E_R_pred']:.2f} E[Rt]={r['E_R_true']:.2f} "
          f"minRt={r['min_R_true']:.2f} {r['emax_source']:26s} emax={r['e_max']:.2f} Rg={r['R_guaranteed']:.2f} "
          f"mean:{r['holds_for_mean']} every:{r['holds_for_every_run']}")
