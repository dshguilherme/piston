"""Frozen-shell-field open-screen BEM power (bempp-cl) with mesh/quadrature convergence.

usage:
  python r0_bem_shells.py run <spherecap|cylpanel> <m> [qbump]   -> json in results/r0_revision/bem_shell_runs
  python r0_bem_shells.py agg                                    -> CSVs + npz + printed tables
Imports validated solver pieces from r0_bem_core (not modified).
"""
import sys, os, json, time, gc
import numpy as np
import scipy.io as sio
import scipy.linalg as sla

ROOT = r"C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course"
sys.path.insert(0, os.path.join(ROOT, "src"))
import bempp_cl.api as bem
import r0_bem_core as core

MESH = os.path.join(ROOT, "results", "r0_revision", "bem_mesh")
RUNS = os.path.join(ROOT, "results", "r0_revision", "bem_shell_runs")
OUT = os.path.join(ROOT, "results", "r0_revision")


def load(geom, m):
    d = sio.loadmat(os.path.join(MESH, f"{geom}_m{m}.mat"))
    labels = [str(x[0]) for x in d["labels"].ravel()]
    return dict(V=d["vertices"].astype(float), T=d["triangles"].astype(np.int64) - 1,
                vn=np.asarray(d["vn"], dtype=complex), N=d["normals"].astype(float),
                labels=labels, freq=float(d["freq"].ravel()[0]), omega=float(d["omega"].ravel()[0]),
                rho0=float(d["rho0"].ravel()[0]), c0=float(d["c0"].ravel()[0]))


def run(geom, m, qbump=0):
    D = load(geom, m)
    q = bem.GLOBAL_PARAMETERS.quadrature
    print(f"[{geom} m={m}] default quadrature regular={q.regular} singular={q.singular}; bump=+{qbump}", flush=True)
    q.regular += qbump
    q.singular += qbump
    print(f"   using regular={q.regular} singular={q.singular}", flush=True)
    omega, rho0, c0 = D["omega"], D["rho0"], D["c0"]
    k = omega / c0
    grid = bem.Grid(np.ascontiguousarray(D["V"].T), np.ascontiguousarray(D["T"].T.astype(np.uint32)))
    nv, nt = grid.number_of_vertices, grid.number_of_elements
    # normal agreement
    en = np.asarray(grid.normals)                       # (nt,3)
    vn_mean = D["N"][D["T"]].mean(axis=1)
    dots = np.einsum("ij,ij->i", en, vn_mean)
    bad = int(np.sum(dots <= 0))
    print(f"   nv={nv} nt={nt} k={k:.5f} normal disagreements = {bad} (min dot {dots.min():.4f})", flush=True)
    if bad:
        print("   STOP: normals disagree"); sys.exit(2)
    edges = np.concatenate([D["V"][D["T"][:, [0, 1]]], D["V"][D["T"][:, [1, 2]]], D["V"][D["T"][:, [2, 0]]]])
    hmax = float(np.linalg.norm(edges[:, 0] - edges[:, 1], axis=1).max())
    bnd = core.boundary_vertices(grid)
    ii = np.where(~bnd)[0]
    n = len(ii)
    print(f"   h_max={hmax:.4f}  interior dofs={n}  boundary vertices={int(bnd.sum())}  est dense matrix GB={n*n*16/1e9:.2f}", flush=True)
    t0 = time.time()
    space, vmap = core.full_p1(grid)
    M = core.mass_matrix_vertex(grid)
    op = bem.operators.boundary.helmholtz.hypersingular(space, space, space, k)
    A = np.asarray(op.weak_form().A)
    di = vmap[ii]
    Aii = A[np.ix_(di, di)]
    del A, op
    gc.collect()
    Aii *= core.SIGN_SCREEN
    t_asm = time.time() - t0
    print(f"   assembled+restricted in {t_asm:.1f}s", flush=True)
    vn = D["vn"]
    rhs = 1j * omega * rho0 * (M @ vn)
    t0 = time.time()
    lu = sla.lu_factor(Aii, overwrite_a=True, check_finite=False)
    del Aii
    gc.collect()
    mus = sla.lu_solve(lu, rhs[ii], check_finite=False)
    del lu
    gc.collect()
    t_sol = time.time() - t0
    print(f"   LU+solve ({vn.shape[1]} rhs) in {t_sol:.1f}s", flush=True)
    res = []
    for j, lab in enumerate(D["labels"]):
        mu = np.zeros(nv, dtype=complex)
        mu[ii] = mus[:, j]
        P = core.SIGN_P * 0.5 * np.real(np.conj(vn[:, j]) @ (M @ mu))
        Pf = core.far_field_power_dl(grid, mu, k, rho0=rho0, c0=c0)
        res.append(dict(design=lab, P_bem=float(P), P_far=float(Pf)))
        print(f"   {lab:14s} P_bem={P:.6e} P_far={Pf:.6e} far/bem={Pf/P:.5f}", flush=True)
    out = dict(geom=geom, m=m, qbump=qbump, regular=int(q.regular), singular=int(q.singular), nv=nv, nt=nt,
               n_dofs=int(n), h_max=hmax, assembly_s=t_asm, solve_s=t_sol, normal_bad=bad, results=res)
    os.makedirs(RUNS, exist_ok=True)
    with open(os.path.join(RUNS, f"{geom}_m{m}_q{qbump}.json"), "w") as f:
        json.dump(out, f, indent=1)
    if m == 2 and qbump == 0:   # check same-mesh sinc with core's function vs export csv
        for j, lab in enumerate(D["labels"]):
            ps = core.sinc_power_mesh(grid, vn[:, j], k, rho0=rho0, c0=c0)
            print(f"   sinc_mesh (core) {lab}: {ps:.8e}", flush=True)


# ---------------------------------------------------------------------------- aggregation
def read_csv(path):
    import csv
    with open(path) as f:
        return list(csv.DictReader(f))


def aggregate():
    import csv, glob
    sinc_rows = read_csv(os.path.join(MESH, "export_sinc_check.csv"))
    sinc = {(r["geometry"], r["design"], int(r["m"])): float(r["P_sinc_mesh"]) for r in sinc_rows}
    ref = read_csv(os.path.join(OUT, "frozen_field_quadrature.csv"))
    refd = {}
    for r in ref:
        if int(r["surface_order"]) == 8 and int(r["angular_order"]) == 64:
            refd[(r["geometry"], r["design"])] = (float(r["legacy_centroid_W"]), float(r["consistent_spline_W"]))
    runs = [json.load(open(p)) for p in sorted(glob.glob(os.path.join(RUNS, "*.json")))]
    levels = []
    for r in runs:
        for x in r["results"]:
            ps = sinc[(r["geom"], x["design"], r["m"])]
            levels.append(dict(geometry=r["geom"], design=x["design"], m=r["m"], n_vertices=r["nv"], n_triangles=r["nt"],
                               n_dofs=r["n_dofs"], h_max=r["h_max"], quad_setting=f"reg{r['regular']}_sing{r['singular']}",
                               qbump=r["qbump"], P_bem=x["P_bem"], P_far=x["P_far"], far_over_bem=x["P_far"] / x["P_bem"],
                               P_sinc_mesh=ps, same_mesh_ratio=ps / x["P_bem"], assembly_s=r["assembly_s"], solve_s=r["solve_s"]))
    cols = ["geometry", "design", "m", "n_vertices", "n_triangles", "n_dofs", "h_max", "quad_setting", "P_bem", "P_far",
            "far_over_bem", "P_sinc_mesh", "same_mesh_ratio", "assembly_s", "solve_s"]
    with open(os.path.join(OUT, "bem_shell_levels.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for l in levels:
            w.writerow(l)
    np.savez(os.path.join(OUT, "bem_shell_levels.npz"),
             **{f"{l['geometry']}|{l['design']}|m{l['m']}|{l['quad_setting']}": np.array([l['P_bem'], l['P_far'], l['P_sinc_mesh'], l['h_max']]) for l in levels})

    def extrap(h, P):
        """returns dict of estimates from base-quad levels (sorted coarse->fine)."""
        out = {}
        h = np.asarray(h); P = np.asarray(P)
        out["finest"] = P[-1]
        if len(P) >= 2:
            out["linear_p1"] = P[-1] + (P[-1] - P[-2]) * h[-1] / (h[-2] - h[-1])
        if len(P) >= 3:
            d1, d2 = P[-2] - P[-3], P[-1] - P[-2]
            rr = (h[-2] / h[-3], h[-1] / h[-2])
            # observed order from three values assuming geometric h ratio (use mean ratio r)
            if d1 * d2 > 0 and abs(d2) < abs(d1) * 10:
                r = np.sqrt(rr[0] * rr[1])
                # solve (h2^p - h3^p)/(h1^p - h2^p) = d2/d1 numerically with actual h
                from scipy.optimize import brentq
                f = lambda p: (h[-2] ** p - h[-1] ** p) / (h[-3] ** p - h[-2] ** p) - d2 / d1
                try:
                    p = brentq(f, 0.05, 8.0)
                    out["order"] = p
                    out["richardson"] = P[-1] + (P[-1] - P[-2]) * h[-1] ** p / (h[-2] ** p - h[-1] ** p)
                except ValueError:
                    pass
            # linear extrapolation fits: lstsq P = a + b h with all points
            c = np.polyfit(h, P, 1)
            out["linfit_all"] = c[1]
            c2 = np.polyfit(h, P, 2)
            out["quadfit_all"] = c2[2]
        return out

    comp = []
    print("\n=== per-level table (base quadrature) ===")
    keys = sorted({(l["geometry"], l["design"]) for l in levels})
    for g, d in keys:
        L = sorted([l for l in levels if l["geometry"] == g and l["design"] == d and l["qbump"] == 0], key=lambda l: l["m"])
        if not L:
            continue
        h = [l["h_max"] for l in L]
        eb = extrap(h, [l["P_bem"] for l in L])
        es = extrap(h, [l["P_sinc_mesh"] for l in L])
        er = extrap(h, [l["same_mesh_ratio"] for l in L])
        Pc, Ps = refd[(g, d)]
        print(f"{g:9s} {d:20s} m={[l['m'] for l in L]} P_bem={['%.5e'%l['P_bem'] for l in L]} ratio_sm={['%.4f'%l['same_mesh_ratio'] for l in L]}")
        print(f"      bem extrap: " + ", ".join(f"{k}={v:.5e}" for k, v in eb.items()))
        print(f"      sinc extrap: " + ", ".join(f"{k}={v:.5e}" for k, v in es.items()) + f"  | consistent_spline={Ps:.5e}")
        print(f"      sinc-extrap vs consistent: " + ", ".join(f"{k}:{v/Ps-1:+.2e}" for k, v in es.items() if k != 'order'))
        print(f"      same-mesh ratio extrap: " + ", ".join(f"{k}={v:.5f}" for k, v in er.items()))
        # choose primary: Richardson if order in [0.5,2.5], else linear p1 from two finest
        cand = {k: v for k, v in eb.items() if k in ("richardson", "linear_p1", "linfit_all", "quadfit_all")}
        if "richardson" in eb and 0.5 <= eb["order"] <= 2.5:
            meth, Pe = "richardson", eb["richardson"]
        elif "linear_p1" in eb:
            meth, Pe = "linear_p1", eb["linear_p1"]
        else:
            meth, Pe = "finest", eb["finest"]
        spread = max([abs(v - eb["finest"]) for v in cand.values()] + [abs(Pe - eb["finest"])] + [abs(v - Pe) for v in cand.values()]) if cand else 0.0
        comp.append(dict(geometry=g, design=d, P_centroid_sinc=Pc, P_consistent_sinc=Ps, P_bem_finest=eb["finest"],
                         P_bem_extrap=Pe, extrap_method=meth, observed_order=eb.get("order", float("nan")),
                         P_bem_uncertainty=spread, ratio_centroid_over_bem=Pc / Pe, ratio_consistent_over_bem=Ps / Pe,
                         ratio_centroid_over_bem_finest=Pc / eb["finest"], ratio_consistent_over_bem_finest=Ps / eb["finest"],
                         same_mesh_ratio_finest=L[-1]["same_mesh_ratio"], same_mesh_ratio_extrap=er.get("richardson", er.get("linear_p1", float("nan")))))
    cc = list(comp[0].keys())
    with open(os.path.join(OUT, "bem_comparison.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cc)
        w.writeheader()
        for c in comp:
            w.writerow(c)
    print("\n=== quadrature sensitivity (m=2) ===")
    for g, d in keys:
        a = [l for l in levels if l["geometry"] == g and l["design"] == d and l["m"] == 2]
        b0 = [l for l in a if l["qbump"] == 0]; b1 = [l for l in a if l["qbump"] > 0]
        if b0 and b1:
            print(f"{g:9s} {d:20s} base {b0[0]['P_bem']:.6e}  raised({b1[0]['quad_setting']}) {b1[0]['P_bem']:.6e}  rel.change {b1[0]['P_bem']/b0[0]['P_bem']-1:+.3e}")


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]) if len(sys.argv) > 4 else 0)
    else:
        aggregate()
