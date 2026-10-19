"""Cap m=2 100 Hz: BEM face decomposition, far field, sinc pattern; saves npz + numbers csv."""
import sys, os, json, time, gc, csv
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import scipy.linalg as sla

ROOT = r"C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course"
sys.path.insert(0, os.path.join(ROOT, "src"))
import bempp_cl.api as bem
import r0_bem_core as core
import r0_bem_shells as shells

VIZ = os.path.join(ROOT, "results", "r0_revision", "viz")
NPZ = os.path.join(VIZ, "cap_m2_100hz_bem.npz")
CSVO = os.path.join(VIZ, "numbers_cap_100hz.csv")
JSONF = os.path.join(ROOT, "results", "r0_revision", "bem_shell_runs", "spherecap_m2_q0.json")
LEVCSV = os.path.join(ROOT, "results", "r0_revision", "bem_shell_levels.csv")
T00 = time.time()

D = shells.load("spherecap", 2)
V, T, vn, N, labels = D["V"], D["T"], D["vn"], D["N"], D["labels"]
omega, rho0, c0 = D["omega"], D["rho0"], D["c0"]
k = omega / c0
print("freq", D["freq"], "k", k, "kR", k * 10.0, labels, flush=True)
grid = bem.Grid(np.ascontiguousarray(V.T), np.ascontiguousarray(T.T.astype(np.uint32)))
nv, nt = grid.number_of_vertices, grid.number_of_elements
edges = np.concatenate([V[T[:, [0, 1]]], V[T[:, [1, 2]]], V[T[:, [2, 0]]]])
elen = np.linalg.norm(edges[:, 0] - edges[:, 1], axis=1)
hmax = float(elen.max())
lam = c0 / D["freq"]
print("nv", nv, "nt", nt, "hmax", hmax, "lambda/hmax", lam / hmax, flush=True)

space, vmap = core.full_p1(grid)
M = core.mass_matrix_vertex(grid)
bnd = core.boundary_vertices(grid)
ii = np.where(~bnd)[0]
rhs = 1j * omega * rho0 * (M @ vn)

# ---- hypersingular solve
CKPT = os.path.join(VIZ, "cap_m2_100hz_mu_checkpoint.npy")
if os.path.exists(CKPT):
    mu = np.load(CKPT); t_H = 784.9   # first run's H assembly time (checkpoint reuse after OOM crash in far field)
    print("loaded mu checkpoint (H assembly 784.9 s in first run)", flush=True)
else:
    t0 = time.time()
    op = bem.operators.boundary.helmholtz.hypersingular(space, space, space, k)
    A = np.asarray(op.weak_form().A)
    di = vmap[ii]
    Aii = A[np.ix_(di, di)] * core.SIGN_SCREEN
    del A, op; gc.collect()
    t_H = time.time() - t0
    print(f"H assembled {t_H:.1f}s", flush=True)
    mus = sla.solve(Aii, rhs[ii], check_finite=False)
    del Aii; gc.collect()
    mu = np.zeros((nv, vn.shape[1]), dtype=complex)
    mu[ii] = mus
    np.save(CKPT, mu)

# ---- double-layer weak form (stage 1 saves Kw@mu and exits to release bempp memory; stage 2 reloads)
KCK = os.path.join(VIZ, "cap_m2_100hz_Kmu_checkpoint.npy")
t_K = 81.3
if os.path.exists(KCK):
    Kw_mu = np.load(KCK)
    print("loaded Kw@mu checkpoint", flush=True)
else:
    t0 = time.time()
    opK = bem.operators.boundary.helmholtz.double_layer(space, space, space, k)
    Kd = np.asarray(opK.weak_form().A)
    mu_dof = np.zeros_like(mu)
    mu_dof[vmap] = mu                      # dof order
    Kmu_dof = Kd @ mu_dof
    Kw_mu = Kmu_dof[vmap]                  # back to vertex order (rows)
    np.save(KCK, Kw_mu)
    print(f"K assembled+applied {time.time()-t0:.1f}s; saved checkpoint; exiting stage 1", flush=True)
    sys.exit(10)
P_bem = 0.5 * np.real(np.einsum("ij,ij->j", np.conj(vn), M @ mu))
Mmu = M @ mu
KvnS = 0.5 * np.real(np.einsum("ij,ij->j", np.conj(vn), Kw_mu))     # 0.5 Re(vn^H Kw mu)
Mvn_mu = 0.5 * np.real(np.einsum("ij,ij->j", np.conj(vn), 0.5 * Mmu))

js = json.load(open(JSONF))
pj = {r["design"]: r for r in js["results"]}
P_json = np.array([pj[l]["P_bem"] for l in labels])
Pfar_json = np.array([pj[l]["P_far"] for l in labels])
print("regression P_bem vs json rel:", P_bem / P_json - 1, flush=True)

# ---- jump / sign self-check
x_nodal_all = np.zeros_like(mu)
lu = spla.splu(M.tocsc())
for j in range(mu.shape[1]):
    x_nodal_all[:, j] = lu.solve(np.real(Kw_mu[:, j])) + 1j * lu.solve(np.imag(Kw_mu[:, j]))

vh = np.zeros(nv); cnt = np.zeros(nv)
for a, b in ((0, 1), (1, 2), (2, 0)):
    l = np.linalg.norm(V[T[:, a]] - V[T[:, b]], axis=1)
    for c in (a, b):
        np.add.at(vh, T[:, c], l); np.add.at(cnt, T[:, c], 1)
hloc = vh / cnt
polar = np.degrees(np.arccos(V[:, 2] / np.linalg.norm(V, axis=1)))
cand = np.where((polar < 20.0) & (~bnd))[0]
rng = np.random.default_rng(0)
sel = np.sort(rng.choice(cand, 20, replace=False))
print("jump check vertices:", len(sel), "polar range", polar[sel].min(), polar[sel].max(), flush=True)
jump_stats = {}
check_designs = ["quadvel-rad", "uniform"]
jump_sign = {}
for lab in check_designs:
    j = labels.index(lab)
    dof = np.zeros(nv, dtype=complex); dof[vmap] = mu[:, j]
    gf = bem.GridFunction(space, coefficients=dof)
    for fac in (0.5, 0.25):
        dl = fac * hloc[sel]
        res = {}
        for sgn in (+1, -1):
            pts = (V[sel] + sgn * dl[:, None] * N[sel]).T
            pot = bem.operators.potential.helmholtz.double_layer(space, np.ascontiguousarray(pts), k)
            res[sgn] = np.asarray(pot.evaluate(gf)).ravel()
        jmp = res[1] - res[-1]
        avg = 0.5 * (res[1] + res[-1])
        mu_s = mu[sel, j]; x_s = x_nodal_all[sel, j]
        e_j = np.linalg.norm(jmp - mu_s) / np.linalg.norm(mu_s)
        e_j_neg = np.linalg.norm(jmp + mu_s) / np.linalg.norm(mu_s)
        e_a = np.linalg.norm(avg - x_s) / np.linalg.norm(x_s)
        coef = np.vdot(mu_s, jmp) / np.vdot(mu_s, mu_s)
        jump_stats[(lab, fac)] = (e_j, e_j_neg, e_a, coef)
        print(f"  {lab} delta={fac}h: |jump-mu|/|mu|={e_j:.4e}  |jump+mu|/|mu|={e_j_neg:.4e}  "
              f"|avg-x|/|x|={e_a:.4e}  LS coef jump/mu={coef:.4f}", flush=True)
coefs = [v[3].real for v in jump_stats.values()]
S = 1.0 if np.mean(coefs) > 0 else -1.0
print("JUMP SIGN s =", S, "(+1 means p(+n)-p(-n)=mu as assumed; -1 means plus/minus SWAPPED)", flush=True)

# ---- face decomposition
Pp = 0.5 * np.real(np.einsum("ij,ij->j", np.conj(vn), Kw_mu + S * 0.5 * Mmu))
Pm = 0.5 * np.real(np.einsum("ij,ij->j", np.conj(vn), S * 0.5 * Mmu - Kw_mu))
print("P_plus+P_minus vs P_bem rel:", (Pp + Pm) / P_bem - 1, flush=True)
p_plus = x_nodal_all + S * mu / 2
p_minus = x_nodal_all - S * mu / 2
I_plus = 0.5 * np.real(p_plus * np.conj(vn))
I_minus = -0.5 * np.real(p_minus * np.conj(vn))

# ---- far field
def bem_far(mu_j, dirs):
    dof = np.zeros(nv, dtype=complex); dof[vmap] = mu_j
    gfn = bem.GridFunction(space, coefficients=dof)
    out = []
    for s in range(0, len(dirs), 800):      # chunked: the far-field operator is dense (ndirs x ndofs)
        ff = bem.operators.far_field.helmholtz.double_layer(space, np.ascontiguousarray(dirs[s:s + 800].T), k)
        out.append(np.asarray(ff.evaluate(gfn)).ravel()); del ff; gc.collect()
    return np.concatenate(out)

Xq, Wq, _ = core.surface_samples(grid, vn[:, 0])
def sinc_far(vj, dirs, chunk=200):
    _, _, vq = core.surface_samples(grid, vj)
    wv = Wq * vq
    out = np.zeros(len(dirs), dtype=complex)
    for s in range(0, len(dirs), chunk):
        ph = np.exp(-1j * k * (dirs[s:s + chunk] @ Xq.T))
        out[s:s + chunk] = ph @ wv
    return rho0 * omega / (2 * np.sqrt(2) * np.pi) * np.abs(out)

L = max(np.ptp(grid.vertices, axis=1))
n = int(max(96, 4 * k * L + 24))
dirsG, WG = core._sphere_dirs(n, 2 * n)
th = np.deg2rad(np.arange(0, 181, 1.0)); ph = np.deg2rad(np.arange(0, 361, 2.0))
TH, PH = np.meshgrid(th, ph, indexing="ij")
dirsP = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], axis=-1).reshape(-1, 3)
north = dirsG[:, 2] > 0
fac = 1.0 / (2 * rho0 * c0)

far = {}
sinc_mesh_all = {}
lev = {(r["design"], r["geometry"], r["m"]): r for r in csv.DictReader(open(LEVCSV))}
for j, lab in enumerate(labels):
    sinc_mesh_all[lab] = core.sinc_power_mesh(grid, vn[:, j], k, rho0=rho0, c0=c0, chunk=250)
    print(f"sinc_mesh {lab}: {sinc_mesh_all[lab]:.8e}", flush=True)
for lab in check_designs:
    j = labels.index(lab)
    t0 = time.time()
    pG = bem_far(mu[:, j], dirsG); pP = bem_far(mu[:, j], dirsP)
    sG = sinc_far(vn[:, j], dirsG); sP = sinc_far(vn[:, j], dirsP)
    d = dict(
        P_int_bem=fac * np.sum(WG * np.abs(pG) ** 2), P_int_sinc=fac * np.sum(WG * sG ** 2),
        P_north_bem=fac * np.sum(WG[north] * np.abs(pG[north]) ** 2), P_south_bem=fac * np.sum(WG[~north] * np.abs(pG[~north]) ** 2),
        P_north_sinc=fac * np.sum(WG[north] * sG[north] ** 2), P_south_sinc=fac * np.sum(WG[~north] * sG[~north] ** 2),
        dPdO_bem=(np.abs(pP) ** 2 * fac).reshape(TH.shape), dPdO_sinc=(sP ** 2 * fac).reshape(TH.shape))
    far[lab] = d
    print(f"far {lab}: P_int_bem={d['P_int_bem']:.6e} (P_bem {P_bem[j]:.6e}, P_far_json {Pfar_json[j]:.6e}) "
          f"P_int_sinc={d['P_int_sinc']:.6e} (sinc_mesh {sinc_mesh_all[lab]:.6e})  {time.time()-t0:.0f}s", flush=True)

# ---- save npz
sv = dict(vertices=V, triangles=T, vn=vn, mu=mu, x=x_nodal_all, p_plus=p_plus, p_minus=p_minus, I_plus=I_plus, I_minus=I_minus,
          labels=np.array(labels), theta_deg=np.arange(0, 181, 1.0), phi_deg=np.arange(0, 361, 2.0),
          P_bem=P_bem, P_bem_json=P_json, P_far_json=Pfar_json, P_plus=Pp, P_minus=Pm,
          P_sinc_mesh=np.array([sinc_mesh_all[l] for l in labels]), jump_sign=S, kR=k * 10.0, k=k, freq=D["freq"],
          h_max=hmax, lam_over_hmax=lam / hmax, rho0=rho0, c0=c0, nv=nv, nt=nt, jump_sel=sel)
for lab in check_designs:
    key = lab.replace("-", "_")
    for kk, vv in far[lab].items():
        sv[f"{kk}__{key}"] = vv
np.savez(NPZ, **sv)

# ---- numbers csv
rows = []
src_npz = "results/r0_revision/viz/cap_m2_100hz_bem.npz"
src_json = "results/r0_revision/bem_shell_runs/spherecap_m2_q0.json"
src_lev = "results/r0_revision/bem_shell_levels.csv"
def add(q, des, val, unit, src=src_npz, fig="fig2/fig3"):
    rows.append([fig, q, des, f"{val:.10e}" if isinstance(val, (float, np.floating)) else val, unit, src])
for j, lab in enumerate(labels):
    ps = sinc_mesh_all[lab]
    lv = lev.get((lab, "spherecap", "2"))
    add("P_bem", lab, P_bem[j], "W"); add("P_bem_json", lab, P_json[j], "W", src_json)
    add("P_bem_regression_rel_err", lab, P_bem[j] / P_json[j] - 1, "1", src_json)
    add("P_plus", lab, Pp[j], "W"); add("P_minus", lab, Pm[j], "W")
    add("P_plus/P_bem", lab, Pp[j] / P_bem[j], "1"); add("P_minus/P_bem", lab, Pm[j] / P_bem[j], "1")
    add("P_plus+P_minus_over_P_bem_minus1", lab, (Pp[j] + Pm[j]) / P_bem[j] - 1, "1")
    add("P_sinc_mesh", lab, ps, "W")
    if lv: add("P_sinc_mesh_csv", lab, float(lv["P_sinc_mesh"]), "W", src_lev)
    add("P_sinc_mesh/P_bem", lab, ps / P_bem[j], "1"); add("2*P_sinc_mesh/P_bem", lab, 2 * ps / P_bem[j], "1")
    add("P_sinc/P_plus", lab, ps / Pp[j], "1"); add("P_sinc/P_minus", lab, ps / Pm[j], "1")
for lab in check_designs:
    j = labels.index(lab); d = far[lab]
    add("P_far_integral_bem", lab, d["P_int_bem"], "W"); add("P_far_integral_bem_rel_err_vs_P_bem", lab, d["P_int_bem"] / P_bem[j] - 1, "1")
    add("P_far_integral_bem_rel_err_vs_P_far_json", lab, d["P_int_bem"] / Pfar_json[j] - 1, "1", src_json)
    add("P_far_integral_sinc", lab, d["P_int_sinc"], "W")
    add("P_far_integral_sinc_rel_err_vs_sinc_mesh", lab, d["P_int_sinc"] / sinc_mesh_all[lab] - 1, "1")
    if lev.get((lab, "spherecap", "2")):
        add("P_far_integral_sinc_rel_err_vs_csv", lab, d["P_int_sinc"] / float(lev[(lab, 'spherecap', '2')]["P_sinc_mesh"]) - 1, "1", src_lev)
    for m_ in ("bem", "sinc"):
        add(f"P_north_theta<90_{m_}", lab, d[f"P_north_{m_}"], "W"); add(f"P_south_theta>90_{m_}", lab, d[f"P_south_{m_}"], "W")
        add(f"north_fraction_{m_}", lab, d[f"P_north_{m_}"] / d[f"P_int_{m_}"], "1")
    for fac_ in (0.5, 0.25):
        e_j, e_jn, e_a, coef = jump_stats[(lab, fac_)]
        add(f"jump_check_rel_err_delta{fac_}h", lab, e_j, "1"); add(f"jump_check_rel_err_oppositesign_delta{fac_}h", lab, e_jn, "1")
        add(f"mean_check_rel_err_delta{fac_}h", lab, e_a, "1"); add(f"jump_LS_coef_delta{fac_}h", lab, float(coef.real), "1")
add("jump_sign_s", "all", float(S), "1"); add("kR", "all", k * 10.0, "1"); add("h_max", "all", hmax, "m")
add("lambda_over_h_max", "all", lam / hmax, "elements/wavelength"); add("nv", "all", float(nv), "1"); add("nt", "all", float(nt), "1")
with open(CSVO, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["figure", "quantity", "design", "value", "unit", "source_file"]); w.writerows(rows)
print(f"DONE total {time.time()-T00:.0f}s (H asm {t_H:.0f}s, K asm {t_K:.0f}s)", flush=True)


