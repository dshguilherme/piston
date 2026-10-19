"""Frozen cap field (m=2 mesh): BEM power vs frequency sweep. Resumable; appends CSV rows incrementally.
usage: python -u r0_viz_cap_freq_sweep.py
"""
import sys, os, time, gc, csv
import numpy as np
import scipy.linalg as sla

ROOT = r"C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course"
sys.path.insert(0, os.path.join(ROOT, "src"))
import bempp_cl.api as bem
import r0_bem_core as core
import r0_bem_shells as sh

VIZ = os.path.join(ROOT, "results", "r0_revision", "viz")
CSV = os.path.join(VIZ, "cap_m2_freq_sweep.csv")
COLS = ["kR", "k", "freq_Hz", "design", "P_bem", "P_far", "far_over_bem", "P_sinc_mesh", "ratio_sinc_over_bem",
        "P_plus", "P_minus", "elements_per_wavelength", "resolved", "assembly_s"]
KRS = [0.1, 0.3, 1, 2, 4, 6, 8, 10, 12, 14, 16, 18.3183, 21]
KW_SET = {1, 4, 10, 16}
RC = 10.0
C0 = 343.0


def done_kRs():
    if not os.path.exists(CSV):
        return set()
    with open(CSV) as f:
        return {round(float(r["kR"]), 6) for r in csv.DictReader(f)}


def main():
    os.makedirs(VIZ, exist_ok=True)
    D = sh.load("spherecap", 2)
    rho0, c0 = D["rho0"], D["c0"]
    vn = D["vn"]
    grid = bem.Grid(np.ascontiguousarray(D["V"].T), np.ascontiguousarray(D["T"].T.astype(np.uint32)))
    nv = grid.number_of_vertices
    edges = np.concatenate([D["V"][D["T"][:, [0, 1]]], D["V"][D["T"][:, [1, 2]]], D["V"][D["T"][:, [2, 0]]]])
    hmax = float(np.linalg.norm(edges[:, 0] - edges[:, 1], axis=1).max())
    print(f"nv={nv} h_max={hmax:.4f}  kR_max(6 el/wl) = {2*np.pi*RC/(6*hmax):.4f}  freq in file={D['freq']}", flush=True)
    bnd = core.boundary_vertices(grid)
    ii = np.where(~bnd)[0]
    M = core.mass_matrix_vertex(grid)
    space, vmap = core.full_p1(grid)
    di = vmap[ii]
    new = not os.path.exists(CSV)
    fh = open(CSV, "a", newline="")
    w = csv.writer(fh)
    if new:
        w.writerow(COLS); fh.flush()
    T0 = time.time()
    for kR in KRS:
        if round(float(kR), 6) in done_kRs():
            print(f"skip kR={kR} (done)", flush=True); continue
        k = kR / RC
        omega = k * c0
        freq = omega / (2 * np.pi)
        epw = (2 * np.pi / k) / hmax
        resolved = bool(epw >= 6)
        print(f"== kR={kR} k={k:.6f} f={freq:.4f} Hz  el/wl={epw:.3f} resolved={resolved}", flush=True)
        t0 = time.time()
        op = bem.operators.boundary.helmholtz.hypersingular(space, space, space, k)
        A = np.asarray(op.weak_form().A)
        Aii = A[np.ix_(di, di)]
        del A, op
        gc.collect()
        Aii *= core.SIGN_SCREEN
        t_asm = time.time() - t0
        print(f"   assembled in {t_asm:.1f}s", flush=True)
        rhs = 1j * omega * rho0 * (M @ vn)
        lu = sla.lu_factor(Aii, overwrite_a=True, check_finite=False)
        del Aii; gc.collect()
        mus = sla.lu_solve(lu, rhs[ii], check_finite=False)
        del lu; gc.collect()
        MU = np.zeros((nv, vn.shape[1]), dtype=complex)
        MU[ii, :] = mus
        Kw = None
        if any(abs(kR - s) < 1e-9 for s in KW_SET):
            op = bem.operators.boundary.helmholtz.double_layer(space, space, space, k)
            Kd = np.asarray(op.weak_form().A)
            Kw = Kd[np.ix_(vmap, vmap)]
            del Kd, op; gc.collect()
        rows = []
        for j, lab in enumerate(D["labels"]):
            mu = MU[:, j]
            Mmu = M @ mu
            P = core.SIGN_P * 0.5 * np.real(np.conj(vn[:, j]) @ Mmu)
            Pf = core.far_field_power_dl(grid, mu, k, rho0=rho0, c0=c0)
            if kR <= 0.5:
                Pf2 = core.far_field_power_dl(grid, mu, k, nth=192, nph=384, rho0=rho0, c0=c0)
                print(f"   [fine-grid far check] {lab:14s} P_far(192x384)={Pf2:.6e} far_fine/bem={Pf2/P:.5f} (default grid P_far={Pf:.6e})", flush=True)
            Ps = core.sinc_power_mesh(grid, vn[:, j], k, rho0=rho0, c0=c0)
            Pp = Pm = ""
            if Kw is not None:
                Kmu = Kw @ mu
                Pp = float(0.5 * np.real(np.conj(vn[:, j]) @ (Kmu + 0.5 * Mmu)))
                Pm = float(0.5 * np.real(np.conj(vn[:, j]) @ (0.5 * Mmu - Kmu)))
                print(f"   {lab:14s} P+={Pp:.6e} P-={Pm:.6e} sum/P_bem={(Pp+Pm)/P:.8f}", flush=True)
            print(f"   {lab:14s} P_bem={P:.6e} P_far={Pf:.6e} far/bem={Pf/P:.5f} sinc/bem={Ps/P:.5f}", flush=True)
            rows.append([kR, k, freq, lab, float(P), float(Pf), float(Pf / P), float(Ps), float(Ps / P), Pp, Pm, epw,
                         int(resolved), t_asm])
        for r in rows:
            w.writerow(r)
        fh.flush()
        del Kw, MU, mus
        gc.collect()
        print(f"   kR done; elapsed total {time.time()-T0:.0f}s", flush=True)
    fh.close()
    print(f"ALL DONE total {time.time()-T0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
