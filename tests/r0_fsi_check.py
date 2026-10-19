"""Experiment 7 (FSI check), step 2: coupled structural-acoustic response of the regime-pilot designs.

Structure: the pilot's model, unchanged (element arrays exported by tests/r0_fsi_export_ctx.m from the same pilot_ctx):
    K(t) = sum_e Bke t^3 + Ske t,  M(t) = sum_e Me t,  K_D = (1 + i w beta) K - w^2 M,  C = beta K   (e^{+i w t})
Fluid: Z_plus = conj(Z_minus), Z_minus from tests/r0_fsi_zop.py (part 'Zc', Richardson L = 4000/8000), scaled by s:
    (K_D + i w s Z_plus) U = F          (dense complex LU on the free dofs)
    s = 0     : in vacuo (pilot)
    s = 1     : the pilot's device, rho' c' = rho0 c0 (air impedance), c' = w R / kR  ->  rho' = rho0 c0 / c'
    s = c'/c0 : physical air density rho0 = 1.21 with the same c' (impedance rho0 c'), the other reading of the device
Radiated power with the pilot's factor:  P = s w^2 ||Fx U||^2  (s = 1 for in vacuo and device rows).
Power balance: P_in = 1/2 Re(F^H iwU), P_diss = 1/2 w^2 beta U^H K U, P_rad(Z) = 1/2 w^2 s Re(U^H Z_plus U).

Designs: uniform t = 0.05 and the final designs tbest of the pilot runs XU, XS (exact), SU, SX (sinc), EU, EX (ERP);
for the anchor c30_k18 the XU/SU/SX/XS runs are the Experiment-5 runs X0/S0/SX/XS (as in the pilot).

Run:  python -u tests/r0_fsi_check.py <cell> [...]
Writes results/r0_revision/fsi/fsi_<cell>.csv, fsi_sweep_<cell>.csv, fsi_modes_<cell>.csv, fsi_<cell>.log
"""
import os, sys, time
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
import numpy as np
import scipy.io as sio
import scipy.sparse as sp
import scipy.linalg as sla
import csv

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from r0_fsi_zop import load_zop, assemble, OUT
from r0_pilot_factor import CELLS, load_samples, R
from r0_nov_capspec import RHO0, C0

ROOT = os.path.dirname(HERE)
PIL = os.path.join(ROOT, "results", "r0_revision", "pilot")
EXP5 = os.path.join(ROOT, "results", "r0_revision", "exact_opt")
RUNS = ["XU", "XS", "SU", "SX", "EU", "EX"]
SWEEP = [1e-3, 1e-2, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0]
NMODES = 40


def runfile(cell, code):
    m = {"XU": "X0", "SU": "S0", "SX": "SX", "XS": "XS"}
    if cell == "c30_k18" and code in m:
        return os.path.join(EXP5, f"run_{m[code]}.mat")
    return os.path.join(PIL, f"run_{cell}_{code}.mat")


class Struct:
    def __init__(self, th_deg):
        c = sio.loadmat(os.path.join(OUT, f"ctx_th{th_deg}.mat"))
        self.Bke, self.Ske, self.Me = c["Bke"], c["Ske"], c["Me"]
        self.rows = (c["Rows"].astype(np.int64) - 1).ravel(); self.cols = (c["Cols"].astype(np.int64) - 1).ravel()
        self.free = c["free"].ravel().astype(np.int64) - 1
        self.F = c["F"].ravel(); self.w = float(c["w"].ravel()[0]); self.beta = float(c["beta"].ravel()[0]); self.ndof = int(c["ndof"].ravel()[0])
        self.Uu = c["Uu"].ravel()

    def KM(self, t):
        t = t.ravel()[:, None]
        K = sp.coo_matrix(((self.Bke * t ** 3 + self.Ske * t).ravel(), (self.rows, self.cols)), (self.ndof, self.ndof)).tocsr()
        Mm = sp.coo_matrix(((self.Me * t).ravel(), (self.rows, self.cols)), (self.ndof, self.ndof)).tocsr()
        return K, Mm


def main(cells):
    for cell in cells:
        t0 = time.time()
        logf = open(os.path.join(OUT, f"fsi_{cell}.log"), "w")

        def log(s):
            print(s, flush=True); logf.write(s + "\n"); logf.flush()
        th_deg, kR = CELLS[cell]
        cp = 2 * np.pi * 100.0 * R / kR
        sphys = cp / C0
        log(f"cell {cell}: theta0 {th_deg}, kR {kR:.4f}, c' {cp:.2f} m/s; device rho' = {RHO0*C0/cp:.5f} kg/m^3 (s=1); "
            f"physical-density variant rho = {RHO0} kg/m^3, impedance x{sphys:.4f} (s={sphys:.4f})")
        st = Struct(th_deg); w = st.w; fr = st.free
        S = load_samples(th_deg); Bs, ws = S["B"], S["w"]
        Fx = sio.loadmat(os.path.join(PIL, f"factor_{cell}.mat"), variable_names=["Fx"])["Fx"]
        d = load_zop(cell)
        Zp = np.conj(assemble(d, "x", "Zc"))[np.ix_(fr, fr)]          # e^{+iwt}
        Zp8 = np.conj(assemble(d, "8000", "Zc"))[np.ix_(fr, fr)]
        Hx = {"x": None}
        Q, offs = d["Q"], d["offs"]
        Hxb = [2 * d["Hb8000"][i] - d["Hb4000"][i] for i in range(len(offs) - 1)]
        H8b = list(d["Hb8000"])
        Ma_unit = Zp.imag / w                                            # added mass per unit s (fixed-frequency)
        log(f"  Z assembled; ||Herm part|| / ||anti-Herm part|| = {np.linalg.norm(Zp.real)/np.linalg.norm(Zp.imag):.3e} "
            f"(e^{{+iwt}}: Re Z = resistance, Im Z = w * added mass); Im-part eigenvalue range "
            f"[{np.linalg.eigvalsh(0.5*(Ma_unit+Ma_unit.T))[[0,-1]]}] kg (unit s)")

        def pquad(blocks, U):
            QU = Q @ U
            return w ** 2 * sum(np.real(np.vdot(QU[offs[i]:offs[i+1]], blocks[i] @ QU[offs[i]:offs[i+1]]))
                                for i in range(len(blocks)))

        def capfrac(U):
            v = 1j * w * U; QV = Q @ v
            return (np.sum(np.abs(QV) ** 2) / (2 * np.pi * R * R)) / np.sum(ws * np.abs(Bs @ v) ** 2)

        # designs
        designs = [("uniform", 0.05 * np.ones(960))]
        for r in RUNS:
            fn = runfile(cell, r)
            if os.path.exists(fn):
                designs.append((r, sio.loadmat(fn, variable_names=["tbest"])["tbest"].ravel()))
        pil = {}
        with open(os.path.join(PIL, "pilot_runs.csv")) as f:
            for row in csv.DictReader(f):
                if row["cell"] == cell:
                    pil[row["run"]] = float(row["P_exact"]); pil["uniform"] = float(row["P_exact_uniform"])
        rows, srows, mrows = [], [], []
        Pref = {}
        for name, t in designs:
            K, Mm = st.KM(t)
            Kf = K[fr][:, fr].toarray()
            Mf = Mm[fr][:, fr].toarray()
            KD = (1 + 1j * w * st.beta) * Kf - w ** 2 * Mf
            Ff = st.F[fr]
            res = {}
            for var, s, Zuse in (("vac", 0.0, Zp), ("dev", 1.0, Zp), ("phys", sphys, Zp), ("dev_L8000", 1.0, Zp8)):
                if var == "phys" and abs(sphys - 1) < 1e-9:
                    res[var] = res["dev"]; continue
                A = KD + 1j * w * s * Zuse if s else KD
                Uf = sla.solve(A, Ff)
                U = np.zeros(st.ndof, complex); U[fr] = Uf
                sp_ = s if var == "phys" else 1.0              # impedance scale for power in this medium
                Pfx = sp_ * w ** 2 * np.sum(np.abs(Fx @ U) ** 2)
                Phx = sp_ * pquad(Hxb, U); Ph8 = sp_ * pquad(H8b, U)
                Pin = 0.5 * np.real(np.vdot(Ff, 1j * w * Uf))
                Pdiss = 0.5 * w ** 2 * st.beta * np.real(np.vdot(Uf, Kf @ Uf))
                PradZ = 0.5 * w ** 2 * s * np.real(np.vdot(Uf, Zuse @ Uf)) if s else 0.0
                Perp = sp_ * 0.5 * RHO0 * C0 * w ** 2 * np.sum(ws * np.abs(Bs @ U) ** 2)
                KE = 0.25 * w ** 2 * np.real(np.vdot(Uf, Mf @ Uf))
                KEa = 0.25 * w ** 2 * np.real(np.vdot(Uf, Ma_unit @ Uf)) * (sp_ if var == "phys" else 1.0)
                res[var] = dict(U=U, Pfx=Pfx, Phx=Phx, Ph8=Ph8, Pin=Pin, Pdiss=Pdiss, PradZ=PradZ, Perp=Perp,
                                drive=st.F @ U, KE=KE, KEa=KEa, cap=capfrac(U))
            if name == "uniform":
                ucheck = np.linalg.norm(res["vac"]["U"] - st.Uu) / np.linalg.norm(st.Uu)
                log(f"  uniform in-vacuo state, Python dense vs MATLAB export: rel diff {ucheck:.2e}")
                Pref = {v: res[v]["Pfx"] for v in res}
            v0 = res["vac"]
            pp = pil.get(name, np.nan)
            for var in ("vac", "dev", "phys", "dev_L8000"):
                r_ = res[var]
                base = res["vac"]["Pfx"] * (sphys if var == "phys" else 1.0)
                rows.append(dict(cell=cell, theta0_deg=th_deg, kR=kR, design=name, variant=var,
                                 s=(0.0 if var == "vac" else (sphys if var == "phys" else 1.0)),
                                 P_rad_Fx=r_["Pfx"], P_rad_Hx_unclipped=r_["Phx"], P_rad_H8000=r_["Ph8"],
                                 P_rad_vac_same_medium=base, ratio_fsi_over_vac=r_["Pfx"] / base,
                                 dB_fsi_over_vac=10 * np.log10(r_["Pfx"] / base),
                                 P_ERP=r_["Perp"], sigma=r_["Pfx"] / r_["Perp"],
                                 dB_vs_uniform_same_variant=10 * np.log10(r_["Pfx"] / Pref[var]),
                                 P_in=r_["Pin"], P_diss=r_["Pdiss"], P_rad_Z=r_["PradZ"],
                                 balance_resid_rel_Pin=(r_["Pin"] - r_["Pdiss"] - r_["PradZ"]) / r_["Pin"],
                                 PradZ_over_PradFx=(r_["PradZ"] / r_["Pfx"] if r_["PradZ"] else np.nan),
                                 drive_point_ratio=abs(r_["drive"] / v0["drive"]),
                                 dU_rel=np.linalg.norm(r_["U"] - v0["U"]) / np.linalg.norm(v0["U"]),
                                 KE_struct=r_["KE"], KE_added_over_KE=r_["KEa"] / r_["KE"],
                                 vel_frac_in_fluid_basis=r_["cap"], P_exact_pilot_csv=pp))
            log(f"  {name:8s} vac P {v0['Pfx']:.4e} (pilot {pp:.4e})  dev P {res['dev']['Pfx']:.4e} "
                f"x{res['dev']['Pfx']/v0['Pfx']:.4g}  phys x{res['phys']['Pfx']/v0['Pfx']:.4g}  "
                f"L8000-coupling x{res['dev_L8000']['Pfx']/v0['Pfx']:.4g}  KEa/KE dev {res['dev']['KEa']/res['dev']['KE']:.3e}  "
                f"cap {v0['cap']:.6f}  bal {rows[-3]['balance_resid_rel_Pin']:.1e}  Hx/Fx {res['dev']['Phx']/res['dev']['Pfx']:.4f}")
            # fluid-loading sweep
            for s in sorted(set(SWEEP + [sphys])):
                Uf = sla.solve(KD + 1j * w * s * Zp, Ff)
                U = np.zeros(st.ndof, complex); U[fr] = Uf
                srows.append(dict(cell=cell, design=name, s=s, ratio_fsi_over_vac=np.sum(np.abs(Fx @ U) ** 2) /
                                  np.sum(np.abs(Fx @ v0["U"]) ** 2),
                                  P_rad_Fx_unitimpedance=w ** 2 * np.sum(np.abs(Fx @ U) ** 2)))
            # resonance shifts (fixed-frequency added mass at the drive kR), lowest NMODES modes
            if name in ("uniform", "XU", "XS", "SU"):
                Ks = 0.5 * (Kf + Kf.T); Ms = 0.5 * (Mf + Mf.T); Ma = 0.5 * (Ma_unit + Ma_unit.T)
                n = Ks.shape[0]
                fl = {}
                for var, s in (("vac", 0.0), ("dev", 1.0), ("phys", sphys)):
                    mu = sla.eigh(Ms + s * Ma, Ks, eigvals_only=True, subset_by_index=[n - NMODES, n - 1])
                    fl[var] = np.sort(1.0 / (2 * np.pi * np.sqrt(mu)))
                for i in range(NMODES):
                    mrows.append(dict(cell=cell, design=name, mode=i + 1, f_vac=fl["vac"][i], f_dev=fl["dev"][i],
                                      f_phys=fl["phys"][i], shift_dev_pct=100 * (fl["dev"][i] / fl["vac"][i] - 1),
                                      shift_phys_pct=100 * (fl["phys"][i] / fl["vac"][i] - 1)))
                sd = [100 * (fl["dev"][i] / fl["vac"][i] - 1) for i in range(NMODES)]
                sph = [100 * (fl["phys"][i] / fl["vac"][i] - 1) for i in range(NMODES)]
                log(f"    modes {name}: f1 vac {fl['vac'][0]:.3f} Hz; shift dev median {np.median(sd):.4f}% "
                    f"[{min(sd):.4f}, {max(sd):.4f}]; phys median {np.median(sph):.4f}% [{min(sph):.4f}, {max(sph):.4f}]")
            del KD, Kf, Mf
        for fn, data in ((f"fsi_{cell}.csv", rows), (f"fsi_sweep_{cell}.csv", srows), (f"fsi_modes_{cell}.csv", mrows)):
            with open(os.path.join(OUT, fn), "w", newline="") as f:
                wr = csv.DictWriter(f, fieldnames=list(data[0].keys())); wr.writeheader(); wr.writerows(data)
        log(f"  runtime {time.time()-t0:.0f} s")
        log("R0_FSI_CHECK_DONE")
        logf.close()


if __name__ == "__main__":
    main(sys.argv[1:])
