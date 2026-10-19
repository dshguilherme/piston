# REPRODUCE-PACKAGE COPY. Do not edit numerics here; edit the original and re-copy.
# Original: tests/r0_opt_cap_checks.py  (sha256 6599bbcba94df8e7)
# Changed: original lines 11-13 (HERE, sys.path.insert(0,HERE), "from r0_opt_exact_factor import OUT, RHO0, C0, load_samples"):
#          sys.path gets C.TESTS and the ORIGINAL r0_opt_exact_factor is imported for RHO0, C0, load_samples only;
#          OUT is replaced by EXO = C.EXO (stored inputs, read only; original lines 25, 33, 34, 44 read via EXO) and
#          W = C.out("smoke","r0_opt_cap_checks") i.e. $REPRO_OUT/smoke/r0_opt_cap_checks (original line 15: cap_reopt_checks.log is written there).
#          Added the config sys.path insert and "import config as C". The log is also printed to stdout, as in the original.
# Notes: load_samples(8) (original module) reads exact_opt/cap_samples_q8.mat from the original module's own OUT = <ROOT>/results/r0_revision/exact_opt
#        (the same folder as C.EXO when tests/ lives under REPRO_ROOT) and writes nothing. Importing r0_opt_exact_factor imports r0_nov_cap_frozen,
#        which executes os.makedirs(results/r0_revision/novelty, exist_ok=True) at import (an existing folder; creates nothing when the stored results are present).
# Reads: results/r0_revision/exact_opt/{run_S0,run_X0,run_XS,run_SX,run_SR,run_XR}.mat, cap_path_powers.npz, cap_path_states.mat (154 MB),
#        cap_exact_factor_q8.mat (160 MB), cap_exact_factor_q12.mat (160 MB), cap_samples_q8.mat
"""Robustness checks for the cap re-optimisation (r0_opt_cap_reopt.m / r0_opt_cap_eval.py).
1. convergence tails of each run (relative change of J over the last 50 / 100 outer iterations);
2. acoustic-quadrature check: exact power of the run optima with the q12 factor vs the q8 factor;
3. exact radiation efficiency sigma = P / (1/2 rho c int |v|^2 dS) of uniform and optima.
Writes results/r0_revision/exact_opt/cap_reopt_checks.log
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir))
import config as C
import numpy as np
import scipy.io as sio

sys.path.insert(0, C.TESTS)
from r0_opt_exact_factor import RHO0, C0, load_samples

EXO = C.EXO
W = C.out("smoke", "r0_opt_cap_checks")
logf = open(os.path.join(W, "cap_reopt_checks.log"), "w")


def log(s=""):
    print(s, flush=True); logf.write(s + "\n")


runs = ["S0", "X0", "XS", "SX", "SR", "XR"]
log("1. convergence tails (J under the run's own operator)")
for r in runs:
    d = sio.loadmat(os.path.join(EXO, f"run_{r}.mat"), squeeze_me=True, struct_as_record=False)
    J = np.atleast_1d(d["hist"].J); n = J.size
    jb = np.minimum.accumulate(J)
    c50 = (jb[-51] - jb[-1]) / jb[-1] if n > 51 else np.nan
    c100 = (jb[-101] - jb[-1]) / jb[-1] if n > 101 else np.nan
    log(f"  {r}: niter {n}, stop {d['meta'].stopReason}, Jbest {jb[-1]:.6e}, best-so-far decrease over last 50 its "
        f"{c50:.2e} ({10*np.log10(1+c50) if np.isfinite(c50) else np.nan:.4f} dB), last 100 its {c100:.2e} "
        f"({10*np.log10(1+c100) if np.isfinite(c100) else np.nan:.4f} dB)")
P = np.load(os.path.join(EXO, "cap_path_powers.npz"))
st = sio.loadmat(os.path.join(EXO, "cap_path_states.mat"), variable_names=["U", "run", "iter"])
run = st["run"].ravel().astype(int); it = st["iter"].ravel().astype(int)
sel = [np.where((run == k + 1) & (it == 0))[0][0] for k in range(len(runs))] + [int(P["iu"]), int(P["ir"])]
names = runs + ["uniform", "paper quadvel-rad"]
U = st["U"][:, sel]
del st
w0 = 2 * np.pi * 100.0
log("\n2. quadrature check, exact power of optima: q12 factor / q8 factor")
P8 = {}
for nq in (8, 12):
    F = sio.loadmat(os.path.join(EXO, f"cap_exact_factor_q{nq}.mat"), variable_names=["Fx"])["Fx"]
    P8[nq] = w0 ** 2 * np.sum(np.abs(F @ U) ** 2, 0)
    del F
for n, a, b in zip(names, P8[8], P8[12]):
    log(f"  {n:18s} q8 {a:.7e}  q12 {b:.7e}  rel {b/a-1:+.2e}")
log("\n3. exact radiation efficiency sigma = P_exact / (1/2 rho c w0^2 int |u_n|^2 dS)")
S = load_samples(8)
un = S["B"] @ U
q = 0.5 * RHO0 * C0 * w0 ** 2 * (S["w"][:, None] * np.abs(un) ** 2).sum(0)
for n, p, qq in zip(names, P8[8], q):
    log(f"  {n:18s} P_exact {p:.5e}  P_ref {qq:.5e}  sigma {p/qq:.4e}")
log("R0_OPT_CAP_CHECKS_DONE")
