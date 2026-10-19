"""Flat-disc limit of the spectral cap solver: uniform v_n on caps of rim radius a = 1 with a/R -> 0, compared with the
independent Hankel-domain disc solution (tests/r0_nov_disc_hankel.py, results/r0_revision/novelty/disc_two_sided.csv).
L is scaled as L ~ 1/theta0 so that L*theta0 is fixed; Richardson in 1/L from (L/2, L).
Output: results/r0_revision/novelty/cap_flat_limit.csv (+ stdout)."""
import os, sys, csv, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from r0_nov_capspec import CapBasis, solve_cap, RHO0, C0

OUT = os.path.join(os.path.dirname(HERE), "results", "r0_revision", "novelty")
rows = []
for a_R in (0.4, 0.2, 0.1, 0.05):
    R = 1.0 / a_R
    th = np.arcsin(a_R)
    L = int(round(1600 / th))
    t = time.time()
    bas = CapBasis(th, 12, L, [0])
    b = R * R * 2 * np.pi * (bas.profile(bas.s, 0) * bas.wu).sum(axis=1)
    for ka in (0.1, 1.0, 3.0, 10.0):
        k = ka  # a = 1
        P1 = solve_cap(bas, R, k, {0: b}, L=L // 2)
        P2 = solve_cap(bas, R, k, {0: b}, L=L)
        Px = 2 * P2 - P1
        dip = 8 / (27 * np.pi) * RHO0 * C0 * ka ** 4
        rows.append(dict(a_over_R=a_R, theta0=th, L=L, ka=ka, P_two_sided=Px, P_over_dipole=Px / dip,
                         L_rel_change=(P2 - P1) / Px, P_over_rhoc_pi_a2=Px / (RHO0 * C0 * np.pi)))
        print(f"a/R={a_R:5.3f} L={L:6d} ka={ka:5.1f}  P2s/P_dipole={Px/dip:.6f}  P2s/(rho c pi a^2)={Px/(RHO0*C0*np.pi):.5f}"
              f"  L-change {(P2-P1)/Px:+.1e}  ({time.time()-t:.0f}s)", flush=True)
with open(os.path.join(OUT, "cap_flat_limit.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
