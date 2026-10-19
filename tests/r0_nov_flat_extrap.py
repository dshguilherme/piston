"""Richardson extrapolation in (a/R)^2 of the spectral-cap flat limit and comparison with the Hankel disc solution.
Reads results/r0_revision/novelty/cap_flat_limit.log (or .csv when present) and disc_two_sided.csv.
Output: results/r0_revision/novelty/cap_flat_extrap.csv (+ stdout)."""
import os, re, csv
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "r0_revision", "novelty")
RHO0, C0 = 1.21, 343.0
pat = re.compile(r"a/R=([\d.]+)\s+L=\s*(\d+)\s+ka=\s*([\d.]+)\s+P2s/P_dipole=([\d.eE+-]+)")
cap = {}
for line in open(os.path.join(OUT, "cap_flat_limit.log")):
    m = pat.search(line)
    if m:
        aR, L, ka, pd = float(m[1]), int(m[2]), float(m[3]), float(m[4])
        cap[(aR, ka)] = pd
disc = {round(float(r["ka"]), 6): float(r["ratio_two_sided_over_dipole_asymptote"])
        for r in csv.DictReader(open(os.path.join(OUT, "disc_two_sided.csv")))}
rows = []
aRs = sorted({k[0] for k in cap}, reverse=True)
for ka in sorted({k[1] for k in cap}):
    d = disc[round(ka, 6)]
    vals = [cap[(a, ka)] for a in aRs if (a, ka) in cap]
    a_used = [a for a in aRs if (a, ka) in cap]
    # extrapolate with the two smallest a/R values, assuming error ~ (a/R)^2
    a1, a2 = a_used[-2], a_used[-1]
    v1, v2 = vals[-2], vals[-1]
    r2 = (a1 / a2) ** 2
    vx = v2 + (v2 - v1) / (r2 - 1)
    rows.append(dict(ka=ka, aR_used=f"{a1},{a2}", cap_values=";".join(f"{a}:{v:.6f}" for a, v in zip(a_used, vals)),
                     cap_extrap=vx, disc_hankel=d, rel_diff=vx / d - 1,
                     diff_ratio_last=(vals[-3] - d) / (v1 - d) if len(vals) >= 3 else (v1 - d) / (v2 - d)))
    print(f"ka={ka:5.1f}  cap(a/R): {rows[-1]['cap_values']}  extrap {vx:.6f}  disc {d:.6f}  rel diff {vx/d-1:+.1e}")
with open(os.path.join(OUT, "cap_flat_extrap.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
