"""Compare the in-flight bempp low-frequency sweep (results/r0_revision/viz/cap_m2_freq_sweep.csv, another agent's
file, READ ONLY) with the spectral two-sided reference (results/r0_revision/novelty/cap_frozen_sweep.csv).
Note: the BEM uses the m=2 mesh trace, the spectral sweep the m=8 trace (at 100 Hz the trace effect is 0.1-0.9 %).
Output: results/r0_revision/novelty/bem_lowk_check.csv (+ stdout)."""
import os, csv

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
bem = list(csv.DictReader(open(os.path.join(ROOT, "results", "r0_revision", "viz", "cap_m2_freq_sweep.csv"))))
spec = {(round(float(r["kR"]), 6), r["design"]): r
        for r in csv.DictReader(open(os.path.join(ROOT, "results", "r0_revision", "novelty", "cap_frozen_sweep.csv")))}
rows = []
for b in bem:
    key = (round(float(b["kR"]), 6), b["design"])
    if key not in spec:
        continue
    s = spec[key]
    P2 = float(s["P_2s"])
    rows.append(dict(kR=key[0], design=key[1], P_bem_surface=float(b["P_bem"]), P_bem_far=float(b["P_far"]),
                     P_spectral=P2, surface_over_spectral=float(b["P_bem"]) / P2, far_over_spectral=float(b["P_far"]) / P2,
                     far_over_surface_bem=float(b["far_over_bem"])))
    r = rows[-1]
    print(f"kR {r['kR']:5.2f} {r['design']:13s} BEM surface/spectral {r['surface_over_spectral']:8.3f}   "
          f"BEM far-field/spectral {r['far_over_spectral']:6.3f}   (BEM far/surface {r['far_over_surface_bem']:.3f})")
with open(os.path.join(ROOT, "results", "r0_revision", "novelty", "bem_lowk_check.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
