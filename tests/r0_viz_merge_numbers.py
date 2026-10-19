"""Merge per-figure number CSVs into results/r0_revision/viz/viz_numbers.csv and add derived
cross-design statistics for diagnostic 2 (face share vs sinc/BEM). Also checks PNG corner pixels."""
import os, csv, glob
import numpy as np
from PIL import Image

VIZ = r"C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course\results\r0_revision\viz"
COLS = ["figure", "quantity", "design", "geometry", "kR", "resolved", "value", "unit", "source_file", "origin_csv"]
out = []
for name in ["numbers_fig1.csv", "numbers_cap_100hz.csv", "numbers_fig4.csv"]:
    p = os.path.join(VIZ, name)
    if not os.path.exists(p):
        print("MISSING", name); continue
    with open(p) as f:
        for r in csv.DictReader(f):
            row = {c: r.get(c, "") for c in COLS}
            row["origin_csv"] = name
            out.append(row)

# diagnostic 2 cross-design statistics from numbers_cap_100hz.csv
d = {}
for r in out:
    if r["origin_csv"] == "numbers_cap_100hz.csv" and r["quantity"] in ("P_plus/P_bem", "P_sinc_mesh/P_bem"):
        d.setdefault(r["design"], {})[r["quantity"]] = float(r["value"])
des = sorted(k for k, v in d.items() if len(v) == 2)
a = np.array([d[k]["P_plus/P_bem"] for k in des]); b = np.array([d[k]["P_sinc_mesh/P_bem"] for k in des])
corr = float(np.corrcoef(a, b)[0, 1])
src = "numbers_cap_100hz.csv (derived)"
for q, v in [("corr(P_plus/P_bem, P_sinc/P_bem) over 6 cap designs", corr),
             ("min P_sinc/P_plus over 6 cap designs", float((b / a).min())),
             ("max P_sinc/P_plus over 6 cap designs", float((b / a).max())),
             ("min P_minus/P_bem over 6 cap designs", float((1 - a).min())),
             ("max P_minus/P_bem over 6 cap designs", float((1 - a).max())),
             ("min 2*P_sinc/P_bem over 6 cap designs", float((2 * b).min())),
             ("max 2*P_sinc/P_bem over 6 cap designs", float((2 * b).max()))]:
    out.append(dict(figure="fig2", quantity=q, design="all6", geometry="spherecap", kR=18.3183, resolved=0,
                    value=v, unit="-", source_file=src, origin_csv="derived"))
    print(f"{q}: {v:.4f}")

for p in sorted(glob.glob(os.path.join(VIZ, "fig*.png"))):
    im = np.asarray(Image.open(p).convert("RGB"))
    c = [tuple(int(x) for x in im[i, j]) for i in (0, -1) for j in (0, -1)]
    ok = all(t == (255, 255, 255) for t in c)
    print(os.path.basename(p), c, "OK" if ok else "NOT WHITE")
    out.append(dict(figure=os.path.basename(p), quantity="corner_pixels_all_white", value=int(ok), unit="bool",
                    source_file=os.path.basename(p), origin_csv="derived"))

with open(os.path.join(VIZ, "viz_numbers.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=COLS)
    w.writeheader()
    for r in out:
        w.writerow({c: r.get(c, "") for c in COLS})
print("rows", len(out))
