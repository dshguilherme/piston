# Visual diagnostics of the two-sided BEM comparison (cap and panel)

All numbers below are in `viz_numbers.csv` (column `origin_csv` names the per-figure file; `source_file` names the raw data). Figures are PNG (200 dpi) + PDF, white backgrounds (corner pixels checked 255,255,255 for all six PNGs; rows `corner_pixels_all_white` in the CSV).

Conventions: cap = sphere segment R = 10 m, half-angle 30 deg, design frequency 100 Hz (kR = 18.32). "+" face = side the shell normal points to = outer/convex face of the cap; "-" face = inner/concave face. mu = p+ - p- is the BEM pressure jump. BEM = bempp-cl two-sided open screen (mu = 0 at the free edge), P1, mesh m=2 unless stated. The m=2 cap mesh has only 5.24 elements per wavelength at 100 Hz (below the usual 6); the independent exact spectral solver (`results/r0_revision/novelty/`) agrees with the m=2 BEM far-field power to 1.1 % at 100 Hz (quadvel-rad 0.9886, uniform 0.9907), so this does not change any conclusion below.

## Fig. 1 - `fig1_field_maps.png/.pdf` (m=2) and `fig1b_field_maps_m3.png/.pdf` (m=3)

Caption. Normal-velocity fields imposed on the BEM, drawn on the exported BEM triangle meshes (component-aware spline trace from `tests/r0_bem_export_fields.m`, used unchanged). Rows: cap quadvel-rad, cap uniform (top view x-y), cylindrical panel uniform, panel long-horizon (unrolled: arc length s = 25 atan2(x,y) vs axial z). Columns: |v_n| (log scale, four decades below the maximum; the clamped edge is ~0), arg v_n as stored in the exported mesh (flat per triangle, cyclic colour map), and a zoom (red box) with every triangle edge drawn, centred on the steepest nodal line. Fig. 1b repeats |v_n| and the zoom on the m=3 meshes.

Key numbers: max |v_n| / area-RMS |v_n| (m/s): cap quadvel-rad 3.25e-6 / 1.27e-6; cap uniform 3.28e-5 / 6.81e-6; panel uniform 2.89e-6 / 9.25e-7; panel long-horizon 1.70e-6 / 5.15e-7. Elements per wavelength (lambda/h_max): cap 5.24 (m=2), 7.85 (m=3); panel 46.5 (m=2), 69.7 (m=3). Max and RMS change by 0.1-1.7 % from m=2 to m=3.

What it shows: the cap fields are axisymmetric ring patterns with antiphase rings separated by thin nodal circles. The nodal circles are crossed by only about 3-6 triangles at m=2, so the acoustic mesh, not the structural field, limits resolution. The panel fields are axial bands with alternating phase. The long-horizon field has slanted nodal lines and many antiphase lobes, which is the cancelling structure that makes its radiated power tiny.

## Fig. 2 - `fig2_cap_face_pressure.png/.pdf`

Caption. Surface pressure on each face of the cap at 100 Hz (BEM m=2), designs quadvel-rad (top) and uniform (bottom). The faces are separated as p+ = (K + 1/2) mu and p- = (K - 1/2) mu, where K is the double-layer direct value (Galerkin weak form, L2-projected to vertices). The columns are |p+| (outer face), |p-| (inner face), |mu|, and the local per-face power density I+ = Re(p+ conj v_n)/2 and I- = -Re(p- conj v_n)/2 (W/m^2, linear, symmetric limits). Text boxes give the face shares and the sinc ratios.

Checks:
- The jump sign was confirmed by evaluating the double-layer potential at x +/- delta n, with delta = h/4. The fitted jump/mu ratio was 0.998 (quadvel-rad) and 1.009 (uniform), against a swapped-sign ratio near -1.
- P+ + P- = P_BEM to 2e-15.
- P_BEM reproduces the committed run (`bem_shell_runs/spherecap_m2_q0.json`) to 3e-15.

Key numbers at 100 Hz, six cap designs:

| design | P+/P_BEM (outer) | P-/P_BEM (inner) | P_sinc/P_BEM | 2 P_sinc/P_BEM | P_sinc/P+ |
|---|---|---|---|---|---|
| uniform | 0.685 | 0.315 | 0.633 | 1.265 | 0.923 |
| quadvel-rad | 0.623 | 0.377 | 0.611 | 1.222 | 0.980 |
| quadvel-I | 0.502 | 0.498 | 0.559 | 1.118 | 1.114 |
| inputpower | 0.498 | 0.502 | 0.537 | 1.074 | 1.080 |
| dyncomp | 0.485 | 0.515 | 0.525 | 1.049 | 1.082 |
| quadvel-ones | 0.427 | 0.573 | 0.489 | 0.978 | 1.144 |

Correlation of P_sinc/P_BEM with P+/P_BEM over the six designs: r = 0.980.

Here P_sinc is the mesh sinc on the same m=2 mesh.

**Verdict on the hypothesis "the factor ~2 is both faces radiating versus one sheet": SUPPORTED at the 100 Hz design point, with two qualifications.**

Evidence for it:
- Both faces radiate a large share at 100 Hz. The inner face carries 31-57 % of the BEM power.
- The single-sheet sinc estimate accounts for 49-63 % of the BEM power.
- That fraction matches the outer-face share to within -8 % / +14 % (P_sinc/P+ = 0.92-1.14), and it tracks the face share from design to design (r = 0.98).

Qualification 1: it is not literally "two equal faces x one sheet". The faces are unequal for two designs (uniform 69/31, quadvel-rad 62/38), and doubling the sinc power overshoots by up to 27 % (2 P_sinc/P_BEM = 0.98-1.27).

Qualification 2: the one-face correspondence holds only near the design point. In the frequency sweep (Fig. 4, bottom panel):
- kR = 1: the inner face absorbs power (P-/P = -4.7, P+/P = +5.7 for quadvel-rad). The two-sided screen cancels, and sinc over-predicts 15x.
- kR = 4: the shares are about 40/60 but sinc/BEM is only 0.25.
- kR = 10: quadvel-rad has P+/P = 0.38 and P_sinc/P = 0.475, so sinc/P+ = 1.26.

The face-sharing picture explains the ratio ~0.5 at high kR in aggregate. It is not an identity between the sinc model and one particular face.

Not independently verified: the face split itself. It rests on the m=2 double-layer weak form and was not cross-checked with the spectral solver or at m=3.

## Fig. 3 - `fig3_cap_directivity.png/.pdf` (quadvel-rad) and `fig3b_cap_directivity_uniform.png/.pdf`

Caption. Far-field directivity at 100 Hz, BEM (two-sided, m=2) against the sinc model.
- Panels (a)-(b) are great-circle cuts phi = 0/180 and 90/270, plotted as 10 log10 of dP/dOmega relative to the BEM maximum (40 dB range). +z is the outer (convex) axis and -z the inner (concave) axis.
- Panels (c)-(d) are theta-phi maps on a common colour scale. Panel (e) is 10 log10(sinc/BEM), masked where both are 40 dB down.
- The sinc pattern is the free-space monopole-sheet far field scaled by sqrt(2). This is the unique monopole-sheet pattern whose integral equals the sinc power, because the sinc prefactor is exactly twice the free-space monopole power (breathing-sphere limit 2.013 at kR = 0.1).

Integral checks:
- quadvel-rad: the BEM pattern integral is 1.16067e-8 W against surface P_BEM 1.16063e-8 W (rel. 3.1e-5). The sinc pattern integral is 7.09228e-9 W against the mesh sinc power 7.09228e-9 W (rel. 3e-14).
- uniform: BEM rel. 1.6e-5, sinc rel. 1.3e-14.

Key numbers:
- quadvel-rad hemisphere split (theta < 90 | theta > 90): BEM 0.506 | 0.494, sinc 0.521 | 0.479.
- uniform hemisphere split: BEM 0.502 | 0.498, sinc 0.666 | 0.334.
- quadvel-rad peak dP/dOmega: BEM 1.022e-8 W/sr, sinc 7.18e-9 W/sr.

What it shows: the deficit is not a missing hemisphere. Both models radiate into both half-spaces. The BEM lobes are higher by a few dB over most directions, and both patterns are almost equal on the -z axis. The sinc pattern has deep nulls (uniform at theta ~ 25 and 125 deg) that the dipole-layer BEM pattern does not have.

## Fig. 4 - `fig4_cap_power_vs_frequency.png/.pdf`

Caption. Radiated power of the frozen cap fields (normal velocity held at the 100 Hz field, only k varied) against kR, for quadvel-rad (blue) and uniform (brown).
- Top panel: sinc spline-consistent and centroid (MATLAB, `cap_sinc_freq_sweep.csv`), sinc on the BEM mesh, the BEM m=2 far-field power (filled = resolved, open = unresolved), and the BEM surface power (faint, not used). Red squares mark points where surface and far-field BEM power differ by more than 2 %. The black/grey dash-dot curves are the independent exact spectral two-sided solution (`novelty/cap_frozen_sweep.csv`, read only).
- Middle panel: sinc / two-sided ratios, with the breathing-sphere sinc/exterior benchmark (analytic and `spherical-exterior-benchmark.csv`) for reference.
- Bottom panel: per-face shares at kR = 1, 4, 10, 16 (symlog axis).
- Grey band: unresolved, kR > 16.01, below 6 elements per wavelength on m=2 (h_max = 0.654 m). kR = 18.32 (100 Hz) and 21 are unresolved.

Which BEM column is plotted: the surface-pairing power is wrong at low kR. far/surface is 0.08-0.40 at kR = 0.1 and 0.72-0.95 at kR = 0.3 for all six designs. A finer far-field grid (192 x 384) gives identical values, and the BEM far-field power agrees with the spectral solution to within 0.3 % there. The far-field column is therefore plotted. From kR = 1 upward, far/surface is within 1 % (kR = 1: 0.99) and within 1e-4 from kR = 4.

BEM far-field power against the spectral solution: 0.97-1.007 over all points. The worst values are at the unresolved kR = 21 (0.971 and 0.970) and at kR = 8-10 for quadvel-rad (0.984 and 0.980).

Consistent sinc / BEM far-field power:

| kR | 0.1 | 0.3 | 1 | 2 | 4 | 6 | 8 | 10 | 12 | 14 | 16 | 18.3* | 21* |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| quadvel-rad | 1859 | 203 | 15.1 | 2.11 | 0.249 | 0.758 | 1.03 | 0.475 | 0.463 | 0.522 | 0.531 | 0.611 | 0.592 |
| uniform | 558 | 61.7 | 5.28 | 1.09 | 0.221 | 0.483 | 0.706 | 0.492 | 0.559 | 0.697 | 0.532 | 0.633 | 0.601 |

\* unresolved. The ratios in the table are mesh-sinc / far-field. The plotted spline-consistent ratios differ by less than 1.5 %; the exact values are in `numbers_fig4.csv`.

What it shows:
- Low kR (kR < 2): the two-sided screen is a dipole layer whose two faces cancel, so its power rises much faster with k than the sinc power. The sinc model over-predicts by orders of magnitude, up to 1859x at kR = 0.1.
- kR ~ 2: the curves cross.
- kR ~ 4: sinc under-predicts by 4x.
- kR >~ 10: the ratio settles around 0.46-0.70, i.e. near the "one sheet against two faces" value of 0.5.
- The factor ~2 at 100 Hz is therefore a high-frequency plateau value, not a constant of the model.
- The breathing-sphere curve shows the opposite limit for a closed body: sinc gives 2x the exterior power at low kR, because the prefactor carries a baffle-image factor of 2.

## Files

- Scripts (tests/):
  - `r0_viz_fig1_field_maps.py`
  - `r0_viz_cap_bem_100hz.py`
  - `r0_viz_fig2_face_pressure.py`
  - `r0_viz_fig3_directivity.py`
  - `r0_viz_cap_freq_sweep.py`
  - `r0_viz_sinc_freq_sweep.m`
  - `r0_viz_fig4_power_vs_frequency.py`
  - `r0_viz_merge_numbers.py`
- Data (viz/):
  - `cap_m2_100hz_bem.npz`
  - `numbers_cap_100hz.csv`
  - `numbers_fig1.csv`
  - `numbers_fig4.csv`
  - `cap_m2_freq_sweep.csv`
  - `cap_sinc_freq_sweep.csv`
  - `viz_numbers.csv` (merged)
  - logs and checkpoint `.npy` files of the restarted runs

The BEM sweep was killed twice by external memory pressure (other sessions' jobs; no traceback) and resumed from its CSV. The kR = 18.3183 sweep row reproduces the committed run.
