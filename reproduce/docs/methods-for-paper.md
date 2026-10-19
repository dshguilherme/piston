# Methods text for the paper (draft, for the authors)

Source rule: every statement below was taken from the code and logs, with the location in brackets so the authors can check it before moving text into the manuscript. Notes in `paper_jsv/supplement/novelty/` were used only where marked "(note)". **[UNCONFIRMED]** marks a statement that could not be confirmed from code or logs. Paths are relative to `harmonic_shell_course/`. The voice follows `paper_jsv/.claude/skills/paper-humanize/references/author-voice.md`; the bracketed sources and this header are for checking and should be removed in the paper.

---

## Computational pipeline

All results come from one pipeline with six steps: a shell model gives the response to a harmonic load, the response is sampled as normal velocity on the mid-surface, a radiation operator turns that velocity into a power, the power in decibels is the objective, a gradient-based optimizer changes the element thicknesses, and every design is finally evaluated under all operators. The same structural model serves every operator, so the operators differ only in the matrix that maps velocity to power.

**Shell model and response.** The cap is the mid-surface of a sphere of radius $R=10$ m, mapped from the unit square by polar angle $\theta=\theta_0 u$ and azimuth $\phi=2\pi v$ [`tests/obstacle_build_ctx.m:48-50,161-164`]. Quadratic $C^1$ B-splines on $40\times24$ elements, periodic in $v$, carry the Kirchhoff--Love displacement, which gives 3024 coefficients [`tests/r0_pilot_runs.m:211-213`; `results/r0_revision/exact_opt/cap_exact_factor_q8.log:1`] (the split $3\times(40+2)\times24$ is our arithmetic; no log prints it [UNCONFIRMED as a printed value]). The element arrays come from the authors' GeoPDEs extension (`op_KL_shells_elements`, `op_u_v_elements`) [`obstacle_build_ctx.m:76-77`], so stiffness and mass are cheap to reassemble for any thickness vector $\mathbf h$:
$\mathbf K(\mathbf h)=\sum_e (h_e^3\mathbf B_e+h_e\mathbf S_e)$, $\mathbf M(\mathbf h)=\sum_e h_e\mathbf M_e$ [`tests/r0_pilot_eval.m:15-16`].
The state equation is $[(1+i\omega\beta)\mathbf K-\omega^2\mathbf M]\mathbf q=\mathbf F$ on the free coefficients, with $\beta=0.1/\omega$ and no mass-proportional damping [`r0_pilot_eval.m:13-21`; `r0_pilot_runs.m:217`]. The rim control points are fixed in all three components (2952 free coefficients remain) [`r0_pilot_runs.m:215`; `tests/r0_fsi_eval.m:6`]; whether this leaves the rotations free in the sense of the continuum model is [UNCONFIRMED] by a run. A force of 1 N along the axis is shared equally by the 24 control points of the apex row [`obstacle_build_ctx.m:131-141`]. One sparse LU factorization per design serves the state and the adjoint solve, because the dynamic stiffness is complex symmetric [`r0_pilot_eval.m:19-37`].

**Velocity sampling.** A sparse matrix $\mathbf B$ maps $\mathbf q$ to the normal displacement at surface quadrature points and returns their surface weights [`src/r0_surface_samples.m:25-38`]; the normal velocity is $i\omega\mathbf B\mathbf q$. The sinc operator uses the $4\times4$ structural Gauss points (15 360 points); the exact operator and ERP use $8\times8$ points (61 440) [`r0_pilot_runs.m:42`; `tests/r0_opt_export_samples.m:27`].

**Radiation operators.** Each operator is stored as a factor $\mathbf F$, and the power is $P=\omega^2\|\mathbf F\mathbf q\|^2$ [`r0_pilot_eval.m:22-23`]. Three factors are used:

- *Sinc.* The code does not form the dense kernel $j_0(k|\mathbf x_q-\mathbf x_r|)$. It writes $j_0$ as an average of plane waves over directions and integrates that average with $32$ Gauss--Legendre nodes in $\cos\theta$ times 64 azimuths, so that $\mathbf L^H\mathbf L=\alpha\,\mathbf B^H\mathbf W\,j_0(k r)\,\mathbf W\mathbf B$ up to the angular quadrature, positive semidefinite by construction [`src/r0_angular_gram.m:24-44`; `src/r0_build_consistent.m:10-13`]. At the design point it agrees with the spectral (Funk--Hecke) sinc power to $10^{-7}$ over 3121 designs [`results/r0_revision/exact_opt/cap_reopt_eval.log:1`]; for other cells the agreement is $2\times10^{-7}$ to $1.5\times10^{-4}$ (note: exp6 §2.3).
- *Exact two-sided.* Per azimuthal order $|m|\le20$, a Galerkin basis of 48 functions $\sin^{|m|}\theta\, s\,P_n(2s^2-1)e^{im\phi}$, with $s^2=(\cos\theta-\cos\theta_0)/(1-\cos\theta_0)$, is whitened (Gram eigenvalues below $10^{-12}$ of the largest dropped), which leaves 35 to 48 functions per order and 1730 rows in total at $\theta_0=30^\circ$ [`tests/r0_nov_capspec.py:135-160`; `cap_exact_factor_q8.log:2-42`]. The radiation blocks are built with the spherical-harmonic series truncated at $L=4000$ and $L=8000$, combined as $\mathbf H=2\mathbf H_{8000}-\mathbf H_{4000}$, eigendecomposed, and their negative eigenvalues set to zero; then $\mathbf F_x=\boldsymbol\Lambda^{1/2}\mathbf V^H\mathbf T\mathbf G_m$ [`tests/r0_opt_exact_factor.py:32,44,57-77`]. Two forms of the block are in use. Experiment 5 (Table 2, Fig. 6) takes $\mathbf H_L=\tfrac{\rho c}{4}(\mathbf A_w^{-1}+\mathbf A_w^{-H})$ [`r0_opt_exact_factor.py:69`]. The regime study and the fluid-loading check take $\mathbf H_L=\tfrac{\rho c}{2}\mathbf A_w^{-1}\,\mathrm{Re}(\mathbf A_w)\,\mathbf A_w^{-H}$, positive semidefinite by construction, because the first form is dominated by round-off at low $kR$ [`tests/r0_pilot_factor.py:65-71`]; at the design point the two agree to $10^{-7}$ (note: exp6 §2.3; gate at `r0_pilot_factor.py:145-154`).
- *ERP.* $\mathbf F_e=\mathrm{diag}(\sqrt{\rho_0c_0w_q/2})\,\mathbf B$, the phase-blind power $\tfrac12\rho_0c_0\omega^2\sum_q w_q|\mathbf B\mathbf q|_q^2$ [`r0_pilot_runs.m:43`].

For the regime cells the medium is scaled, not the structure: the sound speed becomes $c'=\omega R/kR$ and the density is chosen so that $\rho'c'=\rho_0c_0$ [`r0_pilot_runs.m:33-34`].

**Objective.** The optimizer sees $f=100\,(A+10\log_{10}J)/p_0$ with $A=100$ and $p_0$ the value of $A+10\log_{10}J$ at the first iteration [`tests/r0_pilot_track.m:42-44`]. This is a scaled level, not a ratio of powers, so it depends on the units of $J$. In the regime study each operator is therefore normalized so that the uniform design has $J=J_{\rm ref}=5.9464\times10^{-7}$ W, and a warm start is rescaled so that $J({\rm start})=J_{\rm ref}$ [`r0_pilot_runs.m:33,48,89-94`]; the Experiment 5 runs use the unnormalized operators [`tests/r0_opt_cap_reopt.m:20-24`]. How $J_{\rm ref}$ was chosen is not stated in the code; it equals the exact power of the uniform design at the design point to five digits (`cap_reopt_summary.csv`: $5.9464\times10^{-7}$ W) [the intent is UNCONFIRMED]. The constraint keeps the mean element thickness at or below 50 mm, with all elements weighted equally, not by area [`r0_pilot_track.m:21,32,41`].

**Optimizer.** We use GCMMA [`src/gcmma_solver.m`] with at most 20 inner iterations, asymptote growth 1.3 and $\varepsilon_{raa}=\varepsilon_{raa0}=10^{-3}$, which override the library defaults of 10, 1.1 and $10^{-6}$ [`r0_pilot_track.m:35`; `gcmma_solver.m:106-118`]; there is no move limit. Thicknesses stay within 25 to 100 mm and every cold start is the uniform 50 mm [`r0_pilot_runs.m:37`; `r0_pilot_track.m:24-31`]. A run stops after 600 outer iterations (150 for the coupled re-optimization) or when, for three consecutive iterations, the mean scaled objective of the last 10 iterations differs from that of the previous 10 by less than $10^{-4}$ relative, the largest thickness change is below $10^{-4}$ of the bound range, and the constraint holds to $10^{-6}$ [`r0_pilot_track.m:58-68`]. The design returned is the iterate with the lowest $J$, not the last one [`r0_pilot_track.m:50`]; that is why X0 ends at a volume fraction of 0.472 [`exact_opt/cap_reopt_summary.csv`]. The gradient comes from the adjoint, $\mathbf K_D\mathbf z=2\omega^2\overline{\mathbf F^H\mathbf F\mathbf q}$ and $\partial J/\partial h_e=-\mathrm{Re}(\mathbf z_e^T\,\partial\mathbf K_{D,e}\,\mathbf q_e)$ [`r0_pilot_eval.m:25-37`], and was checked against central differences in a seeded random direction (`rng(517)`), with relative errors of $3.5\times10^{-7}$ to $1.6\times10^{-6}$ [`r0_pilot_runs.m:57-64`; `results/r0_revision/pilot/ops.log`]. The regime-study runs used two threads, the coupled re-optimizations four, and the Experiment 5 runs the MATLAB default [`drivers/pilot_queue.ps1`; `r0_fsi_reopt.m:15`].

**Robust objective (Experiment 10; in progress, no result yet).** The code adds a sample-average objective for thickness scatter: $h_{e,i}=h_e(1+\sigma\xi_{e,i})$ with independent standard normal $\xi$, $N$ fixed samples drawn once per run (`rng(4242)`, the first $N$ of 32 columns), and $J_{\rm rob}=\frac1N\sum_iP_i$, or mean plus one or two standard deviations [`tests/r0_robust_eval.m:5,44,71-82`; `tests/r0_robust_runs.m:31-38`]. It averages power, not decibels; the dB scaling and the GCMMA settings are those above [`tests/r0_robust_core.m:31-34`]. Unit tests pass ($\sigma=0$ reproduces the nominal objective to $10^{-15}$) [`results/r0_revision/robust/unit_c30_k18.log`], but no robust optimization run exists on disk, and the cells, $N$, $\sigma$ and mode that will be used are [UNCONFIRMED].

**Evaluation.** Every design is evaluated under ERP, sinc, the spectral sinc and the exact operator (also with $L=8000$ only and with 32 basis functions, as resolution checks) [`r0_pilot_runs.m:141-176`]. The designs of a cell are the uniform thickness, 40 seeded random thickness fields (`rng(1000+`$\theta_0$`)`, cosine-sine series with weights $1/(1+a+b)$, amplitude 5 to 25 mm, projected onto the bounds at mean 50 mm) [`tests/r0_pilot_setup.m:67-82`], every fifth iterate of each run and its returned design. Levels are $10\log_{10}$ of the power relative to the uniform design under the same operator; the gap $s$ of a design is $10\log_{10}$ of its exact power over the lowest exact power among all evaluated designs of the cell, a lower bound on the true gap [`tests/r0_pilot_stats.py:55,78,92`].

**Fluid loading.** The fluid operator is assembled per order from the same Galerkin basis, with its resistive part taken from the positive form above and its reactive part from the direct inverse, both Richardson-combined [`tests/r0_fsi_zop.py:92-113`]; it is conjugated to the structural time convention, symmetrized, and added as $i\omega s\mathbf Z$ to the dynamic stiffness, with $s=1$ for the scaled medium and $s=c'/c_0$ for air density [`tests/r0_fsi_export_z.py:22-28`; `tests/r0_fsi_check.py:73-75,124`]. Input, dissipated and radiated power balance to about $10^{-10}$ (note: exp7 §3; formulas at `r0_fsi_check.py:129-132`).

### Settings

| Item | Value | Source |
|---|---|---|
| Geometry | $R=10$ m; $\theta_0=30^\circ$ (design point), 90°, 120° | `tests/r0_pilot_runs.m:197-216` |
| Discretization | quadratic $C^1$ B-splines, $40\times24$ elements, 3024 coefficients (2952 free) | `r0_pilot_runs.m:211-213`; `cap_exact_factor_q8.log:1`; `r0_fsi_eval.m:6` |
| Material, damping | $E=200$ GPa, $\nu=0.3$, 7850 kg/m³; $\beta=0.1/\omega$, $\alpha=0$ | `r0_pilot_runs.m:212,217` |
| Load, frequency | 1 N axial at the apex row; 100 Hz | `r0_pilot_runs.m:214-215`; `obstacle_build_ctx.m:131-141` |
| Regime cells | $(\theta_0,kR)$ = (30,1), (30,3), (90,1.5), (90,7), (120,1.5), (120,7), (30,18.32); $c'=\omega R/kR$, $\rho'c'=\rho_0c_0$ | `r0_pilot_runs.m:33-34,197-208` |
| Air | $\rho_0=1.21$ kg/m³, $c_0=343$ m/s | `tests/r0_nov_capspec.py:26-27` |
| Sinc operator | $4\times4$ Gauss points; $32\times64$ directions | `r0_pilot_runs.m:42`; `src/r0_build_consistent.m:10-11` |
| Exact operator | 48 functions per order, $\lvert m\rvert\le20$, whitening $10^{-12}$ (1730 rows at 30°), $8\times8$ Gauss points, $L=4000, 8000$, eigenvalue clipping | `tests/r0_opt_exact_factor.py:32,44,70,76`; `r0_nov_capspec.py:149` |
| Objective | $100(A+10\log_{10}J)/p_0$, $A=100$; $J_{\rm ref}=5.9464\times10^{-7}$ W (regime study) | `r0_pilot_track.m:42-44`; `r0_pilot_runs.m:33` |
| Design variables | 960 thicknesses, 25-100 mm, start 50 mm, mean ≤ 50 mm (unweighted) | `r0_pilot_runs.m:37`; `r0_pilot_track.m:21,30` |
| GCMMA | 20 inner iterations, growth 1.3, $\varepsilon_{raa}=10^{-3}$; 600 outer (150 coupled) | `r0_pilot_track.m:35`; `r0_pilot_runs.m:36` |
| Stopping rule | windows of 10, tolerance $10^{-4}$ on objective and thickness change, constraint $10^{-6}$, 3 consecutive | `r0_pilot_runs.m:37-38`; `r0_pilot_track.m:58-68` |
| Random designs, seeds | 40 per cell, `rng(1000+θ0)`; FD checks `rng(517)` | `r0_pilot_setup.m:67-68`; `r0_pilot_runs.m:57` |
| Threads | 2 (regime runs), 4 (coupled re-opt.), default (Exp. 5) | `drivers/pilot_queue.ps1`; `r0_fsi_reopt.m:15` |
| Software | MATLAB R2025a, GeoPDEs 3.2.2, NURBS 1.4.3; Python 3.13 (NumPy 2.5, SciPy 1.18); Python 3.11 with bempp-cl 0.4.2 | `reproduce/env/` |

### Pseudo-code 1: operators (per cap and frequency)

```
input: structural context (B at 4x4 and 8x8 points, weights w), k = omega/c', rho' c'
# sinc (r0_build_consistent -> r0_angular_gram)
directions s_d, weights a_d  <- 32 Gauss-Legendre nodes in cos(theta) x 64 azimuths
L_sinc[d,:] <- sqrt(alpha * a_d) * exp(-i k s_d . x_q) * w_q * B4x4          # P_sinc = omega^2 ||L_sinc q||^2
# exact two-sided (r0_opt_exact_factor.py; r0_pilot_factor.py)
for m = -20..20:
    Psi_m   <- 48 cap functions of order m at the 8x8 points;  T_m <- whitening (drop Gram eig < 1e-12 max)
    G_m     <- Psi_m diag(w e^{-i m phi}) B8x8
    for L in (4000, 8000):
        A_w   <- T_m A_m(L) T_m^T                     # A_m from lambda_l = x^2 j_l'(x) h_l'(x), l <= L
        H_L   <- (rho c/4)(A_w^-1 + A_w^-H)           # Experiment 5
               or (rho c/2) A_w^-1 Re(A_w) A_w^-H     # regime study, fluid loading
    H <- 2 H_8000 - H_4000;  (V, lam) <- eigh(H);  lam <- max(lam, 0)
    Fx_m <- diag(sqrt(lam)) V^H T_m G_m
F_x <- stack(Fx_m)                                     # 1730 x 3024 at theta0 = 30 deg
# ERP
F_e <- diag(sqrt(rho0 c0 w_q / 2)) B8x8
```

### Pseudo-code 2: optimization loop (r0_pilot_track.m; robust variant r0_robust_core.m)

```
input: operator factor F (normalized so that J(uniform) = Jref; warm start rescaled so J(start) = Jref),
       h0 (uniform 50 mm or the returned design of another run), nmax (600; 150 coupled)
h <- h0;  p0 <- empty;  hold <- 0
for it = 1..nmax:
    K_D <- (1 + i omega beta) K(h) - omega^2 M(h);  LU <- lu(K_D)          # free coefficients only
    q   <- LU \ F_load
    J   <- omega^2 ||F q||^2
    # robust variant: for i = 1..N (fixed samples xi_i):  h_i <- h (1 + sigma xi_i); q_i, P_i as above;
    #                 J <- mean_i P_i  (or mean + k std);  gradient = mean of per-sample adjoints x (1 + sigma xi_i)
    z   <- LU \ (2 omega^2 conj(F^H F q));  dJ_e <- -Re(z_e^T dK_D,e q_e)
    f   <- 100 (100 + 10 log10 J) / p0      (p0 set at it = 1)
    g   <- mean(h) / 0.05 - 1
    h_new <- GCMMA(h, f, df, g, dg; inner <= 20, growth 1.3, raa 1e-3, bounds 25..100 mm)
    keep h as the returned design if J is the lowest so far
    if it >= 21 and |mean f(last 10) - mean f(previous 10)| / |mean f(previous 10)| < 1e-4
         and max|h_new - h| / 75 mm < 1e-4 and g <= 1e-6:  hold <- hold + 1  else hold <- 0
    if hold >= 3: stop ("converged")
    h <- h_new
return the lowest-J iterate, the full history of iterates
```

### Pseudo-code 3: evaluation (r0_pilot_runs.m 'eval'; r0_pilot_stats.py)

```
for each cell (theta0, kR):
    designs <- uniform, 40 random fields (rng(1000 + theta0)), every 5th iterate and the returned design of each run
    for each design h: solve the state once; P_op(h) = omega^2 ||F_op q||^2 for op in {ERP, sinc, sinc-spectral, exact,
                       exact L=8000 only, exact 32 functions}
    level_op(h)   <- 10 log10(P_op(h) / P_op(uniform))
    e_red(h)      <- level_own(h) - level_exact(h)                 # predicted minus true change
    s(h)          <- 10 log10(P_exact(h) / min over designs P_exact)   # lower bound on the true gap
    spread        <- |s(XU) - s(XS)|                                # exact runs from two starts
fluid loading: repeat with K_D + i omega s Z (s = 1 scaled medium; s = c'/c0 air) and report the same levels
```
