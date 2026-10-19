# MATLAB requirements

Needed only for the full recomputation (structural model, optimisation runs, BEM mesh export). The smoke tests and all figure and table drivers run without MATLAB.

## Version and products

| Item | Value | Evidence |
|---|---|---|
| MATLAB | R2025a, 25.1.0.2943329 (Student licence on the authors' machine) | `matlab -batch ver`; pilot logs print "Student License" |
| Toolboxes | **none** for every script behind a paper number | static `requiredFilesAndProducts` of the six entry scripts `r0_pilot_runs`, `r0_opt_cap_reopt`, `r0_fsi_reopt`, `r0_long_runs`, `r0_viz_sinc_freq_sweep`, `r0_bem_export_fields` (`paper_jsv/reproduce_plan/evidence/matlab_required.txt`), plus a grep of 45 `.m` files |
| Exception | `prctile` in `tests/r0_robust_evalstudy.m:60` needs the Statistics and Machine Learning Toolbox | Experiment 10 (robust objective, in progress, not in the audited draft) |
| Not used | `parfor`, `gpuArray`, `parpool`, `fmincon`/`optimoptions`, Java, compiled MEX | grep; the Parallel Computing Toolbox is installed but not called by the r0 call graph (the GPU code in `src/harmonic_kd_form.m:27` is outside it) |
| Base functions of note | sparse `lu` (4-output UMFPACK), `decomposition`, `gmres` (robust code), `eigs`, `pinv`, `legendre`, `besselj`/`bessely`, `writetable`/`cell2table`, `validateattributes` | grep |

## Third-party MATLAB code (not in the repository)

| Package | Version | Licence | Notes |
|---|---|---|---|
| GeoPDEs | 3.2.2 per `DESCRIPTION` (2021-04-07); `Contents.m` says 3.2.0 | GPL-3 | zip snapshot without `.git`, so no commit hash. Folders used: `geopdes/inst/{space,operators,msh,utils,geometry}` |
| NURBS toolbox | 1.4.3 (2021-03-29) | GPL-3+ | no compiled MEX on disk, so the pure `.m` versions ran |
| Author files inside the GeoPDEs folder | `topopt/Functions/{op_KL_shells_elements,op_u_v_elements,buildSpaces,scordelis_problem}.m` | (author code extending GeoPDEs; GPL-3 to match) | build the element stiffness and mass arrays (`tests/obstacle_build_ctx.m:76-77`). SHA-256 in `MANIFEST.csv` (rows `EXTERNAL:`). **Must be vendored before release** (audit gap G03) |

Each entry script needs about 180 files, of which about 20 are in this repository; the rest come from GeoPDEs (about 135 files) and NURBS (27 per script).

## Paths

`config.m` replaces the absolute paths: `c = config(true)` adds `GEOPDES_ROOT` (default `<REPRO_ROOT>/../../geopdes-master`), `REPRO_ROOT`, `src/` and `tests/` to the path. Checked once with `run_smoke_tests.py --only matlab_config --matlab` (no solve): `config(true)` resolves `obstacle_build_ctx`, `op_KL_shells_elements`, `op_u_v_elements`, `gcmma_solver`, `r0_build_consistent`, `sp_bspline` and `nrbmak`.

The original scripts still call `addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'))` and similar (16 scripts; `tests/obstacle_build_ctx.m:29-32`). On another machine those folders do not exist; `addpath` then warns and the paths from `config.m` do the work. The originals were not edited. They also write into `results/r0_revision/` directly, so run them against a scratch copy of the repository (`REPRO_ROOT`) if the stored results must stay untouched.

## Threads and seeds (part of the computation)

| Stage | `maxNumCompThreads` | Seeds |
|---|---|---|
| Frozen fields, Experiment 5 runs (Table 2, Fig. 6) | MATLAB default (12 cores) | `rng(517)` finite-difference direction only |
| Regime study runs (Table 3, Fig. 7) | 2 (set by the queue, `drivers/pilot_queue.ps1`) | `rng(1000+theta0)` random designs (`r0_pilot_setup.m:68`); `rng(517)` |
| FSI structural export / coupled re-optimisation | 2 / 4 (`r0_fsi_export_ctx.m`, `r0_fsi_reopt.m:15`) | `rng(517)` |
| Long runs (Exp. 8), robust (Exp. 10) | 2 | robust: `rng(4242)` optimisation samples, `rng(7777)` evaluation samples |
| Unit, seam, gradient, apex tests | default | 921 (twister), 100+d, 413, 7301 |
| Inside the model | | `src/articleB_r0_build.m:97` reseeds the global stream with `rng(12345)` |

Long GCMMA runs are not bitwise portable across thread counts or BLAS builds (sparse LU round-off about 3e-10); README section "Tolerances" gives what a re-run should match.
