# Design letters and run codes (DRAFT)

Status: DRAFT, compiled by reading code and notes; no solver was run. Paths are relative to `harmonic_shell_course/` (called H) unless stated. Every row cites `file:line`. Tags: **[INFERRED]** = deduced from code or notes, not stated in one place; **[UNCONFIRMED]** = could not be confirmed from the files read.

## 0. The one-paragraph answer

The six frozen cap designs of Section 5 are the uniform thickness (U, `t = 0.05`) plus the five optima in `results/articleB_r0_thesisdamp_spherecap_40x24_f100_partial.mat` (struct array `L`, 5 entries, written by `tests/articleB_r0_thesis_damped.m`). Letters were assigned in `paper_jsv/prototype/figs/make_fig_data_d2.py:21-22` (list order `DES`, dict `LAB`). The earlier-draft letter for the sinc design was R; the manuscript renamed it G because R is the cap radius (`make_fig_data_d2.py:22` comment; `manuscript_proto_voice.tex:5`; the old mapping is still live in `tests/r0_d2_paper_figures.py:40`). There is no design A and no design F in any file read **[UNCONFIRMED: why the letters skip A and F]**.

## 1. Frozen legacy designs

Struct order in the .mat is the `jobs` order of the legacy script (`legacy/articleB_r0_thesis_damped.m:8`), verified by reading the file with h5py: L(1)=dyncomp, L(2)=inputpower, L(3)=quadvel-rad, L(4)=quadvel-I, L(5)=quadvel-ones (each `tbest` is 1x960). Each stores `label, metric, R0mode, tbest, hist, meta, meas`. The mat holds only the variable `L` (see README caveat).

| Letter | Legacy label string (`L(j).label`) | Objective it minimised | Stored in (file, struct index) | Cite |
|---|---|---|---|---|
| U | (none; not optimised) uniform `t = 0.05` | not applicable | built in code: `0.05*ones(ctx.ne,1)` | `tests/r0_frozen_field_audit.m:21` (`labels={'uniform'}`) ; `tests/r0_opt_cap_reopt.m:33` |
| B | `quadvel-ones` | `metric='quadvel'`, `R0mode='ones'`: J = w^2 abs(sum_e (C u)_e)^2, net volume velocity (monopole limit, cancellation-seeking) | partial.mat, `L(5)` | legacy `:8`; `src/articleB_r0_eval.m:24-25,78-81`; letter `make_fig_data_d2.py:22`; text `manuscript_proto_voice.tex:236` |
| C | `inputpower` | `metric='inputpower'`: active input power J = -0.5 w Im(F^H u) | partial.mat, `L(2)` | legacy `:8`; `src/articleB_r0_eval.m:54-57`; `src/harmonic_kd_sens.m:24-29`; `manuscript_proto_voice.tex:236` |
| D | `dyncomp` | `metric='dyncomp'`: dynamic compliance J = abs(F.' u) | partial.mat, `L(1)` | legacy `:8`; `src/articleB_r0_eval.m:50-53`; `src/harmonic_kd_sens.m:21,72` |
| E | `quadvel-I` | `metric='quadvel'`, `R0mode='I'`: J = w^2 norm(u)^2, phase-blind sum of squared velocity coefficients | partial.mat, `L(4)` | legacy `:8`; `src/articleB_r0_eval.m:23-24,76-77`; `manuscript_proto_voice.tex:236` |
| G (was R) | `quadvel-rad` | `metric='quadvel'`, `R0mode='rad'`: J = w^2 (C u)^H R (C u), the radiated power with the legacy Nystrom sinc operator (element-centroid form) | partial.mat, `L(3)` | legacy `:8`; `src/articleB_r0_eval.m:25-26,82-86`; `manuscript_proto_voice.tex:236` |

Common setup of the five runs (legacy script): sphere cap, `R=10`, `alpha=30 deg`, `E=2e11`, `NU=0.3`, `RHO=7850`, `deg=[2 2]`, mesh `nsub=[40 24]`, `nquad=[4 4]`, `bc=2`, pole load, frequency 100 Hz (`:10`); damping `alpha=0`, `beta=0.1/omega` (`:15`); bounds `tmin=0.025, tmax=0.10`, `volfrac=0.5`, `nmax=600`, `ftol=xchgtol=1e-4`, `nhold=3`, `window=10`, `objscale='log'`, `Abar=100`, `dbfac=10` (`:10,12,18`); start: uniform `0.5*tmax = 0.05` (`src`-side default in `tests/articleB_r0_opt_one.m:30`, because `opt.tinit` is removed at legacy `:38`). All five stopped `converged` (`results/articleB_r0_thesisdamp_spherecap_40x24_f100.txt:4-8`).

The same script also contains a cylindrical-panel case (`:11`, `freq 6.53547`); its result file is not among the frozen inputs listed here and was not found in `results/` by the name pattern `articleB_r0_thesisdamp*` (only the sphere cap files exist). **[UNCONFIRMED whether the panel branch ever ran from this script]**; the panel field `t3_longhorizon.mat` comes from `legacy/t3_longhorizon.m` (see README).

Letter-to-index used by the figure code: `DES = [uniform, quadvel-ones, inputpower, dyncomp, quadvel-I, quadvel-rad]` -> `[U, B, C, D, E, G]` (`make_fig_data_d2.py:21-22`; same list in `paper_jsv/reproduce_plan/rescued_scratch/fig_ratio_proto.py:33`, which defines no letter dict). `r0_frozen_field_audit.m:26-31` builds the stored list `labels` (uniform, then the L labels in file order: dyncomp, inputpower, quadvel-rad, quadvel-I, quadvel-ones) in `res/spherecap_frozen_states.mat`; that is why `r0_opt_cap_reopt.m:49` selects the legacy G design by the string `'quadvel-rad'`.

## 2. Experiment 5 runs (exact vs sinc re-optimisation, `tests/r0_opt_cap_reopt.m`)

All: cap `theta0=30 deg`, 100 Hz, objective through the `'consistent'` interface J = w^2 norm(factor u)^2, `nmax=600`, same opt struct as legacy (`r0_opt_cap_reopt.m:24-25`); operator S = spline-consistent sinc `r0_build_consistent(ctx,1.21,343,[4 4],32,64)` (`:22`), operator X = exact two-sided factor `Fx` from `cap_exact_factor_q8.mat` (`:20-21`). Stage string = (operator letter, start letter): `switch stage(1)` operator (`:44`), `switch stage(2)` start (`:45-50`). Results saved as `results/r0_revision/exact_opt/run_<stage>.mat` (`:55`).

| Paper name | Data/run name | Objective | Start design | nmax | Notes / cite |
|---|---|---|---|---|---|
| S0 | S0 | sinc | uniform `t=0.05` (`case '0'`, `:46`) | 600 | `:11`; stopped at nmax (exp5 md:46) |
| X0 | X0 | exact | uniform | 600 | `:11`; not converged, still falling 0.148 dB/100 its (exp5 md:47; manuscript `:297`) |
| XS | XS | exact | end (`tbest`) of `run_S0.mat` (`case 'S'`, `:47`) | 600 | `:12` |
| SX | SX | sinc | end of `run_X0.mat` (`case 'X'`, `:48`) | 600 | `:12` |
| SG | SR | sinc | frozen G = `fz.designs(:,strcmp(fz.labels,'quadvel-rad'))` (`case 'R'`, `:49`) | 600 cap; converged at 109 | `:13`; exp5 md:50; manuscript `:5,310`; "SR/XR in data = SG/XG in paper; design R renamed G" (manuscript `:5`, `make_fig_data_d2.py:22`) |
| XG | XR | exact | frozen G (same as above) | 600 | `:13`; exp5 md:51; manuscript `:313` |

Order of columns of the optimisation-path figure: `for run in ("S0","SR","SX","X0","XR","XS")` (`make_fig_data_d2.py:115`). Paper Table 2 rows (S0, SG, SX, X0, XG, XS) are at `manuscript_proto_voice.tex:309-314`; its last row "G (earlier)" is the frozen legacy design itself (`:316`).

The pilot anchor cell `c30_k18` (30 deg cap at kR = 2*pi*100/343*10 = 18.3, the paper case, `r0_pilot_runs.m:205`) reuses these runs instead of re-running: XU = X0, SU = S0, SX = SX, XS = XS (`r0_pilot_runs.m:189` `map=struct('XU','X0','SU','S0','SX','SX','XS','XS')`; `r0_long_runs.m:232`; `r0_robust_runs.m:83`). SR and XR are never used as pilot codes.

## 3. Experiment 6 (regime pilot) and Experiment 7 codes

Code = two letters. **First letter = objective (rung)**: `X` exact two-sided operator (rung R4), `S` spline-consistent sinc (rung R1), `E` equivalent radiated power ERP (rung R0), phase-blind (`r0_pilot_runs.m:10-13,80-82`; exp6 md:74-80). **Second letter = start**: `U` uniform `t=0.05` (no `tinit`), `X` end design of the same cell's `XU` run, `S` end design of the same cell's `SU` run (`:85-88`). Stage dispatch: `case {'XU','SU','EU','SX','EX','XS','SUn'}` (`:78`).

| Code | Objective | Start | Cite / note |
|---|---|---|---|
| XU | exact | uniform | `:20`, `:85` (no tinit) |
| SU | sinc | uniform | `:20` |
| EU | ERP | uniform | `:20` |
| SX | sinc | end of XU | `:21`, `:86` |
| EX | ERP | end of XU | `:21`, `:86` (stage(2)=='X' always loads XU) |
| XS | exact | end of SU | `:22`, `:87` |
| SUn | sinc, normalised | uniform | `:23`: anchor-only check, "the Experiment-5 sinc runs are unnormalised"; it exists only as an option in the runfile map and evaluation list (`:145`); no SUn run file exists in `results/r0_revision/pilot/` |
| (ops, replicate, grad, eval) | not optimisation runs | | `:19,24-26` |

The 7 cells (`cell_def`, `:197-208`): `c30_k1, c30_k3, c90_k1p5, c90_k7, c120_k1p5, c120_k7` (theta0 = 30/90/120 deg, kR = 1, 3, 1.5, 7, 1.5, 7) and the anchor `c30_k18`. Every pilot run: J normalised so J(uniform)=Jref=5.9464e-7 (`:33,48`) and warm starts rescaled so J(start)=Jref (`:92-94`), nmax default 600 (`:36`).

Experiment 7 (fluid loading):
- `r0_fsi_reopt(cell, zfile, startcode, nmax, tag)` (`tests/r0_fsi_reopt.m:1`): `startcode` is `'U'` or a pilot code whose `run_<cell>_<code>.mat` `tbest` is the start (`:24-25`); output `fsi/reopt_<cell>_<tag>.mat` (`:43`).
- `dev_XU`, `phys_XU` are the **tags** of the two runs in cell `c30_k3`: start = pilot XU (`results/r0_revision/fsi/reopt_c30_k3.log`: "start XU"; `RESULT c30_k3 dev_XU`, 150 its, nmax cap), fluid-loading scale `s=1` for `dev` (device, rho' c' = rho0 c0) and `s = c'/c0 = 6.10611` for `phys` (physical density rho = 1.21), taken from the `iwZ` file made by `tests/r0_fsi_export_z.py` (`:6,21-23`; exp7 md:70-71). Tag = `<variant>_<start code>` **[INFERRED: the tag is a free string passed by the caller; the exact call lines were not found, only the log header and file names]**.
- `r0_fsi_reopt_eval.py`: designs evaluated are `uniform`, `XU`, `XS` and `reopt_<tag>` for each tag on the command line (`:37-40`); variants `vac` (s=0), `dev` (s=1), `phys` (s=sphys) (`:47`). These are the exp7 labels "vac / dev / phys" (exp7 md:69-71).

## 4. Objective-name mapping

| Name | What it is | Defined at |
|---|---|---|
| `dyncomp` | dynamic compliance J = abs(g), g = F.' u | `src/articleB_r0_eval.m:50-53`; `src/harmonic_kd_sens.m:21,72` |
| `inputpower` | active input power J = -0.5 w Im(conj(F).' u) | `src/articleB_r0_eval.m:54-57`; `src/harmonic_kd_sens.m:24-29,75,102` |
| `quadvel` (metric) | quadratic-velocity objective J = w^2 Re(u^H R0 u); needs one adjoint solve | `src/articleB_r0_eval.m:58-66`, `:12-14` |
| `quadvel-I` (`R0mode='I'`) | R0 = identity on the FE coefficients: J = w^2 norm(u)^2, phase-blind | `src/articleB_r0_eval.m:23,76-77` |
| `quadvel-ones` (`R0mode='ones'`) | J = w^2 abs(sum_e (C u)_e)^2, net volume velocity (point-source/monopole limit) | `src/articleB_r0_eval.m:24-25,78-81` |
| `quadvel-rad` (`R0mode='rad'`) | J = w^2 (C u)^H R (C u) = W_rad with the legacy Nystrom sinc kernel (`KS.'*KSSinv*KS`) and area weights | `src/articleB_r0_eval.m:26,82-86,95-98` |
| `'consistent'` mode | J = w^2 norm(factor u)^2 with a precomputed factor (sinc spline-consistent, ERP, or exact) | `src/articleB_r0_eval.m:87-89` |
| sinc (paper) vs exact | S rung = `r0_build_consistent` q4, 32x64; X rung = exact two-sided spectral factor `Fx` | `tests/r0_opt_cap_reopt.m:22,20`; exp6 md:74-83 |

`quadvel-rad` (legacy, `rad`) is the earlier Nystrom form; the paper's sinc operator for new runs is `'consistent'` (exp5 md:42: "The paper used `rad` (legacy Nystrom); here both operators use the `consistent` interface"). The manuscript says Design G "minimized the radiated power under the sinc model, in the element-centroid form" (`manuscript_proto_voice.tex:236`).

## 5. Other codes (long runs, robust runs)

`tests/r0_long_runs.m` (Experiment 8), argument string `'<kind> <cell> <name> <target>'`:
- `A <cell> <run> <target>`, run in {XU, XS, SU, SX}: continue the stored pilot optimum with the same rung (rung = `name(1)`, `:47`) to `<target>` total iterations (`:8-13`). Checkpoints `longruns/ck_A_<cell>_<run>.mat`. Present on disk: `A_c30_k3_{SU,XS,XU}`, `A_c90_k1p5_{XS,XU}`.
- `B <cell> <strategy> <budget>` (`strategy_plan`, `:120-144`): `P` plain exact from uniform; `S1` ERP (budget/8) -> sinc (7 budget/24) -> exact; `S2` sinc (budget/4) -> exact; `S3` sinc (budget/4) -> exact with move limit 0.005; `S4` kR-continuation (sinc at c30_k18, then exact; defined for c30 cells only); `S5` sinc until the sinc/exact gradient cosine < 0.9 for 3 iterations, then exact. Stage letters X/S/E as in section 3. On disk: `B_c30_k3_{P,S1,S2,S3,S4,S5}`. Exp9 says the S2/S3 sinc stages duplicated the SU path (`exp9-pass-rule.md:66`).
- `validate`, `R <...>` (refcheck): `:35-36`.

`tests/r0_robust_runs.m` (Experiment 10), argument `'<cell> <rung> <start> <N> <sigma> <mode> <target> [<fsi>]'` (`:5`): rung X or S (`:7`); start U, X (nominal exact optimum: anchor X0, else XU), S (S0 or SU), SX, XS (`:8-9,81-89`); N samples from `rng(4242); randn(ne,32)` (`:31`); mode `mean | meanstd1 | meanstd2` (`:11,33-38`); optional `fsi` = `dev | phys` (`:14-15`). Tag `<cell>_<rung><start>_N<N>_s<100*sigma>_<mode>[_fsi<v>]` (`:27-28`). No robust checkpoint existed in `results/r0_revision/robust/` at the time of reading ([UNCONFIRMED whether runs completed; directory listed only ck_*.mat names, none found]).

## 6. Cross-check against the T3 and panel files (not in the letter scheme)

`results/t3_longhorizon.mat` (`tbest`) is read as the "legacy_long_horizon" panel design (`tests/r0_frozen_field_audit.m:24`) and has no letter. It was produced by `legacy/t3_longhorizon.m` (see `legacy/README.md`).

## 7. Open points

1. Letters A and F do not occur anywhere; whether they were used in an older draft is **[UNCONFIRMED]**.
2. `tests/r0_d2_paper_figures.py:40` still maps quadvel-rad to R (old draft naming) while `paper_jsv/prototype/figs/make_fig_data_d2.py:22` uses G; any release should carry one mapping (`data/inputs/design_labels.csv` proposed in `paper_jsv/reproduce_plan/REPRODUCIBILITY-AUDIT.md` G02).
3. The 42 pilot runs of exp8 (`exp8-long-runs-multistep.md:5`) are 7 cells x 6 codes; only 38 were run in the pilot queue, the other 4 are the exp5 runs (see `drivers/pilot_queue.txt`).
