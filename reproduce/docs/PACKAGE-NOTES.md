# Package notes: local reproducibility package for the JSV paper

Built 2026-10-18 (night of 17-18 Oct) on branch `jsv-direction2`, from `paper_jsv/reproduce_plan/` (the audit of the same night). **No git writes, no push, no deletions, no edits outside `harmonic_shell_course/reproduce/`.** No MATLAB solve, no BEM solve, no optimization was run; MATLAB was started once only to check that `config.m` resolves the model files. The night-shift owner commits.

## 1. What was built

```
reproduce/                       73 files, 0.51 MB (plus out/, 10.6 MB, git-ignored)
  README.md                      honest status by level: re-plot (minutes), smoke (7 min), re-evaluate, full run (hours of MATLAB, 2.9 GB dataset)
  config.py, config.m            the only place paths are set: REPRO_ROOT (default ..), REPRO_OUT (default ./out), GEOPDES_ROOT
  run_smoke_tests.py/.ps1        15 tests; writes docs/smoke-test-log.md
  MANIFEST.csv + tools/make_manifest.py   691 rows: sha256, size, role, paper item, status vs audit, git-tracked, distribution
  CITATION.cff, LICENSE-NOTE.md  drafts (authors G. H. da Silva, A. Lenzi; no DOI; licence undecided)
  .gitignore                     re-includes *.tex, *.sty, *.txt, *.log inside the package (the root .gitignore drops them)
  env/                           requirements + conda files for the analysis and bempp environments, pip freezes, MATLAB and TeX notes
  figures/                       4 extractor copies + make_figures.py (data + PDF for any subset), tikz/ (17 sources + jsvfigs.sty),
                                 extras/ (F19-F22 generators, in progress), legacy_proto/ (first ratio figure)
  tables/                        make_tab_setup/opt/pilot/fsi.py + tabutil.py: compute the four tables, compare with the typed ones
  tests/                         smoke_sphere_modal.py, smoke_bem_stored.py (new); copies/ (4 original checks, output redirected)
  drivers/                       pilot_queue.ps1 + pilot_queue.txt, build_manuscript.ps1, profiling/ (all rescued from scratch)
  legacy/                        articleB_r0_thesis_damped.m (+ fdcheck), t3_longhorizon.m and its two prerequisites, README
  docs/                          design_labels.md, figure_to_command.md, methods-for-paper.md, smoke-test-log.md, this file
```

Copies carry a header with the original path, the original's SHA-256 prefix and the changed lines. Only path handling changed; numerics, formats and output names are untouched.

### Deviations from `REPO-LAYOUT.md`, and why

1. **No copies of `src/` and `tests/` into `matlab/`, `python/`.** The brief puts the package inside fastSIMP (layout option B). Copying about 60 model scripts two folders down would leave two versions to drift apart. The package calls the originals through `config.*`; only scripts whose paths had to change (figure extractors, smoke checks, rescued glue) are copied. For the slim public repository (option A) the copy step is still needed; `MANIFEST.csv` lists what goes.
2. **No `data/` copies of the CSV results.** Same reason: they are already tracked in `results/r0_revision/`. The manifest gives each file's role, paper item and distribution.
3. **A single `config.py`/`config.m` at the package root** instead of `config/repro_paths.py` and friends, as the brief asked.
4. **A nested `.gitignore`.** REPO-LAYOUT puts `pip-freeze-*.txt`, logs and TikZ sources in git, but the root `.gitignore` ignores `*.txt`, `*.log`, `*.tex`; without the nested file `requirements-*.txt`, the TikZ copies and the pilot queue list would have been dropped silently at commit time.
5. **Smoke tests S2, S3 (BEM re-solves) and S7, S8 (MATLAB gradient and unit tests) replaced** by checks on stored data, because the brief excludes BEM and MATLAB solves: BEM validation from `bem-validation.csv` against closed forms, and the sphere modal formulas recomputed in Python against the stored MATLAB CSVs. The exact-operator check `r0_opt_cap_checks.py` (no solve; uses the stored `.mat` assets) was added.
6. **The regime and FSI extractors keep text-mode output (CRLF on Windows).** REPO-LAYOUT proposed LF; that would break byte identity with the stored `.dat` files. The smoke test reports "byte-identical" on Windows and would report "identical after CRLF->LF" on Linux.
7. **No stage drivers (`stage1..6`, `run_all`).** Wrapping the originals without editing them cannot stop them from writing into `results/`; `docs/figure_to_command.md` gives the commands per stage instead (gap G07 stays open).
8. **No `check_numbers.py` over all 233 trace rows.** The four tables are checked cell by cell and Table 1 against the code; the three flagged statements are recomputed. A full number checker remains to be written (gap G09).

## 2. What the smoke tests showed

Final run `run_smoke_tests.ps1 --matlab` (analysis env, Python 3.13.12; log `docs/smoke-test-log.md`, 02:05): **12 PASS, 0 FAIL, 3 INFO**, 7.5 min in total, of which 113 s are the MATLAB start-up for the path check. The default run without `--matlab` took 4.9 min (01:57: 11 PASS, 3 INFO, 1 SKIP). Both are inside the 10-minute target on a machine shared with two MATLAB jobs.

| Test | Result |
|---|---|
| factor two | P_sinc/P_mono = 2.000000000000 at k = 0.3, 1, 3, 8; CSV and log byte-identical to the stored ones |
| sphere modal formulas | stored MATLAB values (spectrum, direct check, exterior benchmark) equal the Python closed forms to 1e-12; Wronskian form used by the factor-two script agrees |
| flat-disc extrapolation | 3.86e-4, 5.43e-4, 3.20e-5, 1.59e-4 (paper 3.9e-4, 5.4e-4, 3.2e-5, 1.6e-4); CSV byte-identical |
| fluid-operator checks | 14/14; `verify.csv` byte-identical; log identical up to timings |
| BEM from stored data | sphere references = closed form; errors fall monotonically with the mesh; finest error at kR = 1 is 1.4e-4; disc reference = dipole asymptote; screen converges monotonically, far/bem = 1 to 1e-8 |
| exact-operator checks (assets) | log identical to the stored `cap_reopt_checks.log` |
| re-plot data | 23 of 23 plotted `.dat` byte-identical |
| ratio figure | PNG byte-identical; PDF identical except `/CreationDate` |
| tables | Table 1: 29 settings found in code + 3 counts; Table 2: 42 cells (row G now computed from `cap_path_powers.csv`; the audit had it from a note); Table 3: 42 cells; Table F.1: 77 cells + 4 footnote values; all equal to the printed digits |
| TikZ compile | 12/12 audited figures compile on the regenerated data; the 5 in-progress figures also compile |
| F19-F22 (INFO) | the three generators run (mechanism 93-120 s); 17 data files equal the stored ones |

## 3. Findings that correct or extend the audit

1. **`results/t3_longhorizon.mat` does have a generator.** `tests/t3_longhorizon.m` at commit f7461da, only on branch `wt/t_233c661d` (worktree `.worktrees/t_233c661d`); output name, header text, iteration count (8210) and the saved variables match. It needs the f7461da version of `articleB_r0_opt_one.m` (deadline/checkpoint options) and `t2_panel_quadvelrad_warm_aip.mat` from `t2_panel_run.m`. All three are rescued byte-exact in `legacy/`. The brief said to document "no generator found"; the evidence says otherwise, so `legacy/README.md` documents the generator instead.
2. **The legacy design file predates its committed generator.** `articleB_r0_thesisdamp_..._partial.mat` (08:52, 12 Oct) holds only `L`, while the b40db56 script (14:38) saves more variables. Labels, iteration counts and `Jbest` agree with the `.txt` header, so the designs are the same, but b40db56 is not guaranteed to reproduce the bytes. Three of its dependencies (`articleB_r0_eval.m`, `articleB_r0_opt_one.m`, `obstacle_build_ctx.m`) changed after b40db56. `_partial` is just the checkpoint name; the file is complete (5 of 5 designs converged).
3. **Pilot queue.** 46 stage starts in `queue.log`, 38 final runs; the "42 runs" of the paper are these 38 plus the 4 anchor runs, which are the Experiment 5 runs X0/S0/SX/XS reused (`r0_pilot_runs.m:189`). Two aborted starts have no recorded `nmax`.
4. **`r0_fsi_reopt` argument order** is `(cell, zfile, startcode, nmax, tag)` (`tests/r0_fsi_reopt.m:1`), not the order written in the audit.
5. **`J_ref = 5.9464e-7 W`** equals the exact power of the uniform design at the design point (`cap_reopt_summary.csv`). The code does not say this is why it was chosen.
6. **Toolbox exception.** `prctile` in `tests/r0_robust_evalstudy.m:60` needs the Statistics and Machine Learning Toolbox (Experiment 10 only; everything behind the audited draft is base MATLAB).
7. **Experiment 10 (robust objective) has no runs yet.** Only unit tests and two Z exports exist in `results/r0_revision/robust/`; no `exp10` note exists. The methods draft describes the objective from the code and marks the study as in progress.

## 4. What remains manual or unverifiable

- **Upstream recomputation is not wrapped.** The original MATLAB and Python scripts still carry absolute paths (16 MATLAB, 12 Python) and write into `results/`. `config.m` puts the right folders on the path, so the absolute `addpath` calls only warn on another machine, but the output locations cannot be redirected without editing the originals.
- **GeoPDEs extension files are outside the repository** (`geopdes-master/topopt/Functions/`, four files; SHA-256 in `MANIFEST.csv` as `EXTERNAL:` rows). They must be vendored (GPL-3) before release; GeoPDEs itself has no commit hash (3.2.2 zip snapshot).
- **Not re-run, therefore unverified by execution:** every MATLAB solve and optimization; every BEM solve (only stored-data checks); exact-factor builds; pilot factors; FSI checks and re-optimizations; the frozen-field audit; the BEM export.
- **Stored logs that cannot be regenerated verbatim:** `paper_jsv/supplement/bem-validation.log` (its sign line contradicts `SIGN_SCREEN=-1` in the current code); `cap_flat_limit.csv` was never written (a/R = 0.05 ran out of memory).
- **Typed values not generated from data:** the four re-optimized points in `fig_fsi.tex:93,108` (recomputed and confirmed by `make_tab_fsi.py`, but still typed), the stop-status notes in `fig_opt_paths.tex:72-84`, the evidence numbers in `fig_decision.tex`, `fig_model_axes.tex`, `fig_bench_cards.tex` ("2 %" and "6 %" of Fig. 3 were not recomputed).
- **The `-List` strings passed to the pilot queue, the call lines of the FSI re-optimizations and the paths of the `iwZ` files** were never recorded; `drivers/pilot_queue.txt` is a reconstruction.
- **Environments pinned but not rebuilt.** The `.yml`/`requirements` files list the versions that ran; no fresh environment was created from them (no download was allowed). The smoke tests ran in the authors' miniconda base interpreter.
- **Memory** was never logged in the original runs.
- **In-progress work in other sessions:** Experiment 8 long runs (11 files in `longruns/` changed during the build; excluded from the package) and the F19-F22 figure sources (three TikZ files changed since the audit). The TikZ copies in `figures/tikz/` are a snapshot of 18-Oct 01:3x; re-copy them when those figures are final and re-run `tools/make_manifest.py`.

## 5. Size and distribution

- The package itself: 73 files, 0.51 MB (largest: `MANIFEST.csv`, about 0.15 MB). No `.mat`, `.npz` or PDF is inside it.
- `MANIFEST.csv` (regenerated at the end of the build) classifies 691 files: 509 unchanged since the audit, 14 changed (all in `longruns/` or the in-progress F19-F21 TikZ sources), 168 new (package files, logs not in the draft, `longruns/`, `passrule/`, `robust/`), none missing. Distribution: git 461 files, 20.5 MB (CSV, JSON, `.dat`, logs, code, the two frozen input `.mat` files of 213 kB and 379 kB); release asset / Zenodo dataset 111 files, 2.9 GB; regenerate-not-shipped 5 files, 185 MB (basis caches, BEM checkpoints); exclude 110 files, 12.6 MB (`paper_d2/`, `longruns/`, `robust/`, superseded files). Not yet listed: `results/t2_panel_quadvelrad_warm_aip.mat` (63.8 kB, input of the rescued `t3_longhorizon.m`; SHA-256 in `legacy/README.md`).
- Small `.mat`/`.pdf`/`.png` files marked `git` in the manifest are ignored by the root `.gitignore` and would need `git add -f` (or moving into the package, where the nested `.gitignore` could re-include them).

## 6. Manuscript items for the author (the manuscript was not edited)

Values recomputed by the smoke test `flagged_manuscript_numbers` from the stored CSVs:

1. **l.250**: "within 0.43 to 0.70 up to $kR=27$" should read **0.39 to 0.70**; the minimum over $10\le kR\le27$ is 0.393 (design G at $kR=23$; C and D at $kR=11$ give 0.396, 0.397) in `results/r0_revision/novelty/cap_frozen_sweep.csv`.
2. **l.530**: "$\rho$ between 0.011 and 0.048" are **per-cap minima** (150°: 0.0106; 120°: 0.0479) in `novelty/cap_modes_pencil.csv`, not a range over cases; "21 to 95 times" is likewise the worst case per cap.
3. **`figs/fig_decision.tex:33`**: "lost 2--30 dB" needs the medium: **2.2-14.5 dB in the scaled medium, 20.4-29.6 dB with the air density** (best exact design, cells with $kR\le3$, `fsi/fsi_summary.csv`; Table F.1 and the card say 2-15 dB scaled).

From the audit, still open: l.253 cites trace row 43 for a value that comes from row 54; rounding that holds but differs (1.03 % vs 1.025 %; 2.5e-15 vs 2.33e-15; "to 4e-7" holds for the best finite-difference step only).

From reading the code for the methods draft (`docs/methods-for-paper.md`), Table 1 and Appendices A/B omit or differ on:
- the objective is $100(A+10\log_{10}J)/p_0$ with $A=100$, a scaled level, not "power in dB scaled by its first value"; the regime study normalizes each operator to $J_{\rm ref}$ and rescales warm starts, Experiment 5 does not;
- GCMMA overrides (growth 1.3, $\varepsilon_{raa}=10^{-3}$) are not stated; only 1 of the 6 Experiment 5 runs met the stopping rule;
- the sinc operator is built from a $32\times64$ direction quadrature (positive Gram factor), not from the dense $j_0$ kernel written in App. A;
- "48 Galerkin functions per order" is before whitening; 35-48 remain, which is where 1730 rows come from;
- the regime study and FSI use the positive-semidefinite form of the exact block, App. B describes only the clipped form;
- the returned design is the lowest-$J$ iterate (X0 ends at volume fraction 0.472);
- the time convention differs between the structural model ($e^{i\omega t}$) and the exact operator ($e^{-i\omega t}$), handled by conjugation in the FSI code.

## 7. Before release (author decisions and steps)

1. Decide licences (`LICENSE-NOTE.md`); then add `LICENSE`, `LICENSE-DATA`, and the `license:` field of `CITATION.cff`.
2. Vendor the four GeoPDEs extension files; record the GeoPDEs 3.2.2 and NURBS 1.4.3 archive hashes.
3. Commit `reproduce/` (the nested `.gitignore` already re-includes the TikZ sources and text files); the manuscript can stay private.
4. Upload the 111 release assets with `MANIFEST.csv` to a Zenodo dataset record; mint the software DOI; fill the placeholders in `README.md`, `CITATION.cff` and the Data availability statement (`ms:507`).
5. When Experiment 10 or the F19-F22 figures enter the paper: copy their generators, extend `figures/make_figures.py` and the smoke tests, re-copy the TikZ sources, re-run `tools/make_manifest.py` and the smoke tests.
