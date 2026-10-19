# legacy/ - generators of the frozen input designs

Everything here is **rescued, byte-exact** from git history and is not part of the current working tree (the first two exist only on branch `wt/jsv-first-commit`, the third only on the worktree branch `wt/t_233c661d`). Paths below are relative to the repository `FastFormation/`; H = `harmonic_shell_course/`.

These scripts produced the **frozen inputs** of the study: the designs U, B, C, D, E, G of Section 5 of the paper (the paper's frozen cap fields), the panel field of Appendix A, and, through `r0_frozen_field_audit.m`, `r0_opt_cap_reopt.m` and `r0_nystrom_rank_sweep.m`, every later warm start. The two result `.mat` files are shipped as data (section 5), not regenerated.

## 1. Files

| File | Size (bytes) | git blob hash (verified equal to `git hash-object`) | Source commit | SHA-256 |
|---|---|---|---|---|
| `articleB_r0_thesis_damped.m` | 5101 | `9aa7fbd66b435d82a527d1803d16001b43359b8e` | `b40db56ef4d83b9d33a8c0557509cb13bd920ea1` (2026-10-12 14:38 -0300, "test(jsv): thesis-damping R0 ladder + T1a instrumentation"), path `harmonic_shell_course/tests/articleB_r0_thesis_damped.m` | `3cbafb9d3d3c9afd06895eee131e8412cfa2e8e1b77786135e85bc6489c3a0c2` |
| `articleB_r0_thesis_damping_fdcheck.m` | 1416 | `e715d060d5165f03347bcbf52fc8b4ffd5e8ea49` | same commit `b40db56`, same folder | `b4f66c722823a2509596922bce56bd4fae2d091e7dd07a51e34cebddcb5e3e1a` |
| `t3_longhorizon.m` | 3949 | `ed26e8d8cf839f695b309e66117db0f821e0152b` | `f7461da82a857140dc3bb9f39ebb764faa7a9145` (2026-10-14 19:52 -0300, "JSV T3-LH: add deadline/ckpt options to articleB_r0_opt_one + t3_longhorizon driver"), path `harmonic_shell_course/tests/t3_longhorizon.m` | `c4084ab58d9ddcb20759a678e163c31c8ecc6932eb60af91b9ed30701ac88992` |

**Verification of bytes.** Each file was written with `cmd /c "git -C <repo> show <commit>:<path> > <dest>"` (no PowerShell re-encoding). `git hash-object <dest>` then equals the blob hash that `git ls-tree <commit> <path>` reports, for all three files. `git log --all --oneline -- <path>` lists exactly one commit for `articleB_r0_thesis_damped.m` (b40db56), one for `articleB_r0_thesis_damping_fdcheck.m` (b40db56) and one for `t3_longhorizon.m` (f7461da).

The task asked for the first two; `t3_longhorizon.m` was added because the search of section 4 found it.

## 2. `articleB_r0_thesis_damped.m` - what it does (line numbers of this file)

- **Output files.** `matfile = RES/articleB_r0_thesisdamp_<name>_<ns1>x<ns2>_f<ftag>_partial.mat` (`:20`) with `name='spherecap'`, `ns=[40 24]`, `ftag` = `100` -> `results\articleB_r0_thesisdamp_spherecap_40x24_f100_partial.mat`. The text file is `articleB_r0_thesisdamp_<name>_<ns1>x<ns2>_f<ftag>.txt` (`:50`), i.e. **without** `_partial`. Both are written in `results\` (`RES`, `:7`). The `.mat` is saved after every optimisation job (`:31,45`) and once more after the loop (`:49`), always under the same `_partial` name.
- **Cases.** `spherecap` (R = 10, alpha = 30 deg, E = 2e11, NU = 0.3, RHO = 7850, t0 = 0.1, deg [2 2], reg [1 1], nsub [40 24], nquad [4 4], bc 2, pole load, 100 Hz, zeta 0.05; `:10`) and `cylpanel` (R = 25, L = 50, th0 = 40 deg, 6.53547 Hz; `:11`). Only the sphere-cap result file exists in `results\`.
- **Optimised designs** (the `jobs` table, `:8`; labels = `L(j).label`): 1 `dyncomp`, 2 `inputpower`, 3 `quadvel-rad`, 4 `quadvel-I`, 5 `quadvel-ones` (metric/R0mode pairs: dyncomp/-, inputpower/-, quadvel/rad, quadvel/I, quadvel/ones). Objective formulas: see `docs/design_labels.md`. Letters: dyncomp = D, inputpower = C, quadvel-rad = G (earlier R), quadvel-I = E, quadvel-ones = B; U is the uniform start.
- **Settings.** mesh 40x24; frequency 100 Hz (`:10`); damping alpha = 0, beta = 0.1/omega, i.e. 5 % equivalent modal damping at the excitation frequency, thesis Eq. 6.118 (`:2-4,15`); bounds tmin = 0.025, tmax = 0.10, volume fraction <= 0.5, nmax = 600, ftol = xchgtol = 1e-4, nhold = 3, window = 10, feastol = 1e-6, `objscale 'log'`, Abar = 100, dbfac = 10 (`:12,18`); acoustic constants rho = 1.21, c = 343, 128 (rank argument of `articleB_r0_build`, `:16`). Start: uniform `0.5*tmax = 0.05` (default of `opt_one`; `opt.tinit` is removed for a fresh job, `:38`). Resume logic (`:21-25,27-40`): saved designs that already `converged` are skipped.
- **Called functions.** `obstacle_build_ctx`, `articleB_r0_build`, `articleB_r0_opt_one`, `articleB_r0_eval` (through `measure`, `:59`), local `measure`; `geopdes` via `addpath(genpath(.../geopdes-master))` (`:6`); `addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation')` (`:6`) adds the repo root (purpose unclear **[UNCONFIRMED]**). Transitively `articleB_build_ctx`, `gcmma_solver`, geopdes and nurbs.
- **Do they exist on the current working tree (H\src, H\tests)?** All four exist: `tests\obstacle_build_ctx.m`, `src\articleB_r0_build.m`, `tests\articleB_r0_opt_one.m`, `src\articleB_r0_eval.m` (plus `src\gcmma_solver.m`, `tests\articleB_build_ctx.m`). None is missing. They also exist at b40db56 (`git ls-tree -r b40db56 --name-only`). **But three are not the same bytes as at b40db56:** `articleB_r0_build.m`, `gcmma_solver.m`, `articleB_build_ctx.m` are identical (same blob); `articleB_r0_eval.m` (b40db56 `9fa2b267...`, HEAD `f55ca38f...`), `articleB_r0_opt_one.m` (`8ea6f430...` vs `9de9b02a...`) and `obstacle_build_ctx.m` (`d2e646e0...` vs `882a87e1...`) were modified by later commits (HEAD diff: eval +/-18 lines, opt_one +/-79 lines, obstacle_build_ctx a comment and the `cylshell` periodicity, which the sphere cap does not use). For an exact rerun, take those three from `b40db56` (`git show b40db56:harmonic_shell_course/<path>`); results with the HEAD versions are not guaranteed bit-identical **[UNCONFIRMED]**.

### Comparison with the header of the stored result

`results\articleB_r0_thesisdamp_spherecap_40x24_f100.txt` (992 bytes, written 2026-10-12 08:52:08, the same second as the .mat) has the header produced by `:51-52`:

```
R0 ladder using thesis Eq. 6.118 stiffness-proportional damping
case=spherecap mesh=40x24 frequency=100 Hz alpha=0 beta=0.1/omega=0.000159154943092 s (zeta=0.05 at f) t=[0.0250,0.1000] volfrac<=0.50 nmax=600
design outer_iter stopReason evals Jbest Wrad WI Wones rho_coh CV vol offdiag_over_total o_over_d wall_s
```

Format strings, column names, `beta = 0.1/(2 pi 100) = 1.59155e-4`, bounds, volume fraction and `nmax` all match the script. Five rows follow in the order dyncomp, inputpower, quadvel-rad, quadvel-I, quadvel-ones, all `converged` (iterations 90, 475, 392, 347, 44; wall times 59.85, 356.09, 614.38, 647.63, 164.88 s, sum 1842.8 s, about 31 minutes in total). Reading the `.mat` with h5py: 5 entries in `L`, same labels, same order, same iteration counts and `Jbest` as the text, `tbest` of size 1x960.

### Why is the file called `_partial`?

Because the script hard-codes that name for the file it saves (`:20`): it is the **checkpoint/resume file** (the script loads it back at `:21-25` and saves after every job at `:31,45`). The name is never changed after the ladder completes, and the final save (`:49`) overwrites the same `_partial.mat`. The `.mat` is therefore complete (5 of 5 designs, all `converged`), despite the name. The `.txt` summary has no `_partial` suffix (`:50`). Two differences to note:

1. The on-disk `.mat` contains only the variable `L` (plus `#refs#`), while the committed script saves `L, ns, f, alpha, beta, volfrac, opt` (`:31,45,49`). The stored file also has `meta.total_iter` and `meta.tailspread10` fields. So the file was written by a slightly earlier or different version of the script than the b40db56 text **[INFERRED]**: the file's timestamp (08:52) precedes the commit (14:38). The designs and numbers agree with the text header, so this does not change the frozen inputs, but it means b40db56 is not guaranteed to reproduce the file's exact bytes.
2. `quadvel-ones` stops at volume fraction 0.4506 (< 0.5, objective 1e-18, 44 iterations), a degenerate cancellation-seeking optimum.

### Can it be re-run today?

Not as a standalone package. It needs MATLAB (R2025a was used for the static scan; the logs print "Student License", `paper_jsv/reproduce_plan/evidence/matlab_required.txt` shows base MATLAB only, no toolbox) and GeoPDEs (`geopdes-master`, with `nurbs-1.4.3` and `topopt\Functions` on the path via `genpath`), at hard-coded absolute paths (`:5-6`). GeoPDEs is unpinned (3.2.2 vs 3.2.0 zip, per the audit). The logged runtime is about 31 minutes for the sphere cap on the original machine (the panel case follows if run). A rerun writes into `results\` and would overwrite the shipped file, so run it into a scratch copy. The `fdcheck` script is a 12x8-mesh panel adjoint-vs-finite-difference check for dyncomp, inputpower and quadvel (rad, I), tolerance 5e-4 (`articleB_r0_thesis_damping_fdcheck.m:12`); it needs the same environment.

## 3. `t3_longhorizon.m`

`function t3_longhorizon(deadline, nmax, ckpt_every, outbase)` (`:1`) continues the V2 panel design (`tvb` in `results\t2_panel_quadvelrad_warm_aip.mat`, loaded at `:26-27`) with `quadvel`/`rad` on the cylindrical panel (`cylpanel`, R = 25, L = 50, th0 = 40 deg, 6.53547 Hz, 40x24, alpha = 0, beta = 0.1/omega; `:16-20`), same optimiser settings as the ladder (`:23-25`) but `nmax = 100000`, `jfloor = 0` (`:28`), deadline and checkpoint options (`:29-30`). It writes `results\<outbase>.mat` (variables `tbest, hist, meta, M`, `:38`) and `results\<outbase>.txt` (`:39-59`, text starts "T3-LH quadvel-rad warm continuation from t2 V2 endpoint"), default `outbase = 't3_longhorizon'` (`:9`).

- **It matches the stored header** `results\t3_longhorizon.txt`: same first line, then `stopReason: converged`, `iterations: 8210`, `walltime_s: 14218.0` (3.95 h), `Jbest: 4.397e-15`, `J(1)`, `Wrad`, `o/d`, "Jbest every 500 its" and "slope log10(Jbest)..." - the exact fields written at `:41-55`. File times: `.txt` 2026-10-15 04:09, i.e. a run of about 4 h started about 00:10; the smoke runs `t3_longhorizon_smoke1/2` (14 Oct 19:50, 3 iterations) used the same driver.
- **Dependencies.** Needs the modified `articleB_r0_opt_one.m` from the same commit f7461da (adds `opt.deadline`, `opt.ckpt_file`, `opt.ckpt_every`; diff of 8 added lines). The current working tree's `tests\articleB_r0_opt_one.m` does **not** contain `deadline` or `ckpt` (grep, no match), so the driver does not run as-is against H\tests. It also needs `t2_panel_quadvelrad_warm_aip.mat` (63808 bytes, SHA-256 `c8c71efbf76b0e8b9d2e48a80d4c90c9275b840dad0e4b51b8c34fc86ceea05c`) whose own generator is `tests/t2_panel_run.m` (line 65 writes it) on the same worktree branch **[not rescued here]**.
- **Branch.** f7461da is contained only in branch `wt/t_233c661d` (`git branch -a --contains f7461da`), checked out as the worktree `.worktrees\t_233c661d\` where `tests\t3_longhorizon.m` exists on disk. It is not on `jsv-direction2`, `main` or any remote branch.

### 3a. Dependencies of `t3_longhorizon.m`, also rescued

Added after the search above, so that the driver and its two prerequisites sit together. Both were written with `cmd /c "git show f7461da:<path> > <dest>"`; `git hash-object` equals the blob in `git ls-tree f7461da <path>`.

| File | git blob | Role |
|---|---|---|
| `f7461da_articleB_r0_opt_one.m` | `8fb1fce5d16b8d2ee9d86cfc7b8c7f1ced39687c` | the optimiser version with `opt.deadline`, `opt.ckpt_file`, `opt.ckpt_every` that `t3_longhorizon.m` calls. To re-run T3, put this file (renamed `articleB_r0_opt_one.m`) ahead of `tests/` on the MATLAB path; do not replace the working-tree file |
| `f7461da_t2_panel_run.m` | `9b42470fbb246f4825add3b5c8f955f69e42c47e` | writer of `results/t2_panel_quadvelrad_warm_aip.mat` (its line 65), the warm start that `t3_longhorizon.m` continues |

The chain behind the panel field of Appendix A is therefore `t2_panel_run.m` -> `t2_panel_quadvelrad_warm_aip.mat` -> `t3_longhorizon.m` (with the f7461da optimiser, about 4 h) -> `t3_longhorizon.mat`. None of these was re-run.

## 4. Search for the generator of `results\t3_longhorizon.mat`

Original premise (audit `paper_jsv/reproduce_plan/REPRODUCIBILITY-AUDIT.md`, G02; `MANIFEST-draft.csv:521`): "NO GENERATOR FOUND". The following read-only searches were run (repo `FastFormation`):

1. `git log --all --oneline -S "t3_longhorizon"` -> `4cc4ac3`, `4a0e56b`, `3ac4a0a`, `4223123`, `f7461da`.
2. `git log --all --oneline -S "T3-LH"` -> `4a0e56b`, `f7461da`.
3. `git log --all --oneline -S "longhorizon"` -> the same five commits.
4. `git grep -n -e "t3_longhorizon" -e "T3-LH" <ref>` over `HEAD`, `main`, `wt/jsv-first-commit`, `wt/jsv-t4-results`, `wt/t_233c661d`, `origin/main`, `origin/wt/jsv-first-commit`, `origin/wt/jsv-t4-results`. Hits that are **writers**: only `wt/t_233c661d:harmonic_shell_course/tests/t3_longhorizon.m` (lines 1, 9, 40, 60). Hits that are **readers**: `tests/r0_consistent_panel_reopt.m:20`, `tests/r0_frozen_field_audit.m:23`, `tests/r0_nystrom_rank_sweep.m:25` (found by `git grep` on `origin/main`, `wt/jsv-t4-results` and `origin/wt/jsv-t4-results`; no hit on the local `main` or `wt/jsv-first-commit`; absent from the `HEAD` grep, although the three files exist on disk in H\tests). Documents only: `results/t4_synthesis.txt`, `paper_jsv/reproduce_plan/*`, `NIGHT_LOG.md`.
5. `git log --all --oneline -- harmonic_shell_course/tests/t3_longhorizon.m` -> `f7461da` only. `git ls-tree f7461da harmonic_shell_course/tests/t3_longhorizon.m` -> blob `ed26e8d8...`.
6. Filesystem: `Test-Path .worktrees\t_233c661d\harmonic_shell_course\tests\t3_longhorizon.m` is True; `H\tests\t3_longhorizon.m` is False. The earlier audit only searched the checked-out branch, so it missed this.

**Conclusion.** A generator **does exist**: `tests/t3_longhorizon.m` at commit `f7461da` on branch `wt/t_233c661d`, rescued here as `t3_longhorizon.m`. It is the writer of `t3_longhorizon.mat/.txt` on the evidence that (a) the output name defaults to `t3_longhorizon`, (b) the header text is the same, (c) iteration count, wall time and slope fields match the format, (d) times are consistent (commit 14 Oct 19:52, smoke runs 19:50, final 15 Oct 04:09). (e) the saved mat (read with h5py) holds exactly the variables `tbest, hist, meta, M` that `:38` saves, and `meta` carries `deadline, source, f, ns` - the four fields this driver adds at `:35`. This is strong evidence, short of a recorded run log. Caveats: the driver loads another frozen input (`t2_panel_quadvelrad_warm_aip.mat`) and depends on a version of `articleB_r0_opt_one.m` that is not on the current branch. A run of about 4 h cannot be expected to be bit-reproducible after the code change; the file stays a frozen input.

## 5. Frozen inputs to be shipped with SHA-256

Computed with `Get-FileHash -Algorithm SHA256` on the files in `harmonic_shell_course\results\`. Both `.mat` and `.txt` are excluded from git by `.gitignore` (`*.mat` line 10, `*.txt` line 13), so they exist only in the OneDrive copy.

| File | Bytes | SHA-256 |
|---|---|---|
| `articleB_r0_thesisdamp_spherecap_40x24_f100_partial.mat` | 213240 | `b422971989562714229978c94bae24eb80e53907baf4cff95766543d75a71a63` |
| `articleB_r0_thesisdamp_spherecap_40x24_f100.txt` (header/summary of the above) | 992 | `c2ad27d74b3f51ea5ae69232436e490df32d2d99b0877242f3e2c57d9c023019` |
| `t3_longhorizon.mat` | 379307 | `d793ca0b1388865466e7f039b0a2cf9d9b3f5d22107d59491a2dcca71838a844` |
| `t3_longhorizon.txt` (header/summary of the above) | 899 | `fb104505984ae81a557e25ed13d165e961ebd9d12cc2977cedac09b112bb5b6f` |
| `t2_panel_quadvelrad_warm_aip.mat` (input of `t3_longhorizon.m`; not requested, recorded for completeness) | 63808 | `c8c71efbf76b0e8b9d2e48a80d4c90c9275b840dad0e4b51b8c34fc86ceea05c` |

The `t3_longhorizon.mat` and `.txt` hashes equal those in `paper_jsv/reproduce_plan/MANIFEST-draft.csv:521-522`.

Treat the two `.mat` as frozen inputs: ship them as data (they are small, 213 kB and 379 kB), verify the SHA-256 on download, and read them with the variables `L` (sphere cap; use `L(j).tbest`, `L(j).label`) and `tbest` (panel). The scripts here document how they were made; they are not required to reproduce the paper's figures.
