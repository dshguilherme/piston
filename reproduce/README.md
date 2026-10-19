# When is a piston good enough? Reproducibility package

Code and data for

> G. H. da Silva and A. Lenzi, *When is a piston good enough? Limits of radiation simplifications for shell optimization and numerical benchmarks*, Journal of Sound and Vibration (submitted).

- Software DOI: `10.5281/zenodo.XXXXXXX` (placeholder, not yet minted)
- Data DOI: `10.5281/zenodo.YYYYYYY` (placeholder, not yet minted)
- Licence: code `<LICENCE PLACEHOLDER, proposed GPL-3.0-or-later>`; data and figures `<PLACEHOLDER, proposed CC-BY-4.0>`; see `LICENSE-NOTE.md`
- Status: release candidate, 2026-10-19. The folder above `reproduce/` is the repository root (`REPRO_ROOT`); see the top-level `README.md` for the layout.

## What can be reproduced, and at what cost

| Level | What it does | Time on the authors' laptop | Needs | Status |
|---|---|---|---|---|
| 1. Re-plot | Regenerate the data of every data figure and its PDF, the four tables (computed and compared cell by cell with the typed tables), from the stored CSV results | about 1 min (+ 3 min to compile all TikZ figures) | Python (analysis env), TeX Live 2025 | **works**: 23 of 23 plotted data files byte-identical, all tables match (see `docs/smoke-test-log.md`) |
| 2. Smoke test | Level 1 plus the analytic and operator checks: factor of two, sphere modal formulas, flat-disc extrapolation, fluid-operator checks, BEM validation from stored data, exact-operator checks on the stored designs | about 7 min | as level 1; the exact-operator check also needs the large `.mat` files (it skips without them) | **works** (`python run_smoke_tests.py`) |
| 3. Re-evaluate stored designs | Exact factors, regime-study evaluation, fluid-loading checks from the stored designs | tens of minutes to hours | MATLAB R2025a + GeoPDEs, Python, the 2.9 GB dataset | commands listed in `docs/figure_to_command.md`; not wrapped |
| 4. Full recomputation | Everything from the structural model, including all optimisation runs and BEM solves | about 13-15 h of MATLAB and Python, plus about 6 h of BEM | MATLAB R2025a + GeoPDEs 3.2.2 + the authors' four GeoPDEs extension files, two Python envs, 16 GB RAM (32 GB recommended) | possible with author help only: the original scripts are called unchanged and still carry absolute paths, and the GeoPDEs extension files are not yet vendored (see `docs/PACKAGE-NOTES.md`) |

The large binary results (`.mat`, `.npz`; 111 files, about 2.9 GB) are **not** in git. They are listed with SHA-256 in `MANIFEST.csv` (distribution `release asset / Zenodo dataset`) and will go to the Zenodo dataset record. Levels 1 and 2 need none of them except for one smoke test.

## Contents

```
README.md               this file
config.py, config.m     the only place paths are set (REPRO_ROOT, REPRO_OUT, GEOPDES_ROOT); no absolute paths
run_smoke_tests.py/.ps1 smoke tests; write docs/smoke-test-log.md
MANIFEST.csv            every input, result, figure source and package file with size, SHA-256, role, paper item, distribution
CITATION.cff            draft citation metadata (no DOI yet)
LICENSE-NOTE.md         proposed licences (GPL-3.0 code, CC-BY-4.0 data) for the authors to decide
env/                    requirements and conda files for the two Python environments, MATLAB and TeX notes, pip freezes
figures/                figure drivers (copies of the original extractors with paths parametrised) + tikz/ (figure sources)
  make_figures.py       one entry point: extract data and compile PDFs for any subset of figures
  extras/               in-progress figures F19-F22 (not in the audited draft)
  legacy_proto/         superseded first version of the ratio figure (rescued from a scratchpad)
tables/                 Table 1 (code-read check) and Tables 2, 3, F.1 computed from the stored CSVs and compared with the paper
tests/                  smoke-test scripts; tests/copies/ holds copies of original check scripts with output redirected
drivers/                pilot_queue.ps1 + pilot_queue.txt (regime-study run queue, reconstructed), build_manuscript.ps1, profiling/
legacy/                 rescued generators of the frozen input designs, byte-exact from git history, with README
tools/make_manifest.py  regenerates MANIFEST.csv
docs/                   design_labels.md, figure_to_command.md, methods-for-paper.md, smoke-test-log.md, PACKAGE-NOTES.md
```

Every copied script starts with a header naming its original, the original's SHA-256 prefix and the exact lines changed. Only path handling was changed; numerical code, formats and output names are as in the originals.

## Requirements

| Component | Version used for the paper | File |
|---|---|---|
| Python, analysis | 3.13.12; numpy 2.5.2, scipy 1.18.1, matplotlib 3.11.2 | `env/environment-analysis.yml`, `env/requirements-analysis.txt` |
| Python, BEM | 3.11.17; bempp-cl 0.4.2, numpy 2.4.6, scipy 1.17.1, numba 0.68.0 (CPU backend) | `env/environment-bempp.yml`, `env/requirements-bempp.txt` |
| MATLAB | R2025a (25.1.0.2943329), base MATLAB only | `env/matlab_requirements.md` |
| GeoPDEs / NURBS | 3.2.2 / 1.4.3 (pure `.m`) + four author extension files | `env/matlab_requirements.md` |
| TeX | TeX Live 2025, pdfTeX 1.40.27, pgfplots 1.18.1 | `env/tex_requirements.md` |
| Machine | Intel i5-12500H, 12 cores, 15.7 GB RAM, Windows 11; shared with other jobs, so times are upper bounds | |

The environment files pin the versions that ran; the environments were not rebuilt from these files.

## Quick start

```powershell
conda env create -f reproduce/env/environment-analysis.yml
conda activate piston-analysis
cd reproduce
python run_smoke_tests.py            # about 7 min; PASS/FAIL per test; writes docs/smoke-test-log.md
python figures/make_figures.py       # all paper figures -> out/figs/pdf
python tables/make_tab_opt.py        # Table 2 (also make_tab_setup.py, make_tab_pilot.py, make_tab_fsi.py)
```

Paths: by default `REPRO_ROOT` is the folder above `reproduce/` and every output goes to `reproduce/out/`. Set `REPRO_ROOT`, `REPRO_OUT`, `GEOPDES_ROOT`, `REPRO_PDFLATEX`, `MATLAB_EXE` to override. In MATLAB: `addpath reproduce; c = config(true);`.

## Paper item -> command

See `docs/figure_to_command.md` for every figure and table with its command, inputs, outputs and expected time, for re-plotting and for recomputation.

## Expected outputs and tolerances

| What | Tolerance |
|---|---|
| Re-plot from stored data | byte-identical `.dat` (the two regime/FSI extractors write CRLF on Windows, as the stored files; on Linux they are identical after CRLF->LF); tables equal to the printed digits |
| Analytic checks (factor two, sphere modes, flat disc, fluid operator) | as stored, byte-identical logs and CSVs under the analysis environment |
| Re-evaluation of stored designs | 1e-6 relative in power |
| BEM vs exact on identical traces | 0.1-0.6 % (fine mesh), 1.03 % (coarse) |
| Re-optimisation, design point (Table 2) | end levels within about 0.1 dB at the same thread count; the 9.5 dB path gap within about 1 dB |
| Regime study (Table 3, Fig. 7) | each cell's verdict unchanged; gaps within about 1 dB, except the 30 deg cap at kR = 3 (reference resolves its exact designs to about +/-2 dB) |
| Fluid loading (Fig. 8, Table F.1) | coupled powers within 1 %; level changes within 0.2 dB, except as above |

Gradient-based runs of 600 iterations are not bitwise portable across machines, BLAS builds or thread counts; the thread count of each run is in `env/matlab_requirements.md`.

## Inputs that are not recomputed

- **Designs U, B, C, D, E, G** (Section 5) are frozen inputs: `results/articleB_r0_thesisdamp_spherecap_40x24_f100_partial.mat` (213 kB). Its generator, `articleB_r0_thesis_damped.m`, is rescued byte-exact in `legacy/`; letters and objectives are in `docs/design_labels.md`.
- **The cylindrical-panel field** of Appendix A, `results/t3_longhorizon.mat` (379 kB), is a frozen input. Its generator `t3_longhorizon.m` and its two prerequisites were found on branch `wt/t_233c661d` (commit f7461da) and are rescued in `legacy/`; the audit had listed it as having no generator.

## Known limitations

- One Windows machine, MATLAB R2025a only. Other releases untested.
- Designs with exact efficiencies of 1e-8 and below sit at the resolution limit of the reference; small numerical differences can move them by more than the tolerances above. The verdicts are what the package guarantees.
- Memory was never logged in the original runs; the 16/32 GB figures are estimates from matrix sizes and two out-of-memory events.

## Citation and contact

See `CITATION.cff`. Please cite the article and the version DOI of the release you used.
Guilherme Henrique da Silva, LVA, Universidade Federal de Santa Catarina. dshguilherme@gmail.com
