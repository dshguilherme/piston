# piston: reproduction package

Code, stored results and figure sources for

> G. H. da Silva and A. Lenzi, *When is a piston good enough? Limits of radiation simplifications for shell optimization and numerical benchmarks*, submitted to the Journal of Sound and Vibration.

The detailed guide is [`reproduce/README.md`](reproduce/README.md): what can be reproduced at which cost, environments, tolerances and known limits. Every figure and table is mapped to its command in [`reproduce/docs/figure_to_command.md`](reproduce/docs/figure_to_command.md) and in the table below.

## Layout

The repository root plays the role of `REPRO_ROOT` in `reproduce/config.py` and `reproduce/config.m`.

```
reproduce/       package: config, smoke tests, figure and table drivers, environments, docs, MANIFEST.csv
src/             MATLAB model code (isogeometric Kirchhoff-Love shell, radiation operators, GCMMA)
tests/           MATLAB and Python scripts that produced every stored result (r0_* for this paper)
results/         stored text results: CSV tables, run logs up to 2 MB, small frozen input designs
paper_jsv/prototype/figs/   TikZ sources and plotted data files of the paper figures
vendor/geopdes_topopt_functions/   the authors' four GeoPDEs extension files (element arrays of the shell)
LARGE_FILES.csv  binary results and large logs that are not in git (about 3.2 GB), with SHA-256
```

## Quick start

```bash
conda env create -f reproduce/env/environment-analysis.yml
```

```bash
conda activate piston-analysis
```

```bash
python reproduce/run_smoke_tests.py
```

```bash
python reproduce/figures/make_figures.py --paper
```

The smoke tests take about 2 minutes (15 tests; on the authors' machine 10 pass, 2 are skipped without the binary dataset and MATLAB, 3 report information) and write `reproduce/docs/smoke-test-log.md`. `make_figures.py --paper` builds the 19 figures of the manuscript into `reproduce/out/figs/pdf` in about 2 minutes; re-plotting needs no binary file. The exact-operator smoke test, any re-evaluation and any re-optimization need the binary dataset listed in `LARGE_FILES.csv`.

## Results of the later sections

| Paper item | Command | Notes |
|---|---|---|
| Guaranteed level: bar figure, pass-rule tables and scatter figure | `python tests/r0_passrule_stats.py`, then `python tests/r0_fig_passrule_v2.py` | reads `results/r0_revision/pilot`, `exact_opt`, `longruns`; writes `results/r0_revision/passrule` and `paper_jsv/prototype/figs/data` |
| Long runs and staged strategies: tables and figure | `python tests/r0_long_summary.py` (needs pandas), then `python tests/r0_fig_longruns.py` | runs: `tests/r0_long_runs.m`, about 13 h with 2 threads |
| Robust design: table of expected levels | `python tests/r0_robust_summary.py <evalstudy CSV>` | runs: `tests/r0_robust_runs.m`, `r0_robust_eval.m`, `r0_robust_evalstudy.m` (2 threads; `prctile` needs the Statistics and Machine Learning Toolbox) |

Re-running `tests/r0_passrule_stats.py` on this repository reproduces every pass-rule number quoted in the paper. It also refreshes the diagnostic rows of the case `c30_k3_long`, because the long runs of that case were extended after the stored CSVs were written; those rows are not quoted in the paper.

## MATLAB setup

1. Download GeoPDEs 3.2.2 and the NURBS toolbox 1.4.3 (both GPL-3) into `external/geopdes`, or set `GEOPDES_ROOT`.
2. In MATLAB: `addpath reproduce; c = config(true);` puts GeoPDEs, `vendor/geopdes_topopt_functions`, `src/` and `tests/` on the path.

The original scripts in `tests/` are kept as they ran. Some still call `addpath` with the absolute paths of the authors' machine; on another machine those calls only warn, and `config(true)` provides the paths. These scripts write into `results/`, so run them on a copy of the repository if the stored results must stay untouched.

## Data

The binary results (`.mat`, `.npz`) and the large logs are listed in `LARGE_FILES.csv` with their SHA-256. They will be deposited as a Zenodo dataset when the paper is submitted; its DOI will be added here and in `CITATION.cff`.

## License

Not yet chosen; see `reproduce/LICENSE-NOTE.md` (proposed: GPL-3.0-or-later for the code, as required by the GeoPDEs extension files, and CC-BY-4.0 for data and figures).

## Contact

Guilherme Henrique da Silva, Laboratório de Vibrações e Acústica, Universidade Federal de Santa Catarina, dshguilherme@gmail.com
