# drivers/

Glue that existed only in a temporary session scratchpad during the study, rescued (audit gap G05) and made path-independent. None of these was run while building the package.

| File | What it is | Caution |
|---|---|---|
| `pilot_queue.ps1` | the queue runner of the regime study (Experiment 6): one fresh MATLAB session per `cell:code:nmax` item, `maxNumCompThreads(2)`, skips runs already saved, appends to `queue.log` | calls the unedited `tests/r0_pilot_runs.m`, which writes into `results/r0_revision/pilot/`; run it with `REPRO_ROOT` pointing at a scratch copy of the repository |
| `pilot_queue.txt` | the item list, reconstructed from `results/r0_revision/pilot/queue.log` and the run logs (38 final runs; superseded attempts commented with the reason) | the `-List` strings actually passed were not recorded |
| `build_manuscript.ps1` | pdflatex, bibtex, pdflatex x2 on the manuscript, with error and warning counts | needs the paper source (`paper_jsv/prototype/`), which is not part of the package |
| `profiling/prof_*.m` | timing helpers for the sparse solves (LU variants, thread counts) | not used for any paper number; load release assets `pilot/factor_c30_k3.mat`, `ops_c30_k3.mat` |

There are no stage drivers that wrap the full recomputation (audit gap G07 remains open): the original MATLAB and Python scripts write into `results/` and carry absolute paths, and wrapping them without editing them was out of scope. `docs/figure_to_command.md` lists the commands per stage in order.
