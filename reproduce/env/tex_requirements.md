# TeX requirements

The figures are standalone TikZ/pgfplots sources (`figures/tikz/fig_*.tex` with `jsvfigs.sty`); they read their data from a `data/` folder beside them. The ratio-versus-kR figure (Fig. 4) is matplotlib, not TeX.

| Item | Version used | Notes |
|---|---|---|
| Distribution | TeX Live 2025 | `C:\texlive\2025` on the authors' machine; any TeX Live 2025 should do |
| Engine | pdfTeX 3.141592653-2.6-1.40.27 | `pdflatex`; no figure needs LuaLaTeX |
| pgf / TikZ | 3.1.10 | libraries `arrows.meta, calc, patterns, decorations.markings, positioning, fit` (`jsvfigs.sty:14-15`); `shapes.geometric` in `fig_decision.tex`, `fig_opt_paths.tex` |
| pgfplots | 1.18.1, `compat=1.18` | `pgfplotstable` in `fig_regime_map`, `fig_opt_paths`, `fig_pilot_regime`, `fig_fsi`; library `patchplots` in `fig_regime_map` |
| Classes | `standalone` (figures); `elsarticle` 3.4c (2025/01/11), options `preprint,12pt` (manuscript) | |
| Manuscript packages | `amsmath, amssymb, bm, booktabs, graphicx, array, placeins [section], hyperref [hidelinks]`; bibliography style `elsarticle-num` | the manuscript is not part of the package |

Commands: `python run_smoke_tests.py --only tikz_compile` compiles all figure sources against the regenerated data in a scratch folder (`$REPRO_OUT/smoke/tex`). By hand: copy `figures/tikz/*` and a `data/` folder into one directory and run `pdflatex -interaction=nonstopmode -halt-on-error <name>.tex`. `drivers/build_manuscript.ps1` builds the manuscript (pdflatex, bibtex, pdflatex twice) when the paper source is present.

Set `REPRO_PDFLATEX` if `pdflatex` is not on `PATH`.
