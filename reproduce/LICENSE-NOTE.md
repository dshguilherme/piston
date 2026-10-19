# Licence note (DRAFT, for the authors to decide)

No licence has been chosen. Nothing in this folder grants rights yet. This note records the proposal from the reproducibility audit and the reasons for it.

| Part | Proposed licence | Reason |
|---|---|---|
| Code (MATLAB, Python, PowerShell in this package and the `src/`, `tests/` scripts it calls) | **GPL-3.0-or-later** | The MATLAB model calls GeoPDEs (GPL-3) and the NURBS toolbox (GPL-3+) directly, and the four author files in `geopdes-master/topopt/Functions/` extend GeoPDEs operators. Distributing them together with code under a more permissive licence would still leave the combined work under the GPL; choosing GPL-3.0-or-later for everything keeps the boundary simple. bempp-cl (MIT), NumPy, SciPy and Matplotlib (BSD-style) are compatible. |
| Data, figures and tables (CSV, JSON, `.dat`, `.mat`, `.npz`, logs, TikZ figure sources, rendered PDFs) | **CC-BY-4.0** | Standard for research data on Zenodo; requires attribution only. |

Points to settle before release:

1. **Confirm GPL-3.0-or-later for the code**, or obtain a different arrangement for the GeoPDEs-dependent parts. If the four GeoPDEs extension files are vendored (audit gap G03), they need a GPL header with the authors' copyright.
2. **Confirm CC-BY-4.0 for the data.** The TikZ sources are figure *sources*; treating them as data (CC-BY) or as code (GPL) is the authors' choice.
3. **MATLAB Student licence.** The MATLAB results were computed under a Student licence. Whether its terms restrict publishing results computed with it was not assessed.
4. **Third-party text.** No third-party code is copied into this package; GeoPDEs and NURBS are referenced, not shipped.
5. When decided, add `LICENSE` (code) and `LICENSE-DATA` (data) files with the full texts and update `README.md` and `CITATION.cff` (`license:` field).
