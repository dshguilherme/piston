"""Paths for the reproducibility package (Python side). No absolute paths anywhere in the package.

Environment variables (all optional):
  REPRO_ROOT  folder that holds results/, src/, tests/ and paper_jsv/ (in fastSIMP: harmonic_shell_course/).
              Default: the parent of the folder that holds this file.
  REPRO_OUT   folder that every re-run writes into. Default: <this folder>/out. Nothing in the package writes into
              results/ or into the paper tree.

Use from a script in a subfolder:
    import os, sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
    import config as C
"""
import os

PKG = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.environ.get("REPRO_ROOT") or os.path.join(PKG, os.pardir))
OUT = os.path.abspath(os.environ.get("REPRO_OUT") or os.path.join(PKG, "out"))

# stored inputs (read only)
RESULTS = os.path.join(ROOT, "results")                 # legacy frozen inputs live here
RES = os.path.join(RESULTS, "r0_revision")              # all r0 result files
NOV = os.path.join(RES, "novelty")
VIZ = os.path.join(RES, "viz")
EXO = os.path.join(RES, "exact_opt")
PIL = os.path.join(RES, "pilot")
FSI = os.path.join(RES, "fsi")
SRC = os.path.join(ROOT, "src")
TESTS = os.path.join(ROOT, "tests")
PAPER = os.path.join(ROOT, "paper_jsv")
FIGS_STORED = os.path.join(PAPER, "prototype", "figs")
FIGDATA_STORED = os.path.join(FIGS_STORED, "data")
SUPPLEMENT = os.path.join(PAPER, "supplement")

# package folders
TIKZ = os.path.join(PKG, "figures", "tikz")
DOCS = os.path.join(PKG, "docs")


def out(*parts):
    """Return a folder under REPRO_OUT, creating it."""
    p = os.path.join(OUT, *parts)
    os.makedirs(p, exist_ok=True)
    return p


def rel(p):
    """Path relative to REPRO_ROOT, with forward slashes (for logs and manifests)."""
    return os.path.relpath(os.path.abspath(p), ROOT).replace(os.sep, "/")


if __name__ == "__main__":
    for k in ("PKG", "ROOT", "OUT", "RES", "SRC", "TESTS", "FIGS_STORED", "TIKZ"):
        v = globals()[k]
        print("%-12s %s %s" % (k, v, "" if os.path.exists(v) else "(missing)"))
