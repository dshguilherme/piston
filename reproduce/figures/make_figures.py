"""Regenerate the paper's figures from the stored results: data extraction, then the PDF.

Usage:  python reproduce/figures/make_figures.py [--paper] [--only fig_fsi fig_opt_paths ...] [--extras] [--no-pdf]

For each figure the driver runs its extractor (copies of the original scripts, paths parametrised) and
compiles the TikZ source with pdflatex in $REPRO_OUT/figs/build. Data go to $REPRO_OUT/figs/data, PDFs to
$REPRO_OUT/figs/pdf. Nothing is written into results/ or into the paper tree. Schematic figures need no data.
--extras also builds the in-progress figures F19-F22 (not in the audited draft).
--paper builds the 17 figures of the submitted manuscript from their sources in paper_jsv/prototype/figs. Their data
come from the extractors below and, for the pass-rule, long-run and fluid-loading figures, from the scripts
tests/r0_fig_passrule_v2.py, tests/r0_fig_longruns.py and tests/r0_fig_fsi_f18b.py, which write into the stored
data folder and are run separately.
"""
import argparse
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir))
import config as C  # noqa: E402

# figure -> (paper item, extractor scripts in figures/, data files it plots)
FIGS = {
    "fig_cap_problem": ("Fig. 1 fig:problem", [], []),
    "fig_factor_two": ("Fig. 2 fig:factortwo", [], []),
    "fig_model_axes": ("Fig. 3 fig:axes", [], []),
    "fig_ratio_kR_proto": ("Fig. 4 fig:ratiokR (matplotlib)", ["make_fig_ratio_kR.py"], []),
    "fig_regime_map": ("Fig. 5 fig:regime", ["make_fig_data_d2.py"], ["f8_*"]),
    "fig_opt_paths": ("Fig. 6 fig:paths", ["make_fig_data_d2.py"], ["f13_*"]),
    "fig_pilot_regime": ("Fig. 7 fig:pilot", ["make_fig_pilot_f17.py"], ["f17_*"]),
    "fig_fsi": ("Fig. 8 fig:fsi", ["make_fig_fsi_f18.py"], ["f18_fsi"]),
    "fig_decision": ("Fig. 9 fig:decision", [], []),
    "fig_sinc_construction": ("Fig. A.1 fig:construction", [], []),
    "fig_verification": ("Fig. B.1 fig:validation", ["make_fig_data_d2.py"], ["f6_*"]),
    "fig_face_power": ("Fig. D.1 fig:faces", ["make_fig_data_d2.py"], ["f9_shares"]),
    "fig_bench_cards": ("Fig. G.1 fig:cards", [], []),
}
EXTRAS = {
    "fig_cap_geometry": ("F19 (in progress)", ["extras/r0_fig_geometry.py"], ["f19_*"]),
    "fig_designs_paper": ("F20 (in progress)", ["extras/r0_fig_designs.py"], ["f20_*"]),
    "fig_designs_pilot": ("F21 (in progress)", ["extras/r0_fig_designs.py"], ["f21_*"]),
    "fig_mechanism": ("F22 (in progress)", ["extras/r0_fig_mechanism.py"], ["f22_*"]),
    "fig_pipeline_eng": ("pipeline schematic (in progress)", [], []),
}

PAPER = {
    "fig_cap_geometry": ("Fig. 1", ["extras/r0_fig_geometry.py"], ["f19_*"]),
    "fig_pipeline_eng": ("pipeline schematic", [], []),
    "fig_opt_paths_v2": ("paths from the uniform start", ["make_fig_data_d2.py"], ["f13_*"]),
    "fig_designs_paper_v2": ("design-point designs", ["extras/r0_fig_designs.py"], ["f20_*"]),
    "fig_mechanism": ("cancellation mechanism", ["extras/r0_fig_mechanism.py"], ["f22_*"]),
    "fig_pilot_regime_v2": ("regime study", ["make_fig_pilot_f17.py"], ["f17_*"]),
    "fig_passbars_v2": ("guaranteed level (data: tests/r0_fig_passrule_v2.py)", [], []),
    "fig_longruns": ("long runs (data: tests/r0_fig_longruns.py)", [], []),
    "fig_fsi_v2": ("fluid loading (data: tests/r0_fig_fsi_f18b.py)", [], []),
    "fig_verification": ("verification", ["make_fig_data_d2.py"], ["f6_*"]),
    "fig_regime_map_v2": ("ratio maps", ["make_fig_data_d2.py"], ["f8_*"]),
    "fig_regime_mean": ("ratio means", ["make_fig_data_d2.py"], ["f8_*"]),
    "fig_face_power_v2": ("face power", ["make_fig_data_d2.py"], ["f9_shares"]),
    "fig_designs_pilot_v2": ("low-kR designs", ["extras/r0_fig_designs.py"], ["f21_*"]),
    "fig_pilot_erp": ("ERP in the regime study", ["make_fig_pilot_f17.py"], ["f17_*"]),
    "fig_bench_cards": ("benchmark cards", [], []),
    "fig_passrule_scatter_v2": ("error against predicted level (data: tests/r0_fig_passrule_v2.py)", [], []),
}


def needs_dataset(stderr):
    return "FileNotFoundError" in stderr and any(x in stderr for x in (".mat'", ".npz'", ".npy'"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=[])
    ap.add_argument("--extras", action="store_true")
    ap.add_argument("--paper", action="store_true")
    ap.add_argument("--no-pdf", action="store_true")
    a = ap.parse_args()
    figs = PAPER if a.paper else dict(FIGS, **(EXTRAS if a.extras else {}))
    only = {x for s in a.only for x in s.split(",") if x}
    names = [n for n in figs if not only or n in only]
    unknown = only - set(figs)
    if unknown:
        sys.exit("unknown figure(s): %s; known: %s" % (", ".join(sorted(unknown)), ", ".join(figs)))

    done = set()
    for n in names:
        for s in figs[n][1]:
            if s in done:
                continue
            t = time.time()
            r = subprocess.run([sys.executable, os.path.join(HERE, s)], cwd=C.out("figs"), capture_output=True, text=True)
            print("extract %-28s exit %d (%.1f s)" % (s, r.returncode, time.time() - t))
            if r.returncode and needs_dataset(r.stderr):
                print("        %-28s needs the binary dataset (LARGE_FILES.csv); the stored data are used" % s)
            elif r.returncode:
                print(r.stderr[-2000:])
                sys.exit(1)
            done.add(s)
    if a.no_pdf:
        return
    tex = shutil.which(os.environ.get("REPRO_PDFLATEX", "pdflatex"))
    if not tex:
        sys.exit("pdflatex not found; set REPRO_PDFLATEX")
    build, pdf = C.out("figs", "build"), C.out("figs", "pdf")
    bdata = os.path.join(build, "data")
    os.makedirs(bdata, exist_ok=True)
    for src in (C.FIGDATA_STORED, os.path.join(C.OUT, "figs", "data")):  # regenerated data overwrite stored copies
        if os.path.isdir(src):
            for f in os.listdir(src):
                shutil.copy2(os.path.join(src, f), os.path.join(bdata, f))
    tikz = C.FIGS_STORED if a.paper else C.TIKZ
    for f in os.listdir(tikz):
        if f.endswith((".tex", ".sty")):
            shutil.copy2(os.path.join(tikz, f), os.path.join(build, f))
    for n in names:
        if n == "fig_ratio_kR_proto":
            shutil.copy2(os.path.join(C.OUT, "figs", n + ".pdf"), os.path.join(pdf, n + ".pdf"))
            print("pdf     %-28s from matplotlib" % n)
            continue
        t = time.time()
        r = subprocess.run([tex, "-interaction=nonstopmode", "-halt-on-error", n + ".tex"], cwd=build,
                           capture_output=True, text=True, errors="replace")
        ok = r.returncode == 0 and os.path.exists(os.path.join(build, n + ".pdf"))
        if ok:
            shutil.copy2(os.path.join(build, n + ".pdf"), os.path.join(pdf, n + ".pdf"))
        print("pdf     %-28s %s (%.1f s)  %s" % (n, "ok" if ok else "FAILED", time.time() - t, figs[n][0]))
    print("PDFs in", pdf)


if __name__ == "__main__":
    main()
