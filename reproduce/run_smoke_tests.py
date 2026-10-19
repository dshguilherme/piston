"""Smoke tests of the reproducibility package. Target: under 10 minutes, no MATLAB solve, no BEM solve.

Usage:  python reproduce/run_smoke_tests.py [--only NAME[,NAME...]] [--matlab] [--no-tex]
        (Python: the analysis environment, env/environment-analysis.yml; the stored logs were made with it.)

Each test prints PASS, FAIL, SKIP or INFO. INFO tests report facts (versions, known manuscript issues,
in-progress figures) and never fail the run. Results and timings go to docs/smoke-test-log.md; all files
written by the tests go under $REPRO_OUT (default reproduce/out). Exit code 1 if any test FAILs.

  --matlab   also start MATLAB once to check that config.m puts the model code on the path
             (no solve; off by default because long MATLAB jobs may be running on the machine)
  --no-tex   skip the pdflatex compile check
"""
import argparse
import csv
import datetime
import filecmp
import hashlib
import os
import platform
import re
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C  # noqa: E402

PY = sys.executable
RESULTS = []


def run(cmd, cwd=None, timeout=600, env=None):
    t = time.time()
    e = dict(os.environ, PYTHONIOENCODING="utf-8", MPLBACKEND="Agg")
    if env:
        e.update(env)
    p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, encoding="utf-8", errors="replace",
                       timeout=timeout, env=e)
    return p.returncode, p.stdout, p.stderr, time.time() - t


def py(script, *args, timeout=600, log=None):
    rc, out, err, dt = run([PY, script, *args], cwd=C.out("smoke"), timeout=timeout)
    if log:
        with open(log, "w", encoding="utf-8", newline="\n") as f:
            f.write(out)
            if err.strip():
                f.write("\n--- stderr ---\n" + err)
    return rc, out, err, dt


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def norm_eol(b):
    return b.replace(b"\r\n", b"\n")


def cmp_text(a, b, rtol=0.0):
    """'byte-identical', 'identical after CRLF->LF', 'numerically equal (rtol)' or None if different."""
    ba, bb = open(a, "rb").read(), open(b, "rb").read()
    if ba == bb:
        return "byte-identical"
    if norm_eol(ba) == norm_eol(bb):
        return "identical after CRLF->LF"
    if rtol > 0:
        ta, tb = norm_eol(ba).decode().split(), norm_eol(bb).decode().split()
        if len(ta) != len(tb):
            return None
        for x, y in zip(ta, tb):
            if x == y:
                continue
            xs, ys = x.strip(",").split(","), y.strip(",").split(",")
            if len(xs) != len(ys):
                return None
            for u, v in zip(xs, ys):
                if u == v:
                    continue
                try:
                    fu, fv = float(u), float(v)
                except ValueError:
                    return None
                if abs(fu - fv) > rtol * max(abs(fu), abs(fv), 1e-300):
                    return None
        return "numerically equal (rtol %.0e)" % rtol
    return None


def norm_log(text):
    """Remove run-dependent parts of a log: timings and absolute 'wrote' paths."""
    out = []
    for ln in text.replace("\r\n", "\n").split("\n"):
        if ln.startswith("wrote "):
            continue
        ln = re.sub(r"\d+(\.\d+)?\s?s\b", "<t>s", ln)
        out.append(ln.rstrip())
    return "\n".join(out).strip()


def record(name, status, seconds, detail):
    RESULTS.append((name, status, seconds, detail))
    print("[%s] %-34s %6.1f s  %s" % (status, name, seconds, detail.split("\n")[0]), flush=True)


def test(name):
    def deco(fn):
        fn.test_name = name
        return fn
    return deco


# ----------------------------------------------------------------------------- tests

@test("env_versions")
def t_env():
    want = {"python": "3.13.12", "numpy": "2.5.2", "scipy": "1.18.1", "matplotlib": "3.11.2"}
    got = {"python": platform.python_version()}
    for m in ("numpy", "scipy", "matplotlib"):
        try:
            got[m] = __import__(m).__version__
        except ImportError:
            got[m] = "missing"
    diff = [k for k in want if got[k] != want[k]]
    tex = shutil.which(os.environ.get("REPRO_PDFLATEX", "pdflatex")) or "not found"
    detail = ", ".join("%s %s" % (k, got[k]) for k in want) + "; pdflatex: " + tex
    if diff:
        detail += "; differs from the pinned analysis env in: " + ", ".join(diff) + " (byte identity of logs may fail)"
    return "INFO", detail


@test("config_paths")
def t_config():
    need = {"REPRO_ROOT": C.ROOT, "results/r0_revision": C.RES, "tests": C.TESTS, "src": C.SRC,
            "stored figure data": C.FIGDATA_STORED, "tikz copies": C.TIKZ}
    miss = [k for k, v in need.items() if not os.path.isdir(v)]
    return ("FAIL" if miss else "PASS"), ("missing: " + ", ".join(miss) if miss else "REPRO_ROOT and stored results found; REPRO_OUT = %s" % C.rel(C.OUT))


@test("factor_two")
def t_factor2():
    d = C.out("smoke", "r0_nov_sphere_factor2")
    rc, out, err, dt = py(os.path.join(C.PKG, "tests", "copies", "r0_nov_sphere_factor2.py"), log=os.path.join(d, "stdout.txt"))
    if rc:
        return "FAIL", "exit %d: %s" % (rc, err.strip()[-300:])
    ratios = [float(x) for x in re.findall(r"ratio=([0-9.]+)", out)]
    worst = max(abs(r - 2) for r in ratios) if ratios else 1
    m = re.findall(r"max \|.*?\|.*?: ([0-9.eE+-]+)", out)
    ok_ratio = len(ratios) == 4 and worst < 1e-11
    csvcmp = cmp_text(os.path.join(d, "sphere_factor2.csv"), os.path.join(C.NOV, "sphere_factor2.csv"), rtol=1e-12)
    logcmp = norm_log(out) == norm_log(open(os.path.join(C.NOV, "sphere_factor2.log"), encoding="utf-8").read())
    ok = ok_ratio and csvcmp is not None
    return ("PASS" if ok else "FAIL"), ("P_sinc/P_mono = 2 to %.1e at k = 0.3, 1, 3, 8; sphere identities %s; CSV vs stored: %s; "
                                        "log vs stored: %s" % (worst, "/".join(m), csvcmp or "DIFFERENT",
                                                               "identical" if logcmp else "differs"))


@test("sphere_modal_formulas")
def t_sphere():
    d = C.out("smoke", "sphere_modal")
    rc, out, err, dt = py(os.path.join(C.PKG, "tests", "smoke_sphere_modal.py"), log=os.path.join(d, "stdout.txt"))
    nok = out.count("\nok ") + out.startswith("ok ")
    return ("PASS" if rc == 0 else "FAIL"), "%d checks ok against sphere-modal-spectrum.csv, direct-check, exterior benchmark%s" % (
        nok, "" if rc == 0 else "; " + (out + err).strip()[-300:])


@test("flat_disc_extrapolation")
def t_flat():
    d = C.out("smoke", "r0_nov_flat_extrap")
    rc, out, err, dt = py(os.path.join(C.PKG, "tests", "copies", "r0_nov_flat_extrap.py"), log=os.path.join(d, "stdout.txt"))
    if rc:
        return "FAIL", "exit %d: %s" % (rc, err.strip()[-300:])
    new = os.path.join(d, "cap_flat_extrap.csv")
    c = cmp_text(new, os.path.join(C.NOV, "cap_flat_extrap.csv"), rtol=1e-12)
    with open(new, newline="") as f:
        rel = [abs(float(r["rel_diff"])) for r in csv.DictReader(f)]
    paper = [3.9e-4, 5.4e-4, 3.2e-5, 1.6e-4]
    ok2 = all(float("%.1e" % a) == b for a, b in zip(rel, paper)) and len(rel) == 4
    return ("PASS" if (c and ok2) else "FAIL"), "|1 - cap/disc| = %s (paper 3.9e-4, 5.4e-4, 3.2e-5, 1.6e-4); CSV vs stored: %s" % (
        ", ".join("%.2e" % a for a in rel), c or "DIFFERENT")


@test("fluid_operator_checks")
def t_fsi():
    d = C.out("smoke", "r0_fsi_verify")
    rc, out, err, dt = py(os.path.join(C.PKG, "tests", "copies", "r0_fsi_verify.py"), timeout=300, log=os.path.join(d, "stdout.txt"))
    if rc:
        return "FAIL", "exit %d: %s" % (rc, err.strip()[-300:])
    m = re.search(r"(\d+)/(\d+) checks passed", out)
    c = cmp_text(os.path.join(d, "verify.csv"), os.path.join(C.FSI, "verify.csv"), rtol=1e-10)
    logcmp = norm_log(out) == norm_log(open(os.path.join(C.FSI, "verify.log"), encoding="utf-8").read())
    ok = bool(m) and m.group(1) == m.group(2) and c is not None
    return ("PASS" if ok else "FAIL"), "%s checks passed; verify.csv vs stored: %s; log vs stored (timings removed): %s" % (
        m.group(0).split()[0] if m else "?", c or "DIFFERENT", "identical" if logcmp else "differs")


@test("bem_from_stored_data")
def t_bem():
    d = C.out("smoke", "bem_stored")
    rc, out, err, dt = py(os.path.join(C.PKG, "tests", "smoke_bem_stored.py"), log=os.path.join(d, "stdout.txt"))
    nok = sum(1 for ln in out.splitlines() if ln.startswith("ok "))
    return ("PASS" if rc == 0 else "FAIL"), "%d checks ok (sphere closed form, disc asymptote, monotone convergence, far/bem); no solve%s" % (
        nok, "" if rc == 0 else "; " + (out + err).strip()[-300:])


@test("exact_operator_checks_assets")
def t_optchecks():
    need = ["run_S0.mat", "run_X0.mat", "run_XS.mat", "run_SX.mat", "run_SR.mat", "run_XR.mat", "cap_path_states.mat",
            "cap_exact_factor_q8.mat", "cap_exact_factor_q12.mat", "cap_samples_q8.mat", "cap_path_powers.npz"]
    miss = [n for n in need if not os.path.exists(os.path.join(C.EXO, n))]
    if miss:
        return "SKIP", "release assets not present: " + ", ".join(miss)
    d = C.out("smoke", "r0_opt_cap_checks")
    rc, out, err, dt = py(os.path.join(C.PKG, "tests", "copies", "r0_opt_cap_checks.py"), timeout=300, log=os.path.join(d, "stdout.txt"))
    if rc:
        return "FAIL", "exit %d: %s" % (rc, err.strip()[-300:])
    produced = os.path.join(d, "cap_reopt_checks.log")
    text = open(produced, encoding="utf-8").read() if os.path.exists(produced) else out
    same = norm_log(text) == norm_log(open(os.path.join(C.EXO, "cap_reopt_checks.log"), encoding="utf-8").read())
    return ("PASS" if same and "R0_OPT_CAP_CHECKS_DONE" in text else "FAIL"), (
        "q8/q12 quadrature, convergence tails, efficiencies from stored designs and F_x; log vs stored: %s" % ("identical" if same else "DIFFERENT"))


DAT23 = ["f13_S0", "f13_SR", "f13_SX", "f13_X0", "f13_XR", "f13_XS", "f17_pt_erp", "f17_pt_exact", "f17_pt_sinc", "f17_regret",
         "f18_fsi", "f6_bem", "f6_flat_ka0.1", "f6_flat_ka1", "f6_flat_ka10", "f6_flat_ka3", "f8_cells_cos0", "f8_cells_cos7",
         "f8_geomean_designs30", "f8_geomean_modes", "f8_rho1_cos0", "f8_rho1_cos7", "f9_shares"]


@test("replot_figure_data")
def t_replot():
    d = C.out("smoke", "replot")
    for s in ("make_fig_data_d2.py", "make_fig_pilot_f17.py", "make_fig_fsi_f18.py"):
        rc, out, err, dt = py(os.path.join(C.PKG, "figures", s), log=os.path.join(d, s + ".txt"))
        if rc:
            return "FAIL", "%s exit %d: %s" % (s, rc, err.strip()[-300:])
    od = os.path.join(C.OUT, "figs", "data")
    kinds, bad = {}, []
    for n in DAT23:
        c = cmp_text(os.path.join(od, n + ".dat"), os.path.join(C.FIGDATA_STORED, n + ".dat"))
        if c is None:
            bad.append(n)
        else:
            kinds[c] = kinds.get(c, 0) + 1
    summary = "; ".join("%d %s" % (v, k) for k, v in sorted(kinds.items()))
    return ("FAIL" if bad else "PASS"), "23 plotted .dat files: %s%s" % (summary, "; DIFFERENT: " + ", ".join(bad) if bad else "")


@test("replot_ratio_figure")
def t_ratio():
    d = C.out("smoke", "replot")
    rc, out, err, dt = py(os.path.join(C.PKG, "figures", "make_fig_ratio_kR.py"), log=os.path.join(d, "make_fig_ratio_kR.txt"))
    if rc:
        return "FAIL", "exit %d: %s" % (rc, err.strip()[-300:])
    png_new, png_old = os.path.join(C.OUT, "figs", "fig_ratio_kR_proto.png"), os.path.join(C.FIGS_STORED, "fig_ratio_kR_proto.png")
    pdf_new, pdf_old = os.path.join(C.OUT, "figs", "fig_ratio_kR_proto.pdf"), os.path.join(C.FIGS_STORED, "fig_ratio_kR_proto.pdf")
    png_same = filecmp.cmp(png_new, png_old, shallow=False)
    strip = lambda b: re.sub(rb"/CreationDate \(D:[^)]*\)", b"", b)
    pdf_same = strip(open(pdf_new, "rb").read()) == strip(open(pdf_old, "rb").read())
    status = "PASS" if png_same else "FAIL"
    return status, "PNG %s; PDF %s" % ("byte-identical" if png_same else "DIFFERENT (font or library version?)",
                                       "identical except /CreationDate" if pdf_same else "differs beyond /CreationDate")


@test("replot_extras_F19_F22")
def t_extras():
    d = C.out("smoke", "replot")
    msgs = []
    for s in ("r0_fig_geometry.py", "r0_fig_designs.py", "r0_fig_mechanism.py"):
        rc, out, err, dt = py(os.path.join(C.PKG, "figures", "extras", s), timeout=400, log=os.path.join(d, s + ".txt"))
        msgs.append("%s %s (%.0f s)" % (s, "ok" if rc == 0 else "exit %d" % rc, dt))
    od = os.path.join(C.OUT, "figs", "data")
    same = diff = 0
    for n in sorted(os.listdir(od)):
        if re.match(r"f(19|20|21|22)_", n) and os.path.exists(os.path.join(C.FIGDATA_STORED, n)):
            if cmp_text(os.path.join(od, n), os.path.join(C.FIGDATA_STORED, n)):
                same += 1
            else:
                diff += 1
    return "INFO", "in-progress figures, not in the audited draft: %s; %d data files equal to stored, %d differ" % ("; ".join(msgs), same, diff)


@test("tables")
def t_tables():
    d = C.out("smoke", "tables")
    msgs, fail = [], False
    for s in ("make_tab_setup.py", "make_tab_opt.py", "make_tab_pilot.py", "make_tab_fsi.py"):
        rc, out, err, dt = py(os.path.join(C.PKG, "tables", s), log=os.path.join(d, s + ".txt"))
        fail |= rc != 0
        m = [x.split("; wrote")[0] for x in re.findall(r"(tab:\w+: [^\n]*)", out)]
        msgs.append((m[0] if m else s) + ("" if rc == 0 else " [exit %d]" % rc))
    return ("FAIL" if fail else "PASS"), "; ".join(msgs)


AUDITED_TIKZ = ["fig_cap_problem", "fig_factor_two", "fig_model_axes", "fig_regime_map", "fig_opt_paths", "fig_pilot_regime",
                "fig_fsi", "fig_decision", "fig_sinc_construction", "fig_verification", "fig_face_power", "fig_bench_cards"]


@test("tikz_compile")
def t_tex(args):
    if args.no_tex:
        return "SKIP", "--no-tex"
    tex = shutil.which(os.environ.get("REPRO_PDFLATEX", "pdflatex"))
    if not tex:
        return "SKIP", "pdflatex not found (set REPRO_PDFLATEX or put TeX Live on PATH)"
    wd = C.out("smoke", "tex")
    dd = os.path.join(wd, "data")
    os.makedirs(dd, exist_ok=True)
    for n in os.listdir(C.FIGDATA_STORED):
        shutil.copy2(os.path.join(C.FIGDATA_STORED, n), os.path.join(dd, n))
    regen = os.path.join(C.OUT, "figs", "data")
    nreg = 0
    if os.path.isdir(regen):
        for n in os.listdir(regen):
            shutil.copy2(os.path.join(regen, n), os.path.join(dd, n))
            nreg += 1
    for n in os.listdir(C.TIKZ):
        shutil.copy2(os.path.join(C.TIKZ, n), os.path.join(wd, n))
    ok, bad, extra = [], [], []
    names = AUDITED_TIKZ + sorted(os.path.splitext(n)[0] for n in os.listdir(C.TIKZ)
                                  if n.endswith(".tex") and os.path.splitext(n)[0] not in AUDITED_TIKZ)
    for name in names:
        rc, out, err, dt = run([tex, "-interaction=nonstopmode", "-halt-on-error", name + ".tex"], cwd=wd, timeout=120)
        good = rc == 0 and os.path.exists(os.path.join(wd, name + ".pdf"))
        if name in AUDITED_TIKZ:
            (ok if good else bad).append(name)
        else:
            extra.append("%s %s" % (name, "ok" if good else "FAILED"))
    status = "FAIL" if bad else "PASS"
    return status, "%d/%d audited TikZ figures compile with regenerated data (%d regenerated files overlaid on stored data)%s; extras: %s" % (
        len(ok), len(AUDITED_TIKZ), nreg, "; FAILED: " + ", ".join(bad) if bad else "", ", ".join(extra) or "none")


@test("flagged_manuscript_numbers")
def t_flagged():
    sw = list(csv.DictReader(open(os.path.join(C.NOV, "cap_frozen_sweep.csv"), newline="")))
    v = [float(r["ratio"]) for r in sw if 10 <= float(r["kR"]) <= 27]
    pe = [r for r in csv.DictReader(open(os.path.join(C.NOV, "cap_modes_pencil.csv"), newline="")) if r["ratio_sinc_optimal"] != "nan"]
    pm = {th: min(float(r["ratio_sinc_optimal"]) for r in pe if r["theta0_deg"] == th) for th in ("120", "150")}
    fs = [r for r in csv.DictReader(open(os.path.join(C.FSI, "fsi_summary.csv"), newline="")) if float(r["kR"]) <= 3.01]
    loss = {"dev": [], "phys": []}
    for c in sorted({r["cell"] for r in fs}):
        d = {r["design"]: r for r in fs if r["cell"] == c}
        best = lambda var: min(float(d["XU"]["dB_vs_uniform_" + var]), float(d["XS"]["dB_vs_uniform_" + var]))
        for var in loss:
            loss[var].append(best(var) - best("vac"))
    return "INFO", ("l.250 'within 0.43 to 0.70 up to kR=27': data give %.3f to %.3f (suggest 0.39 to 0.70); "
                    "l.530 'rho between 0.011 and 0.048': per-cap minima %.4f (150 deg), %.4f (120 deg), not a range; "
                    "fig_decision.tex:33 'lost 2--30 dB': scaled medium %.1f-%.1f dB, air %.1f-%.1f dB (name the medium)") % (
        min(v), max(v), pm["150"], pm["120"], min(loss["dev"]), max(loss["dev"]), min(loss["phys"]), max(loss["phys"]))


@test("matlab_config")
def t_matlab(args):
    if not args.matlab:
        return "SKIP", "not requested (--matlab); no MATLAB started"
    exe = shutil.which(os.environ.get("MATLAB_EXE", "matlab"))
    if not exe:
        return "SKIP", "matlab not found"
    code = ("cd('%s'); c=config(true); fns={'obstacle_build_ctx','op_KL_shells_elements','op_u_v_elements','gcmma_solver',"
            "'r0_build_consistent','sp_bspline','nrbmak'}; for i=1:numel(fns), w=which(fns{i}); fprintf('WHICH %%s|%%s\\n',fns{i},w); end; "
            "fprintf('GEOPDES %%s\\n',c.GEOPDES_ROOT);") % C.PKG.replace("\\", "/")
    rc, out, err, dt = run([exe, "-batch", code], timeout=300)
    found = re.findall(r"WHICH (\w+)\|(.*)", out)
    miss = [f for f, w in found if not w.strip()]
    ok = rc == 0 and found and not miss
    return ("PASS" if ok else "FAIL"), "config(true) resolves %d/%d model functions%s" % (
        len(found) - len(miss), len(found), "; missing: " + ", ".join(miss) if miss else "")


TESTS = [t_env, t_config, t_factor2, t_sphere, t_flat, t_fsi, t_bem, t_optchecks, t_replot, t_ratio, t_extras, t_tables,
         t_tex, t_flagged, t_matlab]


def write_log(total):
    p = os.path.join(C.DOCS, "smoke-test-log.md")
    counts = {}
    for _, s, _, _ in RESULTS:
        counts[s] = counts.get(s, 0) + 1
    with open(p, "w", encoding="utf-8", newline="\n") as f:
        f.write("# Smoke-test log\n\n")
        f.write("Generated by `reproduce/run_smoke_tests.py` on %s. Do not edit by hand; re-run to update.\n\n" %
                datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
        f.write("- Python %s (%s), platform %s\n" % (platform.python_version(), os.path.basename(PY), platform.platform()))
        f.write("- REPRO_ROOT = `%s`, REPRO_OUT = `%s`\n" % (C.rel(C.ROOT) or ".", C.rel(C.OUT)))
        f.write("- Total wall time: %.1f s (%.1f min). Counts: %s\n\n" % (total, total / 60, ", ".join("%s %d" % kv for kv in sorted(counts.items()))))
        f.write("| Test | Status | Time (s) | Result |\n|---|---|---|---|\n")
        for n, s, t, d in RESULTS:
            f.write("| %s | %s | %.1f | %s |\n" % (n, s, t, d.replace("|", "\\|").replace("\n", " ")))
        f.write("\nOutputs of every test are under `$REPRO_OUT/smoke/` (stdout of each script in `stdout.txt`).\n")
    print("wrote", p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=[])
    ap.add_argument("--matlab", action="store_true")
    ap.add_argument("--no-tex", action="store_true")
    args = ap.parse_args()
    only = {x for a in args.only for x in a.split(",") if x}
    t0 = time.time()
    for fn in TESTS:
        if only and fn.test_name not in only:
            continue
        ts = time.time()
        try:
            res = fn(args) if fn.__code__.co_argcount else fn()
        except subprocess.TimeoutExpired as e:
            res = ("FAIL", "timeout after %s s" % e.timeout)
        except Exception as e:  # a crashing test is a failed test, reported with its message
            res = ("FAIL", "%s: %s" % (type(e).__name__, e))
        record(fn.test_name, res[0], time.time() - ts, res[1])
    total = time.time() - t0
    nfail = sum(1 for r in RESULTS if r[1] == "FAIL")
    print("\n%d tests, %d FAIL, total %.1f s" % (len(RESULTS), nfail, total))
    if not only:
        write_log(total)
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
