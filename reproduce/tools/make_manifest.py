"""Regenerate reproduce/MANIFEST.csv for the files that exist now.

Input : paper_jsv/reproduce_plan/MANIFEST-draft.csv (audited draft) + a scan for new files.
Output: reproduce/MANIFEST.csv   (read-only with respect to every data file)

Env: GEOPDES_ROOT = geopdes-master folder (default <ROOT>/../../geopdes-master).
Run : python tools/make_manifest.py
"""
import collections
import csv
import fnmatch
import glob
import hashlib
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
import config as C

ROOT = C.ROOT
PKG = C.PKG
DRAFT = os.path.join(C.PAPER, "reproduce_plan", "MANIFEST-draft.csv")
OUTCSV = os.path.join(PKG, "MANIFEST.csv")
GEO = os.path.abspath(os.environ.get("GEOPDES_ROOT") or os.path.join(ROOT, os.pardir, os.pardir, "geopdes-master"))
EXT_PREFIX = "EXTERNAL:geopdes-master/"
COLS = ["file", "role", "paper_item", "producer", "bytes", "sha256", "status_vs_audit",
        "git_tracked_now", "distribution", "draft_distribution"]
MB = 1024 * 1024
D_GIT, D_REL, D_REG, D_EXC, D_EXT = ("git", "release asset / Zenodo dataset", "regenerate (not shipped)",
                                     "exclude", "external (GeoPDEs folder; vendor)")
FROZEN_LEGACY = {"results/articleB_r0_thesisdamp_spherecap_40x24_f100_partial.mat",
                 "results/articleB_r0_thesisdamp_spherecap_40x24_f100.txt",
                 "results/t3_longhorizon.mat", "results/t3_longhorizon.txt"}


def abspath_of(rel):
    if rel.startswith(EXT_PREFIX):
        return os.path.join(GEO, rel[len(EXT_PREFIX):].replace("/", os.sep))
    return os.path.join(ROOT, rel.replace("/", os.sep))


def sha_size(p):
    try:
        h = hashlib.sha256()
        n = 0
        with open(p, "rb") as f:
            while True:
                b = f.read(4 * MB)
                if not b:
                    break
                h.update(b)
                n += len(b)
        return n, h.hexdigest()
    except OSError:
        return None, None


def git_files():
    out = subprocess.run(["git", "-C", ROOT, "ls-files", "-z"], capture_output=True, check=True).stdout
    return set(x.decode("utf-8", "replace") for x in out.split(b"\0") if x)


def classify_new(rel):
    ext = rel.rsplit(".", 1)[-1].lower() if "." in rel.rsplit("/", 1)[-1] else ""
    if rel.startswith("reproduce/"):
        return "package file"
    if rel.startswith("paper_jsv/prototype/figs/data/"):
        return "figure data"
    if rel.startswith("paper_jsv/prototype/figs/") and ext in ("tex", "sty", "py"):
        return "figure source (tikz)" if ext in ("tex", "sty") else "figure source (script)"
    if ext in ("mat", "npz", "npy"):
        return "large binary"
    if ext in ("log", "err", "txt"):
        return "log"
    if ext in ("csv", "json"):
        return "result data"
    return "other"


def paper_item_new(rel):
    if re.search(r"(^|[/_])f(19|20|21|22)([_.]|$)", rel.rsplit("/", 1)[-1]) and "figs/" in rel:
        return "in-progress figures F19-F22 (not in audited draft)"
    return ""


def distribution(rel, size, draft_dist):
    if rel.startswith(EXT_PREFIX):
        return D_EXT
    low = rel.lower()
    base = low.rsplit("/", 1)[-1]
    ext = base.rsplit(".", 1)[-1] if "." in base else ""
    if "/paper_d2/" in "/" + low or "/longruns/" in low or "/robust/" in low or ".bak" in base:
        return D_EXC
    if draft_dist.startswith("exclude (") or draft_dist == "exclude":
        return D_EXC
    if fnmatch.fnmatch(base, "cap_basis_*.npz") or base.endswith("_checkpoint.npy") or fnmatch.fnmatch(base, "zop_*.npz"):
        return D_REG
    if rel in FROZEN_LEGACY:
        return D_GIT
    if draft_dist.startswith("release asset"):
        return D_REL
    if ext in ("mat", "npz", "npy", "pdf", "png") and size is not None and size > MB:
        return D_REL
    return D_GIT


def scan_new():
    found = set()

    def add_abs(p):
        if os.path.isfile(p):
            found.add(os.path.relpath(p, ROOT).replace(os.sep, "/"))

    for dp, dn, fn in os.walk(C.RES):
        for f in fn:
            add_abs(os.path.join(dp, f))
    for pat in ("articleB_r0_thesisdamp_spherecap_40x24_f100_partial.*", "t3_longhorizon.*"):
        for p in glob.glob(os.path.join(C.RESULTS, pat)):
            add_abs(p)
    for p in glob.glob(os.path.join(C.FIGDATA_STORED, "*")):
        add_abs(p)
    for ext in ("tex", "sty", "py"):
        for p in glob.glob(os.path.join(C.FIGS_STORED, "*." + ext)):
            add_abs(p)
    for ext in ("log", "csv"):
        for p in glob.glob(os.path.join(C.SUPPLEMENT, "*." + ext)):
            add_abs(p)
    for dp, dn, fn in os.walk(PKG):
        dn[:] = [d for d in dn if d != "__pycache__" and not (os.path.abspath(os.path.join(dp, d)) == os.path.abspath(C.OUT))]
        for f in fn:
            p = os.path.join(dp, f)
            if os.path.abspath(p) == os.path.abspath(OUTCSV):
                continue
            add_abs(p)
    return found


def mb(n):
    return "%.2f MB (%.3f GB)" % (n / MB, n / MB / 1024)


def main():
    tracked = git_files()
    with open(DRAFT, newline="", encoding="utf-8") as f:
        draft = list(csv.DictReader(f))
    rows = []
    draft_files = set()
    changed, missing = [], []
    for d in draft:
        rel = d["file"]
        draft_files.add(rel)
        p = abspath_of(rel)
        size, sha = (sha_size(p) if os.path.isfile(p) else (None, None))
        if sha is None:
            status = "missing"
            missing.append(rel)
        elif sha == d["sha256"]:
            status = "unchanged"
        else:
            status = "changed"
            changed.append((rel, d["bytes"], size))
        is_ext = rel.startswith(EXT_PREFIX)
        rows.append({
            "file": rel, "role": d["role"], "paper_item": d["paper_item"], "producer": d["producer"],
            "bytes": "" if size is None else size, "sha256": sha or "", "status_vs_audit": status,
            "git_tracked_now": "no" if is_ext else ("yes" if rel in tracked else "no"),
            "distribution": distribution(rel, size if size is not None else int(d["bytes"] or 0), d["proposed_distribution"]),
            "draft_distribution": d["proposed_distribution"],
        })
    new = sorted(scan_new() - draft_files)
    print("scanning done: %d draft rows, %d new files; hashing..." % (len(draft), len(new)), flush=True)
    for rel in new:
        size, sha = sha_size(abspath_of(rel))
        if sha is None:
            continue
        rows.append({
            "file": rel, "role": classify_new(rel), "paper_item": paper_item_new(rel), "producer": "",
            "bytes": size, "sha256": sha, "status_vs_audit": "new",
            "git_tracked_now": "yes" if rel in tracked else "no",
            "distribution": distribution(rel, size, ""), "draft_distribution": "",
        })
    rows.sort(key=lambda r: r["file"])
    with open(OUTCSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=COLS, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    print("\nSUMMARY  (%d rows -> %s)" % (len(rows), C.rel(OUTCSV)))
    print("status counts:")
    for k, v in sorted(collections.Counter(r["status_vs_audit"] for r in rows).items()):
        print("  %-10s %d" % (k, v))
    print("distribution counts / total size:")
    cnt = collections.Counter(r["distribution"] for r in rows)
    tot = collections.Counter()
    for r in rows:
        tot[r["distribution"]] += int(r["bytes"] or 0)
    for k in sorted(cnt):
        print("  %-36s %5d files  %s" % (k, cnt[k], mb(tot[k])))
    print("git_tracked_now: " + ", ".join("%s=%d" % kv for kv in sorted(collections.Counter(r["git_tracked_now"] for r in rows).items())))
    print("missing files (%d):" % len(missing))
    for m in missing:
        print("  " + m)
    print("changed files (%d), old -> new bytes:" % len(changed))
    for rel, o, n in changed:
        print("  %s  %s -> %s" % (rel, o, n))


if __name__ == "__main__":
    main()
