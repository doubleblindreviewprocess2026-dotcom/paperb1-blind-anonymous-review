#!/usr/bin/env python3
"""Re-run the automatic check of Paper B1 on the published calculation records and compare it with the published report.

    python code/recheck_numbers.py [--verbose] [--keep DIR]

What it does
  1. copies data/csv/results/phase*.csv and the small JSON records (data/records/*.json) into a temporary folder,
  2. runs the check module of the notebooks (code/modules/b1v5_check.py and its helpers, executed in one namespace in the order of the
     notebook cells) on that folder,
  3. compares the new check_report.csv with data/csv/results/check_report.csv, cell by cell.
Expected result: 139 items (PASS 54, WARN 14, FAIL 5, NOT RUN 40, N/A 4, INFO 22) and "identical".
Needs only Python 3 and numpy. Runs in a few seconds. Nothing in the site folder is changed.
"""
import argparse, contextlib, csv, glob, io, os, shutil, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SITE = os.path.dirname(HERE)
MODS = os.path.join(HERE, "modules")
# the order of the notebook cells (common helpers, CPU helpers, geometry pack, stored numbers, literature table, check, figures, design, design check)
ORDER = ["b1v5_common.py", "b1v5_cpu.py", "b1v5_geompack.py", "b1v5_v4_stored.py", "b1v5_lit_table.py", "b1v5_check.py",
         "b1v5_gpu_figs.py", "b1v5_1h.py", "b1v5_1h_check.py"]


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.reader(f))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--verbose", action="store_true", help="print the full output of the check")
    ap.add_argument("--keep", metavar="DIR", help="work in DIR and keep the files (default: a temporary folder)")
    a = ap.parse_args()
    try:
        import numpy  # noqa: F401
    except ImportError:
        sys.exit("numpy is needed:  pip install numpy")

    work = a.keep or tempfile.mkdtemp(prefix="b1_recheck_")
    os.makedirs(work, exist_ok=True)
    for f in glob.glob(os.path.join(SITE, "data", "csv", "results", "phase*.csv")) + glob.glob(os.path.join(SITE, "data", "records", "*.json")):
        shutil.copy(f, work)

    ns = {"__name__": "__main__"}
    out = io.StringIO()
    with contextlib.redirect_stdout(out if not a.verbose else sys.stdout):
        for f in ORDER:
            with open(os.path.join(MODS, f), encoding="utf-8") as fh:
                exec(compile(fh.read(), f, "exec"), ns)
        ns["setup_root"](work, mount=False)
        ns["set_design"](ns["DESIGN_REACTIONS"])
        rep = ns["run_check"](work)

    new = read_csv(os.path.join(work, "check_report.csv"))
    old = read_csv(os.path.join(SITE, "data", "csv", "results", "check_report.csv"))
    diffs = [(i, o, n) for i, (o, n) in enumerate(zip(old, new)) if o != n]
    same = len(old) == len(new) and not diffs
    print("items: %d (published report: %d)" % (len(new) - 1, len(old) - 1))
    print("counts:", ", ".join("%s %d" % (k, v) for k, v in rep["counts"].items()))
    if same:
        print("check_report.csv: identical to the published file")
    else:
        print("check_report.csv: %d differing rows, length %d against %d" % (len(diffs), len(new), len(old)))
        for i, o, n in diffs[:10]:
            print("  row %d\n    published: %s\n    new:       %s" % (i, o, n))
    if not a.keep:
        shutil.rmtree(work, ignore_errors=True)
    sys.exit(0 if same else 1)


if __name__ == "__main__":
    main()
