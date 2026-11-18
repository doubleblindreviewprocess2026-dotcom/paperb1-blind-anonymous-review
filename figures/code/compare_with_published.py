#!/usr/bin/env python3
"""Compare a figure rebuild with the published files.

    python figures/code/build_figures.py            # writes figures_rebuilt/
    python figures/code/compare_with_published.py [--out DIR]

Checks, and prints one line each:
  * every *_data.csv in the rebuild is byte-identical to data/csv/figures/
  * every tables/*.csv in the rebuild is byte-identical to data/csv/tables/ (where the site has the file)
  * every PNG in the rebuild has the same pixels as figures/png/ (needs Pillow and numpy; PNG metadata chunks are ignored)
Exit status 0 when everything agrees. PDF and SVG files are not compared byte by byte (they contain creation dates).
"""
import argparse, filecmp, glob, os, sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ap = argparse.ArgumentParser()
ap.add_argument("--out", default=os.path.join(REPO, "figures_rebuilt"))
out = os.path.abspath(ap.parse_args().out)
bad = 0


def report(kind, name, ok, extra=""):
    global bad
    bad += (not ok)
    print(("ok    " if ok else "DIFF  ") + kind + " " + name + extra)


for f in sorted(glob.glob(os.path.join(out, "*_data.csv"))):
    o = os.path.join(REPO, "data", "csv", "figures", os.path.basename(f))
    report("figure data", os.path.basename(f), os.path.exists(o) and filecmp.cmp(f, o, shallow=False))
for f in sorted(glob.glob(os.path.join(out, "tables", "*.csv"))):
    o = os.path.join(REPO, "data", "csv", "tables", os.path.basename(f))
    if os.path.exists(o):
        report("table", os.path.basename(f), filecmp.cmp(f, o, shallow=False))
try:
    import numpy as np
    from PIL import Image
    for f in sorted(glob.glob(os.path.join(out, "*.png"))):
        b = os.path.basename(f)
        o = os.path.join(REPO, "figures", "png", b)
        if "contact_sheet" in b or not os.path.exists(o):
            continue
        a, c = np.asarray(Image.open(f).convert("RGBA")), np.asarray(Image.open(o).convert("RGBA"))
        report("png pixels", b, a.shape == c.shape and bool((a == c).all()))
except ImportError:
    print("Pillow or numpy missing: PNG pixels not compared")
print("all agree" if not bad else "%d difference(s)" % bad)
sys.exit(1 if bad else 0)
