#!/usr/bin/env python3
"""Remove everything from PNG files that is not needed to draw them (text chunks, time stamps, EXIF, embedded provenance blocks).

    python tools/strip_png_metadata.py figures/png/*.png figures/svg/*.svg

The image data is copied byte for byte, so the pictures are identical. build_site.py --zip applies the same step to the PNG files in the zip,
and SHA256SUMS.txt lists the hashes of the stripped files. Run this on the folder if you upload the files themselves instead of the zip.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sitelib import strip_png, strip_svg

for p in sys.argv[1:]:
    with open(p, "rb") as f:
        old = f.read()
    new = strip_svg(old) if p.lower().endswith(".svg") else strip_png(old)
    if new != old:
        with open(p, "wb") as f:
            f.write(new)
        print("stripped %s: %d -> %d bytes" % (p, len(old), len(new)))
