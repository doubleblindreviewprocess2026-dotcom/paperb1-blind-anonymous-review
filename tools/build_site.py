#!/usr/bin/env python3
"""Rebuild the pages and JSON files of the site from the CSV files in data/.

    python tools/build_site.py                 rebuild data/json/*.json, every page, SHA256SUMS.txt and b1_geometries.zip
    python tools/build_site.py --zip PATH      also write a zip of the whole site to PATH (the folder content at the top level of the zip)
    python tools/build_site.py --sync --writeup W --results R --notebooks N --manuscript M
                                               first copy the sources into this folder (see tools/sync_sources.py), then rebuild

Needs Python 3.8+ and numpy. Reads only files inside this folder (unless --sync). A column that is missing from tools/data_dictionary.py stops the build.
Pages: index.html, results/, demos/, figures/, data/, code/, checks/, about/ (each a folder with index.html, so that the addresses end in a slash).
Set `anonymous: true` in js/config.js and rebuild to add <meta name="robots" content="noindex"> to every page (the script common.js adds it too).
"""
import argparse, hashlib, os, re, sys, zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
SITE = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import build_data, build_pages  # noqa: E402
from sitelib import strip_png, strip_svg  # noqa: E402


def write(rel, text):
    path = os.path.join(SITE, rel)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    return path


def canon_bytes(path):
    """file content as it is published: PNG files without metadata chunks and SVG files without a provenance block (some environments add such a block to every image they write)"""
    with open(path, "rb") as f:
        data = f.read()
    if path.lower().endswith(".png"):
        return strip_png(data)
    if path.lower().endswith(".svg"):
        return strip_svg(data)
    return data


def sha256(path):
    return hashlib.sha256(canon_bytes(path)).hexdigest()


def anonymous_flag():
    s = open(os.path.join(SITE, "js", "config.js"), encoding="utf-8").read()
    return bool(re.search(r"anonymous\s*:\s*true", s))


def zip_geometries():
    gdir = os.path.join(SITE, "data", "geometries")
    out = os.path.join(SITE, "b1_geometries.zip")
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(os.listdir(gdir)):
            if f.endswith(".xyz"):
                zi = zipfile.ZipInfo("b1_geometries/" + f, date_time=(2026, 11, 16, 0, 0, 0))
                zi.compress_type = zipfile.ZIP_DEFLATED
                zi.external_attr = 0o644 << 16
                with open(os.path.join(gdir, f), "rb") as fh:
                    z.writestr(zi, fh.read())
    return out


def write_sums():
    lines = []
    for top in ("data", "code", "figures", "b1_geometries.zip"):
        p = os.path.join(SITE, top)
        if os.path.isfile(p):
            lines.append((top, sha256(p)))
            continue
        for root, dirs, files in os.walk(p):
            dirs.sort()
            for f in sorted(files):
                full = os.path.join(root, f)
                if "__pycache__" in full:
                    continue
                lines.append((os.path.relpath(full, SITE).replace(os.sep, "/"), sha256(full)))
    lines.sort()
    write("SHA256SUMS.txt", "".join("%s  %s\n" % (h, p) for p, h in lines))
    return len(lines)


def make_zip(dest):
    skip_dirs = {"__pycache__", ".git"}
    n = 0
    with zipfile.ZipFile(dest, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for root, dirs, files in os.walk(SITE):
            dirs[:] = sorted(d for d in dirs if d not in skip_dirs)
            for f in sorted(files):
                full = os.path.join(root, f)
                if os.path.abspath(full) == os.path.abspath(dest) or f.endswith(".pyc"):
                    continue
                arc = os.path.relpath(full, SITE).replace(os.sep, "/")
                zi = zipfile.ZipInfo(arc, date_time=(2026, 11, 16, 12, 0, 0))
                zi.compress_type = zipfile.ZIP_DEFLATED
                zi.external_attr = 0o644 << 16
                z.writestr(zi, canon_bytes(full), compresslevel=9)
                n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sync", action="store_true")
    ap.add_argument("--writeup"); ap.add_argument("--results"); ap.add_argument("--notebooks"); ap.add_argument("--manuscript")
    ap.add_argument("--zip", metavar="PATH")
    a = ap.parse_args()
    if a.sync:
        import sync_sources
        sys.argv = ["sync_sources.py", "--writeup", a.writeup, "--results", a.results, "--notebooks", a.notebooks, "--manuscript", a.manuscript]
        sync_sources.main()
    anon = anonymous_flag()
    res = build_data.build_all()
    R = res["results"]
    write("index.html", build_pages.index_page(R, anon))
    write("results/index.html", build_pages.results_page(R, anon))
    write("demos/index.html", build_pages.demos_page(anon))
    write("figures/index.html", build_pages.figures_page(anon))
    write("data/index.html", build_pages.data_page(anon))
    write("code/index.html", build_pages.code_page(anon))
    write("checks/index.html", build_pages.checks_page(anon))
    write("about/index.html", build_pages.about_page(anon))
    open(os.path.join(SITE, ".nojekyll"), "w").close()
    zip_geometries()
    n = write_sums()
    print("pages written (anonymous=%s); %d files in SHA256SUMS.txt" % (anon, n))
    if a.zip:
        k = make_zip(a.zip)
        print("zip: %d files, %.2f MB -> %s" % (k, os.path.getsize(a.zip) / 1e6, a.zip))


if __name__ == "__main__":
    main()
