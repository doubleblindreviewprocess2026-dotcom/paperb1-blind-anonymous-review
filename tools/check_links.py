#!/usr/bin/env python3
"""Static check of the links in the site folder (no browser needed).

    python tools/check_links.py

1. every relative link, image, script and stylesheet in every .html file points to a file (or folder with index.html) that exists,
   and every #anchor points to an id on the target page;
2. no absolute http(s) URLs are loaded by a page (no external scripts, styles, fonts or images); links to other sites are listed.
Exit code 0 = clean.
"""
import os, re, sys
from html.parser import HTMLParser

SITE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []; self.ids = set()

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if "id" in a:
            self.ids.add(a["id"])
        for k in ("href", "src"):
            if k in a and a[k] is not None:
                self.links.append((tag, k, a[k]))


def resolve(base_file, url):
    path = url.split("#")[0].split("?")[0]
    if path == "":
        return base_file
    full = os.path.normpath(os.path.join(os.path.dirname(base_file), path))
    if os.path.isdir(full):
        full = os.path.join(full, "index.html")
    return full


def main():
    problems = []
    html_files = []
    for root, dirs, files in os.walk(SITE):
        dirs[:] = [d for d in dirs if d not in ("__pycache__", ".git")]
        for f in files:
            if f.endswith(".html"):
                html_files.append(os.path.join(root, f))
    ids = {}
    parsed = {}
    for h in html_files:
        p = Links(); p.feed(open(h, encoding="utf-8").read())
        ids[h] = p.ids; parsed[h] = p
    n_links = 0; external = set()
    for h, p in parsed.items():
        rel = os.path.relpath(h, SITE)
        for tag, k, url in p.links:
            if url.startswith(("mailto:", "javascript:", "data:")):
                continue
            n_links += 1
            if re.match(r"^[a-z]+://", url) or url.startswith("//"):
                external.add(url)
                if tag in ("script", "link", "img", "iframe", "source"):
                    problems.append("%s: external resource %s" % (rel, url))
                continue
            if url == "#":
                continue   # placeholder href that common.js fills in (Colab and GitHub buttons)
            target = resolve(h, url)
            if not os.path.exists(target):
                problems.append("%s: broken link %s -> %s" % (rel, url, os.path.relpath(target, SITE)))
                continue
            if "#" in url:
                frag = url.split("#", 1)[1]
                if frag and target.endswith(".html"):
                    if target not in ids:
                        q = Links(); q.feed(open(target, encoding="utf-8").read()); ids[target] = q.ids
                    if frag not in ids[target]:
                        problems.append("%s: anchor #%s not found in %s" % (rel, frag, os.path.relpath(target, SITE)))
    print("%d html files, %d links checked, %d external links (not loaded by pages)" % (len(html_files), n_links, len(external)))
    for e in sorted(external):
        print("  external link:", e)
    if problems:
        print("PROBLEMS:")
        for x in problems:
            print("  ", x)
        sys.exit(1)
    print("clean")


if __name__ == "__main__":
    main()
