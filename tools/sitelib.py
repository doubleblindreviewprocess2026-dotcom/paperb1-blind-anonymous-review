"""Small helpers for build_site.py: escaping, number formats, page layout, tables. Standard library only."""
import html, os, re, struct

MINUS = "−"

NAV = [("", "Home"), ("results", "Results"), ("demos", "Demos"), ("figures", "Figures"), ("data", "Data"), ("code", "Code"), ("checks", "Checks"), ("about", "About")]
SITE_NAME = "Paper B1 companion"


def esc(s):
    return html.escape(str(s), quote=True)


def fmt(x, d=2):
    """number with a true minus sign; empty or non-numeric gives an en dash"""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "–"
    if v != v:
        return "–"
    return ("%.*f" % (d, v)).replace("-", MINUS)


def sgn(x, d=2):
    """signed number: + or true minus; values that round to zero get no sign"""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "–"
    s = "%.*f" % (d, v)
    if float(s) == 0:
        return s.replace("-", "")
    return ("+" if v > 0 else MINUS) + s.replace("-", "")


def sign_class(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return ""
    return "pos" if v > 0 else ("neg" if v < 0 else "")


def human_size(n):
    for u in ("B", "kB", "MB", "GB"):
        if n < 1000 or u == "GB":
            return ("%d %s" % (n, u)) if u == "B" else ("%.1f %s" % (n, u) if n < 100 else "%d %s" % (n, u))
        n /= 1000.0


def png_size(path):
    with open(path, "rb") as f:
        head = f.read(32)
    if head[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    return struct.unpack(">II", head[16:24])


def table(headers, rows, cls="", caption=None, num_cols=(), row_cls=None, label=None):
    """headers: list of str (html allowed); rows: list of lists of str (html allowed). num_cols: indices right aligned."""
    out = ['<div class="tablewrap" tabindex="0"%s role="region" aria-label="%s">' % ("", esc(label or caption or "table"))]
    out.append('<table%s>' % (' class="%s"' % cls if cls else ""))
    if caption:
        out.append("<caption>%s</caption>" % caption)
    out.append("<thead><tr>%s</tr></thead>" % "".join('<th scope="col"%s>%s</th>' % (' class="num"' if i in num_cols else "", h) for i, h in enumerate(headers)))
    out.append("<tbody>")
    for ri, r in enumerate(rows):
        rc = (row_cls(ri, r) if row_cls else "") or ""
        out.append("<tr%s>%s</tr>" % (' class="%s"' % rc if rc else "", "".join(
            "<td%s>%s</td>" % (' class="num"' if i in num_cols else "", c) for i, c in enumerate(r))))
    out.append("</tbody></table></div>")
    return "\n".join(out)


def page(slug, title, body, scripts=(), desc="", anonymous=False):
    """Full HTML page. slug '' is the home page; others are directories one level below the root."""
    pre = "" if slug == "" else "../"
    nav = []
    for s, label in NAV:
        href = (pre + s + "/") if s else (pre or "./")
        cur = ' aria-current="page"' if s == slug else ""
        nav.append('<a href="%s"%s>%s</a>' % (href, cur, label))
    full_title = title if slug == "" else "%s | %s" % (title, SITE_NAME)
    robots = '<meta name="robots" content="noindex">\n' if anonymous else ""
    scripts_html = "\n".join('<script src="%s%s"></script>' % (pre, s) for s in ["js/config.js", "js/common.js"] + list(scripts))
    return """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>%(title)s</title>
<meta name="description" content="%(desc)s">
<meta name="color-scheme" content="light dark">
%(robots)s<link rel="stylesheet" href="%(pre)scss/style.css">
<link rel="icon" href="data:image/svg+xml,%%3Csvg xmlns=%%27http://www.w3.org/2000/svg%%27 viewBox=%%270 0 16 16%%27%%3E%%3Crect width=%%2716%%27 height=%%2716%%27 rx=%%273%%27 fill=%%27%%230b5f78%%27/%%3E%%3Cpath d=%%27M2 12 C5 12 6 3 8 3 S11 12 14 12%%27 fill=%%27none%%27 stroke=%%27white%%27 stroke-width=%%271.6%%27/%%3E%%3C/svg%%3E">
<script>try{var t=localStorage.getItem("b1-theme");if(t==="dark"||t==="light")document.documentElement.setAttribute("data-theme",t);}catch(e){}</script>
</head>
<body>
<a class="skip" href="#main">Skip to the content</a>
<header class="site"><div class="in">
<a class="brand" href="%(home)s">Paper B1 companion</a>
<nav aria-label="Main">%(nav)s</nav>
<button type="button" class="theme" aria-label="Switch theme">Dark theme</button>
</div></header>
<main id="main">
%(body)s
</main>
<footer class="site"><div class="in">
<p>Companion website and data for a study of reaction barriers and reaction paths with Skala-1.1. Code: MIT licence. Data and figures: CC BY 4.0 (see <a href="%(pre)sabout/#licence">About</a>).</p>
<p data-hide-anon>&copy; <span data-cfg="year"></span> <span data-cfg="holder"></span>. Questions: <span data-cfg="contact"></span>.</p>
<p data-show-anon>Anonymous copy for review.</p>
<p>The pages use no external services and set no cookies. The data files are static: every number on this site comes from a CSV file in <a href="%(pre)sdata/">data</a>.</p>
</div></footer>
<noscript><div class="note warn" style="max-width:1088px;margin:0 auto 20px"><p>JavaScript is off. The tables and numbers on the pages still show, but the interactive views (sorting, filters, demos) need JavaScript.</p></div></noscript>
%(scripts)s
</body>
</html>
""" % dict(title=esc(full_title), desc=esc(desc), robots=robots, pre=pre, home=(pre or "./"), nav="".join(nav), body=body, scripts=scripts_html)


def colab_buttons(path, pre):
    """buttons for a notebook path relative to the site root; the Colab and GitHub links are set by common.js from js/config.js"""
    return ('<a class="btn" data-colab="%s" href="#">Open in Colab</a> <a class="btn alt" data-github="%s" href="#">View on GitHub</a> '
            '<a class="btn alt" href="%s" download>Download .ipynb</a>' % (esc(path), esc(path), esc(pre + path)))


def slugify(s):
    return re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")


KEEP_PNG_CHUNKS = {b"IHDR", b"PLTE", b"IDAT", b"IEND", b"pHYs", b"tRNS", b"gAMA", b"sRGB", b"cHRM", b"sBIT"}


def strip_svg(data):
    """SVG bytes without an embedded provenance block (a <metadata><c2pa:manifest> element and its namespace declaration)"""
    data = re.sub(rb"<metadata>\s*<c2pa:manifest>.*?</c2pa:manifest>\s*</metadata>", b"", data, flags=re.S)
    return re.sub(rb'\s+xmlns:c2pa="[^"]*"', b"", data)


def strip_png(data):
    """PNG bytes without the chunks that are not needed to draw the image (text chunks, time, EXIF, embedded provenance blocks).
    The image data (IDAT) is copied byte for byte."""
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError("not a PNG file")
    out = [data[:8]]; i = 8
    while i < len(data):
        n = struct.unpack(">I", data[i:i + 4])[0]
        if data[i + 4:i + 8] in KEEP_PNG_CHUNKS:
            out.append(data[i:i + 12 + n])
        i += 12 + n
    return b"".join(out)
