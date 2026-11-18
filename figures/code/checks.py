"""Grey-scale and deuteranopia simulation, contact sheets. colorspacious is not installed in the build environment,
so the deuteranopia simulation uses the Machado et al. (2009) matrix (severity 1.0) applied in linear RGB."""
import glob, os
import numpy as np
from PIL import Image, ImageDraw
from figcommon import OUT, METH, LOG, note

CHK = OUT + "checks/"
os.makedirs(CHK, exist_ok=True)
M_DEUT = np.array([[0.367322, 0.860646, -0.227968],
                   [0.280085, 0.672501, 0.047413],
                   [-0.011820, 0.042940, 0.968881]])


def srgb2lin(c):
    c = np.asarray(c, float)
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def lin2srgb(c):
    c = np.clip(c, 0, 1)
    return np.where(c <= 0.0031308, c * 12.92, 1.055 * c ** (1 / 2.4) - 0.055)


def deut(arr01):
    return lin2srgb(srgb2lin(arr01) @ M_DEUT.T)


def lab(rgb):
    lin = srgb2lin(rgb)
    M = np.array([[0.4124, 0.3576, 0.1805], [0.2126, 0.7152, 0.0722], [0.0193, 0.1192, 0.9505]])
    xyz = lin @ M.T / np.array([0.95047, 1.0, 1.08883])
    f = np.where(xyz > 0.008856, np.cbrt(xyz), 7.787 * xyz + 16 / 116)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])


def hex2rgb(h):
    return np.array([int(h[i:i + 2], 16) / 255 for i in (1, 3, 5)])


def palette_report():
    names = ["Skala-1.1", "wB97M-V", "B3LYP-D3BJ", "r2SCAN", "PBE"]
    lines = []
    cols = {n: hex2rgb(METH[n]["c"]) for n in names}
    cols["black (CCSD(T))"] = hex2rgb("#000000"); cols["grey (xTB)"] = hex2rgb("#7F7F7F")
    names = list(cols)
    lines.append("| pair | dE76 normal | dE76 deuteranopia | grey L* difference |")
    lines.append("|---|---|---|---|")
    worst = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = cols[names[i]], cols[names[j]]
            de_n = np.linalg.norm(lab(a) - lab(b)); de_d = np.linalg.norm(lab(deut(a)) - lab(deut(b)))
            dl = abs(lab(a)[0] - lab(b)[0])
            lines.append(f"| {names[i]} / {names[j]} | {de_n:.0f} | {de_d:.0f} | {dl:.0f} |")
            if de_d < 20 and dl < 12:
                worst.append((names[i], names[j], de_d, dl))
    return lines, worst


def make_sim(files):
    out = []
    for f in files:
        im = Image.open(f).convert("RGB")
        a = np.asarray(im, float) / 255
        g = (0.2126 * srgb2lin(a[..., 0]) + 0.7152 * srgb2lin(a[..., 1]) + 0.0722 * srgb2lin(a[..., 2]))
        g8 = (lin2srgb(g) * 255).astype(np.uint8)
        Image.fromarray(g8, "L").save(CHK + os.path.basename(f).replace(".png", "_grey.png"))
        d8 = (deut(a) * 255).astype(np.uint8)
        Image.fromarray(d8, "RGB").save(CHK + os.path.basename(f).replace(".png", "_deut.png"))
        out.append(os.path.basename(f))
    return out


def contact_sheet(files, path, ncols=4, tw=700, mode="RGB", label_suffix=""):
    thumbs = []
    for f in files:
        im = Image.open(f).convert(mode)
        r = tw / im.width
        im = im.resize((tw, max(1, int(im.height * r))), Image.LANCZOS)
        thumbs.append((os.path.basename(f), im))
    rowsH = []
    for i in range(0, len(thumbs), ncols):
        rowsH.append(max(t[1].height for t in thumbs[i:i + ncols]) + 22)
    W = ncols * (tw + 12) + 12
    H = sum(rowsH) + 12
    bg = 255
    sheet = Image.new(mode, (W, H), (255, 255, 255) if mode == "RGB" else 255)
    dr = ImageDraw.Draw(sheet)
    y = 12
    for ri, i in enumerate(range(0, len(thumbs), ncols)):
        for k, (n, im) in enumerate(thumbs[i:i + ncols]):
            x = 12 + k * (tw + 12)
            dr.text((x, y), n + label_suffix, fill=(0, 0, 0) if mode == "RGB" else 0)
            sheet.paste(im, (x, y + 14))
            dr.rectangle([x - 1, y + 13, x + tw, y + 14 + im.height], outline=(150, 150, 150) if mode == "RGB" else 150)
        y += rowsH[ri]
    sheet.save(path)
    return sheet.size
