"""Shared style, data loading, helpers for the Paper B1 figures (see VISUALS_PLAN.md section 0)."""
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# All paths are relative to the repository root, which is found from the location of this file (figures/code/ -> two levels up),
# so the script runs from any working directory, for example:  python figures/code/build_figures.py
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) + os.sep
STUDY = REPO + "data/csv/study/"      # master_*, summary_stats, path_profiles, sn2_symscan, derived_*   (was writeup/)
RES = REPO + "data/csv/results/"      # calculation records and tables of the run                         (was results_2026-11-16/)
REC = REPO + "data/records/"          # phase05_hcn_grid.npz, phase13_xtb_barriers.json
GEO = REPO + "data/geometries/"       # xyz files

import argparse
_ap = argparse.ArgumentParser(add_help=False)
_ap.add_argument("--out", default=REPO + "figures_rebuilt")
_args, _ = _ap.parse_known_args()
# the figures are written to a NEW folder (default figures_rebuilt/ in the repository root, ignored by git), never over figures/png|pdf|svg
OUT = os.path.abspath(_args.out) + os.sep
TAB = OUT + "tables/"
os.makedirs(OUT, exist_ok=True)
os.makedirs(TAB, exist_ok=True)

DPI = 300  # 300 dpi as requested (the plan text says 600)
MM = 1 / 25.4
SINGLE, DOUBLE = 85 * MM, 170 * MM

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "mathtext.fontset": "dejavusans",
    "font.size": 8,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.dpi": DPI,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "svg.fonttype": "path",
    "figure.dpi": 100,
    "lines.solid_capstyle": "butt",
})

# ---------------------------------------------------------------- 0.1 method table
ORANGE = "#D55E00"; BLUE = "#0072B2"; GREEN = "#009E73"; PINK = "#CC79A7"; SKY = "#56B4E9"
METH = {
    "Skala-1.1":      dict(label="Skala-1.1", c=ORANGE, m="o", ls="-", lw=1.8, open=False),
    "Skala-1.1-noD3": dict(label="Skala-1.1 without D3", c=ORANGE, m="o", ls="--", lw=1.2, open=True),
    "wB97M-V":        dict(label="ωB97M-V", c=BLUE, m="s", ls="-", lw=1.5, open=False),
    "B3LYP-D3BJ":     dict(label="B3LYP-D3(BJ)", c=GREEN, m="^", ls="-", lw=1.2, open=False),
    "r2SCAN":         dict(label="r²SCAN", c=PINK, m="D", ls="-", lw=1.2, open=False),
    "PBE":            dict(label="PBE", c=SKY, m="v", ls="-", lw=1.2, open=False),
    "CCSD(T)+dMP2":   dict(label="CCSD(T) + MP2 basis correction", c="#000000", m="P", ls="-", lw=1.5, open=False),
    "CCSD(T)-fc":     dict(label="CCSD(T), def2-TZVP", c="#000000", m="x", ls=":", lw=1.0, open=False),
    "GFN2-xTB":       dict(label="GFN2-xTB", c="#7F7F7F", m="*", ls="-.", lw=1.0, open=False),
}
DFT6 = ["Skala-1.1", "Skala-1.1-noD3", "wB97M-V", "B3LYP-D3BJ", "r2SCAN", "PBE"]
DFT5 = ["Skala-1.1", "wB97M-V", "B3LYP-D3BJ", "r2SCAN", "PBE"]
MS = {"o": 5.2, "s": 4.8, "^": 5.6, "D": 4.4, "v": 5.6, "P": 7.0, "x": 5.5, "*": 8.0}

RX = ["hcn", "nh3", "h2co", "sn2_ff", "sn2_clcl"]
RX6 = RX + ["ch3nc"]
RXL = {"hcn": "HCN → HNC", "nh3": "NH₃ inversion", "h2co": "H₂CO → H₂ + CO",
       "sn2_ff": "F⁻ + CH₃F", "sn2_clcl": "Cl⁻ + CH₃Cl", "ch3nc": "CH₃NC → CH₃CN"}
KC = "kcal mol⁻¹"
BAND = "#EEEEEE"
GREY = "#7F7F7F"


def mstyle(meth, ms_scale=1.0, **kw):
    """kwargs for ax.plot marker styling of a method."""
    s = METH[meth]
    d = dict(marker=s["m"], color=s["c"], ms=MS[s["m"]] * ms_scale)
    if s["m"] in ("x", "P", "+"):
        d.update(mew=1.5 if s["m"] == "P" else 1.1, mfc=s["c"], mec=s["c"])
        if s["m"] == "P":
            d.update(marker="+", ms=8, mew=2.0)
    elif s["open"]:
        d.update(mfc="white", mec=s["c"], mew=1.1)
    else:
        d.update(mfc=s["c"], mec="black" if s["c"] != "#000000" else s["c"], mew=0.4)
    d.update(kw)
    return d


def handle(meth, line=True, ms_scale=1.0, **kw):
    s = METH[meth]
    st = mstyle(meth, ms_scale)
    return Line2D([0], [0], ls=s["ls"] if line else "none", lw=s["lw"] if line else 0,
                  label=s["label"], **{**st, **kw})


def legend_handles(methods, line=True):
    return [handle(m, line) for m in methods]


def panel(ax, letter, x=0.02, y=0.97, **kw):
    ax.text(x, y, f"({letter})", transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="left", **kw)


def zero_line(ax, horizontal=True):
    (ax.axhline if horizontal else ax.axvline)(0, color="black", lw=0.8, zorder=2)


def band(ax, w=1.0, horizontal=True):
    if horizontal:
        ax.axhspan(-w, w, color=BAND, zorder=0, lw=0)
    else:
        ax.axvspan(-w, w, color=BAND, zorder=0, lw=0)


def rx_ticklabels(ax, rxs=RX, dagger=True, grey_last=False, rotate=0):
    labs = []
    for r in rxs:
        t = RXL[r]
        if r == "h2co" and dagger:
            t += "†"
        labs.append(t)
    ax.set_xticks(range(len(rxs)))
    tl = ax.set_xticklabels([l.replace(" → ", "\n→ ").replace(" + ", "\n+ ") if rotate == 0 else l for l in labs],
                            rotation=rotate)
    if grey_last and rxs[-1] == "ch3nc":
        tl[-1].set_color(GREY)
    return tl


# ---------------------------------------------------------------- data
def load():
    D = {}
    D["bar"] = pd.read_csv(STUDY + "master_barriers.csv")
    D["err"] = pd.read_csv(STUDY + "master_errors.csv")
    D["sum"] = pd.read_csv(STUDY + "summary_stats.csv")
    D["path"] = pd.read_csv(STUDY + "path_profiles.csv")
    D["sym"] = pd.read_csv(STUDY + "sn2_symscan.csv")
    for n in ["basis_species", "ccsdt_diagnostics", "internal_reference", "leave_one_out",
              "sn2_difference", "sn2_minima", "timing_long", "timing_pairs"]:
        D[n] = pd.read_csv(STUDY + f"derived_{n}.csv")
    return D


# ---------------------------------------------------------------- logging / saving
LOG = []      # anchor checks: (figure, what, expected, got, ok)
NOTES = []    # (figure, text) deviations and observations
DATA_FILES = []


def anchor(fig, what, expected, got, tol=0.006):
    ok = abs(float(expected) - float(got)) <= tol
    LOG.append((fig, what, expected, got, ok))
    print(f"  [{'ok' if ok else 'FAIL'}] {fig}: {what}: expected {expected}, got {got:.4f}")
    return ok


def note(fig, text):
    NOTES.append((fig, text))
    print(f"  NOTE {fig}: {text}")


def save(fig, name, formats=("png", "pdf", "svg")):
    for f in formats:
        try:
            fig.savefig(OUT + f"{name}.{f}", dpi=DPI)
        except Exception as e:  # pragma: no cover
            note(name, f"could not write .{f}: {e}")
    plt.close(fig)
    print(f"saved {name}")


def save_data(df, name):
    p = OUT + f"{name}_data.csv"
    df.to_csv(p, index=False)
    DATA_FILES.append(p)


def md_table(df, floatfmt=None):
    """Markdown table without external dependencies."""
    cols = list(df.columns)
    def fmt(v):
        if v is None or (isinstance(v, float) and np.isnan(v)):
            return ""
        if isinstance(v, (float, np.floating)):
            return f"{v:.4g}" if floatfmt is None else format(v, floatfmt)
        return str(v).replace("|", "/").replace("\n", " ")
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(fmt(r[c]) for c in cols) + " |")
    return "\n".join(lines)


def save_table(df, name, caption=None, floatfmt=None):
    df.to_csv(TAB + f"{name}.csv", index=False)
    with open(TAB + f"{name}.md", "w") as f:
        if caption:
            f.write(f"**{caption}**\n\n")
        f.write(md_table(df, floatfmt) + "\n")
    print(f"table {name}")


def fig_ax(w_mm, h_mm, **kw):
    return plt.subplots(figsize=(w_mm * MM, h_mm * MM), **kw)
