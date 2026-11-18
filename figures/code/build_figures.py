#!/usr/bin/env python3
"""Build all figures, tables, captions, checks and the build log for Paper B1.
Run from the repository root (or from anywhere):

    python figures/code/build_figures.py [--out DIR]

Takes about 1-2 minutes. Reads only the files in data/ of this repository (data/csv/study, data/csv/results, data/records,
data/geometries). Writes PNG, PDF, SVG, the plotted CSV files, tables/, checks/, CAPTIONS.md and FIGURE_BUILD_LOG.md to DIR
(default: figures_rebuilt/ in the repository root; the published files in figures/png, pdf, svg are never overwritten).
Then compare:  python figures/code/compare_with_published.py
Order follows VISUALS_PLAN.md section 5: fallback set (Figures 2, 3, 5, 7) first, then the rest."""
import sys, os, glob, traceback, datetime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figcommon import *
import figcommon
import figs_fallback as FB, figs_rest as FR, figs_extra as FX, tables as TB, checks as CK, captions as CP

D = load()
ERR = []


def run(f):
    try:
        f(D)
    except Exception:
        ERR.append((f.__name__, traceback.format_exc()))
        print(traceback.format_exc())


print("##### fallback set")
for f in (FB.fig02, FB.fig03, FB.fig05, FB.fig07):
    run(f)
print("##### rest of main text")
for f in (FB.fig06, FX.fig01, FR.fig04, FR.fig08, FR.fig09, FR.fig10, FR.fig11, FR.fig12, FX.graphical_abstract):
    run(f)
print("##### supplementary")
for f in (FX.s1, FX.s2, FX.s3, FX.s4, FX.s5, FX.s6, FX.s7, FX.s8, FX.s9):
    run(f)
print("##### tables")
try:
    TB.build_tables(D)
except Exception:
    ERR.append(("tables", traceback.format_exc())); print(traceback.format_exc())

# ---- captions
CP.write_captions(OUT + "CAPTIONS.md")

# ---- colour checks
pngs = sorted(p for p in glob.glob(OUT + "*.png") if "contact_sheet" not in p)
CK.make_sim(pngs)
pal_lines, worst = CK.palette_report()
order = (sorted(p for p in pngs if os.path.basename(p).startswith("fig") and not os.path.basename(p).startswith("figS"))
         + [p for p in pngs if "graphical" in p] + sorted(p for p in pngs if os.path.basename(p).startswith("figS")))
sz = CK.contact_sheet(order, OUT + "contact_sheet.png", ncols=4, tw=650)
CK.contact_sheet([CK.CHK + os.path.basename(p).replace(".png", "_grey.png") for p in order], CK.CHK + "contact_sheet_grey.png", ncols=4, tw=650, mode="L")
CK.contact_sheet([CK.CHK + os.path.basename(p).replace(".png", "_deut.png") for p in order], CK.CHK + "contact_sheet_deuteranopia.png", ncols=4, tw=650)

# ---- log
nfig = len([p for p in pngs])
with open(OUT + "FIGURE_BUILD_LOG.md", "w") as f:
    f.write(f"# Figure build log, Paper B1\n\nBuilt {datetime.datetime.now():%Y-%m-%d %H:%M} by `figures/code/build_figures.py` (matplotlib {matplotlib.__version__}, Agg backend, DejaVu Sans, 300 dpi).\n\n")
    f.write(f"Files: {len([p for p in pngs if 'figS' not in p and 'graphical' not in p])} main figures, 1 graphical abstract, {len([p for p in pngs if 'figS' in p])} supplementary figures, "
            f"{len(glob.glob(TAB + '*.csv'))} table CSV files with Markdown versions in `figures/tables/`.\n\n")
    f.write("## Build errors\n\n")
    f.write("None.\n\n" if not ERR else "".join(f"- {n}:\n```\n{t}\n```\n" for n, t in ERR))
    nf = [x for x in LOG if not x[4]]
    f.write(f"## Anchor and consistency checks (VISUALS_PLAN.md section 5 item 3)\n\n{len(LOG) - len(nf)} of {len(LOG)} checks passed; {len(nf)} failed.\n\n")
    f.write("| figure | check | expected | got | result |\n|---|---|---|---|---|\n")
    for fg, what, ex, got, ok in LOG:
        gv = f"{got:.4f}" if isinstance(got, (float, np.floating)) else str(got)
        f.write(f"| {fg} | {what} | {ex} | {gv} | {'ok' if ok else 'FAIL'} |\n")
    f.write("\n## Deviations from the plan and observations\n\n")
    f.write("- Resolution: PNG at 300 dpi as requested (plan text says 600 dpi). Output folder layout: code in `figures/code/` (`build_figures.py` is the single entry point and imports `figcommon.py`, `figs_fallback.py`, `figs_rest.py`, `figs_extra.py`, `tables.py`, `checks.py`, `captions.py`) instead of `W/figures/make_figures.py`.\n")
    f.write("- Fonts: DejaVu Sans (Arial not installed). Tick labels, legends and annotations 8 pt, axis labels 9 pt, panel letters 9 pt bold, nothing below 8 pt. Labels use Unicode sub/superscripts (NH₃, F⁻, kcal mol⁻¹, r²SCAN) instead of mathtext so that they are a full 8 pt at final size.\n")
    f.write("- Legend colours and markers follow section 0.1 in every figure; Okabe-Ito palette, bars with 0.4 pt black edges.\n")
    f.write("- Single-panel Figures 2 and 9 have no panel letter.\n")
    seen = set()
    for fg, t in NOTES:
        if (fg, t) in seen:
            continue
        seen.add((fg, t)); f.write(f"- {fg}: {t}\n")
    f.write("- Fig 3: height 92 mm (plan 80 mm) so that the rotated tick labels of (c) and the legend fit. Figures 4, 8 and 10 keep the planned sizes; Fig 11 uses rotated tick labels.\n")
    f.write("- Fig 5 and S3: lines of the CCSD(T) + MP2 reference and of the error curves are broken across missing points (no interpolation over the gap), instead of straight lines as in the plan.\n")
    f.write("- Fig 4: reaction order follows section 0.2 (h2co before sn2_ff), the plan text in Fig 4 lists h2co after sn2_clcl; ch3nc is last.\n")
    f.write("- Fig 10: plan wording 'symmetric log or a broken axis is not needed' respected (linear -12 to +28).\n")
    f.write("- Table 2 is written twice: `table2_barriers` (formatted strings, barrier and error in one cell) and `table2_barriers_numeric` (all numbers in separate columns). Table S6 is transcribed by hand from dossier A.1, A.3, A.4 (no CSV source).\n")
    f.write("- Not built, as listed in section 4 of the plan: per-class MAE bars, old Fig. 6 map, CH4 stretch for Skala, bootstrap violins, regression plots, 2D SN2 map, H2CO path in main text, pooled-row MAE bars, Skala-rank bar chart, plots of the H2CO electronic error, radar chart, absolute total energies.\n\n")
    f.write("## Colour-blind and grey-scale check\n\n")
    f.write("`colorspacious` is not installed (and could not be installed in this environment). The deuteranopia simulation uses the Machado et al. (2009) matrix (severity 1.0) applied in linear RGB (`code/checks.py`); simulated copies of every PNG are in `figures/checks/` (`*_grey.png`, `*_deut.png`) with contact sheets `checks/contact_sheet_grey.png` and `checks/contact_sheet_deuteranopia.png`.\n\n")
    f.write("Colour distances (CIE76 delta E, normal vision and simulated deuteranopia) and L* differences (grey scale) between the method colours:\n\n" + "\n".join(pal_lines) + "\n\n")
    if worst:
        f.write("Pairs close in both deuteranopia (dE < 20) and grey level (dL < 12), so that only the marker shape separates them: " + "; ".join(f"{a} / {b} (dE {d:.0f}, dL {l:.0f})" for a, b, d, l in worst) + ". All data points also carry distinct marker shapes (circle, square, triangle up, diamond, triangle down, plus, cross, star), so the methods stay separable; this was confirmed on the simulated contact sheets (see visual notes below).\n\n")
    else:
        f.write("No pair of method colours is close in both the deuteranopia simulation and grey level.\n\n")
    qa = os.path.join(os.path.dirname(os.path.abspath(__file__)), "visual_qa_notes.md")
    f.write("## Visual review of the PNG files\n\n")
    f.write(open(qa).read() if os.path.exists(qa) else "(not yet written)\n")
print("done; contact sheet", sz, "; checks failed:", len([x for x in LOG if not x[4]]), "; build errors:", len(ERR))
