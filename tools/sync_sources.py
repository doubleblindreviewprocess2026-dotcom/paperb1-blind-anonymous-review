#!/usr/bin/env python3
"""Copy the calculation records, tables, figures, geometries and notebooks of Paper B1 into this site folder.

Run it only when the sources change (the site itself is rebuilt from data/ by build_site.py):

    python tools/sync_sources.py --writeup PATH/writeup --results PATH/results --notebooks PATH/notebooks --manuscript PATH/writeup/manuscript

Nothing is read from anywhere else. What is copied:
  data/csv/study/     the study tables (master_*, summary_stats, derived_*, path_profiles, sn2_symscan)
  data/csv/tables/    the paper and supplementary tables as CSV
  data/csv/figures/   the data behind each figure
  data/csv/results/   the calculation records and reports of the run
  data/records/       small JSON records of the run (settings, tiers, timing, xTB barriers, geometry-pack check)
  data/geometries/    every .xyz file used, and one multi-frame .xyz per relaxed path
  figures/{png,pdf,svg}/   the figures; data/json/figure_captions.json holds the captions
  code/notebooks/     the notebooks with every output removed
  code/modules/       the Python modules that the notebooks paste into their cells
  code/geompack/      the geometry pack (JSON) with its SHA-256
Not copied: figures/code/ (the figure-build scripts, edited to read their inputs from data/ of this repository; they live only here).
Needs only the standard library and numpy.
"""
import argparse, csv, glob, hashlib, json, os, re, shutil, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SITE = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from sitelib import strip_png  # noqa: E402

STUDY = ["master_barriers", "master_errors", "summary_stats", "path_profiles", "sn2_symscan", "derived_basis_species",
         "derived_ccsdt_diagnostics", "derived_internal_reference", "derived_leave_one_out", "derived_sn2_difference",
         "derived_sn2_minima", "derived_timing_long", "derived_timing_pairs"]
RECORDS = ["config.json", "tiers.json", "phase10_timing.json", "phase13_xtb_barriers.json", "geompack_installed.json", "phase00_verify_summary.json"]
NOTEBOOKS = ["B1v5_CPU_1h_v2", "B1v5_GPU_1h_v3", "B1v5_collect_zip", "B1v5_GPU_diag_sn2", "B1v5_GPU_diag_sym"]
MODULES = ["b1v5_common", "b1v5_cpu", "b1v5_gpu", "b1v5_gpu_figs", "b1v5_1h", "b1v5_1h_cpu", "b1v5_1h_gpu", "b1v5_1h_check", "b1v5_check",
           "b1v5_lit_table", "b1v5_v4_stored", "b1v5_geompack", "b1v5_diag_sn2",
           "build_cpu_1h", "build_gpu_1h", "build_diag_sn2", "build_diag_sym", "build_collect_zip", "build_geompack",
           "check_test_synthetic", "test_1h_check", "test_1h_notebooks", "run_cells"]
# old file name -> number in the paper (main figures); supplementary figures keep their number
FIG_ORDER = [("fig01_study_design", 1), ("fig04_reference_quality", 2), ("fig02_signed_errors", 3), ("fig03_summary_stats", 4),
             ("fig06_nh3_path", 5), ("fig05_hcn_path", 6), ("fig07_sn2_stretch", 7), ("fig11_sn2_geometry", 8), ("fig08_basis_sensitivity", 9),
             ("fig09_d3_contribution", 10), ("fig10_scout", 11), ("fig12_cost_accuracy", 12)]
SI_FIGS = ["figS%d" % i for i in range(1, 10)]


def mk(p):
    os.makedirs(p, exist_ok=True)
    return p


def cp(src, dst):
    mk(os.path.dirname(dst))
    shutil.copyfile(src, dst)


def clean_png(src, dst):
    """copy a PNG without text, time, EXIF and provenance chunks (the image data is copied byte for byte)"""
    mk(os.path.dirname(dst))
    with open(src, "rb") as f:
        data = strip_png(f.read())
    with open(dst, "wb") as f:
        f.write(data)


def clean_svg(src, dst):
    """copy an SVG without its <metadata> block (creator string and date written by the plotting library)"""
    t = open(src, encoding="utf-8").read()
    t = re.sub(r"\s*<metadata>.*?</metadata>", "", t, flags=re.S)
    mk(os.path.dirname(dst))
    open(dst, "w", encoding="utf-8").write(t)


def clean_dir(path, pattern):
    """remove old generated files of one kind so that a re-sync leaves no stale file behind"""
    for f in glob.glob(os.path.join(path, pattern)):
        os.remove(f)


def strip_notebook(src, dst):
    nb = json.load(open(src, encoding="utf-8"))
    for c in nb.get("cells", []):
        if c.get("cell_type") == "code":
            c["outputs"] = []
            c["execution_count"] = None
        c.get("metadata", {}).pop("outputId", None)
        c.get("metadata", {}).pop("executionInfo", None)
    md = nb.get("metadata", {})
    md.pop("widgets", None)
    colab = md.get("colab")
    if isinstance(colab, dict):
        for k in list(colab):
            if k not in ("provenance", "gpuType", "name", "toc_visible"):
                colab.pop(k)
        colab.pop("provenance", None)
    with open(dst, "w", encoding="utf-8") as f:
        json.dump(nb, f, indent=1, ensure_ascii=False)
        f.write("\n")


def path_to_xyz(js, dst):
    d = json.load(open(js, encoding="utf-8"))
    lines = []
    for p in d["points"]:
        atoms = p["atoms"]
        lines.append(str(len(atoms)))
        lines.append("%s %s=%s point p%02d" % (d["system"], d.get("label", "coordinate").replace(" ", "_"), p["coordinate"], p["index"]))
        for el, xyz in atoms:
            lines.append("%-2s %14.8f %14.8f %14.8f" % (el, xyz[0], xyz[1], xyz[2]))
    with open(dst, "w") as f:
        f.write("\n".join(lines) + "\n")


def parse_captions(md_text, kind):
    """caption paragraphs '**Figure N. Title.** text' followed by '<!-- src: ..., figXX_name -->'"""
    out = {}
    paras = [p.strip() for p in re.split(r"\n\s*\n", md_text)]
    for i, p in enumerate(paras):
        m = re.match(r"\*\*Figure (S?\d+)\. (.+?)\*\*\s*(.*)", p, re.S)
        if not m:
            continue
        num, title, body = m.group(1), m.group(2).rstrip("."), m.group(3).strip()
        body = re.sub(r"\s+", " ", body)
        src = None
        if kind == "main" and i + 1 < len(paras):
            s = re.search(r"figures/CAPTIONS\.md, (fig\d+_\w+)", paras[i + 1])
            src = s.group(1) if s else None
        out[num] = {"title": title, "text": body, "src": src}
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--writeup", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--notebooks", required=True)
    ap.add_argument("--manuscript", required=True, help="folder with the checked captions (used for figure captions only)")
    a = ap.parse_args()
    D = os.path.join(SITE, "data")

    # ---- CSV
    for g in ("study", "tables", "figures", "results"):
        clean_dir(mk(os.path.join(D, "csv", g)), "*.csv")
    for n in STUDY:
        cp(os.path.join(a.writeup, n + ".csv"), os.path.join(D, "csv", "study", n + ".csv"))
    for f in sorted(glob.glob(os.path.join(a.writeup, "figures", "tables", "*.csv"))):
        cp(f, os.path.join(D, "csv", "tables", os.path.basename(f)))
    for f in sorted(glob.glob(os.path.join(a.writeup, "figures", "*_data.csv"))):
        cp(f, os.path.join(D, "csv", "figures", os.path.basename(f)))
    for f in sorted(glob.glob(os.path.join(a.results, "*.csv"))) + [os.path.join(a.results, "diag", "diag_sn2_clcl.csv")]:
        cp(f, os.path.join(D, "csv", "results", os.path.basename(f)))

    # ---- records
    clean_dir(mk(os.path.join(D, "records")), "*")
    for n in RECORDS:
        cp(os.path.join(a.results, n), os.path.join(D, "records", n))
    cp(os.path.join(a.results, "phase05_hcn_grid.npz"), os.path.join(D, "records", "phase05_hcn_grid.npz"))

    # ---- geometries (only the six reactions and CH4; the two extension systems are not part of the paper)
    gdir = mk(os.path.join(D, "geometries"))
    clean_dir(gdir, "*")
    for f in sorted(glob.glob(os.path.join(a.results, "geoms", "*.xyz"))):
        b = os.path.basename(f)
        if b.startswith(("sn2_fcl", "sn2_ohf")):
            continue
        cp(f, os.path.join(gdir, b))
    for rx in ("hcn", "nh3", "h2co", "ch4"):
        path_to_xyz(os.path.join(a.results, "geoms", rx + "_path.json"), os.path.join(gdir, rx + "_path_own_scan.xyz"))
    # path frames for the demo: keep the full JSON of the four paths (atoms and coordinates only)
    paths = {}
    for rx in ("hcn", "nh3", "h2co"):
        d = json.load(open(os.path.join(a.results, "geoms", rx + "_path.json"), encoding="utf-8"))
        paths[rx] = {"label": d["label"], "basis": d["basis"], "level": d["level"],
                     "points": [{"index": p["index"], "coordinate": p["coordinate"], "atoms": p["atoms"], "converged": p["converged"]} for p in d["points"]]}
    with open(os.path.join(D, "records", "path_frames.json"), "w") as f:
        json.dump(paths, f, separators=(",", ":"))

    # ---- check of the T1 grid against the stored array
    try:
        import numpy as np
        z = np.load(os.path.join(a.results, "phase05_hcn_grid.npz"))
        rows = list(csv.DictReader(open(os.path.join(D, "csv", "figures", "figS5_hcn_t1_map_data.csv"))))
        t1 = z["t1"]; th = z["theta"]; rr = z["r"]
        bad = 0
        for r in rows:
            i = int(np.argmin(abs(th - float(r["theta"])))); j = int(np.argmin(abs(rr - float(r["r"]))))
            if abs(t1[i, j] - float(r["t1"])) > 1e-4:
                bad += 1
        print("T1 grid: %d cells, %d differ from the .npz array" % (len(rows), bad))
        assert bad == 0
    except ImportError:
        print("numpy missing: T1 grid not cross-checked")

    # ---- figures
    for k in ("png", "pdf", "svg"):
        clean_dir(mk(os.path.join(SITE, "figures", k)), "*")
    names = [n for n, _ in FIG_ORDER] + SI_FIGS_names(a.writeup)
    for n in names:
        for k in ("png", "pdf", "svg"):
            src = os.path.join(a.writeup, "figures", n + "." + k)
            if os.path.exists(src):
                copy_figure(src, os.path.join(SITE, "figures", k, n + "." + k), k)
    ga = "graphical_abstract"
    for k in ("png", "svg"):
        src = os.path.join(a.writeup, "figures", ga + "." + k)
        if os.path.exists(src):
            copy_figure(src, os.path.join(SITE, "figures", k, ga + "." + k), k)
    main_caps = parse_captions(open(os.path.join(a.manuscript, "manuscript_B1_v5_checked.md"), encoding="utf-8").read().split("## Figure captions")[1].split("## References")[0], "main")
    si_text = open(os.path.join(a.manuscript, "SI_B1_v5_checked.md"), encoding="utf-8").read().split("## Supplementary figure captions")[1].split("## Supplementary tables")[0]
    si_caps = parse_captions(si_text, "si")
    figs = []
    for n, num in FIG_ORDER:
        c = main_caps[str(num)]
        assert c["src"] == n, (num, c["src"], n)
        figs.append({"id": n, "label": "Figure %d" % num, "title": c["title"], "caption": c["text"], "si": False})
    for i, n in enumerate(SI_FIGS):
        c = si_caps["S%d" % (i + 1)]
        figs.append({"id": SI_names(a.writeup)[i], "label": "Figure S%d" % (i + 1), "title": c["title"], "caption": c["text"], "si": True})
    with open(os.path.join(mk(os.path.join(D, "json")), "figure_captions.json"), "w", encoding="utf-8") as f:
        json.dump(figs, f, ensure_ascii=False, indent=1)

    # ---- notebooks and modules
    nd = mk(os.path.join(SITE, "code", "notebooks"))
    clean_dir(nd, "*.ipynb")
    for n in NOTEBOOKS:
        strip_notebook(os.path.join(a.notebooks, n + ".ipynb"), os.path.join(nd, n + ".ipynb"))
    md_ = mk(os.path.join(SITE, "code", "modules"))
    clean_dir(md_, "*.py")
    for n in MODULES:
        cp(os.path.join(a.notebooks, n + ".py"), os.path.join(md_, n + ".py"))
    gd = mk(os.path.join(SITE, "code", "geompack"))
    for n in ("b1v5_geompack_1h.json", "b1v5_geompack_1h.json.sha256", "GEOMPACK_README.md"):
        cp(os.path.join(a.notebooks, "geompack", n), os.path.join(gd, n))
    print("sync done")


def copy_figure(src, dst, kind):
    if kind == "png":
        clean_png(src, dst)
    elif kind == "svg":
        clean_svg(src, dst)
    else:
        cp(src, dst)


def SI_FIGS_names(writeup):
    return SI_names(writeup)


def SI_names(writeup):
    out = []
    for i in range(1, 10):
        g = sorted(glob.glob(os.path.join(writeup, "figures", "figS%d_*.png" % i)))
        out.append(os.path.basename(g[0])[:-4])
    return out


if __name__ == "__main__":
    main()
