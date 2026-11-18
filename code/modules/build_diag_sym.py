"""Builds B1v5_GPU_diag_sym.ipynb: the Cl- + CH3Cl Skala diagnostic (about 25 minutes on a T4).
    python build_diag_sym.py [--out B1v5_GPU_diag_sym.ipynb] [--pack ...]
Uses the cells 1 to 6 of the GPU notebook (install, helpers, pack, settings, G0), then the diagnostic phase, the report and the last cell."""
import os, sys, json, hashlib, argparse
import nbformat as nbf
import build_gpu_1h as B
from build_gpu_1h import md, code, banner, read, HERE

INTRO = """# B1 v5: where is the Skala transition state of Cl- + CH3Cl and F- + CH3F?

**What this is.** The first diagnostic showed that the Skala energy at the B3LYP transition state of Cl- + CH3Cl is a reliable SCF state (density fitting changes it by 0.05 kcal/mol), but my scan was along the wrong coordinate: for a D3h SN2 transition state the symmetric C-Cl stretch has the TS as a **minimum**, not a maximum (the reaction coordinate is the antisymmetric stretch). The profile of Skala was still falling at C-Cl = 2.45 A, while wB97M-V has its minimum near 2.35 A. So the Skala transition state is probably looser, and its barrier lower than the -3.1 error of the main run shows.

This notebook scans the symmetric stretch until the minimum is inside the range, for Cl- + CH3Cl and F- + CH3F, with Skala-1.1 and wB97M-V. The minimum gives the barrier of each method at its own transition state (rigid CH3).

**How to run.** GPU runtime (T4 is enough), **Run all**, click **Allow** in the Drive popup at once. About 20 to 30 minutes. Points already computed by the first diagnostic are reused from the Drive folder `skala-b1-v5`. Results: `diag_sn2_clcl.csv` (extended), `DIAG_SYM_SUMMARY.md`, `figures/diag_sym_stretch.png`. The files of the main run are not changed. The runtime is released at the end if Drive was mounted.
"""

CELL_DIAG = '''
run_gphase_1h("phase18", "SN2 symmetric-stretch scan", run_symscan, GPU_CAPS["phase18"])
'''
CELL_REPORT = '''
try:
    q_save_log()
    report_symscan()
except Exception as ex:
    print("could not write the diagnostic summary:", repr(ex)[:200])
'''


def build(out_name="B1v5_GPU_diag_sym.ipynb", pack_path=None):
    import b1v5_geompack as G
    pack_path = pack_path or os.path.join(HERE, "geompack", "b1v5_geompack_1h.json")
    raw = open(pack_path, "rb").read()
    pack = json.loads(raw.decode())
    sha = hashlib.sha256(raw).hexdigest()
    b64 = G.pack_encode(pack)
    info = f"{len(pack['structures'])} structures (pack of the main run)"
    sub = lambda s: s.replace("{COMMON}", read("b1v5_common.py")).replace("{CPU}", read("b1v5_cpu.py")).replace("{V4}", read("b1v5_v4_stored.py")).replace("{PACK}", read("b1v5_geompack.py")) \
        .replace("{GPU}", read("b1v5_gpu.py")).replace("{LIT}", read("b1v5_lit_table.py")).replace("{CHECK}", read("b1v5_check.py")).replace("{FIGS}", read("b1v5_gpu_figs.py")) \
        .replace("{QUEUE}", read("b1v5_1h.py")).replace("{CHECK1H}", read("b1v5_1h_check.py")).replace("{GPU1H}", read("b1v5_1h_gpu.py"))
    config = B.CELL_CONFIG.replace('os.environ.get("B1V5_DEADLINE_MIN", "50")', 'os.environ.get("B1V5_DEADLINE_MIN", "40")')
    assert config != B.CELL_CONFIG
    cells = [md(INTRO),
             md("**Cell 1: Drive, install, shared helpers.**"), code(sub(B.CELL_SETUP)),
             md("**Cell 2: calculation helpers.** Only defines functions."), code(sub(B.CELL_HELPERS)),
             md("**Cell 3: check module and figure code** (unchanged, needed by the queue code). Only defines functions."), code(sub(B.CELL_CHECK_DEFS)),
             md("**Cell 4: deadline queue and GPU phases.** Only defines functions."), code(sub(B.CELL_QUEUE)),
             md("**Cell 5: the diagnostic code.** Only defines functions."), code(read("b1v5_diag_sn2.py")),
             md("**Cell 6: the geometry pack** (data, with its SHA-256)."), code(B.CELL_PACKDATA.replace("{PACKINFO}", info).replace("{SHA}", sha).replace("{B64}", b64)),
             md("**Cell 7: settings.** Deadline minute 40."), code(config),
             md(B.MD_PACK), code(B.CELL_PACK),
             md(B.MD_G0), code(B.CELL_G0),
             md("**Scan phase.** Skala and wB97M-V at the RC and along the symmetric X-C stretch of the own TS, for Cl- + CH3Cl and F- + CH3F. Cap 40 minutes."), code(CELL_DIAG),
             md("**Report.** Writes `DIAG_SYM_SUMMARY.md` and the profile figure to Drive and prints them."), code(CELL_REPORT),
             md(B.MD_END), code(B.CELL_END)]
    nb = nbf.v4.new_notebook(cells=cells)
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    nb.metadata["colab"] = {"provenance": [], "name": out_name, "gpuType": "T4", "machine_shape": "hm"}
    nb.metadata["accelerator"] = "GPU"
    nbf.write(nb, out_name)
    return out_name, sha


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--pack", default=None); ap.add_argument("--out", default="B1v5_GPU_diag_sym.ipynb"); a = ap.parse_args()
    p, sha = build(a.out, a.pack)
    print("wrote", p, os.path.getsize(p) // 1024, "KB; pack", sha[:12])
