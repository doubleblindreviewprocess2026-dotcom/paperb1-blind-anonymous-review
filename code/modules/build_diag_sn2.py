"""Builds B1v5_GPU_diag_sn2.ipynb: the Cl- + CH3Cl Skala diagnostic (about 25 minutes on a T4).
    python build_diag_sn2.py [--out B1v5_GPU_diag_sn2.ipynb] [--pack ...]
Uses the cells 1 to 6 of the GPU notebook (install, helpers, pack, settings, G0), then the diagnostic phase, the report and the last cell."""
import os, sys, json, hashlib, argparse
import nbformat as nbf
import build_gpu_1h as B
from build_gpu_1h import md, code, banner, read, HERE

INTRO = """# B1 v5: diagnostic of the Cl- + CH3Cl transition state with Skala-1.1

**What this is.** The full run gave a Skala error of -3.1 kcal/mol for Cl- + CH3Cl, but Skala converged at the transition state only through the density-fitting rescue, and at the xTB geometry the same barrier came out at 16.0 (wB97M-V: 13.5 and 13.7 at both). This notebook finds out whether the -3.1 is real: restarts from different initial guesses, with and without density fitting, a profile along the C-Cl stretch, and the same profile with wB97M-V.

**How to run.** GPU runtime (T4 is enough), **Run all**, click **Allow** in the Drive popup at the start. About 25 minutes. Results go to the same Drive folder `skala-b1-v5`: `diag_sn2_clcl.csv`, `DIAG_SN2_SUMMARY.md`, `figures/diag_sn2_stretch.png`. The files of the main run are not changed (`RUN_SUMMARY.md` is NOT overwritten). The runtime is released at the end if Drive was mounted.
"""

CELL_DIAG = '''
run_gphase_1h("phase17", "SN2 diagnostics (Cl- + CH3Cl)", run_sn2_diagnostics, GPU_CAPS["phase17"])
'''
CELL_REPORT = '''
try:
    q_save_log()
    report_sn2_diagnostics()
except Exception as ex:
    print("could not write the diagnostic summary:", repr(ex)[:200])
'''


def build(out_name="B1v5_GPU_diag_sn2.ipynb", pack_path=None):
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
             md("**Diagnostic phase.** Skala (density fitting and not), wB97M-V, at the Cl- + CH3Cl RC and TS of the own and the xTB geometries and along the C-Cl stretch. Cap 30 minutes."), code(CELL_DIAG),
             md("**Report.** Writes `DIAG_SN2_SUMMARY.md` and the profile figure to Drive and prints them."), code(CELL_REPORT),
             md(B.MD_END), code(B.CELL_END)]
    nb = nbf.v4.new_notebook(cells=cells)
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    nb.metadata["colab"] = {"provenance": [], "name": out_name, "gpuType": "T4", "machine_shape": "hm"}
    nb.metadata["accelerator"] = "GPU"
    nbf.write(nb, out_name)
    return out_name, sha


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--pack", default=None); ap.add_argument("--out", default="B1v5_GPU_diag_sn2.ipynb"); a = ap.parse_args()
    p, sha = build(a.out, a.pack)
    print("wrote", p, os.path.getsize(p) // 1024, "KB; pack", sha[:12])
