"""Builds B1v5_CPU_1h.ipynb: the one-hour CPU notebook of PLAN_1hour.md.
    python build_cpu_1h.py [--pack geompack/b1v5_geompack_1h.json] [--out B1v5_CPU_1h.ipynb]      (needs nbformat)
The shared modules are pasted VERBATIM into code cells (edit the files, not the notebook); the geometry pack is embedded (gzip + base64) with its SHA-256."""
import os, sys, json, hashlib, argparse
import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
read = lambda n: open(os.path.join(HERE, n)).read()
md = lambda s: nbf.v4.new_markdown_cell(s.strip())
code = lambda s: nbf.v4.new_code_cell(s.strip("\n"))
banner = lambda n: f"# ===== {n} (pasted verbatim from the file of the same name; edit the file, not this cell) =====\n"

# estimated minutes on Colab (2 vCPU): sandbox measurement x 2, see README_1h.md section 7 (filled in by the builder from timings_1h.json when present)
EST = json.load(open(os.path.join(HERE, "timings_1h.json"))).get("cpu_colab_min", {}) if os.path.exists(os.path.join(HERE, "timings_1h.json")) else {}
e = lambda k, d: EST.get(k, d)

INTRO = f"""
# Paper B1 v5, CPU notebook, ONE HOUR (reference side, six reactions)

**What this notebook does.** The reference numbers for the paper, in a deadline of 50 minutes on a plain Colab CPU runtime (2 vCPU, no GPU): it checks the shipped **geometry pack** (SHA-256 and a gradient test), then CCSD(T) (frozen core, T1 diagnostic, RHF stability) at the stationary points of HCN, NH3, H2CO, CH4, CH3NC and F- + CH3F, along the stored paths of HCN, NH3, H2CO and CH4, at the GFN2-xTB scout geometries, a basis correction by density-fitted MP2 (QZVPP minus TZVP), one all-electron CCSD(T) barrier, then (tail, if time remains) MP2 along the paths, CCSD(T) for Cl- + CH3Cl and the HCN grid, and finally the automated check. Skala is **not** used here (that is the GPU notebook). Run this notebook FIRST, then the GPU notebook. No third pass.

**How to run it.** *Runtime > Run all*, click **Allow** in the Google Drive popup of the first cell (the only click), wait about **{e('core', 28)} minutes for the core and up to 50 minutes with the tail** (estimate for 2 vCPU: sandbox time x 2). Nothing else asks for input. At the end the runtime is released (billing stops); set `AUTO_DISCONNECT = False` in the settings cell to keep it.

**The deadline queue.** One ordered list of calculations. Before every calculation the notebook asks: *will it finish before minute 50?* (elapsed + 1.5 x the running mean time of that kind of calculation must be below 50 min, and the phase must stay inside its cap). A calculation that cannot finish is **not started** and is listed as "not run (time)" in `RUN_SUMMARY.md`. At minute 50 nothing new starts; check, summary and release take the last minutes. A second *Run all* continues with what was not run.

**Geometries.** The optimisations, transition states (one imaginary mode), Hessians, zero-point energies, relaxed scans and xTB scout geometries were computed once in the sandbox with the helper functions of the full notebooks (DF-B3LYP-D3(BJ)/def2-TZVP, def2-TZVPD for the SN2 species; GFN2-xTB) and are stored inside this notebook as a pack with a SHA-256. Cell 5 checks the hash, writes the geometries to the Drive folder `MyDrive/skala-b1-v5/` and cell C0 re-checks them with a B3LYP gradient (RMS below 3e-4 hartree/bohr). `REBUILD_GEOMS = True` (settings cell) rebuilds them here with the same code instead (measured about 1.9 hours in the 4-core sandbox, so about 4 hours on Colab; resumable; not part of the hour).

| Phase | Computes | Est. min | Cap |
|---|---|---|---|
| install, Drive | `pyscf==2.15.0`, `pyscf-dispersion==1.7.0` | {e('install', 3)} | 6 |
| C0 | pack hashes, speed calibration, B3LYP gradient at 12 structures, one stored v4 CCSD(T) energy | {e('C0', 2.3)} | 6 |
| C3a | CCSD(T) at the stationary points, HCN to F- + CH3F | {e('C3a', 11.4)} | 14 |
| C3b | CCSD(T) along the paths (40 points) | {e('C3b', 5.2)} | 7 |
| C3c | MP2 basis correction at the stationary points of HCN, NH3, H2CO, CH3NC, F- + CH3F | {e('C3c', 2.3)} | 12 |
| C4 | CCSD(T) at the xTB minimum and peak, four reactions | {e('C4', 2.9)} | 4 |
| C6 | all-electron CCSD(T), HCN | {e('C6', 0.4)} | 2 |
| tail 1 | MP2 at the 31 path points | {e('tail1', 3.2)} | to minute 50 |
| tail 2 | CCSD(T) + MP2 for Cl- + CH3Cl (Tier C to B) | {e('tail2', 10.8)} | to minute 50 |
| tail 3 | CCSD(T)/def2-SVP on the 425-point HCN grid | {e('tail3', 6.4)} | to minute 50 |
| C7 | check, summary, release | {e('C7', 1)} | 3 |

If the runtime disconnects: everything is saved after every calculation; run all again. Extension sessions (one click each, see the README): `EXTENSION = 'CPU1'` (the tail list again, for whatever the hour did not reach) and `EXTENSION = 'GEOM'` (rebuild the geometry pack here).
"""

CELL_SETUP = '''
# ---- Cell 1: the clock starts here (install time counts against the 50 minutes), install, shared helpers, Drive ----
import time as _t
_T0 = _t.time()
import sys, subprocess, importlib.util, os
try:   # Colab: ask for Google Drive FIRST, so that the Allow popup is the first thing you see. The 50-minute clock starts after the mount.
    from google.colab import drive as _drive
    print("=" * 90); print("CLICK 'Connect to Google Drive' and then Allow in the popup NOW. The 50-minute clock starts after you have clicked."); print("=" * 90)
    for _a in range(3):
        try:
            _drive.mount("/content/drive", force_remount=(_a > 0))
            if os.path.isdir("/content/drive/MyDrive"):
                break
        except Exception as _ex:
            print(f"Drive mount attempt {_a + 1} of 3 failed ({repr(_ex)[:80]}); trying again in 20 s: please click Allow.")
            _t.sleep(20)
    _T0 = _t.time()
except ImportError:
    pass
def _have(mod):
    try:
        return importlib.util.find_spec(mod) is not None
    except Exception:
        return False
_need = [s for m, s in (("pyscf", "pyscf==2.15.0"), ("pyscf.dispersion", "pyscf-dispersion==1.7.0")) if not _have(m)]       # the plan: only these two
if _need:
    print("installing:", _need)
    try:
        _r = subprocess.run([sys.executable, "-m", "pip", "install", "-q"] + _need, capture_output=True, text=True, timeout=6 * 60)    # install cap: 6 minutes
        print("pip:", "ok" if _r.returncode == 0 else "FAILED " + _r.stderr[-300:], f"({_t.time() - _T0:.0f} s)")
    except Exception as _ex:
        print("pip stopped at the 6-minute cap:", repr(_ex)[:100])
else:
    print("all packages already installed")

''' + banner("b1v5_common.py") + "{COMMON}" + '''

# mount Drive (the one click) or fall back to a local folder, make the folders
setup_root()
'''

CELL_HELPERS = banner("b1v5_cpu.py") + "{CPU}" + "\n\n" + banner("b1v5_v4_stored.py") + "{V4}" + "\n\n" + banner("b1v5_geompack.py") + "{PACK}"
CELL_CHECK_DEFS = banner("b1v5_lit_table.py") + "{LIT}" + "\n\n" + banner("b1v5_check.py") + "{CHECK}"
CELL_QUEUE = banner("b1v5_1h.py") + "{QUEUE}" + "\n\n" + banner("b1v5_1h_check.py") + "{CHECK1H}" + "\n\n" + banner("b1v5_1h_cpu.py") + "{CPU1H}"

CELL_PACKDATA = '''
# ---- The geometry pack: {PACKINFO}.  SHA-256 of the pack (canonical JSON): this constant is checked in the next cell, and published with the paper. ----
GEOMPACK_SHA256 = "{SHA}"
GEOMPACK_B64 = "{B64}"
'''

CELL_CONFIG = '''
# ---- Settings (all optional; the defaults are the plan) ----
import os, time, json
SMOKE = os.environ.get("B1V5_SMOKE") == "1"      # quick plumbing test (small bases, three systems); never set this on Colab
DEADLINE_MIN = float(os.environ.get("B1V5_DEADLINE_MIN", "50"))   # no calculation is started after this minute (the clock started in cell 1); 5 more minutes are reserved for check and summary
AUTO_DISCONNECT = True       # release the Colab runtime at the very end so billing stops
EXTENSION = os.environ.get("B1V5_EXTENSION") or None             # (the environment variable is for the tests) None = the one-hour core plus its tail.  'CPU1' = extension session 1 (the tail list again).  'GEOM' = rebuild the geometry pack here
REBUILD_GEOMS = False        # True (or EXTENSION = 'GEOM'): rebuild the geometries with the same code (about 4 hours on Colab, NOT part of the hour; resumable) instead of using the pack
RETRY_FAILED = False         # True: calculations that failed earlier are tried again
FORCE_RERUN = False          # True: ignore DONE markers (finished points are still skipped)
RUN_TAIL = True              # the tail list after the core (MP2 at the path points, HCN grid); it runs only if time remains
MP2_MAX_NBF = 450            # an MP2 basis-correction point with more basis functions is skipped (the reference stays uncorrected there, flagged)
if EXTENSION == "GEOM":
    REBUILD_GEOMS = True
RUN_CORE = EXTENSION is None
SYSTEMS_TO_RUN = None
GRID_STRIDE = (6, 4) if SMOKE else (1, 1)      # HCN grid of the tail: every point (smoke test: every 6th / 4th)
if SMOKE:
    BASIS_OVERRIDE.update({"small": "STO-3G", "main": "def2-SVP", "large": "def2-TZVP"})
os.environ["B1V5_ROOT"] = ROOT                                             # read by the very last cell
os.environ["B1V5_AUTO_DISCONNECT"] = "1" if AUTO_DISCONNECT else "0"
if EXTENSION == "GEOM" and "B1V5_DEADLINE_MIN" not in os.environ:
    DEADLINE_MIN = 150.0     # the rebuild of the geometry pack is long (about 4 hours on Colab); it is resumable, a second Run all continues
q_start(DEADLINE_MIN, t0=_T0)
set_design(DESIGN_REACTIONS)
SYSTEMS_TO_RUN = list(SYSTEM_ORDER)
if SMOKE:
    SYSTEMS_TO_RUN = ["hcn", "nh3", "ch4"]
quiet_geometric()
write_json("config.json", dict(design="1h", smoke=SMOKE, basis_override=dict(BASIS_OVERRIDE), started=now(), reactions=list(DESIGN["reactions"]), headline_family="own",
                               own_geometries_for_all=True, deadline_min=DEADLINE_MIN, extension=EXTENSION, rebuild_geoms=bool(REBUILD_GEOMS)))
print("Hardware:", hw()); print("Versions:", ver()); print("Design reactions:", DESIGN["reactions"], "| systems to run:", SYSTEMS_TO_RUN)
print(f"Deadline: minute {DEADLINE_MIN:.0f} (now minute {q_now_min():.1f}); extension: {EXTENSION}; rebuild geometries: {REBUILD_GEOMS}")
if SMOKE:
    print("*** SMOKE TEST: small bases, three systems. These numbers are NOT the paper numbers. ***")
'''

MD_PACK = "**Cell 5: install the geometry pack** (SHA-256 of the pack and of every structure, then geoms/*.json and the rows of phase01 / phase02 / phase04). With `REBUILD_GEOMS = True` the geometries are rebuilt here instead. Seconds."
CELL_PACK = '''
PACK_INFO = {}
if not REBUILD_GEOMS:
    PACK_INFO = pack_setup(GEOMPACK_B64, GEOMPACK_SHA256, systems=SYSTEMS_TO_RUN)
    print(f"geometry pack: SHA-256 {PACK_INFO.get('sha256')}\\n  expected      {GEOMPACK_SHA256}\\n  {'OK: the hash matches and every structure verified' if PACK_INFO.get('ok') else 'PROBLEM: ' + str(PACK_INFO.get('problems'))}")
    print("  installed:", PACK_INFO.get("installed"))
    if not PACK_INFO.get("ok"):
        print("*" * 100); print("The geometry pack did not verify and was NOT installed. Later phases that need geometries will be skipped and say so."); print("*" * 100)
else:
    pip_install_capped(["geometric==1.1.1", "tblite==0.8.0", "ase==3.29.0"], 10)
    print("REBUILD_GEOMS: rebuilding the geometries with the full-size code (not part of the hour; resumable; a rerun continues).")
    t0_ = time.time()
    with phase("geompack_rebuild", "rebuild the geometries", max_min=600) as ph:
        if not ph.skip:
            PACK_INFO = dict(ok=True, rebuilt=rebuild_geoms(SYSTEMS_TO_RUN))
    print(f"rebuild: {(time.time() - t0_) / 60:.1f} min", PACK_INFO)
if EXTENSION == "GEOM":
    print("EXTENSION = 'GEOM': geometries rebuilt; the other cells find their phases done or run the core on the new geometries (set EXTENSION = None for that).")
'''

MD_C0 = "**C0: verify the geometry pack.** SHA-256 was checked in cell 5; here a B3LYP-D3(BJ) energy and gradient is computed at the 10 unimolecular structures and the 2 SN2 transition states (pass: RMS gradient below 3e-4 hartree/bohr) and one stored v4 CCSD(T) energy is rerun (check I9). The speed of this machine against the sandbox is measured once and scales the time estimates of the queue. Cap 6 minutes."
CELL_C0 = '''
if RUN_CORE:
    with qphase("phase00", "C0 verify geometry pack", CPU_CAPS["phase00"]) as ph:
      if not ph.skip:
        c0_verify(PACK_INFO)
'''
MD_C3A = "**C3a: CCSD(T) at the stationary points** (frozen core, def2-TZVP, def2-TZVPD for F- + CH3F; T1 diagnostic and RHF stability). Order: HCN, NH3, H2CO, CH4, CH3NC, then F- + CH3F (four species, the longest; it becomes Tier C, literature only, if the queue refuses it). Cap 14 minutes."
CELL_C3A = '''
if RUN_CORE:
    with qphase("phase03a", "C3a CCSD(T) stationary points", CPU_CAPS["phase03a"]) as ph:
      if not ph.skip:
        c3a_stationary()
'''
MD_C3B = "**C3b: CCSD(T) along the stored paths** HCN (11 points), NH3 (9), H2CO (11), CH4 (9 stretch points). Cap 7 minutes."
CELL_C3B = '''
if RUN_CORE:
    with qphase("phase03b", "C3b CCSD(T) along the paths", CPU_CAPS["phase03b"]) as ph:
      if not ph.skip:
        c3b_paths()
'''
MD_C3C = "**C3c: MP2 basis correction** of the CCSD(T) reference: DF-MP2/def2-QZVPP minus DF-MP2/def2-TZVP (frozen core) at the stationary points of HCN, NH3, H2CO and CH3NC (12 structures). If it runs slower than expected the cap cuts CH3NC first; the reference stays raw there (check I15 says so). Cap 12 minutes."
CELL_C3C = '''
if RUN_CORE:
    with qphase("phase08", "C3c MP2 basis correction (stationary points)", CPU_CAPS["phase08"]) as ph:
      if not ph.skip:
        c3c_mp2_stationary()
'''
MD_C4 = "**C4: scout reference.** CCSD(T) at the GFN2-xTB minimum and peak of HCN, NH3, H2CO and CH3NC (the geometries come from the pack). Cap 4 minutes."
CELL_C4 = '''
if RUN_CORE:
    with qphase("phase04", "C4 scout reference", CPU_CAPS["phase04"]) as ph:
      if not ph.skip:
        c4_scout_reference()
'''
MD_C6 = "**C6: all-electron CCSD(T)** barrier for HCN only (check I10). Cap 2 minutes."
CELL_C6 = '''
if RUN_CORE:
    with qphase("phase06", "C6 all-electron CCSD(T), HCN", CPU_CAPS["phase06"]) as ph:
      if not ph.skip:
        c6_all_electron()
'''
MD_TAIL = "**Tail list** (only if time remains, same admission test; ends at minute 50), in this order: (1) MP2 basis correction at the 31 path points of HCN, NH3 and H2CO; (2) CCSD(T) for Cl- + CH3Cl (four species, def2-TZVPD) with its MP2 correction: this moves the reaction from Tier C to Tier B; (3) CCSD(T)/def2-SVP on the 425-point rigid HCN grid. Whatever is not reached is done by a second Run all with `EXTENSION = 'CPU1'` (the same list; finished items are skipped)."
CELL_TAIL = '''
if (RUN_CORE and RUN_TAIL) or EXTENSION == "CPU1":
    with qphase("phase08t", "tail 1: MP2 basis correction at the path points", 45 if EXTENSION != "CPU1" else 12) as ph:
      if not ph.skip:
        tail_mp2_paths()
    with qphase("phase03x", "tail 2: CCSD(T) and MP2 for Cl- + CH3Cl", 45 if EXTENSION != "CPU1" else 25) as ph:
      if not ph.skip:
        ext_clcl_ccsdt()
    with qphase("phase05", "tail 3: CCSD(T)/def2-SVP on the HCN grid", 45 if EXTENSION != "CPU1" else 12) as ph:
      if not ph.skip:
        tail_hcn_grid(GRID_STRIDE)
'''
MD_C7 = "**C7: automated check.** Compares everything with the literature table and runs the internal checks I1 to I15, with the six reactions of this design. Items that cannot exist in the one-hour design are printed as NOT RUN (never FAIL). Always runs, whatever the clock says. About 2 minutes."
CELL_C7 = '''
with qphase("phase07x_check", "C7 automated check", CPU_CAPS["check"], final=True, force=True) as ph:
    check_result = run_check()
'''
MD_SUM = "**Summary and release.** Writes `RUN_SUMMARY.md` (what ran, what failed, what was not run for time, minutes against the caps) to the Drive folder."
CELL_SUM = '''
try:
    q_save_log()
    _pk = [f"Geometry pack: SHA-256 {PACK_INFO.get('sha256')} (expected {GEOMPACK_SHA256}), verified: {PACK_INFO.get('ok')}, problems: {PACK_INFO.get('problems')}",
           f"Gradient test (C0): {grad_summary()}", f"Extension: {EXTENSION}; rebuild: {REBUILD_GEOMS}; deadline {DEADLINE_MIN} min"]
    print(write_run_summary(extra_lines=["", "## Geometry pack", ""] + _pk + queue_lines()))
except Exception as ex:
    print("could not write RUN_SUMMARY.md:", repr(ex)[:150])
'''
MD_END = "**Last cell: release the runtime.** Independent of everything above: waits a few seconds so Drive can finish writing, then releases the Colab runtime (stops billing) unless `AUTO_DISCONNECT = False`. It runs even if earlier cells failed (use *Runtime > Run all*)."
CELL_END = '''
# This cell does not use any variable from the cells above.
import os, time, json
print("Final cell: giving Drive a few seconds to finish writing ...")
time.sleep(8)
try:   # a short fallback note, in case the summary cell above did not run
    root = os.environ.get("B1V5_ROOT") or "/content/drive/MyDrive/skala-b1-v5"
    if os.path.isdir(root) and not os.path.exists(os.path.join(root, "RUN_SUMMARY.md")):
        st = json.load(open(os.path.join(root, "phase_status.json"))) if os.path.exists(os.path.join(root, "phase_status.json")) else {}
        open(os.path.join(root, "RUN_SUMMARY.md"), "w").write("# Run summary (fallback)\\n\\n" + "\\n".join(f"- {k}: {v.get('status')}" for k, v in sorted(st.items())) + "\\n")
except Exception as ex:
    print("fallback summary not written:", repr(ex)[:100])
_drive_ok = bool(os.environ.get("B1V5_ROOT")) or os.path.isdir("/content/drive/MyDrive") or not os.path.isdir("/content")
if not _drive_ok:      # Google Drive was not mounted: the results are only in the local folder of this runtime. Do NOT release the runtime.
    try:
        import shutil
        _loc = "/content/skala-b1-v5-local"
        _zip = shutil.make_archive("/content/skala-b1-v5-RESULTS", "zip", _loc) if os.path.isdir(_loc) else None
        print("*" * 90); print("GOOGLE DRIVE WAS NOT MOUNTED. The results are in", _loc, "and in", _zip)
        print("The runtime is NOT released, so the results are not lost. Download the zip (the download starts below, or use the Files panel on the left),")
        print("then disconnect the runtime yourself (Runtime > Disconnect and delete runtime)."); print("*" * 90)
        if _zip:
            from google.colab import files
            files.download(_zip)
    except Exception as ex:
        print("zip or download did not work (", repr(ex)[:100], "): use the Files panel on the left to download /content/skala-b1-v5-local.")
elif os.environ.get("B1V5_AUTO_DISCONNECT", "1") == "1":
    try:
        from google.colab import runtime
        print("Releasing the Colab runtime now (this stops billing).")
        runtime.unassign()
    except Exception as ex:
        print("Not on Colab, or the runtime could not be released (", repr(ex)[:80], "). Nothing else to do.")
else:
    print("AUTO_DISCONNECT is off: the runtime stays connected. Disconnect it yourself to stop billing.")
'''


def build(out_name="B1v5_CPU_1h.ipynb", pack_path=None):
    import b1v5_geompack as G
    pack_path = pack_path or os.path.join(HERE, "geompack", "b1v5_geompack_1h.json")
    raw = open(pack_path, "rb").read()
    pack = json.loads(raw.decode())
    assert G.pack_bytes(pack) == raw, "the pack file is not in canonical form (rebuild it with build_geompack.py)"
    sha = hashlib.sha256(raw).hexdigest()
    b64 = G.pack_encode(pack)
    info = f"{len(pack['structures'])} structures, {len(pack['paths'])} paths, {len(pack['xtb'])} xTB sets, {len(raw) // 1024} KB JSON, {len(b64) // 1024} KB embedded; systems {', '.join(pack['systems'])}; built {pack['created'][:10]} with {pack['builder']['versions'].split(';')[0]}"
    sub = lambda s: s.replace("{COMMON}", read("b1v5_common.py")).replace("{CPU}", read("b1v5_cpu.py")).replace("{V4}", read("b1v5_v4_stored.py")).replace("{PACK}", read("b1v5_geompack.py")) \
        .replace("{LIT}", read("b1v5_lit_table.py")).replace("{CHECK}", read("b1v5_check.py")).replace("{QUEUE}", read("b1v5_1h.py")).replace("{CHECK1H}", read("b1v5_1h_check.py")).replace("{CPU1H}", read("b1v5_1h_cpu.py"))
    cells = [md(INTRO),
             md("**Cell 1: the clock, install, shared helpers, Drive.** The clock of the 50-minute deadline starts in the first line. Installs `pyscf==2.15.0` and `pyscf-dispersion==1.7.0` (cap 6 minutes), defines the shared helper functions and mounts Google Drive (one permission click)."), code(sub(CELL_SETUP)),
             md("**Cell 2: calculation helpers** (PySCF, geomeTRIC, xTB helpers of the full notebooks, the stored v4 numbers, and the geometry-pack tools). Only defines functions; a few seconds."), code(sub(CELL_HELPERS)),
             md("**Cell 3: check module** (literature table and the automated check, pasted unchanged). Only defines functions."), code(sub(CELL_CHECK_DEFS)),
             md("**Cell 4: deadline queue, six-reaction design and CPU phases.** Only defines functions and constants."), code(sub(CELL_QUEUE)),
             md("**Cell 5a: the geometry pack itself** (data, with its SHA-256)."), code(CELL_PACKDATA.replace("{PACKINFO}", info).replace("{SHA}", sha).replace("{B64}", b64)),
             md("**Cell 6: settings.** Deadline, extension sessions, rebuild flag, release of the runtime. Starts the queue clock."), code(CELL_CONFIG),
             md(MD_PACK), code(CELL_PACK),
             md(MD_C0), code(CELL_C0), md(MD_C3A), code(CELL_C3A), md(MD_C3B), code(CELL_C3B), md(MD_C3C), code(CELL_C3C), md(MD_C4), code(CELL_C4), md(MD_C6), code(CELL_C6),
             md(MD_TAIL), code(CELL_TAIL), md(MD_C7), code(CELL_C7), md(MD_SUM), code(CELL_SUM), md(MD_END), code(CELL_END)]
    nb = nbf.v4.new_notebook(cells=cells)
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    nb.metadata["colab"] = {"provenance": [], "name": out_name}
    nbf.validate(nb)
    nbf.write(nb, os.path.join(HERE, out_name) if not os.path.isabs(out_name) else out_name)
    return out_name, sha


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pack", default=None)
    ap.add_argument("--out", default="B1v5_CPU_1h.ipynb")
    a = ap.parse_args()
    p, sha = build(a.out, a.pack)
    full = p if os.path.isabs(p) else os.path.join(HERE, p)
    print("wrote", full, os.path.getsize(full) // 1024, "KB; pack SHA-256", sha)
