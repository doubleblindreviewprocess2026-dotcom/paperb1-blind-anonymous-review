"""Builds B1v5_GPU_1h.ipynb: the one-hour GPU notebook of PLAN_1hour.md.
    python build_gpu_1h.py [--pack geompack/b1v5_geompack_1h.json] [--out B1v5_GPU_1h.ipynb]      (needs nbformat)
The shared modules are pasted VERBATIM into code cells (edit the files, not the notebook); the geometry pack is embedded (gzip + base64) with its SHA-256."""
import os, sys, json, hashlib, argparse
import nbformat as nbf

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
read = lambda n: open(os.path.join(HERE, n)).read()
md = lambda s: nbf.v4.new_markdown_cell(s.strip())
code = lambda s: nbf.v4.new_code_cell(s.strip("\n"))
banner = lambda n: f"# ===== {n} (pasted verbatim from the file of the same name; edit the file, not this cell) =====\n"

INTRO = """
# Paper B1 v5, GPU notebook, ONE HOUR (Skala-1.1 and four standard functionals, six reactions)

**What this notebook does.** The GPU side of the paper in a deadline of 50 minutes on a T4 (L4 and A100 reach further down the same list): Skala-1.1 (with its D3 term, headline) and its "without D3" energy from the same SCF, wB97M-V, B3LYP-D3(BJ), r2SCAN and PBE at every stationary point and fragment of six reactions (HCN, NH3, H2CO, CH3NC, F- + CH3F, Cl- + CH3Cl), Skala, wB97M-V and B3LYP-D3(BJ) along the relaxed paths of HCN, NH3, H2CO and the CH4 stretch (40 points), the same three at the GFN2-xTB scout geometries, a basis check with def2-QZVPP for HCN, NH3, H2CO (Skala and wB97M-V), and finally the figures, the tables, the frozen-criteria outcome table and the automated check on whatever exists. Everything that does not fit is listed as NOT RUN (never FAIL).

**How to run it.** Run the **CPU notebook first**, then this one (order CPU, GPU; no third pass). *Runtime > Change runtime type > GPU* (T4 is enough; L4 or A100 recommended), then *Runtime > Run all*. The only click is the Google Drive permission popup in the first cell. The geometries come from the same shipped, hashed **geometry pack** as the CPU notebook (no geometry work here). At the end the runtime is released (set `AUTO_DISCONNECT = False` to keep it).

**The deadline queue.** One ordered list. Before every calculation: *will it finish before minute 50?* (elapsed + 1.5 x the running mean of that kind of calculation below 50 minutes, and the phase inside its cap). If not, it is not started and is listed as "not run (time)" in `RUN_SUMMARY.md`. After the core comes a tail list (T4: expected to be empty; L4: reactions 7 and 8; A100: items 1 to 4, perhaps 5). Five minutes are reserved for figures, check and summary, which always run. A second *Run all* continues with what was not run.

**Skala with D3.** The headline method `Skala-1.1` is the published functional with its D3 term (package default); `Skala-1.1-noD3` is the same SCF minus the stored D3 energy. D3 is switched off only in the G0 reproduction of the stored v4 energies, which were computed without it.

**If something is not available.** No GPU: only the G0 diagnostics run and the notebook says so. Skala weights cannot be downloaded within the cap: the Skala rows are missing, the other functionals run, the check says NOT RUN. `STANDIN = True` (settings cell) is the dry run on the CPU with cheap functionals in tiny bases (numbers are not results).

| Phase | Computes | Est. min (T4) | Cap |
|---|---|---|---|
| install, weights | pinned pyscf, Skala wheel, gpu4pyscf, Skala weights (Drive cache) | 8 | 15 |
| pack | hash check, geometries to Drive | under 1 | |
| G0 | functional probe, Skala check, stored v4 HCN point 0 (gate, 1e-5 hartree) | 2 | 4 |
| G1a | 20 species x 5 functionals, def2-TZVP(D) | 13 | 18 |
| G1b | 40 path points x 3 functionals | 11 | 12 |
| G3 | scout check: 36 energies | 5 | 6 |
| G2 | def2-QZVPP, 9 species, Skala and wB97M-V | 8 | 8 |
| tail | by GPU class (see above) | T4 0 | to minute 50 |
| G5, check | figures, tables, outcome table, automated check | 3 | 5 |

Extension sessions (one click each, any order, see the README): `EXTENSION = 'GPU1'` (tail 1 to 3: reactions 7 and 8, PBE and r2SCAN on paths, QZVPP(D) for CH3NC and the SN2 systems) and `EXTENSION = 'GPU2'` (tail 4 and 5: PBE and r2SCAN at the scout geometries, Skala-relaxed HCN, CH3NC, F- + CH3F). `EXTENSION = 'GEOM'` is a flag of the CPU notebook (rebuild of the geometry pack).
"""

CELL_SETUP = '''
# ---- Cell 1: the clock starts here (install time counts against the 50 minutes), install with a 15-minute cap, shared helpers, Drive ----
import time as _t
_T0 = _t.time()
import sys, subprocess, importlib.util, shutil, os
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
SKALA_VERSION = "2026.9"      # the version verified on Colab (2026.4 gave wrong forces)
PYSCF_VERSION = "2.12.1"      # PINNED: skala 2026.9 needs pyscf >= 2.11 and < 2.15 (wheel metadata), gpu4pyscf 1.8.1 needs >= 2.8.  Not yet run on Colab.
DISPERSION_VERSION = "1.7.0"  # pyscf-dispersion (the D3 library Skala and B3LYP-D3(BJ) use); the CPU notebook pins the same
INSTALL_CAP_MIN = 15          # install + weights together may not use more than this
def _left_s():
    return max(20.0, 60.0 * INSTALL_CAP_MIN - (_t.time() - _T0))
def _pip(*a):
    try:
        r = subprocess.run([sys.executable, "-m", "pip", "install", "-q", *a], capture_output=True, text=True, timeout=_left_s())
    except Exception as ex:
        print("FAIL (cap or error) " + " ".join(a), repr(ex)[:80]); return False
    print(("ok   " if r.returncode == 0 else "FAIL ") + " ".join(a) + f"  ({_t.time() - _T0:.0f} s since start)")
    if r.returncode:
        print(r.stderr[-500:])
    return r.returncode == 0
def _have(mod):
    try:
        return importlib.util.find_spec(mod) is not None
    except Exception:
        return False
def _gpu_present():
    try:
        out = subprocess.run(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"], capture_output=True, text=True, timeout=15).stdout.strip()
        return bool(out)
    except Exception:
        return False
GPU_PRESENT = _gpu_present()
try:
    if GPU_PRESENT and not (_have("skala") and _have("gpu4pyscf")):
        _pip(f"pyscf=={PYSCF_VERSION}", f"pyscf-dispersion=={DISPERSION_VERSION}")                  # pyscf first, pinned, so that no later package pulls another version
        if not _pip(f"skala-cuda12x=={SKALA_VERSION}"):                                             # GPU wheel for CUDA 12 (the Colab default); needs gpu4pyscf-cuda12x >= 1.8, < 1.9
            _pip(f"skala=={SKALA_VERSION}")
        if not _have("gpu4pyscf"):
            _pip("gpu4pyscf-cuda12x==1.8.1")                                                         # pinned; NO unpinned fallback.  If this fails, the GPU route is off and the notebook says so
    elif not GPU_PRESENT:
        print("NO GPU found in this session: Skala and gpu4pyscf are not installed. Only the G0 diagnostics will run (Runtime > Change runtime type > GPU).")
    need = [spec for mod, spec in {"pyscf": f"pyscf=={PYSCF_VERSION}", "pyscf.dispersion": f"pyscf-dispersion=={DISPERSION_VERSION}", "openpyxl": "openpyxl", "matplotlib": "matplotlib"}.items() if not _have(mod)]
    if need:
        print("installing:", need)
        _pip(*need)
except Exception as ex:
    print("install step problem (the notebook carries on):", repr(ex)[:200])

''' + banner("b1v5_common.py") + "{COMMON}" + '''

# mount Drive (the one click) or fall back to a local folder, make the folders
setup_root()
'''

CELL_HELPERS = (banner("b1v5_cpu.py") + "{CPU}" + "\n\n" + banner("b1v5_v4_stored.py") + "{V4}" + "\n\n" + banner("b1v5_geompack.py") + "{PACK}" + "\n\n" + banner("b1v5_gpu.py") + "{GPU}")
CELL_CHECK_DEFS = banner("b1v5_lit_table.py") + "{LIT}" + "\n\n" + banner("b1v5_check.py") + "{CHECK}" + "\n\n" + banner("b1v5_gpu_figs.py") + "{FIGS}"
CELL_QUEUE = banner("b1v5_1h.py") + "{QUEUE}" + "\n\n" + banner("b1v5_1h_check.py") + "{CHECK1H}" + "\n\n" + banner("b1v5_1h_gpu.py") + "{GPU1H}"

CELL_PACKDATA = '''
# ---- The geometry pack: {PACKINFO}.  SHA-256 of the pack (canonical JSON): this constant is checked in the settings cell. ----
GEOMPACK_SHA256 = "{SHA}"
GEOMPACK_B64 = "{B64}"
'''

CELL_CONFIG = '''
# ---- Settings (all optional; the defaults are the plan) ----
import os, time, json
STANDIN = False              # True: dry run on the CPU with cheap functionals in tiny bases, HCN, NH3 and CH4 only, in a separate folder. NOT results.
STANDIN = STANDIN or os.environ.get("B1V5_STANDIN") == "1"
DEADLINE_MIN = float(os.environ.get("B1V5_DEADLINE_MIN", "50"))   # no calculation is started after this minute (clock: first line of cell 1); 5 more minutes are reserved for figures, check, summary
AUTO_DISCONNECT = True       # release the Colab runtime at the very end so billing stops
EXTENSION = os.environ.get("B1V5_EXTENSION") or None             # (the environment variable is for the tests) None = the one-hour core plus the tail of this GPU class.  'GPU1' = EXT-GPU-1 (tail 1 to 3).  'GPU2' = EXT-GPU-2 (tail 4 and 5)
RETRY_FAILED = True          # True: calculations that failed earlier are tried again (second Run all: the Skala Cl- + CH3Cl TS gets the new PBE-start rescue)
FORCE_RERUN = True           # True: ignore DONE markers (finished points are still skipped), so a finished phase with failed points is looked at again
STOP_ON_G0_FAIL = True       # stop the later phases if G0 cannot reproduce the stored v4 energy (1e-4 hartree)
RUN_CORE = EXTENSION is None
RUN_TAIL = True

# ---- what hardware do we have? ----
GPU_NAME, GPU_MEM_GB, GPU_CLASS = detect_gpu()
NO_GPU = GPU_CLASS == "NONE" and not STANDIN
_plan = GPU_PLAN[GPU_CLASS]
CALC_MAX_SEC, GPU_QZ_MAX_NBF = _plan["calc_max_sec"], _plan["qz_max_nbf"]
SYSTEMS_TO_RUN = list(DESIGN_REACTIONS)
if STANDIN:
    if not ROOT.endswith("-STANDIN"):
        setup_root(ROOT + "-STANDIN", mount=False)
    cfg_ = read_json("config.json", {}) or {}
    BASIS_OVERRIDE.update(cfg_.get("basis_override") or {"small": "STO-3G", "main": "3-21G", "large": "6-31G"})
    SYSTEMS_TO_RUN = ["hcn", "nh3", "ch4"]
    SCAN_STRIDE = 4
    GPU_QZ_MAX_NBF, CALC_MAX_SEC = 400, 600
os.environ["B1V5_ROOT"] = ROOT                                             # read by the very last cell
os.environ["B1V5_AUTO_DISCONNECT"] = "1" if AUTO_DISCONNECT else "0"
q_start(DEADLINE_MIN, t0=_T0)
Q["scale"] = GPU_SPEED.get(GPU_CLASS, 1.0)                                  # the priors are T4 times (plan estimates); L4 and A100 factors are scaled guesses, the running mean takes over after the first calculations
set_design(DESIGN_REACTIONS)
quiet_geometric()
write_json("config_gpu.json", dict(design="1h", standin=STANDIN, gpu=GPU_NAME, gpu_class=GPU_CLASS, basis_override=dict(BASIS_OVERRIDE), started=now(), deadline_min=DEADLINE_MIN, extension=EXTENSION,
                                   tail_items=tail_items_for(EXTENSION, GPU_CLASS)))
print("Hardware:", hw()); print("Versions:", ver()); print("GPU:", GPU_NAME or "none", "| class:", GPU_CLASS, "| tail items by class / extension:", tail_items_for(EXTENSION, GPU_CLASS))
print("Design reactions:", DESIGN["reactions"], "| systems to run:", SYSTEMS_TO_RUN)
print(f"Deadline: minute {DEADLINE_MIN:.0f} (now minute {q_now_min():.1f}); extension: {EXTENSION}")
if STANDIN:
    print("*** STAND-IN DRY RUN: cheap CPU functionals, tiny bases, folder", ROOT, ". These numbers are NOT Skala and NOT results. ***")
if NO_GPU:
    print("*" * 100); print("NO GPU DETECTED. Choose Runtime > Change runtime type > GPU, then Run all again.\\nThis time only the G0 diagnostics will run."); print("*" * 100)
_cpu_cfg = read_json("config.json", {}) or {}
if _cpu_cfg.get("smoke") and not STANDIN:
    print("WARNING: the CPU notebook in this folder ran as a SMOKE TEST (small bases). Its CCSD(T) numbers are not the paper numbers; delete the folder and run the CPU notebook properly.")
if not _cpu_cfg and not STANDIN:
    print("NOTE: no CPU-notebook output found in this folder. Run the CPU notebook first (order: CPU, then GPU). The geometries below come from the pack anyway; the reference numbers will be missing.")
'''

MD_PACK = "**Pack and weights.** Checks the SHA-256 of the geometry pack and of every structure, writes geoms/*.json and the rows of phase01 / phase02 / phase04 (idempotent: the CPU notebook did the same). Then restores the Skala weights from the Drive cache or downloads them once with a wall-clock cap (the rest of the 15 minutes of the install cap), and caches them on Drive."
CELL_PACK = '''
PACK_INFO = pack_setup(GEOMPACK_B64, GEOMPACK_SHA256, systems=None)
print(f"geometry pack: SHA-256 {PACK_INFO.get('sha256')}\\n  expected      {GEOMPACK_SHA256}\\n  {'OK: the hash matches and every structure verified' if PACK_INFO.get('ok') else 'PROBLEM: ' + str(PACK_INFO.get('problems'))}")
if not PACK_INFO.get("ok"):
    print("*" * 100); print("The geometry pack did not verify and was NOT installed. Phases that need geometries will be skipped and say so."); print("*" * 100)
WEIGHTS_INFO = weights_setup(max(1.0, GPU_CAPS["install"] - q_now_min()))
print(f"minute {q_now_min():.1f}: install, pack and weights done")
'''

MD_G0 = "**G0: gate.** Functional probe (which of PBE, r2SCAN, B3LYP, wB97M-V run on the GPU), Skala check, SHA-256 of the weights in use, then the stored v4 HCN point 0 with Skala (D3 off, as in v4), PBE and B3LYP: must agree to 1e-5 hartree, else the later phases stop. Cap 4 minutes."
CELL_G0 = '''
run_gphase_1h("phase10", "G0 gate", g0_1h, GPU_CAPS["phase10"])
load_gstate()
print("G0 state:", {k: GSTATE.get(k) for k in ("routes", "anchor", "gate_ok", "skala_ok", "skala_route", "grad_ok")})
'''
MD_G1A = "**G1a: five functionals at the stationary points.** All species of the six reactions (20), own B3LYP geometries from the pack, def2-TZVP (def2-TZVPD for SN2). Order: Skala (with its noD3 row), wB97M-V, B3LYP-D3(BJ), r2SCAN, PBE: if time runs out PBE and r2SCAN are cut first. Cap 18 minutes."
CELL_G1A = '''
if RUN_CORE:
    run_gphase_1h("phase11", "G1a stationary points", lambda: (g1a_stationary(), g1a_align()), GPU_CAPS["phase11"])
'''
MD_G1B = "**G1b: paths.** Skala, wB97M-V and B3LYP-D3(BJ) along the stored relaxed paths of HCN (11), NH3 (9), H2CO (11) and the CH4 stretch (9): 40 points per functional, one functional at a time. Cap 12 minutes."
CELL_G1B = '''
if RUN_CORE:
    run_gphase_1h("phase11p", "G1b path points", g1b_paths, GPU_CAPS["phase11p"])
'''
MD_G3 = "**G3: scout check.** Skala, wB97M-V and B3LYP-D3(BJ) at the xTB minimum and peak of the six reactions (geometries from the pack), and xTB's own barriers. Cap 6 minutes."
CELL_G3 = '''
if RUN_CORE:
    run_gphase_1h("phase13", "G3 scout check", g3_scout_1h, GPU_CAPS["phase13"])
'''
MD_G2 = "**G2: basis check.** def2-QZVPP at the 9 stationary species of HCN, NH3, H2CO with Skala and wB97M-V, system by system so that a barrier is complete. Cap 8 minutes."
CELL_G2 = '''
if RUN_CORE:
    run_gphase_1h("phase12", "G2 basis check", g2_basis_1h, GPU_CAPS["phase12"])
'''
MD_TAIL = "**Tail list** (same admission test; ends at minute 50). T4: expected to be empty. L4: reactions 7 and 8 (F- + CH3Cl, OH- + CH3F), perhaps PBE and r2SCAN along the paths. A100: items 1 to 4, perhaps the Skala-relaxed HCN. `EXTENSION = 'GPU1'` runs items 1 to 3 and `'GPU2'` items 4 and 5 in a session of their own."
CELL_TAIL = '''
_tail = tail_items_for(EXTENSION, GPU_CLASS) if (RUN_TAIL or EXTENSION) else []
print("tail items to try:", _tail)
run_tail(_tail)
'''
MD_G5 = "**G5: figures and tables.** Reactions 7 and 8 enter the figures and the check only if their stationary points were computed (tail item 1). Regenerates the figures, `numbers_v5.csv`, the tables and `findings_outcome.csv` (criteria frozen on 2026-11-15, five of six reactions) from all CSV files in the Drive folder, on whatever data exists; a figure or table whose data is missing is skipped with a message. Always runs, cap 5 minutes."
CELL_G5 = '''
set_design(design_with_extension_data())
print("design for the figures and the check:", DESIGN["reactions"])
run_gphase_1h("phase15", "G5 figures and tables", lambda: g5_make_all(STANDIN), GPU_CAPS["phase15"], final=True, force=True)
'''
MD_CHECK = "**Automated check.** The same PASS / WARN / FAIL check as C7 of the CPU notebook, now with the GPU results; items that cannot exist in the one-hour design are NOT RUN. Always runs."
CELL_CHECK = '''
def _run_check():
    global check_result
    check_result = run_check()
if not (NO_GPU and not STANDIN):
    with qphase("phase15x_check", "automated check", GPU_CAPS["check"], final=True, force=True) as ph:
        _run_check()
'''
MD_SUM = "**Summary and release.** Writes `RUN_SUMMARY.md` (what ran, what failed, what was not run for time, minutes against the caps, GPU type) to the Drive folder."
CELL_SUM = '''
try:
    q_save_log()
    print(write_run_summary(extra_lines=g_summary_lines_1h() + ["", f"Geometry pack: SHA-256 {PACK_INFO.get('sha256')} (expected {GEOMPACK_SHA256}), verified: {PACK_INFO.get('ok')}, problems: {PACK_INFO.get('problems')}",
                                                                 f"Extension: {EXTENSION}; tail items tried: {_tail}; deadline {DEADLINE_MIN} min; GPU class {GPU_CLASS}"]))
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


def build(out_name="B1v5_GPU_1h.ipynb", pack_path=None):
    import b1v5_geompack as G
    pack_path = pack_path or os.path.join(HERE, "geompack", "b1v5_geompack_1h.json")
    raw = open(pack_path, "rb").read()
    pack = json.loads(raw.decode())
    assert G.pack_bytes(pack) == raw, "the pack file is not in canonical form (rebuild it with build_geompack.py)"
    sha = hashlib.sha256(raw).hexdigest()
    b64 = G.pack_encode(pack)
    info = f"{len(pack['structures'])} structures, {len(pack['paths'])} paths, {len(pack['xtb'])} xTB sets, {len(raw) // 1024} KB JSON, {len(b64) // 1024} KB embedded; systems {', '.join(pack['systems'])}; built {pack['created'][:10]} with {pack['builder']['versions'].split(';')[0]}"
    sub = lambda s: s.replace("{COMMON}", read("b1v5_common.py")).replace("{CPU}", read("b1v5_cpu.py")).replace("{V4}", read("b1v5_v4_stored.py")).replace("{PACK}", read("b1v5_geompack.py")) \
        .replace("{GPU}", read("b1v5_gpu.py")).replace("{LIT}", read("b1v5_lit_table.py")).replace("{CHECK}", read("b1v5_check.py")).replace("{FIGS}", read("b1v5_gpu_figs.py")) \
        .replace("{QUEUE}", read("b1v5_1h.py")).replace("{CHECK1H}", read("b1v5_1h_check.py")).replace("{GPU1H}", read("b1v5_1h_gpu.py"))
    cells = [md(INTRO),
             md("**Cell 1: the clock, install (cap 15 minutes), shared helpers, Drive.** The clock of the 50-minute deadline starts in the first line. Installs the pinned pyscf 2.12.1, pyscf-dispersion, Skala 2026.9 and gpu4pyscf 1.8.1 (if a GPU is present), defines the shared helper functions and mounts Google Drive (one permission click)."), code(sub(CELL_SETUP)),
             md("**Cell 2: calculation helpers** (CPU helpers, stored v4 numbers, geometry-pack tools, GPU helpers). Only defines functions and constants."), code(sub(CELL_HELPERS)),
             md("**Cell 3: check module, literature table and figure code** (pasted unchanged). Only defines functions."), code(sub(CELL_CHECK_DEFS)),
             md("**Cell 4: deadline queue, six-reaction design and GPU phases.** Only defines functions and constants."), code(sub(CELL_QUEUE)),
             md("**Cell 5: the geometry pack itself** (data, with its SHA-256)."), code(CELL_PACKDATA.replace("{PACKINFO}", info).replace("{SHA}", sha).replace("{B64}", b64)),
             md("**Cell 6: settings.** Deadline, extension sessions, dry-run flag, GPU detection. Starts the queue clock."), code(CELL_CONFIG),
             md(MD_PACK), code(CELL_PACK),
             md(MD_G0), code(CELL_G0), md(MD_G1A), code(CELL_G1A), md(MD_G1B), code(CELL_G1B), md(MD_G3), code(CELL_G3), md(MD_G2), code(CELL_G2),
             md(MD_TAIL), code(CELL_TAIL), md(MD_G5), code(CELL_G5), md(MD_CHECK), code(CELL_CHECK), md(MD_SUM), code(CELL_SUM), md(MD_END), code(CELL_END)]
    nb = nbf.v4.new_notebook(cells=cells)
    nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    nb.metadata["colab"] = {"provenance": [], "name": out_name, "gpuType": "T4", "machine_shape": "hm"}
    nb.metadata["accelerator"] = "GPU"
    nbf.validate(nb)
    nbf.write(nb, os.path.join(HERE, out_name) if not os.path.isabs(out_name) else out_name)
    return out_name, sha


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pack", default=None)
    ap.add_argument("--out", default="B1v5_GPU_1h.ipynb")
    a = ap.parse_args()
    p, sha = build(a.out, a.pack)
    full = p if os.path.isabs(p) else os.path.join(HERE, p)
    print("wrote", full, os.path.getsize(full) // 1024, "KB; pack SHA-256", sha)
