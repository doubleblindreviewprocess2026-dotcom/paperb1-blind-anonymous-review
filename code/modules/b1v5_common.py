# ==============================================================================
# b1v5_common.py  --  shared helpers for the Paper B1 v5 notebooks (CPU and GPU).
#
# This file is pasted unchanged into the first code cell of BOTH notebooks, and
# it can also be imported as a normal module.  It only needs numpy.
#
# What is in here, in plain words:
#   1. where the files live (Google Drive folder or a local folder)
#   2. one results table per phase (CSV), written after EVERY calculation,
#      so a disconnect never loses work and a rerun skips what is done
#   3. DONE_phaseXX marker files
#   4. geometry files (geoms/*.json and geoms/*.xyz)
#   5. the system database: the eight reactions plus the CH4 stretch, with
#      start geometries, reference tiers and the species each reaction needs
#   6. unit constants, hardware and package-version recorders
# ==============================================================================
import os, sys, json, csv, io, math, time, platform, datetime, importlib.metadata
import numpy as np

# ------------------------------------------------------------------ constants
HARTREE_KCAL = 627.5094740631          # 1 hartree in kcal/mol
BOHR = 0.529177210903                  # 1 bohr in Angstrom
HARTREE_CM1 = 219474.6313705           # 1 hartree in wavenumbers
EV_HARTREE = 1.0 / 27.211386245988
KCAL = HARTREE_KCAL                    # short alias

# ------------------------------------------------------------------ names used in the CSV files
# Geometry levels (column "geometry_level")
LV_DB = "db_qcisd_mg3"        # published NHTBH38/08 database geometry
LV_OWN = "own_b3lyp"          # our DF-B3LYP-D3(BJ)/def2-TZVP(D) optimisation
LV_SCANPEAK = "own_scanpeak"  # transition state taken from the top of our relaxed scan (fallback)
LV_SCAN = "own_scan"          # a point of our relaxed scan
LV_XTB = "xtb_gfn2"           # GFN2-xTB geometry (scout)
LV_SKALA = "skala_relaxed"    # Skala-relaxed geometry (GPU notebook)
LV_GRID = "rigid_grid"        # v4 rigid HCN grid
LV_NB4 = "v4_nb4"             # v4 NB4 path geometry
# Order in which a stored geometry is picked when the same species exists at several levels (review 2026-11-15, row 4):
# our own B3LYP geometry first (the HEADLINE protocol), the published database geometry last (I11 cross-check only).
LEVEL_PREFERENCE = [LV_OWN, LV_SCANPEAK, LV_DB]
# A barrier or reaction energy needs ALL its species at ONE geometry family; a family is a set of levels that count as "the same protocol".
# own_scanpeak is a transition state taken from the top of our own B3LYP scan: same B3LYP protocol, flagged in the label.
FAMILY_LEVELS = {"own": (LV_OWN, LV_SCANPEAK), "db": (LV_DB,)}
FAMILY_OF = {LV_OWN: "own", LV_SCANPEAK: "own", LV_DB: "db"}
HEADLINE_FAMILY = "own"       # 'own' (default), 'db', or 'any' (first complete family per quantity: own, then db); set from config.json by the check

# Method names (column "method").  Both notebooks must use exactly these strings.
# Skala: "Skala-1.1" is the HEADLINE = the published functional WITH its D3 term (skala package default with_dftd3=True; the D3 energy is
# stored separately in "extra"); "Skala-1.1-noD3" = the same SCF minus that D3 energy (written as a second row).  v4 used no D3 ("v4-skala").
M_SKALA, M_PBE, M_R2SCAN, M_B3LYP, M_WB97MV = "Skala-1.1", "PBE", "r2SCAN", "B3LYP-D3BJ", "wB97M-V"
M_SKALA_NOD3 = "Skala-1.1-noD3"
M_WB97XV = "wB97X-V"                       # fallback anchor (plan section 5)
M_CCSDT_FC, M_CCSDT_AE, M_XTB, M_B2PLYP, M_DFB3LYP = "CCSD(T)-fc", "CCSD(T)-ae", "xTB-GFN2", "B2PLYP", "DF-B3LYP-D3BJ"
M_MP2DF = "DF-MP2-fc"                      # density-fitted MP2, frozen core: only used for the basis correction of the CCSD(T) reference
M_CCSDT_BC = "CCSD(T)+dMP2"                # DERIVED (never stored): CCSD(T)-fc/TZVP(D) + [MP2/QZVPP(D) - MP2/TZVP(D)]
LADDER = [M_SKALA, M_PBE, M_R2SCAN, M_B3LYP, M_WB97MV]          # the five functionals of the paper (Skala = with D3)
LADDER_ALL = [M_SKALA, M_SKALA_NOD3, M_PBE, M_R2SCAN, M_B3LYP, M_WB97MV]
ANCHORS = (M_WB97MV, M_WB97XV)             # anchor functional and its fallback
CORE_METHODS = (M_DFB3LYP, M_CCSDT_FC, M_WB97MV, M_WB97XV)      # methods whose failure can affect a system (checks I1, I6)

COLUMNS = ["system", "point_id", "coordinate", "geometry_level", "method", "basis",
           "energy_hartree", "converged", "t1", "imaginary_modes", "seconds", "hardware", "versions", "extra"]
# "extra" is a small JSON text for things that do not deserve their own column
# (zero-point energy, RHF stability, error messages, frozen-core flag ...).

PHASE_CSV = {  # file names for the phases; GPU names follow the plan
    "phase00": "phase00_timing.json", "phase01": "phase01_stationary.csv", "phase02": "phase02_scans.csv",
    "phase03": "phase03_ccsdt.csv", "phase04": "phase04_scout_ref.csv", "phase05": "phase05_hcn_grid.csv",
    "phase06": "phase06_ccsdt_basis.csv", "phase07": "phase07_optional.csv", "phase08": "phase08_mp2_basis.csv",
    "phase10": "phase10_timing.json", "phase11": "phase11_dft.csv", "phase12": "phase12_basis.csv",
    "phase13": "phase13_scout.csv", "phase14": "phase14_skala_relaxed.csv"}

# ------------------------------------------------------------------ root folder
ROOT = None


def in_colab():
    return "google.colab" in sys.modules


def setup_root(root=None, mount=True):
    """Pick the working folder and create it.  Order: argument, env B1V5_ROOT, Google Drive, local folder."""
    global ROOT
    root = root or os.environ.get("B1V5_ROOT")
    if not root and mount and in_colab():
        from google.colab import drive
        for attempt in range(3):                      # the mount sometimes fails on the first try (popup not answered in time, credentials not yet propagated)
            try:
                drive.mount("/content/drive", force_remount=(attempt > 0))
                if os.path.isdir("/content/drive/MyDrive"):
                    root = "/content/drive/MyDrive/skala-b1-v5"
                    break
            except Exception as ex:
                print(f"Google Drive mount, attempt {attempt + 1} of 3 failed ({repr(ex)[:80]}); " + ("trying again in 20 s: please click Allow in the Drive popup." if attempt < 2 else ""))
            if attempt < 2:
                time.sleep(20)
        if root is None:   # no Drive: carry on locally and say so
            print("*" * 90); print("COULD NOT MOUNT GOOGLE DRIVE -> using a local folder (/content/skala-b1-v5-local). Nothing is saved to Drive.")
            print("The last cell will then NOT release the runtime and will offer a zip to download. Better: stop, Runtime > Restart, Run all again and click Allow."); print("*" * 90)
    root = root or os.path.join(os.getcwd(), "skala-b1-v5-local")
    os.makedirs(os.path.join(root, "geoms"), exist_ok=True)
    ROOT = os.path.abspath(root)
    print("Working folder:", ROOT)
    return ROOT


def path(*parts):
    if ROOT is None:
        setup_root(mount=False)
    return os.path.join(ROOT, *parts)


# ------------------------------------------------------------------ hardware and versions
def hardware_info():
    """One short text line: CPU model, cores, memory, GPU (if any)."""
    cpu = platform.processor() or platform.machine()
    try:
        for line in open("/proc/cpuinfo"):
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip(); break
    except Exception:
        pass
    mem = ""
    try:
        for line in open("/proc/meminfo"):
            if line.startswith("MemTotal"):
                mem = f", {int(line.split()[1]) / 1048576:.0f} GB RAM"; break
    except Exception:
        pass
    gpu = ""
    try:
        import subprocess
        out = subprocess.run(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"], capture_output=True, text=True, timeout=10).stdout.strip()
        gpu = f", GPU {out.splitlines()[0]}" if out else ""
    except Exception:
        pass
    return f"{cpu}, {os.cpu_count()} cores{mem}{gpu}"


def package_versions(names=("pyscf", "geometric", "tblite", "ase", "numpy", "scipy", "pandas", "skala", "skala-cuda12x", "gpu4pyscf",
                            "gpu4pyscf-cuda12x", "cupy-cuda12x", "torch", "huggingface_hub", "pyscf-dispersion")):
    out = []
    for n in names:
        try:
            out.append(f"{n}={importlib.metadata.version(n)}")
        except Exception:
            pass
    return ";".join(out) + f";python={platform.python_version()}"


_HW, _VER = [None], [None]


def hw():
    if _HW[0] is None:
        _HW[0] = hardware_info()
    return _HW[0]


def ver():
    if _VER[0] is None:
        _VER[0] = package_versions()
    return _VER[0]


class Timer:
    """with Timer() as t: ...  then t.seconds."""
    def __enter__(self):
        self.t0 = time.time(); return self

    def __exit__(self, *a):
        self.seconds = time.time() - self.t0


def now():
    return datetime.datetime.now().isoformat(timespec="seconds")


# ------------------------------------------------------------------ results tables (CSV, long format)
def make_row(system, point_id, level, method, basis, energy, converged=True, coordinate="", t1="", imag="",
             seconds=0.0, extra=None):
    return dict(system=system, point_id=point_id, coordinate=coordinate, geometry_level=level, method=method, basis=basis,
                energy_hartree=("" if energy is None else repr(float(energy))), converged=int(bool(converged)), t1=t1,
                imaginary_modes=imag, seconds=f"{seconds:.2f}", hardware=hw(), versions=ver(),
                extra=json.dumps(extra or {}, sort_keys=True))


def row_key(r):
    return (str(r["system"]), str(r["point_id"]), str(r["geometry_level"]), str(r["method"]), str(r["basis"]))


def read_rows(csv_name):
    """All rows of a CSV as a list of dicts ([] when the file does not exist)."""
    p = path(csv_name)
    if not os.path.exists(p):
        return []
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh))


def append_row(csv_name, row):
    """Add one row safely: write a temporary file, then rename it over the old one."""
    p = path(csv_name)
    old = open(p).read() if os.path.exists(p) else ",".join(COLUMNS) + "\n"
    if old and not old.endswith("\n"):
        old += "\n"
    buf = io.StringIO()
    csv.DictWriter(buf, fieldnames=COLUMNS, extrasaction="ignore", lineterminator="\n").writerow(row)
    tmp = p + ".tmp"
    with open(tmp, "w") as fh:
        fh.write(old + buf.getvalue()); fh.flush(); os.fsync(fh.fileno())
    os.replace(tmp, p)


def have(csv_name, system, point_id, level, method, basis, retry_failed=False):
    """True when this exact calculation is already in the CSV (so we skip it)."""
    key = (system, str(point_id), level, method, basis)
    for r in read_rows(csv_name):
        if row_key(r) == key and not (retry_failed and r["converged"] in ("0", "False")):
            return True
    return False


def run_point(csv_name, system, point_id, level, method, basis, fn, coordinate="", retry_failed=False, quiet=False):
    """Do one calculation unless it is already saved.

    fn() must return a dict with 'energy' and optionally 'converged', 't1', 'imag', 'extra'.
    A crash inside fn is stored as a failed row (converged=0) and the run goes on.
    Returns the new row, or None when skipped.
    """
    if have(csv_name, system, point_id, level, method, basis, retry_failed):
        if not quiet:
            print(f"  skip (already done): {system} {point_id} {method} {basis}")
        return None
    check_time()                      # stops the phase cleanly if a time limit has passed
    t0 = time.time()
    try:
        out = fn()
        row = make_row(system, point_id, level, method, basis, out.get("energy"), out.get("converged", True), coordinate,
                       out.get("t1", ""), out.get("imag", ""), time.time() - t0, out.get("extra"))
    except PhaseTimeout:
        raise
    except Exception as ex:
        row = make_row(system, point_id, level, method, basis, float("nan"), False, coordinate, "", "", time.time() - t0,
                       {"error": repr(ex)[:300]})
    append_row(csv_name, row)
    if not quiet:
        e = row["energy_hartree"]
        print(f"  {system} {point_id} {method} {basis}: E={float(e) if e not in ('', 'nan') else float('nan'):.8f} "
              f"conv={row['converged']} {row['seconds']} s")
    return row


def extra_of(row):
    try:
        return json.loads(row.get("extra") or "{}")
    except Exception:
        return {}


# ------------------------------------------------------------------ DONE markers and a small phase ledger
def is_done(phase):
    return os.path.exists(path(f"DONE_{phase}"))


def mark_done(phase, note=""):
    with open(path(f"DONE_{phase}"), "w") as fh:
        fh.write(f"{now()} {note}\n")


def log_phase(phase, seconds, note=""):
    """Keep a running list of what ran and how long it took (phase_log.json)."""
    p = path("phase_log.json")
    log = json.load(open(p)) if os.path.exists(p) else {}
    rec = log.get(phase, {"seconds": 0.0, "runs": 0})
    rec["seconds"] += seconds; rec["runs"] += 1; rec["last"] = now(); rec["note"] = note
    log[phase] = rec
    json.dump(log, open(p + ".tmp", "w"), indent=1); os.replace(p + ".tmp", p)


def write_json(name, obj):
    p = path(name)
    json.dump(obj, open(p + ".tmp", "w"), indent=1, default=float); os.replace(p + ".tmp", p)


def read_json(name, default=None):
    p = path(name)
    return json.load(open(p)) if os.path.exists(p) else default


# ------------------------------------------------------------------ guards: errors never stop the notebook, clocks stop phases cleanly
import traceback, contextlib

MAX_TOTAL_HOURS = 3.2        # whole-notebook budget (the notebook sets its own value)
DEFAULT_PHASE_MIN = 85.0     # per-phase wall-time limit in minutes (plan: no phase over 85 min)
_CLOCK = {"t0": time.time(), "phase": None, "deadline": None, "total_deadline": None}


class PhaseTimeout(Exception):
    """Raised at the next calculation after a time limit has passed; the phase then stops cleanly."""


def start_clock(max_total_hours=None):
    """Call once at the top of the notebook."""
    global MAX_TOTAL_HOURS
    if max_total_hours is not None:
        MAX_TOTAL_HOURS = max_total_hours
    _CLOCK["t0"] = time.time()
    _CLOCK["total_deadline"] = _CLOCK["t0"] + 3600.0 * MAX_TOTAL_HOURS


def elapsed_hours():
    return (time.time() - _CLOCK["t0"]) / 3600.0


def check_time():
    """Called before every single calculation (run_point does it for you)."""
    now_ = time.time()
    if _CLOCK["deadline"] and now_ > _CLOCK["deadline"]:
        raise PhaseTimeout(f"phase time limit reached in {_CLOCK['phase']}")
    if _CLOCK["total_deadline"] and now_ > _CLOCK["total_deadline"]:
        raise PhaseTimeout(f"total time budget of {MAX_TOTAL_HOURS} h reached")


ERRORS_CSV = "phase_errors.csv"


def log_error(phase, system, exc, note=""):
    """Append one line (phase, system, short traceback) to phase_errors.csv. Never raises."""
    try:
        tb = "".join(traceback.format_exception_only(type(exc), exc)).strip()
        where = traceback.extract_tb(exc.__traceback__)[-1:] if exc.__traceback__ else []
        where = "; ".join(f"{os.path.basename(f.filename)}:{f.lineno} {f.name}" for f in where)
        p = path(ERRORS_CSV)
        new = not os.path.exists(p)
        with open(p, "a", newline="") as fh:
            w = csv.writer(fh)
            if new:
                w.writerow(["time", "phase", "system", "error", "where", "note"])
            w.writerow([now(), phase, system, tb[:300], where[:200], note])
        print(f"  !! {phase} {system}: {tb[:200]}  (logged in {ERRORS_CSV})")
    except Exception:
        print("  !! could not log error", repr(exc)[:200])


@contextlib.contextmanager
def item_guard(phase, system=""):
    """Wrap ONE item (a system, a species, a point) of a loop: a crash is logged and the loop carries on.
    A time-limit stop (PhaseTimeout) is passed up so the whole phase ends cleanly."""
    try:
        yield
    except PhaseTimeout:
        raise
    except Exception as ex:
        log_error(phase, system, ex)


def _set_status(phase, **kw):
    st = read_json("phase_status.json", {})
    st.setdefault(phase, {}).update(kw)
    write_json("phase_status.json", st)


@contextlib.contextmanager
def phase(name, title="", max_min=None, force=False):
    """Wrap a whole phase.  Whatever happens inside, the notebook goes on to the next cell.
        with phase("phase01", "C1 geometries", max_min=40) as ph:
            if ph.skip: ...           # (the body is simply skipped when the phase is already DONE)
    The phase is marked DONE only if it ran to the end without error or time-out."""
    class _Ph:
        skip = False
    ph = _Ph()
    print(f"\n===== {name} {title} =====  (notebook time so far {elapsed_hours():.2f} h)")
    if is_done(name) and not force:
        print(f"  {name} already DONE (marker found) -> skipped. Delete DONE_{name} to run it again.")
        _set_status(name, status="already done", title=title)
        ph.skip = True
        yield ph
        return
    t0 = time.time()
    _CLOCK["phase"] = name
    _CLOCK["deadline"] = t0 + 60.0 * (max_min or DEFAULT_PHASE_MIN)
    status, note = "done", ""
    try:
        check_time()
    except PhaseTimeout as ex:          # the total budget is already used up: do not start this phase
        ph.skip = True; status, note = "stopped (time limit)", str(ex)
        print("  ** not started:", ex)
    try:
        yield ph
    except PhaseTimeout as ex:
        status, note = "stopped (time limit)", str(ex)
        print("  ** phase stopped cleanly:", ex, "- rerun the notebook to continue, finished points are kept.")
    except Exception as ex:
        status, note = "FAILED", repr(ex)[:200]
        log_error(name, "(whole phase)", ex)
    finally:
        _CLOCK["deadline"] = None
        secs = time.time() - t0
        n_err = 0
        try:
            n_err = sum(1 for r in csv.DictReader(open(path(ERRORS_CSV))) if r["phase"] == name) if os.path.exists(path(ERRORS_CSV)) else 0
        except Exception:
            pass
        if status == "done" and n_err:
            status = "done with errors"
        if status == "done":
            mark_done(name, title)
        log_phase(name, secs, status)
        _set_status(name, status=status, title=title, seconds=round(secs, 1), errors=n_err, note=note, finished=now())
        print(f"  {name}: {status}, {secs / 60:.1f} min")


def skip_phase(name, title, reason):
    """Record 'skipped: missing input ...' and carry on."""
    print(f"  {name} skipped: {reason}")
    _set_status(name, status="skipped", title=title, note=reason, finished=now())
    log_error(name, "", RuntimeError("skipped: " + reason), note="skipped")


def missing_inputs(*names):
    return [n for n in names if not os.path.exists(path(n))]


def write_run_summary(extra_lines=()):
    """RUN_SUMMARY.md: what ran, what failed, what was skipped, total hours."""
    st = read_json("phase_status.json", {})
    L = ["# Run summary", "", f"Written: {now()}", f"Hardware: {hw()}", f"Versions: {ver()}",
         f"Notebook time: {elapsed_hours():.2f} h (budget {MAX_TOTAL_HOURS} h)", "", "| phase | title | status | minutes | errors | note |", "|---|---|---|---|---|---|"]
    for k in sorted(st):
        v = st[k]
        L.append(f"| {k} | {v.get('title', '')} | {v.get('status', '')} | {v.get('seconds', 0) / 60:.1f} | {v.get('errors', '')} | {str(v.get('note', ''))[:80]} |")
    ran = [k for k in st if st[k].get("status", "").startswith(("done", "already"))]
    bad = [k for k in st if st[k].get("status", "") in ("FAILED",) or "error" in st[k].get("status", "") or "stopped" in st[k].get("status", "")]
    skp = [k for k in st if st[k].get("status") == "skipped"]
    L += ["", f"Ran: {', '.join(sorted(ran)) or 'nothing'}", f"Failed or stopped early: {', '.join(sorted(bad)) or 'none'}",
          f"Skipped: {', '.join(sorted(skp)) or 'none'}", ""]
    L += list(extra_lines)
    if os.path.exists(path(ERRORS_CSV)):
        L += ["", "## Errors (phase_errors.csv)", ""]
        for r in list(csv.DictReader(open(path(ERRORS_CSV))))[:60]:
            L.append(f"- {r['phase']} / {r['system']}: {r['error'][:140]}")
    open(path("RUN_SUMMARY.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L[:30]))
    return "\n".join(L)


# ------------------------------------------------------------------ geometries
# An "atoms" list is [(symbol, (x, y, z)), ...] in Angstrom (same format as fpes.py).
def write_xyz(p, atoms, comment=""):
    lines = [str(len(atoms)), comment] + [f"{s:2s} {x: .8f} {y: .8f} {z: .8f}" for s, (x, y, z) in atoms]
    open(p, "w").write("\n".join(lines) + "\n")


def read_xyz(p):
    L = open(p).read().strip().splitlines()
    n = int(L[0].split()[0])
    return [(l.split()[0], tuple(float(v) for v in l.split()[1:4])) for l in L[2:2 + n]]


def save_geom(name, atoms, charge=0, **meta):
    """Write geoms/<name>.json and geoms/<name>.xyz (atomic)."""
    d = dict(name=name, charge=int(charge), atoms=[[s, list(map(float, c))] for s, c in atoms], saved=now(), **meta)
    pj = path("geoms", name + ".json")
    json.dump(d, open(pj + ".tmp", "w"), indent=1, default=float); os.replace(pj + ".tmp", pj)
    write_xyz(path("geoms", name + ".xyz"), atoms, f"{name} charge={charge} {json.dumps(meta, default=float)[:120]}")


def load_geom(name):
    """Return dict(atoms=[(sym,(x,y,z))...], charge=..., ...) or None."""
    pj = path("geoms", name + ".json")
    if not os.path.exists(pj):
        return None
    d = json.load(open(pj))
    d["atoms"] = [(s, tuple(c)) for s, c in d["atoms"]]
    return d


def geom_name(system, species, level):
    return f"{system}_{species}_{level}".replace("-", "m").replace("+", "p")


def best_geom(system, species):
    """The preferred stored geometry of a stationary species (our own B3LYP first, then the scan-peak TS, then the database).  Returns (level, geom) or (None, None)."""
    for lv in LEVEL_PREFERENCE:
        g = load_geom(geom_name(system, species, lv))
        if g:
            return lv, g
    return None, None


def geometry_of(system, species, level):
    return load_geom(geom_name(system, species, level))


# ------------------------------------------------------------------ geometry builders (from fpes.py; same formulas as v4)
M_C, M_N = 12.0, 14.003074
R_CN_DEFAULT = 1.16


def hcn_geometry(theta_deg, r, R_CN=R_CN_DEFAULT):
    """HCN/HNC rigid grid of v4: C and N on the z axis, H at distance r from the C-N centre of mass.
    theta = 0: H on the C side (HCN); theta = 180: H on the N side (HNC).  (Copy of fpes.hcn_geometry.)"""
    zc = -R_CN * M_N / (M_C + M_N)
    zn = R_CN * M_C / (M_C + M_N)
    t = math.radians(theta_deg)
    return [("C", (0.0, 0.0, zc)), ("N", (0.0, 0.0, zn)), ("H", (r * math.sin(t), 0.0, -r * math.cos(t)))]


def hcn_grid_axes(n_theta=25, n_r=17):
    return np.linspace(0.0, 180.0, n_theta), np.linspace(1.05, 2.30, n_r)


def nh3_geometry(umbrella_deg, r_NH=1.012):
    """NH3 with the three H below the plane by `umbrella_deg` (0 = planar TS, about 22 = minimum).  (Copy of fpes.nh3_geometry.)"""
    t = math.radians(umbrella_deg)
    rho, z = r_NH * math.cos(t), -r_NH * math.sin(t)
    atoms = [("N", (0.0, 0.0, 0.0))]
    for k in range(3):
        p = 2 * math.pi * k / 3
        atoms.append(("H", (rho * math.cos(p), rho * math.sin(p), z)))
    return atoms


def ch4_geometry(stretch, r_CH=1.087):
    """CH4 with one C-H bond stretched; the other three at tetrahedral directions.  (Copy of fpes.ch4_geometry.)"""
    v = 1 / math.sqrt(3)
    dirs = [(v, v, v), (-v, -v, v), (-v, v, -v), (v, -v, -v)]
    atoms = [("C", (0.0, 0.0, 0.0))]
    for k, d in enumerate(dirs):
        L = stretch if k == 0 else r_CH
        atoms.append(("H", (d[0] * L, d[1] * L, d[2] * L)))
    return atoms


CH4_STRETCH = np.linspace(0.90, 1.70, 9)   # the 9 stored v4 points


def _unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


def _ch3(c, axis, r=1.09, ang=109.5, phase=0.0):
    """Three H atoms of a CH3 group on carbon at c.  `axis` points from C to the substituent;
    the H atoms sit at angle `ang` from that axis (109.5 = normal tetrahedral CH3X)."""
    c, a = np.asarray(c, float), _unit(axis)
    t = np.array([1.0, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1.0, 0])
    u = _unit(t - a * np.dot(t, a)); w = np.cross(a, u)
    out = []
    for k in range(3):
        p = phase + 2 * math.pi * k / 3
        d = math.cos(math.radians(ang)) * a + math.sin(math.radians(ang)) * (math.cos(p) * u + math.sin(p) * w)
        out.append(("H", tuple(c + r * d)))
    return out


def _pt(s, x, y, z):
    return (s, (float(x), float(y), float(z)))


# ------------------------------------------------------------------ the system database
# tier A: CCSD(T) on whole path; B: CCSD(T) at stationary points only; C: literature + ladder only.
# zpe: how the literature value treats zero-point energy ("excluded", "included", "unknown").
# reactants / products: the separated species that define the reaction energy and "from reactants" barriers.
# rc / pc: reactant / product complex species (SN2 only).  For symmetric SN2 reactions pc is the same species as rc.
SYSTEMS = {
    "hcn": dict(n=1, name="HCN to HNC", cls="H shift", tier="A", sn2=False, db=True, zpe="excluded",
                reactants=["HCN"], products=["HNC"], rc=None, pc=None, ts="TS"),
    "nh3": dict(n=2, name="NH3 inversion", cls="inversion", tier="A", sn2=False, db=False, zpe="excluded",
                reactants=["NH3"], products=["NH3"], rc=None, pc=None, ts="TS"),
    "h2co": dict(n=3, name="H2CO to H2 + CO", cls="elimination", tier="A", sn2=False, db=False, zpe="included",
                 reactants=["H2CO"], products=["H2", "CO"], rc=None, pc=None, ts="TS"),
    "ch3nc": dict(n=4, name="CH3NC to CH3CN", cls="group shift", tier="B", sn2=False, db=False, zpe="unknown",
                  reactants=["CH3NC"], products=["CH3CN"], rc=None, pc=None, ts="TS"),
    "sn2_ff": dict(n=5, name="F- + CH3F", cls="SN2", tier="B", sn2=True, db=True, zpe="excluded",
                   reactants=["F-", "CH3F"], products=["F-", "CH3F"], rc="RC", pc="RC", ts="TS"),
    "sn2_clcl": dict(n=6, name="Cl- + CH3Cl", cls="SN2", tier="C", sn2=True, db=True, zpe="excluded",
                     reactants=["Cl-", "CH3Cl"], products=["Cl-", "CH3Cl"], rc="RC", pc="RC", ts="TS"),
    "sn2_fcl": dict(n=7, name="F- + CH3Cl", cls="SN2", tier="C", sn2=True, db=True, zpe="excluded",
                    reactants=["F-", "CH3Cl"], products=["CH3F", "Cl-"], rc="RC", pc="PC", ts="TS"),
    "sn2_ohf": dict(n=8, name="OH- + CH3F", cls="SN2", tier="C", sn2=True, db=True, zpe="excluded",
                    reactants=["OH-", "CH3F"], products=["CH3OH", "F-"], rc="RC", pc="PC", ts="TS"),
    "ch4": dict(n=9, name="CH4 C-H stretch", cls="bond stretch", tier="A", sn2=False, db=False, zpe="n/a",
                reactants=["CH4"], products=[], rc=None, pc=None, ts=None),
}
REACTIONS = [k for k in SYSTEMS if k != "ch4"]          # the eight reactions
SYSTEM_ORDER = list(SYSTEMS)
TIER_B = [k for k, v in SYSTEMS.items() if v["tier"] == "B"]

# Which species (stationary points and fragments) each system needs, and what kind they are:
# 'min' = minimum, 'ts' = transition state.
SPECIES_KIND = {
    "hcn": {"HCN": "min", "HNC": "min", "TS": "ts"},
    "nh3": {"NH3": "min", "TS": "ts"},
    "h2co": {"H2CO": "min", "TS": "ts", "H2": "min", "CO": "min"},
    "ch3nc": {"CH3NC": "min", "CH3CN": "min", "TS": "ts"},
    "sn2_ff": {"F-": "min", "CH3F": "min", "RC": "min", "TS": "ts"},
    "sn2_clcl": {"Cl-": "min", "CH3Cl": "min", "RC": "min", "TS": "ts"},
    "sn2_fcl": {"F-": "min", "CH3Cl": "min", "RC": "min", "TS": "ts", "PC": "min", "CH3F": "min", "Cl-": "min"},
    "sn2_ohf": {"OH-": "min", "CH3F": "min", "RC": "min", "TS": "ts", "PC": "min", "CH3OH": "min", "F-": "min"},
    "ch4": {"CH4": "min"},
}
CHARGE = {"F-": -1, "Cl-": -1, "OH-": -1}


def species_charge(system, species):
    """Total charge: -1 for the anion, its complexes and the TS of every SN2 system, else 0."""
    if SYSTEMS[system]["sn2"]:
        return 0 if species in ("CH3F", "CH3Cl", "CH3OH") else -1
    return 0


def basis_for(system, size="main"):
    """Basis set name. size: 'small' | 'main' | 'large'. SN2 systems get the diffuse (D) version.
    BASIS_OVERRIDE can replace the sizes (used only for the quick sandbox test)."""
    base = {"small": "def2-SVP", "main": "def2-TZVP", "large": "def2-QZVPP"}[size]
    if BASIS_OVERRIDE.get(size):
        return BASIS_OVERRIDE[size]
    return base + "D" if SYSTEMS[system]["sn2"] else base


BASIS_OVERRIDE = {}


def _build_start_geometries():
    """Start geometries (Angstrom). Minima are standard experimental-like values; transition-state guesses were refined with GFN2-xTB
    so that geomeTRIC starts close.  Nothing here is a result: C1 re-optimises everything it uses."""
    G = {}
    # --- system 1: HCN -> HNC (linear minima, triangular TS)
    G[("hcn", "HCN")] = [_pt("C", 0, 0, 0), _pt("N", 0, 0, 1.156), _pt("H", 0, 0, -1.065)]
    G[("hcn", "HNC")] = [_pt("C", 0, 0, 0), _pt("N", 0, 0, 1.169), _pt("H", 0, 0, 1.169 + 0.995)]
    th = math.radians(72.0)
    G[("hcn", "TS")] = [_pt("C", 0, 0, 0), _pt("N", 0, 0, 1.18), _pt("H", 1.19 * math.sin(th), 0, 1.19 * math.cos(th))]
    # --- system 2: NH3 inversion
    G[("nh3", "NH3")] = nh3_geometry(22.0)
    G[("nh3", "TS")] = nh3_geometry(0.0, 0.995)
    # --- system 3: H2CO -> H2 + CO
    G[("h2co", "H2CO")] = [_pt("C", 0, 0, 0), _pt("O", 0, 0, 1.205), _pt("H", 0.9435, 0, -0.585), _pt("H", -0.9435, 0, -0.585)]
    # TS found by DF-B3LYP/def2-SVP from three different guesses (all converged to this structure): C-H 1.10, C...H 1.68, H...H 1.31
    G[("h2co", "TS")] = [_pt("C", -0.132, 0, -0.014), _pt("O", -0.228, 0, 1.146), _pt("H", 1.466, 0, -0.542), _pt("H", 0.258, 0, -1.048)]
    G[("h2co", "H2")] = [_pt("H", 0, 0, -0.3705), _pt("H", 0, 0, 0.3705)]
    G[("h2co", "CO")] = [_pt("C", 0, 0, 0), _pt("O", 0, 0, 1.128)]
    # --- system 4: CH3NC -> CH3CN
    G[("ch3nc", "CH3NC")] = [_pt("C", 0, 0, 0), _pt("N", 0, 0, 1.424), _pt("C", 0, 0, 2.594)] + _ch3((0, 0, 0), (0, 0, 1), 1.09, 109.5)
    G[("ch3nc", "CH3CN")] = [_pt("C", 0, 0, 0), _pt("C", 0, 0, 1.46), _pt("N", 0, 0, 2.62)] + _ch3((0, 0, 0), (0, 0, 1), 1.09, 109.5)
    mid = np.array([0, 0, 0.595])
    cme = np.array([1.80, 0, 0.595])
    G[("ch3nc", "TS")] = [_pt("C", 0, 0, 0), _pt("N", 0, 0, 1.19), ("C", tuple(cme))] + _ch3(cme, mid - cme, 1.09, 100.0)
    # --- systems 5-8: SN2 (X- + CH3Y).  Backside geometry: X ... C - Y on the z axis.
    def ch3x(sym, d, r=1.09, sgn=+1, c0=0.0):
        """CH3-sym with carbon at z=c0 and sym on +/- z at distance d."""
        c = np.array([0, 0, c0]); ax = np.array([0, 0, sgn])
        return [_pt("C", *c), _pt(sym, 0, 0, c0 + sgn * d)] + _ch3(c, ax, r, 108.5)
    G[("sn2_ff", "F-")] = [_pt("F", 0, 0, 0)]
    G[("sn2_ff", "CH3F")] = ch3x("F", 1.38)
    G[("sn2_ff", "RC")] = ch3x("F", 1.43, sgn=-1)[:2] + [_pt("F", 0, 0, 2.55)] + ch3x("F", 1.43, sgn=-1)[2:]
    G[("sn2_ff", "TS")] = [_pt("C", 0, 0, 0), _pt("F", 0, 0, 1.83), _pt("F", 0, 0, -1.83)] + [_pt("H", 1.07 * math.cos(a), 1.07 * math.sin(a), 0) for a in (0, 2.0944, 4.1888)]
    G[("sn2_clcl", "Cl-")] = [_pt("Cl", 0, 0, 0)]
    G[("sn2_clcl", "CH3Cl")] = ch3x("Cl", 1.785)
    G[("sn2_clcl", "RC")] = ch3x("Cl", 1.82, sgn=-1)[:2] + [_pt("Cl", 0, 0, 3.15)] + ch3x("Cl", 1.82, sgn=-1)[2:]
    G[("sn2_clcl", "TS")] = [_pt("C", 0, 0, 0), _pt("Cl", 0, 0, 2.33), _pt("Cl", 0, 0, -2.33)] + [_pt("H", 1.07 * math.cos(a), 1.07 * math.sin(a), 0) for a in (0, 2.0944, 4.1888)]
    # F- + CH3Cl -> CH3F + Cl-
    G[("sn2_fcl", "F-")] = [_pt("F", 0, 0, 0)]
    G[("sn2_fcl", "CH3Cl")] = ch3x("Cl", 1.785)
    G[("sn2_fcl", "RC")] = ch3x("Cl", 1.85, sgn=-1)[:2] + [_pt("F", 0, 0, 2.6)] + ch3x("Cl", 1.85, sgn=-1)[2:]
    G[("sn2_fcl", "TS")] = [_pt("C", 0, 0, 0), _pt("F", 0, 0, 1.85), _pt("Cl", 0, 0, -2.20)] + [_pt("H", 1.07 * math.cos(a), 1.07 * math.sin(a), 0.03) for a in (0, 2.0944, 4.1888)]
    G[("sn2_fcl", "PC")] = ch3x("F", 1.43, sgn=-1)[:2] + [_pt("Cl", 0, 0, 3.15)] + ch3x("F", 1.43, sgn=-1)[2:]
    G[("sn2_fcl", "CH3F")] = ch3x("F", 1.38)
    G[("sn2_fcl", "Cl-")] = [_pt("Cl", 0, 0, 0)]
    # OH- + CH3F -> CH3OH + F-
    G[("sn2_ohf", "OH-")] = [_pt("O", 0, 0, 0), _pt("H", 0, 0, 0.96)]
    G[("sn2_ohf", "CH3F")] = ch3x("F", 1.38)
    # reactant complex: hydrogen-bonded HO- ... H-CH2F (O near one methyl H)
    cf = ch3x("F", 1.38)
    hb = np.array(cf[2][1]); u = _unit(hb - np.array(cf[0][1]))
    o = hb + 1.9 * u
    G[("sn2_ohf", "RC")] = cf + [_pt("O", *o), _pt("H", o[0] + 0.96 * u[0], o[1] + 0.96 * u[1], o[2] + 0.96 * u[2])]
    # TS: HO...C...F nearly collinear, O-H bent away
    G[("sn2_ohf", "TS")] = [_pt("C", 0, 0, 0), _pt("F", 0, 0, -1.90), _pt("O", 0, 0, 1.95), _pt("H", 0.90, 0, 2.30)] + [_pt("H", 1.07 * math.cos(a), 1.07 * math.sin(a), 0.03) for a in (0, 2.0944, 4.1888)]
    G[("sn2_ohf", "CH3OH")] = [_pt("C", 0, 0, 0), _pt("O", 0, 0, 1.42), _pt("H", 0.89, 0, 1.73)] + _ch3((0, 0, 0), (0, 0, 1), 1.09, 109.0)
    # product complex: CH3OH ... F-, F- bound to the O-H hydrogen
    G[("sn2_ohf", "PC")] = G[("sn2_ohf", "CH3OH")] + [_pt("F", 2.45, 0, 1.95)]
    G[("sn2_ohf", "F-")] = [_pt("F", 0, 0, 0)]
    # --- CH4
    G[("ch4", "CH4")] = ch4_geometry(1.087)
    return G


START_GEOMS = _build_start_geometries()


def system_species(system):
    return list(SPECIES_KIND[system])


def start_atoms(system, species):
    return START_GEOMS[(system, species)]


def n_electrons(atoms, charge=0):
    Z = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "Cl": 17}
    return sum(Z[s] for s, _ in atoms) - charge


def preview_systems():
    for k, v in SYSTEMS.items():
        print(f"{v['n']:>2} {k:9s} {v['name']:18s} class={v['cls']:12s} tier={v['tier']}  species={system_species(k)}")
