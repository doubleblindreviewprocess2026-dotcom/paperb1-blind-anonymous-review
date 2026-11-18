# ==============================================================================
# b1v5_gpu.py  --  calculation helpers and phase logic for the GPU notebook (G-1 to G4 and the optional phase).
# Pasted into the second code cell of the GPU notebook (after b1v5_common.py, b1v5_cpu.py and b1v5_v4_stored.py,
# which share the same namespace); it is not meant to be imported on its own.
#
# What is in here, in plain words:
#   1. detect the GPU and decide how much work to do
#   2. one function per functional (PBE, r2SCAN, B3LYP-D3(BJ), wB97M-V) and for Skala-1.1, on the GPU or the CPU,
#      with a ladder of retries when the SCF does not converge
#   3. a geometry optimiser (geomeTRIC) that takes energies and gradients from any of them
#   4. the phases: G-1 (rebuild missing CPU inputs), G0 (checks), G1 (all functionals), G2 (basis check),
#      G3 (xTB scout check), G4 (Skala-relaxed geometries), optional (paths of three more SN2 systems)
# Stand-in mode (STANDIN = True) swaps Skala and the large functionals for cheap PySCF functionals in tiny bases
# so that the whole notebook can be tested on a plain CPU.  Its numbers are NOT results.
#
# Skala and its D3 term (review 2026-11-15, rows 1 and 2) -- what was READ in the skala 2026.9 wheel (pip download, not run):
#   * SkalaKS / SkalaRKS (pyscf and gpu4pyscf) and the ASE calculator all default to with_dftd3=True: the PACKAGE DEFAULT IS D3 ON.
#   * The D3 variant is not hard-coded in the package: the dispersion model is built as DFTD3Dispersion(mol, xc.get_d3_settings()), where the
#     string "b3lyp5" is stored inside the weights file (extra file 'expected_d3_settings', identical in skala-1.1.fun and skala-1.1-rev1.fun).
#     pyscf-dispersion then uses its default version 'd3bj' (rational / Becke-Johnson damping) and atm=False (two-body only, no Axilrod-Teller-Muto
#     term).  In pyscf-dispersion 1.7.0 the parameter set "b3lyp5" gives the same energy as "b3lyp".  The D3 energy is added in energy_nuc().
#   * In skala 2026.9 the name "skala-1.1" resolves to the file skala-1.1-rev1.fun (the same trained model, retraced); "skala-1.1-rev0" is the
#     original skala-1.1.fun.  Both weights files are verified against a SHA-256 table inside the package.
#   * skala 2026.9 requires pyscf >= 2.11 and < 2.15 (the GPU notebook pins pyscf 2.12.1; the CPU notebook uses 2.15.0).
# WHAT THIS CODE ASSUMES: Skala-1.1 (headline) = D3 on, the package default; Skala-1.1-noD3 = the same SCF minus the stored D3 energy; the v4 gate
# (G0, I9) = D3 off, as in v4.  The D3 energy, the D3 settings string and the dispersion flag are written into the "extra" of every Skala row.
# Nothing here falls back silently: if with_dftd3 cannot be set, or no D3 energy comes back, the point fails loudly.
# ==============================================================================
import os, re, io, sys, json, csv, glob, math, time, shutil, tarfile, tempfile, subprocess, contextlib, importlib.util
import numpy as np
try:
    b1v5_common_loaded = COLUMNS  # noqa  (names exist when the notebook cells above ran)
except NameError:
    from b1v5_common import *     # noqa
try:
    b1v5_cpu_loaded = ccsdt       # noqa
except NameError:
    from b1v5_cpu import *        # noqa

# ------------------------------------------------------------------ settings (the notebook's settings cell overwrites these)
STANDIN = False               # True: cheap stand-in functionals and tiny bases on CPU (plumbing test, NOT results)
RETRY_FAILED = False
FORCE_RERUN = False
SYSTEMS_TO_RUN = list(SYSTEM_ORDER)
GPU_NAME, GPU_MEM_GB, GPU_CLASS = "", 0.0, "NONE"
NO_GPU = False
CALC_MAX_SEC = 1200           # one single SCF may not run longer than this (seconds); the point is then logged as failed
CONV_TOL = 1e-9               # SCF convergence (hartree), as in v4
GRID_LEVEL = 3                # DFT integration grid level, as in v4
GPU_QZ_MAX_NBF = 700          # G2: skip a large-basis calculation with more basis functions than this
CPU_FALLBACK_MAX_NBF = 220    # a failed GPU calculation is retried on the CPU only below this size
STOP_ON_G0_FAIL = True        # G0 reproduction of stored v4 energies off by more than 1e-4 hartree: stop the later phases
G4_SYSTEMS = ["hcn", "ch3nc", "sn2_ff"]
G4_OPTIONAL = ["nh3", "h2co"]
G4_DELTA = {"angle": 4.0, "distance": 0.04, "dihedral": 4.0}   # spacing of the three constrained points around the TS coordinate
G4_MAX_SHIFT = 3              # how often the three-point window may move if the top is not in the middle
G4_CONV = "GAU_LOOSE"         # geomeTRIC convergence set for the Skala optimisations
G4_MAXITER = 60
OPT_PATH_SYSTEMS = ["sn2_clcl", "sn2_fcl", "sn2_ohf"]
REGEN_PRIORITY = list(SYSTEM_ORDER)

CSV_REPRO = "phase10_repro.csv"
CSV_OPT = "phase16_optional.csv"
LV_SCANGPU = "own_scan_gpu"   # level name for the optional GPU relaxed scans (B3LYP without dispersion)
M_B3LYP_NODISP = "B3LYP-nodisp"   # used only to relax the optional GPU scans

# Method names and the order in which they are run when time is short (the last ones are cut first: PBE, then r2SCAN)
PRIORITY = [M_SKALA, "ANCHOR", M_B3LYP, M_R2SCAN, M_PBE]

# what G0 learns, saved in phase10_timing.json and reloaded by every later phase
GSTATE = dict(routes={}, anchor=M_WB97MV, gate_ok=None, skala_ok=None, skala_route="gpu", skala_mol_mode="cpu",
              grad_ok=None, grad_variant="plain", loaded=False)
SKALA_NAME = "skala-1.1"      # in skala 2026.9 this is the rev1 file; use "skala-1.1-rev0" to load the original weights file (see header)
DM_CACHE = {}


# ------------------------------------------------------------------ 1. GPU detection and scaling
def detect_gpu():
    """(name, memory in GB, class).  class is A100 (also H100), L4, T4, OTHER (any other GPU) or NONE."""
    try:
        out = subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader,nounits"],
                             capture_output=True, text=True, timeout=15).stdout.strip()
    except Exception:
        out = ""
    if not out:
        return "", 0.0, "NONE"
    name, _, mem = out.splitlines()[0].partition(",")
    name = name.strip()
    try:
        mem_gb = float(mem.strip()) / 1024.0
    except Exception:
        mem_gb = 0.0
    up = name.upper()
    if any(k in up for k in ("A100", "H100", "H200", "A800", "GH200")):
        cls = "A100"
    elif re.search(r"\bL4\b", up):
        cls = "L4"
    elif "T4" in up:
        cls = "T4"
    else:
        cls = "OTHER"
    return name, mem_gb, cls


# What each GPU class gets.  T4 (and unknown cards) trim the optional work; point counts of the plan are never reduced.
GPU_PLAN = {
    "A100": dict(optional=True, g4_extra=True, qz_max_nbf=900, calc_max_sec=1500),
    "L4": dict(optional=True, g4_extra=True, qz_max_nbf=700, calc_max_sec=1500),
    "T4": dict(optional=False, g4_extra=False, qz_max_nbf=550, calc_max_sec=1200),
    "OTHER": dict(optional=False, g4_extra=False, qz_max_nbf=550, calc_max_sec=1200),
    "NONE": dict(optional=False, g4_extra=False, qz_max_nbf=300, calc_max_sec=1200),
}


def gpu_available():
    """True when a GPU is present and gpu4pyscf can be imported."""
    return GPU_CLASS != "NONE" and importlib.util.find_spec("gpu4pyscf") is not None


def gpu_cleanup():
    """Free GPU memory between calculations (never raises)."""
    try:
        import gc
        gc.collect()
        try:
            import cupy as cp
            cp.get_default_memory_pool().free_all_blocks()
        except Exception:
            pass
        try:
            import torch
            torch.cuda.empty_cache()
        except Exception:
            pass
    except Exception:
        pass


class CalcTimeout(Exception):
    """One single calculation took longer than CALC_MAX_SEC."""


@contextlib.contextmanager
def calc_time_limit(sec):
    """Stops one calculation after `sec` seconds (only where the operating system allows it; otherwise no limit)."""
    import signal, threading
    if not sec or threading.current_thread() is not threading.main_thread() or not hasattr(signal, "setitimer"):
        yield
        return

    def _h(signum, frame):
        raise CalcTimeout(f"single calculation longer than {sec} s")
    old = signal.signal(signal.SIGALRM, _h)
    signal.setitimer(signal.ITIMER_REAL, float(sec))
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)


# ------------------------------------------------------------------ 2. functionals
# method name -> PySCF functional string (and dispersion).  The same strings are used on the GPU and the CPU.
XC_SPEC = {
    M_PBE: dict(xc="pbe"),
    M_R2SCAN: dict(xc="r2scan"),
    M_B3LYP: dict(xc="b3lyp", disp="d3bj"),
    M_WB97MV: dict(xc="wb97m_v"),
    M_WB97XV: dict(xc="wb97x_v"),
    M_B3LYP_NODISP: dict(xc="b3lyp"),
    "v4-pbe": dict(xc="pbe"),
    "v4-b3lyp": dict(xc="b3lyp"),                 # v4 used plain B3LYP (no dispersion)
}
# cheap stand-ins (STANDIN = True only): Skala-1.1 -> BLYP, wB97M-V -> wB97X (no non-local term)
STANDIN_SPEC = {
    M_SKALA: dict(xc="blyp", disp="d3bj", d3xc="b3lyp5"),        # the stand-in carries a D3 term too, so the with / without D3 plumbing is tested
    M_PBE: dict(xc="pbe"), M_R2SCAN: dict(xc="r2scan"), M_B3LYP: dict(xc="b3lyp", disp="d3bj"),
    M_WB97MV: dict(xc="wb97x"), M_WB97XV: dict(xc="wb97x"), M_B3LYP_NODISP: dict(xc="b3lyp"),
    "v4-pbe": dict(xc="pbe"), "v4-b3lyp": dict(xc="b3lyp"), "v4-skala": dict(xc="blyp"),
}
V4_BASE = {"v4-pbe": M_PBE, "v4-b3lyp": M_B3LYP, "v4-skala": M_SKALA}   # which method's route a v4 method follows
V4_KEY = {"v4-ccsd_t": "ccsd_t", "v4-pbe": "pbe", "v4-b3lyp": "b3lyp", "v4-skala": "skala"}


def is_skala(method):
    return method in (M_SKALA, "v4-skala") and not STANDIN


def spec_of(method):
    return (STANDIN_SPEC if STANDIN else XC_SPEC)[method]


def route_of(method):
    """'gpu', 'cpu', 'standin' or 'none' for a method (from G0; the GPU when G0 has not said anything)."""
    if STANDIN:
        return "standin"
    m = V4_BASE.get(method, method)
    if m == M_SKALA:
        return GSTATE.get("skala_route", "gpu") if GSTATE.get("skala_ok") is not False else "none"
    r = GSTATE["routes"].get(m)
    if r:
        return r
    return "gpu" if gpu_available() else "cpu"


def d3_energy(atoms, charge, basis, xc, version="d3bj"):
    """DFT-D3(BJ) dispersion energy (hartree), computed with the PySCF dispersion package on the CPU (two-body, atm=False: the same call that
    Skala uses; xc is only the NAME of the damping-parameter set, 'b3lyp5' for Skala and 'b3lyp' for B3LYP-D3(BJ) give the same energy)."""
    from pyscf.dispersion.dftd3 import DFTD3Dispersion
    mol = make_mol(atoms, charge, basis)
    return float(DFTD3Dispersion(mol, xc, version=version, atm=False).get_dispersion()["energy"])


def _to_np(x):
    return x.get() if hasattr(x, "get") else np.asarray(x)


def _build_skala(mol, route, df=False, variant="plain", d3=True):
    """Skala-1.1 through its PySCF bindings (GPU: skala.gpu4pyscf, CPU: skala.pyscf).
    d3=True: with its D3 term (the package default and the published functional; headline).  d3=False: dispersion off, ONLY for the v4 gate.
    The keyword with_dftd3 is always passed explicitly, so the setting never depends on a package default; an older signature without it
    raises an error (review row 2: no silent fallback that runs whatever the default is)."""
    os.environ.setdefault("TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD", "1")      # PyTorch >= 2.6 loading fix used in the v7 optimiser notebook
    if route == "gpu":
        from skala.gpu4pyscf import SkalaKS
    else:
        from skala.pyscf import SkalaKS
    m = mol.to_gpu() if (route == "gpu" and (variant == "v7" or GSTATE.get("skala_mol_mode") == "gpu")) else mol
    try:
        mf = SkalaKS(m, xc=SKALA_NAME, with_dftd3=bool(d3))
    except TypeError as ex:
        raise RuntimeError("this skala version does not accept with_dftd3: refusing to guess the dispersion setting "
                           f"(install the pinned skala 2026.9; original error: {ex!r})") from ex
    if d3 and getattr(mf, "with_dftd3", None) is None:
        raise RuntimeError("Skala was asked for its D3 term but built none (no 'expected_d3_settings' in the weights?): refusing to label the row 'with D3'")
    if variant == "v7":      # the settings of the tested v7 optimiser notebook (density fitting, full grid, no small-density cut)
        from pyscf import dft
        mf.grids.level = GRID_LEVEL
        mf.grids.prune = dft.gen_grid.nwchem_prune
        mf.small_rho_cutoff = 0.0
        mf.max_memory = 2000
        mf = mf.density_fit(auxbasis="def2-universal-jkfit")
        mf.init_guess = "minao"
    elif df:
        mf = mf.density_fit()
    return mf


def skala_d3_info(mf):
    """(D3 energy in hartree or None, name of the D3 parameter set or None) of a finished Skala calculation."""
    d3 = None
    try:
        d3 = mf.scf_summary.get("dispersion")
    except Exception:
        pass
    w = getattr(mf, "with_dftd3", None)
    if d3 is None and w:
        try:
            d3 = w.get_energy()
        except Exception:
            d3 = None
    return (None if d3 is None else float(d3)), (getattr(w, "_xc", None) if w else None)


def _build_mf(method, mol, route, df=False, variant="plain"):
    from pyscf import dft
    if is_skala(method):
        return _build_skala(mol, route, df, variant, d3=(method != "v4-skala"))
    mf = dft.RKS(mol, xc=spec_of(method)["xc"])
    mf.grids.level = GRID_LEVEL
    if route == "gpu":
        mf = mf.to_gpu()
    if df:
        mf = mf.density_fit()
    return mf


def _scf_once(method, atoms, charge, basis, route, dm0=None, df=False, level_shift=0.0, max_cycle=100, variant="plain"):
    """One SCF.  Returns dict(energy, converged, dm)."""
    mol = make_mol(atoms, charge, basis)
    mf = _build_mf(method, mol, route, df, variant)
    mf.conv_tol = CONV_TOL
    mf.max_cycle = max_cycle
    if level_shift:
        mf.level_shift = level_shift
    with calc_time_limit(CALC_MAX_SEC):
        e = float(mf.kernel(dm0=dm0)) if dm0 is not None else float(mf.kernel())
    conv = bool(getattr(mf, "converged", True))
    dm = None
    try:
        dm = _to_np(mf.make_rdm1())
    except Exception:
        pass
    d3, d3xc = skala_d3_info(mf) if is_skala(method) else (None, None)
    del mf
    return dict(energy=e, converged=conv, dm=dm, nbf=int(mol.nao), d3=d3, d3xc=d3xc)


def energy_of(method, atoms, charge, basis, dm_key=None, variant="plain", force_df=False):
    """Energy of one geometry with one method, with the retry ladder of plan section 5:
    warm start from the neighbouring density, fresh start, level shift, density fitting, then the CPU (not for Skala).
    force_df=True runs ONLY the density-fitted calculation (used to recompute the partner species of a barrier whose other species needed the
    density-fitting retry, review row 11).  Skala: the energy includes its D3 term; the D3 energy is returned in extra['d3_hartree'].
    Returns dict(energy, converged, extra).  Raises only if every attempt crashed."""
    route = route_of(method)
    if route == "none":
        raise RuntimeError(f"{method} is not available (see G0)")
    dm0 = DM_CACHE.get(dm_key) if dm_key else None
    nao = None
    if dm0 is not None:
        try:
            nao = make_mol(atoms, charge, basis).nao
            if dm0.shape != (nao, nao):
                dm0 = None
        except Exception:
            dm0 = None
    attempts = []
    if force_df:
        attempts.append(("density fitting", dict(df=True)))
        attempts.append(("density fitting + level shift", dict(df=True, level_shift=0.3, max_cycle=300)))
    else:
        if dm0 is not None:
            attempts.append(("warm start", dict(dm0=dm0)))
        attempts.append(("default", {}))
        attempts.append(("level shift", dict(level_shift=0.3, max_cycle=300)))
        if route == "gpu":
            attempts.append(("density fitting", dict(df=True)))
            if not is_skala(method):
                attempts.append(("cpu", dict(route="cpu")))
        if is_skala(method):
            # rescue for a Skala SCF that does not converge (seen on the Cl- + CH3Cl TS in the first Colab run): start from the converged PBE density
            attempts.append(("PBE start + level shift", dict(rescue=True, level_shift=0.3, max_cycle=300)))
            attempts.append(("PBE start + level shift 0.5 + density fitting", dict(rescue=True, df=True, level_shift=0.5, max_cycle=400)))
    best, last_err, tried = None, None, []
    for name, kw in attempts:
        kw = dict(kw)
        r_route = kw.pop("route", route)
        if kw.pop("rescue", False):
            try:
                pr = _scf_once(M_PBE, atoms, charge, basis, r_route)
                if not pr["converged"] or pr["dm"] is None:
                    tried.append(name + " skipped: PBE start did not converge")
                    continue
                kw["dm0"] = pr["dm"]
            except PhaseTimeout:
                raise
            except Exception as ex:
                tried.append(name + " skipped: PBE start crashed " + repr(ex)[:60])
                gpu_cleanup()
                continue
        if r_route == "cpu" and route == "gpu":
            if nao is None:
                nao = make_mol(atoms, charge, basis).nao
            if nao > CPU_FALLBACK_MAX_NBF:
                continue
        check_time()
        try:
            res = _scf_once(method, atoms, charge, basis, r_route, variant=variant, **kw)
            tried.append(name)
            if res["converged"]:
                best = dict(res, attempt=name, route=r_route, df=bool(kw.get("df")))
                break
            if best is None:
                best = dict(res, attempt=name + " (not converged)", route=r_route, df=bool(kw.get("df")))
        except PhaseTimeout:
            raise
        except Exception as ex:
            last_err = ex
            tried.append(name + " crashed: " + repr(ex)[:80])
            gpu_cleanup()
    gpu_cleanup()
    if best is None:
        raise last_err if last_err else RuntimeError("no SCF attempt possible")
    e = best["energy"]
    sp = {} if is_skala(method) else spec_of(method)
    d3_hartree, d3_set = None, None
    if is_skala(method):
        if method != "v4-skala":                                  # D3 on: the energy already contains it; it must have come back
            if best.get("d3") is None:
                raise RuntimeError("Skala returned no D3 energy although with_dftd3=True: refusing to store the row as 'with D3'")
            d3_hartree, d3_set = best["d3"], best.get("d3xc")
    elif sp.get("disp") and np.isfinite(e):
        d3_hartree = d3_energy(atoms, charge, basis, sp.get("d3xc", sp["xc"]), sp["disp"])
        d3_set = sp.get("d3xc", sp["xc"])
        e += d3_hartree
    if dm_key and best["converged"] and best["dm"] is not None:
        DM_CACHE[dm_key] = best["dm"]
    extra = dict(route=best["route"], attempt=best["attempt"], nbf=best["nbf"], df=best["df"])
    if d3_hartree is not None:
        extra.update(d3_hartree=float(d3_hartree), d3_set=d3_set, d3_version="d3bj, two-body (atm=False)")
    if method in (M_SKALA, "v4-skala"):
        extra.update(with_dftd3=bool(method != "v4-skala"), skala_name=SKALA_NAME, weights=GSTATE.get("weights_short", ""))
    if len(tried) > 1:
        extra["tried"] = tried
    if STANDIN:
        extra["standin"] = True
    return dict(energy=float(e), converged=bool(best["converged"]), extra=extra)


def reset_dm(prefix=None):
    """Forget stored densities (call at the start of every new path)."""
    for k in list(DM_CACHE):
        if prefix is None or k[0] == prefix:
            DM_CACHE.pop(k, None)


# ------------------------------------------------------------------ energies and gradients for the optimiser
class EG:
    """Energy and gradient of one method for one molecule (charge and basis fixed); used by geomeTRIC.
    eg(atoms) -> (energy in hartree, gradient in hartree/bohr as a flat array); energy(atoms) -> dict(energy, converged)."""

    def __init__(self, method, charge, basis, variant="plain"):
        self.method, self.charge, self.basis, self.variant = method, charge, basis, variant
        self.dm0, self.n_grad, self.n_energy = None, 0, 0

    def _run(self, atoms, want_grad):
        route = route_of(self.method)
        if route == "none":
            raise RuntimeError(f"{self.method} is not available")
        mol = make_mol(atoms, self.charge, self.basis)
        last = None
        for kw in ({"dm0": self.dm0}, {}, {"level_shift": 0.3, "max_cycle": 300}):
            if kw.get("dm0", 0) is None:
                continue
            mf = _build_mf(self.method, mol, route, False, self.variant)
            mf.conv_tol = CONV_TOL * (0.1 if want_grad else 1.0)
            mf.max_cycle = kw.get("max_cycle", 100)
            if kw.get("level_shift"):
                mf.level_shift = kw["level_shift"]
            with calc_time_limit(CALC_MAX_SEC):
                e = float(mf.kernel(dm0=kw["dm0"]) if "dm0" in kw else mf.kernel())
            if getattr(mf, "converged", True):
                try:
                    self.dm0 = _to_np(mf.make_rdm1())
                except Exception:
                    self.dm0 = None
                g = None
                if want_grad:
                    with calc_time_limit(CALC_MAX_SEC):
                        g = np.asarray(_to_np(mf.nuc_grad_method().kernel()), float).ravel()
                gpu_cleanup()
                return e, g
            last = e
            gpu_cleanup()
        raise RuntimeError(f"SCF did not converge in the gradient engine (last energy {last})")

    def energy(self, atoms):
        check_time()
        self.n_energy += 1
        e, _ = self._run(atoms, False)
        return dict(energy=e, converged=True)

    def eg(self, atoms):
        check_time()
        self.n_grad += 1
        e, g = self._run(atoms, True)
        return e, g


def opt_custom(atoms, charge, eg, constraints=None, convergence_set="GAU_LOOSE", maxiter=60):
    """Geometry optimisation (minimum) with geomeTRIC, energies and gradients from `eg(atoms)`.
    Returns dict(atoms, converged, n_steps, error).  When the optimiser gives up, the last geometry it tried is returned."""
    if len(atoms) == 1:
        return dict(atoms=atoms, converged=True, n_steps=0, error="")
    import geometric.optimize as go
    import geometric.engine as ge, geometric.molecule as gm
    quiet_geometric()

    class _Engine(ge.Engine):
        def __init__(self):
            M = gm.Molecule()
            M.elem = [s for s, _ in atoms]
            M.xyzs = [np.array([c for _, c in atoms])]
            super().__init__(M)
            self.n, self.last = 0, None

        def calc_new(self, coords, dirname):
            check_time()                           # lets the phase watchdog stop an optimisation between two gradients
            self.n += 1
            at = [(s, tuple(c * BOHR)) for (s, _), c in zip(atoms, coords.reshape(-1, 3))]
            self.last = at
            e, g = eg(at)
            return {"energy": e, "gradient": np.asarray(g, float).ravel()}
    eng = _Engine()
    err = ""
    with tempfile.TemporaryDirectory() as d:
        kw = dict(customengine=eng, input=os.path.join(d, "x"), maxiter=maxiter, convergence_set=convergence_set, transition=False)
        if constraints:
            cp = os.path.join(d, "cons.txt")
            open(cp, "w").write(constraint_text(constraints))
            kw["constraints"] = cp
        try:
            M = go.run_optimizer(**kw)
            conv = True
        except PhaseTimeout:
            raise
        except Exception as ex:
            conv, M, err = False, None, repr(ex)[:160]
    if M is not None:
        new = [(s, tuple(map(float, c))) for s, c in zip(M.elem, M.xyzs[-1])]
    else:
        new = eng.last or atoms
    return dict(atoms=new, converged=conv, n_steps=eng.n, error=err)


def fd_gradient_check(method, atoms, charge, basis, variant="plain", h=0.0025, n_comp=2):
    """Analytic gradient against a central finite difference of the energy (same engine), on the two largest components.
    Returns dict(ratio list, net_gradient, ok, status).  ratio = analytic / finite difference."""
    eng = EG(method, charge, basis, variant)
    e0, g = eng.eg(atoms)
    g = g.reshape(-1, 3)
    order = np.argsort(-np.abs(g).ravel())
    picks, seen_atoms = [], set()
    for k in order:                                    # largest components, on different atoms when possible
        i, c = divmod(int(k), 3)
        if i in seen_atoms and len(seen_atoms) < len(atoms):
            continue
        picks.append((i, c)); seen_atoms.add(i)
        if len(picks) == n_comp:
            break
    ratios, rows = [], []
    for i, c in picks:
        ep, em = [], []
        for sgn, store in ((+1, ep), (-1, em)):
            at = [(s, tuple(np.array(x) + (sgn * h if (j == i) else 0.0) * np.eye(3)[c])) for j, (s, x) in enumerate(atoms)]
            store.append(EG(method, charge, basis, variant).energy(at)["energy"])
        fd = (ep[0] - em[0]) / (2 * h / BOHR)          # hartree per bohr
        an = float(g[i, c])
        r = an / fd if abs(fd) > 1e-7 else float("nan")
        ratios.append(r); rows.append(dict(atom=int(i), xyz="xyz"[c], analytic=an, finite_difference=float(fd), ratio=r))
    net = float(np.linalg.norm(g.sum(axis=0)))
    worst = max([abs(r - 1.0) for r in ratios if np.isfinite(r)] or [9.9])
    status = "PASS" if (worst <= 0.05 and net <= 1e-3) else ("WARN" if (worst <= 0.10 and net <= 3e-3) else "FAIL")
    return dict(components=rows, net_gradient=net, worst_ratio_deviation=float(worst), status=status, ok=status in ("PASS", "WARN"),
                energy=float(e0), variant=variant)


# ------------------------------------------------------------------ Skala weights on Drive
def hf_hub_dir():
    home = os.environ.get("HF_HOME") or os.path.join(os.path.expanduser("~"), ".cache", "huggingface")
    return os.environ.get("HF_HUB_CACHE") or os.path.join(home, "hub")


def restore_weights_cache():
    """Unpack the Skala weights from Drive (saved by an earlier session) when this session has none yet."""
    tgz, hub = path("skala_weights_cache.tgz"), hf_hub_dir()
    if os.path.exists(tgz) and not glob.glob(os.path.join(hub, "models--*kala*")):
        try:
            os.makedirs(hub, exist_ok=True)
            with tarfile.open(tgz) as t:
                t.extractall(hub)
            return "restored from Drive"
        except Exception as ex:
            return "could not restore: " + repr(ex)[:80]
    return "nothing to restore" if not os.path.exists(tgz) else "already present locally"


def save_weights_cache():
    """Copy the Skala weights (Hugging Face cache) to Drive once, so later sessions need no download."""
    tgz, hub = path("skala_weights_cache.tgz"), hf_hub_dir()
    dirs = glob.glob(os.path.join(hub, "models--*kala*"))
    if not dirs:
        return "no Hugging Face cache for Skala found (weights may ship inside the package)"
    if os.path.exists(tgz):
        return "already on Drive"
    try:
        with tarfile.open(tgz + ".tmp", "w:gz") as t:
            for d in dirs:
                t.add(d, arcname=os.path.basename(d))
        os.replace(tgz + ".tmp", tgz)
        return f"saved {os.path.getsize(tgz) / 1e6:.0f} MB to Drive"
    except Exception as ex:
        return "could not save: " + repr(ex)[:80]


def skala_weights_info():
    """Review row 10: which weights file is used.  Returns SHA-256 of the file, its name, the Hugging Face revision (the snapshot folder name
    of the download), the SHA-256 the skala package expects, and the D3 parameter-set name stored in the file ('expected_d3_settings')."""
    import hashlib, zipfile
    info = dict(skala_name=SKALA_NAME, standin=STANDIN)
    if STANDIN:
        info["note"] = "stand-in: no Skala weights are used"
        return info
    p, expected = None, None
    try:
        from skala.functional import resolve_functional_artifact
        art = resolve_functional_artifact(SKALA_NAME)
        p, expected = str(art.path), getattr(art, "expected_hash", None)
    except Exception as ex:
        info["resolve_error"] = repr(ex)[:160]
        cand = sorted(glob.glob(os.path.join(hf_hub_dir(), "models--microsoft--skala-*", "snapshots", "*", "*.fun")), key=os.path.getmtime)
        p = cand[-1] if cand else None
    if not p or not os.path.exists(p):
        info["note"] = "weights file not found"
        return info
    info["file"] = os.path.basename(os.path.realpath(p))
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    info["sha256"] = h.hexdigest()
    info["package_expected_sha256"] = expected
    info["hash_matches_package_table"] = (None if expected is None else expected == info["sha256"])
    m = re.search(r"snapshots[/\\]([0-9a-f]{40})", p)
    info["huggingface_revision"] = m.group(1) if m else "unknown (file not in a Hugging Face snapshot folder)"
    info["huggingface_repo"] = "microsoft/" + SKALA_NAME.partition("-rev")[0]
    try:
        with zipfile.ZipFile(p) as z:
            nm = [n for n in z.namelist() if n.endswith("extra/expected_d3_settings")]
            info["expected_d3_settings"] = z.read(nm[0]).decode("utf-8").strip() if nm else "not found in file"
    except Exception as ex:
        info["expected_d3_settings"] = "could not read: " + repr(ex)[:80]
    GSTATE["weights_short"] = info["sha256"][:12]
    return info


# ------------------------------------------------------------------ saved state of G0
def save_gstate(**extra):
    d = dict(routes=GSTATE["routes"], anchor=GSTATE["anchor"], gate_ok=GSTATE["gate_ok"], skala_ok=GSTATE["skala_ok"],
             skala_route=GSTATE["skala_route"], skala_mol_mode=GSTATE["skala_mol_mode"], grad_ok=GSTATE["grad_ok"],
             grad_variant=GSTATE["grad_variant"], weights_short=GSTATE.get("weights_short", ""))
    old = read_json("phase10_timing.json", {}) or {}
    old.update(d)
    old.update(extra)
    write_json("phase10_timing.json", old)


def run_gphase(name, title, fn, max_min, force=None, blocked=None):
    """Run one phase of the notebook: skip it with a clear message when it must not start, otherwise guard it (errors are logged, never raised)."""
    why = blocked if blocked is not None else phase_blocked(name)
    if why:
        skip_phase(name, title, why)
        return
    with phase(name, title, max_min=max_min, force=(FORCE_RERUN if force is None else force)) as ph:
        if not ph.skip:
            fn()


def load_gstate():
    """Read what G0 learned (routes per functional, anchor, gates).  Silent when G0 has not run."""
    st = read_json("phase10_timing.json", None)
    if not st:
        return False
    for k in ("routes", "anchor", "gate_ok", "skala_ok", "skala_route", "skala_mol_mode", "grad_ok", "grad_variant", "weights_short"):
        if k in st and st[k] is not None:
            GSTATE[k] = st[k]
    GSTATE["loaded"] = True
    return True


def anchor_method():
    return GSTATE.get("anchor") or M_WB97MV


def ladder_methods(stage="stationary"):
    """The functionals to run, in priority order.  stage 'path' leaves out CPU-routed ones (CPU is for stationary points only)."""
    out = []
    for m in PRIORITY:
        m = anchor_method() if m == "ANCHOR" else m
        r = route_of(m)
        if r == "none":
            continue
        if stage == "path" and r == "cpu":
            continue
        out.append(m)
    return out


def phase_blocked(name):
    """Reason (text) why a GPU phase must not start, or None.  Used before `with phase(...)`."""
    if name == "phase10":                                    # G0 always runs: with no GPU it is the diagnostics-only phase
        return None
    if NO_GPU and not STANDIN:
        return "no GPU detected: only the G0 diagnostics run (Runtime > Change runtime type > GPU, then Run all again)"
    load_gstate()
    if name in ("phase11", "phase11p", "phase12", "phase13", "phase14", "phase16") and GSTATE.get("gate_ok") is False and STOP_ON_G0_FAIL:
        return "G0 reproduction of the stored v4 energies failed (difference over 1e-4 hartree): old and new numbers must not be mixed. Fix the install, delete DONE_phase10 and run again (or set STOP_ON_G0_FAIL = False)"
    if name == "phase14":
        if GSTATE.get("skala_ok") is False:
            return "Skala is not usable (see G0): G4 skipped"
        if GSTATE.get("grad_ok") is not True:
            return "the Skala gradient check of G0 did not pass (or did not run): G4 skipped and the limit is stated in the summary"
    return None


# ------------------------------------------------------------------ geometry helpers shared by the phases
def g_save_path(system, pts, label="", suffix="path", level=LV_SCAN):
    p = path("geoms", f"{system}_{suffix}.json")
    d = dict(system=system, level=level, basis=basis_for(system, "main"), label=label, saved=now(),
             points=[dict(index=int(q["index"]), coordinate=float(q["coordinate"]), atoms=[[s, list(map(float, c))] for s, c in q["atoms"]],
                          energy=(None if q.get("energy") is None else float(q["energy"])), converged=bool(q.get("converged", True))) for q in pts])
    json.dump(d, open(p + ".tmp", "w"), indent=1)
    os.replace(p + ".tmp", p)


def g_load_path(system, suffix="path"):
    p = path("geoms", f"{system}_{suffix}.json")
    if not os.path.exists(p):
        return []
    try:
        return [dict(q, atoms=[(s, tuple(c)) for s, c in q["atoms"]]) for q in json.load(open(p))["points"]]
    except Exception:
        return []


def g_expected_points(system):
    if system == "ch4":
        return len(CH4_STRETCH)
    n = SCANS[system]["n"]
    return len(sorted(set(range(0, n, SCAN_STRIDE)) | {n - 1}))


def g_scout_species(system):
    S = SYSTEMS[system]
    pmin = S["pc"] or (S["products"][0] if len(S["products"]) == 1 else None)
    return [S["rc"] or S["reactants"][0], "TS", pmin]


def g_levels(system, sp):
    """Geometry levels stored for a species, preferred first."""
    return [lv for lv in LEVEL_PREFERENCE if geometry_of(system, sp, lv)]


# ------------------------------------------------------------------ G-1: rebuild what the CPU notebook did not leave (same helpers, budget-limited)
def g_missing_inputs(systems):
    miss = dict(stationary=[], paths=[], xtb=[], tiers=not os.path.exists(path("tiers.json")))
    for s in systems:
        for sp, kind in SPECIES_KIND[s].items():
            if geometry_of(s, sp, LV_OWN) is None and geometry_of(s, sp, LV_SCANPEAK) is None:      # our own B3LYP geometry is the headline protocol
                miss["stationary"].append((s, sp, kind))
        if s == "ch4" or (s in SCANS and SCANS[s].get("dft", True)):
            if len(g_load_path(s)) < g_expected_points(s):
                miss["paths"].append(s)
        if s in REACTIONS:
            names = [geom_name(s, x, LV_XTB) for x in g_scout_species(s) if x]
            if any(load_geom(n) is None for n in names):
                miss["xtb"].append(s)
    return miss


def g_regen_stationary(system, sp, kind):
    """Same as C1 of the CPU notebook: optimise one species (minimum or transition state), frequencies, CSV row."""
    ch, basis = species_charge(system, sp), basis_for(system, "main")
    name = geom_name(system, sp, LV_OWN)
    g = load_geom(name)
    if g is None:
        t0 = time.time()
        at = start_atoms(system, sp)
        if len(at) > 1:
            r = opt_dft(at, ch, basis, transition=(kind == "ts"), hessian=(kind == "ts"))
            at, conv = r["atoms"], r["converged"]
        else:
            conv = True
        save_geom(name, at, ch, species=sp, kind=kind, level=LV_OWN, opt_converged=bool(conv), basis=basis, seconds=round(time.time() - t0, 1), regenerated_by="GPU notebook G-1")
        g = load_geom(name)

    def calc():
        at, a, retried = refine_stationary(g["atoms"], ch, basis, kind)
        if retried and at is not g["atoms"]:
            save_geom(name, at, ch, species=sp, kind=kind, level=LV_OWN, opt_converged=True, basis=basis, retried=True)
        a["converged"] = bool(a["converged"] and g.get("opt_converged", True))
        a.pop("imag_modes", None)
        a["extra"]["retried"] = bool(retried)
        return a
    run_point(PHASE_CSV["phase01"], system, sp, LV_OWN, M_DFB3LYP, basis, calc, retry_failed=RETRY_FAILED)


def g_regen_scan(system):
    """Same as C2 of the CPU notebook: one relaxed scan (geometry file plus one CSV row per point)."""
    if system == "ch4":
        g_save_path("ch4", [dict(index=i, coordinate=float(s), atoms=ch4_geometry(float(s)), energy=None, converged=True) for i, s in enumerate(CH4_STRETCH)], "C-H stretch (A)")
        return
    sp = SCANS[system]
    g0 = geometry_of(system, sp["start"], LV_OWN) or best_geom(system, sp["start"])[1]
    if g0 is None:
        raise RuntimeError(f"missing start geometry {sp['start']} for the scan of {system}")
    ch, basis = species_charge(system, sp["start"]), basis_for(system, "main")
    store = g_load_path(system)
    try:
        pts = scan_points_dft(system, g0["atoms"], ch, basis, store, lambda pts, system=system, sp=sp: g_save_path(system, pts, sp["label"]))
    finally:
        for q in g_load_path(system):
            run_point(PHASE_CSV["phase02"], system, f"p{q['index']:02d}", LV_SCAN, M_DFB3LYP, basis_for(system, "main"),
                      lambda q=q: dict(energy=q["energy"], converged=q["converged"]), coordinate=f"{q['coordinate']:.4f}", quiet=True)
    if best_geom(system, "TS")[1] is None:           # TS search failed: the top of the scan becomes the TS (plan section 5)
        k = peak_of(pts)
        if k is None:
            raise RuntimeError("no peak along the scan and no TS")
        save_geom(geom_name(system, "TS", LV_SCANPEAK), pts[k]["atoms"], ch, species="TS", kind="ts", level=LV_SCANPEAK, from_scan_point=k)


def g_regen_xtb(system):
    """Same as the geometry part of C4: GFN2-xTB minima, scan peak and product (no CCSD(T) here)."""
    csv_ = PHASE_CSV["phase04"]
    rmin, _, pmin = g_scout_species(system)
    ch = species_charge(system, rmin)
    names = {sp: geom_name(system, sp, LV_XTB) for sp in (rmin, "TS", pmin) if sp}
    lv, g = best_geom(system, rmin)
    a0 = g["atoms"] if g else start_atoms(system, rmin)
    r = opt_xtb(a0, ch)
    save_geom(names[rmin], r["atoms"], ch, species=rmin, level=LV_XTB, converged=r["converged"])
    if pmin:
        lv, g = best_geom(system, pmin)
        r2 = opt_xtb(g["atoms"] if g else start_atoms(system, pmin), species_charge(system, pmin))
        save_geom(names[pmin], r2["atoms"], species_charge(system, pmin), species=pmin, level=LV_XTB, converged=r2["converged"])
    sp_ = SCANS[system]
    if sp_["start"] == rmin:
        g0 = r["atoms"]
    else:
        lv0, gg = best_geom(system, sp_["start"])
        g0 = gg["atoms"] if gg else start_atoms(system, sp_["start"])
    pts = scan_points_xtb(system, g0, species_charge(system, sp_["start"]))
    k = peak_of(pts)
    if k is None:
        save_geom(names["TS"], pts[len(pts) // 2]["atoms"], ch, species="TS", level=LV_XTB, no_peak=True)
        append_row(csv_, make_row(system, "TS", LV_XTB, M_XTB, "GFN2", float("nan"), False, extra={"no_peak": True, "scan_max_index": int(np.nanargmax([q["energy"] for q in pts]))}))
    else:
        save_geom(names["TS"], pts[k]["atoms"], ch, species="TS", level=LV_XTB, scan_point=k)
    for sp, nm in names.items():
        gg = load_geom(nm)
        if gg and not gg.get("no_peak"):
            run_point(csv_, system, sp, LV_XTB, M_XTB, "GFN2", lambda gg=gg, sp=sp: xtb_energy(gg["atoms"], species_charge(system, sp)), retry_failed=RETRY_FAILED)


def g_regen():
    """Phase G-1: find what is missing, rebuild it with the CPU helpers, and say so."""
    systems = [s for s in REGEN_PRIORITY if s in SYSTEMS_TO_RUN]
    miss = g_missing_inputs(systems)
    n = len(miss["stationary"]) + len(miss["paths"]) + len(miss["xtb"])
    print(f"  missing from the CPU notebook: {len(miss['stationary'])} stationary geometries, {len(miss['paths'])} paths, {len(miss['xtb'])} xTB scout sets, tiers.json {'missing' if miss['tiers'] else 'present'}")
    lines = []
    if miss["tiers"]:
        write_json("tiers.json", {k: v["tier"] for k, v in SYSTEMS.items()})
        lines.append("tiers.json was missing: plan defaults written (the CPU notebook C0 timing rule was not applied)")
    if n == 0:
        lines.append("nothing missing: CPU inputs found")
    for s, sp, kind in miss["stationary"]:
        with item_guard("phase09", f"{s}:{sp}"):
            print(f"  regenerating geometry {s} {sp} ({kind})")
            g_regen_stationary(s, sp, kind)
            lines.append(f"regenerated stationary geometry {s}:{sp}")
    for s in miss["paths"]:
        with item_guard("phase09", f"{s}:path"):
            print(f"  regenerating relaxed scan {s}")
            g_regen_scan(s)
            lines.append(f"regenerated path {s}")
    for s in miss["xtb"]:
        with item_guard("phase09", f"{s}:xtb"):
            print(f"  regenerating xTB scout geometries {s}")
            g_regen_xtb(s)
            lines.append(f"regenerated xTB scout geometries {s}")
    left = g_missing_inputs(systems)
    left_n = len(left["stationary"]) + len(left["paths"]) + len(left["xtb"])
    lines.append(f"still missing after G-1: {left_n} item(s)" + ("" if not left_n else " (" + ", ".join([f"{a}:{b}" for a, b, _ in left["stationary"]] + [x + ":path" for x in left["paths"]] + [x + ":xtb" for x in left["xtb"]])[:300] + ")"))
    write_json("phase09_regen.json", dict(missing_before=dict(stationary=[f"{a}:{b}" for a, b, _ in miss["stationary"]], paths=miss["paths"], xtb=miss["xtb"], tiers=miss["tiers"]),
                                           log=lines, finished=now()))
    for l in lines:
        print("  G-1:", l)
    return lines


# ------------------------------------------------------------------ G0
def g_v4_series(method, k_list):
    """The NB4 HCN path points (v4 geometry list) with one v4 method; rows go to phase10_repro.csv."""
    basis = "STO-3G" if STANDIN else "def2-TZVP"
    reset_dm("nb4")
    for k in k_list:
        th, r_ = NB4_HCN_GEOMS[k]
        at = hcn_geometry(th, r_)
        if method == "v4-ccsd_t":
            fn = lambda at=at: ccsdt(at, 0, basis, frozen_core=False, stability=False)
        else:
            fn = lambda at=at: energy_of(method, at, 0, basis, dm_key=("nb4", method, basis))
        run_point(CSV_REPRO, "hcn_nb4", f"nb4_{k:02d}", LV_NB4, method, basis, fn, coordinate=f"{th},{r_}", retry_failed=RETRY_FAILED, quiet=(len(k_list) > 1))


def g_gate():
    """Compare every rerun v4 energy with the stored value.  Returns dict(max_diff, per_method, status)."""
    per, worst_d = {}, 0.0
    for r in read_rows(CSV_REPRO):
        key = V4_KEY.get(r["method"])
        if not key or r["energy_hartree"] in ("", "nan"):
            continue
        k = int(r["point_id"].split("_")[-1])
        d = abs(float(r["energy_hartree"]) - NB4_HCN_E[key][k])
        per[r["method"]] = max(per.get(r["method"], 0.0), d)
        worst_d = max(worst_d, d)
    status = "no data" if not per else ("PASS" if worst_d <= 1e-5 else ("WARN" if worst_d <= 1e-4 else "FAIL"))
    return dict(max_diff=worst_d, per_method=per, status=status)


def g_probe_one(method, atoms, basis):
    """Does this functional run on the GPU, and does it agree with the CPU?  Returns a record; route = 'gpu' | 'cpu' | 'none'."""
    from pyscf import dft
    rec = dict(method=method, basis=basis, gpu_ok=False, cpu_ok=False, gpu_seconds=None, cpu_seconds=None, gpu_cpu_diff=None, error_gpu="", error_cpu="")
    sp = spec_of(method)
    e_cpu = None
    t0 = time.time()
    try:                                                   # CPU reference, native PySCF (with the dispersion switch for B3LYP)
        mf = dft.RKS(make_mol(atoms, 0, basis), xc=sp["xc"])
        mf.grids.level = GRID_LEVEL
        mf.conv_tol = 1e-9
        if sp.get("disp"):
            mf.disp = sp["disp"]
        with calc_time_limit(CALC_MAX_SEC):
            e_cpu = float(mf.kernel())
        rec["cpu_ok"] = bool(mf.converged) and np.isfinite(e_cpu)
        rec["cpu_seconds"] = round(time.time() - t0, 1)
        rec["e_cpu"] = e_cpu
    except Exception as ex:
        rec["error_cpu"] = repr(ex)[:150]
    if gpu_available() and not STANDIN:
        t0 = time.time()
        try:
            res = _scf_once(method, atoms, 0, basis, "gpu")
            e = res["energy"] + (d3_energy(atoms, 0, basis, sp["xc"], sp["disp"]) if sp.get("disp") else 0.0)
            rec["gpu_seconds"] = round(time.time() - t0, 1)
            rec["e_gpu"] = e
            if e_cpu is not None:
                rec["gpu_cpu_diff"] = abs(e - e_cpu)
                rec["gpu_ok"] = bool(res["converged"]) and abs(e - e_cpu) < 5e-5
            else:
                rec["gpu_ok"] = bool(res["converged"])
        except Exception as ex:
            rec["error_gpu"] = repr(ex)[:150]
        gpu_cleanup()
    else:
        rec["error_gpu"] = "no GPU / gpu4pyscf in this session" if not STANDIN else "stand-in: GPU not used"
    rec["route"] = "gpu" if rec["gpu_ok"] else ("cpu" if rec["cpu_ok"] else "none")
    return rec


def g_probe_functionals():
    """Plan section 5: which functionals work on the GPU.  Sets GSTATE['routes'] and GSTATE['anchor']."""
    atoms, basis = start_atoms("hcn", "HCN"), ("STO-3G" if STANDIN else "def2-SVP")
    recs = {}
    for m in (M_PBE, M_R2SCAN, M_B3LYP, M_WB97MV):
        with item_guard("phase10", f"probe {m}"):
            recs[m] = g_probe_one(m, atoms, basis)
            r = recs[m]
            print(f"  probe {m:11s}: GPU {'ok' if r['gpu_ok'] else 'NO'} ({r['gpu_seconds']} s, difference to CPU {r['gpu_cpu_diff']}), CPU {'ok' if r['cpu_ok'] else 'NO'} ({r['cpu_seconds']} s) -> route {r['route']}"
                  + (f"  [{r['error_gpu']}]" if r["error_gpu"] and not r["gpu_ok"] else ""))
    routes = {m: r["route"] for m, r in recs.items()}
    anchor = M_WB97MV
    if routes.get(M_WB97MV, "none") == "none":                # fallback of plan section 5: wB97X-V
        with item_guard("phase10", "probe wB97X-V"):
            recs[M_WB97XV] = g_probe_one(M_WB97XV, atoms, basis)
            routes[M_WB97XV] = recs[M_WB97XV]["route"]
            anchor = M_WB97XV if routes[M_WB97XV] != "none" else M_WB97MV
            print(f"  wB97M-V does not run: the anchor becomes {anchor} (route {routes[M_WB97XV]})")
    elif routes[M_WB97MV] == "cpu" and not STANDIN:
        print("  wB97M-V runs on the CPU only: it will be computed at stationary points only (plan section 5)")
    for m in (M_PBE, M_R2SCAN, M_B3LYP):
        if routes.get(m) == "cpu" and not STANDIN:
            print(f"  {m} runs on the CPU only: stationary points only (plan section 5)")
    if STANDIN:                                               # the stand-in runs everything on the CPU, paths included
        routes = {m: "standin" for m in list(routes) + [M_WB97XV]}
    GSTATE["routes"], GSTATE["anchor"] = routes, anchor
    return recs


def g_skala_check():
    """Does Skala run, on the GPU (else the CPU)?  Sets skala_route / skala_ok (preliminary) and returns a record."""
    atoms, basis = start_atoms("hcn", "HCN"), ("STO-3G" if STANDIN else "def2-SVP")
    rec = dict(basis=basis, tried=[])
    if STANDIN:
        GSTATE["skala_route"], GSTATE["skala_ok"] = "standin", True
        t0 = time.time()
        r = energy_of(M_SKALA, atoms, 0, basis)
        rec.update(route="standin", energy=r["energy"], seconds=round(time.time() - t0, 1), note="stand-in functional, not Skala")
        return rec
    for route, mode in (("gpu", "cpu"), ("gpu", "gpu"), ("cpu", "cpu")):
        if route == "gpu" and not gpu_available():
            continue
        GSTATE["skala_route"], GSTATE["skala_mol_mode"], GSTATE["skala_ok"] = route, mode, True
        t0 = time.time()
        try:
            res = _scf_once(M_SKALA, atoms, 0, basis, route)
            gpu_cleanup()
            if res["converged"] and np.isfinite(res["energy"]):
                rec.update(route=route, mol_mode=mode, energy=res["energy"], seconds=round(time.time() - t0, 1))
                print(f"  Skala runs: route {route} (molecule object: {mode}), E = {res['energy']:.8f}, {rec['seconds']} s")
                return rec
            rec["tried"].append(f"{route}/{mode}: not converged")
        except PhaseTimeout:
            raise
        except Exception as ex:
            rec["tried"].append(f"{route}/{mode}: {repr(ex)[:140]}")
            gpu_cleanup()
    GSTATE["skala_ok"] = False
    print("  Skala does NOT run in this session:", rec["tried"])
    return rec


def g_skala_gradient_check():
    """Plan G0: one Skala gradient against finite differences.  If the plain engine fails, the settings of the tested v7 optimiser notebook are tried."""
    basis = "STO-3G" if STANDIN else basis_for("hcn", "main")
    atoms = hcn_geometry(60.0, 1.40)            # a bent, off-equilibrium HCN with a large gradient
    out = dict(basis=basis, tests=[])
    for variant in (["plain"] if STANDIN else ["plain", "v7"]):
        with item_guard("phase10", f"gradient check {variant}"):
            t0 = time.time()
            r = fd_gradient_check(M_SKALA, atoms, 0, basis, variant=variant)
            r["seconds"] = round(time.time() - t0, 1)
            out["tests"].append(r)
            print(f"  gradient check ({variant}): {r['status']}, ratios {[round(c['ratio'], 4) for c in r['components']]}, net gradient {r['net_gradient']:.1e}")
            if r["ok"]:
                GSTATE["grad_ok"], GSTATE["grad_variant"] = True, variant
                out["chosen"] = variant
                return out
    GSTATE["grad_ok"] = False
    return out


def g0_run():
    """Phase G0 (everything in one place, every step guarded)."""
    info = dict(started=now(), hardware=hw(), versions=ver(), gpu_name=GPU_NAME, gpu_class=GPU_CLASS, gpu_mem_gb=round(GPU_MEM_GB, 1), standin=STANDIN)
    with item_guard("phase10", "weights cache"):
        info["weights_restore"] = restore_weights_cache()
        print("  Skala weights:", info["weights_restore"])
    if NO_GPU and not STANDIN:
        info["note"] = "no GPU: diagnostics only"
        save_gstate(**info)
        print("  ** NO GPU: G0 stops here after the diagnostics. **")
        return info
    with item_guard("phase10", "probe functionals"):
        info["probe"] = g_probe_functionals()
        save_gstate(**info)
    with item_guard("phase10", "Skala check"):
        info["skala"] = g_skala_check()
        save_gstate(**info)
    with item_guard("phase10", "weights cache save"):
        if GSTATE.get("skala_ok") and not STANDIN:
            info["weights_save"] = save_weights_cache()
            print("  Skala weights on Drive:", info["weights_save"])
    with item_guard("phase10", "weights hash"):                       # review row 10: SHA-256 of the weights and the Hugging Face revision
        info["skala_weights"] = skala_weights_info()
        w = info["skala_weights"]
        print("  Skala weights:", {k: w.get(k) for k in ("file", "sha256", "huggingface_revision", "hash_matches_package_table", "expected_d3_settings")})
        save_gstate(**info)
    # reproduce the first stored NB4 HCN point with all four v4 methods, then the other 11 points (DFT methods) for check I9
    for m in ("v4-ccsd_t", "v4-pbe", "v4-b3lyp", "v4-skala"):
        if m == "v4-skala" and not GSTATE.get("skala_ok"):
            continue
        with item_guard("phase10", f"repro {m}"):
            g_v4_series(m, [0])
    first = g_gate()
    print(f"  point 0 of NB4: largest |difference| to the stored v4 energy {first['max_diff']:.1e} hartree ({first['status']}); per method {first['per_method']}")
    for m in ("v4-pbe", "v4-b3lyp", "v4-skala"):
        if m == "v4-skala" and not GSTATE.get("skala_ok"):
            continue
        with item_guard("phase10", f"repro all points {m}"):
            g_v4_series(m, list(range(12)))
    gate = g_gate()
    info["gate"] = gate
    if STANDIN:
        GSTATE["gate_ok"] = True
        print("  stand-in: the gate is not applied (stand-in energies differ from the stored v4 energies by design)")
    elif gate["status"] == "no data":
        GSTATE["gate_ok"] = None
        print("  gate: no reproduced energies (all attempts failed); later phases run, but nothing was verified")
    else:
        GSTATE["gate_ok"] = gate["status"] != "FAIL"
        if gate["status"] == "WARN":
            print("  gate WARN: differences between 1e-5 and 1e-4 hartree: version drift, state it in the paper")
        if gate["status"] == "FAIL":
            print("  ** GATE FAILED: differences over 1e-4 hartree. The later phases will not start (STOP_ON_G0_FAIL). **")
    save_gstate(**info)
    if GSTATE.get("skala_ok"):
        with item_guard("phase10", "Skala gradient"):
            info["gradient"] = g_skala_gradient_check()
    else:
        GSTATE["grad_ok"] = False
    save_gstate(**info)
    return info


# ------------------------------------------------------------------ G1 / G2 / G3 building blocks
def ensure_nod3(csv_name, system, pid, level, basis, coordinate=""):
    """'Skala without D3' (review row 1): the SAME SCF as the headline Skala row minus the D3 energy stored in its extra.
    Writes the second row (method Skala-1.1-noD3) unless it exists; also repairs a run that was cut between the two rows.
    A failed headline row gives a failed noD3 row (so that the two variants always have the same points)."""
    if have(csv_name, system, pid, level, M_SKALA_NOD3, basis, RETRY_FAILED):
        return None
    src = [r for r in read_rows(csv_name) if row_key(r) == (system, str(pid), level, M_SKALA, basis)]
    if not src:
        return None
    src.sort(key=lambda r: r["converged"] in ("1", "True"))            # a repeated headline row (RETRY_FAILED) that converged is used, else the last failed one
    r = src[-1]
    x = extra_of(r)
    ok = r["converged"] in ("1", "True") and r["energy_hartree"] not in ("", "nan") and x.get("d3_hartree") is not None
    e = float(r["energy_hartree"]) - float(x["d3_hartree"]) if ok else float("nan")
    row = make_row(system, pid, level, M_SKALA_NOD3, basis, e, ok, coordinate or r["coordinate"], seconds=0.0,
                   extra=dict(x, derived_from=M_SKALA, same_scf=True, with_dftd3=False))
    append_row(csv_name, row)
    return row


def g_energy_row(csv_name, system, pid, level, method, basis, atoms, charge, coordinate="", dm_key=None):
    """One energy, saved at once (skipped when already in the CSV).  For Skala the second row 'Skala-1.1-noD3' is written too."""
    row = run_point(csv_name, system, pid, level, method, basis,
                    lambda: energy_of(method, atoms, charge, basis, dm_key=dm_key), coordinate=coordinate, retry_failed=RETRY_FAILED, quiet=True)
    if method == M_SKALA:
        ensure_nod3(csv_name, system, pid, level, basis, coordinate)
    return row


def g_align_retry_modes(csv_name, standin_levels=None):
    """Review row 11: when one species of a barrier was computed with the density-fitting retry (extra.df = true) and its partner without,
    recompute the partner(s) with density fitting too.  The new row gets the point id 'SP@df' (same CSV), so nothing is overwritten and the
    check (item I14) uses the aligned energy.  Returns the number of rows added."""
    rows = [r for r in read_rows(csv_name) if "@" not in r["point_id"] and not re.match(r"^p\d+$", r["point_id"]) and r["converged"] in ("1", "True")]
    groups = {}
    for r in rows:
        groups.setdefault((r["system"], r["geometry_level"], r["method"], r["basis"]), []).append(r)
    n = 0
    for (system, level, method, basis), rs in groups.items():
        if method in (M_SKALA_NOD3,) or method.startswith("v4-"):
            continue
        fitted = [r for r in rs if extra_of(r).get("df")]
        if not fitted or len(fitted) == len(rs):
            continue
        need = [r for r in rs if not extra_of(r).get("df")]
        for r in need:
            g = geometry_of(system, r["point_id"], level)
            if g is None:
                continue
            with item_guard("phase11", f"{system}:{r['point_id']}:{method}:align"):
                new = run_point(csv_name, system, r["point_id"] + "@df", level, method, basis,
                                lambda g=g, system=system, r=r, method=method, basis=basis: energy_of(method, g["atoms"], species_charge(system, r["point_id"]), basis, force_df=True),
                                retry_failed=RETRY_FAILED, quiet=True)
                if new is not None:
                    n += 1
                    if method == M_SKALA:
                        ensure_nod3(csv_name, system, r["point_id"] + "@df", level, basis)
    if n:
        print(f"  retry-mode alignment: {n} partner species recomputed with density fitting in {csv_name}")
    return n


def g1_system(system):
    """G1, stationary points of one system, all functionals, at EVERY stored geometry level (review row 4): our own B3LYP geometries first
    (the headline protocol), then the downloaded database geometries wherever the CPU notebook found them (I11 cross-check).
    The headline level is complete for the whole system before the cross-check level is started."""
    basis = basis_for(system, "main")
    present = {sp: g_levels(system, sp) for sp in SPECIES_KIND[system]}
    for sp, lvs in present.items():
        if not lvs:
            log_error("phase11", f"{system}:{sp}", RuntimeError(f"skipped: missing input geometry for {sp} (CPU notebook C1 / G-1)"), note="skipped")
    levels = [lv for lv in LEVEL_PREFERENCE if any(lv in lvs for lvs in present.values())]
    for lv in levels:                                      # own_b3lyp (and own_scanpeak) first, database last
        for sp in SPECIES_KIND[system]:
            if lv not in present[sp]:
                continue
            g = geometry_of(system, sp, lv)
            for m in ladder_methods("stationary"):
                with item_guard("phase11", f"{system}:{sp}:{m}"):
                    g_energy_row(PHASE_CSV["phase11"], system, sp, lv, m, basis, g["atoms"], species_charge(system, sp))


def g1_path(system, method):
    """G1, one functional along the relaxed path of one system (densities of the previous point are used as the start)."""
    pts = g_load_path(system)
    if not pts:
        log_error("phase11p", f"{system}:path", RuntimeError(f"skipped: missing input geoms/{system}_path.json (CPU notebook C2 / G-1)"), note="skipped")
        return
    basis = basis_for(system, "main")
    reset_dm(system)
    for q in pts:
        with item_guard("phase11p", f"{system}:p{q['index']}:{method}"):
            g_energy_row(PHASE_CSV["phase11"], system, f"p{q['index']:02d}", LV_SCAN, method, basis, q["atoms"], species_charge(system, "TS"),
                         coordinate=f"{q['coordinate']:.4f}", dm_key=(system, method, basis))


def g1_dft_stationary():
    """G1, stage 1: stationary points and fragments at every stored geometry level, then the retry-mode alignment."""
    for system in [s for s in SYSTEMS_TO_RUN]:
        with item_guard("phase11", system):
            print(f"  {system}: stationary points")
            g1_system(system)
    with item_guard("phase11", "retry-mode alignment"):
        g_align_retry_modes(PHASE_CSV["phase11"])


def g1_dft_paths():
    """G1, stage 2 (its own phase 'phase11p' so that no phase is longer than 85 minutes): path points, one functional at a time."""
    systems = [s for s in SYSTEMS_TO_RUN]
    skipped = [m for m in ladder_methods("stationary") if m not in ladder_methods("path")]
    if skipped:
        print(f"  note: {', '.join(skipped)} run on the CPU, so they are computed at stationary points only (not along paths)")
    for method in ladder_methods("path"):                    # one functional at a time (cheapest to cut last)
        for system in systems:
            with item_guard("phase11p", f"{system}:path:{method}"):
                print(f"  {system}: path, {method}")
                g1_path(system, method)


def g1_dft():
    g1_dft_stationary()
    g1_dft_paths()


def g2_basis():
    """G2: the same stationary points in the small and the large basis, method by method (PBE and r2SCAN are cut first if time runs out)."""
    csv2 = PHASE_CSV["phase12"]
    species = [(s, sp) for s in SYSTEMS_TO_RUN if s != "ch4" for sp in SPECIES_KIND[s]]
    for m in ladder_methods("stationary"):
        for size in ("large", "small"):
            print(f"  {m}, {size} basis")
            for system, sp in species:
                with item_guard("phase12", f"{system}:{sp}:{m}:{size}"):
                    lv, g = best_geom(system, sp)
                    if g is None:
                        log_error("phase12", f"{system}:{sp}", RuntimeError("skipped: missing input geometry"), note="skipped")
                        continue
                    b, ch = basis_for(system, size), species_charge(system, sp)
                    if size == "large":
                        nbf = make_mol(g["atoms"], ch, b).nao
                        if nbf > GPU_QZ_MAX_NBF:
                            if not have(csv2, system, sp, lv, m, b):
                                log_error("phase12", f"{system}:{sp}", RuntimeError(f"skipped: {b} has {nbf} basis functions, above the limit {GPU_QZ_MAX_NBF} for this GPU"), note="skipped")
                            continue
                    g_energy_row(csv2, system, sp, lv, m, b, g["atoms"], ch)
    with item_guard("phase12", "retry-mode alignment"):
        g_align_retry_modes(csv2)


def g3_scout():
    """G3 (review row 9): ALL five functionals (Skala with and without D3, PBE, r2SCAN, B3LYP-D3(BJ), wB97M-V) at the xTB minimum, peak and product of
    every reaction, so that the geometry penalty can be shown not to be Skala-specific.  xTB's own barrier (from the C4 energies) is recorded too."""
    csv3 = PHASE_CSV["phase13"]
    reactions = [s for s in SYSTEMS_TO_RUN if s in REACTIONS]
    order = sorted(reactions, key=lambda s: (SYSTEMS[s]["tier"] != "C", SYSTEM_ORDER.index(s)))   # Tier C first (it has no CCSD(T))
    for m in ladder_methods("stationary"):
        for system in order:
            basis = basis_for(system, "main")
            for sp in g_scout_species(system):
                if not sp:
                    continue
                with item_guard("phase13", f"{system}:{sp}:{m}"):
                    g = load_geom(geom_name(system, sp, LV_XTB))
                    if g is None:
                        log_error("phase13", f"{system}:{sp}", RuntimeError("skipped: missing input xTB geometry (CPU notebook C4 / G-1)"), note="skipped")
                        continue
                    if g.get("no_peak"):
                        continue                                # recorded by the CPU notebook as a result: nothing to compute
                    g_energy_row(csv3, system, sp, LV_XTB, m, basis, g["atoms"], species_charge(system, sp))
    with item_guard("phase13", "retry-mode alignment"):
        g_align_retry_modes(csv3)
    with item_guard("phase13", "xTB own barriers"):
        out = {}
        e = {(r["system"], r["point_id"]): float(r["energy_hartree"]) for r in read_rows(PHASE_CSV["phase04"])
             if r["method"] == M_XTB and r["geometry_level"] == LV_XTB and r["converged"] in ("1", "True") and r["energy_hartree"] not in ("", "nan")}
        for system in reactions:
            rmin = g_scout_species(system)[0]
            out[system] = None if (system, "TS") not in e or (system, rmin) not in e else round((e[(system, "TS")] - e[(system, rmin)]) * HARTREE_KCAL, 3)
        write_json("phase13_xtb_barriers.json", dict(unit="kcal/mol", barrier_TS_minus_reactant_or_complex=out,
                                                      note="GFN2-xTB energies at GFN2-xTB geometries (from CPU notebook C4); null = no xTB peak or missing"))
        print("  xTB's own barriers (kcal/mol):", out)


# ------------------------------------------------------------------ G4: Skala-relaxed minima and transition states
def g4_minima(system):
    S = SYSTEMS[system]
    out = []
    for sp in S["reactants"] + S["products"] + [x for x in (S["rc"], S["pc"]) if x]:
        if sp not in out and SPECIES_KIND[system].get(sp) == "min":
            out.append(sp)
    return out


def g4_system(system, variant):
    S, spc = SYSTEMS[system], SCANS.get(system)
    basis, csv4 = basis_for(system, "main"), PHASE_CSV["phase14"]

    def relax(atoms, ch, constraints=None):
        eng = EG(M_SKALA, ch, basis, variant)
        return opt_custom(atoms, ch, eng.eg, constraints=constraints, convergence_set=G4_CONV, maxiter=G4_MAXITER)

    def e_std(atoms, ch):                                     # the energy with the same engine as G1 (not the gradient engine)
        return energy_of(M_SKALA, atoms, ch, basis)
    # ---- minima
    for sp in g4_minima(system):
        with item_guard("phase14", f"{system}:{sp}"):
            lv, g = best_geom(system, sp)
            if g is None:
                log_error("phase14", f"{system}:{sp}", RuntimeError("skipped: missing input geometry"), note="skipped")
                continue
            ch, nm = species_charge(system, sp), geom_name(system, sp, LV_SKALA)
            gs = load_geom(nm)
            if gs is None:
                t0 = time.time()
                r = relax(g["atoms"], ch)
                save_geom(nm, r["atoms"], ch, species=sp, kind="min", level=LV_SKALA, opt_converged=bool(r["converged"]), n_gradients=r["n_steps"], seconds=round(time.time() - t0, 1), error=r["error"])
                gs = load_geom(nm)
            run_point(csv4, system, sp, LV_SKALA, M_SKALA, basis,
                      lambda gs=gs, ch=ch: (lambda e: dict(energy=e["energy"], converged=bool(e["converged"] and gs.get("opt_converged", True)),
                                                           extra=dict(e["extra"], n_gradients=gs.get("n_gradients"), opt_error=gs.get("error", ""), start_level=lv)))(e_std(gs["atoms"], ch)),
                      retry_failed=RETRY_FAILED)
            ensure_nod3(csv4, system, sp, LV_SKALA, basis)
    # ---- transition state: constrained optimisations around the TS coordinate, then a parabola through the top three points
    if not S["ts"] or spc is None:
        return
    with item_guard("phase14", f"{system}:TS"):
        if have(csv4, system, "TS", LV_SKALA, M_SKALA, basis) and not RETRY_FAILED:
            return
        lv, gTS = best_geom(system, "TS")
        if gTS is None:
            log_error("phase14", f"{system}:TS", RuntimeError("skipped: missing input TS geometry"), note="skipped")
            return
        ch = species_charge(system, "TS")
        kind, idx = spc["kind"], spc["idx"]
        c_ts, delta = coord_value(kind, idx, gTS["atoms"]), G4_DELTA[kind]
        sp_file = path("geoms", f"{system}_skala_ts_scan.json")
        store = json.load(open(sp_file)) if os.path.exists(sp_file) else {}

        def point(c, start):
            key = f"{c:.5f}"
            if key not in store:
                t0 = time.time()
                r = relax(start, ch, constraints=[(kind, idx, float(c))])
                e = e_std(r["atoms"], ch)
                store[key] = dict(coordinate=float(c), atoms=[[s, list(x)] for s, x in r["atoms"]], energy=e["energy"], converged=bool(r["converged"] and e["converged"]),
                                  d3=e["extra"].get("d3_hartree"), n_gradients=r["n_steps"], seconds=round(time.time() - t0, 1), error=r["error"])
                json.dump(store, open(sp_file + ".tmp", "w"), indent=1)
                os.replace(sp_file + ".tmp", sp_file)
                print(f"   {system} TS window: coordinate {c:.4f}  E = {e['energy']:.8f}  ({r['n_steps']} gradients)")
            return store[key]
        for c in (c_ts - delta, c_ts, c_ts + delta):
            point(c, gTS["atoms"])
        for shift in range(G4_MAX_SHIFT + 1):           # move the window if the top is not in the middle
            keys = sorted(store, key=float)
            k = int(np.nanargmax([store[k_]["energy"] for k_ in keys]))
            if (0 < k < len(keys) - 1) or shift == G4_MAX_SHIFT:
                break
            new_c = float(keys[0]) - delta if k == 0 else float(keys[-1]) + delta
            point(new_c, [(s, tuple(x)) for s, x in store[keys[k]]["atoms"]])
        keys = sorted(store, key=float)
        pts = [dict(coordinate=store[k_]["coordinate"], energy=store[k_]["energy"]) for k_ in keys]
        k = int(np.nanargmax([p["energy"] for p in pts]))
        interior = 0 < k < len(pts) - 1
        if interior:
            xv, ev = parabola_peak(pts, k)
        else:
            xv, ev = pts[k]["coordinate"], pts[k]["energy"]
        top = store[keys[k]]
        save_geom(geom_name(system, "TS", LV_SKALA), [(s, tuple(x)) for s, x in top["atoms"]], ch, species="TS", kind="ts", level=LV_SKALA,
                  vertex_coordinate=xv, interior=bool(interior), method="constrained optimisation + parabola vertex (not a saddle-point search)")
        conv_ts = bool(interior and all(store[x]["converged"] for x in keys[max(0, k - 1):k + 2]))
        d3s = [store[k_].get("d3") for k_ in keys]
        nod3 = None
        if all(v is not None for v in d3s):            # the same window without the D3 term (same SCFs, same geometries)
            pts_n = [dict(coordinate=store[k_]["coordinate"], energy=store[k_]["energy"] - store[k_]["d3"]) for k_ in keys]
            ev_n = parabola_peak(pts_n, k)[1] if interior else pts_n[k]["energy"]
            nod3 = (ev_n, [round(p["energy"], 8) for p in pts_n])
        common = dict(method="constrained optimisation + three-point parabola VERTEX (not a saddle-point search)", interior=bool(interior), window=[round(float(x), 4) for x in keys],
                      vertex_coordinate=xv, reference_ts_coordinate=c_ts, standin=STANDIN)
        if nod3 is not None and not have(csv4, system, "TS", LV_SKALA, M_SKALA_NOD3, basis):      # written first: a cut between the two rows is repaired on rerun
            append_row(csv4, make_row(system, "TS", LV_SKALA, M_SKALA_NOD3, basis, nod3[0], conv_ts, coordinate=f"{xv:.4f}",
                                      extra=dict(common, window_energies=nod3[1], derived_from=M_SKALA, same_scf=True, with_dftd3=False)))
        row = make_row(system, "TS", LV_SKALA, M_SKALA, basis, ev, conv_ts, coordinate=f"{xv:.4f}",
                       extra=dict(common, window_energies=[round(p["energy"], 8) for p in pts], window_d3=[None if v is None else round(v, 8) for v in d3s], with_dftd3=True))
        append_row(csv4, row)
        print(f"   {system} TS: parabola vertex at {xv:.4f} (reference TS {c_ts:.4f}), E = {ev:.8f}, interior maximum: {interior}")


def g4_report(system):
    """Print Skala//Skala against Skala//reference barrier for one system (from the CSV files)."""
    S, basis = SYSTEMS[system], basis_for(system, "main")

    def E(csv_, sp, lv=None):
        for r in read_rows(csv_):
            if r["system"] == system and r["point_id"] == sp and r["method"] == M_SKALA and r["basis"] == basis and r["converged"] in ("1", "True") and r["energy_hartree"] not in ("", "nan") and (lv is None or r["geometry_level"] == lv):
                if lv is not None or r["geometry_level"] in LEVEL_PREFERENCE:
                    return float(r["energy_hartree"])
        return None
    out = {}
    for tag, csv_, lv in (("Skala//Skala", PHASE_CSV["phase14"], LV_SKALA), ("Skala//reference", PHASE_CSV["phase11"], None)):
        ts = E(csv_, "TS", lv)
        parts = [E(csv_, sp, lv) for sp in S["reactants"]]
        out[tag] = None if ts is None or any(p is None for p in parts) else (ts - sum(parts)) * HARTREE_KCAL
    if out["Skala//Skala"] is not None and out["Skala//reference"] is not None:
        print(f"   {system}: barrier Skala//Skala {out['Skala//Skala']:.2f}, Skala//reference {out['Skala//reference']:.2f}, difference {out['Skala//Skala'] - out['Skala//reference']:+.2f} kcal/mol")
    return out


def g4_relaxed(systems=None):
    variant = GSTATE.get("grad_variant", "plain")
    todo = [s for s in (systems or G4_SYSTEMS) if s in SYSTEMS_TO_RUN]
    summary = {}
    for system in todo:
        with item_guard("phase14", system):
            print(f"  {system}: Skala-relaxed minima and TS (gradient settings: {variant})")
            g4_system(system, variant)
            summary[system] = g4_report(system)
    old = read_json("phase14_summary.json", {}) or {}
    old.update(summary)
    write_json("phase14_summary.json", old)


# ------------------------------------------------------------------ optional phase: GPU relaxed paths for further systems, then all functionals along them
def gopt_paths():
    csv_o = CSV_OPT
    for system in [s for s in OPT_PATH_SYSTEMS if s in SCANS and s in SYSTEMS_TO_RUN]:
        with item_guard("phase16", f"{system}:scan"):
            sp = SCANS[system]
            lv, g0 = best_geom(system, sp["start"])
            if g0 is None:
                log_error("phase16", system, RuntimeError(f"skipped: missing start geometry {sp['start']}"), note="skipped")
                continue
            ch, basis = species_charge(system, sp["start"]), basis_for(system, "main")
            eng = EG(M_B3LYP_NODISP, ch, basis)
            print(f"  {system}: relaxed scan with {M_B3LYP_NODISP} on the GPU ({sp['label']})")
            _walk(system, g0["atoms"], ch,
                  lambda at, c, ch=ch: opt_custom(at, ch, eng.eg, constraints=[c], convergence_set=G4_CONV, maxiter=G4_MAXITER),
                  lambda at, eng=eng: eng.energy(at), g_load_path(system, "pathgpu"),
                  lambda pts, system=system, sp=sp: g_save_path(system, pts, sp["label"] + " (GPU B3LYP, no dispersion)", suffix="pathgpu", level=LV_SCANGPU))
    for m in ladder_methods("path"):
        for system in [s for s in OPT_PATH_SYSTEMS if s in SCANS and s in SYSTEMS_TO_RUN]:
            with item_guard("phase16", f"{system}:path:{m}"):
                pts = g_load_path(system, "pathgpu")
                if not pts:
                    continue
                basis = basis_for(system, "main")
                reset_dm(system)
                ch = species_charge(system, SCANS[system]["start"])
                print(f"  {system}: functional {m} along the GPU path")
                for q in pts:
                    with item_guard("phase16", f"{system}:p{q['index']}:{m}"):
                        g_energy_row(csv_o, system, f"p{q['index']:02d}", LV_SCANGPU, m, basis, q["atoms"], ch, coordinate=f"{q['coordinate']:.4f}", dm_key=(system, m, basis))


# ------------------------------------------------------------------ summary helpers
def csv_stats(csv_name):
    rows = read_rows(csv_name)
    if not rows:
        return None
    bad = [r for r in rows if r["converged"] in ("0", "False") or r["energy_hartree"] in ("", "nan")]
    return dict(rows=len(rows), failed=len(bad))


def g_summary_lines():
    L = ["", "## GPU notebook details", "", f"GPU: {GPU_NAME or 'none'} ({GPU_CLASS}, {GPU_MEM_GB:.0f} GB)"]
    if STANDIN:
        L += ["", "**STAND-IN DRY RUN: every number in this folder comes from cheap CPU functionals in tiny bases. None of it is Skala or a result.**"]
    st = read_json("phase10_timing.json", {}) or {}
    L.append(f"Routes per functional: {st.get('routes')}; anchor functional: {st.get('anchor')}")
    L.append(f"G0 gate (stored v4 energies reproduced, Skala WITHOUT D3 as in v4): {st.get('gate_ok')}; Skala usable: {st.get('skala_ok')}; Skala gradient check passed: {st.get('grad_ok')} (settings: {st.get('grad_variant')})")
    w = st.get("skala_weights") or {}
    L.append(f"Skala weights: {w.get('skala_name')} -> file {w.get('file')}, SHA-256 {w.get('sha256')}, Hugging Face revision {w.get('huggingface_revision')}, "
             f"matches the package table: {w.get('hash_matches_package_table')}, D3 parameter set stored in the weights: {w.get('expected_d3_settings')}")
    L.append("Skala-1.1 = with its D3 term (package default, headline); Skala-1.1-noD3 = same SCF minus the stored D3 energy; D3 off only in the v4 gate.")
    if st.get("anchor") and st.get("anchor") != M_WB97MV:
        L.append(f"NOTE: the anchor functional is {st['anchor']}, not wB97M-V: check item R2 of the CPU notebook looks for wB97M-V rows and will say NOT RUN.")
    L += ["", "| table | rows | failed or unconverged |", "|---|---|---|"]
    for name in (CSV_REPRO, PHASE_CSV["phase11"], PHASE_CSV["phase12"], PHASE_CSV["phase13"], PHASE_CSV["phase14"], CSV_OPT):
        s = csv_stats(name)
        L.append(f"| {name} | {s['rows'] if s else 'not written'} | {s['failed'] if s else ''} |")
    cfg = read_json("config.json", {}) or {}
    if cfg.get("smoke") and not STANDIN:
        L.append("WARNING: config.json says the CPU notebook ran as a smoke test (small bases). Its numbers are not the paper numbers.")
    if os.path.exists(path("phase09_regen.json")):
        r = read_json("phase09_regen.json", {})
        L += ["", "G-1 (rebuilt CPU inputs): " + "; ".join(r.get("log", []))[:400]]
    return L
