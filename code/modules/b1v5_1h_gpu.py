# ==============================================================================
# b1v5_1h_gpu.py  --  the phases of the 1-hour GPU notebook (PLAN_1hour.md section 2, "GPU notebook").
# Pasted after b1v5_gpu.py (GPU helpers, energy_of with the Skala D3 handling, G0 helpers), b1v5_gpu_figs.py, b1v5_1h.py (queue) and b1v5_1h_check.py.
# Calculations run through qrun (admission test before EVERY calculation).  Nothing of b1v5_gpu.py is edited; this file calls it.
#
#   install, weights   (cell 1 and weights_setup): pinned skala-cuda12x 2026.9, gpu4pyscf-cuda12x 1.8.1, pyscf; weights from Hugging Face, then the Drive cache
#   G0   gate: stored v4 HCN point 0 with Skala (D3 OFF as in v4), PBE, B3LYP; 1e-5 hartree.  No gradient test (G4 is not in the core).
#   G1a  20 species of the six reactions, five functionals, def2-TZVP(D); Skala with D3, "without D3" is the same SCF minus the stored D3 energy
#   G1b  path points (HCN 11, NH3 9, H2CO 11, CH4 9 = 40): Skala, wB97M-V, B3LYP-D3(BJ)
#   G3   xTB minimum and peak of the six reactions: Skala, wB97M-V, B3LYP-D3(BJ)
#   G2   def2-QZVPP(D) at the 9 stationary species of HCN, NH3, H2CO: Skala and wB97M-V
#   tail (1) F- + CH3Cl and OH- + CH3F stationary points, five functionals (2) PBE and r2SCAN on the paths (3) QZVPP(D) for CH3NC and the SN2 systems, B3LYP added
#        (4) PBE and r2SCAN at the xTB geometries (5) Skala-relaxed HCN (needs the Skala gradient test), CH3NC, F- + CH3F
#   G5   figures, tables, findings_outcome.csv, the automated check: on whatever data exists
#   EXTENSION = 'GPU1' runs tail items 1 to 3, 'GPU2' items 4 and 5, 'GEOM' rebuilds the geometry pack.  GPU class decides how much of the tail is tried by default.
# ==============================================================================
import os, json, time
import numpy as np

GPU_CAPS = dict(install=15, phase10=4, phase11=24, phase11p=12, phase13=12, phase12=8, phase15=5, check=5)       # minutes (plan section 2)
GPU_SPEED = dict(T4=1.0, L4=0.75, A100=0.35, OTHER=1.0, NONE=1.0)           # time per energy relative to the T4 (L4 and A100: SCALED GUESSES, not measured)
GPU_TAIL_BY_CLASS = dict(T4=[], OTHER=[], L4=[1, 2], A100=[1, 2, 3, 4, 5], NONE=[])                              # plan: T4 none expected; L4 (1), perhaps (2); A100 (1) to (4), perhaps (5)
GPU_PRIOR_T4 = {"gpu": 8.0, "gpu:tz": 7.0, "gpu:tz:hcn": 5.0, "gpu:tz:nh3": 5.0, "gpu:tz:h2co": 6.0, "gpu:tz:ch4": 5.0, "gpu:tz:ch3nc": 9.0,
                "gpu:tz:sn2_ff": 12.0, "gpu:tz:sn2_clcl": 16.0, "gpu:tz:sn2_fcl": 14.0, "gpu:tz:sn2_ohf": 14.0,
                "gpu:qz": 27.0, "gpu:qz:ch3nc": 60.0, "gpu:qz:sn2_ff": 90.0, "gpu:qz:sn2_clcl": 120.0, "gpu:qz:sn2_fcl": 110.0, "gpu:qz:sn2_ohf": 110.0,
                "gpu:v4": 15.0, "gpu:g4": 900.0, "gpu:align": 40.0, "gpu:probe": 20.0}
Q["prior"].update(GPU_PRIOR_T4)
CSV11, CSV12, CSV13 = PHASE_CSV["phase11"], PHASE_CSV["phase12"], PHASE_CSV["phase13"]
PATH_METHODS = [M_SKALA, "ANCHOR", M_B3LYP]               # G1b and G3 (core)
BASIS_METHODS = [M_SKALA, "ANCHOR"]                       # G2 (core)
PATH_SYSTEMS = ["hcn", "nh3", "h2co", "ch4"]
BASIS_SYSTEMS = ["hcn", "nh3", "h2co"]
SKALA_PREFETCH_FAILED = False


def anchor_name():
    return anchor_method()


def resolve_methods(names):
    out = []
    for m in names:
        m = anchor_name() if m == "ANCHOR" else m
        if route_of(m) != "none" and m not in out:
            out.append(m)
    return out


def path_ok_methods(names):
    """Methods that may run along paths (a CPU-routed functional is for stationary points only)."""
    return [m for m in resolve_methods(names) if route_of(m) != "cpu" or STANDIN]


# ------------------------------------------------------------------ install, weights
def weights_setup(timeout_min):
    """Restore the Skala weights from the Drive cache, else download them once (HuggingFace) in a SUBPROCESS with a wall-clock cap, then cache them on Drive.
    A slow or failed download costs at most `timeout_min`; the queue then simply reaches less.  Never raises."""
    import subprocess
    global SKALA_PREFETCH_FAILED
    out = dict(restore=restore_weights_cache())
    print("  Skala weights:", out["restore"])
    if STANDIN or NO_GPU:
        out["prefetch"] = "skipped (stand-in or no GPU)"
        return out
    code = ("from skala.functional import resolve_functional_artifact as r; a = r('%s'); print(a.path)" % SKALA_NAME)
    t0 = time.time()
    try:
        r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=max(30.0, 60.0 * timeout_min))
        out["prefetch"] = ("ok: " + r.stdout.strip()[-120:]) if r.returncode == 0 else ("failed: " + (r.stderr.strip()[-200:] or "no message"))
        SKALA_PREFETCH_FAILED = r.returncode != 0
    except Exception as ex:
        out["prefetch"] = "stopped at the cap: " + repr(ex)[:100]
        SKALA_PREFETCH_FAILED = True
    out["seconds"] = round(time.time() - t0, 1)
    if not SKALA_PREFETCH_FAILED:
        out["save"] = save_weights_cache()
    print(f"  Skala weights prefetch: {out['prefetch']} ({out['seconds']} s); cache: {out.get('save')}")
    return out


# ------------------------------------------------------------------ queue wrappers
def q_energy_row(csv_name, system, pid, level, method, basis, atoms, charge, kind, coordinate="", dm_key=None):
    """One energy through the queue; Skala also gets its 'without D3' row from the same SCF (free)."""
    row = qrun(csv_name, system, pid, level, method, basis, lambda: energy_of(method, atoms, charge, basis, dm_key=dm_key), kind,
               coordinate=coordinate, retry_failed=RETRY_FAILED, label=f"{system}:{pid}:{method}:{basis}")
    if method == M_SKALA:
        ensure_nod3(csv_name, system, pid, level, basis, coordinate)
    return row


def own_geom(system, sp):
    lv, g = best_geom(system, sp)
    return lv, g


def blocked_1h(name):
    """Reason (text) why a GPU phase must not start, or None."""
    if name == "phase10":
        return None
    if NO_GPU and not STANDIN:
        return "no GPU detected: only the G0 diagnostics run (Runtime > Change runtime type > GPU, then Run all again)"
    load_gstate()
    if GSTATE.get("gate_ok") is False and STOP_ON_G0_FAIL:
        return "G0 reproduction of the stored v4 energy failed (difference over 1e-4 hartree): old and new numbers must not be mixed. Fix the install, delete DONE_phase10 and run again (or set STOP_ON_G0_FAIL = False)"
    return None


def run_gphase_1h(name, title, fn, cap_min, final=False, force=None):
    why = blocked_1h(name)
    if why:
        skip_phase(name, title, why)
        return
    with qphase(name, title, cap_min, final=final, force=(FORCE_RERUN if force is None else force)) as ph:
        if not ph.skip:
            fn()


# ------------------------------------------------------------------ G0
def g0_1h():
    """G0 of the 1-hour design: weights cache, functional probe, Skala check, weights hash, the stored v4 HCN point 0 with Skala (D3 OFF), PBE, B3LYP."""
    info = dict(started=now(), hardware=hw(), versions=ver(), gpu_name=GPU_NAME, gpu_class=GPU_CLASS, gpu_mem_gb=round(GPU_MEM_GB, 1), standin=STANDIN, design="1h")
    if NO_GPU and not STANDIN:
        info["note"] = "no GPU: diagnostics only"
        save_gstate(**info)
        print("  ** NO GPU: G0 stops here after the diagnostics. **")
        return info
    if SKALA_PREFETCH_FAILED:
        globals()["CALC_MAX_SEC"] = 120                      # the weights are not cached: do not wait 20 minutes for a download inside the Skala check
        print("  weights prefetch failed earlier: the Skala check is limited to 120 s")
    ok, why = q_admit("gpu:probe")
    if ok:
        with item_guard("phase10", "probe functionals"):
            info["probe"] = g_probe_functionals()
            save_gstate(**info)
        with item_guard("phase10", "Skala check"):
            info["skala"] = g_skala_check()
            save_gstate(**info)
    else:
        q_refuse("G0 probe and Skala check", "gpu:probe", why)
    with item_guard("phase10", "weights cache save"):
        if GSTATE.get("skala_ok") and not STANDIN:
            info["weights_save"] = save_weights_cache()
    with item_guard("phase10", "weights hash"):
        info["skala_weights"] = skala_weights_info()
        w = info["skala_weights"]
        print("  Skala weights:", {k: w.get(k) for k in ("file", "sha256", "huggingface_revision", "hash_matches_package_table", "expected_d3_settings")})
        save_gstate(**info)
    th, r_ = NB4_HCN_GEOMS[0]
    at = hcn_geometry(th, r_)
    basis = "STO-3G" if STANDIN else "def2-TZVP"
    for m in ("v4-skala", "v4-pbe", "v4-b3lyp"):
        if m == "v4-skala" and not GSTATE.get("skala_ok"):
            continue
        with item_guard("phase10", f"v4 point 0 {m}"):
            reset_dm("nb4")
            qrun(CSV_REPRO, "hcn_nb4", "nb4_00", LV_NB4, m, basis, lambda m=m: energy_of(m, at, 0, basis, dm_key=("nb4", m, basis)), "gpu:v4",
                 coordinate=f"{th},{r_}", retry_failed=RETRY_FAILED, label=f"stored v4 HCN point 0, {m}")
    gate = g_gate()
    info["gate"] = gate
    print(f"  stored v4 HCN point 0: largest |difference| {gate['max_diff']:.1e} hartree ({gate['status']}); per method {gate['per_method']}")
    if STANDIN:
        GSTATE["gate_ok"] = True
        print("  stand-in: the gate is not applied (stand-in energies differ from the stored v4 energies by design)")
    elif gate["status"] == "no data":
        GSTATE["gate_ok"] = None
        print("  gate: no reproduced energy (all attempts failed or were refused); later phases run, but nothing was verified")
    else:
        GSTATE["gate_ok"] = gate["status"] != "FAIL"
        if gate["status"] == "WARN":
            print("  gate WARN: differences between 1e-5 and 1e-4 hartree: version drift, state it in the paper")
        if gate["status"] == "FAIL":
            print("  ** GATE FAILED: difference over 1e-4 hartree. The later phases will not start (STOP_ON_G0_FAIL). **")
    GSTATE["grad_ok"] = None                              # no gradient test in the core (tail item 5 runs it when it is reached)
    save_gstate(**info)
    return info


# ------------------------------------------------------------------ G1a, G1b, G3, G2
def g1a_stationary(systems=None, methods=None):
    """Every species of the design reactions with the ladder, own geometry, main basis."""
    meths = resolve_methods(methods or ["ANCHOR" if x == "ANCHOR" else x for x in PRIORITY])
    for system in [s for s in (systems or SYSTEMS_TO_RUN) if s != "ch4"]:
        basis = basis_for(system, "main")
        print(f"  {system}: stationary points")
        for sp in SPECIES_KIND[system]:
            lv, g = own_geom(system, sp)
            if g is None:
                log_error("phase11", f"{system}:{sp}", RuntimeError(f"skipped: missing input geometry for {sp} (geometry pack)"), note="skipped")
                continue
            for m in meths:
                with item_guard("phase11", f"{system}:{sp}:{m}"):
                    q_energy_row(CSV11, system, sp, lv, m, basis, g["atoms"], species_charge(system, sp), f"gpu:tz:{system}:{m}")


def g1a_align():
    """Review row 11: partners of a density-fitted retry are recomputed the same way (rare; admitted once)."""
    rows = read_rows(CSV11)
    if not any(extra_of(r).get("df") for r in rows):
        return
    ok, why = q_admit("gpu:align")
    if not ok:
        q_refuse("retry-mode alignment", "gpu:align", why)
        return
    g_align_retry_modes(CSV11)


def g1b_paths(systems=None, methods=None):
    """One functional at a time along the relaxed paths (HCN 11, NH3 9, H2CO 11, CH4 9 points)."""
    meths = path_ok_methods(methods or PATH_METHODS)
    for method in meths:
        for system in [s for s in PATH_SYSTEMS if s in (systems or SYSTEMS_TO_RUN)]:
            pts = g_load_path(system)
            if not pts:
                log_error("phase11p", f"{system}:path", RuntimeError(f"skipped: missing input geoms/{system}_path.json (geometry pack)"), note="skipped")
                continue
            basis, ch = basis_for(system, "main"), species_charge(system, "TS")
            reset_dm(system)
            print(f"  {system}: path, {method}")
            for q in pts:
                with item_guard("phase11p", f"{system}:p{q['index']}:{method}"):
                    q_energy_row(CSV11, system, f"p{q['index']:02d}", LV_SCAN, method, basis, q["atoms"], ch, f"gpu:tz:{system}:{method}",
                                 coordinate=f"{q['coordinate']:.4f}", dm_key=(system, method, basis))


def scout_species(system):
    S = SYSTEMS[system]
    return [S["rc"] or S["reactants"][0], "TS"]


def g3_scout_1h(methods=None, systems=None):
    """The xTB minimum and peak of every design reaction with the functionals of the list (core: Skala, wB97M-V, B3LYP-D3(BJ))."""
    meths = resolve_methods(methods or PATH_METHODS)
    for m in meths:
        for system in [s for s in (systems or SYSTEMS_TO_RUN) if s in REACTIONS]:
            basis = basis_for(system, "main")
            for sp in scout_species(system):
                with item_guard("phase13", f"{system}:{sp}:{m}"):
                    g = load_geom(geom_name(system, sp, LV_XTB))
                    if g is None:
                        log_error("phase13", f"{system}:{sp}", RuntimeError("skipped: missing input xTB geometry (geometry pack)"), note="skipped")
                        continue
                    if g.get("no_peak"):
                        continue
                    q_energy_row(CSV13, system, sp, LV_XTB, m, basis, g["atoms"], species_charge(system, sp), f"gpu:tz:{system}:{m}")
    g3_xtb_barriers()


def g3_xtb_barriers():
    with item_guard("phase13", "xTB own barriers"):
        out = {}
        e = {(r["system"], r["point_id"]): float(r["energy_hartree"]) for r in read_rows(PHASE_CSV["phase04"])
             if r["method"] == M_XTB and r["geometry_level"] == LV_XTB and r["converged"] in ("1", "True") and r["energy_hartree"] not in ("", "nan")}
        for system in [s for s in SYSTEMS_TO_RUN if s in REACTIONS]:
            rmin = scout_species(system)[0]
            out[system] = None if (system, "TS") not in e or (system, rmin) not in e else round((e[(system, "TS")] - e[(system, rmin)]) * HARTREE_KCAL, 3)
        write_json("phase13_xtb_barriers.json", dict(unit="kcal/mol", barrier_TS_minus_reactant_or_complex=out,
                                                      note="GFN2-xTB energies at GFN2-xTB geometries (from the geometry pack); null = no xTB peak or missing"))
        print("  xTB's own barriers (kcal/mol):", out)


def g2_basis_1h(systems=None, methods=None):
    """def2-QZVPP (QZVPPD for SN2) at the stationary species; system by system so that a barrier is complete (core: HCN, NH3, H2CO; Skala and wB97M-V)."""
    meths = resolve_methods(methods or BASIS_METHODS)
    for system in systems or [s for s in BASIS_SYSTEMS if s in SYSTEMS_TO_RUN]:
        b = basis_for(system, "large")
        print(f"  {system}: {b}")
        for m in meths:
            for sp in SPECIES_KIND[system]:
                with item_guard("phase12", f"{system}:{sp}:{m}"):
                    lv, g = own_geom(system, sp)
                    if g is None:
                        log_error("phase12", f"{system}:{sp}", RuntimeError("skipped: missing input geometry"), note="skipped")
                        continue
                    ch = species_charge(system, sp)
                    nbf = make_mol(g["atoms"], ch, b).nao
                    if nbf > GPU_QZ_MAX_NBF:
                        if not have(CSV12, system, sp, lv, m, b):
                            log_error("phase12", f"{system}:{sp}", RuntimeError(f"skipped: {b} has {nbf} basis functions, above the limit {GPU_QZ_MAX_NBF} for this GPU"), note="skipped")
                        continue
                    q_energy_row(CSV12, system, sp, lv, m, b, g["atoms"], ch, f"gpu:qz:{system}:{m}")


# ------------------------------------------------------------------ the tail
def ext78_present():
    return [s for s in EXT_REACTIONS if all(geometry_of(s, sp, LV_OWN) is not None or geometry_of(s, sp, LV_SCANPEAK) is not None for sp in SPECIES_KIND[s])]


def tail1_ext78():
    """(1) F- + CH3Cl and OH- + CH3F stationary points, five functionals (needs the geometries of the extended pack)."""
    ss = ext78_present()
    if not ss:
        log_error("tail1_ext78", "", RuntimeError("skipped: the geometry pack has no geometries for F- + CH3Cl and OH- + CH3F (REBUILD_GEOMS / EXT-geom builds them)"), note="skipped")
        print("  no geometries for reactions 7 and 8 in the pack: tail item 1 skipped")
        return
    g1a_stationary(systems=ss)


def tail2_paths_pbe_r2scan():
    g1b_paths(methods=[M_R2SCAN, M_PBE])


def tail3_qzvpp():
    """(3) QZVPP(D) for CH3NC and the SN2 systems (Skala, wB97M-V, B3LYP) and B3LYP for the three Tier A reactions."""
    g2_basis_1h(systems=["hcn", "nh3", "h2co"], methods=[M_B3LYP])
    g2_basis_1h(systems=[s for s in ["ch3nc", "sn2_ff", "sn2_clcl"] + ext78_present() if s in SYSTEMS_TO_RUN or s in EXT_REACTIONS], methods=[M_SKALA, "ANCHOR", M_B3LYP])


def tail4_scout_more():
    g3_scout_1h(methods=[M_R2SCAN, M_PBE])


def g_gradient_test_lazy():
    """The Skala gradient test of the full G0, run only when the tail reaches the Skala-relaxed geometries."""
    if GSTATE.get("grad_ok") is not None:
        return GSTATE["grad_ok"]
    if not GSTATE.get("skala_ok"):
        return False
    res = g_skala_gradient_check()
    save_gstate(gradient=res)
    return bool(GSTATE.get("grad_ok"))


def tail5_skala_relaxed(systems=("hcn",)):
    """(5) Skala-relaxed minima and transition state (parabola vertex) for HCN, and CH3NC and F- + CH3F in EXT-GPU-2; geomeTRIC is installed on demand."""
    if not _have_mod("geometric"):
        pip_install_capped(["geometric==1.1.1"], 4)
    if not g_gradient_test_lazy():
        log_error("tail5_skala_relaxed", "", RuntimeError("skipped: the Skala gradient test did not pass (or Skala is unusable)"), note="skipped")
        print("  Skala gradient test not passed: tail item 5 skipped")
        return
    variant = GSTATE.get("grad_variant", "plain")
    for system in systems:
        ok, why = q_admit("gpu:g4")
        if not ok:
            q_refuse(f"{system}: Skala-relaxed minima and TS", "gpu:g4", why)
            continue
        with item_guard("tail5_skala_relaxed", system):
            t0 = time.time()
            g4_system(system, variant)
            q_record("gpu:g4", time.time() - t0)
            g4_report(system)


def _have_mod(m):
    import importlib.util
    try:
        return importlib.util.find_spec(m) is not None
    except Exception:
        return False


TAIL_ITEMS = {1: ("tail1_ext78", "tail (1) F- + CH3Cl and OH- + CH3F, five functionals", tail1_ext78),
              2: ("tail2_paths", "tail (2) PBE and r2SCAN along the paths", tail2_paths_pbe_r2scan),
              3: ("tail3_qzvpp", "tail (3) QZVPP(D) for CH3NC and SN2, B3LYP added", tail3_qzvpp),
              4: ("tail4_scout", "tail (4) PBE and r2SCAN at the xTB geometries", tail4_scout_more),
              5: ("tail5_skala_relaxed", "tail (5) Skala-relaxed HCN", tail5_skala_relaxed)}


def tail_items_for(extension, gpu_class):
    """Which tail items to try: the extension session names them; otherwise the GPU class decides (T4: none, as the plan expects)."""
    if extension == "GPU1":
        return [1, 2, 3]
    if extension == "GPU2":
        return [4, 5]
    return list(GPU_TAIL_BY_CLASS.get(gpu_class, []))


def run_tail(items):
    for k in items:
        name, title, fn = TAIL_ITEMS[k]
        if k == 5 and EXTENSION == "GPU2":
            fn = lambda: tail5_skala_relaxed(("hcn", "ch3nc", "sn2_ff"))
        if q_remaining_min() < 1.0:
            q_refuse(title, "gpu:tail", f"deadline: {q_now_min():.1f} min used")
            continue
        run_gphase_1h(name, title, fn, max(1.0, q_remaining_min()))


def design_with_extension_data():
    """After the queue: reactions 7 and 8 enter the check and the figures if their stationary points were computed (tail item 1)."""
    have11 = {r["system"] for r in read_rows(CSV11) if r["system"] in EXT_REACTIONS}
    return DESIGN_REACTIONS + [s for s in EXT_REACTIONS if s in have11]


def g_summary_lines_1h():
    return g_summary_lines() + [""] + queue_lines()
