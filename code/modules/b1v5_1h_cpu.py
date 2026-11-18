# ==============================================================================
# b1v5_1h_cpu.py  --  the phases of the 1-hour CPU notebook (PLAN_1hour.md section 2, "CPU notebook").
# Pasted after b1v5_1h.py (queue) and b1v5_geompack.py.  Calculations use the helpers of b1v5_cpu.py (ccsdt, dfmp2, make_ks, ...) and run through
# qrun, i.e. with the admission test of the deadline queue in front of EVERY calculation.
#
#   C0   verify the geometry pack: SHA-256 (done at installation), B3LYP gradient at the 10 unimolecular structures and the 2 SN2 transition
#        states (pass: RMS gradient below 3e-4 hartree/bohr), one stored v4 CCSD(T) energy (check I9)
#   C3a  CCSD(T), frozen core, def2-TZVP(D), T1, RHF stability at the stationary points: HCN, NH3, H2CO, CH4, CH3NC, then F- + CH3F
#   C3b  CCSD(T) along the stored paths: HCN 11, NH3 9, H2CO 11, CH4 9 points
#   C3c  DF-MP2/def2-QZVPP minus DF-MP2/def2-TZVP at the stationary points of HCN, NH3, H2CO, CH3NC (basis correction of the reference)
#   C4   CCSD(T) at the xTB minimum and peak: HCN, NH3, H2CO, CH3NC
#   C6   all-electron CCSD(T) barrier, HCN only (check I10)
#   tail (only if time remains, in this order): MP2 correction at the 31 path points; CCSD(T) + MP2 for Cl- + CH3Cl (Tier C -> B); CCSD(T)/def2-SVP on the 425-point HCN grid
#   EXT-CPU-1 (EXTENSION = 'CPU1'): the same tail list in a session of its own (finished items are skipped)
# ==============================================================================
import os, json, time
import numpy as np

CPU_CAPS = dict(install=6, phase00=6, phase03a=14, phase03b=7, phase08=12, phase04=4, phase06=2, check=3)       # minutes (plan section 2)
CPU_STATIONARY_ORDER = ["hcn", "nh3", "h2co", "ch4", "ch3nc", "sn2_ff"]                                        # C3a order of the plan
CPU_PATH_SYSTEMS = ["hcn", "nh3", "h2co", "ch4"]
CPU_MP2_SYSTEMS = ["hcn", "nh3", "h2co", "ch3nc", "sn2_ff"]                                                    # measured: the whole list takes about 1.2 minutes in the sandbox
CPU_SCOUT_REF_SYSTEMS = ["hcn", "nh3", "h2co", "ch3nc"]
MP2_MAX_NBF = 450                  # a point above this is skipped and flagged (the reference stays uncorrected there)
RETRY_FAILED = False

# Seconds measured in the 4-core sandbox at the plan's basis sets (README_1h.md section 7); the speed factor
# of the machine (Colab 2 vCPU: about 2) is applied on top of these before the first calculation of a kind.
CPU_PRIOR_SANDBOX = {                       # MEASURED in the 4-core sandbox on the real pack (README_1h.md section 7); seconds per calculation
    "grad": 5.0, "ccsdt": 6.0, "ccsdt:hcn": 5.0, "ccsdt:nh3": 2.5, "ccsdt:h2co": 6.0, "ccsdt:ch4": 2.5, "ccsdt:ch3nc": 30.0,
    "ccsdt:sn2_ff": 20.0, "ccsdt:sn2_ff:F-": 3.0, "ccsdt:sn2_ff:RC": 95.0, "ccsdt:sn2_ff:TS": 115.0,
    "ccsdt:sn2_clcl": 25.0, "ccsdt:sn2_clcl:Cl-": 3.0, "ccsdt:sn2_clcl:RC": 135.0, "ccsdt:sn2_clcl:TS": 140.0,
    "ccsdt:hcn_nb4": 5.0, "ccsdt:hcn_grid": 0.5, "ccsdt:hcn:ae": 5.5,
    "mp2": 3.0, "mp2:hcn:main": 0.8, "mp2:hcn:large": 2.0, "mp2:nh3:main": 0.3, "mp2:nh3:large": 1.7, "mp2:h2co:main": 1.1, "mp2:h2co:large": 3.3,
    "mp2:ch3nc:main": 3.0, "mp2:ch3nc:large": 4.7, "mp2:sn2_ff:main": 2.5, "mp2:sn2_ff:large": 5.5, "mp2:sn2_clcl:main": 3.0, "mp2:sn2_clcl:large": 7.5,
}
CAL_REF_SANDBOX_S = 2.7            # B3LYP energy + gradient of HCN at def2-TZVP in the sandbox (calibration of the speed factor)
Q["prior"].update(CPU_PRIOR_SANDBOX)


def cpu_calibrate():
    """Speed factor of THIS machine against the sandbox (one 2-second calculation); scales every prior.  Never raises."""
    try:
        t0 = time.time()
        mf = make_ks(start_atoms("hcn", "HCN"), 0, "def2-TZVP")
        mf.kernel(); mf.nuc_grad_method().kernel()
        sec = time.time() - t0
        Q["scale"] = float(min(8.0, max(0.7, sec / CAL_REF_SANDBOX_S)))
        print(f"  speed calibration: {sec:.1f} s for the reference calculation (sandbox {CAL_REF_SANDBOX_S} s) -> priors x {Q['scale']:.2f}")
    except Exception as ex:
        Q["scale"] = 2.0
        print("  speed calibration failed, priors x 2:", repr(ex)[:100])


def cc_row(csv_name, system, pid, level, atoms, ch, basis, kind, coordinate="", method=None, fc=True, stab=True, label=None):
    method = method or (M_CCSDT_FC if fc else M_CCSDT_AE)
    return qrun(csv_name, system, pid, level, method, basis, lambda: ccsdt(atoms, ch, basis, frozen_core=fc, stability=stab), kind,
                coordinate=coordinate, retry_failed=RETRY_FAILED, label=label or f"{system}:{pid}:CCSD(T){'' if fc else '-ae'}:{basis}")


def mp2_row(csv_name, system, pid, level, atoms, ch, size, kind, coordinate=""):
    b = basis_for(system, size)
    nbf = make_mol(atoms, ch, b).nao
    if nbf > MP2_MAX_NBF:
        log_error(Q["phase"] or "phase08", f"{system}:{pid}", RuntimeError(f"skipped: {b} has {nbf} basis functions, above MP2_MAX_NBF={MP2_MAX_NBF}"), note="skipped")
        return None
    return qrun(csv_name, system, pid, level, M_DFMP2_NAME, b, lambda: dfmp2(atoms, ch, b), kind, coordinate=coordinate, retry_failed=RETRY_FAILED,
                label=f"{system}:{pid}:DF-MP2:{b}")


M_DFMP2_NAME = M_MP2DF          # noqa  (b1v5_common)


def stationary_geom(system, sp):
    lv, g = best_geom(system, sp)
    return lv, g


def load_path_pts(system):
    p = path("geoms", f"{system}_path.json")
    if not os.path.exists(p):
        return []
    return [dict(q, atoms=[(s, tuple(c)) for s, c in q["atoms"]]) for q in json.load(open(p))["points"]]


def update_tiers():
    """A system whose CCSD(T) stationary rows are not all there (the queue refused them) becomes Tier C: the plan says so for F- + CH3F."""
    tiers = {k: v["tier"] for k, v in SYSTEMS.items()}
    rows = read_rows(PHASE_CSV["phase03"])
    for s in DESIGN["reactions"] + ["ch4"]:
        if tiers[s] not in ("A", "B"):
            continue
        have_ = {r["point_id"] for r in rows if r["system"] == s and r["method"] == M_CCSDT_FC and r["geometry_level"] in (LV_OWN, LV_SCANPEAK) and r["converged"] in ("1", "True")}
        if not set(SPECIES_KIND[s]) <= have_:
            if tiers[s] == "B" or s in ("hcn", "nh3", "h2co"):
                print(f"  {s}: CCSD(T) at the stationary points is incomplete ({sorted(set(SPECIES_KIND[s]) - have_)} missing): tier {tiers[s]} -> C")
                tiers[s] = "C"
    write_json("tiers.json", tiers)
    return tiers


# ------------------------------------------------------------------ C0
def c0_verify(pack_info):
    """The gradient test and one stored v4 CCSD(T) energy.  The SHA-256 test was done when the pack was installed (pack_setup)."""
    print(f"  pack: sha256 {pack_info.get('sha256')}  verified: {pack_info.get('ok')}  structures/paths checked: {pack_info.get('n_verified')}  problems: {pack_info.get('problems')}")
    cpu_calibrate()
    pack = pack_info.get("_pack")
    if pack is not None:
        gradient_test(pack, systems=SYSTEMS_TO_RUN, run=qrun)
    g = grad_summary()
    print(f"  gradient test: {g['passed']} of {g['n']} structures pass (RMS gradient < {GRAD_RMS_TOL:g} hartree/bohr); worst {g['worst_rms']:.1e} at {g['worst_where']}")
    th, r_ = NB4_HCN_GEOMS[0]
    qrun("phase00_repro.csv", "hcn_nb4", "nb4_00", LV_NB4, "v4-ccsd_t", "def2-TZVP", lambda: ccsdt(hcn_geometry(th, r_), 0, "def2-TZVP", frozen_core=False, stability=False),
         "ccsdt:hcn_nb4", coordinate=f"{th},{r_}", label="stored v4 CCSD(T) point 0 (I9)")
    tiers = {k: v["tier"] for k, v in SYSTEMS.items()}
    write_json("tiers.json", tiers)
    write_json("phase00_verify_summary.json", dict(gradient=g, pack_sha256=pack_info.get("sha256"), pack_ok=pack_info.get("ok"), speed_factor=Q["scale"], hardware=hw(), versions=ver(), at=now()))
    return g


# ------------------------------------------------------------------ C3a, C3b, C3c, C4, C6
def c3a_stationary(systems=None):
    csv_ = PHASE_CSV["phase03"]
    for system in [s for s in CPU_STATIONARY_ORDER if s in (systems or SYSTEMS_TO_RUN)]:
        print(f"  {system}: stationary points")
        for sp in SPECIES_KIND[system]:
            lv, g = stationary_geom(system, sp)
            if g is None:
                log_error("phase03a", f"{system}:{sp}", RuntimeError("skipped: missing input geometry (geometry pack not installed?)"), note="skipped")
                continue
            with item_guard("phase03a", f"{system}:{sp}"):
                cc_row(csv_, system, sp, lv, g["atoms"], species_charge(system, sp), basis_for(system, "main"), f"ccsdt:{system}:{sp}")
    return update_tiers()


def c3b_paths(systems=None):
    csv_ = PHASE_CSV["phase03"]
    for system in [s for s in CPU_PATH_SYSTEMS if s in (systems or SYSTEMS_TO_RUN)]:
        pts = load_path_pts(system)
        if not pts:
            log_error("phase03b", system, RuntimeError(f"skipped path: missing input geoms/{system}_path.json (geometry pack)"), note="skipped")
            continue
        print(f"  {system}: {len(pts)} path points")
        for q in pts:
            with item_guard("phase03b", f"{system}:p{q['index']}"):
                cc_row(csv_, system, f"p{q['index']:02d}", LV_SCAN, q["atoms"], 0, basis_for(system, "main"), f"ccsdt:{system}:path", coordinate=f"{q['coordinate']:.4f}")


def c3c_mp2_stationary(systems=None):
    csv_ = PHASE_CSV["phase08"]
    for system in [s for s in CPU_MP2_SYSTEMS if s in (systems or SYSTEMS_TO_RUN)]:
        print(f"  {system}: MP2 basis correction at the stationary points")
        for sp in SPECIES_KIND[system]:
            lv, g = stationary_geom(system, sp)
            if g is None:
                continue
            for size in ("main", "large"):
                with item_guard("phase08", f"{system}:{sp}:{size}"):
                    mp2_row(csv_, system, sp, lv, g["atoms"], species_charge(system, sp), size, f"mp2:{system}:{size}")


def c4_scout_reference(systems=None):
    csv_ = PHASE_CSV["phase04"]
    for system in [s for s in CPU_SCOUT_REF_SYSTEMS if s in (systems or SYSTEMS_TO_RUN)]:
        S = SYSTEMS[system]
        rmin = S["rc"] or S["reactants"][0]
        for sp in (rmin, "TS"):
            g = geometry_of(system, sp, LV_XTB)
            if g is None or g.get("no_peak"):
                continue
            with item_guard("phase04", f"{system}:{sp}"):
                cc_row(csv_, system, sp, LV_XTB, g["atoms"], species_charge(system, sp), basis_for(system, "main"), f"ccsdt:{system}:{sp}")


def c6_all_electron(systems=None):
    csv_ = PHASE_CSV["phase06"]
    if "hcn" not in (systems or SYSTEMS_TO_RUN):
        return
    for sp in ("HCN", "TS"):
        lv, g = best_geom("hcn", sp)
        if g is not None:
            with item_guard("phase06", f"hcn:{sp}:ae"):
                cc_row(csv_, "hcn", sp, lv, g["atoms"], 0, basis_for("hcn", "main"), "ccsdt:hcn:ae", method=M_CCSDT_AE, fc=False, stab=False)


# ------------------------------------------------------------------ the tail (only if time remains)
def tail_mp2_paths(systems=None):
    """MP2 basis correction at the 31 path points of HCN, NH3, H2CO (CH4: stays raw)."""
    csv_ = PHASE_CSV["phase08"]
    for system in [s for s in ("hcn", "nh3", "h2co") if s in (systems or SYSTEMS_TO_RUN)]:
        pts = load_path_pts(system)
        print(f"  {system}: MP2 correction at {len(pts)} path points")
        for q in pts:
            for size in ("main", "large"):
                with item_guard("phase08", f"{system}:p{q['index']}:{size}"):
                    mp2_row(csv_, system, f"p{q['index']:02d}", LV_SCAN, q["atoms"], 0, size, f"mp2:{system}:path:{size}", coordinate=f"{q['coordinate']:.4f}")


def grid_npz(theta, r, basis):
    E = np.full((len(theta), len(r)), np.nan); conv = np.zeros((len(theta), len(r)), bool); t1 = np.full((len(theta), len(r)), np.nan)
    for x in read_rows(PHASE_CSV["phase05"]):
        _, i, j = x["point_id"].split("_")
        if x["energy_hartree"] not in ("", "nan"):
            E[int(i), int(j)] = float(x["energy_hartree"]); conv[int(i), int(j)] = x["converged"] in ("1", "True"); t1[int(i), int(j)] = float(x["t1"] or "nan")
    np.savez(path("phase05_hcn_grid.npz"), E=E, converged=conv, t1=t1, theta=theta, r=r, basis=basis, method=M_CCSDT_FC)
    return E


def tail_hcn_grid(stride=(1, 1)):
    """CCSD(T)/def2-SVP on the 25 x 17 rigid HCN grid of v4 (425 points, about 0.4 s each in the sandbox)."""
    theta, r = hcn_grid_axes(25, 17)
    basis = "def2-SVP"
    order = [(i, j) for i in range(0, 25, stride[0]) for j in range(0, 17, stride[1])]
    try:
        for n, (i, j) in enumerate(order):
            with item_guard("phase05", f"g_{i}_{j}"):
                qrun(PHASE_CSV["phase05"], "hcn_grid", f"g_{i}_{j}", LV_GRID, M_CCSDT_FC, basis,
                     lambda i=i, j=j: ccsdt(hcn_geometry(theta[i], r[j]), 0, basis, stability=False), "ccsdt:hcn_grid",
                     coordinate=f"{theta[i]:.1f},{r[j]:.4f}", retry_failed=RETRY_FAILED, label=f"hcn grid {i},{j}", quiet=(n % 50 != 0))
    finally:
        E = grid_npz(theta, r, basis)
        print(f"  grid points with an energy: {int(np.isfinite(E).sum())} of {E.size}")


def ext_clcl_ccsdt():
    """EXT-CPU-1: CCSD(T) for the Cl- + CH3Cl stationary points (def2-TZVPD, frozen core)."""
    for sp in SPECIES_KIND["sn2_clcl"]:
        lv, g = stationary_geom("sn2_clcl", sp)
        if g is not None:
            with item_guard("phase03x", f"sn2_clcl:{sp}"):
                cc_row(PHASE_CSV["phase03"], "sn2_clcl", sp, lv, g["atoms"], species_charge("sn2_clcl", sp), basis_for("sn2_clcl", "main"), f"ccsdt:sn2_clcl:{sp}")
    for sp in SPECIES_KIND["sn2_clcl"]:                              # the MP2 basis correction of this reference (cheap: about 0.5 minute in the sandbox)
        lv, g = stationary_geom("sn2_clcl", sp)
        if g is not None:
            for size in ("main", "large"):
                with item_guard("phase03x", f"sn2_clcl:{sp}:mp2:{size}"):
                    mp2_row(PHASE_CSV["phase08"], "sn2_clcl", sp, lv, g["atoms"], species_charge("sn2_clcl", sp), size, f"mp2:sn2_clcl:{size}")
    rows = read_rows(PHASE_CSV["phase03"])
    have_ = {r["point_id"] for r in rows if r["system"] == "sn2_clcl" and r["method"] == M_CCSDT_FC and r["converged"] in ("1", "True")}
    if set(SPECIES_KIND["sn2_clcl"]) <= have_:                      # the reference now exists: Cl- + CH3Cl moves from Tier C to Tier B (check R1 / R2 and the headline set see it)
        tiers = read_json("tiers.json", {}) or {k: v["tier"] for k, v in SYSTEMS.items()}
        tiers["sn2_clcl"] = "B"
        write_json("tiers.json", tiers)
        print("  Cl- + CH3Cl: CCSD(T) complete at the stationary points: tier C -> B")

