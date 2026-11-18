"""Test of the automated check (b1v5_check.py) on SYNTHETIC CSV files.

It writes made-up energies into a temporary folder, designed so that some checks PASS, some WARN and some FAIL, runs the
check, and asserts the expected statuses.  Nothing here is a chemistry result.
Run:  python check_test_synthetic.py        (needs only numpy)

After the review of 2026-11-15 the headline geometry level is our own B3LYP geometry (LV_OWN); the database level (LV_DB) is only a cross-check.
Scenarios: A (mixed), B (second set of failures), C (empty folder), D (rules of the review: one geometry level per quantity, density-fitting retry
alignment, small-reference tolerance, unverified rows, I6 / I1 not dropping a system, MP2 basis correction, Table 2 statistics).
"""
import os, sys, tempfile, json, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import b1v5_common as C
import b1v5_check as K
import b1v5_v4_stored as V4
import b1v5_lit_table as LT
import numpy as np

KC = C.HARTREE_KCAL
TZ, QZ = {"A": "def2-TZVP", "D": "def2-TZVPD"}, {"A": "def2-QZVPP", "D": "def2-QZVPPD"}


def bas(system, size):
    return C.basis_for(system, size)


class Synth:
    """Writes rows into the phase CSV files of a temporary root."""
    def __init__(self, root):
        C.setup_root(root, mount=False)
        self.n = 0

    def put(self, csv_name, system, pid, level, method, basis, kcal, conv=True, t1="", imag="", extra=None):
        S_ = C.SYSTEMS.get(system)
        frag = S_ and any(pid in lst and len(lst) > 1 for lst in (S_["reactants"], S_["products"]))
        base = -50.0 if frag else -100.0           # two separated fragments add up to the same baseline as one molecule
        C.append_row(csv_name, C.make_row(system, pid, level, method, basis, base + kcal / KC, conv, "", t1, imag, 0.1, extra))
        self.n += 1

    def block(self, csv_name, system, method, basis, level, energies, **kw):
        for sp, e in energies.items():
            self.put(csv_name, system, sp, level, method, basis, e, **kw)

    def mp2(self, system, level, deltas):
        """DF-MP2 rows in the main and the large basis; deltas = {point: (QZ minus TZ energy in kcal/mol)}  (review row 5)."""
        for pid, d in deltas.items():
            self.put("phase08_mp2_basis.csv", system, pid, level, C.M_MP2DF, bas(system, "main"), 0.0)
            self.put("phase08_mp2_basis.csv", system, pid, level, C.M_MP2DF, bas(system, "large"), d)


def build(root, scen="A"):
    B = scen == "B"
    S = Synth(root)
    m = C
    tiers = {k: v["tier"] for k, v in C.SYSTEMS.items()}
    C.write_json("tiers.json", tiers)
    P3, P11, P12, P4, P6, P1 = "phase03_ccsdt.csv", "phase11_dft.csv", "phase12_basis.csv", "phase04_scout_ref.csv", "phase06_ccsdt_basis.csv", "phase01_stationary.csv"
    FC, WB = m.M_CCSDT_FC, m.M_WB97MV
    # ---------------- HCN (tier A): everything fine, except in scenario B (unconverged point, own vs database geometry 1.2, nb4 rerun off)
    hcn = lambda ts, p: {"HCN": 0.0, "TS": ts, "HNC": p}
    S.block(P3, "hcn", FC, TZ["A"], m.LV_DB, hcn(48.9, 15.0), t1=0.015, extra={"stable": 1})
    S.block(P3, "hcn", FC, TZ["A"], m.LV_OWN, hcn(48.9 + (1.2 if B else 0.1), 15.0), t1=0.015, extra={"stable": 1})
    S.block(P6, "hcn", FC, QZ["A"], m.LV_OWN, hcn(49.2, 15.2))
    S.block(P6, "hcn", m.M_CCSDT_AE, TZ["A"], m.LV_OWN, {"HCN": 0.0, "TS": 49.0 + (0.5 if B else 0.0)})
    S.block(P11, "hcn", WB, TZ["A"], m.LV_OWN, hcn(48.2, 15.2)); S.block(P12, "hcn", WB, QZ["A"], m.LV_OWN, hcn(48.5, 15.2))
    for meth, ts in ((m.M_SKALA, 48.0), (m.M_SKALA_NOD3, 48.3), (m.M_PBE, 46.7), (m.M_R2SCAN, 47.5), (m.M_B3LYP, 47.8)):
        S.block(P11, "hcn", meth, TZ["A"], m.LV_OWN, hcn(ts, 15.0))
    S.block(P11, "hcn", WB, TZ["A"], m.LV_DB, hcn(48.3, 15.2))                       # database-geometry cross-check rows
    S.mp2("hcn", m.LV_OWN, {"HCN": 0.0, "TS": -0.4, "HNC": -0.1}); S.mp2("hcn", m.LV_DB, {"HCN": 0.0, "TS": -0.4, "HNC": -0.1})
    S.block(P1, "hcn", m.M_DFB3LYP, TZ["A"], m.LV_OWN, {"HCN": 0.0, "HNC": 15.0}, imag=0)
    S.put(P1, "hcn", "TS", m.LV_OWN, m.M_DFB3LYP, TZ["A"], 48.0, imag=1)
    if B:
        S.put(P11, "hcn", "HNC", m.LV_OWN, m.M_PBE, TZ["A"], 15.0, conv=False)              # unconverged PBE point: removes PBE's entry only (review row 3)
        S.put(P1, "hcn", "H", m.LV_OWN, m.M_DFB3LYP, TZ["A"], 0.0, conv=False)              # unconverged DF-B3LYP (geometry) point: a core method fails the system
    for i, e in enumerate([0, 5, 15, 30, 44, 49.0, 44, 30, 20, 17, 15]):
        S.put(P3, "hcn", f"p{i:02d}", m.LV_SCAN, FC, TZ["A"], e, t1=0.015, extra={"stable": 1})
    S.block(P4, "hcn", FC, TZ["A"], m.LV_XTB, {"HCN": 0.0, "TS": 48.9 + 0.4}, t1=0.015, extra={"stable": 1})
    S.mp2("hcn", m.LV_SCAN, {f"p{i:02d}": 0.0 for i in range(11)})
    S.put("phase00_repro.csv", "hcn_nb4", "nb4_00", m.LV_NB4, "v4-ccsd_t", "def2-TZVP", 0.0)
    # replace the dummy energy by the stored v4 energy (+3e-6, or +2e-4 in scenario B)
    rows = open(C.path("phase00_repro.csv")).read().splitlines()
    cells = rows[1].split(","); cells[6] = repr(V4.NB4_HCN_E["ccsd_t"][0] + (2e-4 if B else 3e-6)); rows[1] = ",".join(cells)
    open(C.path("phase00_repro.csv"), "w").write("\n".join(rows) + "\n")
    # ---------------- NH3 (tier A): CCSD(T) barrier 7.0 (9.0 in B), T1 0.025, double-peaked path, basis shift 0.8
    nh = lambda ts: {"NH3": 0.0, "TS": ts}
    S.block(P3, "nh3", FC, TZ["A"], m.LV_OWN, nh(9.0 if B else 7.0), t1=0.025, extra={"stable": 1})
    S.mp2("nh3", m.LV_OWN, {"NH3": 0.0, "TS": -0.8})
    S.block(P11, "nh3", WB, TZ["A"], m.LV_OWN, nh(5.2)); S.block(P12, "nh3", WB, QZ["A"], m.LV_OWN, nh(6.0))
    S.put(P1, "nh3", "NH3", m.LV_OWN, m.M_DFB3LYP, TZ["A"], 0.0, imag=0); S.put(P1, "nh3", "TS", m.LV_OWN, m.M_DFB3LYP, TZ["A"], 5.0, imag=1)
    for i, e in enumerate([0, 3, 7, 4, 2, 5, 7.2, 5, 3, 1, 0]):
        S.put(P3, "nh3", f"p{i:02d}", m.LV_SCAN, FC, TZ["A"], e, t1=0.015, extra={"stable": 1})
    S.block(P4, "nh3", FC, TZ["A"], m.LV_XTB, {"NH3": 0.0, "TS": (9.0 if B else 7.0) + 2.0}, t1=0.015, extra={"stable": 1})
    # ---------------- H2CO (tier A): zero-point INCLUDED in the literature value
    zp = lambda sp: {"H2CO": 0.0, "TS": -3.0, "H2": 0.0, "CO": 0.0}[sp]
    for sp, e in {"H2CO": 0.0, "TS": 85.0, "H2": 0.0, "CO": 0.0}.items():
        S.put(P3, "h2co", sp, m.LV_OWN, FC, TZ["A"], e, t1=0.012, extra={"stable": 1})
        S.put(P1, "h2co", sp, m.LV_OWN, m.M_DFB3LYP, TZ["A"], e, imag=(1 if sp == "TS" else 0), extra={"zpe_hartree": 0.026 + zp(sp) / KC})
    S.block(P11, "h2co", WB, TZ["A"], m.LV_OWN, {"H2CO": 0.0, "TS": (-5.0 if B else 86.6), "H2": 0.0, "CO": 0.0}); S.block(P12, "h2co", WB, QZ["A"], m.LV_OWN, {"H2CO": 0.0, "TS": 88.0, "H2": 0.0, "CO": 0.0})
    S.mp2("h2co", m.LV_OWN, {"H2CO": 0.0, "TS": -0.5, "H2": 0.0, "CO": 0.0})
    # ---------------- CH3NC (tier B): zero-point convention unknown; TS with TWO imaginary modes; scout penalty 3.5
    for sp, e in {"CH3NC": 0.0, "TS": 38.9, "CH3CN": -20.0}.items():
        S.put(P3, "ch3nc", sp, m.LV_OWN, FC, TZ["A"], e, t1=0.015, extra={"stable": 1})
        S.put(P1, "ch3nc", sp, m.LV_OWN, m.M_DFB3LYP, TZ["A"], e, imag=(2 if sp == "TS" else 0), extra={"zpe_hartree": 0.05 + (-1.0 / KC if sp == "TS" else 0.0)})
    S.block(P11, "ch3nc", WB, TZ["A"], m.LV_OWN, {"CH3NC": 0.0, "TS": 40.5, "CH3CN": -19.0}); S.block(P12, "ch3nc", WB, QZ["A"], m.LV_OWN, {"CH3NC": 0.0, "TS": 41.0, "CH3CN": -19.0})
    S.mp2("ch3nc", m.LV_OWN, {"CH3NC": 0.0, "TS": -0.3, "CH3CN": 0.0})
    S.block(P4, "ch3nc", FC, TZ["A"], m.LV_XTB, {"CH3NC": 0.0, "TS": 38.9 + 3.5}, t1=0.015, extra={"stable": 1})
    # ---------------- F- + CH3F (tier B): CCSD(T) fine, wB97M-V barrier too high (R2 FAIL), Skala lowest of ladder (I7 WARN)
    sn = lambda ts: {"F-": 0.0, "CH3F": 0.0, "RC": -13.0, "TS": ts}
    S.block(P3, "sn2_ff", FC, TZ["D"], m.LV_OWN, sn(0.5), t1=0.014, extra={"stable": 1})
    S.mp2("sn2_ff", m.LV_OWN, {"F-": 0.0, "CH3F": 0.0, "RC": 0.0, "TS": -0.4})
    for meth, ts in ((m.M_SKALA, 10.0), (m.M_PBE, 12.0), (m.M_R2SCAN, 12.5), (m.M_B3LYP, 13.5)):
        S.block(P11, "sn2_ff", meth, TZ["D"], m.LV_OWN, sn(ts - 13.0))
    S.block(P11, "sn2_ff", WB, TZ["D"], m.LV_OWN, sn(3.8)); S.block(P12, "sn2_ff", WB, QZ["D"], m.LV_OWN, sn(4.0))
    # ---------------- Cl- + CH3Cl (tier C): wB97M-V + ladder; PBE lowest (I7 PASS); xTB scan without a peak (I12 WARN)
    cl = lambda ts: {"Cl-": 0.0, "CH3Cl": 0.0, "RC": -9.0, "TS": ts}
    for meth, ts in ((m.M_SKALA, 2.5), (m.M_PBE, -12.0 if B else 1.0), (m.M_R2SCAN, 2.0), (m.M_B3LYP, 3.0), (WB, 3.5)):
        S.block(P11, "sn2_clcl", meth, TZ["D"], m.LV_OWN, cl(ts))
    S.block(P12, "sn2_clcl", WB, QZ["D"], m.LV_OWN, cl(3.3))
    if B:
        S.put(P11, "sn2_clcl", "RC", m.LV_OWN, m.M_R2SCAN, TZ["D"], -9.0, conv=False)           # unconverged r2SCAN point: only r2SCAN's entry goes
    for sp, e, im in (("Cl-", 0.0, 0), ("CH3Cl", 0.0, 0), ("RC", -9.0, 0), ("TS", 3.0, 1)):
        S.put(P1, "sn2_clcl", sp, m.LV_OWN, m.M_DFB3LYP, TZ["D"], e, imag=im)
    S.put(P1, "sn2_ff", "TS", m.LV_DB, m.M_DFB3LYP, TZ["D"], 0.5, imag=1)                      # database geometry only: I4 WARN (not tested)
    C.append_row(P4, C.make_row("sn2_clcl", "TS", m.LV_XTB, m.M_XTB, "GFN2", float("nan"), False, extra={"no_peak": True}))
    # ---------------- CH4: T1 0.04 and unstable RHF at the last point (I2 FAIL, I3 FAIL)
    for i in range(9):
        bad = i == 8
        S.put(P3, "ch4", f"p{i:02d}", m.LV_SCAN, FC, TZ["A"], 10.0 * i, t1=(0.04 if bad else 0.01), extra={"stable": 0 if bad else 1})
    return S


def st(rep, item, scope_start):
    return [x["status"] for x in rep["items"] if x["item"] == item and x["scope"].startswith(scope_start)]


def expect(rep, item, scope_start, status):
    got = st(rep, item, scope_start)
    ok = got == [status]
    print(("  ok  " if ok else "  BAD ") + f"{item} {scope_start}: expected {status}, got {got}")
    return ok


def expect_set(rep, item, scope_start, statuses):
    got = sorted(st(rep, item, scope_start))
    ok = got == sorted(statuses)
    print(("  ok  " if ok else "  BAD ") + f"{item} {scope_start}: expected {sorted(statuses)}, got {got}")
    return ok


def check(cond, text):
    print(("  ok  " if cond else "  BAD ") + text)
    return bool(cond)


def scenario_d():
    """The rules of the review of 2026-11-15, tested on the library functions and on the Table 2 code."""
    ok = True
    # ---- D1: one geometry level per quantity (row 4)
    tmp = tempfile.mkdtemp(prefix="b1v5_synth_D1_")
    S = Synth(tmp)
    C.write_json("tiers.json", {k: v["tier"] for k, v in C.SYSTEMS.items()})
    sysid = "sn2_ohf"; SS = C.SYSTEMS[sysid]
    bD = bas(sysid, "main")
    sps = list(dict.fromkeys(SS["reactants"] + [SS["rc"], SS["ts"]]))
    for sp in sps:
        S.put("phase11_dft.csv", sysid, sp, C.LV_OWN if sp != SS["ts"] else C.LV_DB, C.M_WB97MV, bD, 0.0 if sp != SS["ts"] else 5.0)
    D = K.Data()
    ok &= check(K.quantity(D, sysid, "barrier_fwd_cplx", C.M_WB97MV, bD) is None, "TS at the database level, minima at the own level: barrier is None (NOT RUN), not mixed")
    S.put("phase11_dft.csv", sysid, SS["ts"], C.LV_OWN, C.M_WB97MV, bD, 5.0)
    D = K.Data()
    v = K.quantity(D, sysid, "barrier_fwd_cplx", C.M_WB97MV, bD)
    ok &= check(v is not None and abs(v - 5.0) < 1e-6, f"after the own-level TS is added the barrier is {v}")
    i = K.quantity_info(D, sysid, "barrier_fwd_cplx", C.M_WB97MV, bD)
    ok &= check(i["level"] == C.LV_OWN, f"level written with the value: {i['level']}")
    ok &= check(K.quantity(D, sysid, "barrier_fwd_cplx", C.M_WB97MV, bD, level="db") is None, "database family requested explicitly: not complete, None")
    # ---- D2: density-fitting retry alignment (row 11)
    sysid = "sn2_fcl"; SS = C.SYSTEMS[sysid]; bD = bas(sysid, "main")
    for sp in list(dict.fromkeys(SS["reactants"] + [SS["rc"], SS["ts"]])):
        S.put("phase11_dft.csv", sysid, sp, C.LV_OWN, C.M_WB97MV, bD, 0.0 if sp != SS["ts"] else 6.0, extra={"route": "gpu", "df": bool(sp == SS["ts"])})
    D = K.Data()
    ok &= check(K.quantity(D, sysid, "barrier_fwd_cplx", C.M_WB97MV, bD) is None, "fitted TS with unfitted partners and no recomputed partner: barrier dropped")
    R = K.Report(); K.i14_retry_mode(D, R)
    ok &= check(any(x["item"] == "I14" and x["status"] == K.WARN for x in R.items), "I14 WARN names the mixed barrier")
    for sp in list(dict.fromkeys(SS["reactants"] + [SS["rc"]])):
        S.put("phase11_dft.csv", sysid, sp + "@df", C.LV_OWN, C.M_WB97MV, bD, 0.2, extra={"route": "gpu", "df": True})
    D = K.Data()
    v = K.quantity(D, sysid, "barrier_fwd_cplx", C.M_WB97MV, bD)
    ok &= check(v is not None and abs(v - 5.8) < 1e-6, f"with the partners recomputed with density fitting the barrier uses them: {v} (6.0 - 0.2 expected)")
    R = K.Report(); K.i14_retry_mode(D, R)
    ok &= check(any(x["item"] == "I14" and x["status"] == K.INFO and "recomputed" in x["message"] for x in R.items), "I14 INFO records the alignment")
    shutil.rmtree(tmp, ignore_errors=True)
    # ---- D3: tolerance for small references (row 7)
    t = (1.0, 2.0)
    for diff, ref, want in ((0.6, 5.12, K.WARN), (0.4, 5.12, K.PASS), (1.0, 40.0, K.PASS), (0.6, 3.0, K.WARN), (1.5, 40.0, K.WARN), (2.5, 40.0, K.FAIL)):
        g = K.grade(diff, t, ref)
        ok &= check(g == want, f"grade(diff {diff}, ref {ref}) = {g}, expected {want}")
    # ---- D4: literature rows (row 6) and Table 2 statistics (row 8) on scenario A data
    tmp = tempfile.mkdtemp(prefix="b1v5_synth_D4_")
    build(tmp, "A")
    ns = dict(K.__dict__)
    exec(compile(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "b1v5_gpu_figs.py")).read(), "b1v5_gpu_figs.py", "exec"), ns)
    D = ns["Data"]()
    L = {l["id"]: l for l in LT.LIT}
    ok &= check(L["L6"]["verified"] == "no" and not LT.lit_in_stats(L["L6"]), "L6 (CH3NC 38.25): verified no, excluded from statistics")
    ok &= check("L19" in L and L["L19"]["value"] is None and not LT.lit_in_stats(L["L19"]), "L19 placeholder exists: no value, excluded from statistics")
    ok &= check(not LT.lit_in_stats(L["L3"]) and LT.lit_scorable(L["L3"]), "L3 (L1 minus L2) is scored in the supplement but not in pooled statistics")
    for m in (C.M_SKALA, C.M_B3LYP, C.M_CCSDT_FC, C.M_CCSDT_BC):
        ids = [r["lit_id"] for r in ns["pooled_rows"](D, m)]
        ok &= check("L3" not in ids and "L6" not in ids and "L19" not in ids, f"pooled rows of {m} contain neither L3, L6 nor L19: {ids}")
    hs = ns["headline_rows"](D, C.M_SKALA)
    systems = [r["system"] for r in hs]
    ok &= check(len(systems) == len(set(systems)) and "ch3nc" not in systems and set(systems) <= {"hcn", "nh3", "h2co", "sn2_ff", "sn2_clcl", "sn2_fcl", "sn2_ohf"},
                f"headline set has one barrier per reaction, no CH3NC: {systems}")
    ok &= check(all(r["level"].startswith("own_b3lyp") for r in hs), "headline rows carry the geometry level 'own_b3lyp'")
    line = ns["_stat_line"](C.M_SKALA, "main", "headline", hs)
    ok &= check(line["n_reactions"] == len(systems) and line["n_rows"] == len(hs), f"n printed beside the mean: n_reactions {line['n_reactions']}, n_rows {line['n_rows']}")
    ok &= check(ns["TIE"] == 0.5, "ties are defined as within 0.5 kcal/mol")
    e_bc = {r["system"]: r["error"] for r in ns["headline_rows"](D, C.M_CCSDT_BC)}
    e_raw = {r["system"]: r["error"] for r in ns["headline_rows"](D, C.M_CCSDT_FC)}
    ok &= check(abs((e_raw["hcn"] - e_bc["hcn"]) - 0.4) < 1e-6, f"basis-corrected CCSD(T) error differs from the raw one by the MP2 shift: {e_raw['hcn']:.2f} vs {e_bc['hcn']:.2f}")
    # the raw reference stays as a column of check R1
    rep = K.run_check(tmp, write=False)
    r1 = [x for x in rep["items"] if x["item"] == "R1" and x["scope"].startswith("L1 ")][0]
    ok &= check(r1["raw_value"] != "" and abs(float(r1["raw_value"]) - float(r1["value"]) - 0.4) < 1e-6, f"R1 L1 value {r1['value']} (corrected) and raw column {r1['raw_value']}")
    shutil.rmtree(tmp, ignore_errors=True)
    return ok


def main():
    allok = True
    # ------------------------------------------------ scenario A
    tmp = tempfile.mkdtemp(prefix="b1v5_synth_A_")
    S = build(tmp, "A")
    print(f"\n##### scenario A: {S.n} synthetic rows in {tmp}\n")
    rep = K.run_check(tmp)
    print("\n##### assertions, scenario A")
    exp = [("R1", "L1 ", "PASS"), ("R1", "L2 ", "PASS"), ("R1", "L3 ", "PASS"), ("R1", "L4 ", "WARN"), ("R1", "L5 ", "PASS"), ("R1", "L6 ", "WARN"),
           ("R1", "L7 ", "PASS"), ("R1", "L8 ", "PASS"), ("R1", "L9 ", "N/A"),
           ("R2", "L1 ", "PASS"), ("R2", "L4 ", "WARN"), ("R2", "L5 ", "WARN"), ("R2", "L6 ", "WARN"), ("R2", "L7 ", "FAIL"), ("R2", "L8 ", "WARN"),
           ("R2", "L9 ", "PASS"), ("R2", "L10 ", "PASS"), ("R2", "L11 ", "NOT RUN"), ("R3", "L6 ", "WARN"), ("R4", "L1 ", "INFO"),
           ("I1", "hcn", "PASS"), ("I1", "sn2_fcl", "NOT RUN"), ("I2", "hcn", "PASS"), ("I2", "nh3", "WARN"), ("I2", "ch4", "FAIL"), ("I2", "sn2_clcl", "N/A"),
           ("I3", "ch4", "FAIL"), ("I3", "hcn", "PASS"),
           ("I4", "hcn", "PASS"), ("I4", "ch3nc", "FAIL"), ("I4", "sn2_ff", "WARN"), ("I5", "hcn", "PASS"), ("I5", "nh3", "FAIL"), ("I5", "h2co", "NOT RUN"),
           ("I6", "hcn", "PASS"), ("I6", "h2co", "PASS"), ("I7", "hcn", "PASS"), ("I7", "sn2_ff", "WARN"), ("I7", "sn2_clcl", "PASS"),
           ("I8", "hcn", "PASS"), ("I8", "nh3", "WARN"), ("I8", "h2co", "FAIL"), ("I9", "v4 energies", "PASS"), ("I9", "v4 scout", "PASS"), ("I9", "v4 path", "NOT RUN"),
           ("I10", "hcn", "PASS"), ("I11", "hcn own", "PASS"), ("I12", "hcn", "PASS"), ("I12", "nh3", "WARN"), ("I12", "ch3nc", "FAIL"), ("I12", "sn2_clcl", "WARN"),
           ("I13", "hcn", "INFO")]
    for e in exp:
        allok &= expect(rep, *e)
    dec = {k: v[0] for k, v in rep["decisions"].items()}
    want = dict(hcn="MAIN", nh3="MAIN", h2co="MAIN", ch3nc="SUPPLEMENT", sn2_ff="SUPPLEMENT", sn2_clcl="MAIN", sn2_fcl="INCOMPLETE", sn2_ohf="INCOMPLETE")
    for k, v in want.items():
        ok = dec[k] == v; allok &= ok
        print(("  ok  " if ok else "  BAD ") + f"decision {k}: expected {v}, got {dec[k]}")
    ok = rep["paper"].startswith("4 of 8") and "PROVISIONAL" in rep["paper"]; allok &= ok
    print(("  ok  " if ok else "  BAD ") + "paper level:", rep["paper"])
    for f in ("check_report.csv", "check_report.md"):
        ok = os.path.exists(os.path.join(tmp, f)); allok &= ok; print(("  ok  " if ok else "  BAD ") + f, "written")
    # ------------------------------------------------ scenario B (a second set of failures)
    tmp2 = tempfile.mkdtemp(prefix="b1v5_synth_B_")
    build(tmp2, "B")
    print("\n##### scenario B (more failures)\n")
    repB = K.run_check(tmp2)
    print("\n##### assertions, scenario B")
    for e in [("R1", "L4 ", "FAIL"), ("I1", "hcn", "FAIL"), ("I9", "v4 energies", "FAIL"), ("I10", "hcn", "WARN"), ("I11", "hcn own", "FAIL")]:
        allok &= expect(repB, *e)
    # review row 3: an unconverged PBE point or a PBE barrier below zero removes that method's entry, never the system
    allok &= expect_set(repB, "I6", "sn2_clcl", ["PASS", "INFO"])
    allok &= check(any(x["item"] == "I1" and x["scope"].startswith("sn2_clcl") and x["status"] == "WARN" for x in repB["items"]), "I1 sn2_clcl: the unconverged r2SCAN point is a WARN (that method entry is dropped), not a FAIL of the system")
    decB = {k: v[0] for k, v in repB["decisions"].items()}
    for k, v in dict(hcn="SUPPLEMENT", nh3="MAIN, CCSD(T) role withdrawn", sn2_clcl="MAIN").items():
        ok = decB[k] == v; allok &= ok
        print(("  ok  " if ok else "  BAD ") + f"decision {k}: expected {v}, got {decB[k]}")
    # ------------------------------------------------ scenario C: an empty folder must not crash
    tmp3 = tempfile.mkdtemp(prefix="b1v5_synth_C_")
    print("\n##### scenario C: empty folder (nothing run yet)\n")
    repC = K.run_check(tmp3)
    ok = repC["counts"]["FAIL"] == 0 and repC["counts"]["PASS"] <= 3 and all(v[0] == "INCOMPLETE" for k, v in repC["decisions"].items() if k != "ch4")
    allok &= ok
    print(("  ok  " if ok else "  BAD ") + f"empty folder: counts {repC['counts']}, no crash, every system INCOMPLETE")
    print("\n##### scenario D: rules of the review (levels, retry alignment, tolerance, literature rows, Table 2)\n")
    allok &= scenario_d()
    print("\nRESULT:", "ALL ASSERTIONS PASSED" if allok else "SOME ASSERTIONS FAILED")
    for t in (tmp, tmp2, tmp3):
        shutil.rmtree(t, ignore_errors=True)
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
