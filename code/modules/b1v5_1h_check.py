# ==============================================================================
# b1v5_1h_check.py  --  the automated check and the frozen criteria for the SIX-REACTION design of PLAN_1hour.md.
# Pasted after b1v5_check.py (and, in the GPU notebook, after b1v5_gpu_figs.py), in the same namespace.  It does not edit the full-size files:
# it REDEFINES a few functions (i11_geom, decide, paper_level, run_check, findings_outcome) and the FROZEN table.
#
# What changes for six reactions (plan section 5: "F1 and F5 thresholds become 5 of 6 where they were 6 of 8; fix that before the run"):
#   * the check and the figures see only the design (set_design in b1v5_1h.py rebinds REACTIONS, SYSTEM_ORDER and LIT);
#   * every item that cannot exist in the reduced design is printed as NOT RUN (never FAIL) and says why:
#       - literature rows L11 to L18 (F- + CH3Cl, OH- + CH3F): "outside the 1-hour core design (extension EXT-GPU-1)";
#       - I11 own versus database geometry: "dropped by design (decision D2: one geometry protocol, no database download)";
#       - I11 Skala//Skala: "Skala-relaxed geometries are not part of the core (extension EXT-GPU-2)";
#       - I8 / F3 where no QZVPP energies exist (CH3NC and the SN2 systems in the core), I9 path errors (only the first stored v4 point is rerun),
#         I10 (all-electron CCSD(T) for HCN only), CCSD(T) references that the queue did not reach (tier falls back to C): all NOT RUN;
#   * decisions: the two extension reactions get "NOT RUN (extension)", never INCOMPLETE or SUPPLEMENT;
#   * paper-level rule: n of the 6 design reactions (rule: at least 6 reactions in at least 3 classes: the core has exactly six, so one failure turns
#     the full paper into a short communication, as the plan says in section 6);
#   * FROZEN: thresholds are the same 75% / 50% rules, written for the design (F1, F5: 5 of 6 reactions when all six have a verified reference,
#     4 of 5 while CH3NC L6 is unverified; F4: 5 of 6; F3: at least 3 reactions instead of 4, because the core has the basis check for three).
# ==============================================================================
import math
import numpy as np

if "_decide_orig" not in globals():
    _decide_orig = decide                    # noqa  (the full-size rule 4.4, kept)
    _run_check_orig = run_check              # noqa

if "_load_rows_orig" not in globals():
    _load_rows_orig = load_rows              # noqa


def load_rows():
    """As b1v5_check.load_rows, without the gradient-test rows of phase00_verify.csv (method '...-grad'): they are a test, not an energy."""
    rows, files = _load_rows_orig()
    return [r for r in rows if not str(r["method"]).endswith("-grad")], files


NOT_IN_DESIGN = "outside the 1-hour core design (extension session EXT-GPU-1 adds F- + CH3Cl and OH- + CH3F)"


def pack_items(R):
    """P1: SHA-256 of the geometry pack; P2: the gradient test of C0 (RMS gradient below 3e-4 hartree/bohr at the stored geometries)."""
    info = read_json("geompack_installed.json", None)
    if not info:
        R.add("P1", "geometry pack", NOTRUN, msg="geompack_installed.json not found (the pack was not installed, or the geometries were rebuilt: REBUILD_GEOMS)")
    else:
        ok = bool(info.get("ok"))
        R.add("P1", "geometry pack", PASS if ok else FAIL, str(info.get("sha256", ""))[:16], str(info.get("expected_sha256", ""))[:16], "", "SHA-256 equal",
              ("pack and all %s structure/path hashes verified" % info.get("n_verified")) if ok else "pack did not verify: " + "; ".join(info.get("problems", []))[:200])
    try:
        g = grad_summary()
    except Exception:
        g = dict(n=0)
    if not g["n"]:
        R.add("P2", "gradient test", NOTRUN, msg="phase00_verify.csv has no gradient rows (C0 not run, or geometries rebuilt)")
    else:
        bad = g["failed"]
        R.add("P2", "gradient test", PASS if not bad else FAIL, "%.1e" % g["worst_rms"], "", "", "RMS < %g" % g["tolerance"],
              f"{g['passed']} of {g['n']} stored structures pass; worst {g['worst_rms']:.1e} at {g['worst_where']}" + ("" if not bad else "; the stored geometries are not stationary on this machine: treat every barrier with caution (EXT-geom rebuilds them)"))


def design_notes(D, R):
    """NOT RUN lines for everything the six-reaction design cannot do (review: NOT RUN, never FAIL)."""
    pack_items(R)
    for lit in lit_out_of_design():
        tag = f"{lit['id']} {lit['system']}:{lit['quantity']}"
        R.add("R1", tag, NOTRUN, msg=NOT_IN_DESIGN)
        R.add("R2", tag, NOTRUN, msg=NOT_IN_DESIGN)
    for s in ALL_REACTIONS:
        if s not in DESIGN["reactions"]:
            for item in ("I1", "I4", "I5", "I6", "I8", "I12"):
                R.add(item, s, NOTRUN, msg=NOT_IN_DESIGN)


def i11_geom_1h(D, R):
    """I11 in the 1-hour design: the database comparison is dropped by decision D2; Skala//Skala only if the extension G4 ran."""
    R.add("I11", "own vs database geometry", NOTRUN, msg="dropped by design (decision D2 of PLAN_1hour: one geometry protocol, own B3LYP geometries; the database download is untested and blocked in the sandbox)")
    any_ = False
    for s in REACTIONS:
        mb = basis_for(s, "main")
        q = "barrier_fwd_cplx" if SYSTEMS[s]["sn2"] else "barrier_fwd"
        sk = quantity(D, s, q, M_SKALA, mb, level=LV_SKALA)
        rf = quantity(D, s, q, M_SKALA, mb)
        if sk is not None and rf is not None:
            any_ = True
            st = grade(sk - rf, TOL["I11_skala"])
            R.add("I11", f"{s} Skala//Skala", st, _fnum(sk), _fnum(rf), _fnum(sk - rf), f"{TOL['I11_skala'][0]}/{TOL['I11_skala'][1]}",
                  "Skala//Skala (TS = parabola vertex) minus Skala//reference barrier" + ("; report both barriers" if st == FAIL else ""))
    if not any_:
        R.add("I11", "Skala//Skala", NOTRUN, msg="Skala-relaxed geometries are not part of the 1-hour core (extension EXT-GPU-2, tail item 5)")


def decide_1h(D, R):
    out = _decide_orig(D, R)
    for s in ALL_REACTIONS:
        if s not in out:
            out[s] = ("NOT RUN (extension)", NOT_IN_DESIGN)
    return out


def paper_level_1h(dec):
    main = [s for s in DESIGN["reactions"] if dec[s][0].startswith("MAIN")]
    inc = [s for s in DESIGN["reactions"] if dec[s][0] == "INCOMPLETE"]
    classes = {SYSTEMS[s]["cls"] for s in main}
    n = len(main)
    if n >= 6 and len(classes) >= 3:
        v = "FULL PAPER as defined (the core has exactly six reactions: no margin)"
    elif n >= 4:
        v = "add verified extras (extension EXT-GPU-1: F- + CH3Cl, OH- + CH3F), or submit as a short communication titled as a case study"
    else:
        v = "DO NOT SUBMIT: fix the protocol"
    if inc:
        v = f"PROVISIONAL ({len(inc)} system(s) still incomplete: {', '.join(inc)}); so far: " + v
    return f"{n} of {len(DESIGN['reactions'])} design reactions in the main text, {len(classes)} classes: {v}"


def run_check_1h(root=None, write=True):
    """As run_check of b1v5_check.py, for the design: design notes first, I11 replaced, decisions and the paper-level rule for six reactions."""
    if root:
        setup_root(root, mount=False)
    print("=" * 100); print("Paper B1 v5 automated check (1-hour design, %d reactions)" % len(DESIGN["reactions"])); print("=" * 100)
    D = Data()
    print(f"CSV files found: {', '.join(D.files) or 'none'}  ({len(D.rows)} rows)\n")
    R = Report()
    for name, fn in (("literature checks R1-R4", literature_checks), ("design notes", design_notes), ("I1", i1_scf), ("I2/I3", i2_i3_cc), ("I4", i4_imag),
                     ("I5/I6", i5_i6), ("I7", i7_rank), ("I8", i8_basis), ("I9", i9_v4), ("I10", i10_fc), ("I11", i11_geom_1h), ("I12", i12_scout),
                     ("I13", i13_spread), ("I14", i14_retry_mode), ("I15", i15_basis_correction)):
        try:
            fn(D, R)
        except Exception as ex:                       # a bug in one check must never stop the others
            R.add(name.split()[0], "check code", WARN, msg=f"check crashed: {repr(ex)[:150]}")
    dec = decide_1h(D, R)
    paper = paper_level_1h(dec)
    c = R.counts()
    print("\n" + "-" * 100); print("Decision per system (rule 4.4)")
    for s in ALL_SYSTEM_ORDER:
        if s in dec:
            print(f"  {SYSTEMS[s]['n']} {SYSTEMS[s]['name']:18s} -> {dec[s][0]:30s} {dec[s][1]}")
    print("\nPaper-level rule:", paper)
    print("\nCounts: " + "  ".join(f"{k} {v}" for k, v in c.items()))
    if write:
        with open(path("check_report.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["item", "scope", "status", "value", "reference", "difference", "tolerance", "message", "raw_value", "raw_difference"]); w.writeheader(); w.writerows(R.items)
        md = ["# Check report, Paper B1 v5 (1-hour design: six reactions)", "", f"Written {now()}.  Tolerances in kcal/mol unless stated.  Items that cannot exist in this design are NOT RUN, never FAIL.", "",
              "Counts: " + ", ".join(f"{k} {v}" for k, v in c.items()), "", "## Decision per system", "", "| # | System | Tier | Decision | Why |", "|---|---|---|---|---|"]
        for s in ALL_SYSTEM_ORDER:
            if s in dec:
                md.append(f"| {SYSTEMS[s]['n']} | {SYSTEMS[s]['name']} | {D.tiers.get(s, SYSTEMS[s]['tier'])} | {dec[s][0]} | {dec[s][1]} |")
        md += ["", f"**Paper-level rule:** {paper}", "", "## All items", "", "| status | item | scope | value | reference | difference | tolerance | message | raw CCSD(T) value | raw difference |", "|---|---|---|---|---|---|---|---|---|---|"]
        for x in R.items:
            md.append(f"| {x['status']} | {x['item']} | {x['scope']} | {x['value']} | {x['reference']} | {x['difference']} | {x['tolerance']} | {str(x['message']).replace('|', '/')} | {x.get('raw_value', '')} | {x.get('raw_difference', '')} |")
        open(path("check_report.md"), "w").write("\n".join(md) + "\n")
        print("\nWrote check_report.csv and check_report.md in", path(""))
    return dict(items=R.items, counts=c, decisions=dec, paper=paper)


# redefine the names the notebooks call (the full-size names stay available as _decide_orig, _run_check_orig)
run_check = run_check_1h
decide = decide_1h
paper_level = paper_level_1h


# ------------------------------------------------------------------ frozen criteria for the six-reaction design (figures module)
FROZEN_DATE = "2026-11-15 (1-hour design, written before any real number exists)"
F3_MIN_REACTIONS = 3        # the core has the basis check (def2-QZVPP) for HCN, NH3 and H2CO only; the full design required 4
FROZEN = {
    "F1": dict(finding="Barriers at relaxed stationary points, five functionals, headline set (one primary barrier per verified reaction, SN2 from the complex; n = 5 of the 6 design reactions while the CH3NC value L6 is unverified, 6 once it is confirmed)",
               supports="Skala-1.1 (with D3) mean absolute error is at most 2.0 kcal/mol AND at most the wB97M-V (anchor) mean absolute error + 0.5",
               negative="Skala mean absolute error exceeds the better of B3LYP-D3(BJ) and r2SCAN by more than 0.5 kcal/mol",
               inconclusive="everything else (a tie within 0.5 with B3LYP or r2SCAN, or Skala above 2.0 without losing)"),
    "F2": dict(finding="Error along the relaxed Tier A paths (HCN, NH3, H2CO, CH4) against the CCSD(T)/def2-TZVP reference (basis-corrected where the MP2 correction exists)",
               supports="range (max minus min) of the signed Skala error is at most 1.0 kcal/mol on at least 3 of the 4 paths",
               negative="maximum absolute Skala error is at least 2.0 kcal/mol on at least 2 of the 4 paths",
               inconclusive="everything else; always inconclusive when fewer than 3 paths have a usable reference (R1 failed or rows missing)"),
    "F3": dict(finding="Basis-set stability of the primary barrier, def2-TZVP(D) to def2-QZVPP(D), HCN, NH3, H2CO (the core); more after EXT-GPU-1",
               supports="Skala shift at most 0.5 kcal/mol AND Skala's rank among the five functionals unchanged, in at least 75% of the reactions with data (at least 3 in the 1-hour design)",
               negative="Skala shift above 1.0 kcal/mol in at least 2 reactions",
               inconclusive="everything else, or fewer than 3 reactions with both bases"),
    "F4": dict(finding="Scout check: barriers at xTB geometries against the barrier at our own B3LYP geometries, every functional of G3 (Skala with D3, wB97M-V, B3LYP-D3(BJ) in the core)",
               supports="geometry penalty below 1.0 kcal/mol in at least 75% of the reactions (5 of 6) for EVERY functional computed",
               negative="penalty above 3.0 kcal/mol in at least 3 reactions, or in every reaction of one class (SN2: both F- + CH3F and Cl- + CH3Cl)",
               inconclusive="everything else, or fewer than 4 reactions with xTB-geometry energies"),
    "F5": dict(finding="Robustness of Skala's rank among the five functionals (by absolute error; errors within 0.5 kcal/mol of each other tie)",
               supports="the most frequent rank (ties share a rank) occurs in at least 75% of the reactions with a verified reference (5 of 6; 4 of 5 while L6 is unverified)",
               negative="the most frequent rank occurs in fewer than 50% of the reactions",
               inconclusive="everything else, or fewer than 4 reactions"),
    "F6": dict(finding="Negative results kept: CH4 C-H stretch (and the rigid HCN map against CCSD(T)/def2-SVP when the optional CPU tail ran)",
               supports="T1 below 0.02 (and RHF stable at the CH4 points) at at least 90% of the CH4 points (and of the HCN grid cells in the 120 kcal/mol window, if computed)",
               negative="that fraction is below 50%: the reference is truncated at the distance where it fails, and this is said in the paper",
               inconclusive="everything else"),
}


def findings_outcome_1h(D):
    """findings_outcome of b1v5_gpu_figs.py for the design: same arithmetic, the SN2 clause of F4 uses the SN2 systems of the design, and F3 needs 3 reactions."""
    rows = []
    fn = _functionals(D)
    anc = fn[1]
    sn2 = [s for s in REACTIONS if SYSTEMS[s]["sn2"]]

    def add(fid, measured, n, verdict, sentence, status=None):
        f = FROZEN[fid]
        rows.append(dict(finding=fid, description=f["finding"], criterion_supports=f["supports"], criterion_negative=f["negative"], criterion_inconclusive=f["inconclusive"],
                         frozen_on=FROZEN_DATE, measured=measured, n=n, verdict=verdict, sentence_for_the_paper=sentence))
    # F1
    mae = {}
    for m in fn + [M_SKALA_NOD3]:
        h = headline_rows(D, m)
        if h:
            mae[m] = (float(np.mean([abs(r["error"]) for r in h])), len(h))
    if M_SKALA in mae and anc in mae:
        S, W = mae[M_SKALA][0], mae[anc][0]
        B = min([mae[m][0] for m in (M_B3LYP, M_R2SCAN) if m in mae] or [float("nan")])
        v = "supports" if (S <= 2.0 and S <= W + 0.5) else ("negative" if (B == B and S > B + 0.5) else "inconclusive")
        add("F1", "MAE (kcal/mol, n): " + ", ".join(f"{m} {a:.2f} ({n})" for m, (a, n) in mae.items()), mae[M_SKALA][1], v,
            f"On {mae[M_SKALA][1]} reactions with a verified reference Skala-1.1 (with D3) has a mean absolute barrier error of {S:.2f} kcal/mol against {W:.2f} for {anc}"
            + (f" and {B:.2f} for the better of B3LYP-D3(BJ) and r2SCAN" if B == B else "") + ".")
    else:
        add("F1", "Skala or anchor headline rows missing", 0, "NOT RUN (data missing)", "")
    # F2
    sk_paths, n_valid = [], 0
    for s in ("hcn", "nh3", "h2co", "ch4"):
        r = _path_errors(D, s, M_SKALA, M_CCSDT_BC) or _path_errors(D, s, M_SKALA, M_CCSDT_FC)
        if r is not None:
            n_valid += 1
            sk_paths.append((s, float(r["err"].max() - r["err"].min()), float(np.max(np.abs(r["err"])))))
    if n_valid:
        flat, big = sum(1 for _, rg, _m in sk_paths if rg <= 1.0), sum(1 for _, _r, mx in sk_paths if mx >= 2.0)
        v = "inconclusive" if n_valid < 3 else ("supports" if flat >= 3 else ("negative" if big >= 2 else "inconclusive"))
        add("F2", "; ".join(f"{s}: range {rg:.2f}, max |error| {mx:.2f}" for s, rg, mx in sk_paths), n_valid, v,
            f"Along {n_valid} relaxed paths the Skala error against CCSD(T) has a range of {min(r for _, r, _ in sk_paths):.2f} to {max(r for _, r, _ in sk_paths):.2f} kcal/mol.")
    else:
        add("F2", "no Skala and CCSD(T) path rows", 0, "NOT RUN (data missing)", "")
    # F3
    shifts, rank_same, n_b = [], 0, 0
    for r in basis_table(D):
        if r["method"] == M_SKALA and r["shift_main_to_large"] is not None:
            n_b += 1
            shifts.append((r["key"], r["shift_main_to_large"]))
            if abs(r["shift_main_to_large"]) <= 0.5 and r["rank_changed"] in ("no", None):
                rank_same += 1
    if n_b:
        big = sum(1 for _, x in shifts if abs(x) > 1.0)
        v = "inconclusive" if n_b < F3_MIN_REACTIONS else ("supports" if rank_same >= math.ceil(0.75 * n_b) else ("negative" if big >= 2 else "inconclusive"))
        add("F3", "Skala shifts main to large basis: " + ", ".join(f"{k} {x:+.2f}" for k, x in shifts), n_b, v,
            f"The Skala barrier moves by at most {max(abs(x) for _, x in shifts):.2f} kcal/mol from the triple-zeta to the quadruple-zeta basis in {n_b} reactions.")
    else:
        add("F3", "no basis-check rows (G2 not reached)", 0, "NOT RUN (data missing)", "")
    # F4
    pen, n_react = {}, set()
    for r in table3_scout_cost(D, record=False)[0]:
        if r["geometry_penalty_xtb"] is not None and r["method"] in fn:
            pen.setdefault(r["method"], {})[r["key"]] = r["geometry_penalty_xtb"]
            n_react.add(r["key"])
    if pen:
        thr = math.ceil(0.75 * len(n_react))
        small_all = all(sum(1 for x in v.values() if abs(x) < 1.0) >= thr for v in pen.values())
        allp = {k: x for v in pen.values() for k, x in v.items()}
        n_big = len({k for v in pen.values() for k, x in v.items() if abs(x) > 3.0})
        sn2_all = bool(sn2) and all(any(abs(pen[m].get(k, 0)) > 3.0 for m in pen) for k in sn2)
        v = "inconclusive" if len(n_react) < 4 else ("supports" if small_all else ("negative" if (n_big >= 3 or sn2_all) else "inconclusive"))
        add("F4", "geometry penalty (kcal/mol): " + "; ".join(f"{m}: " + ", ".join(f"{k} {x:+.2f}" for k, x in sorted(vv.items())) for m, vv in pen.items()), len(n_react), v,
            f"Barriers at GFN2-xTB geometries differ from those at B3LYP geometries by at most {max(abs(x) for x in allp.values()):.2f} kcal/mol over {len(n_react)} reactions and {len(pen)} functionals.")
    else:
        add("F4", "no xTB-geometry energies (G3 not reached)", 0, "NOT RUN (data missing)", "")
    # F5
    ranks = []
    for s in REACTIONS:
        l = lit_primary(s)
        if l is None:
            continue
        err = {}
        for m in fn:
            c = compare_to_lit(D, l, m, basis_for(s, "main"))
            if c is not None and c.get("value") is not None:
                err[m] = c["diff"]
        if len(err) == len(fn):
            ranks.append(_rank_with_ties(err)[M_SKALA])
    if ranks:
        from collections import Counter
        modal, cnt = Counter(ranks).most_common(1)[0]
        v = "inconclusive" if len(ranks) < 4 else ("supports" if cnt >= math.ceil(0.75 * len(ranks)) else ("negative" if cnt < 0.5 * len(ranks) else "inconclusive"))
        add("F5", f"Skala rank per reaction (1 = best, ties within {TIE}): {ranks}; most frequent rank {modal} in {cnt} of {len(ranks)}", len(ranks), v,
            f"Skala's rank among the five functionals is {modal} in {cnt} of {len(ranks)} reactions.")
    else:
        add("F5", "needs all five functionals at the stationary points", 0, "NOT RUN (data missing)", "")
    # F6
    fr = _t1_fraction(D)
    if fr:
        good, tot = sum(a for a, _ in fr.values()), sum(b for _, b in fr.values())
        frac = good / tot if tot else float("nan")
        v = "supports" if frac >= 0.9 else ("negative" if frac < 0.5 else "inconclusive")
        sk = _path_errors(D, "ch4", M_SKALA, M_CCSDT_BC) or _path_errors(D, "ch4", M_SKALA, M_CCSDT_FC)
        loc = "" if sk is None else f" The largest Skala error along the CH4 stretch is {np.max(np.abs(sk['err'])):.2f} kcal/mol at C-H = {sk['coord'][int(np.argmax(np.abs(sk['err'])))]:.2f} A."
        add("F6", "reference reliable (T1 < 0.02) at: " + ", ".join(f"{k} {a}/{b}" for k, (a, b) in fr.items()), tot, v,
            f"CCSD(T) is reliable (T1 below 0.02) at {good} of {tot} points of the CH4 stretch" + (" and the HCN map." if "hcn_grid" in fr else ".") + loc)
    else:
        add("F6", "no CH4 CCSD(T) rows", 0, "NOT RUN (data missing)", "")
    return rows


try:                                   # the figures module (GPU notebook only): use the design's criteria and the design's findings function
    findings_outcome = findings_outcome_1h
except Exception:
    pass
