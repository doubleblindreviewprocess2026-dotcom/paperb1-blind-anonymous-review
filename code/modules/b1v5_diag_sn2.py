# ==============================================================================
# b1v5_diag_sn2.py -- diagnostic of the Cl- + CH3Cl transition state with Skala-1.1 (2026-11-16).
# Question: is the Skala error of -3.1 kcal/mol (own B3LYP TS geometry, density-fitted rescue) real, or an SCF / density-fitting / geometry artefact?
# In the full run Skala at the xTB TS (C-Cl 2.20 A) gave a barrier of 16.0, at the B3LYP TS (C-Cl 2.35 A) only 10.3, while wB97M-V changed by 0.2.
# Pasted after b1v5_gpu.py, b1v5_1h.py, b1v5_1h_gpu.py.  Rows go to diag_sn2_clcl.csv (not read by the check).  Nothing else is edited.
#   D1 Skala (with D3, density fitting) at RC and TS of the own B3LYP geometries, and at the xTB RC and TS
#   D2 the same TS without density fitting, started from (a) default guess, (b) the density-fitted density, (c) the PBE density
#   D3 rigid stretch of both C-Cl bonds around the own TS (CH3 fixed): Skala and wB97M-V energy profiles, the Skala maximum (parabola vertex)
#   D4 fragments Cl- and CH3Cl (density fitting) for the separated-reactants barrier
# ==============================================================================
import os, json, time
import numpy as np

DIAG_CSV = "diag_sn2_clcl.csv"
DIAG_SYS = "sn2_clcl"
DIAG_LEVEL = "diag"
DIAG_CAP_MIN = 30
STRETCH_SKALA = [-0.30, -0.20, -0.15, -0.10, -0.05, 0.10]          # offsets of both C-Cl distances from the own TS (2.348 A), in angstrom; 0.00 is the TS itself (D1)
STRETCH_WB97MV = [-0.20, -0.10, 0.10]
DIAG_DM = {}
GPU_CAPS["phase17"] = DIAG_CAP_MIN
GPU_CAPS["phase18"] = 40


def _atoms_of(sp, level=None):
    if level:
        g = geometry_of(DIAG_SYS, sp, level)
    else:
        _, g = best_geom(DIAG_SYS, sp)
    return None if g is None else g["atoms"]


def _stretched(atoms, d, halogens=("Cl",)):
    """Move both halogen atoms along their C-X bond by d angstrom (CH3 untouched)."""
    ic = [i for i, a in enumerate(atoms) if a[0] == "C"][0]
    c = np.array(atoms[ic][1], float)
    out = []
    for a in atoms:
        if a[0] in halogens:
            v = np.array(a[1], float) - c
            n = np.linalg.norm(v)
            out.append([a[0], [float(x) for x in c + v * (1.0 + d / n)]])
        else:
            out.append([a[0], [float(x) for x in a[1]]])
    return out


def _cl_c(atoms):
    ic = [i for i, a in enumerate(atoms) if a[0] == "C"][0]
    c = np.array(atoms[ic][1], float)
    return sorted(round(float(np.linalg.norm(np.array(a[1], float) - c)), 4) for a in atoms if a[0] == "Cl")


def _diag_scf(method, atoms, charge, basis, df, start_key=None, store_key=None, plan=None):
    """One energy with an explicit ladder.  plan = list of (name, kwargs); default ladder: warm start from DIAG_DM[start_key], default guess,
    PBE start + level shift 0.5.  Returns dict(energy, converged, extra) like energy_of."""
    route = route_of(method)
    nao = make_mol(atoms, charge, basis).nao
    if plan is None:
        plan = []
        dm = DIAG_DM.get(start_key) if start_key else None
        if dm is not None and dm.shape == (nao, nao):
            plan.append(("warm start", dict(dm0=dm)))
        plan += [("default", {}), ("PBE start + level shift 0.5", dict(rescue=True, level_shift=0.5, max_cycle=400))]
    best, tried = None, []
    for name, kw in plan:
        kw = dict(kw)
        if kw.pop("rescue", False):
            pr = _scf_once(M_PBE, atoms, charge, basis, route_of(M_PBE), df=df)
            if not pr["converged"] or pr["dm"] is None:
                tried.append(name + " skipped (PBE start failed)"); continue
            kw["dm0"] = pr["dm"]
        check_time()
        try:
            res = _scf_once(method, atoms, charge, basis, route, df=df, **kw)
        except PhaseTimeout:
            raise
        except Exception as ex:
            tried.append(name + " crashed " + repr(ex)[:60]); gpu_cleanup(); continue
        tried.append(name + (" ok" if res["converged"] else " not converged"))
        if res["converged"]:
            best = dict(res, attempt=name); break
        if best is None:
            best = dict(res, attempt=name + " (not converged)")
    gpu_cleanup()
    if best is None:
        raise RuntimeError("no attempt possible: " + "; ".join(tried))
    if store_key and best["converged"] and best["dm"] is not None:
        DIAG_DM[store_key] = best["dm"]
    ex = dict(df=bool(df), attempt=best["attempt"], tried=tried, nbf=best["nbf"])
    if best.get("d3") is not None:
        ex["d3_hartree"] = float(best["d3"])
    return dict(energy=float(best["energy"]), converged=bool(best["converged"]), extra=ex)


def _drow(pid, method, atoms, charge, basis, df, kind, start_key=None, store_key=None, plan=None, coordinate=""):
    return qrun(DIAG_CSV, DIAG_SYS, pid, DIAG_LEVEL, method, basis,
                lambda: _diag_scf(method, atoms, charge, basis, df, start_key, store_key, plan), kind, coordinate=coordinate, retry_failed=True,
                label=f"{pid}:{method}")


def run_sn2_diagnostics():
    basis = basis_for(DIAG_SYS, "main")
    ts_own, rc_own = _atoms_of("TS", LV_OWN), _atoms_of("RC", LV_OWN)
    ts_xtb, rc_xtb = (geometry_of(DIAG_SYS, "TS", LV_XTB) or {}).get("atoms"), (geometry_of(DIAG_SYS, "RC", LV_XTB) or {}).get("atoms")
    if ts_own is None or rc_own is None:
        raise RuntimeError("the own-geometry TS or RC of Cl- + CH3Cl is missing (geometry pack)")
    chg_ts, chg_rc = species_charge(DIAG_SYS, "TS"), species_charge(DIAG_SYS, "RC")
    print(f"  C-Cl distances at the own TS: {_cl_c(ts_own)}, at the xTB TS: {_cl_c(ts_xtb) if ts_xtb else 'missing'}")
    S, W = M_SKALA, M_WB97MV
    kS, kW = "gpu:tz:sn2_clcl:" + S, "gpu:tz:sn2_clcl:" + W
    # D1: Skala, density fitting, own geometries
    with item_guard("phase17", "D1 RC"):
        _drow("RC_own_DF", S, rc_own, chg_rc, basis, True, kS, store_key="rc_own")
    with item_guard("phase17", "D1 TS"):
        _drow("TS_own_DF", S, ts_own, chg_ts, basis, True, kS, store_key="ts_own_df")
    # D2: TS without density fitting
    with item_guard("phase17", "D2 TS non-DF default"):
        _drow("TS_own_noDF_default", S, ts_own, chg_ts, basis, False, kS, plan=[("default, 60 cycles", dict(max_cycle=60))])
    if DIAG_DM.get("ts_own_df") is not None:
        with item_guard("phase17", "D2 TS non-DF from DF density"):
            _drow("TS_own_noDF_fromDF", S, ts_own, chg_ts, basis, False, kS, plan=[("start from the density-fitted density", dict(dm0=DIAG_DM["ts_own_df"], max_cycle=100))])
    with item_guard("phase17", "D2 TS non-DF PBE start"):
        _drow("TS_own_noDF_PBEstart", S, ts_own, chg_ts, basis, False, kS, plan=[("PBE start + level shift 0.5", dict(rescue=True, level_shift=0.5, max_cycle=300))])
    # D1 continued: xTB geometries
    if ts_xtb and rc_xtb:
        with item_guard("phase17", "D1 xTB RC"):
            _drow("RC_xtb_DF", S, rc_xtb, chg_rc, basis, True, kS)
        with item_guard("phase17", "D1 xTB TS"):
            _drow("TS_xtb_DF", S, ts_xtb, chg_ts, basis, True, kS)
    # D3: stretch profile, Skala (warm start from the neighbouring point, from TS outward and inward)
    order = sorted([d for d in STRETCH_SKALA if d < 0], reverse=True) + sorted([d for d in STRETCH_SKALA if d > 0])
    for d in order:
        with item_guard("phase17", f"D3 Skala {d:+.2f}"):
            at = _stretched(ts_own, d)
            _drow(f"TS_stretch{d:+.2f}_DF", S, at, chg_ts, basis, True, kS, start_key="ts_own_df" if abs(d) > 0 else None, coordinate=f"{_cl_c(at)[0]:.3f}")
    # D3: wB97M-V profile and RC, density fitting
    with item_guard("phase17", "D3 wB97M-V RC"):
        _drow("RC_own_DF", W, rc_own, chg_rc, basis, True, kW)
    with item_guard("phase17", "D3 wB97M-V TS"):
        _drow("TS_own_DF", W, ts_own, chg_ts, basis, True, kW)
    for d in STRETCH_WB97MV:
        with item_guard("phase17", f"D3 wB97M-V {d:+.2f}"):
            at = _stretched(ts_own, d)
            _drow(f"TS_stretch{d:+.2f}_DF", W, at, chg_ts, basis, True, kW, coordinate=f"{_cl_c(at)[0]:.3f}")
    # D4: fragments (separated reactants), Skala density fitting
    for sp in ("Cl-", "CH3Cl"):
        with item_guard("phase17", f"D4 {sp}"):
            at = _atoms_of(sp, LV_OWN)
            if at:
                _drow(f"{sp}_own_DF", S, at, species_charge(DIAG_SYS, sp), basis, True, kS)


def _drows():
    out = {}
    for r in read_rows(DIAG_CSV):
        if r["converged"] in ("1", "True") and r["energy_hartree"] not in ("", "nan"):
            out[(r["point_id"], r["method"])] = (float(r["energy_hartree"]), r)
    return out


def report_sn2_diagnostics():
    """Prints the findings and writes DIAG_SN2_SUMMARY.md and diag_sn2_stretch.png."""
    H = HARTREE_KCAL
    R = _drows()
    S, W = M_SKALA, M_WB97MV
    allrows = read_rows(DIAG_CSV)
    L = ["# Diagnostic: Cl- + CH3Cl with Skala-1.1", "", f"Written {now()}.  Basis {basis_for(DIAG_SYS, 'main')}, density fitting where marked DF.  Energies in kcal/mol relative to the stated reference.", ""]
    L += ["## Every calculation", "", "| point | method | converged | energy (hartree) | attempt | seconds |", "|---|---|---|---|---|---|"]
    for r in allrows:
        x = extra_of(r)
        L.append(f"| {r['point_id']} | {r['method']} | {r['converged']} | {r['energy_hartree']} | {x.get('attempt', '')} | {r['seconds']} |")
    L.append("")

    def E(pid, m):
        return R[(pid, m)][0] if (pid, m) in R else None

    def barrier(ts, rc, m):
        a, b = E(ts, m), E(rc, m)
        return None if a is None or b is None else (a - b) * H
    L += ["## Barrier from the complex (TS minus RC), density fitting", ""]
    for m in (S, W):
        for lab, ts, rc in (("own B3LYP geometries", "TS_own_DF", "RC_own_DF"), ("xTB geometries", "TS_xtb_DF", "RC_xtb_DF")):
            b = barrier(ts, rc, m)
            L.append(f"- {m}, {lab}: " + (f"{b:.2f}" if b is not None else "not available"))
    L.append("- literature (CCSD(T) + basis correction of this study): 13.41 / 13.61; the full run gave Skala 10.29 (own, rescue) and 15.9 (xTB), wB97M-V 13.50 and 13.71")
    L.append("")
    L += ["## Is the density-fitted TS the same SCF state as without density fitting?", ""]
    e_df = E("TS_own_DF", S)
    for pid in ("TS_own_noDF_default", "TS_own_noDF_fromDF", "TS_own_noDF_PBEstart"):
        e = E(pid, S)
        L.append(f"- {pid}: " + ("did not converge" if e is None else f"{e:.8f} hartree, {'' if e_df is None else f'{(e - e_df) * H:+.2f} kcal/mol relative to the density-fitted TS'}"))
    L.append("")
    ds_s = sorted([(float(k[0].split('stretch')[1].split('_')[0]), v[0]) for k, v in R.items() if k[1] == S and k[0].startswith("TS_stretch")] + ([(0.0, e_df)] if e_df is not None else []))
    ds_w = sorted([(float(k[0].split('stretch')[1].split('_')[0]), v[0]) for k, v in R.items() if k[1] == W and k[0].startswith("TS_stretch")] + ([(0.0, E('TS_own_DF', W))] if E('TS_own_DF', W) is not None else []))
    ts_dist = _cl_c(_atoms_of("TS", LV_OWN))[0]
    L += ["## Stretch profile (both C-Cl bonds moved together, CH3 fixed), energy relative to the RC (own geometry)", "", "| C-Cl distance (A) | offset | Skala-1.1 | wB97M-V |", "|---|---|---|---|"]
    offs = sorted({d for d, _ in ds_s} | {d for d, _ in ds_w})
    rcS, rcW = E("RC_own_DF", S), E("RC_own_DF", W)
    dS, dW = dict(ds_s), dict(ds_w)
    for d in offs:
        a = f"{(dS[d] - rcS) * H:.2f}" if d in dS and rcS is not None else ""
        b = f"{(dW[d] - rcW) * H:.2f}" if d in dW and rcW is not None else ""
        L.append(f"| {ts_dist + d:.3f} | {d:+.2f} | {a} | {b} |")
    L.append("")
    verdict = []
    for lab, ds, rc in (("Skala-1.1", ds_s, rcS), ("wB97M-V", ds_w, rcW)):
        if len(ds) >= 3 and rc is not None:
            i = int(np.argmax([e for _, e in ds]))
            lo, hi = max(0, i - 1), min(len(ds), i + 2)
            seg = ds[lo:hi]
            if len(seg) == 3:
                x = np.array([d for d, _ in seg]); y = np.array([(e - rc) * H for _, e in seg])
                c = np.polyfit(x, y, 2)
                xv = -c[1] / (2 * c[0]) if c[0] < 0 else x[1]
                yv = float(np.polyval(c, xv))
                verdict.append(f"{lab}: highest point at offset {ds[i][0]:+.2f} A (C-Cl {ts_dist + ds[i][0]:.3f} A), {(ds[i][1] - rc) * H:.2f} kcal/mol; parabola vertex offset {xv:+.2f} A, {yv:.2f} kcal/mol"
                               + ("" if 0 < i < len(ds) - 1 else "  ** the maximum is at the edge of the scanned range: the profile has no top inside it **"))
            else:
                verdict.append(f"{lab}: highest point at offset {ds[i][0]:+.2f} A, {(ds[i][1] - rc) * H:.2f} kcal/mol (fewer than three points around it)")
    L += ["## Maximum along the stretch", ""] + [f"- {v}" for v in verdict] + [""]
    sf = barrier("TS_own_DF", "RC_own_DF", S)
    L += ["## Reading (automatic, check it)", "",
          "- If the Skala maximum lies at a shorter C-Cl distance than the B3LYP TS (2.35 A) and its height is near 16 kcal/mol, the -3.1 error of the full run came from evaluating Skala at a geometry that is not its transition state; "
          "the Skala barrier of this reaction is then the maximum of the profile (about +2.5 against the reference), and a Skala-relaxed TS is the proper next step.",
          "- If the profile is irregular (jumps between neighbouring points) or the non-DF TS energy differs from the DF one by more than 0.5 kcal/mol, the SCF solution is not reliable and the Skala value for this reaction should be reported as not converged.",
          "- If the profile is smooth with its top near 2.35 A and about 10 kcal/mol, the -3.1 is a property of Skala for this reaction."]
    open(path("DIAG_SN2_SUMMARY.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5.5, 3.6))
        for lab, ds, rc, col in (("Skala-1.1", ds_s, rcS, "#c0392b"), ("wB97M-V", ds_w, rcW, "#2c3e50")):
            if ds and rc is not None:
                ax.plot([ts_dist + d for d, _ in ds], [(e - rc) * H for _, e in ds], "o-", color=col, label=lab)
        ax.axhline(13.41, color="grey", ls="--", lw=1, label="CCSD(T) barrier 13.4"); ax.axvline(ts_dist, color="grey", ls=":", lw=1)
        ax.set_xlabel("C-Cl distance (angstrom), both bonds"); ax.set_ylabel("energy relative to RC (kcal/mol)"); ax.legend(fontsize=8); fig.tight_layout()
        os.makedirs(path("figures"), exist_ok=True); fig.savefig(path("figures", "diag_sn2_stretch.png"), dpi=150); plt.close(fig)
        print("wrote figures/diag_sn2_stretch.png")
    except Exception as ex:
        print("figure not written:", repr(ex)[:100])


# ==============================================================================
# Second diagnostic (2026-11-16): the symmetric stretch of a D3h SN2 transition state is NOT the reaction coordinate: the TS is a MINIMUM along it
# (the reaction coordinate is the antisymmetric stretch).  For a method, the TS energy is therefore the MINIMUM of E(d) along the symmetric stretch
# (rigid CH3), and the method's TS has the C-X distance of that minimum.  Scan until the minimum is inside the range.
# ==============================================================================
SYM_REFS = {"sn2_clcl": dict(name="Cl- + CH3Cl", cc=13.61, lit=13.41), "sn2_ff": dict(name="F- + CH3F", cc=13.60, lit=13.38)}
SYM_PLAN = {   # system -> {method: offsets (A) from the own B3LYP TS}
    "sn2_clcl": {M_SKALA: [0.0, 0.10, 0.20, 0.30, 0.40, -0.05, -0.10], M_WB97MV: [0.0, -0.05, 0.05, -0.10, 0.10]},
    "sn2_ff":   {M_SKALA: [0.0, 0.10, 0.20, -0.10, -0.05, 0.05, 0.15], M_WB97MV: [0.0, -0.10, -0.05, 0.05, 0.10, -0.15]},
}


def run_symscan():
    for system, plan in SYM_PLAN.items():
        basis = basis_for(system, "main")
        ts = (geometry_of(system, "TS", LV_OWN) or {}).get("atoms")
        rc = (geometry_of(system, "RC", LV_OWN) or {}).get("atoms")
        if ts is None or rc is None:
            print(f"  {system}: own TS or RC missing, skipped"); continue
        chg_ts, chg_rc = species_charge(system, "TS"), species_charge(system, "RC")
        print(f"  {system}: X-C distances at the own TS: {_cl_c_any(ts)}")
        for method, offsets in plan.items():
            kind = f"gpu:tz:{system}:{method}"
            with item_guard("phase18", f"{system} RC {method}"):
                qrun(DIAG_CSV, system, "RC_own_DF", DIAG_LEVEL, method, basis,
                     lambda m=method: _diag_scf(m, rc, chg_rc, basis, True, None, None, None), kind, retry_failed=True, label=f"{system}:RC:{method}")
            key = f"{system}_{method}_ts"
            for d in offsets:
                with item_guard("phase18", f"{system} {method} {d:+.2f}"):
                    at = ts if d == 0 else _stretched(ts, d, ("Cl", "F"))
                    pid = "TS_own_DF" if d == 0 else f"TS_stretch{d:+.2f}_DF"
                    qrun(DIAG_CSV, system, pid, DIAG_LEVEL, method, basis,
                         lambda m=method, at=at, key=key, d=d: _diag_scf(m, at, chg_ts, basis, True, key if d != 0 else None, key if d == 0 else None, None), kind,
                         coordinate=f"{_cl_c_any(at)[0]:.3f}", retry_failed=True, label=f"{system}:{pid}:{method}")


def _cl_c_any(atoms):
    ic = [i for i, a in enumerate(atoms) if a[0] == "C"][0]
    c = np.array(atoms[ic][1], float)
    return sorted(round(float(np.linalg.norm(np.array(a[1], float) - c)), 4) for a in atoms if a[0] in ("Cl", "F"))


def report_symscan():
    H = HARTREE_KCAL
    rows = [r for r in read_rows(DIAG_CSV) if r["converged"] in ("1", "True") and r["energy_hartree"] not in ("", "nan")]
    L = ["# Symmetric-stretch scan: the transition state is the MINIMUM along this coordinate", "", f"Written {now()}. Rigid CH3 (B3LYP geometry), both X-C bonds moved together, density fitting, def2-TZVPD. Energies relative to the reaction complex (RC) of the same method, kcal/mol.", ""]
    ok_fig = []
    for system, ref in SYM_REFS.items():
        ts_atoms = (geometry_of(system, "TS", LV_OWN) or {}).get("atoms")
        if ts_atoms is None:
            continue
        d0 = _cl_c_any(ts_atoms)[0]
        L += [f"## {ref['name']}", "", f"B3LYP TS: X-C = {d0:.3f} A. Reference barrier from the complex: CCSD(T) + basis correction {ref['cc']:.2f}, literature {ref['lit']:.2f}.", ""]
        table, res = {}, {}
        for method in (M_SKALA, M_WB97MV):
            rc = [float(r["energy_hartree"]) for r in rows if r["system"] == system and r["point_id"] == "RC_own_DF" and r["method"] == method]
            pts = []
            for r in rows:
                if r["system"] == system and r["method"] == method and (r["point_id"] == "TS_own_DF" or r["point_id"].startswith("TS_stretch")):
                    d = 0.0 if r["point_id"] == "TS_own_DF" else float(r["point_id"].split("stretch")[1].split("_")[0])
                    pts.append((d, float(r["energy_hartree"])))
            pts = sorted(set(pts))
            if not rc or len(pts) < 3:
                res[method] = None; continue
            table[method] = {d: (e - rc[0]) * H for d, e in pts}
            ys = [(e - rc[0]) * H for _, e in pts]; xs = [d for d, _ in pts]
            i = int(np.argmin(ys))
            if 0 < i < len(xs) - 1:
                c = np.polyfit(xs[i - 1:i + 2], ys[i - 1:i + 2], 2)
                xv = -c[1] / (2 * c[0]) if c[0] > 0 else xs[i]
                yv = float(np.polyval(c, xv)) if c[0] > 0 else ys[i]
                res[method] = (xv, yv, True, ys[i])
            else:
                res[method] = (xs[i], ys[i], False, ys[i])
        offs = sorted({d for t in table.values() for d in t})
        L += ["| X-C distance (A) | offset | " + " | ".join(table) + " |", "|---|---|" + "---|" * len(table)]
        for d in offs:
            L.append(f"| {d0 + d:.3f} | {d:+.2f} | " + " | ".join((f"{table[m][d]:.2f}" if d in table[m] else "") for m in table) + " |")
        L.append("")
        for m, v in res.items():
            if v is None:
                L.append(f"- {m}: not enough points"); continue
            xv, yv, inner, ymin = v
            err = yv - ref["cc"]
            L.append(f"- **{m}**: minimum along the symmetric stretch at offset {xv:+.2f} A (X-C {d0 + xv:.3f} A), TS barrier **{yv:.2f}** kcal/mol, error against CCSD(T) {err:+.2f}"
                     + ("" if inner else "  ** the minimum is at the edge of the scanned range: the true value is lower still, scan further **"))
        L.append("")
        ok_fig.append((system, ref, d0, table))
    L += ["## Reading", "",
          "- The Skala value at the B3LYP geometry is an upper bound of Skala's TS energy; the minimum along the symmetric stretch is the better estimate (the CH3 internal coordinates are not relaxed, a small effect).",
          "- A minimum at a clearly longer X-C distance than that of wB97M-V means Skala prefers a looser transition state: report the Skala barrier at its own TS, and the barrier at B3LYP geometries as a geometry-sensitivity result.",
          "- Points on the repulsive side (short X-C, such as the xTB TS of the full run) are not barriers of any functional and must not be read as one."]
    open(path("DIAG_SYM_SUMMARY.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, max(1, len(ok_fig)), figsize=(5.2 * max(1, len(ok_fig)), 3.6), squeeze=False)
        for ax, (system, ref, d0, table) in zip(axs[0], ok_fig):
            for m, col in ((M_SKALA, "#c0392b"), (M_WB97MV, "#2c3e50")):
                if m in table:
                    xs = sorted(table[m]); ax.plot([d0 + d for d in xs], [table[m][d] for d in xs], "o-", color=col, label=m)
            ax.axhline(ref["cc"], color="grey", ls="--", lw=1, label=f"CCSD(T) {ref['cc']:.1f}"); ax.axvline(d0, color="grey", ls=":", lw=1)
            ax.set_title(ref["name"], fontsize=9); ax.set_xlabel("X-C distance (A), both bonds"); ax.set_ylabel("energy relative to RC (kcal/mol)"); ax.legend(fontsize=7)
        fig.tight_layout(); os.makedirs(path("figures"), exist_ok=True); fig.savefig(path("figures", "diag_sym_stretch.png"), dpi=150); plt.close(fig)
        print("wrote figures/diag_sym_stretch.png")
    except Exception as ex:
        print("figure not written:", repr(ex)[:100])
