"""Build the JSON files in data/json/ from the CSV files in data/csv/ (and data/records/). Used by build_site.py.

Every number in a JSON file is read from a CSV file; nothing is typed in here. The JSON files are what the interactive pages load.
Needs numpy for the parabola fits and the 2-D views of the path geometries.
"""
import csv, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SITE = os.path.dirname(HERE)
D = os.path.join(SITE, "data")

RX_ORDER = ["hcn", "nh3", "h2co", "sn2_ff", "sn2_clcl", "ch3nc"]
RX_LABEL = {"hcn": "HCN → HNC", "nh3": "NH₃ inversion", "h2co": "H₂CO → H₂ + CO", "ch3nc": "CH₃NC → CH₃CN",
            "sn2_ff": "F⁻ + CH₃F", "sn2_clcl": "Cl⁻ + CH₃Cl"}
M_ORDER = ["Skala-1.1", "Skala-1.1-noD3", "wB97M-V", "B3LYP-D3BJ", "r2SCAN", "PBE", "CCSD(T)-fc", "CCSD(T)+dMP2"]
M_LABEL = {"Skala-1.1": "Skala-1.1", "Skala-1.1-noD3": "Skala-1.1 without D3", "wB97M-V": "ωB97M-V", "B3LYP-D3BJ": "B3LYP-D3(BJ)", "r2SCAN": "r²SCAN",
           "PBE": "PBE", "CCSD(T)-fc": "CCSD(T), frozen core", "CCSD(T)+dMP2": "CCSD(T)+ΔMP2", "GFN2-xTB": "GFN2-xTB"}
SRC_SHORT = [("NHTBH38", "NHTBH38/08 database page"), ("Lin et al", "Lin et al. 2002"), ("Feller", "Feller et al. 2000"), ("Nguyen and Stanton", "Nguyen and Stanton 2018 (abstract)")]
PATH_VIEW = {"hcn": (0, 1), "nh3": (0, 2), "h2co": (0, 1)}


def rd(rel):
    with open(os.path.join(D, "csv", rel), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fl(x):
    try:
        return float(x) if x not in ("", None) else None
    except ValueError:
        return None


def r4(x):
    return None if x is None else round(x, 4)


def short_source(s):
    for k, v in SRC_SHORT:
        if s.startswith(k):
            return v
    return s


def build_results():
    me = rd("study/master_errors.csv")
    intr = {(r["reaction"], r["method"]): r for r in rd("study/derived_internal_reference.csv")}
    rx = []
    for k in RX_ORDER:
        rows = [r for r in me if r["reaction"] == k and r["is_primary_barrier"] == "True" and r["basis_kind"] == "main"]
        sk = [r for r in rows if r["method"] == "Skala-1.1"][0]
        ir = intr[(k, "Skala-1.1")]
        rx.append({"key": k, "label": RX_LABEL[k], "cls": sk["cls"], "tier": sk["tier"], "quantity": sk["quantity"], "lit_id": sk["lit_id"],
                   "ref": fl(sk["reference_value_kcal"]), "ref_unc": fl(sk["reference_unc"]), "ref_source": short_source(sk["reference_source"]), "ref_source_full": sk["reference_source"],
                   "zp": sk["reference_zero_point"], "verified": sk["reference_verified"] == "yes", "headline": sk["in_headline"] == "True",
                   "zpe_added": fl(sk["zpe_added_kcal"]), "own_ref": fl(ir["own_ccsdt_dmp2_kcal"]), "own_usable": ir["own_reference_usable"] == "True", "note": sk["note"],
                   "basis": sk["basis"]})
    rows = []
    for r in me:
        if r["is_primary_barrier"] != "True":
            continue
        k, m, b = r["reaction"], r["method"], r["basis_kind"]
        own = None
        if b == "main" and (k, m) in intr and m in M_ORDER:
            own = fl(intr[(k, m)]["signed_error_vs_own_ccsdt_kcal"])
        rows.append({"r": k, "m": m, "b": b, "basis": r["basis"], "comp": r4(fl(r["computed_kcal"])), "elec": r4(fl(r["electronic_kcal"])), "zpe": r4(fl(r["zpe_added_kcal"])),
                     "err": r4(fl(r["signed_error_kcal"])), "own": r4(own), "wd": r["method_withdrawn_by_check"] == "True", "src": r["source_files"]})
    ss = rd("study/summary_stats.csv")
    pub = [{"m": r["method"], "b": r["basis_kind"], "n": int(r["n_rows"]), "mae": fl(r["MAE"]), "mse": fl(r["MSE"]), "max": fl(r["MaxAE"]), "lo": fl(r["MAE_ci95_low"]), "hi": fl(r["MAE_ci95_high"])}
           for r in ss if r["set"] == "HEADLINE"]
    loo = [{"removed": r["removed"], "m": r["method"], "n": int(r["n"]), "mae": fl(r["MAE"]), "mse": fl(r["MSE"])} for r in rd("study/derived_leave_one_out.csv")]
    t3 = rd("tables/table3_statistics.csv")
    inv = {v: k for k, v in M_LABEL.items()}
    inv["CCSD(T) + MP2 corr."] = "CCSD(T)+dMP2"
    internal = [{"m": inv[r["method"]], "n": int(r["n"]), "mae": fl(r["MAE"]), "mse": fl(r["MSE"]), "max": fl(r["MaxAE"])} for r in t3 if r["block"].startswith("(ii)")]
    return {"reactions": rx, "methods": [{"key": m, "label": M_LABEL[m]} for m in M_ORDER], "rows": rows, "published": pub, "leave_one_out": loo, "published_internal": internal,
            "sources": {"errors": "data/csv/study/master_errors.csv", "internal": "data/csv/study/derived_internal_reference.csv", "stats": "data/csv/study/summary_stats.csv",
                        "loo": "data/csv/study/derived_leave_one_out.csv"}}


def kabsch(P, Q):
    """rotation that best maps P onto Q (both centred, N x 3)"""
    H = P.T @ Q
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    Dm = np.diag([1, 1, d])
    return Vt.T @ Dm @ U.T


def build_paths():
    pp = rd("study/path_profiles.csv")
    frames = json.load(open(os.path.join(D, "records", "path_frames.json"), encoding="utf-8"))
    out = {}
    for k in ("hcn", "nh3", "h2co"):
        pts = frames[k]["points"]
        X = [np.array([a[1] for a in p["atoms"]], float) for p in pts]
        X = [x - x.mean(axis=0) for x in X]
        ref = X[len(X) // 2]
        aligned = []
        for x in X:
            R = kabsch(x, ref)
            aligned.append(x @ R.T)
        allx = np.vstack(aligned)
        # principal axes of the union of all aligned frames
        w, v = np.linalg.eigh(np.cov(allx.T))
        axes = v[:, np.argsort(w)[::-1]]
        ia, ib = PATH_VIEW[k]
        ic = ({0, 1, 2} - {ia, ib}).pop()
        proj = [a @ axes for a in aligned]
        atoms = [p["atoms"][i][0] for p in pts for i in range(1)] and [a[0] for a in pts[0]["atoms"]]
        pdata = []
        for p, pr in zip(pts, proj):
            pdata.append({"id": "p%02d" % p["index"], "coord": r4(p["coordinate"]),
                          "xyz": [[round(float(pr[i, ia]), 3), round(float(pr[i, ib]), 3), round(float(pr[i, ic]), 3)] for i in range(len(atoms))]})
        series = {}
        for r in pp:
            if r["reaction"] != k:
                continue
            ref_name = r["reference"]
            s = series.setdefault(ref_name, {"ref": {}, "err": {}, "used": {}})
            s["ref"][r["point_id"]] = fl(r["ref_relative_kcal"])
            s["err"].setdefault(r["method"], {})[r["point_id"]] = fl(r["signed_error_kcal"])
            s["used"][r["point_id"]] = r["used_in_reference_stats"] == "True"
        out[k] = {"label": RX_LABEL[k], "coord_label": frames[k]["label"], "basis": frames[k]["basis"], "elements": atoms, "points": pdata, "series": series,
                  "sources": "data/csv/study/path_profiles.csv; geometries: data/geometries/%s_path_own_scan.xyz" % k}
    return out


def build_sn2():
    sc = rd("study/sn2_symscan.csv")
    mins = rd("study/derived_sn2_minima.csv")
    out = {}
    for rxk in ("sn2_ff", "sn2_clcl"):
        entry = {"label": RX_LABEL[rxk], "methods": {}}
        for m in ("Skala-1.1", "wB97M-V"):
            pts = []
            for r in sc:
                if r["reaction"] != rxk or r["method"] != m or r["kind"] not in ("symmetric_stretch_scan", "own_B3LYP_TS"):
                    continue
                if r["density_fitting"] != "True":
                    continue
                pts.append({"x": fl(r["cx_distance_A"]), "y": r4(fl(r["rel_to_RC_B3LYPgeom_kcal"])), "y_noD3": r4(fl(r["rel_to_RC_noD3_kcal"])), "conv": r["converged"] == "1",
                            "id": r["point_id"], "attempt": r["attempt"], "ts": r["kind"] == "own_B3LYP_TS"})
            pts.sort(key=lambda p: (p["x"], not p["conv"]))
            mrow = [x for x in mins if x["reaction"] == rxk and x["method"] == m][0]
            fit_x = [float(s) for s in mrow["fit_points_A"].split(";")]
            fit_y = []
            for fx in fit_x:
                cand = [p for p in pts if abs(p["x"] - fx) < 1e-3 and p["conv"]]
                fit_y.append(cand[0]["y"])
            c = np.polyfit(fit_x, fit_y, 2)
            vx = -c[1] / (2 * c[0]); vy = np.polyval(c, vx)
            assert abs(vx - float(mrow["vertex_cx_A"])) < 2e-3 and abs(vy - float(mrow["vertex_energy_kcal"])) < 2e-3, (rxk, m, vx, vy, mrow)
            entry["methods"][m] = {"points": pts, "fit_x": fit_x, "coef": [float(x) for x in c],
                                   "vertex": {"x": fl(mrow["vertex_cx_A"]), "y": fl(mrow["vertex_energy_kcal"])}, "energy_at_ts": fl(mrow["energy_at_b3lyp_ts_kcal"]),
                                   "err_vs_ccsdt": fl(mrow["error_vertex_vs_ccsdt"]), "err_vs_lit": fl(mrow["error_vertex_vs_lit"])}
            entry["ts_x"] = fl(mrow["b3lyp_ts_cx_A"]); entry["ref_ccsdt"] = fl(mrow["ref_ccsdt_dmp2_kcal"]); entry["ref_lit"] = fl(mrow["literature_kcal"])
        entry["sources"] = "data/csv/study/sn2_symscan.csv, data/csv/study/derived_sn2_minima.csv"
        out[rxk] = entry
    return out


def build_t1():
    rows = rd("figures/figS5_hcn_t1_map_data.csv")
    th = sorted(set(float(r["theta"]) for r in rows)); rr = sorted(set(float(r["r"]) for r in rows))
    t1 = [[None] * len(rr) for _ in th]; er = [[None] * len(rr) for _ in th]; win = [[0] * len(rr) for _ in th]
    for r in rows:
        i = th.index(float(r["theta"])); j = rr.index(float(r["r"]))
        t1[i][j] = round(float(r["t1"]), 5); er[i][j] = round(float(r["E_rel_kcal"]), 2); win[i][j] = 1 if r["in_window"] == "True" else 0
    return {"theta": th, "r": [round(x, 4) for x in rr], "t1": t1, "erel": er, "win": win, "window_kcal": 120,
            "sources": "data/csv/figures/figS5_hcn_t1_map_data.csv (values of data/csv/results/phase05_hcn_grid.csv)"}


def build_cost():
    pairs = rd("study/derived_timing_pairs.csv")
    keep = ["B3LYP-D3BJ", "PBE", "Skala-1.1", "r2SCAN", "wB97M-V"]
    rows = [{"s": r["system"], "p": r["point_id"], "basis": r["basis"], "nbf": int(float(r["nbf"])), **{m: fl(r[m]) for m in keep}} for r in pairs]
    ss = {r["method"]: r for r in rd("study/summary_stats.csv") if r["set"] == "HEADLINE" and r["basis_kind"] == "main"}
    mae = {m: {"mae": fl(ss[m]["MAE"]), "n": int(ss[m]["n_rows"])} for m in keep}
    f12 = {r["method"]: r for r in rd("figures/fig12_cost_accuracy_data.csv") if r["panel"] == "b"}
    pub = {m: {"median": fl(f12[m]["median_seconds"]), "loo_min": fl(f12[m]["loo_min"]), "loo_max": fl(f12[m]["loo_max"]), "n_common": fl(f12[m]["n_common"])} for m in keep}
    return {"methods": keep, "rows": rows, "mae": mae, "published": pub,
            "sources": "data/csv/study/derived_timing_pairs.csv, data/csv/study/summary_stats.csv, data/csv/figures/fig12_cost_accuracy_data.csv"}


def build_basis():
    mb = rd("study/master_barriers.csv")
    keep = ("Skala-1.1", "Skala-1.1-noD3", "wB97M-V")
    basis = []
    d3 = []
    for r in mb:
        if r["method"] in keep and r["basis_large"]:
            basis.append({"r": r["reaction"], "m": r["method"], "main": r4(fl(r["barrier_main_used_kcal"])), "large": r4(fl(r["barrier_large_used_kcal"])), "shift": r4(fl(r["shift_main_to_large_kcal"])),
                          "err_main": r4(fl(r["signed_error_main_kcal"])), "err_large": r4(fl(r["signed_error_large_kcal"])), "main_basis": r["basis_main"], "large_basis": r["basis_large"]})
    by = {(r["reaction"], r["method"]): r for r in mb}
    for k in RX_ORDER:
        a = by[(k, "Skala-1.1")]; b = by[(k, "Skala-1.1-noD3")]
        d3.append({"r": k, "with": r4(fl(a["barrier_main_electronic_kcal"])), "without": r4(fl(b["barrier_main_electronic_kcal"])),
                   "d3": r4(fl(a["barrier_main_electronic_kcal"]) - fl(b["barrier_main_electronic_kcal"])), "err_with": r4(fl(a["signed_error_main_kcal"])), "err_without": r4(fl(b["signed_error_main_kcal"])),
                   "headline": a["in_headline"] == "True"})
    sp = [{"system": r["system"], "species": r["species"], "m": r["method"], "lowering": fl(r["lowering_TZVP_to_QZVPP_kcal"])} for r in rd("study/derived_basis_species.csv")]
    return {"basis": basis, "d3": d3, "species": sp,
            "sources": "data/csv/study/master_barriers.csv, data/csv/study/derived_basis_species.csv"}


def build_all():
    os.makedirs(os.path.join(D, "json"), exist_ok=True)
    res = {"results": build_results(), "paths": build_paths(), "sn2": build_sn2(), "t1grid": build_t1(), "cost": build_cost(), "basis": build_basis()}
    for k, v in res.items():
        with open(os.path.join(D, "json", k + ".json"), "w", encoding="utf-8") as f:
            json.dump(v, f, ensure_ascii=False, separators=(",", ":"))
    return res
