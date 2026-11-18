# ==============================================================================
# b1v5_gpu_figs.py  --  phase G5: figures and tables for Paper B1 v5, regenerated from the CSV files.
# Pasted into the G5 cell of the GPU notebook, after the check module (b1v5_check.py), whose classes
# (Data, quantity, compare_to_lit, ...) it reuses so that errors are computed exactly as the automated check does.
#
# Writes into the Drive folder:
#   figures/fig1_design.png ... fig6_hcn_map_ch4.png      the six main figures of plan section 1
#   numbers_v5.csv                                        every number used in a figure or table, one per line
#   table1_systems.csv, table2_errors.csv, table2_detail.csv, table3_scout_cost.csv   and  summary_tables.xlsx (if openpyxl is present)
#   and the analysis outputs of review section 4: findings_outcome.csv (F1 to F6 against the FROZEN criteria), path_statistics.csv,
#   basis_table.csv, scout_table.csv, geometry_table.csv, skala_d3_table.csv, reference_quality.csv, excluded_and_failed.csv
# A figure whose data is missing is skipped with a message; the other figures are still made.
# ==============================================================================
import os, csv, json, math
import numpy as np
try:
    COLUMNS  # noqa
except NameError:
    from b1v5_common import *  # noqa

NUMS = []                       # rows of numbers_v5.csv
FIG_STATUS = {}                 # figure/table name -> "made" or the reason it was skipped
# one colour per method, fixed (colour-blind-safe set, in the order of the reference palette); markers differ as well
METHOD_STYLE = {
    M_SKALA: ("#2a78d6", "o"), M_SKALA_NOD3: ("#7fb3ee", "h"), M_PBE: ("#eb6834", "v"), M_R2SCAN: ("#1baf7a", "^"), M_B3LYP: ("#e87ba4", "s"),
    M_WB97MV: ("#6250d6", "D"), "wB97X-V": ("#6250d6", "D"), M_CCSDT_FC: ("#0b0b0b", "P"), M_CCSDT_BC: ("#0b0b0b", "P"),
}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e1e0d9"


def _plt():
    import matplotlib
    try:
        matplotlib.use("Agg", force=False)
    except Exception:
        pass
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white", "font.size": 9,
                         "axes.edgecolor": "#898781", "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
                         "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                         "legend.frameon": False, "figure.dpi": 100, "savefig.dpi": 200, "axes.axisbelow": True})
    return plt


def _methods_present(D, basis_of=None):
    """Ladder methods that have any row (anchor shown under its real name)."""
    have_ = {r["method"] for r in D.rows}
    out = [m for m in (M_SKALA, M_SKALA_NOD3, M_WB97MV, M_WB97XV, M_B3LYP, M_R2SCAN, M_PBE) if m in have_]
    return out


def _save(fig, name, standin=False):
    os.makedirs(path("figures"), exist_ok=True)
    if standin:
        fig.text(0.5, 0.5, "STAND-IN DRY RUN\nnot Skala, not results", ha="center", va="center", fontsize=26, color="#c3c2b7", alpha=0.45, rotation=20)
    p = path("figures", name)
    fig.savefig(p, bbox_inches="tight")
    try:
        from IPython.display import display, Image
        display(Image(filename=p, width=760))
    except Exception:
        pass
    _plt().close(fig)
    print("  wrote", p)
    return p


def _nofig(plt, msg):
    fig, ax = plt.subplots(figsize=(5, 2))
    ax.text(0.5, 0.5, msg, ha="center", va="center", color=INK2); ax.set_axis_off()
    return fig


def num(table, system, quantity, method, basis, level, value, ref=None, note=""):
    NUMS.append(dict(table=table, system=system, quantity=quantity, method=method, basis=basis, geometry_level=level,
                     value="" if value is None else round(float(value), 4), reference="" if ref is None else round(float(ref), 4),
                     error="" if (value is None or ref is None) else round(float(value) - float(ref), 4), unit="kcal/mol", note=note))


def primary_q(system):
    return "barrier_fwd_cplx" if SYSTEMS[system]["sn2"] else "barrier_fwd"


def lit_primary(system, verified_only=True):
    """The primary literature barrier of a reaction (SN2: from the complex).  verified_only=True: None for unverified rows (L6)."""
    for l in LIT:
        if l["system"] == system and l["quantity"] == primary_q(system) and lit_scorable(l):
            if verified_only and not lit_in_stats(l):
                return None
            return l
    return None


def boot_ci(x, n=2000, seed=12345):
    """95% percentile bootstrap interval of the mean absolute error, resampling single values (kept for simple lists)."""
    x = np.abs(np.asarray(x, float))
    if len(x) < 3:
        return (float("nan"), float("nan"))
    rng = np.random.RandomState(seed)
    m = x[rng.randint(0, len(x), (n, len(x)))].mean(axis=1)
    return tuple(float(v) for v in np.percentile(m, [2.5, 97.5]))


def boot_by_reaction(groups, n=2000, seed=12345):
    """95% bootstrap interval of the mean absolute error, resampling REACTIONS (not rows): `groups` = {reaction: [signed errors]}.
    The rows of one reaction are dependent (L3 = L1 - L2, forward and reverse barriers, SN2 from reactants and from the complex), so
    they stay together.  Fewer than 3 reactions: no interval."""
    keys = [k for k, v in groups.items() if len(v)]
    if len(keys) < 3:
        return (float("nan"), float("nan"))
    rng = np.random.RandomState(seed)
    arr = [np.abs(np.asarray(groups[k], float)) for k in keys]
    means = []
    for _ in range(n):
        pick = rng.randint(0, len(keys), len(keys))
        means.append(float(np.concatenate([arr[i] for i in pick]).mean()))
    return tuple(float(v) for v in np.percentile(means, [2.5, 97.5]))


TIE = 0.5                       # kcal/mol: two mean absolute errors (or per-reaction errors) within this are a TIE (review row 8, F5)
CLASS_ORDER = ["H shift", "inversion", "elimination", "group shift", "SN2"]


def stat_methods(D):
    """Methods that get a line in Table 2: the five functionals (Skala with and without D3), the CCSD(T) reference (corrected and raw)."""
    have_ = {k[2] for k in D.idx}
    return [m for m in (M_SKALA, M_SKALA_NOD3, M_WB97MV, M_WB97XV, M_B3LYP, M_R2SCAN, M_PBE, M_CCSDT_BC, M_CCSDT_FC) if m in have_]


def headline_rows(D, method, kind="main"):
    """ONE primary barrier per reaction with a verified reference (SN2 from the complex).  Returns a list of dicts."""
    out = []
    for s in REACTIONS:
        l = lit_primary(s)
        if l is None:
            continue
        c = compare_to_lit(D, l, method, basis_for(s, kind))
        if c is not None and c.get("value") is not None:
            out.append(dict(system=s, cls=SYSTEMS[s]["cls"], lit_id=l["id"], computed=c["value"], literature=l["value"], error=c["diff"], level=c["level"]))
    return out


def pooled_rows(D, method, kind="main"):
    """Every verified literature row except the dependent L3 (supplement statistics)."""
    out = []
    for l in LIT:
        if not lit_in_stats(l):
            continue
        c = compare_to_lit(D, l, method, basis_for(l["system"], kind))
        if c is not None and c.get("value") is not None:
            out.append(dict(system=l["system"], cls=SYSTEMS[l["system"]]["cls"], lit_id=l["id"], computed=c["value"], literature=l["value"], error=c["diff"], level=c["level"]))
    return out


def _stat_line(method, kind, label, rows, by="system"):
    if not rows:
        return None
    e = np.array([r["error"] for r in rows], float)
    groups = {}
    for r in rows:
        groups.setdefault(r[by], []).append(r["error"])
    lo, hi = boot_by_reaction(groups)
    return dict(method=method, basis_kind=kind, set=label, n_rows=len(e), n_reactions=len(groups), MAE=round(float(np.mean(np.abs(e))), 3), MSE=round(float(np.mean(e)), 3),
                MaxAE=round(float(np.max(np.abs(e))), 3), MAE_ci95_low=round(lo, 3), MAE_ci95_high=round(hi, 3),
                geometry_level="; ".join(sorted({r["level"] for r in rows})), reactions=", ".join(sorted(groups)))


# ------------------------------------------------------------------ tables
def table1_systems(D):
    rows = []
    for s in SYSTEM_ORDER:
        S = SYSTEMS[s]
        lits = [l for l in LIT if l["system"] == s]
        fmt = lambda l: f"{l['id']} {l['quantity']} " + ("(value not yet read)" if l["value"] is None else f"{l['value']:+.2f}") + (f" +/- {l['unc']}" if l["unc"] else "") + ("" if l["verified"] == "yes" else " [UNVERIFIED, not in statistics]")
        rows.append(dict(number=S["n"], system=S["name"], key=s, cls=S["cls"], tier=D.tiers.get(s, S["tier"]),
                         species=", ".join(SPECIES_KIND[s]), literature=("; ".join(fmt(l) for l in lits) or "none (negative result)"),
                         source=("; ".join(sorted({l["source"][:70] for l in lits})) or ""), zero_point=S["zpe"], verified=", ".join(sorted({l["verified"] for l in lits})) or ""))
    return rows


def table2_errors(D):
    """Table 2 (review row 8).  Summary lines per method and basis: (a) HEADLINE = one primary barrier per verified reaction (SN2 from the complex),
    (b) per class, (c) pooled verified rows without the dependent L3 (supplement).  n is printed with every mean; the bootstrap resamples reactions.
    The detail table has every literature row with signed error, reference, convention, geometry level and basis (supplement)."""
    out, detail = [], []
    for m in stat_methods(D):
        for kind in ("main", "large"):
            if m in (M_CCSDT_BC, M_CCSDT_FC) and kind == "large":
                continue
            for l in LIT:
                if not lit_scorable(l):
                    continue
                c = compare_to_lit(D, l, m, basis_for(l["system"], kind))
                if c is None or c.get("value") is None:
                    continue
                prim = lit_primary(l["system"]) is l
                detail.append(dict(method=m, basis_kind=kind, basis=basis_for(l["system"], kind), lit_id=l["id"], system=l["system"], quantity=l["quantity"],
                                   computed=round(c["value"], 3), literature=l["value"], error=round(c["diff"], 3), zero_point_convention=l["zp"], geometry_level=c["level"],
                                   verified=l["verified"], in_headline=bool(prim and lit_in_stats(l)), in_pooled=lit_in_stats(l), note=c["note"][:80]))
                num("table2", l["system"], l["quantity"] + " " + l["id"], m, basis_for(l["system"], kind), c["level"], c["value"], l["value"], c["note"][:60])
            h = headline_rows(D, m, kind)
            line = _stat_line(m, kind, "HEADLINE: one primary barrier per reaction", h)
            if line:
                out.append(line)
            for cls in CLASS_ORDER:
                cl = _stat_line(m, kind, f"class {cls}", [r for r in h if r["cls"] == cls])
                if cl:
                    out.append(cl)
            pl = _stat_line(m, kind, "pooled verified rows, L3 excluded (supplement)", pooled_rows(D, m, kind))
            if pl:
                out.append(pl)
    return out, detail


def skala_d3_table(D):
    """Skala with and without D3, per reaction (review section 4, item 7)."""
    rows = []
    for s in REACTIONS:
        q, mb = primary_q(s), basis_for(s, "main")
        a, b = quantity_info(D, s, q, M_SKALA, mb), quantity_info(D, s, q, M_SKALA_NOD3, mb)
        if a is None and b is None:
            continue
        l = lit_primary(s, verified_only=False)
        lit_ok = l is not None and lit_in_stats(l)
        ref = reference_quantity(D, s, q, mb)
        rows.append(dict(system=SYSTEMS[s]["name"], key=s, quantity=q, geometry_level=(a or b)["level"],
                         skala_with_D3=None if a is None else round(a["value"], 3), skala_without_D3=None if b is None else round(b["value"], 3),
                         D3_contribution=None if a is None or b is None else round(a["value"] - b["value"], 3),
                         literature=None if l is None or l["value"] is None else l["value"], literature_verified=bool(lit_ok),
                         error_with_D3=None if (a is None or not lit_ok) else round(a["value"] - l["value"], 3),
                         error_without_D3=None if (b is None or not lit_ok) else round(b["value"] - l["value"], 3),
                         CCSDT_reference=None if ref is None else round(ref["value"], 3), reference_is_basis_corrected=None if ref is None else ref["corrected"]))
        for tag, v in (("with D3", a), ("without D3", b)):
            if v is not None:
                num("skala_d3", s, q, M_SKALA if tag == "with D3" else M_SKALA_NOD3, mb, v["level"], v["value"], l["value"] if lit_ok else None)
    return rows


def table3_scout_cost(D, record=True):
    """Scout check per reaction (review row 9): ALL functionals at the reference (own) geometry, the xTB geometries and, for Skala, the Skala-relaxed
    geometry; xTB's own barrier; the split of the error into geometry penalty and functional error; measured cost per point."""
    scout = []
    for s in REACTIONS:
        q, mb = primary_q(s), basis_for(s, "main")
        lit = lit_primary(s)
        ref = reference_quantity(D, s, q, mb)
        xtb = quantity(D, s, q, M_XTB, "GFN2", level=LV_XTB)
        target = ref["value"] if ref else (lit["value"] if lit else None)
        tname = ("CCSD(T)+dMP2" if ref and ref["corrected"] else ("CCSD(T) raw" if ref else ("literature" if lit else "none")))
        for m in [x for x in (M_SKALA, M_SKALA_NOD3, M_WB97MV, M_WB97XV, M_B3LYP, M_R2SCAN, M_PBE) if any(k[2] == x for k in D.idx if k[0] == s)]:
            a = quantity(D, s, q, m, mb)
            b = quantity(D, s, q, m, mb, level=LV_XTB)
            c = quantity(D, s, q, m, mb, level=LV_SKALA) if m in (M_SKALA, M_SKALA_NOD3) else None
            if all(v is None for v in (a, b, c)):
                continue
            scout.append(dict(system=SYSTEMS[s]["name"], key=s, method=m, quantity=q, literature=(lit["value"] if lit else None), reference_used=tname, reference=target,
                              CCSDT_ref_geometry=(ref["value"] if ref else None), ref_geometry=a, xtb_geometry=b, relaxed_geometry=c, xtb_own_barrier=xtb,
                              geometry_penalty_xtb=(None if a is None or b is None else b - a),
                              functional_error_at_ref_geometry=(None if a is None or target is None else a - target),
                              total_error_at_xtb_geometry=(None if b is None or target is None else b - target),
                              skala_relaxed_minus_reference_geometry=(None if a is None or c is None else c - a)))
            for tag, v in (("reference geometry", a), ("xTB geometry", b), ("own-method relaxed geometry (TS = parabola vertex)", c)):
                if v is not None and record:
                    num("table3", s, q, m, mb, tag, v, lit["value"] if lit else None)
    cost = {}
    for fn in ("phase11_dft.csv", "phase12_basis.csv", "phase13_scout.csv", "phase14_skala_relaxed.csv", "phase04_scout_ref.csv"):
        for r in read_rows(fn):
            try:
                sec = float(r["seconds"])
            except Exception:
                continue
            if r["converged"] in ("1", "True") and float(r["seconds"] or 0) > 0:
                cost.setdefault((r["method"], r["hardware"].split(", GPU")[-1] if "GPU" in r["hardware"] else "CPU"), []).append(sec)
    cost_rows = [dict(method=m, hardware=h, points=len(v), mean_seconds=round(float(np.mean(v)), 2), median_seconds=round(float(np.median(v)), 2)) for (m, h), v in sorted(cost.items())]
    for c in (cost_rows if record else []):
        num("table3", "", "seconds per point", c["method"], "", c["hardware"], c["mean_seconds"], None, f"{c['points']} points, unit seconds")
    return scout, cost_rows


# ------------------------------------------------------------------ figures
def fig1_design(D, plt, standin):
    cols = [("stationary\nenergies", lambda r: r["file"] == "phase11_dft.csv" and not is_path_id(r["pid"])),
            ("path\nenergies", lambda r: r["file"] in ("phase11_dft.csv", "phase16_optional.csv") and is_path_id(r["pid"])),
            ("CCSD(T)\nrows", lambda r: r["method"].startswith("CCSD(T)")),
            ("MP2 basis\ncorrection", lambda r: r["file"] == "phase08_mp2_basis.csv"),
            ("basis\ncheck", lambda r: r["file"] == "phase12_basis.csv"),
            ("xTB scout\n(G3)", lambda r: r["file"] == "phase13_scout.csv"),
            ("Skala-\nrelaxed (G4)", lambda r: r["file"] == "phase14_skala_relaxed.csv")]
    M = np.array([[sum(1 for r in D.rows if r["system"] == s and f(r) and r["E"] is not None) for _, f in cols] for s in SYSTEM_ORDER])
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    ax.imshow((M > 0).astype(float), cmap="Blues", vmin=0, vmax=2.2, aspect="auto")
    ax.grid(False)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([c for c, _ in cols])
    ax.set_yticks(range(len(SYSTEM_ORDER)))
    ax.set_yticklabels([f"{SYSTEMS[s]['n']} {SYSTEMS[s]['name']}  ({SYSTEMS[s]['cls']}, tier {D.tiers.get(s, SYSTEMS[s]['tier'])})" for s in SYSTEM_ORDER])
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, str(M[i, j]) if M[i, j] else "-", ha="center", va="center", color=INK, fontsize=8)
    ax.xaxis.tick_top(); ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title("Figure 1. Study design: number of energies saved per system and calculation type", loc="left", fontsize=9, pad=28)
    for i, s in enumerate(SYSTEM_ORDER):
        for j, (c, _) in enumerate(cols):
            num("fig1", s, "rows saved: " + c.replace("\n", " "), "", "", "", M[i, j])
    return _save(fig, "fig1_design.png", standin)


def fig2_barrier_errors(D, plt, standin):
    methods = [m for m in _methods_present(D) if m not in (M_CCSDT_FC, M_CCSDT_BC)]
    systems = [s for s in REACTIONS if lit_primary(s)]                      # verified references only: CH3NC (L6) is left out
    data = {m: [] for m in methods + [M_CCSDT_BC]}
    for s in systems:
        l = lit_primary(s)
        for m in data:
            c = compare_to_lit(D, l, m, basis_for(s, "main"))
            if c is None and m == M_CCSDT_BC:                               # no MP2 basis correction for this system: raw CCSD(T)
                c = compare_to_lit(D, l, M_CCSDT_FC, basis_for(s, "main"))
            v = None if c is None or c.get("value") is None else c["diff"]
            data[m].append(v)
            if v is not None:
                num("fig2", s, l["quantity"] + " error vs " + l["id"], m, basis_for(s, "main"), "preferred", c["value"], l["value"])
    if not any(v is not None for m in methods for v in data[m]):
        FIG_STATUS["fig2"] = "skipped: no functional energies at stationary points yet"
        return _save(_nofig(plt, "Figure 2: no functional energies yet"), "fig2_barrier_errors.png", standin)
    fig, ax = plt.subplots(figsize=(8.6, 4.0))
    w = 0.8 / len(methods)
    x = np.arange(len(systems))
    for k, m in enumerate(methods):
        v = np.array([np.nan if y is None else y for y in data[m]], float)
        ax.bar(x + (k - (len(methods) - 1) / 2) * w, v, w * 0.9, color=METHOD_STYLE[m][0], label=m, edgecolor="white", linewidth=0.5)
    v = np.array([np.nan if y is None else y for y in data[M_CCSDT_BC]], float)
    if np.isfinite(v).any():
        ax.plot(x, v, linestyle="none", marker="P", color=INK, markersize=7, label="CCSD(T) + MP2 basis correction (tiers A, B)")
    ax.axhline(0, color=INK, lw=0.8)
    ax.axhspan(-1.0, 1.0, color="#f0efec", zorder=0)
    ax.set_xticks(x); ax.set_xticklabels([f"{SYSTEMS[s]['n']} {SYSTEMS[s]['name']}" for s in systems], rotation=20, ha="right")
    ax.set_ylabel("signed error vs literature (kcal/mol)")
    ax.set_title(f"Figure 2. Barrier errors by functional, own B3LYP geometries (def2-TZVP, def2-TZVPD for SN2; shaded: +/- 1 kcal/mol; n = {len(systems)} reactions with a verified reference; "
                 "Skala-1.1 = with D3)", loc="left", fontsize=8)
    ax.legend(ncol=3, loc="best", fontsize=8)
    return _save(fig, "fig2_barrier_errors.png", standin)


def _path_series(D, system, method, basis):
    return {r["pid"]: (r["E"], r["coord"]) for r in _path(D, system, method, basis)}


def fig3_paths(D, plt, standin):
    systems = [s for s in ("hcn", "nh3", "h2co", "ch4") if s in SYSTEMS]
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.2))
    made = 0
    for ax, s in zip(axes.ravel(), systems):
        mb = basis_for(s, "main")
        ref = _path_series(D, s, M_CCSDT_BC, mb)                              # basis-corrected reference (review row 5)
        raw_note = ""
        if len(ref) < 3:
            ref, raw_note = _path_series(D, s, M_CCSDT_FC, mb), " (raw CCSD(T): no MP2 correction)"
        ax.set_title(f"{SYSTEMS[s]['n']} {SYSTEMS[s]['name']}{raw_note}", loc="left", fontsize=8)
        label = SCANS[s]["label"] if s in SCANS else "C-H distance (A)"
        ax.set_xlabel(label); ax.set_ylabel("error vs CCSD(T)+dMP2 (kcal/mol)")
        if len(ref) < 3:
            ax.text(0.5, 0.5, "no CCSD(T) path rows", transform=ax.transAxes, ha="center", color=INK2); continue
        for m in _methods_present(D):
            if m in (M_CCSDT_FC, M_CCSDT_BC):
                continue
            ser = _path_series(D, s, m, mb)
            common = [p for p in sorted(ref) if p in ser]
            if len(common) < 3:
                continue
            e = np.array([(ser[p][0] - ser[common[0]][0]) - (ref[p][0] - ref[common[0]][0]) for p in common]) * HARTREE_KCAL
            xs = [ref[p][1] if ref[p][1] is not None else i for i, p in enumerate(common)]
            col, mk = METHOD_STYLE[m]
            ax.plot(xs, e, marker=mk, color=col, lw=1.2, ms=4, label=m)
            made += 1
            for p, x_, v in zip(common, xs, e):
                num("fig3", s, "path error " + p, m, mb, LV_SCAN, v, 0.0, f"coordinate {x_:.3f}; relative to first common point")
        ax.axhline(0, color=INK, lw=0.8)
    if made == 0:
        plt.close(fig)
        FIG_STATUS["fig3"] = "skipped: needs CCSD(T) path rows (CPU notebook C3) and path energies (G1)"
        return _save(_nofig(plt, "Figure 3: needs CCSD(T) paths and G1 paths"), "fig3_path_errors.png", standin)
    axes.ravel()[0].legend(fontsize=7, ncol=2)
    fig.suptitle("Figure 3. Error along the relaxed paths against the basis-corrected CCSD(T) (energy relative to the first path point)", x=0.01, ha="left", fontsize=9)
    fig.tight_layout()
    return _save(fig, "fig3_path_errors.png", standin)


def fig4_basis(D, plt, standin):
    methods = [m for m in _methods_present(D) if m != M_CCSDT_FC]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    made = 0
    for ax, (a_kind, b_kind, title) in zip(axes, (("main", "large", "main to large basis (TZVP to QZVPP)"), ("small", "main", "small to main basis (SVP to TZVP)"))):
        systems = [s for s in REACTIONS]
        w = 0.8 / max(1, len(methods))
        x = np.arange(len(systems))
        for k, m in enumerate(methods):
            v = []
            for s in systems:
                q = primary_q(s)
                a, b = quantity(D, s, q, m, basis_for(s, a_kind)), quantity(D, s, q, m, basis_for(s, b_kind))
                v.append(np.nan if a is None or b is None else b - a)
                if a is not None and b is not None:
                    num("fig4", s, f"basis shift {a_kind}->{b_kind} of {q}", m, f"{basis_for(s, a_kind)}->{basis_for(s, b_kind)}", "preferred", b - a)
            v = np.array(v, float)
            if np.isfinite(v).any():
                made += 1
            ax.bar(x + (k - (len(methods) - 1) / 2) * w, v, w * 0.9, color=METHOD_STYLE[m][0], label=m, edgecolor="white", linewidth=0.5)
        for yv in (-0.5, 0.5):
            ax.axhline(yv, color=INK2, lw=0.8, ls="--")
        ax.axhline(0, color=INK, lw=0.8)
        ax.set_xticks(x); ax.set_xticklabels([f"{SYSTEMS[s]['n']}" for s in systems])
        ax.set_xlabel("reaction number"); ax.set_title(title, loc="left", fontsize=9)
    axes[0].set_ylabel("barrier shift (kcal/mol)")
    axes[0].legend(fontsize=7, ncol=2)
    fig.suptitle("Figure 4. Basis-set shifts of the barrier (dashed: +/- 0.5 kcal/mol)", x=0.01, ha="left", fontsize=9)
    if made == 0:
        FIG_STATUS["fig4"] = "skipped: no basis-check rows yet (G2)"
        plt.close(fig)
        return _save(_nofig(plt, "Figure 4: no basis-check rows yet (G2)"), "fig4_basis_shifts.png", standin)
    fig.tight_layout()
    return _save(fig, "fig4_basis_shifts.png", standin)


def fig5_scout(D, plt, standin):
    systems = [s for s in REACTIONS]
    x = np.arange(len(systems))
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0), gridspec_kw=dict(width_ratios=[1.5, 1]))
    ax, ax2 = axes
    made = 0
    anc = anchor_method() if "anchor_method" in globals() else M_WB97MV
    for key, (lvl, mk, lab) in {"ref": (None, "o", "Skala (D3), own B3LYP geometry"), "xtb": (LV_XTB, "s", "Skala, xTB geometries (scout)"), "skala": (LV_SKALA, "^", "Skala, Skala-relaxed (TS = parabola vertex)")}.items():
        v = []
        for s in systems:
            val = quantity(D, s, primary_q(s), M_SKALA, basis_for(s, "main"), **({"level": lvl} if lvl else {}))
            v.append(np.nan if val is None else val)
            if val is not None:
                num("fig5", s, primary_q(s), M_SKALA, basis_for(s, "main"), lvl or "preferred", val, (lit_primary(s) or {}).get("value"))
        v = np.array(v, float)
        if np.isfinite(v).any():
            made += 1
        ax.plot(x, v, linestyle="none", marker=mk, color=METHOD_STYLE[M_SKALA][0], ms=7, mfc=(METHOD_STYLE[M_SKALA][0] if key == "ref" else "white"), label=lab)
    lv = [(lit_primary(s, verified_only=False) or {}).get("value", np.nan) for s in systems]
    lv = [np.nan if v is None else v for v in lv]
    ax.plot(x, lv, linestyle="none", marker="_", color=INK, ms=16, mew=2, label="literature")
    cc = [(reference_quantity(D, s, primary_q(s), basis_for(s, "main")) or {}).get("value") for s in systems]
    ax.plot(x, [np.nan if c is None else c for c in cc], linestyle="none", marker="P", color=INK, mfc="white", ms=7, label="CCSD(T) + dMP2")
    ax.set_xticks(x); ax.set_xticklabels([str(SYSTEMS[s]["n"]) for s in systems]); ax.set_xlabel("reaction number"); ax.set_ylabel("barrier (kcal/mol)")
    ax.legend(fontsize=7, loc="best")
    ax.set_title("three routes to the barrier", loc="left", fontsize=9)
    w = 0.38
    meths_pen = [m for m in (M_SKALA, anc, M_B3LYP, M_R2SCAN, M_PBE) if any(r["method"] == m and r["level"] == LV_XTB for r in D.rows)] or [M_SKALA, anc]
    w = 0.8 / max(1, len(meths_pen))
    for k, m in enumerate(meths_pen):
        pen = []
        for s in systems:
            a = quantity(D, s, primary_q(s), m, basis_for(s, "main"))
            b = quantity(D, s, primary_q(s), m, basis_for(s, "main"), level=LV_XTB)
            pen.append(np.nan if a is None or b is None else b - a)
            if a is not None and b is not None:
                num("fig5", s, "geometry penalty (xTB minus reference geometry)", m, basis_for(s, "main"), LV_XTB, b - a)
        ax2.bar(x + (k - (len(meths_pen) - 1) / 2) * w, pen, w * 0.9, color=METHOD_STYLE.get(m, ("#898781", "o"))[0], label=m, edgecolor="white")
    ax2.axhspan(-1, 1, color="#f0efec", zorder=0); ax2.axhline(0, color=INK, lw=0.8)
    ax2.set_xticks(x); ax2.set_xticklabels([str(SYSTEMS[s]["n"]) for s in systems]); ax2.set_xlabel("reaction number")
    ax2.set_ylabel("xTB minus reference geometry (kcal/mol)"); ax2.set_title("geometry penalty (shaded: +/- 1)", loc="left", fontsize=9)
    ax2.legend(fontsize=7)
    fig.suptitle("Figure 5. Scout check: Skala energies at xTB geometries against Skala's own barrier", x=0.01, ha="left", fontsize=9)
    if made == 0:
        FIG_STATUS["fig5"] = "skipped: no Skala barriers yet (G1/G3/G4)"
    fig.tight_layout()
    return _save(fig, "fig5_scout_check.png", standin)


def _v4_grid(a):
    return np.array([[np.nan if v is None else v for v in row] for row in a], float)


def fig6_hcn_map_ch4(D, plt, standin):
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.9), gridspec_kw=dict(width_ratios=[1, 1, 1.1]))
    cc_npz = path("phase05_hcn_grid.npz")
    have_grid = os.path.exists(cc_npz)
    if have_grid:
        z = np.load(cc_npz, allow_pickle=True)
        Ecc, theta, rr = z["E"], np.asarray(z["theta"], float), np.asarray(z["r"], float)

        def rel(E):
            hc = np.where(theta < 45)[0]
            i = np.unravel_index(np.nanargmin(np.where(np.isfinite(E[hc]), E[hc], np.nan)), E[hc].shape)
            return (E - E[hc][i]) * HARTREE_KCAL
        rc = rel(Ecc)
        for ax, (nm, G) in zip(axes[:2], (("Skala (stored v4 grid)", HCN_GRID_SKALA), ("PBE (stored v4 grid)", HCN_GRID_PBE))):
            try:
                Eg = _v4_grid(G)
                err = rel(Eg) - rc
                err = np.where(rc <= 120.0, err, np.nan)
                lim = np.nanmax(np.abs(err)) if np.isfinite(err).any() else 1.0
                lim = max(1.0, min(lim, 8.0))
                im = ax.imshow(err.T, origin="lower", aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim,
                               extent=[theta[0], theta[-1], rr[0], rr[-1]])
                ax.grid(False)
                fig.colorbar(im, ax=ax, label="error vs CCSD(T)/def2-SVP (kcal/mol)")
                ax.set_xlabel("theta (deg; 0 = HCN, 180 = HNC)"); ax.set_ylabel("r (A)"); ax.set_title(nm, loc="left", fontsize=9)
                n_ok = int(np.isfinite(err).sum())
                for (i, j), v in np.ndenumerate(err):
                    if np.isfinite(v):
                        num("fig6", "hcn_grid", "map error", nm.split()[0], "def2-SVP", LV_GRID, v, 0.0, f"theta {theta[i]:.1f}, r {rr[j]:.3f}")
                ax.text(0.02, 0.97, f"{n_ok} cells", transform=ax.transAxes, va="top", fontsize=7, color=INK2)
            except Exception as ex:
                ax.text(0.5, 0.5, "map not available:\n" + repr(ex)[:60], transform=ax.transAxes, ha="center", fontsize=7)
    else:
        for ax in axes[:2]:
            ax.text(0.5, 0.5, "no phase05_hcn_grid.npz\n(CPU notebook C5)", transform=ax.transAxes, ha="center", color=INK2); ax.set_axis_off()
        FIG_STATUS["fig6a"] = "skipped: no CCSD(T) HCN grid (C5)"
    ax = axes[2]
    mb = basis_for("ch4", "main")
    ref = _path_series(D, "ch4", M_CCSDT_FC, mb)
    made = 0
    if len(ref) >= 3:
        for m in _methods_present(D):
            if m == M_CCSDT_FC:
                continue
            ser = _path_series(D, "ch4", m, mb)
            common = [p for p in sorted(ref) if p in ser]
            if len(common) < 3:
                continue
            e = np.array([(ser[p][0] - ser[common[0]][0]) - (ref[p][0] - ref[common[0]][0]) for p in common]) * HARTREE_KCAL
            xs = [ref[p][1] if ref[p][1] is not None else i for i, p in enumerate(common)]
            col, mk = METHOD_STYLE[m]
            ax.plot(xs, e, marker=mk, color=col, lw=1.2, ms=4, label=m)
            made += 1
            for p, x_, v in zip(common, xs, e):
                num("fig6", "ch4", "stretch error " + p, m, mb, LV_SCAN, v, 0.0, f"C-H {x_:.3f} A")
        ax.axhline(0, color=INK, lw=0.8); ax.legend(fontsize=7, ncol=2)
    if made == 0:
        ax.text(0.5, 0.5, "no CH4 CCSD(T) or functional rows yet", transform=ax.transAxes, ha="center", color=INK2)
        FIG_STATUS["fig6b"] = "skipped: no CH4 stretch rows"
    ax.set_xlabel("C-H distance (A)"); ax.set_ylabel("error vs CCSD(T) (kcal/mol)"); ax.set_title("CH4 C-H stretch", loc="left", fontsize=9)
    fig.suptitle("Figure 6. HCN rigid-map error (left, middle) and CH4 stretch (right)", x=0.01, ha="left", fontsize=9)
    fig.tight_layout()
    return _save(fig, "fig6_hcn_map_ch4.png", standin)


# ------------------------------------------------------------------ analysis outputs of review section 4
# THE FROZEN OUTCOME CRITERIA.  Dated 2026-11-15, written before any real number exists; the same wording is in PLAN_experiments.md section 7.
# "supports" / "negative" / "inconclusive" are decided by exactly these rules.  Do not change them after a real number has been seen.
FROZEN_DATE = "2026-11-15"
FROZEN = {
    "F1": dict(finding="Barriers at relaxed stationary points, five functionals, headline set (one primary barrier per verified reaction, SN2 from the complex)",
               supports="Skala-1.1 (with D3) mean absolute error is at most 2.0 kcal/mol AND at most the wB97M-V (anchor) mean absolute error + 0.5",
               negative="Skala mean absolute error exceeds the better of B3LYP-D3(BJ) and r2SCAN by more than 0.5 kcal/mol",
               inconclusive="everything else (a tie within 0.5 with B3LYP or r2SCAN, or Skala above 2.0 without losing)"),
    "F2": dict(finding="Error along the relaxed Tier A paths (HCN, NH3, H2CO, CH4) against the basis-corrected CCSD(T)",
               supports="range (max minus min) of the signed Skala error is at most 1.0 kcal/mol on at least 3 of the 4 paths",
               negative="maximum absolute Skala error is at least 2.0 kcal/mol on at least 2 of the 4 paths",
               inconclusive="everything else; always inconclusive when fewer than 3 paths have a usable reference (R1 failed or rows missing)"),
    "F3": dict(finding="Basis-set stability of the primary barrier, def2-TZVP(D) to def2-QZVPP(D)",
               supports="Skala shift at most 0.5 kcal/mol AND Skala's rank among the five functionals unchanged, in at least 75% of the reactions with data (at least 4 reactions)",
               negative="Skala shift above 1.0 kcal/mol in at least 2 reactions",
               inconclusive="everything else, or fewer than 4 reactions with both bases"),
    "F4": dict(finding="Scout check: barriers at xTB geometries against the barrier at our own B3LYP geometries, every functional of G3",
               supports="geometry penalty below 1.0 kcal/mol in at least 75% of the reactions (6 of 8) for EVERY functional computed",
               negative="penalty above 3.0 kcal/mol in at least 3 reactions, or in every reaction of one class (SN2: all four)",
               inconclusive="everything else, or fewer than 4 reactions with xTB-geometry energies"),
    "F5": dict(finding="Robustness of Skala's rank among the five functionals (by absolute error; errors within 0.5 kcal/mol of each other tie)",
               supports="the most frequent rank (ties share a rank) occurs in at least 75% of the reactions (6 of 8)",
               negative="the most frequent rank occurs in fewer than 50% of the reactions",
               inconclusive="everything else, or fewer than 4 reactions"),
    "F6": dict(finding="Negative results kept: CH4 C-H stretch and the rigid HCN map against CCSD(T)/def2-SVP",
               supports="T1 below 0.02 (and RHF stable at the CH4 points) at at least 90% of the CH4 points and HCN grid cells in the 120 kcal/mol window",
               negative="that fraction is below 50%: the reference is truncated at the distance where it fails, and this is said in the paper",
               inconclusive="everything else"),
}


def _rank_with_ties(errs):
    """Rank (1 = smallest absolute error) of each method; errors within TIE of the better neighbour share the rank."""
    items = sorted(errs.items(), key=lambda kv: abs(kv[1]))
    ranks, cur, lead = {}, 0, None
    for m, e in items:
        if lead is None or abs(e) - lead > TIE:
            cur, lead = cur + 1, abs(e)
        ranks[m] = cur
    return ranks


def _functionals(D):
    """The five functionals of the paper: Skala (with D3), the anchor actually present, B3LYP, r2SCAN, PBE."""
    have_ = {k[2] for k in D.idx}
    anc = M_WB97MV if M_WB97MV in have_ else (M_WB97XV if M_WB97XV in have_ else M_WB97MV)
    return [M_SKALA, anc, M_B3LYP, M_R2SCAN, M_PBE]


def _path_errors(D, system, method, ref_method):
    mb = basis_for(system, "main")
    ref, ser = _path_series(D, system, ref_method, mb), _path_series(D, system, method, mb)
    common = [p for p in sorted(ref) if p in ser]
    if len(common) < 3:
        return None
    e = np.array([(ser[p][0] - ser[common[0]][0]) - (ref[p][0] - ref[common[0]][0]) for p in common]) * HARTREE_KCAL
    xs = [ref[p][1] if ref[p][1] is not None else float(i) for i, p in enumerate(common)]
    peak = int(np.argmax([ref[p][0] for p in common]))
    return dict(points=common, err=e, coord=xs, ref_peak_index=peak)


def path_statistics(D):
    """Review section 4, item 3: per Tier A system and method, error statistics along the path against the basis-corrected and the raw CCSD(T)."""
    rows = []
    for s in ("hcn", "nh3", "h2co", "ch4"):
        for ref_m, tag in ((M_CCSDT_BC, "CCSD(T)+dMP2"), (M_CCSDT_FC, "CCSD(T) raw")):
            for m in [x for x in (M_SKALA, M_SKALA_NOD3, M_WB97MV, M_WB97XV, M_B3LYP, M_R2SCAN, M_PBE) if any(k[2] == x and k[0] == s for k in D.idx)]:
                r = _path_errors(D, s, m, ref_m)
                if r is None:
                    continue
                e, k = r["err"], int(np.argmax(np.abs(r["err"])))
                where = "no reaction peak" if s == "ch4" else ("at the TS (+/- 1 point)" if abs(k - r["ref_peak_index"]) <= 1 else ("before the TS" if k < r["ref_peak_index"] else "after the TS"))
                rows.append(dict(system=SYSTEMS[s]["name"], key=s, method=m, reference=tag, n_points=len(e), mean_abs_error=round(float(np.mean(np.abs(e))), 3),
                                 max_abs_error=round(float(np.max(np.abs(e))), 3), max_at_coordinate=round(float(r["coord"][k]), 4), max_location=where,
                                 error_min=round(float(e.min()), 3), error_max=round(float(e.max()), 3), error_range=round(float(e.max() - e.min()), 3)))
    return rows


def basis_table(D):
    """Review section 4, item 4: barrier shifts SVP, TZVP, QZVPP per method; rank change yes or no (Skala's rank among the functionals, main against large basis)."""
    rows, fn = [], _functionals(D)
    for s in REACTIONS:
        q = primary_q(s)
        b = {k: basis_for(s, k) for k in ("small", "main", "large")}
        l = lit_primary(s)
        err = {kind: {} for kind in ("main", "large")}
        for m in fn:
            for kind in ("main", "large"):
                c = compare_to_lit(D, l, m, b[kind]) if l else None
                if c is not None and c.get("value") is not None:
                    err[kind][m] = c["diff"]
        rk = {kind: (_rank_with_ties(err[kind]) if len(err[kind]) == len(fn) else None) for kind in err}
        for m in [x for x in (M_SKALA, M_SKALA_NOD3) + tuple(fn[1:]) if any(k[2] == x and k[0] == s for k in D.idx)]:
            v = {k: quantity(D, s, q, m, b[k]) for k in b}
            if all(x is None for x in v.values()):
                continue
            sh = lambda a_, b_: None if v[a_] is None or v[b_] is None else round(v[b_] - v[a_], 3)
            changed = None
            if rk["main"] and rk["large"] and m in rk["main"] and m in rk["large"]:
                changed = "yes" if rk["main"][m] != rk["large"][m] else "no"
            rows.append(dict(system=SYSTEMS[s]["name"], key=s, method=m, quantity=q, small_basis=b["small"], barrier_small=None if v["small"] is None else round(v["small"], 3),
                             main_basis=b["main"], barrier_main=None if v["main"] is None else round(v["main"], 3), large_basis=b["large"],
                             barrier_large=None if v["large"] is None else round(v["large"], 3), shift_small_to_main=sh("small", "main"), shift_main_to_large=sh("main", "large"),
                             rank_main=None if not rk["main"] or m not in rk["main"] else rk["main"][m], rank_large=None if not rk["large"] or m not in rk["large"] else rk["large"][m],
                             rank_changed=changed))
    return rows


def scout_table(D):
    """Review section 4, item 5 (same content as Table 3, one line per reaction and functional)."""
    sc, _ = table3_scout_cost(D, record=False)
    return sc


def geometry_table(D):
    """Review section 4, item 6: the driving coordinate of the transition state (the key distance or angle) at every stored geometry level."""
    rows = []
    for s in REACTIONS:
        sp = SCANS.get(s)
        if not sp:
            continue
        r = dict(system=SYSTEMS[s]["name"], key=s, coordinate=sp["label"])
        for lv in (LV_OWN, LV_SCANPEAK, LV_DB, LV_XTB, LV_SKALA):
            g = geometry_of(s, "TS", lv)
            r[lv] = None if g is None else round(coord_value(sp["kind"], sp["idx"], g["atoms"]), 4)
        if r.get(LV_SKALA) is not None and r.get(LV_OWN) is not None:
            r["skala_minus_b3lyp"] = round(r[LV_SKALA] - r[LV_OWN], 4)
        if r.get(LV_XTB) is not None and r.get(LV_OWN) is not None:
            r["xtb_minus_b3lyp"] = round(r[LV_XTB] - r[LV_OWN], 4)
        if r.get(LV_DB) is not None and r.get(LV_OWN) is not None:
            r["database_minus_b3lyp"] = round(r[LV_DB] - r[LV_OWN], 4)
        if any(r.get(lv) is not None for lv in (LV_OWN, LV_SCANPEAK, LV_DB, LV_XTB, LV_SKALA)):
            rows.append(r)
    return rows


def reference_quality(D):
    """Review section 4, item 8: T1, RHF stability, frozen core, CCSD(T) against the literature (raw and basis-corrected)."""
    rows = []
    for s in REACTIONS:
        if D.tiers.get(s, SYSTEMS[s]["tier"]) not in ("A", "B"):
            continue
        mb = basis_for(s, "main"); q = primary_q(s)
        cc = [r for r in D.rows if r["system"] == s and r["method"] == M_CCSDT_FC and r["basis"] == mb and not is_path_id(r["pid"])]
        t1 = [r["t1"] for r in cc if r["t1"] is not None]
        stab = [r["x"].get("stable") for r in cc if r["x"].get("stable") not in (None, "")]
        l = lit_primary(s, verified_only=False)
        ref = reference_quantity(D, s, q, mb)
        ae, fc = quantity(D, s, "barrier_fwd", M_CCSDT_AE, mb), quantity(D, s, "barrier_fwd", M_CCSDT_FC, mb)
        rows.append(dict(system=SYSTEMS[s]["name"], key=s, tier=D.tiers.get(s, SYSTEMS[s]["tier"]), quantity=q, max_T1=None if not t1 else round(max(t1), 4),
                         rhf_stable_all=None if not stab else all(str(v) not in ("0", "False") for v in stab), frozen_core=True,
                         all_electron_minus_frozen_core=None if ae is None or fc is None else round(ae - fc, 3),
                         CCSDT_raw=None if ref is None else ref["raw"], CCSDT_plus_dMP2=None if ref is None or not ref["corrected"] else round(ref["value"], 3),
                         literature=None if l is None else l["value"], literature_verified=None if l is None else l["verified"],
                         raw_minus_literature=None if ref is None or ref["raw"] is None or l is None or l["value"] is None else round(ref["raw"] - l["value"], 3),
                         corrected_minus_literature=None if ref is None or not ref["corrected"] or l is None or l["value"] is None else round(ref["value"] - l["value"], 3),
                         geometry_level=None if ref is None else ref["level"]))
    return rows


def excluded_and_failed(D):
    """Review section 4, item 9: every point that is left out or failed, with the reason."""
    rows = []
    for r in D.rows:
        why = None
        if not r["conv"] or r["E"] is None:
            why = "SCF / calculation not converged" + (f": {r['x'].get('error')}" if r["x"].get("error") else "")
        elif (r["system"], r["pid"], r["level"], r["basis"], r["method"]) in D.excluded:
            why = f"excluded from reference statistics: T1 {r['t1']} or unstable RHF"
        elif r["x"].get("interior") is False:
            why = "Skala TS: no interior maximum in the window (parabola vertex not usable)"
        if why:
            rows.append(dict(system=r["system"], point=r["pid"], level=r["level"], method=r["method"], basis=r["basis"], file=r["file"], reason=why))
    for key, v in D.mixed.items():
        if not v["aligned"]:
            rows.append(dict(system=key[0], point=",".join(v["partners_missing"]), level="", method=key[1], basis=key[2], file="",
                             reason="barrier would mix density-fitted and unfitted energies; partner not recomputed"))
    p = path("phase_errors.csv")
    if os.path.exists(p):
        for r in csv.DictReader(open(p)):
            rows.append(dict(system=r["system"], point="", level="", method="", basis="", file="phase_errors.csv", reason=f"{r['phase']}: {r['error'][:140]}"))
    return rows


def _t1_fraction(D):
    """Fraction of CH4 stretch points with T1 < 0.02 and RHF stable, and of HCN grid cells (<= 120 kcal/mol) with T1 < 0.02."""
    out = {}
    ch4 = [r for r in D.rows if r["system"] == "ch4" and r["method"] == M_CCSDT_FC and is_path_id(r["pid"]) and r["conv"]]
    if ch4:
        good = [r for r in ch4 if (r["t1"] is not None and r["t1"] < TOL["I2"][0]) and r["x"].get("stable") not in (0, "0", False)]
        out["ch4"] = (len(good), len(ch4))
    z = path("phase05_hcn_grid.npz")
    if os.path.exists(z):
        try:
            d = np.load(z, allow_pickle=True)
            E, t1, theta = d["E"], d["t1"], np.asarray(d["theta"], float)
            ok = np.isfinite(E)
            hc = np.where(theta < 45)[0]
            emin = np.nanmin(E[hc])
            win = ok & ((E - emin) * HARTREE_KCAL <= 120.0)
            out["hcn_grid"] = (int(np.sum(win & np.isfinite(t1) & (t1 < TOL["I2"][0]))), int(np.sum(win)))
        except Exception:
            pass
    return out


def findings_outcome(D):
    """Review section 4, item 1: F1 to F6 against the frozen criteria: measured value, n, verdict, one sentence for the paper.
    A verdict is only as good as the data behind it: with missing data it says 'inconclusive (data missing)'."""
    rows = []
    fn = _functionals(D)
    anc = fn[1]

    def add(fid, measured, n, verdict, sentence):
        f = FROZEN[fid]
        rows.append(dict(finding=fid, description=f["finding"], criterion_supports=f["supports"], criterion_negative=f["negative"], criterion_inconclusive=f["inconclusive"],
                         frozen_on=FROZEN_DATE, measured=measured, n=n, verdict=verdict, sentence_for_the_paper=sentence))

    # ---- F1
    mae = {}
    for m in fn + [M_SKALA_NOD3]:
        h = headline_rows(D, m)
        if h:
            mae[m] = (float(np.mean([abs(r["error"]) for r in h])), len(h))
    if M_SKALA in mae and anc in mae:
        S, W = mae[M_SKALA][0], mae[anc][0]
        B = min([mae[m][0] for m in (M_B3LYP, M_R2SCAN) if m in mae] or [float("nan")])
        if S <= 2.0 and S <= W + 0.5:
            v = "supports"
        elif B == B and S > B + 0.5:
            v = "negative"
        else:
            v = "inconclusive"
        add("F1", "MAE (kcal/mol, n): " + ", ".join(f"{m} {a:.2f} ({n})" for m, (a, n) in mae.items()), mae[M_SKALA][1], v,
            f"On {mae[M_SKALA][1]} reactions with a verified reference Skala-1.1 (with D3) has a mean absolute barrier error of {S:.2f} kcal/mol against {W:.2f} for {anc}"
            + (f" and {B:.2f} for the better of B3LYP-D3(BJ) and r2SCAN" if B == B else "") + ".")
    else:
        add("F1", "Skala or anchor headline rows missing", 0, "inconclusive (data missing)", "")
    # ---- F2
    sk_paths, n_valid = [], 0
    for s in ("hcn", "nh3", "h2co", "ch4"):
        r = _path_errors(D, s, M_SKALA, M_CCSDT_BC)
        if r is not None:
            n_valid += 1
            sk_paths.append((s, float(r["err"].max() - r["err"].min()), float(np.max(np.abs(r["err"])))))
    if n_valid:
        flat = sum(1 for _, rg, _m in sk_paths if rg <= 1.0)
        big = sum(1 for _, _r, mx in sk_paths if mx >= 2.0)
        v = "inconclusive" if n_valid < 3 else ("supports" if flat >= 3 else ("negative" if big >= 2 else "inconclusive"))
        add("F2", "; ".join(f"{s}: range {rg:.2f}, max |error| {mx:.2f}" for s, rg, mx in sk_paths), n_valid, v,
            f"Along {n_valid} relaxed paths the Skala error against basis-corrected CCSD(T) has a range of {min(r for _, r, _ in sk_paths):.2f} to {max(r for _, r, _ in sk_paths):.2f} kcal/mol.")
    else:
        add("F2", "no Skala and CCSD(T) path rows", 0, "inconclusive (data missing)", "")
    # ---- F3
    shifts, rank_same, n_b = [], 0, 0
    for r in basis_table(D):
        if r["method"] == M_SKALA and r["shift_main_to_large"] is not None:
            n_b += 1
            shifts.append((r["key"], r["shift_main_to_large"]))
            if abs(r["shift_main_to_large"]) <= 0.5 and r["rank_changed"] in ("no", None):
                rank_same += 1
    if n_b:
        big = sum(1 for _, x in shifts if abs(x) > 1.0)
        v = "inconclusive" if n_b < 4 else ("supports" if rank_same >= math.ceil(0.75 * n_b) else ("negative" if big >= 2 else "inconclusive"))
        add("F3", "Skala shifts main to large basis: " + ", ".join(f"{k} {x:+.2f}" for k, x in shifts), n_b, v,
            f"The Skala barrier moves by at most {max(abs(x) for _, x in shifts):.2f} kcal/mol from the triple-zeta to the quadruple-zeta basis in {n_b} reactions.")
    else:
        add("F3", "no basis-check rows", 0, "inconclusive (data missing)", "")
    # ---- F4
    pen, n_react = {}, set()
    for r in table3_scout_cost(D, record=False)[0]:
        if r["geometry_penalty_xtb"] is not None and r["method"] in fn:
            pen.setdefault(r["method"], {})[r["key"]] = r["geometry_penalty_xtb"]
            n_react.add(r["key"])
    if pen:
        n_min = min(len(v) for v in pen.values())
        thr = math.ceil(0.75 * len(n_react))
        small_all = all(sum(1 for x in v.values() if abs(x) < 1.0) >= thr for v in pen.values())
        allp = {k: x for v in pen.values() for k, x in v.items()}
        n_big = len({k for v in pen.values() for k, x in v.items() if abs(x) > 3.0})
        sn2_all = all(any(abs(pen[m].get(k, 0)) > 3.0 for m in pen) for k in ("sn2_ff", "sn2_clcl", "sn2_fcl", "sn2_ohf"))
        v = "inconclusive" if len(n_react) < 4 else ("supports" if small_all else ("negative" if (n_big >= 3 or sn2_all) else "inconclusive"))
        add("F4", "geometry penalty (kcal/mol): " + "; ".join(f"{m}: " + ", ".join(f"{k} {x:+.2f}" for k, x in sorted(vv.items())) for m, vv in pen.items()), len(n_react), v,
            f"Barriers at GFN2-xTB geometries differ from those at B3LYP geometries by at most {max(abs(x) for x in allp.values()):.2f} kcal/mol over {len(n_react)} reactions and {len(pen)} functionals.")
    else:
        add("F4", "no xTB-geometry energies", 0, "inconclusive (data missing)", "")
    # ---- F5
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
        add("F5", "needs all five functionals at the stationary points", 0, "inconclusive (data missing)", "")
    # ---- F6
    fr = _t1_fraction(D)
    if fr:
        good, tot = sum(a for a, _ in fr.values()), sum(b for _, b in fr.values())
        frac = good / tot if tot else float("nan")
        v = "supports" if frac >= 0.9 else ("negative" if frac < 0.5 else "inconclusive")
        sk = _path_errors(D, "ch4", M_SKALA, M_CCSDT_BC) or _path_errors(D, "ch4", M_SKALA, M_CCSDT_FC)
        loc = "" if sk is None else f" The largest Skala error along the CH4 stretch is {np.max(np.abs(sk['err'])):.2f} kcal/mol at C-H = {sk['coord'][int(np.argmax(np.abs(sk['err'])))]:.2f} A."
        add("F6", "reference reliable (T1 < 0.02) at: " + ", ".join(f"{k} {a}/{b}" for k, (a, b) in fr.items()), tot, v,
            f"CCSD(T) is reliable (T1 below 0.02) at {good} of {tot} points of the CH4 stretch and the HCN map." + loc)
    else:
        add("F6", "no CH4 / HCN grid CCSD(T) rows", 0, "inconclusive (data missing)", "")
    return rows


# ------------------------------------------------------------------ driver
def g5_make_all(standin=False):
    """Make every figure and table.  Each one is guarded: a failure is logged and the next one still runs."""
    plt = _plt()
    D = Data()
    print(f"  G5 reads {len(D.rows)} rows from {len(D.files)} CSV files")
    NUMS.clear(); FIG_STATUS.clear()
    tables = {}
    with item_guard("phase15", "table1"):
        tables["table1_systems"] = table1_systems(D); FIG_STATUS["table1"] = "made"
    with item_guard("phase15", "table2"):
        t2, det = table2_errors(D)
        tables["table2_errors"], tables["table2_detail"] = t2, det; FIG_STATUS["table2"] = "made" if t2 else "empty (no functional rows yet)"
    with item_guard("phase15", "table3"):
        sc, cost = table3_scout_cost(D)
        tables["table3_scout"], tables["table3_cost"] = sc, cost; FIG_STATUS["table3"] = "made" if sc or cost else "empty"
    for key, fn in (("skala_d3_table", skala_d3_table), ("path_statistics", path_statistics), ("basis_table", basis_table), ("scout_table", scout_table),
                    ("geometry_table", geometry_table), ("reference_quality", reference_quality), ("excluded_and_failed", excluded_and_failed),
                    ("findings_outcome", findings_outcome)):
        with item_guard("phase15", key):
            tables[key] = fn(D); FIG_STATUS[key] = "made" if tables[key] else "empty (data missing)"
    for name, fn in (("fig1", fig1_design), ("fig2", fig2_barrier_errors), ("fig3", fig3_paths), ("fig4", fig4_basis), ("fig5", fig5_scout), ("fig6", fig6_hcn_map_ch4)):
        with item_guard("phase15", name):
            fn(D, plt, standin)
            FIG_STATUS.setdefault(name, "made")
    # CSV files
    def _write(name, rows):
        if not rows:
            return
        keys = []
        for r in rows:
            for k in r:
                if k not in keys:
                    keys.append(k)
        with open(path(name + ".csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    for k, rows in tables.items():
        _write(k, rows)
    keys = ["table", "system", "quantity", "method", "basis", "geometry_level", "value", "reference", "error", "unit", "note"]
    with open(path("numbers_v5.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(NUMS)
    print(f"  numbers_v5.csv: {len(NUMS)} numbers")
    try:                                                      # one workbook with every table (needs openpyxl; the CSV files are the fallback)
        import pandas as pd
        with pd.ExcelWriter(path("summary_tables.xlsx"), engine="openpyxl") as xw:
            for k, rows in tables.items():
                if rows:
                    pd.DataFrame(rows).to_excel(xw, sheet_name=k[:31], index=False)
            pd.DataFrame(NUMS).to_excel(xw, sheet_name="numbers", index=False)
        print("  wrote summary_tables.xlsx")
    except Exception as ex:
        print("  summary_tables.xlsx not written (", repr(ex)[:80], ") - the table CSV files hold the same numbers")
    write_json("phase15_figures_status.json", FIG_STATUS)
    for k, v in FIG_STATUS.items():
        print(f"   {k}: {v}")
    return FIG_STATUS
