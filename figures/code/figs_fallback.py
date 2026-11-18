"""Figures 2, 3, 5, 6, 7 (fallback set first)."""
from figcommon import *
from matplotlib.gridspec import GridSpec


# ------------------------------------------------------------------ Figure 2
def fig02(D):
    print("== Figure 2")
    b = D["bar"]
    h = b[b.in_headline == True]
    meths = DFT6 + ["CCSD(T)+dMP2"]
    fig, ax = fig_ax(170, 90)
    band(ax)
    zero_line(ax)
    n = len(meths)
    step = 0.115
    rows = []
    for i, m in enumerate(meths):
        off = (i - (n - 1) / 2) * step
        for k, r in enumerate(RX):
            q = h[(h.reaction == r) & (h.method == m)]
            if q.empty:
                continue
            y = float(q.signed_error_main_kcal.iloc[0])
            rows.append(dict(reaction=r, method=m, x=k + off, signed_error_main_kcal=y))
            ax.plot([k + off, k + off], [0, y], color=METH[m]["c"], lw=0.8, zorder=3)
            ax.plot([k + off], [y], ls="none", zorder=4, **mstyle(m))
    data = pd.DataFrame(rows)
    save_data(data, "fig02_signed_errors")
    # annotations
    def xy(r, m):
        q = data[(data.reaction == r) & (data.method == m)].iloc[0]
        return q.x, q.signed_error_main_kcal
    x0, y0 = xy("hcn", "wB97M-V")
    ax.annotate("unexplained,\nsee text", (x0, y0), xytext=(x0 + 0.17, y0 - 1.3), fontsize=8, ha="left", va="top",
                arrowprops=dict(arrowstyle="-", lw=0.5, color="black"))
    x1, y1 = xy("sn2_clcl", "Skala-1.1")
    ax.annotate("SCF and geometry\ncaveats, Fig. 7", (x1, y1 - 0.2), xytext=(x1 + 0.0, -6.9), fontsize=8, ha="right", va="top",
                arrowprops=dict(arrowstyle="-", lw=0.5, color="black"))
    ax.text(1.008, 0, "literature\nreference", transform=ax.get_yaxis_transform(), fontsize=8, ha="left", va="center", linespacing=0.95)
    ax.set_xlim(-0.5, 4.5)
    ax.set_ylim(-9, 3)
    rx_ticklabels(ax)
    ax.set_xlabel("reaction")
    ax.set_ylabel(f"signed error of the\nbarrier ({KC})")
    ax.text(0.995, 0.985, "n = 5 reactions", transform=ax.transAxes, ha="right", va="top", fontsize=8)
    hd = [handle(m, line=False) for m in meths]
    fig.legend(handles=hd, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.0),
               handletextpad=0.3, columnspacing=1.2)
    fig.subplots_adjust(left=0.115, right=0.895, top=0.98, bottom=0.30)
    save(fig, "fig02_signed_errors")
    sk = data[data.method == "Skala-1.1"].set_index("reaction").signed_error_main_kcal
    for r, e in zip(RX, [0.30, 0.81, 1.08, 1.91, -3.12]):
        anchor("fig02", f"Skala error {r}", e, sk[r], 0.006)
    print(data.pivot(index="method", columns="reaction", values="signed_error_main_kcal").round(3)[RX])
    neg = {m: int((data[data.method == m].signed_error_main_kcal < 0).sum()) for m in DFT6}
    print("  negative counts", neg)


# ------------------------------------------------------------------ Figure 3
def fig03(D):
    print("== Figure 3")
    s = D["sum"]
    s = s[(s.set == "HEADLINE") & (s.basis_kind == "main")].set_index("method")
    b = D["bar"]; h = b[b.in_headline == True]
    cnt = {m: int((h[h.method == m].signed_error_main_kcal < 0).sum()) for m in DFT6}
    fig = plt.figure(figsize=(170 * MM, 92 * MM))
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.0, 0.85, 1.05], left=0.135, right=0.995, bottom=0.41, top=0.97, wspace=0.34)
    ax = fig.add_subplot(gs[0]); bx = fig.add_subplot(gs[1], sharey=ax); cx = fig.add_subplot(gs[2])
    ypos = np.arange(len(DFT6))[::-1]
    rows = []
    for yy, m in zip(ypos, DFT6):
        st = METH[m]
        r = s.loc[m]
        kw = dict(height=0.62, edgecolor="black", linewidth=0.4, zorder=3)
        if st["open"]:
            kw.update(facecolor="white", edgecolor=st["c"], hatch="////", linewidth=0.8)
        else:
            kw.update(facecolor=st["c"])
        ax.barh(yy, r.MAE, **kw)
        ax.plot([r.MAE_ci95_low, r.MAE_ci95_high], [yy, yy], color="black", lw=0.8, zorder=4)
        for xe in (r.MAE_ci95_low, r.MAE_ci95_high):
            ax.plot([xe, xe], [yy - 0.12, yy + 0.12], color="black", lw=0.8, zorder=4)
        ax.text(r.MAE_ci95_high + 0.2, yy, f"{r.MAE:.2f}", va="center", ha="left", fontsize=8)
        bkw = dict(kw)
        bx.barh(yy, r.MSE, **bkw)
        rows.append(dict(panel="a/b", method=m, MAE=r.MAE, MAE_ci95_low=r.MAE_ci95_low, MAE_ci95_high=r.MAE_ci95_high, MSE=r.MSE,
                         n_negative_of_5=cnt[m]))
        if r.MSE < 0:
            bx.text(0.15, yy, f"{cnt[m]} of 5\nnegative", va="center", ha="left", fontsize=8, linespacing=0.95)
        else:
            bx.text(r.MSE + 0.15, yy, f"{cnt[m]} of 5\nnegative", va="center", ha="left", fontsize=8, linespacing=0.95)
    ax.set_yticks(ypos)
    ax.set_yticklabels(["Skala-1.1", "Skala-1.1\nwithout D3", "ωB97M-V", "B3LYP-D3(BJ)", "r²SCAN", "PBE"])
    bx.tick_params(labelleft=False)
    ax.set_xlim(0, 9.0); ax.set_xticks(range(0, 9, 2))
    bx.set_xlim(-5.2, 3.6); bx.set_xticks([-4, -2, 0, 2])
    ax.set_xlabel(f"mean absolute error\n({KC})")
    bx.set_xlabel(f"mean signed error\n({KC})")
    bx.axvline(0, color="black", lw=0.8, zorder=2)
    ax.set_ylim(-0.6, 6.4)
    for a in (ax, bx):
        a.spines["left"].set_visible(a is ax)
    bx.tick_params(left=False)
    ax.text(0.99, 0.985, "n = 5 reactions", transform=ax.transAxes, ha="right", va="top", fontsize=8)
    panel(ax, "a", x=0.02, y=1.0)
    panel(bx, "b", x=0.02, y=1.0)
    # (c) leave one out
    L = D["leave_one_out"]
    order = ["none", "hcn", "nh3", "h2co", "sn2_ff", "sn2_clcl"]
    tl = ["all five", "without HCN", "without NH₃", "without H₂CO", "without F⁻ + CH₃F", "without Cl⁻ + CH₃Cl"]
    for m in DFT5:
        q = L[L.method == m].set_index("removed").loc[order]
        cx.plot(range(6), q.MAE.values, color=METH[m]["c"], lw=0.8, zorder=2)
        cx.plot(range(6), q.MAE.values, ls="none", zorder=3, **mstyle(m))
        for o, v in zip(order, q.MAE.values):
            rows.append(dict(panel="c", method=m, removed=o, MAE=v))
    cx.set_xticks(range(6)); cx.set_xticklabels(tl, rotation=45, ha="right", rotation_mode="anchor")
    cx.set_ylim(0, 6); cx.set_xlim(-0.4, 5.4)
    cx.set_xlabel("reaction left out")
    cx.set_ylabel(f"mean absolute error ({KC})")
    cx.text(0.99, 0.03, "n = 5 or 4 reactions", transform=cx.transAxes, ha="right", va="bottom", fontsize=8)
    panel(cx, "c", x=0.02, y=1.0)
    hd = legend_handles(DFT6)
    fig.legend(handles=hd, loc="lower center", ncol=6, frameon=False, bbox_to_anchor=(0.5, -0.005), handletextpad=0.2,
               columnspacing=0.7, fontsize=8)
    data = pd.DataFrame(rows)
    save_data(data, "fig03_summary_stats")
    save(fig, "fig03_summary_stats")
    print(data[data.panel == "a/b"].round(3).to_string())
    for m, v in zip(DFT5, [1.44, 1.08, 2.26, 3.29, 4.52]):  # Skala, wB, B3LYP, r2SCAN, PBE
        anchor("fig03", f"MAE {m}", v, s.loc[m].MAE, 0.006)
    print("  counts negative", cnt)
    # claims: ordering survives leave-one-out
    ok_order = True
    for o in order:
        q = L[L.removed == o].set_index("method").MAE
        ok_order &= bool(q["Skala-1.1"] < q["B3LYP-D3BJ"] < q["r2SCAN"] < q["PBE"])
    LOG.append(("fig03", "caption claim: Skala < B3LYP < r2SCAN < PBE in every subset", True, ok_order, ok_order))
    sw = [(o, bool(L[(L.removed == o) & (L.method == 'Skala-1.1')].MAE.iloc[0] < L[(L.removed == o) & (L.method == 'wB97M-V')].MAE.iloc[0])) for o in order]
    print("  Skala < wB97M-V per subset:", sw)
    note("fig03", "Skala-1.1 < wB97M-V holds only when Cl-+CH3Cl is removed (and never otherwise); caption claim that the order 'depends on HCN and on Cl-+CH3Cl' checked: "
         + str(sw))


# ------------------------------------------------------------------ Figure 5 / 6 helpers
def plot_runs(ax, q, xcol, ycol, **kw):
    """Line through points, broken where consecutive point ids are missing (no interpolation across a gap)."""
    q = q.sort_values("point_id")
    ids = [int(v[1:]) for v in q.point_id]
    seg = [[]]
    for i, (k, (_, r)) in zip(ids, enumerate(q.iterrows())):
        if seg[-1] and i - seg[-1][-1][0] > 1:
            seg.append([])
        seg[-1].append((i, r[xcol], r[ycol]))
    for sg in seg:
        if len(sg) > 1:
            ax.plot([a[1] for a in sg], [a[2] for a in sg], **kw)


def path_panels(ax_a, ax_b, P, rx, coord_label, ref_main="CCSD(T)+dMP2", flagged=None, noD3=True):
    p = P[P.reaction == rx]
    dm = p[p.reference == "CCSD(T)+dMP2"]
    raw = p[p.reference == "CCSD(T) raw"]
    rows = []
    refd = dm[dm.method == "Skala-1.1"].sort_values("point_id")
    refr = raw[raw.method == "Skala-1.1"].sort_values("point_id")
    if rx in ("hcn", "h2co"):
        # (a) reference curves
        ax_a.plot(refr.coordinate_deg_or_A, refr.ref_relative_kcal, color="black", ls=":", lw=1.0, marker="x", ms=5.5, mew=1.1, zorder=5)
        plot_runs(ax_a, refd, "coordinate_deg_or_A", "ref_relative_kcal", color="black", ls="-", lw=1.5, zorder=6)
        ax_a.plot(refd.coordinate_deg_or_A, refd.ref_relative_kcal, ls="none", marker="+", ms=8, mew=2.0, color="black", zorder=6)
        use = raw
        meths_a = ["Skala-1.1", "wB97M-V", "B3LYP-D3BJ"]
    else:
        ax_a.plot(refd.coordinate_deg_or_A, refd.ref_relative_kcal, color="black", ls="-", lw=1.5, marker="+", ms=8, mew=2.0, zorder=6)
        use = dm
        meths_a = ["Skala-1.1", "wB97M-V", "B3LYP-D3BJ"]
    for m in meths_a:
        q = use[use.method == m].sort_values("point_id")
        y = q.ref_relative_kcal + q.signed_error_kcal
        ax_a.plot(q.coordinate_deg_or_A, y, color=METH[m]["c"], lw=METH[m]["lw"], zorder=4, ls=METH[m]["ls"])
        ax_a.plot(q.coordinate_deg_or_A, y, ls="none", zorder=5, **mstyle(m))
        for c, v in zip(q.coordinate_deg_or_A, y):
            rows.append(dict(panel="a", method=m, coordinate=c, energy_rel_kcal=v))
    for c, v in zip(refd.coordinate_deg_or_A, refd.ref_relative_kcal):
        rows.append(dict(panel="a", method="CCSD(T)+dMP2", coordinate=c, energy_rel_kcal=v))
    for c, v in zip(refr.coordinate_deg_or_A, refr.ref_relative_kcal):
        rows.append(dict(panel="a", method="CCSD(T) raw", coordinate=c, energy_rel_kcal=v))
    # (b)
    band(ax_b); zero_line(ax_b)
    meths_b = ["Skala-1.1", "Skala-1.1-noD3", "wB97M-V", "B3LYP-D3BJ"] if noD3 else ["Skala-1.1", "wB97M-V", "B3LYP-D3BJ"]
    for m in meths_b:
        q = dm[dm.method == m].sort_values("point_id")
        st = METH[m]
        if m == "Skala-1.1-noD3" and rx == "nh3":
            ax_b.plot(q.coordinate_deg_or_A, q.signed_error_kcal, ls="none", zorder=7, **mstyle(m, 1.6, mfc="none"))
        else:
            plot_runs(ax_b, q, "coordinate_deg_or_A", "signed_error_kcal", color=st["c"], lw=st["lw"], ls=st["ls"], zorder=4)
            ax_b.plot(q.coordinate_deg_or_A, q.signed_error_kcal, ls="none", zorder=5, **mstyle(m))
        for c, v in zip(q.coordinate_deg_or_A, q.signed_error_kcal):
            rows.append(dict(panel="b", method=m, reference="CCSD(T)+dMP2", coordinate=c, signed_error_kcal=v))
    if rx in ("hcn", "h2co"):
        only_raw = sorted(set(raw.point_id) - set(dm.point_id))
        for m in ["Skala-1.1", "wB97M-V", "B3LYP-D3BJ"]:
            q = raw[(raw.method == m) & (raw.point_id.isin(only_raw))]
            st = METH[m]
            ax_b.plot(q.coordinate_deg_or_A, q.signed_error_kcal, ls="none", marker=st["m"], ms=MS[st["m"]] * 0.8, mfc="white",
                      mec=st["c"], mew=1.0, zorder=5)
            for c, v, pid in zip(q.coordinate_deg_or_A, q.signed_error_kcal, q.point_id):
                rows.append(dict(panel="b", method=m, reference="CCSD(T) raw", point_id=pid, coordinate=c, signed_error_kcal=v))
        return pd.DataFrame(rows), only_raw
    return pd.DataFrame(rows), []


def fig05(D):
    print("== Figure 5")
    P = D["path"]
    fig, (a, b) = plt.subplots(2, 1, figsize=(170 * MM, 110 * MM), sharex=True, gridspec_kw=dict(hspace=0.08))
    data, only_raw = path_panels(a, b, P, "hcn", "")
    print("  raw-only points:", only_raw)
    for ax in (a, b):
        ax.set_xlim(184, 0)
        ax.axvline(71.518, color="black", lw=0.6, ls="--", zorder=1)
        ax.axvspan(88, 130, color="black", alpha=0.045, zorder=1, lw=0)
    a.text(109, 1.0, "reference\nflagged", ha="center", va="bottom", fontsize=8, linespacing=0.95)
    a.set_ylim(-3, 58)
    b.set_ylim(-5, 2)
    a.set_ylabel(f"energy relative to\nHCN ({KC})")
    b.set_ylabel(f"error relative to\nCCSD(T) ({KC})")
    b.set_xlabel("H-C-N angle (degrees)")
    for x, t in [(180, "HCN"), (71.518, "TS"), (6, "HNC")]:
        a.text(x, 1.02, t, transform=a.get_xaxis_transform(), ha="center", va="bottom", fontsize=8)
    panel(a, "a", x=0.015, y=0.97); panel(b, "b", x=0.015, y=0.97)
    b.text(0.99, 0.03, "n = 3 methods + variant, 8 points (filled)", transform=b.transAxes, ha="right", va="bottom", fontsize=8, color="white")  # placeholder removed below
    b.texts[-1].remove()
    hd = [handle("Skala-1.1"), handle("Skala-1.1-noD3"), handle("wB97M-V"), handle("B3LYP-D3BJ"),
          Line2D([0], [0], color="black", lw=1.5, marker="+", ms=8, mew=2.0, label="CCSD(T) + MP2 basis correction"),
          Line2D([0], [0], color="black", lw=1.0, ls=":", marker="x", ms=5.5, mew=1.1, label="CCSD(T), def2-TZVP (flagged points open in b)")]
    fig.legend(handles=hd, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0), handletextpad=0.3, columnspacing=1.2)
    fig.subplots_adjust(left=0.10, right=0.985, top=0.93, bottom=0.215)
    save_data(data, "fig05_hcn_path")
    save(fig, "fig05_hcn_path")
    pa = data[data.panel == "b"]
    for m in ["Skala-1.1", "wB97M-V", "B3LYP-D3BJ"]:
        q = pa[pa.method == m]
        print(f"  {m} error range (all b points): {q.signed_error_kcal.min():.2f} .. {q.signed_error_kcal.max():.2f}")
    q = pa[(pa.method == 'wB97M-V') & (pa.reference == 'CCSD(T) raw')].sort_values('coordinate', ascending=False)
    print("  wB97M-V vs raw by angle:\n", q[['coordinate', 'signed_error_kcal']].to_string())
    dd = P[(P.reaction == "hcn") & (P.reference == "CCSD(T) raw") & (P.method.isin(["wB97M-V", "Skala-1.1", "B3LYP-D3BJ"]))]
    mx = {m: dd[dd.method == m].signed_error_kcal.abs().max() for m in ["Skala-1.1", "B3LYP-D3BJ"]}
    LOG.append(("fig05", "caption/plan claim: Skala and B3LYP within about 1 kcal/mol of raw CCSD(T) along path (max |err|)", "~1", f"Skala {mx['Skala-1.1']:.2f}, B3LYP {mx['B3LYP-D3BJ']:.2f}", True))


def fig06(D):
    print("== Figure 6")
    P = D["path"]
    fig, (a, b) = plt.subplots(2, 1, figsize=(85 * MM, 100 * MM), sharex=True, gridspec_kw=dict(hspace=0.08))
    data, _ = path_panels(a, b, P, "nh3", "")
    for ax in (a, b):
        ax.set_xlim(-42, 42)
        ax.axvline(0, color="black", lw=0.6, ls="--", zorder=1)
    a.set_ylim(-0.4, 7.2)
    b.set_ylim(-1.2, 0.8)
    a.set_ylabel(f"energy relative\nto NH₃ minimum\n({KC})")
    b.set_ylabel(f"error relative\nto CCSD(T)\n({KC})")
    b.set_xlabel("umbrella dihedral angle (degrees)")
    panel(a, "a", x=0.02, y=0.97); panel(b, "b", x=0.02, y=0.97)
    hd = [handle("Skala-1.1"), handle("Skala-1.1-noD3", line=False), handle("wB97M-V"), handle("B3LYP-D3BJ"),
          Line2D([0], [0], color="black", lw=1.5, marker="+", ms=8, mew=2.0, label="CCSD(T) + MP2 corr.")]
    fig.legend(handles=hd, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.0), handletextpad=0.3, columnspacing=0.8, fontsize=8)
    fig.subplots_adjust(left=0.27, right=0.985, top=0.98, bottom=0.27)
    a.yaxis.set_label_coords(-0.19, 0.5); b.yaxis.set_label_coords(-0.19, 0.5)
    save_data(data, "fig06_nh3_path")
    save(fig, "fig06_nh3_path")
    q = data[data.panel == "b"]
    for m, exp in [("Skala-1.1", 0.48), ("wB97M-V", -0.84), ("B3LYP-D3BJ", -0.93)]:
        v = q[q.method == m].signed_error_kcal
        got = v.loc[v.abs().idxmax()]
        anchor("fig06", f"NH3 path extreme error {m} (plan text)", exp, got, 0.006)


# ------------------------------------------------------------------ Figure 7
def parabola(xs, ys):
    c = np.polyfit(xs, ys, 2)
    xv = -c[1] / (2 * c[0])
    return c, xv, np.polyval(c, xv)


def fig07(D):
    print("== Figure 7")
    S = D["sym"]
    mins = D["sn2_minima"]
    dif = D["sn2_difference"]
    kinds = ["own_B3LYP_TS", "symmetric_stretch_scan"]
    fig = plt.figure(figsize=(170 * MM, 130 * MM))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.55, 1.0], hspace=0.10, wspace=0.27, left=0.095, right=0.995, top=0.985, bottom=0.215)
    axes = {}
    rows = []
    cfg = {"sn2_clcl": dict(col=0, xlim=(2.00, 2.80), ylim=(0, 35), ts=2.3483, name="Cl⁻ + CH₃Cl", let=("a", "c"), ref=13.41),
           "sn2_ff": dict(col=1, xlim=(1.68, 2.08), ylim=(10, 27), ts=1.8563, name="F⁻ + CH₃F", let=("b", "d"), ref=13.38)}
    for rx, c in cfg.items():
        a = fig.add_subplot(gs[0, c["col"]]); d = fig.add_subplot(gs[1, c["col"]], sharex=a)
        axes[rx] = (a, d)
        s = S[(S.reaction == rx) & S.kind.isin(kinds) & (S.density_fitting == True)]
        for m in ["Skala-1.1", "wB97M-V"]:
            q = s[(s.method == m) & (s.converged == 1)].sort_values("cx_distance_A")
            st = METH[m]
            a.plot(q.cx_distance_A, q.rel_to_RC_B3LYPgeom_kcal, color=st["c"], lw=1.0, zorder=3)
            a.plot(q.cx_distance_A, q.rel_to_RC_B3LYPgeom_kcal, ls="none", zorder=4, **mstyle(m))
            for x, y in zip(q.cx_distance_A, q.rel_to_RC_B3LYPgeom_kcal):
                rows.append(dict(reaction=rx, method=m, cx_distance_A=x, rel_to_RC_kcal=y, converged=1))
            if m == "Skala-1.1":
                a.plot(q.cx_distance_A, q.rel_to_RC_noD3_kcal, color=st["c"], lw=1.0, ls="--", zorder=3)
                a.plot(q.cx_distance_A, q.rel_to_RC_noD3_kcal, ls="none", zorder=4, **mstyle("Skala-1.1-noD3", 0.9))
                for x, y in zip(q.cx_distance_A, q.rel_to_RC_noD3_kcal):
                    rows.append(dict(reaction=rx, method="Skala-1.1-noD3", cx_distance_A=x, rel_to_RC_kcal=y, converged=1))
            # parabola from the fit points (recomputed from the scan rows) -> compare to derived_sn2_minima
            mm = mins[(mins.reaction == rx) & (mins.method == m)].iloc[0]
            fp = [float(v) for v in mm.fit_points_A.split(";")]
            pts = s[(s.method == m) & (s.converged == 1)]
            ys = [float(pts[np.isclose(pts.cx_distance_A, f, atol=6e-4)].rel_to_RC_B3LYPgeom_kcal.iloc[0]) for f in fp]
            co, xv, yv = parabola(fp, ys)
            xx = np.linspace(min(fp) - 0.02, max(fp) + 0.02, 60)
            a.plot(xx, np.polyval(co, xx), color=st["c"], lw=0.8, alpha=0.4, zorder=2)
            a.plot([mm.vertex_cx_A, mm.vertex_cx_A], [mm.vertex_energy_kcal - (c["ylim"][1] - c["ylim"][0]) * 0.03,
                                                       mm.vertex_energy_kcal + (c["ylim"][1] - c["ylim"][0]) * 0.03],
                   color=st["c"], lw=1.0, zorder=5)
            a.plot([mm.vertex_cx_A], [mm.vertex_energy_kcal], marker="*", ms=11, mfc=st["c"], mec="black", mew=0.5, ls="none", zorder=6)
            rows.append(dict(reaction=rx, method=m, kind="parabola_vertex", cx_distance_A=mm.vertex_cx_A, rel_to_RC_kcal=mm.vertex_energy_kcal))
            anchor("fig07", f"{rx} {m} vertex distance (recomputed vs csv)", mm.vertex_cx_A, xv, 0.0006)
            anchor("fig07", f"{rx} {m} vertex energy (recomputed vs csv)", mm.vertex_energy_kcal, yv, 0.0015)
            lab = f"{mm.vertex_cx_A:.3f} Å, {mm.vertex_energy_kcal:.1f}"
            if rx == "sn2_clcl":
                dx, dy, ha = (0.0, 5.2, "center") if m == "Skala-1.1" else (-0.03, 5.4, "right") if False else (0.0, 0, "left")
                if m == "Skala-1.1":
                    a.annotate(lab, (mm.vertex_cx_A, mm.vertex_energy_kcal), xytext=(2.52, 4.3), fontsize=8, ha="left", va="center", color="black",
                               arrowprops=dict(arrowstyle="-", lw=0.5, color=st["c"]))
                else:
                    a.annotate(lab, (mm.vertex_cx_A, mm.vertex_energy_kcal), xytext=(2.03, 6.0), fontsize=8, ha="left", va="center", color="black",
                               arrowprops=dict(arrowstyle="-", lw=0.5, color=st["c"]))
            else:
                if m == "Skala-1.1":
                    a.annotate(lab, (mm.vertex_cx_A, mm.vertex_energy_kcal), xytext=(1.69, 19.5), fontsize=8, ha="left", va="center",
                               arrowprops=dict(arrowstyle="-", lw=0.5, color=st["c"]))
                else:
                    a.annotate(lab, (mm.vertex_cx_A, mm.vertex_energy_kcal), xytext=(1.88, 10.9), fontsize=8, ha="left", va="center",
                               arrowprops=dict(arrowstyle="-", lw=0.5, color=st["c"]))
        # unconverged point
        if rx == "sn2_clcl":
            u = s[(s.method == "Skala-1.1") & (s.converged == 0)]
            y = float(u.rel_to_RC_B3LYPgeom_kcal.iloc[0])
            a.plot([2.298], [y], marker="o", ms=7, mfc="white", mec=ORANGE, mew=1.1, ls="none", zorder=7)
            a.plot([2.298], [y], marker="x", ms=7, mec=ORANGE, mew=1.0, ls="none", zorder=8)
            rows.append(dict(reaction=rx, method="Skala-1.1", cx_distance_A=2.298, rel_to_RC_kcal=y, converged=0))
            anchor("fig07", "unconverged Skala Cl y at 2.298", 15.75, y, 0.006)
            # xTB TS
            xt = S[(S.reaction == rx) & (S.kind == "xtb_TS_geometry") & (S.method == "Skala-1.1")]
            yx = float(xt.rel_to_RC_B3LYPgeom_kcal.iloc[0])
            a.plot([2.201, 2.236], [yx, yx], color=GREY, lw=1.0, zorder=5)
            a.plot([2.201, 2.236], [yx, yx], ls="none", marker="|", ms=5, color=GREY, zorder=5)
            a.plot([2.219], [yx], marker="*", ms=10, mfc=GREY, mec="black", mew=0.4, ls="none", zorder=9)
            rows.append(dict(reaction=rx, method="GFN2-xTB TS geometry (Skala-1.1 energy)", cx_distance_A=2.219, rel_to_RC_kcal=yx))
            anchor("fig07", "xTB TS point energy", 16.96, yx, 0.006)
        # B3LYP TS line, ref lines
        for ax_ in (a, d):
            ax_.axvline(c["ts"], color=GREY, lw=0.8, ls=":", zorder=1)
        a.text(c["ts"] + 0.004 * (c["xlim"][1] - c["xlim"][0]) / 0.8, c["ylim"][1] - 0.5, "B3LYP TS", fontsize=8, color="#555555", ha="left", va="top", rotation=0)
        refc = float(s.ref_ccsdt_dmp2_kcal.iloc[0])
        a.axhline(refc, color="black", lw=1.0, ls="--", zorder=2)
        a.axhline(c["ref"], color=GREY, lw=0.6, zorder=1)
        a.set_xlim(*c["xlim"]); a.set_ylim(*c["ylim"])
        a.tick_params(labelbottom=False)
        a.set_ylabel(f"energy relative to the\nion-molecule complex ({KC})")
        panel(a, c["let"][0], x=0.02, y=0.985)
        # difference
        dd = dif[(dif.reaction == rx) & dif.skala_minus_wB97MV_kcal.notna()].sort_values("cx_distance_A")
        zero_line(d)
        d.plot(dd.cx_distance_A, dd.skala_minus_wB97MV_kcal, color="black", lw=1.0, zorder=3)
        d.plot(dd.cx_distance_A, dd.skala_minus_wB97MV_kcal, ls="none", marker="o", ms=5, mfc="black", mec="black", zorder=4)
        for x, y in zip(dd.cx_distance_A, dd.skala_minus_wB97MV_kcal):
            rows.append(dict(reaction=rx, method="Skala-1.1 minus wB97M-V", cx_distance_A=x, rel_to_RC_kcal=y))
        d.set_ylim(-7, 4)
        d.set_xlabel("C-X distance, both bonds (Å)")
        d.set_ylabel(f"Skala-1.1 minus\nωB97M-V ({KC})")
        panel(d, c["let"][1], x=0.02, y=0.97)
        if rx == "sn2_clcl":
            span = dd.skala_minus_wB97MV_kcal.max() - dd.skala_minus_wB97MV_kcal.min()
            dist = dd.cx_distance_A.max() - dd.cx_distance_A.min()
            print(f"  Cl difference changes by {span:.2f} over {dist:.2f} A")
            anchor("fig07", "Cl difference span (caption: 7.9)", 7.9, span, 0.06)
        else:
            print("  F difference range", dd.skala_minus_wB97MV_kcal.min(), dd.skala_minus_wB97MV_kcal.max())
    hd = [handle("Skala-1.1"), handle("Skala-1.1-noD3"), handle("wB97M-V"),
          Line2D([0], [0], marker="*", ls="none", mfc="white", mec="black", ms=10, label="parabola minimum"),
          Line2D([0], [0], marker="$\\otimes$", ls="none", color=ORANGE, ms=9, label="SCF not converged"),
          Line2D([0], [0], color="black", ls="--", lw=1.0, label="CCSD(T) + MP2 corr., B3LYP TS"),
          Line2D([0], [0], color=GREY, lw=0.6, label="literature barrier"),
          Line2D([0], [0], color=GREY, lw=0.8, ls=":", label="B3LYP TS distance"),
          Line2D([0], [0], marker="*", ls="-", color=GREY, mfc=GREY, mec="black", mew=0.4, ms=9, label="GFN2-xTB TS geometry")]
    fig.legend(handles=hd, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0), handletextpad=0.3, columnspacing=0.8)
    data = pd.DataFrame(rows)
    save_data(data, "fig07_sn2_stretch")
    save(fig, "fig07_sn2_stretch")
    # extra: F panel Skala above wB97M-V throughout
    f = dif[(dif.reaction == "sn2_ff") & dif.skala_minus_wB97MV_kcal.notna()]
    LOG.append(("fig07", "claim: F panel Skala above wB97M-V at every common point", True, bool((f.skala_minus_wB97MV_kcal > 0).all()), bool((f.skala_minus_wB97MV_kcal > 0).all())))
    mm = mins.set_index(["reaction", "method"])
    anchor("fig07", "Skala Cl min distance", 2.423, mm.loc[("sn2_clcl", "Skala-1.1")].vertex_cx_A, 0.0006)
    anchor("fig07", "Skala Cl min energy", 9.53, mm.loc[("sn2_clcl", "Skala-1.1")].vertex_energy_kcal, 0.006)
    anchor("fig07", "wB97M-V Cl min distance", 2.316, mm.loc[("sn2_clcl", "wB97M-V")].vertex_cx_A, 0.0006)
    anchor("fig07", "wB97M-V Cl min energy", 13.34, mm.loc[("sn2_clcl", "wB97M-V")].vertex_energy_kcal, 0.006)
