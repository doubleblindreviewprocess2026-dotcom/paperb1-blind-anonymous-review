"""Figure 1 (schematic), graphical abstract, supplementary figures S1-S9."""
from figcommon import *
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch
from matplotlib.lines import Line2D
from figs_fallback import path_panels, parabola
from figs_rest import err_row

DARK = "#444444"
LG = "#9A9A9A"


# ------------------------------------------------------------------ Figure 1
def fig01(D):
    print("== Figure 1")
    fig = plt.figure(figsize=(170 * MM, 85 * MM))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 170); ax.set_ylim(0, 85); ax.axis("off")
    # block titles
    ax.text(1, 81.8, "Reactions (reference, kcal mol⁻¹)", fontsize=9, fontweight="bold", va="center")
    ax.text(68, 81.8, "Geometries", fontsize=9, fontweight="bold", va="center")
    ax.text(111, 81.8, "Energies", fontsize=9, fontweight="bold", va="center")
    rows = [("HCN → HNC", "H shift, 48.07", "NHTBH38/08 database", False),
            ("NH₃ inversion", "inversion, 5.12", "Lin et al.", False),
            ("H₂CO → H₂ + CO", "elimination, 81.9", "Feller et al., zero-point incl.", False),
            ("F⁻ + CH₃F", "SN2, 13.38 from complex", "NHTBH38/08 database", False),
            ("Cl⁻ + CH₃Cl", "SN2, 13.41 from complex", "NHTBH38/08 database", False),
            ("CH₃NC → CH₃CN", "literature reference unverified,", "internal reference only", True)]
    yc = [71.0, 59.0, 47.0, 35.0, 23.0, 10.5]
    for (n, l2, l3, g), y in zip(rows, yc):
        col = LG if g else "black"
        ec = LG if g else DARK
        ax.add_patch(Rectangle((1, y - 5.3), 64, 10.6, fc="white", ec=ec, lw=0.6, ls="--" if g else "-"))
        ax.text(2.5, y + 3.3, n, fontsize=8, fontweight="bold", color=col, va="center")
        ax.text(2.5, y - 0.1, l2, fontsize=8, color=col, va="center")
        ax.text(2.5, y - 3.5, l3, fontsize=8, color=col, va="center")
    # geometries
    boxes = [(62.0, 79.0, "B3LYP-D3(BJ)/def2-TZVP\n(def2-TZVPD for SN2),\nminima and transition\nstates, Hessian checked"),
             (43.0, 55.0, "GFN2-xTB geometries\n(scout test)"),
             (24.5, 37.5, "relaxed scans: HCN 11,\nNH₃ 9, H₂CO 11 points"),
             (5.0, 19.0, "SN2 symmetric stretch,\nrigid CH₃, density\nfitting")]
    for y0, y1, t in boxes:
        ax.add_patch(Rectangle((69, y0), 38, y1 - y0, fc="white", ec=DARK, lw=0.6))
        ax.text(88, (y0 + y1) / 2, t, fontsize=8, ha="center", va="center", linespacing=1.05)
    ax.plot([65.5, 67.0, 67.0], [70.5, 70.5, 12.0], color=DARK, lw=0.6)
    ax.add_patch(FancyArrowPatch((67.0, 70.5), (68.8, 70.5), arrowstyle="-|>", mutation_scale=6, lw=0.6, color=DARK, shrinkA=0, shrinkB=0))
    for y0, y1, _ in boxes[1:]:
        ax.add_patch(FancyArrowPatch((67.0, (y0 + y1) / 2), (68.8, (y0 + y1) / 2), arrowstyle="-|>", mutation_scale=6, lw=0.6, color=DARK, shrinkA=0, shrinkB=0))
    ax.add_patch(FancyArrowPatch((107.3, 47), (110.3, 47), arrowstyle="-|>", mutation_scale=7, lw=0.8, color=DARK, shrinkA=0, shrinkB=0))
    # energies matrix
    cols = ["stationary points, TZVP", "QZVPP (HCN, NH₃, H₂CO)", "paths", "xTB geometries", "SN2 stretch"]
    meths = ["Skala-1.1", "Skala-1.1-noD3", "wB97M-V", "B3LYP-D3BJ", "r2SCAN", "PBE", "CCSD(T)+dMP2"]
    lab = ["Skala-1.1", "Skala-1.1 no D3", "ωB97M-V", "B3LYP-D3(BJ)", "r²SCAN", "PBE", "CCSD(T)+MP2 corr."]
    # 1 = computed, 0.5 = partial, 0 = none
    fill = {"Skala-1.1": [1, 1, 1, 1, 1], "Skala-1.1-noD3": [1, 1, 1, 1, 1], "wB97M-V": [1, 1, 1, 1, 1],
            "B3LYP-D3BJ": [1, 0, 1, 1, 0], "r2SCAN": [1, 0, 0, 0, 0], "PBE": [1, 0, 0, 0, 0], "CCSD(T)+dMP2": [1, 0, 0.5, 0.5, 0]}
    x0, cw, rh, ytop = 143.5, 5.3, 5.0, 44.5
    for j, c in enumerate(cols):
        ax.text(x0 + j * cw + cw / 2, 80.5, c, rotation=90, fontsize=8, ha="center", va="top")
    for i, (m, l) in enumerate(zip(meths, lab)):
        y = ytop - i * rh - rh / 2
        st = mstyle(m, 0.9)
        ax.plot([112.5], [y], ls="none", **st)
        ax.text(114.6, y, l, fontsize=8, va="center")
        for j, f in enumerate(fill[m]):
            r = Rectangle((x0 + j * cw + 0.5, y - rh / 2 + 0.5), cw - 1.0, rh - 1.0, ec=DARK, lw=0.5,
                          fc=DARK if f == 1 else "white", hatch="////" if f == 0.5 else None)
            ax.add_patch(r)
    # legend of cell fills
    ly = ytop - 7 * rh - 2.6
    ax.add_patch(Rectangle((108, ly - 1.4), 3, 2.8, ec=DARK, lw=0.5, fc=DARK)); ax.text(111.8, ly, "computed", fontsize=8, va="center")
    ax.add_patch(Rectangle((127, ly - 1.4), 3, 2.8, ec=DARK, lw=0.5, fc="white", hatch="////")); ax.text(130.8, ly, "part of the set", fontsize=8, va="center")
    ax.add_patch(Rectangle((153, ly - 1.4), 3, 2.8, ec=DARK, lw=0.5, fc="white")); ax.text(156.8, ly, "not run", fontsize=8, va="center")
    save_data(pd.DataFrame([dict(method=m, **{c: fill[m][j] for j, c in enumerate(cols)}) for m in meths]), "fig01_design_matrix")
    save(fig, "fig01_study_design")
    note("fig01", "Matrix cells for CCSD(T)+MP2 on paths and xTB geometries are hatched (partial: paths exclude 7 flagged points; xTB geometries only HCN, NH3, CH3NC), following dossier A.1/A.5; the plan lists them as computed.")


# ------------------------------------------------------------------ Graphical abstract
def graphical_abstract(D):
    print("== Graphical abstract")
    b = D["bar"]; h = b[b.in_headline == True]
    S = D["sym"]; mins = D["sn2_minima"]
    fig = plt.figure(figsize=(6, 3))
    plt.rcParams.update({"font.size": 10})
    fig.text(0.5, 0.955, "Skala-1.1 barriers: no systematic underestimation", ha="center", va="center", fontsize=11, fontweight="bold")
    a = fig.add_axes([0.105, 0.20, 0.54, 0.64])
    band(a); zero_line(a)
    meths = ["Skala-1.1", "wB97M-V", "PBE"]
    for j, m in enumerate(meths):
        off = (j - 1) * 0.24
        for k, r in enumerate(RX):
            y = float(h[(h.reaction == r) & (h.method == m)].signed_error_main_kcal.iloc[0])
            a.plot([k + off] * 2, [0, y], color=METH[m]["c"], lw=1.2, zorder=3)
            a.plot([k + off], [y], ls="none", zorder=4, **mstyle(m, 1.2))
    a.set_xlim(-0.5, 4.5); a.set_ylim(-9, 5.2)
    a.set_xticks(range(5)); a.set_xticklabels(["HCN", "NH₃", "H₂CO", "F⁻+\nCH₃F", "Cl⁻+\nCH₃Cl"], fontsize=10)
    a.set_yticks([-8, -4, 0, 4]); a.tick_params(axis="y", labelsize=10)
    a.set_ylabel("signed error (kcal mol⁻¹)", fontsize=10)
    a.legend(handles=[handle(m, line=False, ms_scale=1.2) for m in meths], loc="upper left", ncol=3, frameon=False, fontsize=10,
             handletextpad=0.1, columnspacing=0.8, bbox_to_anchor=(-0.02, 1.02))
    # miniature of Fig 7a
    c = fig.add_axes([0.72, 0.20, 0.26, 0.55])
    rx = "sn2_clcl"
    s = S[(S.reaction == rx) & S.kind.isin(["own_B3LYP_TS", "symmetric_stretch_scan"]) & (S.density_fitting == True)]
    for m in ["Skala-1.1", "wB97M-V"]:
        q = s[(s.method == m) & (s.converged == 1)].sort_values("cx_distance_A")
        c.plot(q.cx_distance_A, q.rel_to_RC_B3LYPgeom_kcal, color=METH[m]["c"], lw=1.2, zorder=3)
        c.plot(q.cx_distance_A, q.rel_to_RC_B3LYPgeom_kcal, ls="none", zorder=4, **mstyle(m, 0.8))
        mm = mins[(mins.reaction == rx) & (mins.method == m)].iloc[0]
        c.plot([mm.vertex_cx_A], [mm.vertex_energy_kcal], marker="*", ms=13, mfc=METH[m]["c"], mec="black", mew=0.5, ls="none", zorder=6)
    c.axhline(13.61, color="black", lw=1.0, ls="--", zorder=2)
    c.set_xlim(2.0, 2.8); c.set_ylim(0, 33)
    c.set_xticks([]); c.set_yticks([])
    c.spines["left"].set_visible(True); c.spines["bottom"].set_visible(True)
    fig.text(0.855, 0.115, "one located failure", ha="center", va="center", fontsize=10, fontweight="bold", color=ORANGE)
    for f in ("png", "svg"):
        fig.savefig(OUT + f"graphical_abstract.{f}", dpi=200)
    plt.close(fig)
    plt.rcParams.update({"font.size": 8})
    from PIL import Image
    im = Image.open(OUT + "graphical_abstract.png"); print("  size px", im.size)
    note("graphical_abstract", "Built at 6 x 3 in (1200 x 600 px at 200 dpi) so that 10 pt labels fit; the plan's '5 x 2.5 cm' cannot carry 10 pt labels. Axes of the miniature carry no tick labels (plan: no other text); left panel has a y label and short reaction labels.")


# ------------------------------------------------------------------ S1
def s1(D):
    print("== S1")
    ir = D["internal_reference"]
    q = ir[(ir.own_reference_usable == True) & (ir.method != "Skala-1.1-noD3")]
    fig, ax = fig_ax(85, 85)
    ax.plot([-9, 5], [-9, 5], color="black", lw=0.6, ls="-", zorder=1)
    rows = []
    for m in DFT5:
        for _, r in q[q.method == m].iterrows():
            st = mstyle(m, 0.9)
            if r.reaction == "ch3nc":
                st.update(mfc="white", mec=METH[m]["c"], mew=1.1)
            ax.plot([r.signed_error_vs_literature_kcal], [r.signed_error_vs_own_ccsdt_kcal], ls="none", zorder=3, **st)
            rows.append(dict(reaction=r.reaction, method=m, x_vs_literature=r.signed_error_vs_literature_kcal, y_vs_own_ccsdt=r.signed_error_vs_own_ccsdt_kcal))
    init = {"hcn": "H", "nh3": "N", "ch3nc": "M", "sn2_ff": "F", "sn2_clcl": "C"}
    for rr, ll in init.items():
        z = q[(q.reaction == rr) & (q.method == "Skala-1.1")].iloc[0]
        ax.text(z.signed_error_vs_literature_kcal + 0.3, z.signed_error_vs_own_ccsdt_kcal + 0.15, ll, fontsize=8, fontweight="bold", va="bottom", ha="left")
    ax.set_xlim(-7, 4); ax.set_ylim(-7, 4); ax.set_aspect("equal")
    ax.set_xlabel(f"error against literature ({KC})")
    ax.set_ylabel(f"error against own CCSD(T) + MP2\ncorrection ({KC})")
    ax.text(0.98, 0.02, "H HCN, N NH₃, M CH₃NC,\nF F⁻ + CH₃F, C Cl⁻ + CH₃Cl\n(labels on Skala-1.1)\nopen symbols (M): literature\nvalue unverified", transform=ax.transAxes, fontsize=8, ha="right", va="bottom", linespacing=1.0)
    hd = legend_handles(DFT5, line=False)
    fig.legend(handles=hd, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0), handletextpad=0.2, columnspacing=0.6, fontsize=8)
    fig.subplots_adjust(left=0.17, right=0.98, top=0.98, bottom=0.30)
    save_data(pd.DataFrame(rows), "figS1_two_references")
    save(fig, "figS1_two_references")
    note("figS1", "Reaction initials are printed only next to the Skala-1.1 points (labelling all 25 points overlaps along the diagonal); H2CO excluded because no usable own reference.")


# ------------------------------------------------------------------ S2
def internal_stats(D):
    ir = D["internal_reference"]
    q = ir[ir.own_reference_usable == True]
    out = []
    for m in DFT6:
        e = q[q.method == m].signed_error_vs_own_ccsdt_kcal
        out.append(dict(method=m, n=len(e), MAE=e.abs().mean(), MSE=e.mean(), MaxAE=e.abs().max(), n_negative=int((e < 0).sum())))
    return pd.DataFrame(out)


def s2(D):
    print("== S2")
    st = internal_stats(D).set_index("method")
    fig = plt.figure(figsize=(170 * MM, 62 * MM))
    gs = fig.add_gridspec(1, 2, left=0.17, right=0.995, top=0.97, bottom=0.22, wspace=0.12, width_ratios=[1, 1])
    a = fig.add_subplot(gs[0]); b = fig.add_subplot(gs[1], sharey=a)
    ypos = np.arange(6)[::-1]
    for y, m in zip(ypos, DFT6):
        s_ = METH[m]
        kw = dict(height=0.62, edgecolor="black", linewidth=0.4, zorder=3)
        if s_["open"]:
            kw.update(facecolor="white", edgecolor=s_["c"], hatch="////", linewidth=0.8)
        else:
            kw.update(facecolor=s_["c"])
        a.barh(y, st.loc[m].MAE, **kw); b.barh(y, st.loc[m].MSE, **kw)
        a.text(st.loc[m].MAE + 0.1, y, f"{st.loc[m].MAE:.2f}", va="center", fontsize=8)
        b.text(max(st.loc[m].MSE, 0) + 0.1, y, f"{st.loc[m].MSE:+.2f}".replace("-", "−"), va="center", fontsize=8)
    a.set_yticks(ypos); a.set_yticklabels(["Skala-1.1", "Skala-1.1\nwithout D3", "ωB97M-V", "B3LYP-D3(BJ)", "r²SCAN", "PBE"])
    b.tick_params(labelleft=False, left=False); b.spines["left"].set_visible(False)
    b.axvline(0, color="black", lw=0.8)
    a.set_xlim(0, 5.2); b.set_xlim(-4.3, 1.5)
    a.set_xlabel(f"mean absolute error ({KC})"); b.set_xlabel(f"mean signed error ({KC})")
    a.text(0.99, 0.93, "n = 5 reactions\nCH₃NC included,\nH₂CO excluded", transform=a.transAxes, ha="right", va="top", fontsize=8, linespacing=1.0)
    a.set_ylim(-0.6, 6.0)
    panel(a, "a", x=0.02, y=0.99) ; panel(b, "b", x=0.02, y=0.99)
    a.texts[-1].set_position((0.02, 0.99))
    save_data(st.reset_index(), "figS2_internal_reference_stats")
    save(fig, "figS2_internal_reference_stats")
    print(st.round(3))
    anchor("figS2", "Skala MAE internal ref (plan 1.29)", 1.29, st.loc["Skala-1.1"].MAE, 0.006)
    anchor("figS2", "wB97M-V MAE internal ref (plan 1.29)", 1.29, st.loc["wB97M-V"].MAE, 0.006)


# ------------------------------------------------------------------ S3
def s3(D):
    print("== S3")
    P = D["path"]
    fig, (a, b) = plt.subplots(2, 1, figsize=(170 * MM, 110 * MM), sharex=True, gridspec_kw=dict(hspace=0.08))
    data, only_raw = path_panels(a, b, P, "h2co", "")
    print("  raw-only points", only_raw)
    for ax in (a, b):
        ax.set_xlim(1.92, 0.73)
        ax.axvspan(1.115, 1.43, color="black", alpha=0.045, zorder=1, lw=0)
        ax.axvline(1.3114, color="black", lw=0.6, ls="--", zorder=1)
    a.text(1.27, 2.0, "no reliable\nreference", ha="center", va="bottom", fontsize=8, linespacing=0.95)
    a.text(1.3114, 1.01, "TS", transform=a.get_xaxis_transform(), ha="center", va="bottom", fontsize=8)
    b.set_xlabel("H···H distance (Å)")
    a.set_ylabel(f"energy relative to\nH₂CO ({KC})"); b.set_ylabel(f"error relative to\nCCSD(T) ({KC})")
    a.set_ylim(-4, 100)
    panel(a, "a", x=0.015, y=0.97); panel(b, "b", x=0.015, y=0.97)
    lo, hi = b.get_ylim(); b.set_ylim(min(lo, -4.5), max(hi, 6.5))
    hd = [handle("Skala-1.1"), handle("Skala-1.1-noD3"), handle("wB97M-V"), handle("B3LYP-D3BJ"),
          Line2D([0], [0], color="black", lw=1.5, marker="+", ms=8, mew=2.0, label="CCSD(T) + MP2 basis correction"),
          Line2D([0], [0], color="black", lw=1.0, ls=":", marker="x", ms=5.5, mew=1.1, label="CCSD(T), def2-TZVP (flagged points open in b)")]
    fig.legend(handles=hd, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0), handletextpad=0.3, columnspacing=1.2)
    fig.subplots_adjust(left=0.10, right=0.985, top=0.95, bottom=0.215)
    save_data(data, "figS3_h2co_path")
    save(fig, "figS3_h2co_path")
    r = P[(P.reaction == "h2co") & (P.reference == "CCSD(T) raw") & (P.method == "Skala-1.1")].sort_values("point_id")
    print(r[["point_id", "coordinate_deg_or_A", "ref_relative_kcal"]].round(2).T)


# ------------------------------------------------------------------ S4
def s4(D):
    print("== S4")
    P = D["path"]
    q = P[(P.reaction == "hcn") & (P.reference == "CCSD(T) raw")].pivot(index=["point_id", "coordinate_deg_or_A"], columns="method", values="signed_error_kcal").reset_index()
    q["wB97MV_minus_B3LYP"] = q["wB97M-V"] - q["B3LYP-D3BJ"]
    q["wB97MV_minus_Skala"] = q["wB97M-V"] - q["Skala-1.1"]
    fig, ax = fig_ax(85, 60)
    zero_line(ax)
    ax.plot(q.coordinate_deg_or_A, q.wB97MV_minus_B3LYP, color=GREEN, lw=1.2, marker="^", ms=5.6, mec="black", mew=0.4, zorder=3, label="ωB97M-V minus B3LYP-D3(BJ)")
    ax.plot(q.coordinate_deg_or_A, q.wB97MV_minus_Skala, color=ORANGE, lw=1.8, marker="o", ms=5.2, mec="black", mew=0.4, zorder=3, label="ωB97M-V minus Skala-1.1")
    ax.set_xlim(184, 0)
    ax.set_xlabel("H-C-N angle (degrees)"); ax.set_ylabel(f"difference of relative\nenergies ({KC})")
    ax.legend(frameon=False, loc="lower left", fontsize=8, bbox_to_anchor=(0.0, 0.0), handletextpad=0.3)
    ax.set_ylim(-7.5, 3)
    fig.subplots_adjust(left=0.2, right=0.98, top=0.98, bottom=0.2)
    save_data(q, "figS4_wb97mv_deviation")
    save(fig, "figS4_wb97mv_deviation")
    print(q.round(2).to_string())


# ------------------------------------------------------------------ S5
def s5(D):
    print("== S5")
    d = np.load(REC + "phase05_hcn_grid.npz", allow_pickle=True)
    th, r, t1, E = d["theta"], d["r"], d["t1"], d["E"]
    Er = (E - E.min()) * 627.5094740631
    win = Er < 120
    n_win = int(win.sum()); n_ok = int((win & (t1 < 0.02)).sum())
    anchor("figS5", "cells in window (caption: 369)", 369, n_win, 0.1)
    anchor("figS5", "cells in window with T1<0.02 (caption: 111)", 111, n_ok, 0.1)
    g = pd.read_csv(RES + "phase05_hcn_grid.csv")
    LOG.append(("figS5", "axis meaning: phase05_hcn_grid.csv coordinate = 'theta,r' and npz theta (25), r (17)", "theta, r", f"first coordinate {g.coordinate.iloc[0]}, theta[0]={th[0]}, r[0]={r[0]}", True))
    fig, ax = fig_ax(170, 80)
    dth, dr = th[1] - th[0], r[1] - r[0]
    X, Y = np.meshgrid(np.append(th - dth / 2, th[-1] + dth / 2), np.append(r - dr / 2, r[-1] + dr / 2))
    T = np.where(win, t1, np.nan).T
    pm = ax.pcolormesh(X, Y, T, cmap="cividis", vmin=0.01, vmax=0.06, shading="flat")
    Xc, Yc = np.meshgrid(th, r)
    ax.pcolormesh(X, Y, np.where(win, np.nan, 1.0).T, cmap="Greys", vmin=0, vmax=3, shading="flat")
    cs = ax.contour(Xc, Yc, t1.T, levels=[0.02, 0.03], colors=["white", "black"], linewidths=[1.0, 1.0], linestyles=["-", "--"])
    ax.contour(Xc, Yc, Er.T, levels=[120], colors=["#CC79A7"], linewidths=1.4)
    cb = fig.colorbar(pm, ax=ax, pad=0.015, extend="max")
    cb.set_label("T1 diagnostic"); cb.ax.tick_params(labelsize=8)
    ax.set_xlabel("θ (degrees; 0 = H on the C side, 180 = H on the N side)")
    ax.set_ylabel("r (Å, H to the C-N centre of mass)")
    ax.set_xlim(-dth / 2, 180 + dth / 2); ax.set_ylim(r[0] - dr / 2, r[-1] + dr / 2)
    hd = [Line2D([0], [0], color="white", lw=1.0, label="T1 = 0.02", path_effects=None),
          Line2D([0], [0], color="black", lw=1.0, ls="--", label="T1 = 0.03"),
          Line2D([0], [0], color="#CC79A7", lw=1.4, label="120 kcal mol⁻¹ above the minimum"),
          Patch(fc="#E0E0E0", ec="none", label="outside the window")]
    leg = ax.legend(handles=hd, loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=4, frameon=True, facecolor="#777777", edgecolor="none", fontsize=8,
                    labelcolor="white", handletextpad=0.4, columnspacing=1.2)
    fig.subplots_adjust(left=0.075, right=0.97, top=0.98, bottom=0.27)
    save_data(pd.DataFrame(dict(theta=np.repeat(th, len(r)), r=np.tile(r, len(th)), t1=t1.ravel(), E_rel_kcal=Er.ravel(), in_window=win.ravel())), "figS5_hcn_t1_map")
    save(fig, "figS5_hcn_t1_map")
    print("  window cells", n_win, "T1<0.02:", n_ok, "total cells", t1.size)


# ------------------------------------------------------------------ S6
def ranks(D):
    b = D["bar"]; h = b[(b.in_headline == True) & b.method.isin(DFT5)]
    tab = h.pivot(index="reaction", columns="method", values="signed_error_main_kcal").loc[RX, DFT5].abs()
    tie, strict = {}, {}
    for r in RX:
        v = tab.loc[r].sort_values()
        rem = list(v.index); k = 1; rk = {}
        while rem:
            best = v[rem[0]]
            grp = [m for m in rem if v[m] - best <= 0.5]
            for m in grp:
                rk[m] = k
            rem = [m for m in rem if m not in grp]; k += 1
        tie[r] = rk
        strict[r] = {m: int(i + 1) for i, m in enumerate(v.index)}
    return tab, pd.DataFrame(tie).T.loc[RX, DFT5], pd.DataFrame(strict).T.loc[RX, DFT5]


def s6(D):
    print("== S6")
    tab, tie, strict = ranks(D)
    fig, axs = plt.subplots(1, 2, figsize=(170 * MM, 62 * MM), sharey=True)
    for ax, t, let, ttl in zip(axs, [tie, strict], "ab", ["ties within 0.5 kcal mol⁻¹", "strict"]):
        ax.imshow(t.values, cmap="Greys_r", vmin=0, vmax=7.5, aspect="auto")
        for i in range(5):
            for j in range(5):
                v = t.values[i, j]
                ax.text(j, i, str(v), ha="center", va="center", fontsize=9, fontweight="bold", color="white" if v <= 3 else "black")
        ax.set_xticks(range(5)); ax.set_xticklabels(["Skala-\n1.1", "ωB97M-\nV", "B3LYP-\nD3(BJ)", "r²SCAN", "PBE"], fontsize=8)
        ax.xaxis.tick_top()
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.tick_params(length=0)
        ax.set_xlabel(ttl, fontsize=9, labelpad=4)
    axs[0].set_yticks(range(5)); axs[0].set_yticklabels([RXL[r].replace(" → ", " →\n") if r == "h2co" else RXL[r] for r in RX], fontsize=8)
    fig.subplots_adjust(left=0.19, right=0.99, top=0.74, bottom=0.12, wspace=0.07)
    for ax, let in zip(axs, "ab"):
        fig.text(ax.get_position().x0 - 0.005, 0.97, f"({let})", fontsize=9, fontweight="bold", va="top")
    save_data(pd.concat([tie.add_prefix("tie_"), strict.add_prefix("strict_")], axis=1).reset_index(), "figS6_ranks")
    save(fig, "figS6_ranks")
    print(tie, strict, sep="\n")
    sk = tie["Skala-1.1"].tolist(); ss = strict["Skala-1.1"].tolist()
    LOG.append(("figS6", "Skala ranks with ties [1,2,2,2,2]", "[1, 2, 2, 2, 2]", str(sk), sk == [1, 2, 2, 2, 2]))
    LOG.append(("figS6", "Skala ranks strict [1,5,2,2,2]", "[1, 5, 2, 2, 2]", str(ss), ss == [1, 5, 2, 2, 2]))
    exp = {"hcn": [1, 4, 1, 2, 3], "nh3": [2, 1, 1, 1, 1], "h2co": [2, 1, 3, 4, 5], "sn2_ff": [2, 1, 3, 4, 5], "sn2_clcl": [2, 1, 3, 4, 5]}
    ok = all(tie.loc[r].tolist() == v for r, v in exp.items())
    LOG.append(("figS6", "all 25 tie ranks equal dossier B.5 table", True, ok, ok))


# ------------------------------------------------------------------ S7
def s7(D):
    print("== S7")
    e = D["err"]
    q = e[(e.basis_kind == "main") & e.lit_id.isin(["L2", "L3", "L7", "L9"]) & e.method.isin(DFT6)]
    items = [("L2", "HCN → HNC\nreverse"), ("L3", "HCN → HNC\nreaction energy"), ("L7", "F⁻ + CH₃F\nfrom reactants"), ("L9", "Cl⁻ + CH₃Cl\nfrom reactants")]
    fig, ax = fig_ax(170, 80)
    band(ax); zero_line(ax)
    n = 6; step = 0.13
    rows = []
    for i, m in enumerate(DFT6):
        off = (i - (n - 1) / 2) * step
        for k, (lid, _) in enumerate(items):
            y = float(q[(q.lit_id == lid) & (q.method == m)].signed_error_kcal.iloc[0])
            ax.plot([k + off] * 2, [0, y], color=METH[m]["c"], lw=0.8, zorder=3)
            ax.plot([k + off], [y], ls="none", zorder=4, **mstyle(m))
            rows.append(dict(lit_id=lid, method=m, signed_error_kcal=y))
    ax.set_xticks(range(4)); ax.set_xticklabels([t for _, t in items])
    ax.set_xlim(-0.5, 3.5); ax.set_ylim(-7.5, 3.5)
    ax.set_xlabel("quantity"); ax.set_ylabel(f"signed error ({KC})")
    ax.text(0.99, 0.98, "n = 4 quantities (3 reactions)", transform=ax.transAxes, ha="right", va="top", fontsize=8)
    fig.legend(handles=[handle(m, line=False) for m in DFT6], loc="lower center", ncol=6, frameon=False, bbox_to_anchor=(0.5, 0.0), handletextpad=0.2, columnspacing=0.8, fontsize=8)
    fig.subplots_adjust(left=0.08, right=0.995, top=0.98, bottom=0.30)
    save_data(pd.DataFrame(rows), "figS7_other_quantities")
    save(fig, "figS7_other_quantities")
    print(pd.DataFrame(rows).pivot(index="method", columns="lit_id", values="signed_error_kcal").round(2))


# ------------------------------------------------------------------ S8
def s8(D):
    print("== S8")
    tp = D["timing_pairs"]
    cols = [("ratio_wB97MV_over_Skala", "ωB97M-V / Skala-1.1"), ("ratio_B3LYP_over_Skala", "B3LYP-D3(BJ) / Skala-1.1"), ("ratio_Skala_over_PBE", "Skala-1.1 / PBE")]
    shade = {"def2-TZVP": "#BBBBBB", "def2-TZVPD": "#666666", "def2-QZVPP": "#000000"}
    mk = {"def2-TZVP": "o", "def2-TZVPD": "s", "def2-QZVPP": "^"}
    fig, axs = plt.subplots(1, 3, figsize=(170 * MM, 62 * MM), sharex=True)
    rows = []
    tz = tp[tp.basis.isin(["def2-TZVP", "def2-TZVPD"])]
    for ax, (c, ttl), let in zip(axs, cols, "abc"):
        for bs in shade:
            q = tp[tp.basis == bs].dropna(subset=[c])
            ax.plot(q.nbf, q[c], ls="none", marker=mk[bs], ms=4.5, mfc=shade[bs] if bs != "def2-TZVP" else "white", mec="black", mew=0.6, label=bs)
            for _, r in q.iterrows():
                rows.append(dict(ratio=c, system=r.system, point_id=r.point_id, basis=bs, nbf=r.nbf, value=r[c]))
        ax.axhline(1, color="black", lw=0.8)
        med = tz[c].median()
        ax.axhline(med, color=GREY, lw=1.0, ls="--")
        ax.text(0.98, 0.04 if c != "ratio_Skala_over_PBE" else 0.04, f"median {med:.1f}", transform=ax.transAxes, ha="right", va="bottom", fontsize=8)
        ax.set_yscale("log"); ax.set_xscale("log")
        ax.set_xlabel("number of basis functions", fontsize=8); ax.set_ylabel(f"ratio of wall times\n{ttl}", fontsize=8)
        ax.set_xticks([10, 30, 100]); ax.set_xticklabels(["10", "30", "100"])
        ax.xaxis.set_minor_formatter(plt.NullFormatter()); ax.yaxis.set_minor_formatter(plt.NullFormatter())
        ax.set_yticks([0.2, 0.5, 1, 2, 5, 10]); ax.set_yticklabels(["0.2", "0.5", "1", "2", "5", "10"])
        panel(ax, let, x=0.03, y=0.99)
        ax.set_ylim(0.12, 22)
    axs[1].set_yticks([0.2, 0.5, 1, 2, 5, 10]); axs[1].set_ylim(0.12, 22)
    hd = [Line2D([0], [0], marker=mk[b_], ls="none", mfc=shade[b_] if b_ != "def2-TZVP" else "white", mec="black", ms=5, label=b_) for b_ in shade]
    hd += [Line2D([0], [0], color="black", lw=0.8, label="ratio = 1"), Line2D([0], [0], color=GREY, lw=1.0, ls="--", label="median, TZVP and TZVPD rows")]
    fig.legend(handles=hd, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.0), handletextpad=0.3, columnspacing=0.8, fontsize=8)
    fig.subplots_adjust(left=0.095, right=0.995, top=0.98, bottom=0.30, wspace=0.50)
    save_data(pd.DataFrame(rows), "figS8_timing_ratios")
    save(fig, "figS8_timing_ratios")
    anchor("figS8", "median wB97M-V/Skala ratio (TZVP + TZVPD rows)", 3.2, tz.ratio_wB97MV_over_Skala.median(), 0.05)


# ------------------------------------------------------------------ S9
def s9(D):
    print("== S9")
    b = D["bar"]
    fig, ax = fig_ax(85, 60)
    band(ax); zero_line(ax)
    rows = []
    for k, r in enumerate(["hcn", "nh3", "h2co"]):
        for j, m in enumerate(["Skala-1.1", "wB97M-V"]):
            q = err_row(b, r, m)
            x = k + (j - 0.5) * 0.3
            y0, y1 = q.signed_error_main_kcal, q.signed_error_large_kcal
            ax.plot([x], [y0], ls="none", zorder=4, **mstyle(m))
            ax.annotate("", xy=(x, y1), xytext=(x, y0), arrowprops=dict(arrowstyle="-|>", lw=0.9, color=METH[m]["c"], shrinkA=3, shrinkB=0), zorder=3)
            ax.plot([x], [y1], ls="none", marker=METH[m]["m"], ms=MS[METH[m]["m"]], mfc="white", mec=METH[m]["c"], mew=1.0, zorder=5)
            rows.append(dict(reaction=r, method=m, error_TZVP=y0, error_QZVPP=y1))
    ax.set_xticks(range(3)); ax.set_xticklabels(["HCN →\nHNC", "NH₃\ninversion", "H₂CO →\nH₂ + CO"])
    ax.set_xlim(-0.55, 2.55); ax.set_ylim(-4.5, 2.5)
    ax.set_ylabel(f"signed error\n({KC})")
    ax.text(0.99, 0.02, "n = 3 reactions", transform=ax.transAxes, ha="right", va="bottom", fontsize=8)
    hd = [handle("Skala-1.1", line=False), handle("wB97M-V", line=False),
          Line2D([0], [0], marker="o", ls="none", mfc="white", mec="black", ms=5, label="def2-QZVPP (filled: def2-TZVP)")]
    fig.legend(handles=hd, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.0), handletextpad=0.2, columnspacing=0.8, fontsize=8)
    fig.subplots_adjust(left=0.2, right=0.99, top=0.98, bottom=0.40)
    save_data(pd.DataFrame(rows), "figS9_large_basis_errors")
    save(fig, "figS9_large_basis_errors")
    print(pd.DataFrame(rows).round(2))
