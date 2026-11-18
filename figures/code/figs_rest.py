"""Figures 4, 8, 9, 10, 11, 12."""
from figcommon import *
from matplotlib.lines import Line2D


def err_row(b, r, m):
    q = b[(b.reaction == r) & (b.method == m)]
    return q.iloc[0] if len(q) else None


# ------------------------------------------------------------------ Figure 4
def fig04(D):
    print("== Figure 4")
    b = D["bar"]; dg = D["ccsdt_diagnostics"]
    fig, (a, c) = plt.subplots(1, 2, figsize=(170 * MM, 85 * MM), gridspec_kw=dict(width_ratios=[1, 1.05]))
    band(a, 0.5); zero_line(a)
    rows = []
    rxs = RX6
    for k, r in enumerate(rxs):
        raw = err_row(b, r, "CCSD(T)-fc"); cor = err_row(b, r, "CCSD(T)+dMP2")
        flagged = r in ("h2co", "ch3nc")
        col = GREY if flagged else "black"
        a.plot([k, k], [raw.signed_error_main_kcal, cor.signed_error_main_kcal], color=col, lw=0.8, zorder=3)
        if flagged:
            a.plot([k], [raw.signed_error_main_kcal], ls="none", marker="x", ms=6, mew=1.1, mfc="none", mec=GREY, zorder=4)
            a.plot([k], [cor.signed_error_main_kcal], ls="none", marker="+", ms=9, mew=2.0, mec=GREY, zorder=4)
        else:
            a.plot([k], [raw.signed_error_main_kcal], ls="none", marker="x", ms=6, mew=1.2, mec="black", zorder=4)
            a.plot([k], [cor.signed_error_main_kcal], ls="none", marker="+", ms=9, mew=2.0, mec="black", zorder=4)
        rows.append(dict(panel="a", reaction=r, raw=raw.signed_error_main_kcal, corrected=cor.signed_error_main_kcal, flagged=flagged))
    a.text(2.0, -0.62, "TS reference\nunstable", ha="center", va="top", fontsize=8, color="#555555", linespacing=0.95)
    a.text(5.45, 1.95, "literature value\nunverified", ha="right", va="bottom", fontsize=8, color="#555555", linespacing=0.95)
    a.text(5.45, -0.56, "0.5: tie threshold\nused for ranks", ha="right", va="top", fontsize=8, color="#555555", linespacing=0.95)
    a.set_xlim(-0.5, 5.5); a.set_ylim(-1.2, 2.5)
    a.set_xticks(range(6))
    tl = a.set_xticklabels([RXL[r] for r in rxs], rotation=45, ha="right", rotation_mode="anchor")
    tl[2].set_color(GREY); tl[5].set_color(GREY)
    a.set_xlabel("reaction")
    a.set_ylabel(f"CCSD(T) minus literature\n({KC})")
    panel(a, "a", x=0.02, y=0.99)
    # (b)
    order = ["hcn", "nh3", "h2co", "ch3nc", "sn2_ff", "sn2_clcl", "ch4"]
    rng = np.random.RandomState(7)
    c.axhline(0.02, color="black", lw=0.8, ls="--", zorder=1); c.axhline(0.03, color="black", lw=0.8, ls="--", zorder=1)
    nex = 0
    for k, sname in enumerate(order):
        q = dg[dg.system == sname]
        for i, (_, r) in enumerate(q.iterrows()):
            x = k + (rng.rand() - 0.5) * 0.5
            if r.excluded_from_reference:
                nex += 1
                c.plot([x], [r.t1], ls="none", marker="o", ms=6, mfc="white", mec="black", mew=0.9, zorder=4)
                c.plot([x], [r.t1], ls="none", marker="x", ms=6, mec="black", mew=0.9, zorder=5)
            else:
                c.plot([x], [r.t1], ls="none", marker="o", ms=4.5, mfc="black", mec="black", zorder=4)
            if r.stable == 0:
                c.text(x + 0.09, r.t1, "u", fontsize=8, va="center", ha="left")
            rows.append(dict(panel="b", system=sname, point_id=r.point_id, t1=r.t1, stable=r.stable, excluded=bool(r.excluded_from_reference)))
    c.text(6.6, 0.0205, "0.02", fontsize=8, ha="right", va="bottom"); c.text(6.6, 0.0305, "0.03", fontsize=8, ha="right", va="bottom")
    c.set_xlim(-0.5, 6.6); c.set_ylim(0, 0.037)
    c.set_xticks(range(7))
    tl = c.set_xticklabels(["HCN", "NH₃", "H₂CO", "CH₃NC", "F⁻ + CH₃F", "Cl⁻ + CH₃Cl", "CH₄"], rotation=45, ha="right", rotation_mode="anchor")
    c.set_xlabel("system"); c.set_ylabel("T1 diagnostic")
    panel(c, "b", x=0.02, y=0.99)
    c.text(0.99, 0.03, f"{len(dg)} points", transform=c.transAxes, ha="right", va="bottom", fontsize=8)
    hd = [Line2D([0], [0], marker="x", ls="none", color="black", ms=6, mew=1.2, label="CCSD(T), def2-TZVP"),
          Line2D([0], [0], marker="+", ls="none", color="black", ms=9, mew=2.0, label="CCSD(T) + MP2 basis correction"),
          Line2D([0], [0], marker="o", ls="none", color="black", ms=4.5, label="point used (b)"),
          Line2D([0], [0], marker="$\\otimes$", ls="none", color="black", ms=9, label="point removed (b)")]
    fig.legend(handles=hd, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.0), handletextpad=0.2, columnspacing=1.0)
    # the removed legend glyph with cross: overlay
    fig.subplots_adjust(left=0.115, right=0.995, top=0.98, bottom=0.39, wspace=0.27)
    data = pd.DataFrame(rows)
    save_data(data, "fig04_reference_quality")
    save(fig, "fig04_reference_quality")
    print(data[data.panel == "a"].round(3).to_string())
    LOG.append(("fig04", "points in panel (b) / excluded", "61 / 7", f"{len(dg)} / {nex}", len(dg) == 61 and nex == 7))
    ok = data[(data.panel == "a") & (data.flagged == False)]
    mx = ok.corrected.abs().max()
    print(f"  max |CCSD(T)+dMP2 - lit| (4 usable reactions) {mx:.3f}")
    LOG.append(("fig04", "claim: usable CCSD(T)+dMP2 within about 0.4 of literature (max abs)", "~0.4", f"{mx:.3f}", mx < 0.45))


# ------------------------------------------------------------------ Figure 8
def fig08(D):
    print("== Figure 8")
    b = D["bar"]; bs = D["basis_species"]
    fig = plt.figure(figsize=(170 * MM, 85 * MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.85, 1.15], left=0.115, right=0.995, top=0.975, bottom=0.30, wspace=0.30)
    a = fig.add_subplot(gs[0]); c = fig.add_subplot(gs[1])
    rxs = ["hcn", "nh3", "h2co"]
    band(a, 0.5); zero_line(a)
    rows = []
    wid = 0.30
    for k, r in enumerate(rxs):
        for j, m in enumerate(["Skala-1.1", "wB97M-V"]):
            q = err_row(b, r, m)
            x = k + (j - 0.5) * (wid + 0.04)
            sh = q.shift_main_to_large_kcal
            a.bar(x, sh, width=wid, color=METH[m]["c"], edgecolor="black", linewidth=0.4, zorder=3)
            txt = f"{q.signed_error_main_kcal:+.2f} → {q.signed_error_large_kcal:+.2f}".replace("-", "−")
            a.text(x, min(sh, 0) - 0.14 - (0.25 if sh > 0 else 0), txt, rotation=90, ha="center", va="top", fontsize=8)
            rows.append(dict(panel="a", reaction=r, method=m, shift_kcal=sh, err_main=q.signed_error_main_kcal, err_large=q.signed_error_large_kcal))
            if m == "Skala-1.1":
                qn = err_row(b, r, "Skala-1.1-noD3")
                a.plot([x], [qn.shift_main_to_large_kcal], ls="none", zorder=5, **mstyle("Skala-1.1-noD3", 1.0, marker="o"))
                rows.append(dict(panel="a", reaction=r, method="Skala-1.1-noD3", shift_kcal=qn.shift_main_to_large_kcal,
                                 err_main=qn.signed_error_main_kcal, err_large=qn.signed_error_large_kcal))
    a.set_xticks(range(3)); a.set_xticklabels([RXL[r].replace(" → ", "\n→ ").replace(" + ", "\n+ ") for r in rxs])
    a.set_xlim(-0.55, 2.55); a.set_ylim(-5.8, 1.4)
    a.set_ylabel(f"barrier change,\ndef2-TZVP → def2-QZVPP\n({KC})", fontsize=9)
    a.set_xlabel("reaction")
    panel(a, "a", x=0.02, y=0.99)
    a.text(0.99, 0.99, "n = 3 reactions", transform=a.transAxes, ha="right", va="top", fontsize=8)
    # (b)
    groups = [("hcn", ["HCN", "TS", "HNC"], "HCN → HNC"), ("nh3", ["NH3", "TS"], "NH₃ inversion"), ("h2co", ["H2CO", "TS", "H2", "CO"], "H₂CO → H₂ + CO")]
    nice = {"HCN": "HCN", "HNC": "HNC", "TS": "TS", "NH3": "NH₃", "H2CO": "H₂CO", "H2": "H₂", "CO": "CO"}
    x = 0; xt = []; xl = []; gl = []
    diffs = []
    for sysn, sp, gname in groups:
        x0 = x
        for s_ in sp:
            q = bs[(bs.system == sysn) & (bs.species == s_)].set_index("method")
            ys, yw = q.loc["Skala-1.1"].lowering_TZVP_to_QZVPP_kcal, q.loc["wB97M-V"].lowering_TZVP_to_QZVPP_kcal
            c.plot([x, x], [ys, yw], color="#AAAAAA", lw=1.0, zorder=2)
            c.plot([x], [yw], ls="none", zorder=4, **mstyle("wB97M-V"))
            c.plot([x], [ys], ls="none", zorder=5, **mstyle("Skala-1.1"))
            rows.append(dict(panel="b", system=sysn, species=s_, Skala=ys, wB97MV=yw, diff=ys - yw))
            diffs.append((sysn, s_, ys - yw))
            xt.append(x); xl.append(nice[s_])
            if sysn == "h2co" and s_ == "TS":
                xts = x; yts = (ys, yw)
            x += 1
        gl.append(((x0 + x - 1) / 2, gname))
        x += 0.8
    c.set_xticks(xt); c.set_xticklabels(xl)
    for gx, gname in gl:
        c.text(gx, -0.20, gname, transform=c.get_xaxis_transform(), ha="center", va="top", fontsize=8, fontweight="bold")
    c.set_xlim(-0.6, x - 0.8 + 0.1); c.set_ylim(-8.0, 0.0)
    c.annotate("", xy=(xts + 0.18, yts[0]), xytext=(xts + 0.18, yts[1]), arrowprops=dict(arrowstyle="<->", lw=0.8, color="black", shrinkA=0, shrinkB=0))
    c.text(xts + 0.32, (yts[0] + yts[1]) / 2, "2.2 larger\nfor Skala-1.1", ha="left", va="center", fontsize=8, linespacing=0.95)
    c.set_ylabel(f"total-energy change,\ndef2-TZVP → def2-QZVPP\n({KC})", fontsize=9)
    panel(c, "b", x=0.02, y=0.99)
    hd = [handle("Skala-1.1", line=False), handle("Skala-1.1-noD3", line=False), handle("wB97M-V", line=False)]
    fig.legend(handles=hd, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0), handletextpad=0.2, columnspacing=1.5)
    data = pd.DataFrame(rows)
    save_data(data, "fig08_basis_sensitivity")
    save(fig, "fig08_basis_sensitivity")
    h2 = b[(b.reaction == "h2co") & b.method.isin(["Skala-1.1", "wB97M-V"])].set_index("method").shift_main_to_large_kcal
    anchor("fig08", "H2CO shift Skala", -2.32, h2["Skala-1.1"], 0.006)
    anchor("fig08", "H2CO shift wB97M-V", -0.17, h2["wB97M-V"], 0.006)
    d = pd.DataFrame(diffs, columns=["system", "species", "diff"])
    print(d.round(3).to_string())
    others = d[~((d.system == "h2co") & (d.species == "TS"))]["diff"].abs()
    LOG.append(("fig08", "caption claim: only H2CO TS differs by >0.35 (max |diff| among other species)", "<=0.35", f"{others.max():.3f} (NH3 TS)", others.max() <= 0.35))
    LOG.append(("fig08", "corrected caption claim: only H2CO TS differs by more than 0.5 (max |diff| among other species)", "<=0.5", f"{others.max():.3f}", others.max() <= 0.5))
    note("fig08", f"Plan caption says H2CO TS is the only species where functionals differ by more than 0.35; data give {others.max():.2f} for NH3 TS. Caption reworded to 'more than 0.5'.")
    anchor("fig08", "H2CO TS difference (arrow label 2.2)", 2.2, abs(yts[0] - yts[1]), 0.06)
    sk = b[(b.reaction == "h2co") & (b.method == "Skala-1.1")].iloc[0]
    anchor("fig08", "Skala H2CO error main", 1.08, sk.signed_error_main_kcal, 0.006)
    anchor("fig08", "Skala H2CO error large", -1.24, sk.signed_error_large_kcal, 0.006)


# ------------------------------------------------------------------ Figure 9
def d3_table(D):
    """D3 contribution to each barrier from the UNROUNDED energies in results_2026-11-16/phase11_dft.csv
    (E(Skala-1.1) - E(Skala-1.1-noD3) at TS minus the same at the reactant/complex, converted to kcal/mol).
    master_barriers.csv holds barriers rounded to 4 decimals, so subtracting those rounds the Cl- case to -0.3345 -> -0.335."""
    HK = 627.5094740631
    d = pd.read_csv(RES + "phase11_dft.csv")
    d = d[(d.geometry_level == "own_b3lyp") & (d.converged == 1) & d.method.isin(["Skala-1.1", "Skala-1.1-noD3"])]
    pts = {"hcn": ("HCN", "TS"), "nh3": ("NH3", "TS"), "h2co": ("H2CO", "TS"), "ch3nc": ("CH3NC", "TS"),
           "sn2_ff": ("RC", "TS"), "sn2_clcl": ("RC@df", "TS")}   # Cl- + CH3Cl: density-fitted rows, as in master_barriers.csv
    b = D["bar"]

    def E(sysn, pid, m):
        q = d[(d.system == sysn) & (d.point_id == pid) & (d.method == m)]
        if sysn == "sn2_clcl" and pid == "TS":
            q = q[q.extra.str.contains('"df": true')]
        assert len(q) == 1, (sysn, pid, m, len(q))
        return float(q.energy_hartree.iloc[0])
    rows = []
    for r in RX6:
        a, ts = pts[r]
        bw = (E(r, ts, "Skala-1.1") - E(r, a, "Skala-1.1")) * HK
        bn = (E(r, ts, "Skala-1.1-noD3") - E(r, a, "Skala-1.1-noD3")) * HK
        s = err_row(b, r, "Skala-1.1").barrier_main_electronic_kcal
        anchor("fig09", f"raw-energy barrier {r} vs master_barriers.csv", s, bw, 6e-5)
        rows.append(dict(reaction=r, with_D3=bw, without_D3=bn, D3_contribution=bw - bn))
    return pd.DataFrame(rows)


def fig09(D):
    print("== Figure 9")
    t = d3_table(D)
    ref = pd.read_csv(RES + "skala_d3_table.csv").set_index("key").D3_contribution
    fig, ax = fig_ax(85, 60)
    ypos = np.arange(6)[::-1]
    for y, (_, r) in zip(ypos, t.iterrows()):
        ax.barh(y, r.D3_contribution, height=0.6, color=ORANGE, edgecolor="black", linewidth=0.4, zorder=3)
        v = r.D3_contribution
        ax.text(max(v, 0) + 0.02, y, f"{v:+.3f}".replace("-", "−"), va="center", ha="left", fontsize=8)
        anchor("fig09", f"D3 {r.reaction} vs R/skala_d3_table.csv", ref[r.reaction], v, 0.0006)
    ax.axvline(0, color="black", lw=0.8)
    for x in (-0.5, 0.5):
        ax.axvline(x, color=GREY, lw=0.6)
    ax.set_xlim(-0.55, 0.55); ax.set_ylim(-0.6, 5.6)
    ax.set_yticks(ypos)
    tl = ax.set_yticklabels([RXL[r] for r in RX6])
    tl[-1].set_color(GREY)
    ax.set_xlabel(f"D3 contribution to the\nbarrier ({KC})")
    ax.text(0.01, 0.0, "reference\nunverified", transform=ax.get_yaxis_transform(), fontsize=8, color="#555555", va="center", ha="left")
    ax.texts[-1].remove()
    fig.tight_layout()
    save_data(t, "fig09_d3_contribution")
    save(fig, "fig09_d3_contribution")
    print(t.round(4).to_string())
    LOG.append(("fig09", "max |D3 contribution| <= 0.33 (plan) ", "0.33", f"{t.D3_contribution.abs().max():.4f} (0.33 to 2 decimals)", round(t.D3_contribution.abs().max(), 2) == 0.33))
    print("  NOTE max |D3| = %.4f (Cl-+CH3Cl); plan quotes 0.33 (rounded down from 0.334)" % t.D3_contribution.abs().max())


# ------------------------------------------------------------------ Figure 10
def fig10(D):
    print("== Figure 10")
    b = D["bar"]; ir = D["internal_reference"]
    xtb = json.load(open(REC + "phase13_xtb_barriers.json"))["barrier_TS_minus_reactant_or_complex"]
    own = ir[ir.method == "Skala-1.1"].set_index("reaction").own_ccsdt_dmp2_kcal
    fig, (a, c) = plt.subplots(1, 2, figsize=(170 * MM, 85 * MM), gridspec_kw=dict(width_ratios=[1.0, 1.0]))
    meths = ["Skala-1.1", "wB97M-V", "B3LYP-D3BJ"]
    band(a); zero_line(a)
    rows = []
    for j, m in enumerate(meths):
        off = (j - 1) * 0.22
        for k, r in enumerate(RX6):
            q = err_row(b, r, m)
            y = q.xtb_geometry_penalty_kcal
            st = mstyle(m)
            if r == "h2co":
                st.update(mfc="white", mec=METH[m]["c"], mew=1.1)
            a.plot([k + off, k + off], [0, y], color=METH[m]["c"], lw=0.8, zorder=3)
            a.plot([k + off], [y], ls="none", zorder=4, **st)
            rows.append(dict(panel="a", reaction=r, method=m, x=k + off, xtb_geometry_penalty_kcal=y))
    q = err_row(b, "sn2_clcl", "Skala-1.1")
    a.annotate("+5.7, see Fig. 7", (4 - 0.22, q.xtb_geometry_penalty_kcal), xytext=(3.45, 5.7), fontsize=8, ha="right", va="center",
               arrowprops=dict(arrowstyle="-", lw=0.5, color="black"))
    a.text(2, 1.15, "no xTB transition state\nfound; not a test", ha="center", va="bottom", fontsize=8, linespacing=0.95)
    a.set_xlim(-0.5, 5.5); a.set_ylim(-2.5, 6.8)
    a.set_xticks(range(6))
    a.set_xticklabels([RXL[r] for r in RX6], fontsize=8, rotation=45, ha="right", rotation_mode="anchor")
    a.set_xlabel("reaction")
    a.set_ylabel(f"barrier at xTB minus\nbarrier at B3LYP\ngeometries ({KC})", fontsize=9)
    panel(a, "a", x=0.02, y=0.99)
    # (b)
    band(c); zero_line(c)
    for k, r in enumerate(RX6):
        ref = err_row(b, r, "Skala-1.1").reference_kcal
        openm = r == "ch3nc"
        if r == "ch3nc":
            ref = own[r]
        xv = xtb[r]
        if xv is not None:
            y = xv - ref
            c.plot([k - 0.15, k - 0.15], [0, y], color=GREY, lw=0.8, zorder=3)
            st = mstyle("GFN2-xTB");
            if openm:
                st.update(mfc="white", mec=GREY, mew=1.0)
            c.plot([k - 0.15], [y], ls="none", zorder=4, **st)
            rows.append(dict(panel="b", reaction=r, series="GFN2-xTB barrier", x=k - 0.15, error=y, reference=ref))
        else:
            c.text(k - 0.15, 0.8, "none", ha="center", va="bottom", fontsize=8, color="#555555")
        if r != "h2co":
            yb = err_row(b, r, "Skala-1.1").barrier_at_xtb_geometry_electronic_kcal - ref
            c.plot([k + 0.15, k + 0.15], [0, yb], color=ORANGE, lw=0.8, zorder=3)
            st = mstyle("Skala-1.1")
            if openm:
                st.update(mfc="white", mec=ORANGE, mew=1.1)
            c.plot([k + 0.15], [yb], ls="none", zorder=4, **st)
            rows.append(dict(panel="b", reaction=r, series="Skala-1.1 at xTB geometry", x=k + 0.15, error=yb, reference=ref))
    c.set_xlim(-0.5, 5.5); c.set_ylim(-12, 28)
    c.set_xticks(range(6))
    tl = c.set_xticklabels([RXL[r] for r in RX6], fontsize=8, rotation=45, ha="right", rotation_mode="anchor")
    tl[-1].set_color(GREY)
    c.text(4.7, 18.5, "reference\nunverified", ha="right", va="center", fontsize=8, color="#555555", linespacing=0.95)
    c.set_xlabel("reaction")
    c.set_ylabel(f"signed error of the\nbarrier ({KC})")
    panel(c, "b", x=0.02, y=0.99)
    hd = [handle("Skala-1.1", line=False), handle("wB97M-V", line=False), handle("B3LYP-D3BJ", line=False),
          Line2D([0], [0], marker="*", ls="none", mfc="#7F7F7F", mec="#7F7F7F", ms=8, label="GFN2-xTB barrier (b)"),
          Line2D([0], [0], marker="o", ls="none", mfc="white", mec="black", ms=5, label="open: H₂CO (a), CH₃NC (b)")]
    fig.legend(handles=hd, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.0), handletextpad=0.2, columnspacing=0.9, fontsize=8)
    fig.subplots_adjust(left=0.115, right=0.995, top=0.98, bottom=0.39, wspace=0.30)
    data = pd.DataFrame(rows)
    save_data(data, "fig10_scout")
    save(fig, "fig10_scout")
    print(data.round(3).to_string())
    anchor("fig10", "Skala Cl penalty", 5.68, err_row(b, "sn2_clcl", "Skala-1.1").xtb_geometry_penalty_kcal, 0.006)
    ymax = max(abs(v) for v in data[data.panel == "b"].error)
    print("  max |error| panel b", ymax)


# ------------------------------------------------------------------ Figure 11
def xtb_cx(rx):
    """C-X distances (A) in the xTB TS geometry, from results geoms."""
    f = GEO + f"{rx}_TS_xtb_gfn2.xyz"
    L = open(f).read().splitlines()[2:]
    at = [(l.split()[0], np.array([float(v) for v in l.split()[1:4]])) for l in L if l.strip()]
    c = [p for s, p in at if s == "C"][0]
    return sorted(np.linalg.norm(p - c) for s, p in at if s in ("F", "Cl"))


def fig11(D):
    print("== Figure 11")
    ir = D["internal_reference"]; mins = D["sn2_minima"]; b = D["bar"]
    own = ir[ir.method == "Skala-1.1"].set_index("reaction").own_ccsdt_dmp2_kcal
    dist = {rx: xtb_cx(rx) for rx in ("sn2_ff", "sn2_clcl")}
    print("  xTB C-X distances from xyz:", {k: [round(v, 3) for v in x] for k, x in dist.items()})
    fig, axs = plt.subplots(1, 2, figsize=(85 * MM, 85 * MM), sharey=True)
    rows = []
    cats = ["xTB", "B3LYP", "min"]
    for ax, rx, tt in zip(axs, ["sn2_ff", "sn2_clcl"], ["F⁻ + CH₃F", "Cl⁻ + CH₃Cl"]):
        band(ax); zero_line(ax)
        for m in ["Skala-1.1", "wB97M-V"]:
            xtbv = err_row(b, rx, m).barrier_at_xtb_geometry_electronic_kcal - own[rx]
            b3 = ir[(ir.reaction == rx) & (ir.method == m)].signed_error_vs_own_ccsdt_kcal.iloc[0]
            mm = mins[(mins.reaction == rx) & (mins.method == m)].iloc[0].error_vertex_vs_ccsdt
            ys = [xtbv, b3, mm]
            ax.plot(range(3), ys, color=METH[m]["c"], lw=1.0, zorder=2)
            base = mstyle(m)
            ax.plot([0], [ys[0]], ls="none", zorder=4, **{**base, "fillstyle": "left", "mfc": METH[m]["c"], "mfcalt": "white", "mec": METH[m]["c"], "mew": 1.0})
            ax.plot([1], [ys[1]], ls="none", zorder=4, **base)
            ax.plot([2], [ys[2]], ls="none", zorder=4, **{**base, "mfc": "white", "mec": METH[m]["c"], "mew": 1.1})
            for cat, y in zip(cats, ys):
                rows.append(dict(reaction=rx, method=m, geometry=cat, error_vs_ccsdt_dmp2_kcal=y))
        d = dist[rx]
        mnS = mins[(mins.reaction == rx) & (mins.method == "Skala-1.1")].vertex_cx_A.iloc[0]
        mnW = mins[(mins.reaction == rx) & (mins.method == "wB97M-V")].vertex_cx_A.iloc[0]
        b3d = mins[(mins.reaction == rx)].b3lyp_ts_cx_A.iloc[0]
        labs = [f"xTB TS\n{d[0]:.3f} / {d[1]:.3f}", f"B3LYP TS\n{b3d:.3f}", f"own minimum\n{mnS:.3f} / {mnW:.3f}"]
        ax.set_xticks(range(3)); ax.set_xticklabels(labs, rotation=90, ha="center", va="top", fontsize=8, linespacing=1.0)
        ax.set_xlim(-0.4, 2.4)
        ax.set_ylim(-5.0, 3.9)
        ax.text(0.5, 1.0, tt, transform=ax.transAxes, ha="center", va="top", fontsize=8, fontweight="bold")
        rows.append(dict(reaction=rx, method="distances", geometry=f"xTB {d[0]:.3f}/{d[1]:.3f}; B3LYP {b3d:.3f}; min Skala {mnS:.3f}, wB97M-V {mnW:.3f}", error_vs_ccsdt_dmp2_kcal=np.nan))
    axs[0].set_ylabel(f"barrier minus CCSD(T)\nbarrier at the B3LYP TS\n({KC})")
    axs[1].tick_params(labelleft=False, left=False); axs[1].spines["left"].set_visible(False)
    panel(axs[0], "a", x=0.02, y=0.99); panel(axs[1], "b", x=0.02, y=0.99)
    hd = [handle("Skala-1.1", line=False), handle("wB97M-V", line=False),
          Line2D([0], [0], marker="o", ls="none", fillstyle="left", mfc="black", mfcalt="white", mec="black", ms=5, label="GFN2-xTB TS"),
          Line2D([0], [0], marker="o", ls="none", mfc="black", mec="black", ms=5, label="B3LYP TS"),
          Line2D([0], [0], marker="o", ls="none", mfc="white", mec="black", ms=5, label="own minimum")]
    fig.legend(handles=hd, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0), handletextpad=0.2, columnspacing=0.6, fontsize=8)
    fig.subplots_adjust(left=0.23, right=0.99, top=0.985, bottom=0.46, wspace=0.05)
    data = pd.DataFrame(rows)
    save_data(data, "fig11_sn2_geometry")
    save(fig, "fig11_sn2_geometry")
    print(data.round(3).to_string())
    g = lambda rx, m, geo: data[(data.reaction == rx) & (data.method == m) & (data.geometry == geo)].error_vs_ccsdt_dmp2_kcal.iloc[0]
    for (rx, m, geo, e) in [("sn2_ff", "Skala-1.1", "B3LYP", 1.68), ("sn2_clcl", "Skala-1.1", "B3LYP", -3.32), ("sn2_ff", "wB97M-V", "B3LYP", -1.05),
                            ("sn2_clcl", "wB97M-V", "B3LYP", -0.11), ("sn2_ff", "Skala-1.1", "min", 1.28), ("sn2_clcl", "Skala-1.1", "min", -4.08),
                            ("sn2_ff", "wB97M-V", "min", -1.26), ("sn2_clcl", "wB97M-V", "min", -0.27),
                            ("sn2_ff", "Skala-1.1", "xTB", 0.34), ("sn2_clcl", "Skala-1.1", "xTB", 2.36), ("sn2_clcl", "wB97M-V", "xTB", 0.10)]:
        anchor("fig11", f"{rx} {m} {geo}", e, g(rx, m, geo), 0.006)
    anchor("fig11", "sn2_ff wB97M-V xTB (plan -2.24 from rounded inputs)", -2.24, g("sn2_ff", "wB97M-V", "xTB"), 0.011)
    for rx, ex in [("sn2_ff", [1.757, 1.896]), ("sn2_clcl", [2.201, 2.236])]:
        for got, e in zip(dist[rx], ex):
            anchor("fig11", f"xTB C-X distance {rx}", e, got, 0.0006)


# ------------------------------------------------------------------ Figure 12
def fig12(D):
    print("== Figure 12")
    tp = D["timing_pairs"]; tl = D["timing_long"]; s = D["sum"]; L = D["leave_one_out"]
    fig, (a, c) = plt.subplots(1, 2, figsize=(170 * MM, 85 * MM), gridspec_kw=dict(width_ratios=[1, 1]))
    t = tp[tp.basis.isin(["def2-TZVP", "def2-TZVPD"])]
    rows = []
    cols = {"Skala-1.1": "Skala-1.1", "wB97M-V": "wB97M-V", "B3LYP-D3BJ": "B3LYP-D3BJ", "r2SCAN": "r2SCAN", "PBE": "PBE"}
    for m in DFT5:
        q = t.dropna(subset=[m])
        a.plot(q.nbf, q[m], ls="none", alpha=0.9, zorder=3, **mstyle(m, 0.85))
        for _, r in q.iterrows():
            rows.append(dict(panel="a", method=m, system=r.system, point_id=r.point_id, nbf=r.nbf, seconds=r[m]))
    x = tl[(tl.system == "sn2_clcl") & (tl.point_id == "TS") & (tl.method == "Skala-1.1") & (tl.seconds == 436.9)].iloc[0]
    a.plot([x.nbf], [x.seconds], ls="none", marker="o", ms=8, mfc="white", mec=ORANGE, mew=1.2, zorder=5)
    a.annotate("SCF retries", (x.nbf, x.seconds), xytext=(70, 300), fontsize=8, ha="center", va="center", arrowprops=dict(arrowstyle="-", lw=0.5))
    rows.append(dict(panel="a", method="Skala-1.1 (SCF retries, density fitting)", system="sn2_clcl", point_id="TS", nbf=x.nbf, seconds=x.seconds))
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlim(9, 220); a.set_ylim(0.3, 800)
    from matplotlib.ticker import NullFormatter
    a.set_xticks([10, 20, 50, 100, 200]); a.set_xticklabels(["10", "20", "50", "100", "200"])
    a.xaxis.set_minor_formatter(NullFormatter()); a.yaxis.set_minor_formatter(NullFormatter())
    a.set_yticks([1, 10, 100]); a.set_yticklabels(["1", "10", "100"])
    a.set_xlabel("number of basis functions")
    a.set_ylabel("wall time per single point (s),\nTesla T4")
    panel(a, "a", x=0.02, y=0.99)
    print("  nbf range plotted:", t.nbf.min(), t.nbf.max(), "; rows", len(t))
    # (b)
    com = tp.dropna(subset=["PBE"])
    print("  common structures:", len(com))
    med = {m: com[m].median() for m in DFT5}
    hl = s[(s.set == "HEADLINE") & (s.basis_kind == "main")].set_index("method").MAE
    off = {"Skala-1.1": (6, 4, "left", "bottom"), "wB97M-V": (-6, 6, "right", "bottom"), "B3LYP-D3BJ": (6, 6, "left", "bottom"),
           "r2SCAN": (-6, 6, "right", "bottom"), "PBE": (6, 6, "left", "bottom")}
    for m in DFT5:
        loo = L[(L.method == m) & (L.removed != "none")].MAE
        c.plot([med[m], med[m]], [loo.min(), loo.max()], color=METH[m]["c"], lw=1.0, zorder=2)
        for yv in (loo.min(), loo.max()):
            c.plot([med[m] * 0.97, med[m] * 1.03], [yv, yv], color=METH[m]["c"], lw=1.0, zorder=2)
        c.plot([med[m]], [hl[m]], ls="none", zorder=4, **mstyle(m, 1.2))
        rows.append(dict(panel="b", method=m, median_seconds=med[m], MAE=hl[m], loo_min=loo.min(), loo_max=loo.max(), n_common=len(com)))
    lab = {"Skala-1.1": ("Skala-1.1", (1.1, 0.0), "left", "center"), "wB97M-V": ("ωB97M-V", (1.1, 0.0), "left", "center"),
           "B3LYP-D3BJ": ("B3LYP-D3(BJ)", (1.1, 0.0), "left", "center"), "r2SCAN": ("r²SCAN", (1.1, 0.0), "left", "center"), "PBE": ("PBE", (1.1, 0.0), "left", "center")}
    for m in DFT5:
        txt, (dx, dy), ha, va = lab[m]
        c.text(med[m] * (1.18 if dx > 0 else 1 / 1.18), hl[m] + dy, txt, ha=ha, va=va, fontsize=8)
    from matplotlib.ticker import NullFormatter
    c.set_xscale("log"); c.set_xlim(1.5, 80); c.set_ylim(0, 6)
    c.xaxis.set_minor_formatter(NullFormatter())
    c.set_xticks([2, 5, 10, 20, 40]); c.set_xticklabels(["2", "5", "10", "20", "40"])
    c.set_xlabel("median wall time per single point (s)")
    c.set_ylabel(f"mean absolute error, n = 5 ({KC})")
    panel(c, "b", x=0.02, y=0.99)
    hd = legend_handles(DFT5, line=False)
    hd.append(Line2D([0], [0], marker="o", ls="none", mfc="white", mec=ORANGE, ms=7, label="Cl⁻ + CH₃Cl TS, Skala-1.1"))
    fig.legend(handles=hd, loc="lower center", ncol=6, frameon=False, bbox_to_anchor=(0.5, 0.0), handletextpad=0.1, columnspacing=0.7, fontsize=8)
    fig.subplots_adjust(left=0.095, right=0.995, top=0.98, bottom=0.25, wspace=0.28)
    data = pd.DataFrame(rows)
    save_data(data, "fig12_cost_accuracy")
    save(fig, "fig12_cost_accuracy")
    print("  medians (19 common structures):", {k: round(v, 2) for k, v in med.items()})
    r = t.ratio_wB97MV_over_Skala.median()
    anchor("fig12", "median wB97M-V/Skala ratio (TZVP + TZVPD rows, n=%d)" % len(t), 3.2, r, 0.05)
    print("  ratio medians: all rows %.3f, 19 common %.3f" % (tp.ratio_wB97MV_over_Skala.median(), com.ratio_wB97MV_over_Skala.median()))
    note("fig12", "Plan expects medians near Skala 6, wB97M-V 20, B3LYP 12, r2SCAN 2.5, PBE 2.3 s; recomputed over the 19 common structures: "
         + ", ".join(f"{k} {v:.2f}" for k, v in med.items()) + ". Skala (7.4) and wB97M-V (22.7) differ from the plan's rough figures; recomputed values are plotted.")
    note("fig12", "Median ratio wB97M-V/Skala is 3.22 over the 50 TZVP/TZVPD structures, 3.32 over all 59 rows including QZVPP, 2.86 over the 19 structures common to all five methods.")
