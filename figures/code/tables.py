"""Tables 1-5 and S1-S8 (CSV + Markdown), VISUALS_PLAN.md section 3. No colour, no shading."""
from figcommon import *
from figs_extra import internal_stats

MLAB = {"Skala-1.1": "Skala-1.1", "Skala-1.1-noD3": "Skala-1.1 without D3", "wB97M-V": "ωB97M-V", "B3LYP-D3BJ": "B3LYP-D3(BJ)",
        "r2SCAN": "r²SCAN", "PBE": "PBE", "CCSD(T)-fc": "CCSD(T) raw", "CCSD(T)+dMP2": "CCSD(T) + MP2 corr."}


def f2(v, signed=False):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return ""
    return (f"{v:+.2f}" if signed else f"{v:.2f}").replace("-", "−")


def build_tables(D):
    b = D["bar"]
    # ---------------------------------------------------------------- Table 1
    t1 = pd.read_csv(RES + "table1_systems.csv")
    t1 = t1[t1.key != "ch4"].copy()
    bh76 = {"hcn": "yes", "sn2_ff": "yes", "sn2_clcl": "yes"}
    t1["in_BH76"] = t1.key.map(lambda k: bh76.get(k, "no"))
    t1 = t1.rename(columns={"system": "reaction", "cls": "class", "zero_point": "zero_point_in_reference", "verified": "reference_verified"})
    t1 = t1[["number", "reaction", "class", "species", "literature", "source", "zero_point_in_reference", "reference_verified", "in_BH76"]]
    save_table(t1, "table1_reactions", "Table 1. Reactions, reference values (kcal/mol), sources, zero-point convention of the reference, verification status, BH76 membership. CH3NC is computed but in no statistic. Source: R/table1_systems.csv; BH76 membership from the plan.")
    # ---------------------------------------------------------------- Table 2
    meths = ["Skala-1.1", "Skala-1.1-noD3", "wB97M-V", "B3LYP-D3BJ", "r2SCAN", "PBE", "CCSD(T)-fc", "CCSD(T)+dMP2"]
    rows_n, rows_s = [], []
    for grp, rxs in [("headline", RX), ("not in statistics", ["ch3nc"])]:
        for r in rxs:
            q = b[b.reaction == r]
            ref = q.reference_kcal.iloc[0]
            dn = dict(set=grp, reaction=RXL[r].replace("₂", "2").replace("₃", "3").replace("⁻", "-"), reference_kcal=ref)
            ds = dict(set=grp, reaction=RXL[r], reference_kcal=f"{ref:.2f}")
            for m in meths:
                x = q[q.method == m].iloc[0]
                dn[f"{m}_barrier"] = x.barrier_main_used_kcal; dn[f"{m}_error"] = x.signed_error_main_kcal
                flag = "*" if (r == "h2co" and m.startswith("CCSD")) else ""
                ds[MLAB[m]] = f"{x.barrier_main_used_kcal:.2f} ({f2(x.signed_error_main_kcal, True)}){flag}"
            rows_n.append(dn); rows_s.append(ds)
    t2 = pd.DataFrame(rows_n)
    save_table(t2, "table2_barriers_numeric", None)
    t2s = pd.DataFrame(rows_s)
    t2s.to_csv(TAB + "table2_barriers.csv", index=False)
    with open(TAB + "table2_barriers.md", "w") as f:
        f.write("**Table 2. Computed forward barriers (kcal/mol, used value) with signed error against the reference in parentheses (computed minus reference). def2-TZVP; def2-TZVPD for the two SN2 reactions; B3LYP-D3(BJ) geometries. H2CO values include the B3LYP zero-point difference of -4.988. "
                "* H2CO CCSD(T): withdrawn by the check (TS Hartree-Fock unstable), information only, in no statistic. CH3NC: reference unverified, in no statistic. CCSD(T) + MP2 corr. = CCSD(T)/def2-TZVP plus DF-MP2 (def2-QZVPP minus def2-TZVP).**\n\n")
        f.write(md_table(t2s) + "\n")
    print("table table2_barriers")
    # ---------------------------------------------------------------- Table 3
    s = D["sum"]
    rows = []
    for m in DFT6 + ["CCSD(T)+dMP2"]:
        r = s[(s.set == "HEADLINE") & (s.basis_kind == "main") & (s.method == m)]
        if len(r):
            r = r.iloc[0]
            rows.append(dict(block="(i) literature reference, def2-TZVP(D)", method=MLAB[m], n=int(r.n_reactions), MAE=r.MAE, MSE=r.MSE, MaxAE=r.MaxAE))
    st = internal_stats(D)
    for _, r in st.iterrows():
        rows.append(dict(block="(ii) internal CCSD(T)+MP2 reference, with CH3NC, H2CO excluded", method=MLAB[r.method], n=int(r.n), MAE=r.MAE, MSE=r.MSE, MaxAE=r.MaxAE))
    for m in ["Skala-1.1", "Skala-1.1-noD3", "wB97M-V"]:
        r = s[(s.set == "HEADLINE") & (s.basis_kind == "large") & (s.method == m)].iloc[0]
        rows.append(dict(block="(iii) def2-QZVPP, HCN, NH3, H2CO, literature reference", method=MLAB[m], n=int(r.n_reactions), MAE=r.MAE, MSE=r.MSE, MaxAE=r.MaxAE))
    t3 = pd.DataFrame(rows)
    save_table(t3, "table3_statistics", "Table 3. Mean absolute error (MAE), mean signed error (MSE) and maximum absolute error (MaxAE) in kcal/mol. (i) n = 5 reactions (CCSD(T) rows: 4 reactions, H2CO has no usable CCSD(T)); (ii) internal reference; (iii) def2-QZVPP, n = 3. Source: W/summary_stats.csv (HEADLINE), W/derived_internal_reference.csv.", ".2f")
    # ---------------------------------------------------------------- Table 4
    fo = pd.read_csv(RES + "findings_outcome.csv")
    sens = {
        "F1": "Margin to the threshold is 0.14 kcal/mol (1.443 against 1.08 + 0.5 = 1.58); the verdict needs HCN: without HCN Skala 1.73 and wB97M-V 0.46, F1 not supported. With CH3NC (n = 6) the criterion is not met (1.52 > 1.46).",
        "F2": "Three paths, not four (CH4 has no DFT rows). Against raw CCSD(T) the ranges are 1.06, 0.46, 4.81; the HCN and H2CO CCSD(T) values are flagged by the check, so only NH3 is clean.",
        "F3": "The rank-unchanged part cannot be evaluated (QZVPP only for Skala, Skala without D3 and wB97M-V). Only H2CO (-2.32) exceeds 1.0 kcal/mol.",
        "F4": "H2CO is not a TS test (xTB found no peak). Without CH3NC: Skala 3 of 5, wB97M-V 4 of 5, B3LYP-D3(BJ) 5 of 5 below 1.0.",
        "F5": "Needs the 0.5 tie rule: strict ranks of Skala are 1, 5, 2, 2, 2 (3 of 5 = 60 %, below 75 %); verdict margin is one reaction.",
        "F6": "Mixes two systems: CH4 8 of 9 (88.9 %), HCN grid 111 of 369 (30.1 %); grid stability flag not stored.",
    }
    t4 = pd.DataFrame(dict(finding=fo.finding, criterion_supports=fo.criterion_supports, criterion_negative=fo.criterion_negative,
                           measured=fo.measured, n=fo.n, verdict=fo.verdict, sensitivity_note=fo.finding.map(sens)))
    save_table(t4, "table4_findings", "Table 4. The six criteria fixed before the calculations (F1 to F6), measured value, n, verdict and one sensitivity note each. Source: R/findings_outcome.csv; sensitivity notes from dossier section I.")
    # ---------------------------------------------------------------- Table 5
    m = D["sn2_minima"].copy()
    m["reaction"] = m.reaction.map(RXL)
    m["method"] = m.method.map(MLAB)
    t5 = m[["reaction", "method", "vertex_cx_A", "vertex_energy_kcal", "error_vertex_vs_ccsdt", "error_vertex_vs_lit", "b3lyp_ts_cx_A", "energy_at_b3lyp_ts_kcal", "ref_ccsdt_dmp2_kcal", "literature_kcal"]]
    t5.columns = ["reaction", "functional", "minimum_C-X_distance_A", "energy_at_minimum_kcal", "error_vs_CCSD(T)+MP2_kcal", "error_vs_literature_kcal",
                  "B3LYP_TS_C-X_distance_A", "energy_at_B3LYP_TS_kcal", "CCSD(T)+MP2_barrier_kcal", "literature_barrier_kcal"]
    save_table(t5, "table5_sn2_minima", "Table 5. Minimum of the rigid symmetric stretch (three-point parabola), energy relative to the ion-molecule complex, errors against both references, and the energy at the B3LYP transition state. Upper bounds for relaxed transition-state energies. Source: W/derived_sn2_minima.csv.")
    # ---------------------------------------------------------------- S1
    e = D["err"]
    cols = ["reaction", "lit_id", "quantity", "is_primary_barrier", "in_headline", "in_pooled_stats", "method", "basis_kind", "basis", "geometry_level", "electronic_kcal", "zpe_added_kcal",
            "computed_kcal", "reference_value_kcal", "reference_verified", "signed_error_kcal", "abs_error_kcal", "method_withdrawn_by_check", "note"]
    save_table(e[cols], "tableS1_all_scored_rows", "Table S1. All scored rows (L1 to L10), kcal/mol. Source: W/master_errors.csv.")
    # ---------------------------------------------------------------- S2
    P = D["path"]; dg = D["ccsdt_diagnostics"]
    pv = P.pivot_table(index=["reaction", "reference", "point_id", "coordinate_deg_or_A", "ref_relative_kcal"], columns="method", values="signed_error_kcal").reset_index()
    pv.columns.name = None
    pv = pv.rename(columns={c: f"error_{c}" for c in DFT6 if c in pv.columns})
    d2 = dg[dg.geometry_level == "own_scan"][["system", "point_id", "t1", "stable", "excluded_from_reference"]].rename(columns={"system": "reaction"})
    ts2 = pv.merge(d2, on=["reaction", "point_id"], how="left")
    ts2 = ts2[ts2.reaction.isin(["hcn", "nh3", "h2co"])]
    save_table(ts2, "tableS2_path_points", "Table S2. Path points: coordinate (degrees or Angstrom), reference energy relative to the first point, signed errors (kcal/mol), T1, RHF stability (1 stable, 0 unstable), excluded from the reference or not. Sources: W/path_profiles.csv, W/derived_ccsdt_diagnostics.csv.")
    # ---------------------------------------------------------------- S3
    sy = D["sym"][["reaction", "method", "point_id", "kind", "cx_distance_A", "density_fitting", "converged", "energy_hartree", "d3_hartree", "rel_to_RC_B3LYPgeom_kcal", "rel_to_RC_noD3_kcal", "attempt", "seconds"]]
    save_table(sy, "tableS3_sn2_stretch", "Table S3. SN2 symmetric-stretch energies (kcal/mol relative to the ion-molecule complex), convergence and wall time (s). Source: W/sn2_symscan.csv.")
    # ---------------------------------------------------------------- S4
    bs = D["basis_species"]
    save_table(bs, "tableS4a_basis_species", "Table S4a. Total energies (hartree) of each species at def2-TZVP and def2-QZVPP and the change in kcal/mol. Source: W/derived_basis_species.csv.")
    bb = b[b.reaction.isin(["hcn", "nh3", "h2co"]) & b.method.isin(["Skala-1.1", "Skala-1.1-noD3", "wB97M-V"])][
        ["reaction", "method", "basis_main", "barrier_main_electronic_kcal", "zpe_added_kcal", "barrier_main_used_kcal", "basis_large", "barrier_large_used_kcal", "shift_main_to_large_kcal", "signed_error_main_kcal", "signed_error_large_kcal"]]
    save_table(bb, "tableS4b_basis_barriers", "Table S4b. Barriers at the two basis sets and the shift (kcal/mol). Source: W/master_barriers.csv.")
    # ---------------------------------------------------------------- S5
    g = pd.read_csv(RES + "geometry_table.csv")
    xtbj = json.load(open(REC + "phase13_xtb_barriers.json"))["barrier_TS_minus_reactant_or_complex"]
    rows = []
    for _, r in g.iterrows():
        k = r.key
        d = dict(reaction=RXL[k], key_coordinate=r.coordinate, B3LYP_TS=r.own_b3lyp, xTB_TS=r.xtb_gfn2, xTB_minus_B3LYP=r.xtb_minus_b3lyp)
        for m in ["Skala-1.1", "wB97M-V", "B3LYP-D3BJ"]:
            d[f"penalty_{m}_kcal"] = b[(b.reaction == k) & (b.method == m)].xtb_geometry_penalty_kcal.iloc[0]
        d["GFN2-xTB_own_barrier_kcal"] = xtbj.get(k)
        rows.append(d)
    save_table(pd.DataFrame(rows), "tableS5_scout", "Table S5. Scout test: key coordinate at the B3LYP and GFN2-xTB transition states, geometry penalties (kcal/mol) and the GFN2-xTB barrier itself (empty: none found). H2CO: xTB found no peak, the xTB TS equals the B3LYP TS. Sources: R/geometry_table.csv, W/master_barriers.csv, R/phase13_xtb_barriers.json.")
    # ---------------------------------------------------------------- S6 (transcribed from dossier A.3 and A.4)
    s6 = pd.DataFrame([
        ("Skala", "Skala-1.1 (microsoft/skala-1.1, HF revision 99b5ed87e5f69d9216e1f9e30148b922eaea1241); with_dftd3=True; D3 set b3lyp5, D3(BJ), two-body (atm=False); GPU route", "dossier A.3"),
        ("Skala weights SHA-256", "7f3e8622e1eb520ccd88a55464c3e359ac4d7e5ccbd1fb77a26afa1e1c20a5cd", "dossier A.3"),
        ("Skala without D3", "same SCF, stored D3 energy removed (not an independent calculation)", "dossier A.3"),
        ("Other DFT", "PBE, r2SCAN, B3LYP-D3(BJ) (pyscf d3bj), wB97M-V (VV10 default nonlocal grid); grid level 3, SCF tolerance 1e-9, no density fitting unless stated", "dossier A.3"),
        ("CCSD(T)", "frozen core, def2-TZVP (def2-TZVPD SN2), RHF reference; T1 and RHF stability stored", "dossier A.3"),
        ("CCSD(T) + MP2 correction", "CCSD(T) + [DF-MP2/def2-QZVPP(D) minus DF-MP2/def2-TZVP(D)], frozen core, same geometry", "dossier A.3"),
        ("Geometries", "DF-B3LYP-D3(BJ)/def2-TZVP(D), geomeTRIC, analytic Hessian; geometry pack SHA-256 81e4533319f5fc593cc3e81daf93e6e4db4998df01eea275b3c813fc87dafae6", "dossier A.1, A.3"),
        ("xTB", "GFN2-xTB (tblite 0.8.0)", "dossier A.3"),
        ("Skala density fitting", "Cl- + CH3Cl TS did not converge without density fitting; TS and partners are density-fitted rows", "dossier A.3"),
        ("Hardware, geometry pack build", "Intel Xeon 2.10 GHz, 4 cores, 16 GB", "dossier A.4"),
        ("Software, geometry pack build", "pyscf 2.15.0, geometric 1.1.1, tblite 0.8.0, ase 3.29.0, numpy 2.5.4, scipy 1.18.1, pandas 3.1.1, pyscf-dispersion 1.7.0, python 3.13.18", "dossier A.4"),
        ("Hardware, CPU notebook (CCSD(T), MP2)", "Intel Xeon 2.20 GHz, 2 cores, 13 GB", "dossier A.4"),
        ("Software, CPU notebook", "pyscf 2.15.0, numpy 2.4.6, scipy 1.16.3, pandas 2.2.3, torch 2.14.1+cpu, pyscf-dispersion 1.7.0, python 3.13.18", "dossier A.4"),
        ("Hardware, GPU notebook (all DFT)", "Intel Xeon 2.00 GHz, 8 cores, 51 GB, Tesla T4 (15 GB)", "dossier A.4"),
        ("Software, GPU notebook", "pyscf 2.12.1, geometric 1.1.1, ase 3.29.0, skala 2026.9, skala-cuda12x 2026.9, gpu4pyscf-cuda12x 1.8.1, cupy-cuda12x 14.2.0, torch 2.14.1+cu130, huggingface_hub 1.33.0, pyscf-dispersion 1.7.0, python 3.13.16", "dossier A.4"),
    ], columns=["item", "setting", "source"])
    save_table(s6, "tableS6_software_hardware", "Table S6. Software versions, hardware, Skala weights hash, D3 parameters (transcribed from dossier sections A.1, A.3, A.4).")
    # ---------------------------------------------------------------- S7
    ex = pd.read_csv(RES + "excluded_and_failed.csv")
    ex = ex[~ex.system.astype(str).str.contains("grid")]
    save_table(ex, "tableS7_excluded_failed", f"Table S7. Excluded and failed calculations (rigid-grid rows omitted; {len(ex)} rows). Source: R/excluded_and_failed.csv.")
    # ---------------------------------------------------------------- S8
    tl = D["timing_long"]
    save_table(tl, "tableS8_wall_times", "Table S8. Wall time (s) per structure and method on a Tesla T4. Source: W/derived_timing_long.csv.")
    return dict(t1=t1, t2=t2, t3=t3, t4=t4, t5=t5)
