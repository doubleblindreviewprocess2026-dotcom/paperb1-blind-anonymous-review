# ==============================================================================
# b1v5_lit_table.py  --  the literature table LIT for Paper B1 v5 (plan section 4.1).
# Pasted into the check cell (C7) of the notebooks; also importable.
#
# Every value is in kcal/mol.  Columns:
#   id         L1 ... L19 as in the plan
#   system     key of the system database (b1v5_common.SYSTEMS)
#   quantity   what is compared (see QUANTITIES below)
#   value, unc literature value and its stated uncertainty (None = none stated)
#   source     where it comes from
#   zp         zero-point convention of the literature number: "excluded", "included" or "unknown"
#   verified   "yes", "partial" (value yes, convention not) or "no".  ONLY verified == "yes" rows enter error statistics
#              (review 2026-11-15, row 6: L6 CH3NC is "no" until the full paper is read; L19 is a placeholder without a value).
#              L3 (= L1 minus L2) is kept as a row but is not used in pooled statistics (STATS_EXCLUDED_IDS).
#
# QUANTITIES (how the check computes the computed counterpart):
#   barrier_fwd        E(TS) - E(separated reactants)           (single reactant for systems 1 to 4)
#   barrier_rev        E(TS) - E(separated products)
#   barrier_fwd_cplx   E(TS) - E(reactant complex)              (SN2 only)
#   barrier_rev_cplx   E(TS) - E(product complex)               (SN2 only)
#   reaction_energy    E(separated products) - E(separated reactants)
# ==============================================================================
_NHT = "NHTBH38/08 database page, comp.chem.umn.edu/db/dbs/nhtbh38.html (QCISD/MG3 geometries)"

LIT = [
    dict(id="L1", system="hcn", quantity="barrier_fwd", value=48.07, unc=None, source=_NHT, zp="excluded", verified="yes"),
    dict(id="L2", system="hcn", quantity="barrier_rev", value=32.82, unc=None, source=_NHT, zp="excluded", verified="yes"),
    dict(id="L3", system="hcn", quantity="reaction_energy", value=15.25, unc=None, source=_NHT + "; difference of L1 and L2", zp="excluded", verified="yes"),
    dict(id="L4", system="nh3", quantity="barrier_fwd", value=5.12, unc=None,
         source="Lin et al., J. Chem. Phys. 117, 11265 (2002): 1790 cm-1, CCSD(T)/CBS with corrections", zp="excluded", verified="yes"),
    dict(id="L5", system="h2co", quantity="barrier_fwd", value=81.9, unc=0.3,
         source="Feller, Dupuis, Garrett, J. Chem. Phys. 113, 218 (2000), theory (abstract read)", zp="included", verified="yes"),
    dict(id="L6", system="ch3nc", quantity="barrier_fwd", value=38.25, unc=0.25,
         source="Nguyen and Stanton, Sanibel 2018 abstract, HEAT-345Q (experiment 38.2 to 38.4); the abstract calls it an 'activation energy': "
                "may be an Arrhenius-type or 0 K quantity, not an electronic barrier. Full paper (Nguyen, Thorpe, Bross, Ruscic, Stanton, "
                "J. Phys. Chem. Lett. 2018) not yet read", zp="unknown", verified="no"),
    dict(id="L7", system="sn2_ff", quantity="barrier_fwd", value=-0.34, unc=None, source=_NHT, zp="excluded", verified="yes"),
    dict(id="L8", system="sn2_ff", quantity="barrier_fwd_cplx", value=13.38, unc=None, source=_NHT, zp="excluded", verified="yes"),
    dict(id="L9", system="sn2_clcl", quantity="barrier_fwd", value=3.10, unc=None, source=_NHT, zp="excluded", verified="yes"),
    dict(id="L10", system="sn2_clcl", quantity="barrier_fwd_cplx", value=13.41, unc=None, source=_NHT, zp="excluded", verified="yes"),
    dict(id="L11", system="sn2_fcl", quantity="barrier_fwd", value=-12.54, unc=None, source=_NHT, zp="excluded", verified="yes"),
    dict(id="L12", system="sn2_fcl", quantity="barrier_rev", value=20.11, unc=None, source=_NHT, zp="excluded", verified="yes"),
    dict(id="L13", system="sn2_fcl", quantity="barrier_fwd_cplx", value=3.44, unc=None, source=_NHT, zp="excluded", verified="yes"),
    dict(id="L14", system="sn2_fcl", quantity="barrier_rev_cplx", value=29.42, unc=None, source=_NHT, zp="excluded", verified="yes"),
    dict(id="L15", system="sn2_ohf", quantity="barrier_fwd", value=-2.44, unc=None, source=_NHT, zp="excluded", verified="yes"),
    dict(id="L16", system="sn2_ohf", quantity="barrier_rev", value=17.66, unc=None, source=_NHT, zp="excluded", verified="yes"),
    dict(id="L17", system="sn2_ohf", quantity="barrier_fwd_cplx", value=10.96, unc=None, source=_NHT, zp="excluded", verified="yes"),
    dict(id="L18", system="sn2_ohf", quantity="barrier_rev_cplx", value=47.20, unc=None, source=_NHT, zp="excluded", verified="yes"),
    # PLACEHOLDER: the CH3NC -> CH3CN reaction enthalpy is reported in the same paper (abstract only seen); VALUE NOT YET READ.
    dict(id="L19", system="ch3nc", quantity="reaction_energy", value=None, unc=None,
         source="Nguyen, Thorpe, Bross, Ruscic, Stanton, J. Phys. Chem. Lett. 2018 (OSTI 1461519): reaction enthalpy, value not yet read; "
                "convention (0 K or 298 K, enthalpy or energy) unknown", zp="unknown", verified="no"),
]

# Information only: never scored, never joined to a PASS/WARN/FAIL.
LIT_INFO = [
    dict(id="I-a", system="hcn", quantity="HNC minus HCN (0 K)", value=14.89, unc=None,
         source="ATcT, from the v4 literature check (not re-read)", zp="included", verified="no"),
    dict(id="I-b", system="h2co", quantity="barrier_fwd (experiment)", value=79.2, unc=0.8,
         source="Feller, Dupuis, Garrett 2000, experiment", zp="included", verified="yes"),
]

LIT_BY_ID = {r["id"]: r for r in LIT}
STATS_EXCLUDED_IDS = {"L3"}            # L3 = L1 - L2: not independent, never in pooled statistics (only in the supplement table)


def lit_scorable(l):
    """True when the row has a number to compare with (L19 is a placeholder)."""
    return l.get("value") is not None


def lit_in_stats(l):
    """True when the row may enter ANY error statistic: has a value, is verified, and is not the dependent row L3."""
    return lit_scorable(l) and l.get("verified") == "yes" and l["id"] not in STATS_EXCLUDED_IDS
