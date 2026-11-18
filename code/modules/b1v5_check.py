# ==============================================================================
# b1v5_check.py  --  the automated check (cell C7) of Paper B1 v5, plan section 4.
# Pasted into the C7 cell of the notebooks (after b1v5_common, b1v5_lit_table and b1v5_v4_stored); also importable.
#
# What it does, in plain words:
#   1. reads every phase*.csv in the working folder,
#   2. compares our numbers with the literature table LIT (checks R1 to R4),
#   3. runs internal sanity checks I1 to I13,
#   4. prints one line per item: PASS / WARN / FAIL (or NOT RUN when the data is missing),
#   5. writes check_report.csv and check_report.md,
#   6. decides, system by system, where it goes in the paper (rule 4.4) and whether the paper has enough systems.
# It never crashes on missing files: a missing input just gives "NOT RUN".
# All energy differences are in kcal/mol unless written otherwise.
#
# Changes after the final paper-strategy review (2026-11-15), in short:
#   row 3  I6 is scored on wB97M-V (or its fallback) and CCSD(T) only; other functionals give INFO ("no barrier" is a finding).
#          I1 removes a METHOD's entry, not the system, unless a core method (DF-B3LYP, CCSD(T), anchor) fails.
#   row 4  one geometry level per quantity: all species of a barrier must come from the same geometry family, else NOT RUN.
#          Headline family = our own B3LYP geometries; the database geometries are the I11 cross-check.
#   row 5  R1 and the reference statistics use CCSD(T) + [MP2/QZVPP(D) - MP2/TZVP(D)]; the raw CCSD(T) is kept as a column.
#   row 6  unverified literature rows (L6) and placeholders (L19) are never scored PASS or FAIL and never enter statistics.
#   row 7  PASS tolerance for references below 10 kcal/mol is max(0.5, 10% of the reference).
#   row 11 a barrier must not mix density-fitted and unfitted energies (I14).
# ==============================================================================
import os, re, csv, json, glob, math
import numpy as np

try:
    COLUMNS  # noqa  (defined when b1v5_common was pasted above)
except NameError:
    from b1v5_common import *  # noqa
try:
    LIT  # noqa
except NameError:
    from b1v5_lit_table import *  # noqa
try:
    scout_barrier_kcal  # noqa
    HAVE_V4 = True
except NameError:
    try:
        from b1v5_v4_stored import *  # noqa
        HAVE_V4 = True
    except ImportError:
        HAVE_V4 = False

# ------------------------------------------------------------------ tolerances (kcal/mol unless stated); the author may tighten them
TOL = dict(R1=(1.5, 3.0), R2=(2.0, 4.0), R4_flag=5.0, I1_path=(0.98, 0.95), I2=(0.02, 0.03), I5_floor=0.1, I5_peak=(0.5, 1.5),
           I8=(0.5, 1.0), I9_energy=(1e-5, 1e-4), I9_mae=0.02, I9_scout=0.01, I10=0.3, I11_geom=(0.3, 1.0), I11_skala=(0.5, 1.0), I12=(1.0, 3.0))
V4_MAE = {"v4-pbe": 1.242, "v4-b3lyp": 1.374, "v4-skala": 0.564}
V4_METHODS = {"v4-ccsd_t": "ccsd_t", "v4-pbe": "pbe", "v4-b3lyp": "b3lyp", "v4-skala": "skala"}   # method name in rows -> key in b1v5_v4_stored.NB4_HCN_E
SMALL_REF = 10.0           # kcal/mol: for |reference| below this the PASS tolerance is max(REL_FLOOR, REL_FRAC * |reference|) (review row 7)
REL_FLOOR, REL_FRAC = 0.5, 0.10
I6_SCORED = (M_CCSDT_FC,) + tuple(ANCHORS)       # only these methods can FAIL item I6 (review row 3)
NOT_A_FUNCTIONAL = (M_CCSDT_BC, M_MP2DF, M_SKALA_NOD3, M_DFB3LYP, M_CCSDT_AE, M_XTB)   # skipped where "every functional" is looped
PASS, WARN, FAIL, NOTRUN, NA, INFO = "PASS", "WARN", "FAIL", "NOT RUN", "N/A", "INFO"
_RANK = {PASS: 0, NA: 0, INFO: 0, NOTRUN: 1, WARN: 2, FAIL: 3}


def worst(statuses):
    s = [x for x in statuses if x]
    return max(s, key=lambda x: _RANK[x]) if s else NOTRUN


def cap_warn(status):
    return WARN if status == PASS else status


# ------------------------------------------------------------------ loading
_PATH_ID = re.compile(r"^(p\d+|g_\d+_\d+|nb4_\d+)$")


def is_path_id(pid):
    return bool(_PATH_ID.match(str(pid)))


def _f(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except Exception:
        return None


def load_rows():
    """All rows of all phase*.csv, parsed.  Returns (rows, files_found)."""
    rows, files = [], []
    for fn in sorted(glob.glob(path("phase[0-9]*.csv"))):
        files.append(os.path.basename(fn))
        try:
            for r in csv.DictReader(open(fn, newline="")):
                x = extra_of(r)
                rows.append(dict(file=os.path.basename(fn), system=r["system"], pid=r["point_id"], level=r["geometry_level"],
                                 method=r["method"], basis=r["basis"], E=_f(r["energy_hartree"]), conv=r["converged"] in ("1", "True", "true"),
                                 t1=_f(r.get("t1")), imag=(int(float(r["imaginary_modes"])) if _f(r.get("imaginary_modes")) is not None else None),
                                 coord=_f(r.get("coordinate")), x=x))
        except Exception as ex:
            print("  could not read", fn, repr(ex)[:100])
    return rows, files


def level_label(rows):
    """Text for the geometry level of the rows that make up one quantity."""
    lv = sorted({r["level"] for r in rows})
    if lv == [LV_OWN]:
        return LV_OWN
    if set(lv) <= {LV_OWN, LV_SCANPEAK}:
        return LV_OWN + "+scanpeak_TS"
    return "+".join(lv)


class Data:
    """Index of the CSV rows with the lookups the checks need."""
    def __init__(self):
        self.rows, self.files = load_rows()
        cfg = read_json("config.json", {}) or {}
        for k, v in (cfg.get("basis_override") or {}).items():
            BASIS_OVERRIDE[k] = v
        # headline geometry family ('own' by default; 'db' or 'any' only if the CPU notebook was run with OWN_GEOMETRIES_FOR_ALL = False)
        self.headline = cfg.get("headline_family") or "own"
        self.idx = {}
        for r in self.rows:
            self.idx.setdefault((r["system"], r["pid"], r["method"], r["basis"]), []).append(r)
        self.tiers = read_json("tiers.json", {}) or {k: v["tier"] for k, v in SYSTEMS.items()}
        self.excluded = set()          # (system, point, level, basis, method) excluded by I2 / I3
        for r in self.rows:
            if r["method"].startswith("CCSD(T)") and ((r["t1"] or 0) > TOL["I2"][1] or r["x"].get("stable") in (0, "0", False)):
                self.excluded.add((r["system"], r["pid"], r["level"], r["basis"], r["method"]))
        self.mixed = {}                # (system, method, basis) -> dict about density-fitted / unfitted mixing (item I14)
        self.n_derived = 0
        self._build_derived()

    def _build_derived(self):
        """CCSD(T)/TZVP(D) + [DF-MP2/QZVPP(D) - DF-MP2/TZVP(D)] at every geometry where all three energies exist (review row 5).
        The derived rows live only in the index (never in self.rows), so the other checks do not see them."""
        mp2 = {}
        for r in self.rows:
            if r["method"] == M_MP2DF and r["conv"] and r["E"] is not None:
                mp2.setdefault((r["system"], r["pid"], r["level"]), {})[r["basis"]] = r["E"]
        for r in self.rows:
            if r["method"] != M_CCSDT_FC or not r["conv"] or r["E"] is None or r["system"] not in SYSTEMS:
                continue
            mb, lb = basis_for(r["system"], "main"), basis_for(r["system"], "large")
            if r["basis"] != mb or (r["system"], r["pid"], r["level"], r["basis"], r["method"]) in self.excluded:
                continue
            m = mp2.get((r["system"], r["pid"], r["level"]), {})
            if mb in m and lb in m:
                d = m[lb] - m[mb]
                row = dict(r, method=M_CCSDT_BC, E=r["E"] + d, file="derived", x=dict(r["x"], raw_E=r["E"], dmp2_hartree=d, derived=True))
                self.idx.setdefault((r["system"], r["pid"], M_CCSDT_BC, r["basis"]), []).append(row)
                self.n_derived += 1

    def rows_for(self, system, point, method, basis):
        return self.idx.get((system, point, method, basis), [])

    def _levels(self, level):
        """Tuple of geometry levels for a request: None = the headline family, 'own' / 'db' = a family, anything else = that exact level."""
        if level is None:
            return [FAMILY_LEVELS[f] for f in (("own", "db") if self.headline == "any" else (self.headline,))]
        return [FAMILY_LEVELS[level]] if level in FAMILY_LEVELS else [(level,)]

    def _ok_rows(self, system, sp, method, basis, levels, allow_excluded=False):
        rs = [r for r in self.rows_for(system, sp, method, basis) if r["conv"] and r["E"] is not None and r["x"].get("interior") is not False
              and r["level"] in levels]
        if not allow_excluded:
            rs = [r for r in rs if (r["system"], r["pid"], r["level"], r["basis"], r["method"]) not in self.excluded]
        return sorted(rs, key=lambda r: levels.index(r["level"]))

    def energies(self, system, species, method, basis, level=None, allow_excluded=False):
        """All `species` of one quantity at ONE geometry family (review row 4).  Returns (dict species -> row, level label) or None
        when no single family has every species (the quantity is then NOT RUN, never silently mixed)."""
        for levels in self._levels(level):
            got = {}
            for sp in dict.fromkeys(species):
                rs = self._ok_rows(system, sp, method, basis, levels, allow_excluded)
                if not rs:
                    break
                got[sp] = rs[0]
            else:
                return got, level_label(list(got.values()))
        return None

    def aligned(self, system, sp, method, basis, level, df=True):
        """Partner energy recomputed with the same retry mode (pid 'SP@df'), or None."""
        rs = self._ok_rows(system, sp + ("@df" if df else ""), method, basis, (level,))
        return rs[0] if rs else None

    def E(self, system, sp, method, basis, level=None, allow_excluded=False):
        """Energy (hartree) of ONE stationary species (own geometry first, then the database one); None if absent.
        Anything that combines species must use energies() or quantity(), which keep the geometry level fixed."""
        lv = self._levels(level) if level else [tuple(LEVEL_PREFERENCE)]
        for levels in lv:
            rs = self._ok_rows(system, sp, method, basis, levels, allow_excluded)
            if rs:
                return rs[0]["E"]
        return None

    def zpe(self, system, sp):
        for lv in (LV_OWN, LV_SCANPEAK):
            for r in self.rows_for(system, sp, M_DFB3LYP, basis_for(system, "main")):
                if r["level"] == lv and r["x"].get("zpe_hartree") is not None:
                    return float(r["x"]["zpe_hartree"])
        return None

    def methods_with(self, system, basis):
        return sorted({k[2] for k in self.idx if k[0] == system and k[3] == basis and not is_path_id(k[1]) and "@" not in k[1]})


def _quantity_species(system, q):
    """(plus species, minus species) of a quantity; None when the quantity does not exist for this system."""
    S = SYSTEMS[system]
    if q == "barrier_fwd":
        return [S["ts"]] if S["ts"] else None, S["reactants"]
    if q == "barrier_rev":
        return [S["ts"]] if S["ts"] else None, S["products"]
    if q == "barrier_fwd_cplx":
        return ([S["ts"]] if S["ts"] else None), ([S["rc"]] if S["rc"] else [])
    if q == "barrier_rev_cplx":
        return ([S["ts"]] if S["ts"] else None), ([S["pc"]] if S["pc"] else [])
    if q == "reaction_energy":
        return S["products"], S["reactants"]
    return None, None


def quantity_info(D, system, q, method, basis, zpe=False, level=None):
    """Computed counterpart of a LIT quantity.  Returns dict(value in kcal/mol, level label, aligned) or None when anything is missing,
    a species comes from another geometry family, or fitted and unfitted energies would be mixed without a recomputed partner."""
    plus, minus = _quantity_species(system, q)
    if not plus or not minus:
        return None
    if zpe:
        z = [D.zpe(system, sp) for sp in plus + minus]
        if any(v is None for v in z):
            return None
        return dict(value=(sum(z[:len(plus)]) - sum(z[len(plus):])) * HARTREE_KCAL, level="own_b3lyp", aligned=False)
    got = D.energies(system, plus + minus, method, basis, level)
    if got is None:
        return None
    rows, lab = got
    aligned = False
    dfs = {sp: bool(r["x"].get("df")) for sp, r in rows.items()}
    if len(set(dfs.values())) > 1:                       # review row 11: fitted and unfitted energies in one barrier
        fit_levels = {sp: rows[sp]["level"] for sp in rows}
        new = {}
        for sp, f in dfs.items():
            if not f:
                a = D.aligned(system, sp, method, basis, fit_levels[sp], df=True)
                if a is None:
                    D.mixed[(system, method, basis)] = dict(species=sorted(sp for sp, f2 in dfs.items() if f2), partners_missing=sorted(sp for sp, f2 in dfs.items() if not f2), aligned=False)
                    return None
                new[sp] = a
        rows = dict(rows, **new); aligned = True
        D.mixed[(system, method, basis)] = dict(species=sorted(sp for sp, f2 in dfs.items() if f2), partners_missing=[], aligned=True)
    e = sum(rows[sp]["E"] for sp in plus) - sum(rows[sp]["E"] for sp in minus)
    return dict(value=e * HARTREE_KCAL, level=lab, aligned=aligned)


def quantity(D, system, q, method, basis, zpe=False, level=None, **kw):
    """Value (kcal/mol) of a quantity or None.  See quantity_info."""
    i = quantity_info(D, system, q, method, basis, zpe=zpe, level=level)
    return None if i is None else i["value"]


def reference_quantity(D, system, q, basis):
    """The CCSD(T) reference of a quantity: basis-corrected when the MP2 rows exist (review row 5), else raw.
    Returns dict(value, raw, corrected: bool, level) or None."""
    c = quantity_info(D, system, q, M_CCSDT_BC, basis)
    r = quantity_info(D, system, q, M_CCSDT_FC, basis)
    if c is not None:
        return dict(value=c["value"], raw=None if r is None else r["value"], corrected=True, level=c["level"])
    if r is not None:
        return dict(value=r["value"], raw=r["value"], corrected=False, level=r["level"])
    return None


def compare_to_lit(D, lit, method, basis, level=None):
    """Computed value after the zero-point rule of plan 4.1.  Returns dict(value, diff, note, forced_warn, level) or None."""
    if lit.get("value") is None:
        return None
    qi = quantity_info(D, lit["system"], lit["quantity"], method, basis, level=level)
    if qi is None:
        return None
    el, lab = qi["value"], qi["level"]
    zp = lit["zp"]
    dz = quantity(D, lit["system"], lit["quantity"], method, basis, zpe=True)
    base = dict(level=lab, aligned=qi["aligned"])
    if zp == "excluded":
        return dict(base, value=el, diff=el - lit["value"], note="electronic energy vs zero-point-excluded reference", forced_warn=False)
    if zp == "included":
        if dz is None:
            return dict(base, value=None, diff=None, note="needs the B3LYP zero-point energies from C1", forced_warn=False, missing_zpe=True)
        v = el + dz
        return dict(base, value=v, diff=v - lit["value"], note=f"zero-point difference {dz:+.2f} added", forced_warn=False)
    # unknown convention (review row 6): the electronic value is compared, nothing is picked to look better; never better than WARN
    return dict(base, value=el, diff=el - lit["value"], note="zero-point convention unknown: electronic value shown, not scored", forced_warn=True)


def pass_tolerance(tol, ref):
    """PASS limit; for small references (|ref| below 10 kcal/mol) at most max(0.5, 10% of |ref|)  (review row 7)."""
    if ref is not None and abs(ref) < SMALL_REF:
        return min(tol[0], max(REL_FLOOR, REL_FRAC * abs(ref)))
    return tol[0]


def grade(diff, tol, ref=None):
    d = abs(diff)
    return PASS if d <= pass_tolerance(tol, ref) else (WARN if d <= tol[1] else FAIL)


# ------------------------------------------------------------------ report container
class Report:
    def __init__(self):
        self.items = []

    def add(self, item, scope, status, value="", ref="", diff="", tol="", msg="", raw_value="", raw_difference=""):
        self.items.append(dict(item=item, scope=scope, status=status, value=value, reference=ref, difference=diff, tolerance=tol, message=msg,
                               raw_value=raw_value, raw_difference=raw_difference))
        fmt = lambda v: f"{v:.3f}" if isinstance(v, float) else str(v)
        print(f"{status:8s} {item:4s} {scope:22s} value={fmt(value):>9s} ref={fmt(ref):>9s} diff={fmt(diff):>8s} tol={tol!s:>9s}  {msg}")

    def status(self, item, scope):
        return [x["status"] for x in self.items if x["item"] == item and sys_of(x["scope"]) == scope]

    def counts(self):
        c = {k: 0 for k in (PASS, WARN, FAIL, NOTRUN, NA, INFO)}
        for x in self.items:
            c[x["status"]] += 1
        return c


def sys_of(scope):
    """The system key inside an item scope such as 'L1 hcn:barrier_fwd' or 'nh3 CCSD(T)'."""
    for tok in re.split(r"[ :]", scope):
        if tok in SYSTEMS:
            return tok
    return ""


def _fnum(x, n=2):
    return "" if x is None else round(float(x), n)


# ------------------------------------------------------------------ literature checks
R2_CAUSE = {"sn2": "missing diffuse functions, or wrong fragment reference (complex vs separated)", "h2co": "zero-point handling",
            "ch3nc": "zero-point handling (convention unknown)", "other": "geometry level or functional string"}


def literature_checks(D, R):
    for lit in LIT:
        s = lit["system"]; tier = D.tiers.get(s, SYSTEMS[s]["tier"]); mb, lb = basis_for(s, "main"), basis_for(s, "large")
        tag = f"{lit['id']} {s}:{lit['quantity']}"
        if not lit_scorable(lit):                     # placeholder row (L19): nothing to compare with yet
            R.add("R3", tag, WARN, "", "", "", "never PASS", f"PLACEHOLDER, value not yet read, excluded from all statistics. Source: {lit['source'][:90]}")
            continue
        unverified = lit["verified"] != "yes"
        # --- R1 protocol: our CCSD(T) (basis-corrected) vs literature (tier A and B only)
        if tier not in ("A", "B"):
            R.add("R1", tag, NA, msg="tier C: no CCSD(T) planned for this system")
        else:
            c = compare_to_lit(D, lit, M_CCSDT_BC, mb); corrected = True
            raw = compare_to_lit(D, lit, M_CCSDT_FC, mb)
            if c is None or c.get("value") is None:
                c, corrected = raw, False
            if c is None or c.get("value") is None:
                R.add("R1", tag, NOTRUN, msg="no usable CCSD(T) rows at one geometry level" + (" (zero-point energies missing)" if c is not None else " (excluded by I2/I3, missing, or species at different levels)"))
            else:
                st = grade(c["diff"], TOL["R1"], lit["value"])
                if unverified or c["forced_warn"]:
                    st = WARN                                    # never PASS, and never FAIL, on an unverified reference
                elif not corrected and st == PASS:
                    st = WARN                                    # no PASS without the MP2 basis correction
                msg = c["note"] + f"; geometry level {c['level']}"
                msg += ("; reference = CCSD(T) + [MP2/" + basis_for(s, "large") + " - MP2/" + mb + "]") if corrected else "; MP2 basis correction NOT available: raw CCSD(T) used (no PASS possible)"
                if unverified:
                    msg += "; reference unverified: not scored"
                if st == WARN and not unverified:
                    msg += "; likely cause: remaining basis incompleteness, zero-point convention, or geometry level"
                if st == FAIL:
                    msg += "; RISK: a referee will doubt every reference number. FALLBACK: literature becomes the only reference for this system, CCSD(T) path errors go to the supplement"
                R.add("R1", tag, st, _fnum(c["value"]), lit["value"], _fnum(c["diff"]), f"{TOL['R1'][0]}/{TOL['R1'][1]}" + (f" (PASS {pass_tolerance(TOL['R1'], lit['value']):.2f})" if pass_tolerance(TOL['R1'], lit['value']) != TOL['R1'][0] else ""), msg,
                      raw_value=_fnum(raw["value"]) if raw and raw.get("value") is not None else "", raw_difference=_fnum(raw["diff"]) if raw and raw.get("value") is not None else "")
        # --- R2 anchor: wB97M-V (fallback wB97X-V) in the large basis vs literature (all systems)
        c, used, anchor = None, lb, M_WB97MV
        for anchor in ANCHORS:
            for b in (lb, mb):
                c = compare_to_lit(D, lit, anchor, b)
                if c is not None and c.get("value") is not None:
                    used = b; break
            if c is not None and c.get("value") is not None:
                break
        if c is None or c.get("value") is None:
            R.add("R2", tag, NOTRUN, msg="no wB97M-V (or fallback wB97X-V) rows at one geometry level (GPU notebook G1/G2 not run yet)")
        else:
            st = grade(c["diff"], TOL["R2"], lit["value"])
            if unverified or c["forced_warn"]:
                st = WARN if st in (PASS, FAIL) else st
            msg = c["note"] + f"; geometry level {c['level']}" + ("" if used == lb else f"; large basis {lb} not available, scored at {used}")
            if anchor != M_WB97MV:
                msg += f"; FALLBACK ANCHOR {anchor} used (wB97M-V does not run)"
            if unverified:
                msg += "; reference unverified: not scored"
            if st == WARN and not unverified:
                msg += "; likely cause: " + R2_CAUSE["sn2" if SYSTEMS[s]["sn2"] else (s if s in R2_CAUSE else "other")]
            if st == FAIL:
                msg += "; geometry or setup error suspected. FALLBACK: system leaves the main text (rule 4.4)"
            R.add("R2", tag, st, _fnum(c["value"]), lit["value"], _fnum(c["diff"]), f"{TOL['R2'][0]}/{TOL['R2'][1]}" + (f" (PASS {pass_tolerance(TOL['R2'], lit['value']):.2f})" if pass_tolerance(TOL['R2'], lit['value']) != TOL['R2'][0] else ""), msg)
        # --- R3 unverified reference
        if unverified:
            R.add("R3", tag, WARN, lit["value"], lit["value"], "", "never PASS",
                  f"verified={lit['verified']} ({lit['zp']} zero-point): EXCLUDED from all error statistics; confirm source before submission. Source: {lit['source'][:90]}")
        # --- R4 Skala (with and without D3) and the other functionals vs literature (reported, not scored)
        errs = {}
        for m in LADDER_ALL:
            c = compare_to_lit(D, lit, m, mb)
            if c is not None and c.get("value") is not None:
                errs[m] = c["diff"]
        if not errs:
            R.add("R4", tag, NOTRUN, msg="no functional rows (GPU notebook G1 not run yet)")
        else:
            sk = errs.get(M_SKALA)
            msg = "signed errors (method minus literature): " + ", ".join(f"{m} {e:+.2f}" for m, e in errs.items()) + ("; reference unverified, not in statistics" if unverified else "")
            if sk is not None and abs(sk) > TOL["R4_flag"]:
                R.add("R4", tag, WARN, _fnum(sk), lit["value"], _fnum(sk), f">{TOL['R4_flag']}", "Skala error is large: re-run that point from a fresh guess. " + msg)
            else:
                R.add("R4", tag, INFO, _fnum(sk), lit["value"], _fnum(sk), "not scored", msg)


# ------------------------------------------------------------------ internal checks
def _path(D, system, method, basis):
    rs = [r for k, v in D.idx.items() if k[0] == system and k[2] == method and k[3] == basis and is_path_id(k[1]) for r in v]
    rs = sorted([r for r in rs if r["conv"] and r["E"] is not None], key=lambda r: int(re.sub(r"\D", "", r["pid"]) or 0))
    return rs


def i1_scf(D, R):
    """SCF convergence.  A failed point removes THAT METHOD's entry (it is already left out of every quantity); it fails the SYSTEM only when
    a core method (DF-B3LYP geometries, CCSD(T), the anchor functional) is affected (review row 3)."""
    for s in SYSTEM_ORDER:
        rs = [r for r in D.rows if r["system"] == s and r["method"] != M_XTB and "@" not in r["pid"]]
        if not rs:
            R.add("I1", s, NOTRUN, msg="no rows"); continue
        good = {(r["pid"], r["method"], r["basis"], r["level"]) for r in rs if r["conv"]}          # a failed row that was later repeated and converged does not count
        rs = [r for r in rs if r["conv"] or (r["pid"], r["method"], r["basis"], r["level"]) not in good]
        stat = [r for r in rs if not is_path_id(r["pid"])]; pth = [r for r in rs if is_path_id(r["pid"])]
        bad_core = [f"{r['pid']}/{r['method']}" for r in stat if not r["conv"] and r["method"] in CORE_METHODS]
        bad_other = [f"{r['pid']}/{r['method']}" for r in stat if not r["conv"] and r["method"] not in CORE_METHODS]
        frac_core, frac_other = {}, {}
        for m in sorted({r["method"] for r in pth}):
            xs = [r["conv"] for r in pth if r["method"] == m]
            (frac_core if m in CORE_METHODS else frac_other)[m] = sum(xs) / len(xs)
        worst_core = min(frac_core.values()) if frac_core else None
        weak_other = {m: f for m, f in frac_other.items() if f < TOL["I1_path"][0]}
        frac = worst_core
        if bad_core:
            R.add("I1", s, FAIL, len(bad_core), "", "", "0 unconverged (core methods)", "unconverged stationary points of a core method: " + ", ".join(bad_core[:6]))
        elif worst_core is not None and worst_core < TOL["I1_path"][1]:
            R.add("I1", s, FAIL, _fnum(worst_core, 3), "", "", ">=0.95", f"only {worst_core:.1%} of path points converged for a core method")
        elif bad_other or weak_other or (worst_core is not None and worst_core < TOL["I1_path"][0]):
            msg = []
            if bad_other:
                msg.append("method entries removed (unconverged stationary points): " + ", ".join(bad_other[:8]))
            if weak_other:
                msg.append("path convergence below 98% for " + ", ".join(f"{m} {f:.0%}" for m, f in weak_other.items()))
            if worst_core is not None and worst_core < TOL["I1_path"][0]:
                msg.append(f"core method path convergence {worst_core:.1%}")
            R.add("I1", s, WARN, _fnum(frac, 3) if frac is not None else len(bad_other), "", "", ">=0.98 path, all stationary",
                  "; ".join(msg) + ". Only these method entries are dropped, the system stays (rule 4.4)")
        else:
            R.add("I1", s, PASS, _fnum(frac, 3) if frac is not None else "", "", "", ">=0.98 path, all stationary",
                  f"{len(stat)} stationary rows" + (f", {len(pth)} path rows" if pth else ""))


def i2_i3_cc(D, R):
    for s in SYSTEM_ORDER:
        rs = [r for r in D.rows if r["system"] == s and r["method"].startswith("CCSD(T)") and r["file"] != "phase06_ccsdt_basis.csv" and r["level"] != LV_XTB]
        if not rs:
            R.add("I2", s, NA if D.tiers.get(s) == "C" else NOTRUN, msg="no CCSD(T) rows" + (" (tier C)" if D.tiers.get(s) == "C" else ""))
            R.add("I3", s, NA if D.tiers.get(s) == "C" else NOTRUN, msg="no CCSD(T) rows")
            continue
        t1 = [(r["t1"], r["pid"]) for r in rs if r["t1"] is not None]
        if not t1:
            R.add("I2", s, NOTRUN, msg="no T1 values")
        else:
            mx = max(t1)
            flagged = [f"{p}({v:.3f})" for v, p in t1 if v >= TOL["I2"][0]]
            st = PASS if mx[0] < TOL["I2"][0] else (WARN if mx[0] <= TOL["I2"][1] else FAIL)
            msg = f"max T1 {mx[0]:.4f} at {mx[1]}" + (f"; flagged: {', '.join(flagged[:6])}" if flagged else "")
            if st == FAIL:
                msg += "; point(s) excluded from reference statistics"
            R.add("I2", s, st, _fnum(mx[0], 4), "", "", f"<{TOL['I2'][0]} / {TOL['I2'][1]}", msg)
        stab = [(r["x"].get("stable"), r["pid"]) for r in rs if r["x"].get("stable") not in (None, "")]
        if not stab:
            R.add("I3", s, NOTRUN, msg="no stability results")
        else:
            bad = [p for v, p in stab if str(v) in ("0", "False")]
            R.add("I3", s, FAIL if bad else PASS, len(bad), "", "", "all stable", ("unstable RHF, points excluded: " + ", ".join(bad[:6])) if bad else f"{len(stab)} points stable")


def i4_imag(D, R):
    for s in SYSTEM_ORDER:
        sts, msgs, n_own, n_db = [], [], 0, 0
        for sp, kind in SPECIES_KIND[s].items():
            own = [r for r in D.rows if r["system"] == s and r["pid"] == sp and r["method"] == M_DFB3LYP and r["level"] in (LV_OWN, LV_SCANPEAK) and r["imag"] is not None]
            if own:
                n_own += 1
                want = 1 if kind == "ts" else 0
                if own[0]["imag"] != want:
                    sts.append(FAIL); msgs.append(f"{sp}: {own[0]['imag']} imaginary modes, expected {want}")
            elif any(r["system"] == s and r["pid"] == sp and r["level"] == LV_DB for r in D.rows):
                n_db += 1
        if n_own == 0 and n_db == 0:
            R.add("I4", s, NOTRUN if s != "ch4" else NA, msg="no frequency data" if s != "ch4" else "single minimum; not needed"); continue
        if not sts:
            if n_db and not n_own:
                sts.append(WARN); msgs.append("database geometry, not tested")
            else:
                sts.append(PASS); msgs.append(f"{n_own} own stationary points tested" + (f"; {n_db} database geometries not tested" if n_db else ""))
        R.add("I4", s, worst(sts), "", "", "", "TS 1, minima 0", "; ".join(msgs))


def _e_pair(D, s, sps, m, basis, level=None):
    """Energies of the species `sps` at one geometry family: dict species -> hartree, or None."""
    g = D.energies(s, sps, m, basis, level)
    return None if g is None else {k: r["E"] for k, r in g[0].items()}


def i5_i6(D, R):
    for s in SYSTEM_ORDER:
        S = SYSTEMS[s]
        mb = basis_for(s, "main")
        meths = [m for m in D.methods_with(s, mb) if m not in (M_CCSDT_BC, M_MP2DF, M_SKALA_NOD3, M_XTB)]
        # ---- I6: transition state above both minima.  SCORED on CCSD(T) and the anchor functional only (review row 3);
        #      for every other functional "no barrier" is a finding (INFO), never a reason to drop the system.
        if S["ts"]:
            lows = [x for x in dict.fromkeys(([S["rc"], S["pc"]] if S["sn2"] else S["reactants"][:1] + (S["products"] if len(S["products"]) == 1 and S["products"] != S["reactants"] else []))) if x]
            sum_prod = (not S["sn2"]) and len(S["products"]) > 1
            bad, ok, info = [], 0, []
            for m in meths:
                comps = [([S["ts"], sp], sp, lambda e, sp=sp: e[sp]) for sp in lows]
                if sum_prod:
                    comps.append(([S["ts"]] + S["products"], "sum of products", lambda e: sum(e[q] for q in S["products"])))
                for need, label, f in comps:
                    e = _e_pair(D, s, need, m, mb)
                    if e is None:
                        continue
                    below = e[S["ts"]] <= f(e)
                    if m in I6_SCORED:
                        ok += 1
                        if below:
                            bad.append(f"{m} (TS {(e[S['ts']] - f(e)) * HARTREE_KCAL:+.1f} vs {label})")
                    elif below:
                        info.append((m, label, (e[S["ts"]] - f(e)) * HARTREE_KCAL))
            R.add("I6", s, FAIL if bad else (PASS if ok else NOTRUN), len(bad), "", "", "TS above minima (CCSD(T), anchor)",
                  ("TS not above minimum: " + "; ".join(bad[:4])) if bad else (f"{ok} CCSD(T)/anchor comparisons all fine" if ok else "no CCSD(T) or anchor energies yet at one geometry level"))
            for m, label, d in info:
                R.add("I6", f"{s} {m}", INFO, _fnum(d), "", "", "not scored", f"no barrier with {m}: TS {d:+.1f} kcal/mol vs {label}; a finding for the paper, the system stays")
        else:
            R.add("I6", s, NA, msg="no transition state (bond stretch)")
        # ---- I5: path single-peaked, peak close to the TS energy
        if not S["ts"]:
            R.add("I5", s, NA, msg="no reaction path peak expected (bond stretch)"); continue
        done = False
        for m in [M_CCSDT_FC, M_WB97MV, M_B3LYP, M_DFB3LYP, M_PBE]:
            pts = _path(D, s, m, mb)
            if len(pts) < 5:
                continue
            r0 = S["reactants"][0] if not S["sn2"] else S["rc"]
            ee = _e_pair(D, s, [r0, S["ts"]], m, mb)
            e0 = ee[r0] if ee else (_e_pair(D, s, [r0], m, mb) or {}).get(r0)
            e0 = e0 if e0 is not None else pts[0]["E"]
            E = (np.array([r["E"] for r in pts]) - e0) * HARTREE_KCAL      # relative to the reactant minimum
            d = np.diff(E); d[np.abs(d) < TOL["I5_floor"]] = 0
            sg = np.sign(d[d != 0]); n_max = int(np.sum((sg[:-1] > 0) & (sg[1:] < 0))) if len(sg) > 1 else 0
            ets = ee[S["ts"]] if ee else None
            diff = None if ets is None else float(E.max() - (ets - e0) * HARTREE_KCAL)
            if n_max != 1:
                st, msg = FAIL, f"{n_max} maxima along {len(pts)} {m} path points (need exactly one)"
            elif diff is None:
                st, msg = WARN, "one maximum, but no TS/minimum energy to compare with"
            else:
                a = abs(diff); st = PASS if a <= TOL["I5_peak"][0] else (WARN if a <= TOL["I5_peak"][1] else FAIL)
                msg = f"one maximum ({len(pts)} {m} points); path peak minus TS energy {diff:+.2f}"
            R.add("I5", s, st, _fnum(diff), "", _fnum(diff), f"{TOL['I5_peak'][0]}/{TOL['I5_peak'][1]}", msg)
            done = True
            break
        if not done:
            R.add("I5", s, NOTRUN, msg="fewer than 5 path points in any method")


def i7_rank(D, R):
    mb = basis_for("hcn", "main")
    b = {m: quantity(D, "hcn", "barrier_fwd", m, mb) for m in (M_PBE, M_B3LYP, M_CCSDT_FC)}
    if any(v is None for v in b.values()):
        R.add("I7", "hcn", NOTRUN, msg="needs PBE, B3LYP and CCSD(T) at the HCN stationary points")
    else:
        ok = b[M_PBE] < b[M_B3LYP] < b[M_CCSDT_FC]
        R.add("I7", "hcn", PASS if ok else WARN, "", "", "", "PBE < B3LYP < CCSD(T)", f"barriers PBE {b[M_PBE]:.2f}, B3LYP {b[M_B3LYP]:.2f}, CCSD(T) {b[M_CCSDT_FC]:.2f}" + ("" if ok else "; order differs: check functional strings and dispersion flags"))
    for s in [k for k in SYSTEMS if SYSTEMS[k]["sn2"]]:
        mb = basis_for(s, "main")
        v = {m: quantity(D, s, "barrier_fwd_cplx", m, mb) for m in LADDER}
        if any(x is None for x in v.values()):
            R.add("I7", s, NOTRUN, msg="needs all five functionals"); continue
        ok = min(v, key=v.get) == M_PBE
        R.add("I7", s, PASS if ok else WARN, "", "", "", "PBE lowest", "central barriers: " + ", ".join(f"{m} {x:.2f}" for m, x in v.items()) + ("" if ok else "; PBE is not lowest: check functional strings and dispersion flags"))


def i8_basis(D, R):
    for s in REACTIONS:
        q = "barrier_fwd_cplx" if SYSTEMS[s]["sn2"] else "barrier_fwd"
        mb, lb = basis_for(s, "main"), basis_for(s, "large")
        sh = {}
        for m in sorted({k[2] for k in D.idx if k[0] == s} - {M_CCSDT_BC, M_MP2DF, M_XTB}):
            a, b = quantity(D, s, q, m, mb), quantity(D, s, q, m, lb)
            if a is not None and b is not None:
                sh[m] = b - a
        if not sh:
            R.add("I8", s, NOTRUN, msg=f"no method has the barrier in both {mb} and {lb}"); continue
        w = max(sh.values(), key=abs)
        st = grade(w, TOL["I8"])
        msg = f"{mb} to {lb} shifts: " + ", ".join(f"{m} {x:+.2f}" for m, x in sh.items())
        if st == WARN:
            msg += "; quote the larger basis as the headline"
        if st == FAIL:
            msg += "; this system is reported at the large basis only"
        R.add("I8", s, st, _fnum(w), "", "", f"{TOL['I8'][0]}/{TOL['I8'][1]}", msg)


def i9_v4(D, R):
    if not HAVE_V4:
        R.add("I9", "v4", NOTRUN, msg="stored v4 numbers (b1v5_v4_stored) not available"); return
    # (a) energies of NB4 point(s) vs stored
    worst_d, n, per = 0.0, 0, []
    for m, key in V4_METHODS.items():
        for r in D.rows:
            if r["system"] == "hcn_nb4" and r["method"] == m and r["E"] is not None:
                k = int(r["pid"].split("_")[-1])
                d = abs(r["E"] - NB4_HCN_E[key][k]); n += 1; worst_d = max(worst_d, d); per.append(f"{m}[{k}] {d:.1e}")
    if n == 0:
        R.add("I9", "v4 energies", NOTRUN, msg="no rerun of stored NB4 points yet (C0 does CCSD(T); GPU G0 does the DFT methods)")
    else:
        st = PASS if worst_d <= TOL["I9_energy"][0] else (WARN if worst_d <= TOL["I9_energy"][1] else FAIL)
        msg = f"{n} stored energies rerun, largest |difference| {worst_d:.1e} hartree"
        if st == WARN:
            msg += "; version drift, state it"
        if st == FAIL:
            msg += "; STOP: do not mix old and new numbers"
        R.add("I9", "v4 energies", st, f"{worst_d:.1e}", "", "", f"{TOL['I9_energy'][0]:g}/{TOL['I9_energy'][1]:g} hartree", msg)
    # (b) path mean errors recomputed from all 12 points
    ref = np.array(NB4_HCN_E["ccsd_t"]); rel = lambda E: (np.array(E) - E[0]) * HARTREE_KCAL
    res = []
    for m, mae in V4_MAE.items():
        pts = {int(r["pid"].split("_")[-1]): r["E"] for r in D.rows if r["system"] == "hcn_nb4" and r["method"] == m and r["E"] is not None}
        if len(pts) == 12:
            E = [pts[k] for k in range(12)]
            res.append((m, float(np.mean(np.abs(rel(E) - rel(ref)))), mae))
    if not res:
        R.add("I9", "v4 path errors", NOTRUN, msg="needs all 12 NB4 HCN path points for PBE, B3LYP and Skala (GPU G0 extended)")
    else:
        st = worst([PASS if abs(a - b) <= TOL["I9_mae"] else FAIL for _, a, b in res])
        R.add("I9", "v4 path errors", st, "", "", "", f"{TOL['I9_mae']}", ", ".join(f"{m} {a:.3f} vs {b:.3f}" for m, a, b in res))
    # (c) scout values from the stored arrays
    out = []
    for s, stored in SCOUT_STORED.items():
        v = scout_barrier_kcal(s); out.append((s, v, stored))
    st = worst([PASS if abs(v - t) <= TOL["I9_scout"] else FAIL for _, v, t in out])
    R.add("I9", "v4 scout values", st, "", "", "", f"{TOL['I9_scout']}", ", ".join(f"{s} {v:.3f} vs {t:.3f}" for s, v, t in out))


def i10_fc(D, R):
    out = []
    for s in REACTIONS:
        mb = basis_for(s, "main")
        a, b = quantity(D, s, "barrier_fwd", M_CCSDT_AE, mb), quantity(D, s, "barrier_fwd", M_CCSDT_FC, mb)
        if a is not None and b is not None:
            out.append((s, a - b))
    if not out:
        R.add("I10", "all", NOTRUN, msg="no system has both all-electron and frozen-core CCSD(T) (phase C6)")
    for s, d in out:
        st = PASS if abs(d) <= TOL["I10"] else WARN
        R.add("I10", s, st, _fnum(d), "", "", f"{TOL['I10']}", f"all-electron minus frozen-core barrier {d:+.2f}" + ("" if st == PASS else "; larger: report both"))


def i11_geom(D, R):
    """Own B3LYP geometries (headline) against the published database geometries (cross-check), for every reaction that has both
    (review row 4), then Skala//Skala against Skala//reference."""
    for s in REACTIONS:
        q = "barrier_fwd_cplx" if SYSTEMS[s]["sn2"] else "barrier_fwd"
        mb = basis_for(s, "main")
        if not SYSTEMS[s]["db"]:
            continue
        res = None
        for m in (M_CCSDT_FC,) + ANCHORS:
            basis = mb
            a, b = quantity(D, s, q, m, basis, level="own"), quantity(D, s, q, m, basis, level="db")
            if a is not None and b is not None:
                res = (m, a, b); break
        if res is None:
            R.add("I11", f"{s} own vs database", NOTRUN, msg="needs the barrier at both our own and the downloaded database geometries (database not downloaded, or G1 not run at both levels?)")
        else:
            m, a, b = res
            st = grade(a - b, (TOL["I11_geom"][0], TOL["I11_geom"][1]))
            R.add("I11", f"{s} own vs database", st, _fnum(a), _fnum(b), _fnum(a - b), f"{TOL['I11_geom'][0]}/{TOL['I11_geom'][1]}",
                  f"own minus database {m} barrier" + ("; report both barriers" if st == FAIL else ""))
    any_ = False
    for s in REACTIONS:
        mb = basis_for(s, "main")
        q = "barrier_fwd_cplx" if SYSTEMS[s]["sn2"] else "barrier_fwd"
        sk = quantity(D, s, q, M_SKALA, mb, level=LV_SKALA)
        rf = quantity(D, s, q, M_SKALA, mb)
        if sk is not None and rf is not None:
            any_ = True
            st = grade(sk - rf, TOL["I11_skala"])
            R.add("I11", f"{s} Skala//Skala", st, _fnum(sk), _fnum(rf), _fnum(sk - rf), f"{TOL['I11_skala'][0]}/{TOL['I11_skala'][1]}", "Skala//Skala (TS = parabola vertex) minus Skala//reference barrier" + ("; report both barriers" if st == FAIL else ""))
    if not any_:
        R.add("I11", "Skala//Skala", NOTRUN, msg="GPU phase G4 not run yet")


def i12_scout(D, R):
    for s in REACTIONS:
        S = SYSTEMS[s]; mb = basis_for(s, "main"); rmin = S["rc"] or S["reactants"][0]
        xrows = [r for r in D.rows if r["system"] == s and r["level"] == LV_XTB and r["pid"] == S["ts"]]
        if any((not r["conv"]) or r["x"].get("no_peak") for r in xrows if r["method"] == M_XTB):
            R.add("I12", s, WARN, msg="xTB found no peak along the scan: recorded as a result (outside the recipe's validity)"); continue
        for m, tag in ((M_CCSDT_FC, "CCSD(T)"), (M_WB97MV, "wB97M-V"), (M_WB97XV, "wB97X-V")):
            x = _e_pair(D, s, [S["ts"], rmin], m, mb, level=LV_XTB)
            r = _e_pair(D, s, [S["ts"], rmin], m, mb)
            if x is None or r is None:
                continue
            xb, rb = (x[S["ts"]] - x[rmin]) * HARTREE_KCAL, (r[S["ts"]] - r[rmin]) * HARTREE_KCAL
            pen = xb - rb
            st = grade(pen, TOL["I12"])
            R.add("I12", f"{s} {tag}", st, _fnum(xb), _fnum(rb), _fnum(pen), f"{TOL['I12'][0]}/{TOL['I12'][1]}",
                  "scout (xTB geometries) minus reference-geometry barrier" + ("; system outside the recipe's validity (a result)" if st == FAIL else ""))
            break
        else:
            R.add("I12", s, NOTRUN, msg="no CCSD(T) or wB97M-V energies at xTB geometries (C4 / G3)")


def i13_spread(D, R):
    for s in REACTIONS:
        q = "barrier_fwd_cplx" if SYSTEMS[s]["sn2"] else "barrier_fwd"; mb = basis_for(s, "main")
        v = {m: quantity(D, s, q, m, mb) for m in LADDER}
        std = [v[m] for m in (M_PBE, M_R2SCAN, M_B3LYP, M_WB97MV) if v[m] is not None]
        if len(std) < 2:
            R.add("I13", s, NOTRUN, msg="needs at least two standard functionals"); continue
        lo, hi = min(std), max(std); sk = v[M_SKALA]
        where = "no Skala value" if sk is None else ("inside" if lo <= sk <= hi else "outside")
        R.add("I13", s, INFO, _fnum(hi - lo), "", "", "information only", f"standard functionals {lo:.2f} to {hi:.2f} (range {hi - lo:.2f}); Skala {'' if sk is None else f'{sk:.2f} '}{where}")


def i14_retry_mode(D, R):
    """Review row 11: a barrier must not mix density-fitted (retry) and unfitted energies.  The GPU notebook recomputes the partner species
    with density fitting (pid 'SP@df'); the quantity then uses that partner.  Without it the barrier is dropped from every table."""
    n_rows = sum(1 for r in D.rows if r["file"] != "derived" and r["x"].get("route") is not None)
    if not n_rows:
        R.add("I14", "all", NOTRUN, msg="no GPU rows yet"); return
    n_df = sum(1 for r in D.rows if r["x"].get("df") and "@" not in r["pid"])
    for s in REACTIONS:
        mb = basis_for(s, "main")
        for m in LADDER_ALL + [M_WB97XV]:
            for q in {"barrier_fwd", "barrier_fwd_cplx"}:
                quantity_info(D, s, q, m, mb)                 # fills D.mixed
    shown = 0
    for (s, m, b), v in sorted(D.mixed.items()):
        shown += 1
        if v["aligned"]:
            R.add("I14", f"{s} {m}", INFO, "", "", "", "same retry mode", f"density-fitted retry for {', '.join(v['species'])}: partner(s) recomputed with density fitting, barrier uses them")
        else:
            R.add("I14", f"{s} {m}", WARN, "", "", "", "same retry mode",
                  f"barrier would mix density-fitted ({', '.join(v['species'])}) and unfitted ({', '.join(v['partners_missing'])}) energies: dropped from all tables until the partner is recomputed (rerun G1 of the GPU notebook)")
    if not shown:
        R.add("I14", "all", PASS, n_df, "", "", "no mixed retry modes", f"{n_df} density-fitted retry point(s) among {n_rows} GPU rows; no barrier mixes fitted and unfitted energies")
    n_path_df = sum(1 for r in D.rows if r["x"].get("df") and is_path_id(r["pid"]))
    if n_path_df:
        R.add("I14", "paths", INFO, n_path_df, "", "", "information", f"{n_path_df} path point(s) were computed with the density-fitting retry (error of order 0.1 kcal/mol in a path profile)")


def i15_basis_correction(D, R):
    """Size of the MP2 basis correction of the CCSD(T) barrier (review row 5)."""
    any_ = False
    for s in REACTIONS:
        if D.tiers.get(s, SYSTEMS[s]["tier"]) not in ("A", "B"):
            continue
        q = "barrier_fwd_cplx" if SYSTEMS[s]["sn2"] else "barrier_fwd"; mb = basis_for(s, "main")
        ref = reference_quantity(D, s, q, mb)
        if ref is None:
            continue
        if not ref["corrected"]:
            R.add("I15", s, NOTRUN, msg="no MP2 basis correction rows (C3b) for the barrier species: raw CCSD(T) used"); continue
        any_ = True
        d = ref["value"] - ref["raw"]
        R.add("I15", s, WARN if abs(d) > 3.0 else INFO, _fnum(d), "", "", "information (WARN above 3)",
              f"MP2/{basis_for(s, 'large')} minus MP2/{mb} shifts the CCSD(T) {q} by {d:+.2f} (raw {ref['raw']:.2f} -> corrected {ref['value']:.2f}), level {ref['level']}")
    if not any_ and not [x for x in R.items if x["item"] == "I15"]:
        R.add("I15", "all", NOTRUN, msg="C3b (MP2 basis correction) has not run yet")


# ------------------------------------------------------------------ decision (rule 4.4) and paper-level rule
def decide(D, R):
    out = {}
    for s in SYSTEM_ORDER:
        S = SYSTEMS[s]
        if s == "ch4":
            out[s] = ("negative result (F6)", "bond stretch, not a reaction: T1 and stability decide whether the reference is used"); continue
        sts = lambda item: [x["status"] for x in R.items if x["item"] == item and sys_of(x["scope"]) == s]
        i1, i4, i6 = worst(sts("I1")), worst(sts("I4")), worst(sts("I6"))
        r2 = sts("R2"); r1 = [x for x in sts("R1") if x != NA]
        i2i3 = sts("I2") + sts("I3")
        if FAIL in (i1, i4, i6):
            why = [n for n, v in (("I1", i1), ("I4", i4), ("I6", i6)) if v == FAIL]
            out[s] = ("SUPPLEMENT", f"{'/'.join(why)} failed: retry once from a fresh initial guess; if it fails again, supplement with the failure stated in one sentence in the main text"); continue
        if FAIL in r2:
            out[s] = ("SUPPLEMENT", "R2 failed (wB97M-V far from literature): geometry or setup error suspected; failure stated in the main text"); continue
        missing = [n for n, v in (("I1", i1), ("I4", i4), ("I6", i6)) if v == NOTRUN]
        if r2 == [] or all(x == NOTRUN for x in r2):
            missing.append("R2")
        if D.tiers.get(s, S["tier"]) in ("A", "B") and (not r1 or all(x == NOTRUN for x in r1)) and FAIL not in i2i3:
            missing.append("R1")
        if missing:
            out[s] = ("INCOMPLETE", "not run yet: " + ", ".join(missing)); continue
        if FAIL in r1 or FAIL in i2i3:
            out[s] = ("MAIN, CCSD(T) role withdrawn", "R1 failed or I2/I3 exclude the transition state: kept on the literature reference"); continue
        out[s] = ("MAIN", "I1 (core methods), I4, I6 (CCSD(T) and anchor) not failed; R2 PASS or WARN" + (" (R1 too)" if r1 else ""))
    return out


def paper_level(dec):
    main = [s for s in REACTIONS if dec[s][0].startswith("MAIN")]
    inc = [s for s in REACTIONS if dec[s][0] == "INCOMPLETE"]
    classes = {SYSTEMS[s]["cls"] for s in main}
    n = len(main)
    if n >= 6 and len(classes) >= 3:
        v = "FULL PAPER as defined"
    elif n >= 4:
        v = "add verified extras, or submit as a short communication titled as a case study"
    else:
        v = "DO NOT SUBMIT: fix the protocol"
    if inc:
        v = f"PROVISIONAL ({len(inc)} system(s) still incomplete: {', '.join(inc)}); so far: " + v
    return f"{n} of 8 reactions in the main text, {len(classes)} classes: {v}"


# ------------------------------------------------------------------ main entry
def run_check(root=None, write=True):
    """Run everything.  Returns dict(items, counts, decisions, paper)."""
    if root:
        setup_root(root, mount=False)
    print("=" * 100); print("Paper B1 v5 automated check (C7)"); print("=" * 100)
    D = Data()
    print(f"CSV files found: {', '.join(D.files) or 'none'}  ({len(D.rows)} rows)\n")
    R = Report()
    for name, fn in (("literature checks R1-R4", literature_checks), ("I1", i1_scf), ("I2/I3", i2_i3_cc), ("I4", i4_imag), ("I5/I6", i5_i6),
                     ("I7", i7_rank), ("I8", i8_basis), ("I9", i9_v4), ("I10", i10_fc), ("I11", i11_geom), ("I12", i12_scout), ("I13", i13_spread),
                     ("I14", i14_retry_mode), ("I15", i15_basis_correction)):
        try:
            fn(D, R)
        except Exception as ex:                       # a bug in one check must never stop the others
            R.add(name.split()[0], "check code", WARN, msg=f"check crashed: {repr(ex)[:150]}")
    dec = decide(D, R)
    paper = paper_level(dec)
    c = R.counts()
    print("\n" + "-" * 100); print("Decision per system (rule 4.4)")
    for s in SYSTEM_ORDER:
        print(f"  {SYSTEMS[s]['n']} {SYSTEMS[s]['name']:18s} -> {dec[s][0]:30s} {dec[s][1]}")
    print("\nPaper-level rule:", paper)
    print("\nCounts: " + "  ".join(f"{k} {v}" for k, v in c.items()))
    if write:
        with open(path("check_report.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["item", "scope", "status", "value", "reference", "difference", "tolerance", "message", "raw_value", "raw_difference"]); w.writeheader(); w.writerows(R.items)
        md = ["# Check report, Paper B1 v5", "", f"Written {now()}.  Tolerances in kcal/mol unless stated.", "",
              "Counts: " + ", ".join(f"{k} {v}" for k, v in c.items()), "", "## Decision per system", "", "| # | System | Tier | Decision | Why |", "|---|---|---|---|---|"]
        for s in SYSTEM_ORDER:
            md.append(f"| {SYSTEMS[s]['n']} | {SYSTEMS[s]['name']} | {D.tiers.get(s, SYSTEMS[s]['tier'])} | {dec[s][0]} | {dec[s][1]} |")
        md += ["", f"**Paper-level rule:** {paper}", "", "## All items", "", "| status | item | scope | value | reference | difference | tolerance | message | raw CCSD(T) value | raw difference |", "|---|---|---|---|---|---|---|---|---|---|"]
        for x in R.items:
            md.append(f"| {x['status']} | {x['item']} | {x['scope']} | {x['value']} | {x['reference']} | {x['difference']} | {x['tolerance']} | {str(x['message']).replace('|', '/')} | {x.get('raw_value', '')} | {x.get('raw_difference', '')} |")
        open(path("check_report.md"), "w").write("\n".join(md) + "\n")
        print("\nWrote check_report.csv and check_report.md in", path(""))
    return dict(items=R.items, counts=c, decisions=dec, paper=paper)
