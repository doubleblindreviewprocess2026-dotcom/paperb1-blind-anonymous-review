# ==============================================================================
# b1v5_cpu.py  --  calculation helpers for the CPU notebook (PySCF, geomeTRIC, xTB).
# Pasted into the second code cell of the CPU notebook; also importable.
# Needs b1v5_common (already loaded in the notebook).  Plain functions, no hidden state.
#
# Levels of theory used:
#   geometry / hessian : density-fitted B3LYP-D3(BJ), def2-TZVP (def2-TZVPD for SN2 species)
#   reference energy   : CCSD(T), frozen core, RHF reference, def2-TZVP(D)
#   scout geometry     : GFN2-xTB via tblite
# ==============================================================================
import os, time, tempfile, logging, math, json
import numpy as np
try:
    b1v5_common_loaded = COLUMNS  # noqa  (names exist when the notebook cell above ran)
except NameError:
    from b1v5_common import *     # noqa

Z_OF = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "Cl": 17}
MEM_MB = 9000   # PySCF memory limit; the plan uses 9 GB on a 12 GB Colab machine


def quiet_geometric():
    """geomeTRIC is very chatty; keep only warnings."""
    for name in list(logging.root.manager.loggerDict) + ["geometric", "geometric.nifty", ""]:
        if name == "" or name.startswith("geometric"):
            lg = logging.getLogger(name)
            if name:
                lg.setLevel(logging.ERROR)
    try:
        import geometric.nifty as nf
        nf.logger.setLevel(logging.ERROR)
        for h in list(nf.logger.handlers):
            h.setLevel(logging.ERROR)
    except Exception:
        pass


def atoms_str(atoms):
    return "; ".join(f"{s} {x:.8f} {y:.8f} {z:.8f}" for s, (x, y, z) in atoms)


def make_mol(atoms, charge, basis):
    from pyscf import gto
    return gto.M(atom=atoms_str(atoms), basis=basis, charge=int(charge), spin=0, unit="Angstrom", verbose=0, max_memory=MEM_MB)


def mol_atoms(mol):
    xyz = mol.atom_coords() * BOHR
    return [(mol.atom_pure_symbol(i), tuple(map(float, xyz[i]))) for i in range(mol.natm)]


# ------------------------------------------------------------------ DFT (geometries, Hessians, zero-point energies)
def make_ks(atoms, charge, basis, xc="b3lyp", disp="d3bj", df=True, grid_level=3):
    from pyscf import dft
    mol = make_mol(atoms, charge, basis)
    mf = dft.RKS(mol, xc=xc)
    if df:
        mf = mf.density_fit()
    if disp:
        mf.disp = disp
    mf.grids.level = grid_level
    mf.conv_tol = 1e-9
    mf.max_cycle = 150
    return mf


def dft_energy(atoms, charge, basis, **kw):
    mf = make_ks(atoms, charge, basis, **kw)
    e = mf.kernel()
    if not mf.converged:                       # one more try with Newton steps
        mf = mf.newton(); e = mf.kernel()
    return dict(energy=float(e), converged=bool(mf.converged))


def frequencies(mf):
    """Harmonic analysis from the analytic Hessian.
    Returns (frequencies in cm-1 with imaginary ones negative, zero-point energy in hartree, displacement vectors of the imaginary modes)."""
    from pyscf.hessian import thermo
    mol = mf.mol
    if mol.natm == 1:
        return [], 0.0, []
    h = mf.Hessian().kernel()
    r = thermo.harmonic_analysis(mol, h)
    # PySCF returns imaginary modes as complex numbers (0 + 1900j); we turn them into negative numbers
    fq = [(-abs(f.imag) if abs(np.imag(f)) > 1e-6 else float(np.real(f))) for f in np.array(r["freq_wavenumber"], dtype=complex)]
    zpe = 0.5 * sum(f for f in fq if f > 0) / HARTREE_CM1
    order = [i for i in np.argsort(fq) if fq[i] < -20.0]
    return fq, zpe, [np.array(r["norm_mode"][i]) for i in order]


def analyse_stationary(atoms, charge, basis, **kw):
    """Energy + harmonic frequencies + zero-point energy at a fixed geometry."""
    mf = make_ks(atoms, charge, basis, **kw)
    e = mf.kernel()
    ok = bool(mf.converged)
    fq, zpe, modes = frequencies(mf)
    n_imag = int(sum(1 for f in fq if f < -20.0))      # -20 cm-1 guards against numerical noise near zero
    return dict(energy=float(e), converged=ok, imag=n_imag, imag_modes=modes,
                extra=dict(zpe_hartree=zpe, lowest_freq_cm1=(min(fq) if fq else None), freqs=[round(float(f), 1) for f in fq][:12]))


def refine_stationary(atoms, charge, basis, kind, **kw):
    """Frequencies at `atoms`; if there are MORE imaginary modes than wanted (minimum 0, transition state 1), push the geometry along
    the extra mode(s) and optimise once more (rule 4.4: one retry from a fresh guess).  Returns (atoms, analysis, retried)."""
    want = 1 if kind == "ts" else 0
    a = analyse_stationary(atoms, charge, basis, **kw)
    if a["imag"] <= want or len(atoms) < 3:
        return atoms, a, False
    extra_modes = a["imag_modes"][want:]
    d = sum(extra_modes)
    d = 0.25 * d / np.abs(d).max()                     # largest atom moves 0.25 A
    pushed = [(s, tuple(np.array(c) + d[i])) for i, (s, c) in enumerate(atoms)]
    r = opt_dft(pushed, charge, basis, transition=(kind == "ts"), hessian=(kind == "ts"), **kw)
    a2 = analyse_stationary(r["atoms"], charge, basis, **kw)
    if abs(a2["imag"] - want) < abs(a["imag"] - want):
        a2["converged"] = bool(a2["converged"] and r["converged"])
        return r["atoms"], a2, True
    return atoms, a, True


def constraint_text(items):
    """items: list of ('distance'|'angle'|'dihedral', atom indices (0-based), value) -> geomeTRIC constraint file text."""
    lines = ["$set"]
    for kind, idx, val in items:
        lines.append(f"{kind} " + " ".join(str(i + 1) for i in idx) + f" {val:.6f}")
    return "\n".join(lines) + "\n"


def opt_dft(atoms, charge, basis, constraints=None, transition=False, hessian=False, maxsteps=120, **kw):
    """Geometry optimisation (minimum or transition state) with geomeTRIC.  `constraints` as for constraint_text.
    hessian=True gives geomeTRIC the analytic Hessian at the first step (needed for TS searches)."""
    from pyscf.geomopt import geometric_solver
    quiet_geometric()
    mf = make_ks(atoms, charge, basis, **kw)
    opts = dict(convergence_set="GAU", maxsteps=maxsteps)
    if transition:
        opts["transition"] = True
    if hessian:
        opts["hessian"] = True
    with tempfile.TemporaryDirectory() as d:
        if constraints:
            cp = os.path.join(d, "cons.txt"); open(cp, "w").write(constraint_text(constraints)); opts["constraints"] = cp
        conv, mol = geometric_solver.kernel(mf, assert_convergence=False, **opts)
    return dict(atoms=mol_atoms(mol), converged=bool(conv))


# ------------------------------------------------------------------ CCSD(T)
def ccsdt(atoms, charge, basis, frozen_core=True, stability=True):
    """CCSD(T) with RHF reference.  Returns energy, convergence, T1 diagnostic and RHF stability in 'extra'."""
    from pyscf import scf, cc
    from pyscf.data.elements import chemcore
    mol = make_mol(atoms, charge, basis)
    mf = scf.RHF(mol); mf.conv_tol = 1e-10; mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:
        mf = mf.newton(); mf.kernel()
    ncore = chemcore(mol) if frozen_core else 0
    extra = dict(frozen_core=bool(frozen_core), n_frozen=int(ncore), e_scf=float(mf.e_tot), nbf=int(mol.nao))
    if stability and mol.nelectron > 2:
        try:
            _, _, st_i, st_e = mf.stability(internal=True, external=True, return_status=True)
            extra["stable"] = int(bool(st_i) and bool(st_e)); extra["stable_internal"] = int(bool(st_i)); extra["stable_external"] = int(bool(st_e))
        except Exception as ex:
            extra["stable"] = ""; extra["stability_error"] = repr(ex)[:100]
    mycc = cc.CCSD(mf); mycc.frozen = ncore; mycc.conv_tol = 1e-8; mycc.max_cycle = 100
    mycc.kernel()
    et = mycc.ccsd_t()
    try:
        t1 = float(mycc.get_t1_diagnostic())
    except Exception:                                  # older PySCF without the helper
        t1 = float(np.sqrt(np.linalg.norm(mycc.t1) ** 2 / (2 * mycc.t1.shape[0])))
    extra.update(e_ccsd_corr=float(mycc.e_corr), e_t=float(et))
    return dict(energy=float(mf.e_tot + mycc.e_corr + et), converged=bool(mf.converged and mycc.converged), t1=round(t1, 5), extra=extra)


# ------------------------------------------------------------------ DF-MP2 for the basis correction of the CCSD(T) reference (review row 5)
def dfmp2(atoms, charge, basis, frozen_core=True):
    """Total DF-MP2 energy (RHF with density fitting + DF-MP2 correlation, frozen core), used ONLY for the difference
    MP2/def2-QZVPP(D) - MP2/def2-TZVP(D) that is added to the CCSD(T) energy of the same geometry.
    The fitting error of the two bases largely cancels in that difference.  RHF: auxiliary basis def2-universal-jkfit; MP2: PySCF default RI basis."""
    from pyscf import scf, mp
    from pyscf.data.elements import chemcore
    mol = make_mol(atoms, charge, basis)
    mf = scf.RHF(mol).density_fit(auxbasis="def2-universal-jkfit"); mf.conv_tol = 1e-9; mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:
        mf = scf.RHF(mol).density_fit(auxbasis="def2-universal-jkfit").newton(); mf.conv_tol = 1e-9; mf.kernel()
    ncore = chemcore(mol) if frozen_core else 0
    pt = mp.dfmp2.DFMP2(mf); pt.frozen = ncore
    pt.kernel()
    return dict(energy=float(mf.e_tot + pt.e_corr), converged=bool(mf.converged),
                extra=dict(frozen_core=bool(frozen_core), n_frozen=int(ncore), e_scf=float(mf.e_tot), e_mp2_corr=float(pt.e_corr), nbf=int(mol.nao),
                           auxbasis_jk="def2-universal-jkfit", auxbasis_mp2="pyscf default RI"))


# ------------------------------------------------------------------ B2PLYP (optional phase): hybrid + scaled DF-MP2
def b2plyp(atoms, charge, basis):
    """B2PLYP = 0.53 HF + 0.47 B88 exchange + 0.73 LYP + 0.27 * MP2 correlation (Grimme 2006); frozen core, density fitting.
    NOT yet checked against a published energy (see README)."""
    from pyscf import dft, mp
    from pyscf.data.elements import chemcore
    mol = make_mol(atoms, charge, basis)
    mf = dft.RKS(mol); mf.xc = "0.53*HF + 0.47*B88, 0.73*LYP"; mf = mf.density_fit()
    mf.grids.level = 3; mf.conv_tol = 1e-9; mf.kernel()
    pt = mp.dfmp2.DFMP2(mf); pt.frozen = chemcore(mol)
    pt.kernel()
    return dict(energy=float(mf.e_tot + 0.27 * pt.e_corr), converged=bool(mf.converged),
                extra=dict(e_hybrid=float(mf.e_tot), e_mp2_corr=float(pt.e_corr)))


# ------------------------------------------------------------------ xTB (tblite) with a geomeTRIC engine
def xtb_eg(atoms, charge, method="GFN2-xTB"):
    """Energy (hartree) and gradient (hartree/bohr, flat) for GFN2-xTB."""
    from tblite.interface import Calculator
    nums = np.array([Z_OF[s] for s, _ in atoms])
    pos = np.array([c for _, c in atoms]) / BOHR
    calc = Calculator(method, nums, pos, charge=float(charge))
    calc.set("verbosity", 0)
    res = calc.singlepoint()
    return float(res.get("energy")), np.array(res.get("gradient")).ravel()


def xtb_energy(atoms, charge):
    try:
        e, _ = xtb_eg(atoms, charge)
        return dict(energy=e, converged=True)
    except Exception:
        return dict(energy=float("nan"), converged=False)


def _xtb_engine(atoms, charge):
    import geometric.engine as ge, geometric.molecule as gm

    class XtbEngine(ge.Engine):
        def __init__(self):
            M = gm.Molecule(); M.elem = [s for s, _ in atoms]; M.xyzs = [np.array([c for _, c in atoms])]
            super().__init__(M)
            self.n = 0

        def calc_new(self, coords, dirname):
            self.n += 1
            at = [(s, tuple(c * BOHR)) for (s, _), c in zip(atoms, coords.reshape(-1, 3))]
            e, g = xtb_eg(at, charge)
            return {"energy": e, "gradient": g}
    return XtbEngine()


def opt_xtb(atoms, charge, constraints=None, transition=False, hessian=None, maxiter=300):
    """Geometry optimisation with GFN2-xTB.  A single atom is returned unchanged."""
    if len(atoms) == 1:
        return dict(atoms=atoms, converged=True, energy=xtb_eg(atoms, charge)[0])
    import geometric.optimize as go
    quiet_geometric()
    eng = _xtb_engine(atoms, charge)
    with tempfile.TemporaryDirectory() as d:
        kw = dict(customengine=eng, input=os.path.join(d, "x"), maxiter=maxiter, convergence_set="GAU", transition=transition)
        if constraints:
            cp = os.path.join(d, "cons.txt"); open(cp, "w").write(constraint_text(constraints)); kw["constraints"] = cp
        if hessian is not None:
            hp = os.path.join(d, "hess.txt"); np.savetxt(hp, hessian); kw["hessian"] = "file:" + hp
        try:
            M = go.run_optimizer(**kw); conv = True
        except Exception:                      # not converged, or xTB SCF failure: report 'not converged', never crash
            conv = False; M = None
    if M is None:
        return dict(atoms=atoms, converged=False, energy=float("nan"))
    new = [(s, tuple(map(float, c))) for s, c in zip(M.elem, M.xyzs[-1])]
    return dict(atoms=new, converged=conv, energy=xtb_eg(new, charge)[0])


# ------------------------------------------------------------------ relaxed scans
def coord_value(kind, idx, atoms):
    """Value of a distance (A), angle (deg) or dihedral (deg) in `atoms`."""
    p = [np.array(atoms[i][1]) for i in idx]
    if kind == "distance":
        return float(np.linalg.norm(p[0] - p[1]))
    if kind == "angle":
        a, b = p[0] - p[1], p[2] - p[1]
        return float(np.degrees(np.arccos(np.clip(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b), -1, 1))))
    b0, b1, b2 = p[0] - p[1], p[2] - p[1], p[3] - p[2]
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1; w = b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


# One driving coordinate per reaction (atom numbers 0-based, in the order of the species' start geometry in b1v5_common).
# 'start' = species whose geometry is point zero of the walk; values(g0) = the n driven values in path order (reactant side first).
# The walk starts at the value closest to the start geometry and goes out in both directions, each step starting from the previous one.
SCANS = {
    "hcn": dict(kind="angle", idx=(2, 0, 1), start="HCN", n=21, label="H-C-N angle (deg)",
                values=lambda g0: np.array([180.0] + list(np.linspace(171.0, 6.0, 20)))),
    "nh3": dict(kind="dihedral", idx=(0, 1, 2, 3), start="NH3", n=21, label="N-H1-H2-H3 dihedral (deg)",
                values=lambda g0: np.linspace(coord_value("dihedral", (0, 1, 2, 3), g0), -coord_value("dihedral", (0, 1, 2, 3), g0), 21)),
    "h2co": dict(kind="distance", idx=(2, 3), start="TS", n=21, label="H...H distance (A)",
                 values=lambda g0: np.linspace(1.85, 0.80, 21)),
    "ch3nc": dict(kind="angle", idx=(0, 1, 2), start="CH3NC", n=11, label="C-N-C angle (deg)",
                  values=lambda g0: np.array([180.0] + list(np.linspace(171.0, 9.0, 10)))),
    "sn2_ff": dict(kind="distance", idx=(0, 2), start="RC", n=11, label="C...F(attacking) distance (A)",
                   values=lambda g0: np.linspace(coord_value("distance", (0, 2), g0), 1.43, 11)),
    # the three systems below are scanned with xTB only (C4); DFT paths for them are an optional GPU extra
    "sn2_clcl": dict(kind="distance", idx=(0, 2), start="RC", n=11, label="C...Cl(attacking) distance (A)", dft=False,
                     values=lambda g0: np.linspace(coord_value("distance", (0, 2), g0), 1.80, 11)),
    "sn2_fcl": dict(kind="distance", idx=(0, 2), start="RC", n=11, label="C...F(attacking) distance (A)", dft=False,
                    values=lambda g0: np.linspace(coord_value("distance", (0, 2), g0), 1.40, 11)),
    "sn2_ohf": dict(kind="distance", idx=(0, 5), start="RC", n=11, label="C...O(attacking) distance (A)", dft=False,
                    values=lambda g0: np.linspace(coord_value("distance", (0, 5), g0), 1.42, 11)),
}
SCAN_STRIDE = 1      # the sandbox test sets 3: every third point only
# CH4 has no scan: the 9 stored v4 stretch points are reused (b1v5_common.CH4_STRETCH).


def _walk(system, g0, charge, relax, energy, store, save):
    """Generic relaxed walk.  relax(atoms, constraint) -> dict(atoms, converged); energy(atoms) -> dict(energy, converged)."""
    sp = SCANS[system]
    vals = np.asarray(sp["values"](g0), float)
    keep = sorted(set(range(0, len(vals), SCAN_STRIDE)) | {len(vals) - 1})
    vals = vals[keep]
    n = len(vals)
    c0 = coord_value(sp["kind"], sp["idx"], g0)
    i0 = int(np.argmin(np.abs(vals - c0)))
    pts = {int(p["index"]): p for p in store if int(p["index"]) < n}
    def one(k, prev_atoms, free=False):
        t0 = time.time()
        if free:
            r = dict(atoms=prev_atoms, converged=True)
        else:
            r = relax(prev_atoms, (sp["kind"], sp["idx"], float(vals[k])))
        e = energy(r["atoms"])
        pts[k] = dict(index=k, coordinate=float(vals[k]) if not free else c0, atoms=r["atoms"], energy=e["energy"],
                      converged=bool(r["converged"] and e["converged"]), seconds=time.time() - t0)
        save([pts[j] for j in sorted(pts)])
        print(f"   scan {system} point {k + 1}/{n}  coordinate {pts[k]['coordinate']:.3f}  E={e['energy']:.6f}  converged={pts[k]['converged']}  {time.time() - t0:.0f} s")
    if i0 not in pts:
        one(i0, g0, free=True)
    for k in range(i0 + 1, n):                # towards the end of the list
        if k not in pts:
            check_time(); one(k, pts[k - 1]["atoms"])
    for k in range(i0 - 1, -1, -1):           # towards the beginning
        if k not in pts:
            check_time(); one(k, pts[k + 1]["atoms"])
    return [pts[j] for j in sorted(pts)]


def scan_points_dft(system, g0, charge, basis, store, save, maxsteps=100):
    """Relaxed scan with DF-B3LYP-D3(BJ); `store` = finished points, `save(points)` is called after every point."""
    return _walk(system, g0, charge,
                 lambda at, c: opt_dft(at, charge, basis, constraints=[c], maxsteps=maxsteps),
                 lambda at: dft_energy(at, charge, basis), store, save)


def scan_points_xtb(system, g0, charge):
    """The same walk with GFN2-xTB (fast; no resume needed)."""
    return _walk(system, g0, charge, lambda at, c: opt_xtb(at, charge, constraints=[c]),
                 lambda at: xtb_energy(at, charge), [], lambda pts: None)


def peak_of(points):
    """Index of the highest-energy scan point, or None when the highest point is an end point (no peak found)."""
    E = np.array([p["energy"] for p in points])
    k = int(np.nanargmax(E))
    return None if k in (0, len(E) - 1) else k


def parabola_peak(points, k):
    """Three-point parabola through the highest point and its two neighbours; returns (coordinate, energy) of the vertex."""
    x = np.array([points[j]["coordinate"] for j in (k - 1, k, k + 1)]); y = np.array([points[j]["energy"] for j in (k - 1, k, k + 1)])
    a, b, c = np.polyfit(x, y, 2)
    if a >= 0:
        return float(x[1]), float(y[1])
    xv = -b / (2 * a)
    if not (x.min() <= xv <= x.max()):
        return float(x[1]), float(y[1])
    return float(xv), float(np.polyval([a, b, c], xv))
