# ==============================================================================
# b1v5_geompack.py  --  the GEOMETRY PACK of the 1-hour design (PLAN_1hour.md, decision D1 and section 7.1).
#
# What it is, in plain words:
#   The geometry work (optimisations, transition states with ONE imaginary mode, Hessians, zero-point energies, relaxed scans, GFN2-xTB
#   minima and peaks) is deterministic CPU work.  It was done ONCE, in the sandbox, with the same helper functions the full notebooks use
#   (b1v5_cpu.py: opt_dft, refine_stationary, scan_points_dft, opt_xtb, ...), and stored in ONE JSON file with a SHA-256.  Both 1-hour
#   notebooks carry that file inside them (gzip + base64 in a cell), check its SHA-256, unpack it into geoms/ and phase01 / phase02 /
#   phase04 rows, and the CPU notebook re-verifies the geometries with a gradient test (phase C0).
#
# This file holds:   (1) the builder (also used by REBUILD_GEOMS = True on Colab and by build_geompack.py in the sandbox),
#                    (2) pack loading, hashing, verification and installation into the working folder.
# Pasted unchanged into both 1-hour notebooks (after b1v5_common.py and b1v5_cpu.py, which share the namespace); also importable.
# ==============================================================================
import os, sys, io, json, gzip, base64, hashlib, math, time, platform, datetime
import numpy as np
try:
    b1v5_common_loaded = COLUMNS      # noqa  (names exist when the cells above ran)
except NameError:
    from b1v5_common import *         # noqa
try:
    b1v5_cpu_loaded = ccsdt           # noqa
except NameError:
    from b1v5_cpu import *            # noqa

PACK_FORMAT = "b1v5-geompack-1"
CORE_SYSTEMS = ["hcn", "nh3", "h2co", "ch3nc", "sn2_ff", "sn2_clcl", "ch4"]      # six reactions + the CH4 stretch (PLAN_1hour D3)
EXT78_SYSTEMS = ["sn2_fcl", "sn2_ohf"]                                           # F- + CH3Cl and OH- + CH3F: first extension items
DFT_PATH_SYSTEMS = ["hcn", "nh3", "h2co"]                                        # relaxed DFT scans of the core (CH4: 9 stored rigid points)
PATH_N = {"hcn": 11, "nh3": 9, "h2co": 11}                                       # plan section 2: C3b 'HCN 11 points, NH3 9, H2CO 11, CH4 9'
XTB_N = 21                                                                       # xTB scans are cheap: 21 points for a finer peak
GRAD_RMS_TOL = 3.0e-4                                                            # pass criterion of the gradient test (hartree/bohr)
GRAD_RMS_BUILD_TARGET = 1.5e-4                                                   # the builder re-optimises tighter above this
LEVELS_TEXT = ("geometry: DF-B3LYP-D3(BJ) (PySCF density fitting, grid level 3, conv_tol 1e-9), def2-TZVP (def2-TZVPD for every SN2 species), "
               "geomeTRIC convergence set GAU; transition states by geomeTRIC eigenvector following with the analytic Hessian at step 0; "
               "frequencies and zero-point energy from the analytic DF-B3LYP-D3(BJ) Hessian; xTB: GFN2-xTB (tblite) with a geomeTRIC engine")


# ------------------------------------------------------------------ hashing
def canon_atoms(atoms):
    return [[s, [round(float(v), 8) for v in c]] for s, c in atoms]


def struct_sha(system, species, level, charge, atoms):
    blob = json.dumps([system, species, level, int(charge), canon_atoms(atoms)], separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()


def path_sha(system, pts):
    blob = json.dumps([system, [[int(p["index"]), round(float(p["coordinate"]), 8), canon_atoms([(s, tuple(c)) for s, c in p["atoms"]])] for p in pts]], separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()


def pack_bytes(pack):
    """Canonical bytes of a pack (sorted keys, no spaces): what the SHA-256 of the pack is computed over."""
    return json.dumps(pack, sort_keys=True, separators=(",", ":")).encode()


def pack_sha256(pack):
    return hashlib.sha256(pack_bytes(pack)).hexdigest()


def pack_encode(pack):
    """gzip + base64 text of a pack (what the notebooks embed)."""
    return base64.b64encode(gzip.compress(pack_bytes(pack), 9)).decode()


def pack_decode(text):
    return json.loads(gzip.decompress(base64.b64decode(text)).decode())


def key_of(system, species, level):
    return f"{system}/{species}/{level}"


# ------------------------------------------------------------------ the relaxed-scan specifications of the 1-hour design
def _scan_specs(ts_angle_hcn=72.0, n_xtb=None):
    """One driving coordinate per reaction, as in b1v5_cpu.SCANS, but with 11 / 9 / 11 points for the DFT paths (HCN, NH3, H2CO) and, if n_xtb is
    given, n_xtb points for the xTB scans (all six reactions).  The HCN grid contains the angle of our own transition state, so the path
    has a point at the saddle (constrained optimisation there = the transition state)."""
    n = n_xtb
    X = {}
    X["hcn"] = dict(kind="angle", idx=(2, 0, 1), start="HCN", label="H-C-N angle (deg)",
                    values=(lambda g0: np.array([180.0] + list(np.linspace(171.0, 6.0, n - 1)))) if n else
                    (lambda g0: np.array(sorted({180.0, round(float(ts_angle_hcn), 3)} | set(np.round(np.linspace(171.0, 6.0, 9), 3)), reverse=True))))
    X["nh3"] = dict(kind="dihedral", idx=(0, 1, 2, 3), start="NH3", label="N-H1-H2-H3 dihedral (deg)",
                    values=lambda g0, k=(n or PATH_N["nh3"]): np.linspace(coord_value("dihedral", (0, 1, 2, 3), g0), -coord_value("dihedral", (0, 1, 2, 3), g0), k))
    X["h2co"] = dict(kind="distance", idx=(2, 3), start="TS", label="H...H distance (A)",
                     values=lambda g0, k=(n or PATH_N["h2co"]): np.linspace(1.85, 0.80, k))
    X["ch3nc"] = dict(kind="angle", idx=(0, 1, 2), start="CH3NC", label="C-N-C angle (deg)",
                      values=lambda g0, k=(n or 11): np.array([180.0] + list(np.linspace(171.0, 9.0, k - 1))))
    X["sn2_ff"] = dict(kind="distance", idx=(0, 2), start="RC", label="C...F(attacking) distance (A)",
                       values=lambda g0, k=(n or 11): np.linspace(coord_value("distance", (0, 2), g0), 1.43, k))
    X["sn2_clcl"] = dict(kind="distance", idx=(0, 2), start="RC", label="C...Cl(attacking) distance (A)",
                         values=lambda g0, k=(n or 11): np.linspace(coord_value("distance", (0, 2), g0), 1.80, k))
    X["sn2_fcl"] = dict(kind="distance", idx=(0, 2), start="RC", label="C...F(attacking) distance (A)",
                        values=lambda g0, k=(n or 11): np.linspace(coord_value("distance", (0, 2), g0), 1.40, k))
    X["sn2_ohf"] = dict(kind="distance", idx=(0, 5), start="RC", label="C...O(attacking) distance (A)",
                        values=lambda g0, k=(n or 11): np.linspace(coord_value("distance", (0, 5), g0), 1.42, k))
    return X


def _walk_spec(sp, g0, relax, energy, store, save, label=""):
    """The relaxed walk of b1v5_cpu._walk with an explicit scan specification (so that DFT and xTB scans can differ)."""
    vals = np.asarray(sp["values"](g0), float)
    n = len(vals)
    c0 = coord_value(sp["kind"], sp["idx"], g0)
    i0 = int(np.argmin(np.abs(vals - c0)))
    pts = {int(p["index"]): p for p in store if int(p["index"]) < n}

    def one(k, prev_atoms, free=False):
        t0 = time.time()
        r = dict(atoms=prev_atoms, converged=True) if free else relax(prev_atoms, (sp["kind"], sp["idx"], float(vals[k])))
        e = energy(r["atoms"])
        pts[k] = dict(index=k, coordinate=float(vals[k]) if not free else c0, atoms=r["atoms"], energy=e["energy"],
                      converged=bool(r["converged"] and e["converged"]), seconds=time.time() - t0)
        save([pts[j] for j in sorted(pts)])
        print(f"   {label} point {k + 1}/{n}  coordinate {pts[k]['coordinate']:.3f}  E={e['energy']:.6f}  converged={pts[k]['converged']}  {time.time() - t0:.0f} s", flush=True)
    if i0 not in pts:
        one(i0, g0, free=True)
    for k in range(i0 + 1, n):
        if k not in pts:
            one(k, pts[k - 1]["atoms"])
    for k in range(i0 - 1, -1, -1):
        if k not in pts:
            one(k, pts[k + 1]["atoms"])
    return [pts[j] for j in sorted(pts)]


# ------------------------------------------------------------------ the builder
class PackBuilder:
    """Builds the pack step by step in a work folder (resumable: every finished step is written to pack_work.json).
        PackBuilder(work).build(systems)   then   .assemble(path_out)
    Every step records its seconds in work['timings']."""

    def __init__(self, work, hardware=None, log=print):
        self.work = os.path.abspath(work)
        os.makedirs(self.work, exist_ok=True)
        self.fn = os.path.join(self.work, "pack_work.json")
        self.W = json.load(open(self.fn)) if os.path.exists(self.fn) else dict(structures={}, paths={}, xtb={}, timings=[], notes=[])
        self.log = log
        self.hardware = hardware

    def save(self):
        json.dump(self.W, open(self.fn + ".tmp", "w"), indent=1, default=float)
        os.replace(self.fn + ".tmp", self.fn)

    def tick(self, step, system, species, seconds, note=""):
        self.W["timings"].append(dict(step=step, system=system, species=species, seconds=round(float(seconds), 1), note=note))
        self.log(f"  [{system}:{species}] {step}: {seconds:.1f} s {note}")
        self.save()

    # ---- one stationary structure: optimise, frequencies, gradient check
    def stationary(self, system, sp, kind):
        key = key_of(system, sp, LV_OWN)
        if key in self.W["structures"]:
            return
        ch, basis = species_charge(system, sp), basis_for(system, "main")
        at = start_atoms(system, sp)
        t0 = time.time()
        conv = True
        if len(at) > 1:
            r = opt_dft(at, ch, basis, transition=(kind == "ts"), hessian=(kind == "ts"))
            at, conv = r["atoms"], r["converged"]
        self.tick("optimisation", system, sp, time.time() - t0, f"kind={kind} converged={conv}")
        t0 = time.time()
        at, a, retried = refine_stationary(at, ch, basis, kind)
        self.tick("frequencies + zero-point energy" + (" (+ retry along the extra mode)" if retried else ""), system, sp, time.time() - t0, f"imag={a['imag']}")
        a = self.regrid_analysis(at, ch, basis, kind, a, system, sp)
        g = self.gradient(at, ch, basis, system, sp)
        if g["rms"] > GRAD_RMS_BUILD_TARGET and len(at) > 1:            # re-optimise tighter once so that the gradient test on Colab has margin
            t0 = time.time()
            at2, a2, g2 = self.tighten(at, ch, basis, kind)
            self.tick("tighter re-optimisation (GAU_TIGHT)", system, sp, time.time() - t0, f"rms {g['rms']:.1e} -> {g2['rms']:.1e}, imag {a['imag']} -> {a2['imag']}")
            if g2["rms"] < g["rms"] and abs(a2["imag"] - (1 if kind == "ts" else 0)) <= abs(a["imag"] - (1 if kind == "ts" else 0)):
                at, a, g = at2, a2, g2
        self.W["structures"][key] = dict(system=system, species=sp, level=LV_OWN, kind=kind, charge=int(ch), basis=basis, atoms=canon_atoms(at),
                                         energy_hartree=a["energy"], scf_converged=bool(a["converged"]), opt_converged=bool(conv), imag=int(a["imag"]),
                                         zpe_hartree=a["extra"].get("zpe_hartree"), lowest_freq_cm1=a["extra"].get("lowest_freq_cm1"), freqs=a["extra"].get("freqs"),
                                         retried=bool(retried), rms_grad=g["rms"], max_grad=g["max"], grad_seconds=g["seconds"], hessian_grid=int(a["extra"].get("hessian_grid", 3)),
                                         sha256=struct_sha(system, sp, LV_OWN, ch, at))
        self.save()

    def regrid_analysis(self, at, ch, basis, kind, a, system, sp):
        """Found after the first sandbox build (Cl- + CH3Cl, def2-TZVPD): at integration grid level 3 the very flat modes of an anion with diffuse functions can come out
        imaginary (TS: -134 cm-1, at level 5 the same geometry gives +206 cm-1 and exactly one imaginary mode).  If the number of imaginary modes is not the wanted one,
        the Hessian is repeated ONCE at grid level 5 at the stored geometry and taken if it is closer to the wanted number.  The geometry, the energy and the gradient test
        stay at level 3 (the level the notebooks use); only frequencies and zero-point energy come from level 5 (recorded as hessian_grid = 5)."""
        want = 1 if kind == "ts" else 0
        if len(at) < 2 or a["imag"] == want:
            return a
        t0 = time.time()
        a5 = analyse_stationary(at, ch, basis, grid_level=5)
        self.tick("Hessian repeated at grid level 5", system, sp, time.time() - t0, f"imag {a['imag']} -> {a5['imag']}")
        if abs(a5["imag"] - want) < abs(a["imag"] - want):
            a5["energy"], a5["converged"] = a["energy"], a["converged"]
            a5["extra"]["hessian_grid"] = 5
            return a5
        return a

    def regrid(self, system, sp):
        """Apply regrid_analysis to a structure that is already in the work folder (used for the structures built before this step existed)."""
        key = key_of(system, sp, LV_OWN)
        v = self.W["structures"][key]
        a = dict(energy=v["energy_hartree"], converged=v["scf_converged"], imag=v["imag"], extra=dict(zpe_hartree=v["zpe_hartree"], lowest_freq_cm1=v["lowest_freq_cm1"], freqs=v["freqs"]))
        a2 = self.regrid_analysis([(s, tuple(c)) for s, c in v["atoms"]], v["charge"], v["basis"], v["kind"], a, system, sp)
        if a2 is not a:
            v.update(imag=int(a2["imag"]), zpe_hartree=a2["extra"]["zpe_hartree"], lowest_freq_cm1=a2["extra"]["lowest_freq_cm1"], freqs=a2["extra"]["freqs"], hessian_grid=5)
            fb = key_of(system, sp, LV_SCANPEAK)
            if v["kind"] == "ts" and v["imag"] == 1 and fb in self.W["structures"]:
                del self.W["structures"][fb]
                self.W["notes"].append(f"{system}: the scan-peak fallback TS was dropped: the own TS has one imaginary mode at grid level 5")
            self.W["notes"].append(f"{system}:{sp}: frequencies and zero-point energy from a grid level 5 Hessian (level 3 gave {a['imag']} imaginary modes)")
            self.save()

    def gradient(self, atoms, ch, basis, system="", sp=""):
        t0 = time.time()
        mf = make_ks(atoms, ch, basis)
        e = mf.kernel()
        g = np.asarray(mf.nuc_grad_method().kernel(), float) if len(atoms) > 0 else np.zeros((1, 3))
        sec = time.time() - t0
        rec = dict(rms=float(np.sqrt(np.mean(g ** 2))), max=float(np.abs(g).max()), seconds=round(sec, 1), energy=float(e))
        if system:
            self.tick("gradient at the stored geometry", system, sp, sec, f"rms={rec['rms']:.2e} max={rec['max']:.2e}")
        return rec

    def tighten(self, atoms, ch, basis, kind):
        from pyscf.geomopt import geometric_solver
        quiet_geometric()
        mf = make_ks(atoms, ch, basis)
        opts = dict(convergence_set="GAU_TIGHT", maxsteps=120)
        if kind == "ts":
            opts["transition"] = True; opts["hessian"] = True
        conv, mol = geometric_solver.kernel(mf, assert_convergence=False, **opts)
        new = mol_atoms(mol)
        a2 = analyse_stationary(new, ch, basis)
        return new, a2, self.gradient(new, ch, basis)

    # ---- relaxed DFT scan of one system (HCN, NH3, H2CO), 11 / 9 / 11 points
    def dft_path(self, system):
        if system in self.W["paths"] and self.W["paths"][system].get("complete"):
            return
        ts = self.W["structures"].get(key_of(system, "TS", LV_OWN))
        ts_angle = 72.0
        if system == "hcn" and ts:
            ts_angle = coord_value("angle", (2, 0, 1), [(s, tuple(c)) for s, c in ts["atoms"]])
        spec = _scan_specs(ts_angle)[system]
        start = self.W["structures"][key_of(system, spec["start"], LV_OWN)]
        g0 = [(s, tuple(c)) for s, c in start["atoms"]]
        ch, basis = species_charge(system, spec["start"]), basis_for(system, "main")
        old = self.W["paths"].get(system, {}).get("points", [])
        store = [dict(p, atoms=[(s, tuple(c)) for s, c in p["atoms"]]) for p in old]

        def save(pts):
            self.W["paths"][system] = dict(system=system, label=spec["label"], basis=basis, level=LV_SCAN, complete=False,
                                           points=[dict(index=int(q["index"]), coordinate=float(q["coordinate"]), atoms=canon_atoms(q["atoms"]), energy=float(q["energy"]),
                                                        converged=bool(q["converged"]), seconds=round(float(q.get("seconds", 0.0)), 1)) for q in pts])
            self.save()
        t0 = time.time()
        pts = _walk_spec(spec, g0, lambda at, c: opt_dft(at, ch, basis, constraints=[c], maxsteps=100), lambda at: dft_energy(at, ch, basis), store, save, label=f"{system} DFT scan")
        P = self.W["paths"][system]
        P["complete"] = True; P["sha256"] = path_sha(system, P["points"])
        k = peak_of(pts)
        P["peak_index"] = k
        self.tick("relaxed DFT scan", system, "path", time.time() - t0, f"{len(pts)} points, peak index {k}")

    # ---- the scan-peak fallback for a TS search that failed
    def scanpeak_fallback(self, system):
        spec = _scan_specs()[system]
        start = self.W["structures"].get(key_of(system, spec["start"], LV_OWN))
        if not start:
            return
        g0 = [(s, tuple(c)) for s, c in start["atoms"]]
        ch, basis = species_charge(system, spec["start"]), basis_for(system, "main")
        t0 = time.time()
        pts = _walk_spec(spec, g0, lambda at, c: opt_dft(at, ch, basis, constraints=[c], maxsteps=100), lambda at: dft_energy(at, ch, basis), [], lambda p: None, label=f"{system} fallback scan")
        k = peak_of(pts)
        if k is None:
            self.W["notes"].append(f"{system}: TS search failed and the fallback scan has no peak")
            return
        a = analyse_stationary(pts[k]["atoms"], species_charge(system, "TS"), basis)
        self.W["structures"][key_of(system, "TS", LV_SCANPEAK)] = dict(system=system, species="TS", level=LV_SCANPEAK, kind="ts", charge=int(species_charge(system, "TS")), basis=basis,
                                                                    atoms=canon_atoms(pts[k]["atoms"]), energy_hartree=a["energy"], scf_converged=bool(a["converged"]), opt_converged=True,
                                                                    imag=int(a["imag"]), zpe_hartree=a["extra"].get("zpe_hartree"), lowest_freq_cm1=a["extra"].get("lowest_freq_cm1"),
                                                                    freqs=a["extra"].get("freqs"), retried=False, rms_grad=None, max_grad=None, grad_seconds=None,
                                                                    sha256=struct_sha(system, "TS", LV_SCANPEAK, species_charge(system, "TS"), pts[k]["atoms"]), from_scan_point=int(k))
        self.tick("scan-peak fallback for the TS", system, "TS", time.time() - t0, f"point {k}")

    # ---- the CH4 stretch: the nine stored rigid points (geometry only, energies from C3b / G1b)
    def ch4_path(self):
        if "ch4" in self.W["paths"]:
            return
        pts = [dict(index=i, coordinate=float(s), atoms=canon_atoms(ch4_geometry(float(s))), energy=None, converged=True, seconds=0.0) for i, s in enumerate(CH4_STRETCH)]
        self.W["paths"]["ch4"] = dict(system="ch4", label="C-H stretch (A)", basis=basis_for("ch4", "main"), level=LV_SCAN, complete=True, points=pts, sha256=path_sha("ch4", pts))
        self.save()

    # ---- GFN2-xTB scout geometries: minimum (reactant or complex) and the peak of a relaxed xTB scan
    def xtb(self, system):
        if system in self.W["xtb"]:
            return
        S = SYSTEMS[system]
        rmin = S["rc"] or S["reactants"][0]
        ch = species_charge(system, rmin)
        own = self.W["structures"].get(key_of(system, rmin, LV_OWN))
        a0 = [(s, tuple(c)) for s, c in own["atoms"]] if own else start_atoms(system, rmin)
        t0 = time.time()
        r = opt_xtb(a0, ch)
        e_min = xtb_energy(r["atoms"], ch)["energy"]
        self.tick("xTB minimum", system, rmin, time.time() - t0, f"converged={r['converged']}")
        spec = _scan_specs(n_xtb=XTB_N)[system]
        if spec["start"] == rmin:
            g0 = r["atoms"]
        else:
            st = self.W["structures"].get(key_of(system, spec["start"], LV_OWN))
            g0 = [(s, tuple(c)) for s, c in st["atoms"]] if st else start_atoms(system, spec["start"])
        t0 = time.time()
        pts = _walk_spec(spec, g0, lambda at, c: opt_xtb(at, species_charge(system, spec["start"]), constraints=[c]),
                         lambda at: xtb_energy(at, species_charge(system, spec["start"])), [], lambda p: None, label=f"{system} xTB scan")
        k = peak_of(pts)
        ec = species_charge(system, "TS")
        rec = dict(system=system, reactant_species=rmin, charge=int(ch), method="GFN2-xTB (tblite)", scan=dict(label=spec["label"], n_points=len(pts)),
                   minimum=dict(atoms=canon_atoms(r["atoms"]), energy_hartree=float(e_min), converged=bool(r["converged"]), sha256=struct_sha(system, rmin, LV_XTB, ch, r["atoms"])))
        if k is None:
            mid = pts[len(pts) // 2]
            rec["peak"] = dict(atoms=canon_atoms(mid["atoms"]), energy_hartree=None, converged=False, no_peak=True, scan_max_index=int(np.nanargmax([q["energy"] for q in pts])),
                               sha256=struct_sha(system, "TS", LV_XTB, ec, mid["atoms"]))
        else:
            rec["peak"] = dict(atoms=canon_atoms(pts[k]["atoms"]), energy_hartree=float(pts[k]["energy"]), converged=bool(pts[k]["converged"]), no_peak=False, scan_point=int(k),
                               coordinate=float(pts[k]["coordinate"]), sha256=struct_sha(system, "TS", LV_XTB, ec, pts[k]["atoms"]))
        self.W["xtb"][system] = rec
        self.tick("xTB relaxed scan + peak", system, "TS", time.time() - t0, f"{len(pts)} points, peak {k}")

    # ---- everything for a list of systems
    def build(self, systems):
        for system in systems:
            for sp, kind in SPECIES_KIND[system].items():
                with item_guard("geompack", f"{system}:{sp}"):
                    self.stationary(system, sp, kind)
            ts = self.W["structures"].get(key_of(system, "TS", LV_OWN))
            if SYSTEMS[system]["ts"] and (ts is None or ts["imag"] != 1):
                with item_guard("geompack", f"{system}:fallback"):
                    self.scanpeak_fallback(system)
            if system in DFT_PATH_SYSTEMS:
                with item_guard("geompack", f"{system}:path"):
                    self.dft_path(system)
            if system == "ch4":
                self.ch4_path()
        for system in systems:
            if system in REACTIONS:
                with item_guard("geompack", f"{system}:xtb"):
                    self.xtb(system)

    # ---- the pack file
    def assemble(self, systems, design, extra_note=""):
        W = self.W
        pack = dict(format=PACK_FORMAT, design=design, systems=list(systems), created=now(), levels=LEVELS_TEXT, note=extra_note,
                    builder=dict(script="build_geompack.py + b1v5_geompack.py + b1v5_cpu.py", versions=ver(), hardware=self.hardware or hw(), python=platform.python_version()),
                    structures={k: v for k, v in sorted(W["structures"].items()) if v["system"] in systems},
                    paths={k: v for k, v in sorted(W["paths"].items()) if k in systems},
                    xtb={k: v for k, v in sorted(W["xtb"].items()) if k in systems},
                    timings=[t for t in W["timings"] if t["system"] in systems], notes=W["notes"])
        return pack


# ------------------------------------------------------------------ loading and verification
def load_pack_text(text, expected_sha256=None):
    """Decode an embedded pack; returns (pack, sha256 of its canonical bytes, matches).  Never raises on a hash mismatch (the caller decides)."""
    pack = pack_decode(text)
    sha = pack_sha256(pack)
    return pack, sha, (expected_sha256 is None or sha == expected_sha256)


def verify_hashes(pack):
    """Re-compute the SHA-256 of every structure and path.  Returns (n_ok, list of problems)."""
    ok, bad = 0, []
    for k, v in pack["structures"].items():
        h = struct_sha(v["system"], v["species"], v["level"], v["charge"], [(s, tuple(c)) for s, c in v["atoms"]])
        if h == v["sha256"]:
            ok += 1
        else:
            bad.append("structure " + k)
    for k, v in pack["paths"].items():
        if path_sha(k, v["points"]) == v.get("sha256"):
            ok += 1
        else:
            bad.append("path " + k)
    for k, v in pack["xtb"].items():
        for part, sp in (("minimum", v["reactant_species"]), ("peak", "TS")):
            p = v[part]
            if struct_sha(k, sp, LV_XTB, p.get("charge", v["charge"] if part == "minimum" else species_charge(k, "TS")), [(s, tuple(c)) for s, c in p["atoms"]]) == p["sha256"]:
                ok += 1
            else:
                bad.append(f"xtb {k} {part}")
    return ok, bad


def pack_geoms(pack, systems=None):
    """Yield (system, species, level, kind, charge, atoms) for every structure that is a geometry (own, scan-peak TS, xTB)."""
    for k, v in pack["structures"].items():
        if systems is None or v["system"] in systems:
            yield v["system"], v["species"], v["level"], v["kind"], v["charge"], [(s, tuple(c)) for s, c in v["atoms"]]
    for sysname, v in pack["xtb"].items():
        if systems is not None and sysname not in systems:
            continue
        yield sysname, v["reactant_species"], LV_XTB, "min", v["charge"], [(s, tuple(c)) for s, c in v["minimum"]["atoms"]]
        p = v["peak"]
        yield sysname, "TS", LV_XTB, "ts", species_charge(sysname, "TS"), [(s, tuple(c)) for s, c in p["atoms"]]


def install_pack(pack, systems=None, sha=None, hardware_label="geometry pack (sandbox build)"):
    """Write the pack into the working folder: geoms/*.json + *.xyz, geoms/<system>_path.json, and the rows of phase01_stationary.csv (DF-B3LYP-D3(BJ)
    energy, imaginary modes, zero-point energy), phase02_scans.csv (path energies) and phase04_scout_ref.csv (xTB energies).
    Idempotent: a geometry file with the same atoms is left alone; rows that exist are not duplicated.  Returns a dict of counts."""
    n = dict(geoms=0, geoms_kept=0, paths=0, rows=0)
    for system, sp, lv, kind, ch, at in pack_geoms(pack, systems):
        name = geom_name(system, sp, lv)
        g = load_geom(name)
        if g is not None and canon_atoms(g["atoms"]) == canon_atoms(at):
            n["geoms_kept"] += 1
            continue
        save_geom(name, at, ch, species=sp, kind=kind, level=lv, opt_converged=True, basis=basis_for(system, "main"), from_pack=True, pack_sha256=sha or "")
        n["geoms"] += 1
    for system, v in pack["paths"].items():
        if systems is not None and system not in systems:
            continue
        p = path("geoms", f"{system}_path.json")
        d = dict(system=system, level=LV_SCAN, basis=v["basis"], label=v["label"], saved=now(), from_pack=True, pack_sha256=sha or "",
                 points=[dict(index=int(q["index"]), coordinate=float(q["coordinate"]), atoms=[[s, list(map(float, c))] for s, c in q["atoms"]],
                              energy=(None if q.get("energy") is None else float(q["energy"])), converged=bool(q.get("converged", True))) for q in v["points"]])
        json.dump(d, open(p + ".tmp", "w"), indent=1); os.replace(p + ".tmp", p)
        n["paths"] += 1
    def put(csv_name, system, pid, level, method, basis, e, conv, imag="", coordinate="", extra=None):
        if have(csv_name, system, pid, level, method, basis):
            return
        r = make_row(system, pid, level, method, basis, e, conv, coordinate, "", imag, 0.0, dict(extra or {}, from_pack=True))
        r["hardware"] = hardware_label
        append_row(csv_name, r)
        n["rows"] += 1
    for k, v in pack["structures"].items():
        if systems is not None and v["system"] not in systems:
            continue
        put(PHASE_CSV["phase01"], v["system"], v["species"], v["level"], M_DFB3LYP, v["basis"], v["energy_hartree"], v["scf_converged"] and v["opt_converged"], imag=v["imag"],
            extra=dict(zpe_hartree=v["zpe_hartree"], lowest_freq_cm1=v["lowest_freq_cm1"], freqs=v["freqs"], retried=v["retried"], rms_grad_build=v["rms_grad"]))
    for system, v in pack["paths"].items():
        if systems is not None and system not in systems or system == "ch4":
            continue
        for q in v["points"]:
            put(PHASE_CSV["phase02"], system, f"p{int(q['index']):02d}", LV_SCAN, M_DFB3LYP, v["basis"], q["energy"], q["converged"], coordinate=f"{q['coordinate']:.4f}")
    for system, v in pack["xtb"].items():
        if systems is not None and system not in systems:
            continue
        for part, sp in (("minimum", v["reactant_species"]), ("peak", "TS")):
            p = v[part]
            if p.get("no_peak"):
                put(PHASE_CSV["phase04"], system, "TS", LV_XTB, M_XTB, "GFN2", float("nan"), False, extra={"no_peak": True, "scan_max_index": p.get("scan_max_index")})
            else:
                put(PHASE_CSV["phase04"], system, sp, LV_XTB, M_XTB, "GFN2", p["energy_hartree"], p["converged"])
    return n


def gradient_test(pack, systems=None, tol=GRAD_RMS_TOL, which=None, quiet=False, csv_name="phase00_verify.csv", admit=None, run=None):
    """The gradient-norm test of phase C0: DF-B3LYP-D3(BJ) energy + gradient at every stored structure in `which` (a list of (system, species) pairs;
    default every unimolecular structure and the SN2 transition states).  PASS: RMS gradient below `tol` hartree/bohr.  One row per structure in
    `csv_name` (energy = the energy recomputed here; extra holds RMS, max, the stored pack energy and the difference).
    `run(csv, system, pid, level, method, basis, fn, kind)` lets the caller route every calculation through the deadline queue."""
    out = []
    todo = which if which is not None else [(v["system"], v["species"]) for k, v in pack["structures"].items()
                                           if v["level"] in (LV_OWN, LV_SCANPEAK) and (not SYSTEMS[v["system"]]["sn2"] or v["species"] == "TS") and v["system"] != "ch4"
                                           and len(v["atoms"]) > 2]          # polyatomic structures only (no H2, CO, atoms)
    for system, sp in todo:
        if systems is not None and system not in systems:
            continue
        v = pack["structures"].get(key_of(system, sp, LV_OWN)) or pack["structures"].get(key_of(system, sp, LV_SCANPEAK))
        if v is None:
            continue
        at = [(s, tuple(c)) for s, c in v["atoms"]]
        def fn(at=at, v=v):
            mf = make_ks(at, v["charge"], v["basis"])
            e = float(mf.kernel())
            g = np.asarray(mf.nuc_grad_method().kernel(), float)
            rms, mx = float(np.sqrt(np.mean(g ** 2))), float(np.abs(g).max())
            return dict(energy=e, converged=bool(mf.converged), imag=v["imag"],
                        extra=dict(rms_grad=rms, max_grad=mx, tolerance=tol, passed=bool(rms < tol), pack_energy=v["energy_hartree"], dE_vs_pack=e - v["energy_hartree"],
                                   rms_grad_at_build=v["rms_grad"]))
        pid = sp
        if run is not None:
            row = run(csv_name, system, pid, v["level"], "DF-B3LYP-D3BJ-grad", v["basis"], fn, f"grad:{system}")
        else:
            row = run_point(csv_name, system, pid, v["level"], "DF-B3LYP-D3BJ-grad", v["basis"], fn, quiet=quiet)
        out.append((system, sp, row))
    return out


def grad_summary(csv_name="phase00_verify.csv", tol=GRAD_RMS_TOL):
    rows = [r for r in read_rows(csv_name) if r["method"] == "DF-B3LYP-D3BJ-grad"]
    n_pass = n_fail = 0
    worst = (0.0, "")
    for r in rows:
        x = extra_of(r)
        if x.get("rms_grad") is None or r["converged"] not in ("1", "True"):
            n_fail += 1; continue
        if x["rms_grad"] < tol:
            n_pass += 1
        else:
            n_fail += 1
        if x["rms_grad"] > worst[0]:
            worst = (x["rms_grad"], f"{r['system']}:{r['point_id']}")
    return dict(n=len(rows), passed=n_pass, failed=n_fail, worst_rms=worst[0], worst_where=worst[1], tolerance=tol)


# ------------------------------------------------------------------ REBUILD_GEOMS: the same builder, on Colab, into the working folder
def rebuild_geoms(systems, work_subdir="geompack_rebuild", pack_sha=None):
    """Not part of the 1-hour budget (extension EXT-geom, about 1.5 to 2 hours, resumable).  Builds every structure with the full-size machinery,
    writes the new geometries over the pack's (the pack versions are copied to geoms_pack/ first) and replaces the phase01 / phase02 / phase04 rows."""
    import shutil
    work = path(work_subdir)
    b = PackBuilder(work)
    b.build(systems)
    pack = b.assemble(systems, "rebuilt on this machine")
    d_old = path("geoms_pack")
    os.makedirs(d_old, exist_ok=True)
    for f in os.listdir(path("geoms")):
        if f.endswith((".json", ".xyz")) and not os.path.exists(os.path.join(d_old, f)):
            shutil.copy2(path("geoms", f), os.path.join(d_old, f))
    for name in (PHASE_CSV["phase01"], PHASE_CSV["phase02"], PHASE_CSV["phase04"]):
        if os.path.exists(path(name)):
            os.replace(path(name), path(name.replace(".csv", "_from_pack.csv")))
    for f in list(os.listdir(path("geoms"))):
        if f.endswith((".json", ".xyz")):
            os.remove(path("geoms", f))
    n = install_pack(pack, systems, sha=pack_sha_of(pack), hardware_label=hw())
    json.dump(pack, open(path("geompack_rebuilt.json"), "w"), indent=1, default=float)
    return n


def pack_sha_of(pack):
    return pack_sha256(pack)
