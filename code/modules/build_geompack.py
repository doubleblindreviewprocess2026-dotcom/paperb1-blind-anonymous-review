"""Builds the geometry pack of the 1-hour design in the sandbox (or anywhere with pyscf, geometric, tblite).

    python build_geompack.py --work DIR [--stage core|ext78] [--out geompack/b1v5_geompack_1h.json] [--smoke]

  --stage core   six reactions + CH4 (the pack the 1-hour notebooks ship)
  --stage ext78  F- + CH3Cl and OH- + CH3F (first extension items); added to the same work folder, pack written with --assemble
  --assemble     only write the pack file from the work folder (no calculation)
  --smoke        tiny bases (STO-3G / def2-SVP), HCN, NH3 and CH4 only: tests the machinery, the numbers are NOT the pack

Resumable: every finished step is in DIR/pack_work.json; run the same command again after an interruption.
Needs: pyscf 2.15.0, pyscf-dispersion 1.7.0, geometric 1.1.1, tblite 0.8.0 (the versions of the shipped pack)."""
import os, sys, json, argparse, hashlib, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ap = argparse.ArgumentParser()
ap.add_argument("--work", required=True)
ap.add_argument("--stage", default="core", choices=["core", "ext78", "none"])
ap.add_argument("--out", default=os.path.join(HERE, "geompack", "b1v5_geompack_1h.json"))
ap.add_argument("--assemble", action="store_true")
ap.add_argument("--smoke", action="store_true")
ap.add_argument("--systems", default="")
ap.add_argument("--regrid", default="", help="system:species,... repeat the Hessian at grid level 5 for structures already in the work folder")
a = ap.parse_args()
os.environ["B1V5_ROOT"] = os.path.abspath(a.work)
import b1v5_common as C
C.setup_root(os.path.abspath(a.work), mount=False)
C.start_clock(100)
import b1v5_cpu as P
import b1v5_geompack as G
for n in dir(C):
    if not n.startswith("__"):
        setattr(G, n, getattr(C, n)) if not hasattr(G, n) else None
for n in dir(P):
    if not n.startswith("__"):
        setattr(G, n, getattr(P, n)) if not hasattr(G, n) else None
if a.smoke:
    C.BASIS_OVERRIDE.update({"small": "STO-3G", "main": "def2-SVP", "large": "def2-TZVP"})
    G.CORE_SYSTEMS = ["hcn", "nh3", "ch4"]
P.quiet_geometric()
systems = a.systems.split(",") if a.systems else (G.CORE_SYSTEMS if a.stage == "core" else G.EXT78_SYSTEMS)
b = G.PackBuilder(a.work, hardware=C.hw())
for item in [x for x in a.regrid.split(",") if x]:
    b.regrid(*item.split(":"))
if not a.assemble and a.stage != "none":
    t0 = time.time()
    b.build(systems)
    print("build finished in %.1f min" % ((time.time() - t0) / 60))
allsys = list(G.CORE_SYSTEMS) + [s for s in G.EXT78_SYSTEMS if any(v["system"] == s for v in b.W["structures"].values())]
design = "1h core: six reactions + CH4" + (" + ext78 (F- + CH3Cl, OH- + CH3F)" if len(allsys) > len(G.CORE_SYSTEMS) else "")
pack = b.assemble(allsys, design)
os.makedirs(os.path.dirname(a.out), exist_ok=True)
data = G.pack_bytes(pack)
open(a.out, "wb").write(data)
sha = hashlib.sha256(data).hexdigest()
open(a.out + ".sha256", "w").write(f"{sha}  {os.path.basename(a.out)}\n")
print("wrote", a.out, len(data), "bytes, sha256", sha)
