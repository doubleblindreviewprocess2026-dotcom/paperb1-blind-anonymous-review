"""Tests of the two 1-hour notebooks that can be run WITHOUT a GPU and without Colab.

    python test_1h_notebooks.py [names ...]          (needs pyscf, geometric, tblite, ase, matplotlib, pandas, openpyxl, nbformat)
    B1V5_SMOKE_PACK=path/to/smoke_pack.json           a pack built with  build_geompack.py --smoke  (HCN, NH3, CH4 in tiny bases; built here if the variable is empty)

The notebooks are rebuilt from the smoke pack into a temporary folder (the real notebooks are not touched) and run cell by cell in one namespace
(run_cells.py), exactly as Colab does.  CPU notebook: B1V5_SMOKE=1 (small bases, HCN, NH3, CH4).  GPU notebook: B1V5_STANDIN=1 (cheap CPU functionals).

CPU tests
  cpu_full      one Run all: every phase done, pack verified, gradient test, queue section in RUN_SUMMARY.md, check ran with six reactions
  cpu_resume    second Run all recomputes nothing; a cut phase is completed without duplicates
  cpu_crash     exceptions injected into one system and into a whole phase: logged, later phases run; a hard stop (KeyboardInterrupt) in the middle of a
                phase, a fresh namespace, and the rerun finishes it
  cpu_deadline  a slow calculation kind and a small deadline: the slow phase stops at its cap, later phases still run, the notebook ends before the deadline
  cpu_admit     the admission test itself (deadline test, cap test, running mean, priors, refusal list, PhaseTimeout at the end of a phase with refusals)
  cpu_ext       EXTENSION = 'CPU1' on the finished folder; EXTENSION = 'GEOM' (REBUILD_GEOMS of one tiny system)
  cpu_final     the last cell with a mocked Colab runtime
GPU tests (stand-in mode)
  gpu_full      one Run all: phases, Skala rows with and without D3, G1a/G1b/G3/G2 rows, figures, tables, check, summary
  gpu_resume    nothing recomputed; cut phase completed
  gpu_crash     whole-phase and one-system crash injected: later phases (G5, check, summary) still run
  gpu_deadline  slow energies, small deadline: later phases refused ("not run (time)"), G5 and check still run, ends before the deadline
  gpu_tail      GPU class A100 (tail items 1 to 5), EXTENSION GPU1 and GPU2: items run or say why they were skipped (no geometries for reactions 7 and 8)
  gpu_nogpu     no GPU, no stand-in: only G0 diagnostics, everything else skipped with a message, summary written
  gpu_gate      failed G0 gate stops G1a to G2 but G5, the check and the summary run
  gpu_routes    functional-probe fallbacks with a faked GPU
  gpu_skalad3   Skala with and without D3 (explicit with_dftd3, loud failure of an old signature, noD3 row = same SCF minus D3, weights hash)
  gpu_align     density-fitting retry alignment through g1a_align
  gpu_fake      the non-stand-in code path with a FAKE GPU against the stored v4 energies (needs B1V5_CPU_SMOKE_ROOT = folder of a finished CPU smoke run)
  gpu_final     the last cell
"""
import os, sys, io, csv, json, glob, shutil, tempfile, types, time, hashlib, contextlib, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault("MPLBACKEND", "Agg")
from run_cells import run_cells, load_cells

RESULTS = []
NBDIR = None
CPU_NB = GPU_NB = None


def check(name, cond, detail=""):
    RESULTS.append((name, bool(cond)))
    print(("  PASS  " if cond else "  FAIL  ") + name + (f"   [{detail}]" if detail and not cond else ""), flush=True)


def fresh_root(tag):
    return tempfile.mkdtemp(prefix=f"b1v5_1h_{tag}_")


def prepare():
    """Build both notebooks from the smoke pack into a temporary folder."""
    global NBDIR, CPU_NB, GPU_NB
    if NBDIR:
        return
    pack = os.environ.get("B1V5_SMOKE_PACK")
    if not pack or not os.path.exists(pack):
        work = tempfile.mkdtemp(prefix="b1v5_smokepack_")
        pack = os.path.join(work, "pack.json")
        print("building the smoke pack (about 5 minutes) ...", flush=True)
        subprocess.run([sys.executable, os.path.join(HERE, "build_geompack.py"), "--work", work, "--smoke", "--out", pack], check=True, capture_output=True)
    import build_cpu_1h, build_gpu_1h
    NBDIR = tempfile.mkdtemp(prefix="b1v5_1h_nb_")
    CPU_NB, GPU_NB = os.path.join(NBDIR, "cpu.ipynb"), os.path.join(NBDIR, "gpu.ipynb")
    build_cpu_1h.build(CPU_NB, pack)
    build_gpu_1h.build(GPU_NB, pack)


def status(root):
    p = os.path.join(root, "phase_status.json")
    return json.load(open(p)) if os.path.exists(p) else {}


def rows(root, name):
    p = os.path.join(root, name)
    return list(csv.DictReader(open(p))) if os.path.exists(p) else []


def digest(root, names):
    return {n: hashlib.md5(open(os.path.join(root, n), "rb").read()).hexdigest() for n in names if os.path.exists(os.path.join(root, n))}


def no_dups(rs):
    keys = [(r["system"], r["point_id"], r["geometry_level"], r["method"], r["basis"]) for r in rs]
    return len(keys) == len(set(keys))


def done(st, ph):
    return st.get(ph, {}).get("status", "").startswith("done")


def run_cpu(root, hooks=None, only=None, ns=None, env=None):
    os.environ["B1V5_ROOT"] = root
    os.environ["B1V5_SMOKE"] = "1"
    os.environ.pop("B1V5_EXTENSION", None)
    for k, v in (env or {}).items():
        os.environ[k] = v
    return run_cells(CPU_NB, hooks=hooks, only=only, ns=ns, log=root + "_log.txt")


def run_gpu(root, standin=True, hooks=None, only=None, ns=None, env=None):
    os.environ["B1V5_ROOT"] = root
    os.environ["B1V5_STANDIN"] = "1" if standin else "0"
    os.environ.pop("B1V5_EXTENSION", None)
    for k, v in (env or {}).items():
        os.environ[k] = v
    return run_cells(GPU_NB, hooks=hooks, only=only, ns=ns, log=root + ("-STANDIN" if standin else "") + "_log.txt")


def groot(root, standin=True):
    return root + "-STANDIN" if standin else root


CPU_CSVS = ["phase00_repro.csv", "phase00_verify.csv", "phase01_stationary.csv", "phase02_scans.csv", "phase03_ccsdt.csv", "phase04_scout_ref.csv",
            "phase05_hcn_grid.csv", "phase06_ccsdt_basis.csv", "phase08_mp2_basis.csv"]
GPU_CSVS = ["phase10_repro.csv", "phase11_dft.csv", "phase12_basis.csv", "phase13_scout.csv", "phase14_skala_relaxed.csv"]
# code cells (0-based): 0 setup, 1 helpers, 2 check defs, 3 queue, 4 pack data, 5 settings, 6 pack, then the phases
CPU_CELL = dict(settings=5, pack=6, c0=7, c3a=8, c3b=9, c3c=10, c4=11, c6=12, tail=13, c7=14, summary=15, end=16)
GPU_CELL = dict(settings=5, pack=6, g0=7, g1a=8, g1b=9, g3=10, g2=11, tail=12, g5=13, check=14, summary=15, end=16)


# =============================================================================================== CPU
def test_cpu_full():
    print("== cpu_full ==")
    base = fresh_root("cpu_full")
    t0 = time.time()
    ns, tm = run_cpu(base)
    print(f"  wall time {time.time() - t0:.0f} s")
    st = status(base)
    check("no notebook cell raised an exception", all(s == "ok" for _, s, _ in tm), str([t for t in tm if t[1] != "ok"]))
    for ph in ("phase00", "phase03a", "phase03b", "phase08", "phase04", "phase06", "phase08t", "phase05", "phase07x_check"):
        check(f"{ph} finished", done(st, ph), st.get(ph, {}).get("status"))
    for f in CPU_CSVS:
        if os.path.exists(os.path.join(base, f)):
            check(f"{f}: no duplicate rows", no_dups(rows(base, f)))
    info = json.load(open(os.path.join(base, "geompack_installed.json")))
    check("pack hash verified and installed", info.get("ok") is True and info.get("sha256") == info.get("expected_sha256"), str(info.get("problems")))
    ver = json.load(open(os.path.join(base, "phase00_verify_summary.json")))
    check("gradient test passed for all structures", ver["gradient"]["n"] >= 4 and ver["gradient"]["failed"] == 0, str(ver["gradient"]))
    check("speed factor measured", 0.7 <= ver["speed_factor"] <= 8)
    gr = [r for r in rows(base, "phase00_verify.csv")]
    check("gradient rows are in the long CSV format with the gradient in extra", gr and all("rms" in r["extra"] for r in gr))
    check("stored v4 CCSD(T) point 0 recomputed (phase00_repro.csv)", any(r["method"] == "v4-ccsd_t" for r in rows(base, "phase00_repro.csv")))
    cc = rows(base, "phase03_ccsdt.csv")
    check("CCSD(T) rows with T1 at the stationary points and along the paths", any(r["point_id"] == "TS" and r["t1"] for r in cc) and any(r["point_id"].startswith("p") for r in cc))
    check("MP2 rows in two bases", len({r["basis"] for r in rows(base, "phase08_mp2_basis.csv")}) >= 2)
    summ = open(os.path.join(base, "RUN_SUMMARY.md")).read()
    check("RUN_SUMMARY has the queue section, the pack hash and the gradient test", "## Deadline queue" in summ and "Geometry pack" in summ and "Gradient test" in summ)
    check("RUN_SUMMARY budget line shows the 55 minute window, not 3.2 h", "budget 0.92 h" in summ, summ[:400])
    rep = open(os.path.join(base, "check_report.md")).read()
    check("check report is the 1-hour version with NOT RUN for the extension reactions", "1-hour design" in rep and "NOT RUN" in rep and "FAIL | R1 | L11" not in rep)
    check("check items P1 and P2 present and PASS", "| PASS | P1 |" in rep and "| PASS | P2 |" in rep)
    check("queue files written", os.path.exists(os.path.join(base, "queue_log.json")) and os.path.exists(os.path.join(base, "queue_not_run.csv")))
    ql = json.load(open(os.path.join(base, "queue_log.json")))["sessions"][-1]
    check("the running means per kind of calculation were recorded (admission test uses them)", len(ql["means"]) >= 5 and all(v > 0 for v in ql["means"].values()), str(ql["means"])[:200])
    check("tiers.json written", os.path.exists(os.path.join(base, "tiers.json")))
    check("config.json says design 1h and smoke", json.load(open(os.path.join(base, "config.json"))).get("design") == "1h")
    return base


def test_cpu_resume(base):
    print("== cpu_resume ==")
    before = digest(base, CPU_CSVS)
    ns, tm = run_cpu(base)
    after = digest(base, CPU_CSVS)
    st = status(base)
    check("second Run all: every core phase 'already done'", all(st[p]["status"] == "already done" for p in ("phase00", "phase03a", "phase03b", "phase08", "phase04", "phase06", "phase08t", "phase05")), str({k: v["status"] for k, v in st.items()}))
    check("second Run all: CSV files byte-identical", before == after, str([k for k in before if before[k] != after.get(k)]))
    check("second Run all: check and summary rerun without error", all(s == "ok" for _, s, _ in tm))
    # cut phase 3b in the middle: remove the DONE marker and the last 12 rows
    p = os.path.join(base, "phase03_ccsdt.csv")
    lines = open(p).read().splitlines()
    open(p, "w").write("\n".join(lines[:-12]) + "\n")
    os.remove(os.path.join(base, "DONE_phase03b"))
    ns, tm = run_cpu(base)
    now_rows = rows(base, "phase03_ccsdt.csv")
    check("cut phase: missing rows recomputed, none duplicated", len(now_rows) == len(lines) - 1 and no_dups(now_rows), f"{len(now_rows)} vs {len(lines) - 1}")
    check("cut phase: kept rows byte-identical", open(p).read().splitlines()[:len(lines) - 12] == lines[:-12])
    # the geometry pack is re-installed idempotently: no duplicate phase01 rows
    check("rerun: pack install did not duplicate rows", no_dups(rows(base, "phase01_stationary.csv")) and no_dups(rows(base, "phase02_scans.csv")))


def test_cpu_crash():
    print("== cpu_crash ==")
    base = fresh_root("cpu_crash")

    def patch(ns):
        orig_cc = ns["ccsdt"]

        def ccsdt(atoms, ch, basis, **k):                # one system's CCSD(T) crashes
            if len(atoms) == 4 and ch == 0 and basis == ns["basis_for"]("nh3", "main") and not k.get("frozen_core") is False and any(a[0] == "N" for a in atoms):
                raise RuntimeError("INJECTED crash in the NH3 CCSD(T)")
            return orig_cc(atoms, ch, basis, **k)
        ns["ccsdt"] = ccsdt

        def boom(*a, **k):
            raise RuntimeError("INJECTED crash of the whole C4 phase")
        ns["c4_scout_reference"] = boom
    ns, tm = run_cpu(base, hooks={CPU_CELL["c0"]: patch})
    st = status(base)
    errs = rows(base, "phase_errors.csv")
    check("all cells ran (no exception escaped)", all(s == "ok" for _, s, _ in tm), str([t for t in tm if t[1] != "ok"]))
    check("system crash stored as failed rows and logged for NH3", any("INJECTED" in r["extra"] and r["system"] == "nh3" for r in rows(base, "phase03_ccsdt.csv")) or any("nh3" in e["system"] and "INJECTED" in e["error"] for e in errs))
    check("whole-phase crash logged for phase04", any(e["phase"] == "phase04" and "INJECTED" in e["error"] for e in errs))
    check("phase04 FAILED", st["phase04"]["status"] == "FAILED", st["phase04"]["status"])
    check("HCN CCSD(T) rows survived", any(r["system"] == "hcn" and r["converged"] in ("1", "True") for r in rows(base, "phase03_ccsdt.csv")))
    for ph in ("phase06", "phase08t", "phase05", "phase07x_check"):
        check(f"later phase {ph} still ran", done(st, ph), st.get(ph, {}).get("status"))
    summ = open(os.path.join(base, "RUN_SUMMARY.md")).read()
    check("RUN_SUMMARY lists the failed phase and quotes the errors", "FAILED" in summ and "INJECTED" in summ)
    check("check ran on the partial data", os.path.exists(os.path.join(base, "check_report.md")))
    # a second run without the injected crash and with RETRY_FAILED finishes everything
    def retry(ns):
        ns["RETRY_FAILED"] = True
        ns["FORCE_RERUN"] = True            # DONE markers are ignored (finished points are still skipped)
    ns, tm = run_cpu(base, hooks={CPU_CELL["c0"]: retry})
    st = status(base)
    check("rerun without the crash: phase04 done", done(st, "phase04"), st["phase04"]["status"])
    bad = [r for r in rows(base, "phase03_ccsdt.csv") if r["system"] == "nh3" and r["converged"] not in ("1", "True")]
    check("rerun with RETRY_FAILED: the NH3 points are recomputed", not bad or all(any(x["point_id"] == r["point_id"] and x["converged"] in ("1", "True") for x in rows(base, "phase03_ccsdt.csv")) for r in bad), str(len(bad)))
    # hard stop in the middle of a phase (runtime killed): the next session is a FRESH namespace
    base2 = fresh_root("cpu_kill")
    count = {"n": 0}

    def killer(ns):
        orig = ns["ccsdt"]

        def ccsdt(*a, **k):
            count["n"] += 1
            if count["n"] == 6:
                raise KeyboardInterrupt("INJECTED hard stop (runtime killed)")
            return orig(*a, **k)
        ns["ccsdt"] = ccsdt
    ns, tm = run_cpu(base2, hooks={CPU_CELL["c0"]: killer}, only=list(range(0, CPU_CELL["c3b"] + 1)))
    n1 = len(rows(base2, "phase03_ccsdt.csv"))
    check("hard stop: the cell stopped, rows written so far are kept", any(s != "ok" for _, s, _ in tm) and n1 >= 4, f"{n1} rows; {tm[-3:]}")
    st1 = status(base2)
    check("hard stop: the interrupted phase is recorded as 'interrupted' and has no DONE marker (not a Exception: the original phase() marked it DONE)",
          any(v["status"] == "interrupted" for v in st1.values()) and not all(os.path.exists(os.path.join(base2, "DONE_" + k)) for k, v in st1.items() if k in ("phase03a",)), str({k: v["status"] for k, v in st1.items()}))
    ns, tm = run_cpu(base2)           # fresh namespace
    st = status(base2)
    check("fresh session after the hard stop: everything finishes", all(s == "ok" for _, s, _ in tm) and done(st, "phase03a") and st["phase03b"]["status"] in ("done", "already done") and done(st, "phase07x_check"), str({k: v["status"] for k, v in st.items()}))
    check("fresh session: no duplicate rows, rows of the first session kept", no_dups(rows(base2, "phase03_ccsdt.csv")) and len(rows(base2, "phase03_ccsdt.csv")) > n1)


def test_cpu_admit():
    print("== cpu_admit (admission test of the queue) ==")
    base = fresh_root("cpu_admit")
    ns, tm = run_cpu(base, only=list(range(0, CPU_CELL["pack"] + 1)))
    Q, q_admit, q_start, q_record, qrun, qphase = ns["Q"], ns["q_admit"], ns["q_start"], ns["q_record"], ns["qrun"], ns["qphase"]
    q_start(10.0, t0=time.time() - 60 * 7.0, final_min=1.0)          # 7 of 10 minutes used
    Q["means"].clear(); Q["counts"].clear(); Q["scale"] = 1.0; Q["prior"]["t_a"] = 60.0; Q["prior"]["t_b"] = 100.0
    ok, why = q_admit("t_a")
    check("deadline test: 7 min used + 1.5 x 60 s = 8.5 min < 10: admitted", ok, why)
    ok, why = q_admit("t_b")
    check("deadline test: 7 min used + 1.5 x 100 s = 9.5 min < 10: admitted", ok, why)
    Q["prior"]["t_c"] = 130.0
    ok, why = q_admit("t_c")
    check("deadline test: 7 min + 1.5 x 130 s = 10.25 min > 10: refused, with the reason", (not ok) and "deadline" in why, why)
    Q["scale"] = 2.0
    ok, why = q_admit("t_b")
    check("the speed factor scales the prior: 1.5 x 200 s: refused", not ok)
    Q["scale"] = 1.0
    q_record("t_c", 30.0); q_record("t_c", 50.0)
    check("a measured running mean replaces the prior (t_c now 40 s)", abs(ns["q_est"]("t_c") - 40.0) < 1e-9 and q_admit("t_c")[0])
    check("kinds fall back to the prefix before ':' (prior of 't_a' for 't_a:hcn:TS')", abs(ns["q_est"]("t_a:hcn:TS") - 60.0) < 1e-9)
    # cap test and refusal bookkeeping, phase left not DONE
    q_start(50.0, t0=time.time(), final_min=5.0)
    Q["refused"].clear()
    calls = []
    try:
        with qphase("phase_tst", "admission test phase", 0.05) as ph:         # cap 3 s
            Q["prior"]["t_slow"] = 10.0
            r = qrun("t_tst.csv", "x", "p0", "lv", "m", "b", lambda: calls.append(1) or 1.0, "t_slow")
    except ns["PhaseTimeout"] as ex:
        check("a phase with refused calculations ends with PhaseTimeout (not DONE)", "not run (time)" in str(ex), str(ex))
    st = status(base).get("phase_tst", {})
    check("the refused calculation was not started", calls == [] and len(Q["refused"]) == 1 and "phase cap" in Q["refused"][0]["reason"], str(Q["refused"]))
    check("phase status 'stopped', no DONE marker", st.get("status", "").startswith("stopped") and not os.path.exists(os.path.join(base, "DONE_phase_tst")), str(st))
    ns["q_save_log"]()
    check("queue_not_run.csv lists the refusal", "t_slow" in open(os.path.join(base, "queue_not_run.csv")).read())
    check("RUN_SUMMARY queue lines mention the refusal", any("not run" in l.lower() for l in ns["queue_lines"]()))
    # hard stop: after deadline + final minutes check_time raises, except in a final phase
    q_start(0.01, t0=time.time() - 60.0, final_min=0.01)
    try:
        ns["check_time"](); raised = False
    except ns["PhaseTimeout"]:
        raised = True
    check("check_time raises after the hard stop", raised)
    Q["final"] = True
    try:
        ns["check_time"](); ok = True
    except Exception:
        ok = False
    Q["final"] = False
    check("check_time does not raise in a final phase (figures, check, summary always run)", ok)
    # q_time_limit interrupts one calculation
    t0 = time.time()
    try:
        with ns["q_time_limit"](1.0):
            while True:
                sum(range(10000))
        hit = False
    except ns["QTimeout"]:
        hit = True
    check("q_time_limit interrupts a long calculation after about 1 s", hit and time.time() - t0 < 4, f"{time.time() - t0:.1f} s")


def test_cpu_deadline():
    print("== cpu_deadline (a slow phase must not take the later phases with it) ==")
    base = fresh_root("cpu_dl")
    t0 = time.time()
    DEADLINE, FINAL, SLOW = 3.5, 0.5, 10.0

    def patch(ns):
        ns["q_start"](DEADLINE, t0=ns["_T0"], final_min=FINAL)
        ns["CPU_CAPS"]["phase03a"] = 0.4
        ns["CPU_CAPS"]["phase03b"] = 0.3
        orig = ns["ccsdt"]

        def ccsdt(*a, **k):                                  # a slow CCSD(T): every call takes 14 s more
            time.sleep(SLOW)
            return orig(*a, **k)
        ns["ccsdt"] = ccsdt
    ns, tm = run_cpu(base, hooks={CPU_CELL["c0"]: patch})
    wall = time.time() - t0
    st = status(base)
    nr = rows(base, "queue_not_run.csv")
    check("all cells ran", all(s == "ok" for _, s, _ in tm), str([t for t in tm if t[1] != "ok"]))
    check("C3a stopped at its cap and is not DONE", st["phase03a"]["status"].startswith("stopped") and not os.path.exists(os.path.join(base, "DONE_phase03a")), st["phase03a"]["status"])
    check("calculations were refused and listed (queue_not_run.csv)", len(nr) > 0 and any(r["phase"] == "phase03a" for r in nr), str(len(nr)))
    check("CCSD(T) rows of the first systems exist", any(r["system"] == "hcn" for r in rows(base, "phase03_ccsdt.csv")))
    check("MP2 phase (a different, cheap kind) still ran after the slow phase", st.get("phase08", {}).get("status", "").startswith(("done", "stopped")) and len(rows(base, "phase08_mp2_basis.csv")) > 0, str(st.get("phase08")))
    check("check and summary ran after the deadline", done(st, "phase07x_check") and os.path.exists(os.path.join(base, "RUN_SUMMARY.md")))
    summ = open(os.path.join(base, "RUN_SUMMARY.md")).read()
    check("RUN_SUMMARY says how many calculations were not run (time)", "Not run (time) in this session" in summ and "none: every admitted" not in summ)
    # the notebook (without the final 8 s sleep of the release cell) ends before deadline + final + the longest admitted calculation
    elapsed_min = (time.time() - ns["_T0"]) / 60.0
    check(f"ended at minute {elapsed_min:.1f}, before deadline {DEADLINE} + final {FINAL} + one slow calculation", elapsed_min < DEADLINE + FINAL + (SLOW + 25.0) / 60.0, f"{elapsed_min:.2f}")
    check("no calculation was started after the deadline", all(float(r["minute"]) <= DEADLINE + 0.5 for r in nr if r["minute"]))
    # a later full-time session completes what was not run
    ns, tm = run_cpu(base)
    st = status(base)
    check("second Run all with the normal clock completes C3a and the rest", done(st, "phase03a") and done(st, "phase03b"), str({k: v["status"] for k, v in st.items()}))
    check("no duplicate rows after completion", no_dups(rows(base, "phase03_ccsdt.csv")))


def test_cpu_ext(base=None):
    print("== cpu_ext ==")
    base = base or test_cpu_full_quiet()
    st0 = status(base)
    ns, tm = run_cpu(base, env={"B1V5_EXTENSION": "CPU1"})
    st = status(base)
    check("EXTENSION CPU1: all cells ran", all(s == "ok" for _, s, _ in tm), str([t for t in tm if t[1] != "ok"]))
    check("EXTENSION CPU1: the core phases were skipped (not rerun), check and summary ran", done(st, "phase07x_check") and ns["RUN_CORE"] is False)
    errs = open(os.path.join(base, "phase_errors.csv")).read() if os.path.exists(os.path.join(base, "phase_errors.csv")) else ""
    check("EXTENSION CPU1: Cl- + CH3Cl without geometries (smoke pack) says why instead of crashing", "phase03x" in st and (done(st, "phase03x") or "skipped" in errs or st["phase03x"]["status"] != "FAILED"), str(st.get("phase03x")))
    # EXTENSION GEOM with a tiny system
    base2 = fresh_root("cpu_geom")

    def tiny(ns):
        ns["SYSTEMS_TO_RUN"] = ["hcn"]
    os.environ["B1V5_ROOT"] = base2
    ns, tm = run_cpu(base2, hooks={CPU_CELL["pack"]: tiny}, env={"B1V5_EXTENSION": "GEOM"}, only=list(range(0, CPU_CELL["pack"] + 1)))
    check("EXTENSION GEOM: REBUILD_GEOMS set, cells ran", all(s == "ok" for _, s, _ in tm) and ns["REBUILD_GEOMS"] is True, str([t for t in tm if t[1] != "ok"]))
    check("EXTENSION GEOM: geometries rebuilt on this machine (geoms/hcn_HCN_own_b3lyp.json)", os.path.exists(os.path.join(base2, "geoms", "hcn_HCN_own_b3lyp.json")) and os.path.exists(os.path.join(base2, "geompack_rebuilt.json")))
    check("EXTENSION GEOM: no pack-installed record", not os.path.exists(os.path.join(base2, "geompack_installed.json")))


_FULL_CACHE = {}


def test_cpu_full_quiet():
    if "base" not in _FULL_CACHE:
        _FULL_CACHE["base"] = test_cpu_full()
    return _FULL_CACHE["base"]


def test_final_cell(nb, tag):
    print(f"== {tag}_final ==")
    cells = load_cells(nb)
    src = cells[-1][0]
    for case in ("on", "off", "crash"):
        root = fresh_root("final")
        called = []
        mod_colab = types.ModuleType("google.colab"); mod_google = types.ModuleType("google")
        rt = types.ModuleType("google.colab.runtime")
        if case == "crash":
            rt.unassign = lambda: (_ for _ in ()).throw(RuntimeError("mock: cannot unassign"))
        else:
            rt.unassign = lambda: called.append(time.time())
        mod_colab.runtime = rt
        saved = {k: sys.modules.get(k) for k in ("google", "google.colab", "google.colab.runtime")}
        sys.modules.update({"google": mod_google, "google.colab": mod_colab, "google.colab.runtime": rt})
        os.environ["B1V5_ROOT"] = root
        os.environ["B1V5_AUTO_DISCONNECT"] = "0" if case == "off" else "1"
        json.dump({"phase11": {"status": "done"}}, open(os.path.join(root, "phase_status.json"), "w"))
        t0 = time.time(); buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                exec(compile(src, "<final>", "exec"), {"__name__": "__main__"})        # EMPTY namespace: proves independence of the earlier cells
            ok = True
        except Exception:
            ok = False
        dt = time.time() - t0
        for k, v in saved.items():
            if v is None:
                sys.modules.pop(k, None)
            else:
                sys.modules[k] = v
        if case == "on":
            check(f"{tag}: AUTO_DISCONNECT on: runtime.unassign called exactly once, in an empty namespace", ok and len(called) == 1, f"calls {len(called)}")
            check(f"{tag}: it waited a few seconds first (Drive flush)", dt >= 7.5, f"{dt:.1f} s")
            check(f"{tag}: fallback RUN_SUMMARY.md written when none exists", os.path.exists(os.path.join(root, "RUN_SUMMARY.md")))
        elif case == "off":
            check(f"{tag}: AUTO_DISCONNECT off: not called", ok and len(called) == 0)
        else:
            check(f"{tag}: unassign raising an error does not raise out of the cell", ok and "could not be released" in buf.getvalue())
        shutil.rmtree(root, ignore_errors=True)


# =============================================================================================== GPU (stand-in)
def test_gpu_full():
    print("== gpu_full ==")
    base = fresh_root("gpu_full")
    t0 = time.time()
    ns, tm = run_gpu(base)
    print(f"  wall time {time.time() - t0:.0f} s")
    R = groot(base)
    st = status(R)
    print("  ", {k: v["status"] for k, v in sorted(st.items())})
    check("no notebook cell raised an exception", all(s == "ok" for _, s, _ in tm), str([t for t in tm if t[1] != "ok"]))
    for ph in ("phase10", "phase11", "phase11p", "phase13", "phase12", "phase15", "phase15x_check"):
        check(f"{ph} finished", done(st, ph), st.get(ph, {}).get("status"))
    for f in GPU_CSVS[:4]:
        check(f"{f} has rows, no duplicates", len(rows(R, f)) > 0 and no_dups(rows(R, f)))
    r11 = rows(R, "phase11_dft.csv")
    meth = sorted({r["method"] for r in r11})
    check("G1 methods: Skala (with D3), noD3, anchor, B3LYP, r2SCAN, PBE", len(meth) == 6 and "Skala-1.1-noD3" in meth and "Skala-1.1" in meth, str(meth))
    sk = {(r["system"], r["point_id"], r["geometry_level"], r["basis"]): r for r in r11 if r["method"] == "Skala-1.1"}
    nd = {(r["system"], r["point_id"], r["geometry_level"], r["basis"]): r for r in r11 if r["method"] == "Skala-1.1-noD3"}
    check("every Skala row has a noD3 partner", set(sk) == set(nd) and len(sk) > 0, f"{len(sk)} vs {len(nd)}")
    bad = [k for k in sk if sk[k]["converged"] in ("1", "True") and abs(float(sk[k]["energy_hartree"]) - float(json.loads(sk[k]["extra"])["d3_hartree"]) - float(nd[k]["energy_hartree"])) > 1e-12]
    check("noD3 = Skala minus the stored D3 energy (same SCF)", not bad, str(bad[:2]))
    # path rows: three functionals only (the core), PBE and r2SCAN not along the paths
    paths = [r for r in r11 if r["point_id"].startswith("p")]
    pm = sorted({r["method"] for r in paths})
    check("paths: Skala, noD3, wB97M-V and B3LYP only in the core (PBE and r2SCAN are tail item 2)", pm == ["B3LYP-D3BJ", "Skala-1.1", "Skala-1.1-noD3", "wB97M-V"], str(pm))
    check("path systems: hcn, nh3, ch4", {r["system"] for r in paths} == {"hcn", "nh3", "ch4"}, str({r["system"] for r in paths}))
    stat = [r for r in r11 if not r["point_id"].startswith("p")]
    check("stationary points: all five functionals", len({r["method"] for r in stat}) == 6)
    sc = rows(R, "phase13_scout.csv")
    check("G3 core: three functionals (+ noD3) at the xTB geometries", sorted({r["method"] for r in sc}) == ["B3LYP-D3BJ", "Skala-1.1", "Skala-1.1-noD3", "wB97M-V"], str(sorted({r["method"] for r in sc})))
    b12 = rows(R, "phase12_basis.csv")
    check("G2 core: Skala and anchor only", sorted({r["method"] for r in b12}) == ["Skala-1.1", "Skala-1.1-noD3", "wB97M-V"], str(sorted({r["method"] for r in b12})))
    check("G0: stored v4 HCN point 0 rerun with three methods", len(rows(R, "phase10_repro.csv")) == 3)
    ql = json.load(open(os.path.join(R, "queue_log.json")))["sessions"][-1]
    check("GPU: running means per kind recorded", len(ql["means"]) >= 4, str(ql["means"])[:200])
    check("no gradient test and no G4 in the core", not os.path.exists(os.path.join(R, "phase14_skala_relaxed.csv")))
    for f in ("fig1_design.png", "fig2_barrier_errors.png", "fig4_basis_shifts.png", "fig5_scout_check.png", "fig6_hcn_map_ch4.png"):
        check(f"figures/{f} written", os.path.exists(os.path.join(R, "figures", f)) and os.path.getsize(os.path.join(R, "figures", f)) > 3000)
    for f in ("numbers_v5.csv", "summary_tables.xlsx", "table1_systems.csv", "table2_errors.csv", "check_report.md", "RUN_SUMMARY.md", "phase10_timing.json", "findings_outcome.csv",
              "basis_table.csv", "scout_table.csv", "geometry_table.csv", "skala_d3_table.csv", "reference_quality.csv", "phase13_xtb_barriers.json", "queue_log.json"):
        check(f"{f} written", os.path.exists(os.path.join(R, f)))
    fo = rows(R, "findings_outcome.csv")
    check("findings_outcome.csv: F1 to F6, frozen date, verdict (NOT RUN where the data is missing)", [r["finding"] for r in fo] == ["F1", "F2", "F3", "F4", "F5", "F6"] and all(r["frozen_on"].startswith("2026-11-15") and r["verdict"] for r in fo), str([r["verdict"] for r in fo]))
    rep = open(os.path.join(R, "check_report.md")).read()
    check("check: 1-hour design, extension items NOT RUN, I11 dropped by D2", "1-hour design" in rep and "dropped by design" in rep and "| NOT RUN | R1 | L11" in rep)
    summ = open(os.path.join(R, "RUN_SUMMARY.md")).read()
    check("RUN_SUMMARY: stand-in marked, GPU info, queue section, pack hash", "STAND-IN" in summ and "## Deadline queue" in summ and "Geometry pack" in summ)
    return base


def test_gpu_resume(base):
    print("== gpu_resume ==")
    R = groot(base)
    before = digest(R, GPU_CSVS)
    ns, tm = run_gpu(base)
    after = digest(R, GPU_CSVS)
    st = status(R)
    check("second Run all: core phases 'already done' (or re-entered and 'done' with nothing recomputed, FORCE_RERUN is on in the second-run notebook)", all(st[p]["status"] in ("already done", "done") for p in ("phase10", "phase11", "phase11p", "phase13", "phase12")), str({k: v["status"] for k, v in st.items()}))
    check("second Run all: CSV files byte-identical", before == after, str([k for k in before if before[k] != after.get(k)]))
    p = os.path.join(R, "phase11_dft.csv")
    lines = open(p).read().splitlines()
    open(p, "w").write("\n".join(lines[:-30]) + "\n")
    for ph in ("phase11", "phase11p"):
        if os.path.exists(os.path.join(R, "DONE_" + ph)):
            os.remove(os.path.join(R, "DONE_" + ph))
    ns, tm = run_gpu(base)
    now_rows = rows(R, "phase11_dft.csv")
    check("cut phase: rows recomputed, none duplicated", len(now_rows) == len(lines) - 1 and no_dups(now_rows), f"{len(now_rows)} vs {len(lines) - 1}")
    check("cut phase: kept rows byte-identical", open(p).read().splitlines()[:len(lines) - 30] == lines[:-30])


def test_gpu_crash():
    print("== gpu_crash ==")
    base = fresh_root("gpu_crash")

    def patch(ns):
        def boom(*a, **k):
            raise RuntimeError("INJECTED crash of the whole G2 phase")
        ns["g2_basis_1h"] = boom
        orig = ns["energy_of"]

        def energy_of(method, atoms, charge, basis, **k):
            if len(atoms) == 4 and any(a[0] == "N" for a in atoms) and method == "PBE":      # NH3 (and its TS), not HCN
                raise RuntimeError("INJECTED crash in NH3 PBE")
            return orig(method, atoms, charge, basis, **k)
        ns["energy_of"] = energy_of
    ns, tm = run_gpu(base, hooks={GPU_CELL["g0"]: patch})
    R = groot(base)
    st = status(R)
    errs = rows(R, "phase_errors.csv")
    check("all cells ran (no exception escaped)", all(s == "ok" for _, s, _ in tm), str([t for t in tm if t[1] != "ok"]))
    check("whole-phase crash logged for phase12 and the phase is FAILED", any(e["phase"] == "phase12" and "INJECTED" in e["error"] for e in errs) and st["phase12"]["status"] == "FAILED", st["phase12"]["status"])
    check("single crash stored for NH3 PBE", any("INJECTED" in r["extra"] and r["system"] == "nh3" and r["method"] == "PBE" for r in rows(R, "phase11_dft.csv")))
    check("HCN PBE rows survived", any(r["system"] == "hcn" and r["method"] == "PBE" and r["converged"] in ("1", "True") for r in rows(R, "phase11_dft.csv")))
    for ph in ("phase11p", "phase13", "phase15", "phase15x_check"):
        check(f"later phase {ph} still ran", done(st, ph) or "error" in st.get(ph, {}).get("status", ""), st.get(ph, {}).get("status"))
    summ = open(os.path.join(R, "RUN_SUMMARY.md")).read()
    check("RUN_SUMMARY lists the failed phase and quotes the error", "FAILED" in summ and "INJECTED" in summ)


def test_gpu_deadline():
    print("== gpu_deadline ==")
    base = fresh_root("gpu_dl")
    DEADLINE, FINAL, SLOW = 1.6, 0.5, 6.0

    def patch(ns):
        ns["q_start"](DEADLINE, t0=ns["_T0"], final_min=FINAL)
        ns["GPU_CAPS"]["phase11"] = 0.7
        orig = ns["energy_of"]

        def energy_of(*a, **k):
            time.sleep(SLOW)
            return orig(*a, **k)
        ns["energy_of"] = energy_of
    ns, tm = run_gpu(base, hooks={GPU_CELL["g0"]: patch})
    R = groot(base)
    st = status(R)
    nr = rows(R, "queue_not_run.csv")
    check("all cells ran", all(s == "ok" for _, s, _ in tm), str([t for t in tm if t[1] != "ok"]))
    check("G1a stopped at its cap, not DONE", st["phase11"]["status"].startswith("stopped") and not os.path.exists(os.path.join(R, "DONE_phase11")), st["phase11"]["status"])
    check("later phases refused whole items and listed them", len(nr) > 0 and {r["phase"] for r in nr} & {"phase11", "phase11p", "phase13", "phase12"}, str(len(nr)))
    check("G5, the check and the summary still ran", done(st, "phase15") and done(st, "phase15x_check") and os.path.exists(os.path.join(R, "RUN_SUMMARY.md")), str({k: v["status"] for k, v in st.items()}))
    check("figures and the outcome table were made on the partial data", os.path.exists(os.path.join(R, "findings_outcome.csv")) and os.path.exists(os.path.join(R, "figures", "fig1_design.png")))
    summ = open(os.path.join(R, "RUN_SUMMARY.md")).read()
    check("RUN_SUMMARY states how many calculations were not run", "Not run (time) in this session" in summ and "none: every admitted" not in summ)
    elapsed = (time.time() - ns["_T0"]) / 60.0
    check(f"ended at minute {elapsed:.1f}, before deadline {DEADLINE} + final {FINAL} + one slow calculation", elapsed < DEADLINE + FINAL + (SLOW + 30) / 60.0, f"{elapsed:.2f}")
    fo = rows(R, "findings_outcome.csv")
    check("findings with missing data say NOT RUN, never FAIL", all(r["verdict"] in ("supports", "negative", "inconclusive") or r["verdict"].startswith("NOT RUN") for r in fo), str([r["verdict"] for r in fo]))
    ns, tm = run_gpu(base)
    st = status(R)
    check("a second Run all completes the phases", done(st, "phase11") and done(st, "phase11p") and done(st, "phase12"), str({k: v["status"] for k, v in st.items()}))
    check("no duplicate rows after completion", no_dups(rows(R, "phase11_dft.csv")))


def test_gpu_tail():
    print("== gpu_tail ==")
    # (a) A100 class on the stand-in: tail items 1 to 5 are tried
    base = fresh_root("gpu_tail")

    def a100(ns):
        ns["GPU_CLASS"] = "A100"
    ns, tm = run_gpu(base, hooks={GPU_CELL["tail"]: a100})
    R = groot(base)
    st = status(R)
    check("A100: tail items 1 to 5 chosen by class", ns["_tail"] == [1, 2, 3, 4, 5], str(ns["_tail"]))
    check("all cells ran", all(s == "ok" for _, s, _ in tm), str([t for t in tm if t[1] != "ok"]))
    errs = rows(R, "phase_errors.csv")
    check("tail 1: reactions 7 and 8 without geometries in the pack: skipped with a message, no crash", any(e["phase"] == "tail1_ext78" and "skipped" in e["error"] for e in errs) and st.get("tail1_ext78", {}).get("status", "") != "FAILED", str(st.get("tail1_ext78")))
    r11 = rows(R, "phase11_dft.csv")
    paths_m = {r["method"] for r in r11 if r["point_id"].startswith("p")}
    check("tail 2: PBE and r2SCAN along the paths", {"PBE", "r2SCAN"} <= paths_m, str(paths_m))
    sc = {r["method"] for r in rows(R, "phase13_scout.csv")}
    check("tail 4: PBE and r2SCAN at the xTB geometries", {"PBE", "r2SCAN"} <= sc, str(sc))
    b12 = {(r["system"], r["method"]) for r in rows(R, "phase12_basis.csv")}
    check("tail 3: B3LYP added to the basis check", any(m.startswith("B3LYP") for _, m in b12), str(sorted(b12)[:6]))
    g4 = rows(R, "phase14_skala_relaxed.csv")
    check("tail 5: the Skala gradient test ran and Skala-relaxed HCN was computed", ns["GSTATE"].get("grad_ok") is True and any(r["system"] == "hcn" for r in g4), str(ns["GSTATE"].get("grad_ok")))
    check("G5 and the check ran after the tail", done(st, "phase15") and done(st, "phase15x_check"))
    check("RUN_SUMMARY names the tail items", "tail items tried: [1, 2, 3, 4, 5]" in open(os.path.join(R, "RUN_SUMMARY.md")).read())
    # (b) T4 class: nothing in the tail
    check("T4, L4, A100 tail lists as in the plan", ns["tail_items_for"](None, "T4") == [] and ns["tail_items_for"](None, "L4") == [1, 2] and ns["tail_items_for"](None, "A100") == [1, 2, 3, 4, 5]
          and ns["tail_items_for"]("GPU1", "T4") == [1, 2, 3] and ns["tail_items_for"]("GPU2", "T4") == [4, 5])
    # (c) extension sessions on a finished core folder: the core is not rerun
    base2 = fresh_root("gpu_ext")
    ns, tm = run_gpu(base2)
    R2 = groot(base2)
    n_before = len(rows(R2, "phase11_dft.csv"))
    ns, tm = run_gpu(base2, env={"B1V5_EXTENSION": "GPU1"})
    check("EXT-GPU-1: items 1 to 3 run, the core is not rerun, all cells ok", ns["_tail"] == [1, 2, 3] and ns["RUN_CORE"] is False and all(s == "ok" for _, s, _ in tm), str(ns["_tail"]))
    r11 = rows(R2, "phase11_dft.csv")
    check("EXT-GPU-1: PBE and r2SCAN path rows added, no duplicates", {"PBE", "r2SCAN"} <= {r["method"] for r in r11 if r["point_id"].startswith("p")} and no_dups(r11) and len(r11) > n_before)
    ns, tm = run_gpu(base2, env={"B1V5_EXTENSION": "GPU2"})
    check("EXT-GPU-2: items 4 and 5 run", ns["_tail"] == [4, 5] and all(s == "ok" for _, s, _ in tm) and any(r["system"] == "hcn" for r in rows(R2, "phase14_skala_relaxed.csv")), str(ns["_tail"]))
    st2 = status(R2)
    check("extension sessions end with G5, check and summary", done(st2, "phase15") and done(st2, "phase15x_check"))
    # (d) extension geometries present: reactions 7 and 8 get five functionals (tail 1) and enter the design of the figures and the check
    base3 = fresh_root("gpu_ext78")

    def fake_ext(ns):                       # the starting geometries of the species play the stored reaction 7 and 8 geometries (stand-in numbers only)
        for sysname in ("sn2_fcl", "sn2_ohf"):
            for sp in ns["SPECIES_KIND"][sysname]:
                ns["save_geom"](ns["geom_name"](sysname, sp, ns["LV_OWN"]), ns["start_atoms"](sysname, sp), charge=ns["species_charge"](sysname, sp))
    ns, tm = run_gpu(base3, hooks={GPU_CELL["pack"]: fake_ext}, env={"B1V5_EXTENSION": "GPU1"})
    R3 = groot(base3)
    check("EXT-GPU-1 with reaction 7 and 8 geometries: all cells ran", all(s == "ok" for _, s, _ in tm), str([t for t in tm if t[1] != "ok"]))
    r11 = rows(R3, "phase11_dft.csv")
    for sysname in ("sn2_fcl", "sn2_ohf"):
        ms = {r["method"] for r in r11 if r["system"] == sysname}
        check(f"tail 1: {sysname} has all five functionals (+ noD3)", len(ms) == 6, str(sorted(ms)))
    check("the figures and the check use eight reactions after tail 1", ns["DESIGN"]["reactions"] == ns["DESIGN_REACTIONS"] + ["sn2_fcl", "sn2_ohf"], str(ns["DESIGN"]["reactions"]))
    rep = open(os.path.join(R3, "check_report.md")).read()
    check("check report: L11 to L18 are no longer NOT RUN 'outside the design' once reactions 7 and 8 are in", "outside the 1-hour core design" not in rep.split("## All items")[1].split("L11")[1][:400] if "L11" in rep else True)


def test_gpu_nogpu():
    print("== gpu_nogpu ==")
    base = fresh_root("gpu_nogpu")
    ns, tm = run_gpu(base, standin=False)
    R = groot(base, standin=False)
    st = status(R)
    print("  ", {k: v["status"] for k, v in sorted(st.items())})
    check("all cells ran", all(s == "ok" for _, s, _ in tm), str([t for t in tm if t[1] != "ok"]))
    check("G0 ran (diagnostics)", done(st, "phase10"))
    for ph in ("phase11", "phase11p", "phase13", "phase12", "phase15"):
        check(f"{ph} skipped with a message", st.get(ph, {}).get("status") == "skipped" and "no GPU" in st[ph].get("note", ""), str(st.get(ph)))
    tim = json.load(open(os.path.join(R, "phase10_timing.json")))
    check("diagnostics recorded", tim.get("gpu_class") == "NONE" and "no GPU" in tim.get("note", ""))
    check("no energy CSV written", not any(os.path.exists(os.path.join(R, f)) for f in ("phase11_dft.csv", "phase12_basis.csv")))
    check("RUN_SUMMARY written", os.path.exists(os.path.join(R, "RUN_SUMMARY.md")))


def test_gpu_gate():
    print("== gpu_gate ==")
    base = fresh_root("gpu_gate")
    ns, tm = run_gpu(base, only=list(range(0, GPU_CELL["g0"] + 1)))
    R = groot(base)
    ns["NB4_HCN_E"]["pbe"] = [e + 5e-4 for e in ns["NB4_HCN_E"]["pbe"]]
    g = ns["g_gate"]()
    check("g_gate reports FAIL for a 5e-4 hartree difference", g["status"] == "FAIL", str(g))
    ns["GSTATE"]["gate_ok"] = False; ns["save_gstate"]()
    ns["STANDIN"] = False; ns["NO_GPU"] = False
    why = ns["blocked_1h"]("phase11")
    check("blocked_1h stops G1a after a failed gate", why and "G0 reproduction" in why, str(why))
    check("... and G1b, G3, G2", all(ns["blocked_1h"](p) for p in ("phase11p", "phase13", "phase12")))
    check("G0 itself is never blocked", ns["blocked_1h"]("phase10") is None)
    ns["STOP_ON_G0_FAIL"] = False
    check("STOP_ON_G0_FAIL = False lets them start", ns["blocked_1h"]("phase11") is None)
    ns["STANDIN"] = True; ns["STOP_ON_G0_FAIL"] = True
    # an end-to-end run with the gate failing: G1a..G2 skipped, G5/check/summary run
    base2 = fresh_root("gpu_gate2")

    def fail_gate(ns):
        ns["NB4_HCN_E"]["pbe"] = [e + 5e-4 for e in ns["NB4_HCN_E"]["pbe"]]
        orig = ns["g0_1h"]

        def g0(*a, **k):
            r = orig(*a, **k)
            ns["GSTATE"]["gate_ok"] = False
            ns["save_gstate"]()
            return r
        ns["g0_1h"] = g0
    ns, tm = run_gpu(base2, hooks={GPU_CELL["g0"]: fail_gate})
    R2 = groot(base2)
    st = status(R2)
    check("failed gate: G1a, G1b, G3, G2 skipped with the reason", all(st.get(p, {}).get("status") == "skipped" and "G0" in st[p].get("note", "") for p in ("phase11", "phase11p", "phase13", "phase12")), str({p: st.get(p) for p in ("phase11", "phase12")}))
    check("failed gate: the summary and the check still ran", os.path.exists(os.path.join(R2, "RUN_SUMMARY.md")) and done(st, "phase15x_check"), str(st.get("phase15x_check")))


def test_gpu_routes():
    print("== gpu_routes (faked GPU) ==")
    base = fresh_root("gpu_routes")
    ns, tm = run_gpu(base, only=list(range(0, GPU_CELL["pack"] + 1)))
    ns["STANDIN"] = False
    ns["GPU_CLASS"] = "T4"
    ns["gpu_available"] = lambda: True
    ns["GSTATE"].update(routes={"PBE": "gpu", "r2SCAN": "cpu", "B3LYP-D3BJ": "gpu", "wB97M-V": "none", "wB97X-V": "cpu"}, anchor="wB97X-V", skala_ok=True, skala_route="gpu", gate_ok=True)
    st = ns["resolve_methods"](["Skala-1.1", "ANCHOR", "B3LYP-D3BJ", "r2SCAN", "PBE"])
    pa = ns["path_ok_methods"](["Skala-1.1", "ANCHOR", "B3LYP-D3BJ", "r2SCAN", "PBE"])
    check("anchor falls back to wB97X-V", "wB97X-V" in st and "wB97M-V" not in st, str(st))
    check("CPU-routed functionals are stationary-only (not along paths)", "r2SCAN" in st and "r2SCAN" not in pa and "wB97X-V" not in pa, str(pa))
    check("order Skala, anchor, B3LYP, r2SCAN, PBE", st == ["Skala-1.1", "wB97X-V", "B3LYP-D3BJ", "r2SCAN", "PBE"], str(st))
    ns["GSTATE"]["skala_ok"] = False
    check("Skala dropped when unusable", "Skala-1.1" not in ns["resolve_methods"](["Skala-1.1", "ANCHOR"]))
    ns["STANDIN"] = True


def test_gpu_align():
    print("== gpu_align ==")
    base = fresh_root("gpu_align")
    ns, tm = run_gpu(base, only=list(range(0, GPU_CELL["g1a"] + 1)))
    R = groot(base)
    p = os.path.join(R, "phase11_dft.csv")
    rs = list(csv.DictReader(open(p)))
    n_df = 0
    for r in rs:
        if r["system"] == "nh3" and r["point_id"] == "TS" and r["method"] == "PBE" and r["geometry_level"] == "own_b3lyp":
            x = json.loads(r["extra"]); x["df"] = True; r["extra"] = json.dumps(x); n_df += 1
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rs[0].keys())); w.writeheader(); w.writerows(rs)
    ns["g1a_align"]()
    rs = rows(R, "phase11_dft.csv")
    al = [r for r in rs if r["system"] == "nh3" and r["point_id"].endswith("@df") and r["method"] == "PBE"]
    check("a density-fitted TS gets its partner recomputed with density fitting (NH3@df)", n_df == 1 and any(r["point_id"] == "NH3@df" for r in al), f"{[r['point_id'] for r in al]}")
    check("the recomputed partner is flagged df = true; the unfitted row is kept", all(json.loads(r["extra"]).get("df") is True for r in al) and any(r["system"] == "nh3" and r["point_id"] == "NH3" and r["method"] == "PBE" for r in rs))
    n0 = len(rs)
    ns["g1a_align"]()
    check("running the alignment again adds nothing", len(rows(R, "phase11_dft.csv")) == n0 and no_dups(rows(R, "phase11_dft.csv")))


FAKE_SKALA_CALLS = []


def install_fake_gpu(bad_plain_gradient=False, accept_with_dftd3=True):
    import test_gpu_notebook as T
    T.FAKE_SKALA_CALLS = FAKE_SKALA_CALLS
    return T.install_fake_gpu(bad_plain_gradient=bad_plain_gradient, accept_with_dftd3=accept_with_dftd3)


def remove_fake_gpu():
    import test_gpu_notebook as T
    T.remove_fake_gpu()


def test_gpu_skalad3():
    print("== gpu_skalad3 ==")
    base = fresh_root("gpu_d3")
    ns, tm = run_gpu(base, only=list(range(0, GPU_CELL["pack"] + 1)))
    ns["STANDIN"] = False
    ns["GSTATE"].update(routes={}, skala_ok=True, skala_route="cpu", gate_ok=True)
    at = ns["start_atoms"]("hcn", "HCN")
    FAKE_SKALA_CALLS.clear()
    install_fake_gpu()
    mol = ns["make_mol"](at, 0, "STO-3G")
    mf1 = ns["_build_skala"](mol, "cpu", d3=True); mf0 = ns["_build_skala"](mol, "cpu", d3=False)
    check("D3 object present for the headline, absent for the v4 variant; with_dftd3 passed explicitly", mf1.with_dftd3 is not None and mf0.with_dftd3 is None and FAKE_SKALA_CALLS == [True, False], str(FAKE_SKALA_CALLS))
    r1 = ns["energy_of"]("Skala-1.1", at, 0, "STO-3G"); r0 = ns["energy_of"]("v4-skala", at, 0, "STO-3G")
    x1, x0 = r1["extra"], r0["extra"]
    check("headline: with_dftd3 true, D3 energy and parameter set stored", x1.get("with_dftd3") is True and x1.get("d3_hartree") is not None and x1.get("d3_set") == "b3lyp5", str(x1))
    check("v4 variant: D3 off", x0.get("with_dftd3") is False and "d3_hartree" not in x0)
    check("energy minus stored D3 equals the energy without D3 (same SCF)", abs(r1["energy"] - x1["d3_hartree"] - r0["energy"]) < 1e-8)
    ns["setup_root"](base + "-x", mount=False)
    ns["q_energy_row"]("t_d3.csv", "hcn", "HCN", "own_b3lyp", "Skala-1.1", "STO-3G", at, 0, "gpu:tz:hcn")
    rs = rows(base + "-x", "t_d3.csv")
    check("q_energy_row writes the headline row and the noD3 row through the queue", sorted(r["method"] for r in rs) == ["Skala-1.1", "Skala-1.1-noD3"], str([r["method"] for r in rs]))
    remove_fake_gpu(); FAKE_SKALA_CALLS.clear(); install_fake_gpu(accept_with_dftd3=False)
    try:
        ns["energy_of"]("Skala-1.1", at, 0, "STO-3G"); loud = False
    except RuntimeError as ex:
        loud = "with_dftd3" in str(ex)
    check("a Skala without the with_dftd3 keyword fails loudly; nothing was run with a guessed setting", loud and FAKE_SKALA_CALLS == [])
    remove_fake_gpu()


def test_gpu_rescue():
    print("== gpu_rescue (Skala SCF that fails until it starts from the PBE density; RETRY_FAILED repeats a failed row) ==")
    import numpy as np
    base = fresh_root("gpu_rescue")
    ns, tm = run_gpu(base, only=list(range(0, GPU_CELL["pack"] + 1)))
    ns["STANDIN"] = False
    ns["GSTATE"].update(routes={"PBE": "cpu"}, skala_ok=True, skala_route="cpu", gate_ok=True)
    at = ns["start_atoms"]("hcn", "HCN")
    calls = []
    d3 = -0.001

    def fake_scf(method, atoms, charge, basis, route, dm0=None, df=False, level_shift=0.0, max_cycle=100, variant="plain"):
        calls.append((method, dm0 is not None, df, level_shift))
        if method == "PBE":
            return dict(energy=-93.0, converged=True, dm=np.eye(3), nbf=3, d3=None, d3xc=None)
        ok = dm0 is not None and level_shift >= 0.3
        return dict(energy=-93.1, converged=ok, dm=np.eye(3) if ok else None, nbf=3, d3=d3, d3xc="b3lyp5")
    ns["_scf_once"] = fake_scf
    ns["is_skala"] = lambda m: m in ("Skala-1.1", "v4-skala")
    r = ns["energy_of"]("Skala-1.1", at, 0, "STO-3G")
    check("the rescue attempt (PBE start + level shift) converged where the plain attempts did not", r["converged"] and "PBE start" in r["extra"]["attempt"], str(r["extra"]))
    check("the PBE density was computed first and used as the start", ("PBE", False, False, 0.0) in calls and any(c[0] == "Skala-1.1" and c[1] for c in calls), str(calls))
    ns["setup_root"](base + "-x", mount=False)
    ns["RETRY_FAILED"] = False
    ns["append_row"]("t_rs.csv", ns["make_row"]("hcn", "TS", "own_b3lyp", "Skala-1.1", "STO-3G", float("nan"), False, "", extra=dict(error="x", d3_hartree=None)))
    ns["ensure_nod3"]("t_rs.csv", "hcn", "TS", "own_b3lyp", "STO-3G")
    ns["RETRY_FAILED"] = True
    ns["q_energy_row"]("t_rs.csv", "hcn", "TS", "own_b3lyp", "Skala-1.1", "STO-3G", at, 0, "gpu:tz:hcn")
    rs = rows(base + "-x", "t_rs.csv")
    sk = [x for x in rs if x["method"] == "Skala-1.1"]; nd = [x for x in rs if x["method"] == "Skala-1.1-noD3"]
    check("RETRY_FAILED: the failed headline row is repeated (two rows, the last converged) and the noD3 row follows the converged one",
          len(sk) == 2 and sk[-1]["converged"] in ("1", "True") and nd and nd[-1]["converged"] in ("1", "True"), str([(x["method"], x["converged"]) for x in rs]))


def test_gpu_weights():
    print("== gpu_weights (install and weights caps) ==")
    base = fresh_root("gpu_w")
    ns, tm = run_gpu(base, only=list(range(0, GPU_CELL["settings"] + 1)))
    # a hanging download must stop at the cap: fake a subprocess that sleeps
    ns["STANDIN"] = False; ns["NO_GPU"] = False
    import subprocess as sp
    real = sp.run
    seen = {}

    def fake_run(cmd, *a, **k):
        seen["timeout"] = k.get("timeout")
        raise sp.TimeoutExpired(cmd, k.get("timeout") or 1)
    sp.run = fake_run
    ns["subprocess"].run = fake_run if "subprocess" in ns else sp.run
    t0 = time.time()
    try:
        out = ns["weights_setup"](0.5)
    finally:
        sp.run = real
        if "subprocess" in ns:
            ns["subprocess"].run = real
    check("a hanging weights download is stopped by the cap and reported, nothing raised", "stopped at the cap" in out.get("prefetch", "") and ns["SKALA_PREFETCH_FAILED"] is True, str(out))
    check("the cap passed to the download is 30 s for 0.5 min", seen.get("timeout") == 30.0, str(seen))
    ns["STANDIN"] = True
    ns["GSTATE"].update(skala_ok=True)
    ns["SKALA_PREFETCH_FAILED"] = False


def test_gpu_fake():
    print("== gpu_fake ==")
    src = os.environ.get("B1V5_CPU_SMOKE_ROOT")
    if not src or not os.path.isdir(src):
        print("  skipped: set B1V5_CPU_SMOKE_ROOT to the folder of a finished CPU smoke run (test_cpu_full leaves one: its path is printed)")
        return
    base = fresh_root("gpu_fake")
    shutil.rmtree(base)
    shutil.copytree(src, base)
    install_fake_gpu()
    smoke_basis = json.load(open(os.path.join(base, "config.json")))["basis_override"]

    def after_settings(ns):
        ns["BASIS_OVERRIDE"].update(smoke_basis)
        ns["NO_GPU"], ns["GPU_CLASS"], ns["GPU_NAME"] = False, "T4", "Fake T4"
        ns["SYSTEMS_TO_RUN"] = ["hcn", "nh3", "ch4"]
        ns["XC_SPEC"]["wB97M-V"] = dict(xc="wb97x")
        ns["STOP_ON_G0_FAIL"] = False
    ns, tm = run_gpu(base, standin=False, hooks={GPU_CELL["pack"]: after_settings}, only=list(range(0, GPU_CELL["g1a"] + 1)))
    tim = json.load(open(os.path.join(base, "phase10_timing.json")))
    gate = tim["gate"]
    print("  G0 gate per method:", {k: f"{v:.1e}" for k, v in gate["per_method"].items()}, "| routes", tim["routes"])
    check("stored v4 energies reproduced to 1e-5 hartree: PBE (non-stand-in path, fake GPU)", gate["per_method"].get("v4-pbe", 1) < 1e-5, str(gate["per_method"]))
    check("stored v4 energies reproduced to 1e-5 hartree: B3LYP", gate["per_method"].get("v4-b3lyp", 1) < 1e-5, str(gate["per_method"]))
    check("fake Skala (PBE) does NOT match the stored Skala: the gate fails as it must", gate["status"] == "FAIL" and tim.get("gate_ok") is False)
    check("no gradient test in the core G0", tim.get("grad_ok") is None)
    r11 = rows(base, "phase11_dft.csv")
    check("G1a ran after the gate failure because STOP_ON_G0_FAIL was off, with GPU-route rows", len(r11) > 0 and any(json.loads(r["extra"]).get("route") == "gpu" for r in r11))
    sk = [r for r in r11 if r["method"] == "Skala-1.1"]
    check("Skala rows: with_dftd3 true, D3 energy, D3 set; every construction explicit", len(sk) > 0 and all(isinstance(c, bool) for c in FAKE_SKALA_CALLS) and all(json.loads(r["extra"]).get("with_dftd3") is True and json.loads(r["extra"]).get("d3_set") == "b3lyp5" for r in sk if r["converged"] in ("1", "True")))
    remove_fake_gpu()
    shutil.rmtree(base, ignore_errors=True)


CPU_TESTS = ["cpu_full", "cpu_resume", "cpu_admit", "cpu_crash", "cpu_deadline", "cpu_ext", "cpu_final"]
GPU_TESTS = ["gpu_full", "gpu_resume", "gpu_crash", "gpu_deadline", "gpu_tail", "gpu_nogpu", "gpu_gate", "gpu_routes", "gpu_skalad3", "gpu_rescue", "gpu_weights", "gpu_align", "gpu_fake", "gpu_final"]

if __name__ == "__main__":
    prepare()
    want = sys.argv[1:] or (CPU_TESTS + GPU_TESTS)
    cpu_base = gpu_base = None
    for w in want:
        print(f"\n######## {w}", flush=True)
        if w == "cpu_full":
            cpu_base = test_cpu_full_quiet(); print("  CPU smoke folder (use it as B1V5_CPU_SMOKE_ROOT):", cpu_base)
        elif w == "cpu_resume":
            test_cpu_resume(test_cpu_full_quiet())
        elif w == "cpu_ext":
            test_cpu_ext()
        elif w == "cpu_final":
            test_final_cell(CPU_NB, "cpu")
        elif w == "gpu_full":
            gpu_base = test_gpu_full()
        elif w == "gpu_resume":
            gpu_base = gpu_base or test_gpu_full(); test_gpu_resume(gpu_base)
        elif w == "gpu_final":
            test_final_cell(GPU_NB, "gpu")
        else:
            globals()["test_" + w]()
        if cpu_base and not os.environ.get("B1V5_CPU_SMOKE_ROOT"):
            os.environ["B1V5_CPU_SMOKE_ROOT"] = cpu_base
    bad = [n for n, ok in RESULTS if not ok]
    print(f"\n{len(RESULTS) - len(bad)} of {len(RESULTS)} checks passed")
    if bad:
        print("FAILED:", *bad, sep="\n  ")
        sys.exit(1)
