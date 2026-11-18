# ==============================================================================
# b1v5_1h.py  --  the DEADLINE QUEUE and the six-reaction DESIGN of the 1-hour notebooks (PLAN_1hour.md, decisions D3, D4, D5).
# Pasted into BOTH 1-hour notebooks after b1v5_common.py / b1v5_cpu.py (and b1v5_gpu.py in the GPU notebook), in the same namespace.
# It REDEFINES a few names of b1v5_common (check_time) and adds the queue; nothing in the full-size files was changed.
#
# In plain words:
#   * ONE ordered list of calculations per notebook.  Before EVERY calculation `qrun` asks `q_admit`:
#         start only if   elapsed + 1.5 x (running mean time of this kind of calculation) < DEADLINE_MIN        (the deadline test)
#         and             phase elapsed + mean time <= the phase cap                                           (the cap test)
#     The mean is the running mean of the calculations of that kind done in THIS session; before the first one, a prior (measured in the
#     sandbox, scaled by a speed factor measured at the start) is used.  A calculation that fails the test is NOT started and is listed as
#     "not run (time)" in RUN_SUMMARY.md and queue_not_run.csv.  A phase with refused calculations is not marked DONE, so a second Run all
#     continues exactly there.
#   * After the core list comes a TAIL list (same admission test).  At DEADLINE_MIN nothing new is started; figures, check and summary
#     (FINAL_MIN) always run, whatever the clock says.
#   * Nothing here can start a calculation after the deadline; the worst overrun is the single calculation that was admitted (about 4 min).
# ==============================================================================
import os, sys, json, csv, time, contextlib, signal, threading
import numpy as np
try:
    b1v5_common_loaded = COLUMNS      # noqa
except NameError:
    from b1v5_common import *         # noqa

DEADLINE_MIN = 50.0          # no calculation is started after this many minutes (counted from the first line of the first cell)
FINAL_MIN = 5.0              # reserved for figures, check, summary
ADMIT_FACTOR = 1.5           # a calculation is admitted only if elapsed + ADMIT_FACTOR * mean < deadline
CALC_MAX_SEC_BASE = 300      # one SCF attempt may not run longer than this (the plan: no core calculation over about 4 min)
CALC_TOTAL_SEC = 420         # all retry attempts of one calculation together (GPU notebook): then the point is stored as failed and the queue goes on

# ------------------------------------------------------------------ the six-reaction design
DESIGN_REACTIONS = ["hcn", "nh3", "h2co", "ch3nc", "sn2_ff", "sn2_clcl"]
EXT_REACTIONS = ["sn2_fcl", "sn2_ohf"]
LIT_ALL = list(LIT) if "LIT" in globals() else None
if "ALL_REACTIONS" not in globals():
    ALL_REACTIONS = list(REACTIONS)          # the eight of the full design (b1v5_common)
    ALL_SYSTEM_ORDER = list(SYSTEM_ORDER)
DESIGN = dict(reactions=list(DESIGN_REACTIONS))


def set_design(reactions):
    """Restrict the notebook to these reactions (plus the CH4 stretch).  Rebinds REACTIONS, SYSTEM_ORDER, TIER_B and LIT in the notebook namespace,
    so that the check, the figures and the tables of the full-size modules see only this design."""
    g = globals()
    reactions = [r for r in ALL_REACTIONS if r in set(reactions)]
    DESIGN["reactions"] = list(reactions)
    g["REACTIONS"] = list(reactions)
    g["SYSTEM_ORDER"] = list(reactions) + ["ch4"]
    g["TIER_B"] = [k for k in reactions if SYSTEMS[k]["tier"] == "B"]
    if "LIT_ALL" in g and g["LIT_ALL"] is None:
        g["LIT_ALL"] = list(g["LIT"])
    if g.get("LIT_ALL"):
        g["LIT"][:] = [l for l in g["LIT_ALL"] if l["system"] in reactions]
    return list(reactions)


def lit_out_of_design():
    return [l for l in (LIT_ALL or []) if l["system"] not in DESIGN["reactions"]]


# ------------------------------------------------------------------ the queue
Q = dict(t0=None, deadline_min=DEADLINE_MIN, final_min=FINAL_MIN, phase=None, cap_t=None, final=False, means={}, counts={}, prior={}, scale=1.0,
         refused=[], caps={}, started=[], tail_names=[], log=[])
PRIOR_DEFAULT_S = 20.0


class CalcBudget(Exception):
    """One calculation (all its retry attempts together) used more than CALC_TOTAL_SEC: the point is stored as failed and the queue goes on."""


class QTimeout(PhaseTimeout):
    """The hard stop (deadline + final minutes) arrived inside one calculation."""


def q_start(deadline_min=None, t0=None, final_min=None):
    """Call once at the top of the notebook; t0 = time.time() taken in the very first line of the first cell (install time counts)."""
    if deadline_min is not None:
        Q["deadline_min"] = float(deadline_min)
    if final_min is not None:
        Q["final_min"] = float(final_min)
    Q["t0"] = t0 if t0 else time.time()
    _CLOCK["t0"] = Q["t0"]
    _CLOCK["total_deadline"] = None            # the queue has its own hard stop (check_time below)
    globals()["MAX_TOTAL_HOURS"] = round((Q["deadline_min"] + Q["final_min"]) / 60.0, 2)      # only the number printed in RUN_SUMMARY.md
    return Q["t0"]


def q_now_min():
    return (time.time() - Q["t0"]) / 60.0 if Q["t0"] else 0.0


def q_remaining_min():
    return Q["deadline_min"] - q_now_min()


def q_hard_stop_t():
    return Q["t0"] + 60.0 * (Q["deadline_min"] + Q["final_min"])


def check_time():
    """REPLACES b1v5_common.check_time (same name, so run_point and the optimisers pick it up): stops the phase at its limit and everything at the
    hard stop (deadline + final minutes), except in the final phases (figures, check, summary), which must always run."""
    now_ = time.time()
    if _CLOCK.get("deadline") and now_ > _CLOCK["deadline"]:
        raise PhaseTimeout(f"phase time limit reached in {_CLOCK['phase']}")
    if Q["t0"] and not Q["final"] and now_ > q_hard_stop_t():
        raise PhaseTimeout(f"hard stop ({Q['deadline_min'] + Q['final_min']:.0f} min) reached")
    if Q.get("calc_end") and now_ > Q["calc_end"]:
        raise CalcBudget(f"one calculation used more than {CALC_TOTAL_SEC} s with all its retries")


@contextlib.contextmanager
def phase(name, title="", max_min=None, force=False):
    """REPLACES b1v5_common.phase.  Same behaviour, plus one fix found by the 1-hour tests: a hard stop of the cell (KeyboardInterrupt, e.g. the Colab
    'interrupt execution' button, or SystemExit) is not an Exception, so the original marked the interrupted phase DONE and a rerun skipped its remaining
    points.  Here it is recorded as 'interrupted', gets no DONE marker, and the interrupt is passed on."""
    class _Ph:
        skip = False
    ph = _Ph()
    print(f"\n===== {name} {title} =====  (notebook time so far {elapsed_hours():.2f} h)")
    if is_done(name) and not force:
        print(f"  {name} already DONE (marker found) -> skipped. Delete DONE_{name} to run it again.")
        _set_status(name, status="already done", title=title)
        ph.skip = True
        yield ph
        return
    t0 = time.time()
    _CLOCK["phase"] = name
    _CLOCK["deadline"] = t0 + 60.0 * (max_min or DEFAULT_PHASE_MIN)
    status, note, hard = "done", "", None
    try:
        check_time()
    except PhaseTimeout as ex:
        ph.skip = True; status, note = "stopped (time limit)", str(ex)
        print("  ** not started:", ex)
    try:
        yield ph
    except PhaseTimeout as ex:
        status, note = "stopped (time limit)", str(ex)
        print("  ** phase stopped cleanly:", ex, "- rerun the notebook to continue, finished points are kept.")
    except Exception as ex:
        status, note = "FAILED", repr(ex)[:200]
        log_error(name, "(whole phase)", ex)
    except BaseException as ex:                      # KeyboardInterrupt, SystemExit: not DONE, and the interrupt goes on
        status, note, hard = "interrupted", repr(ex)[:200], ex
    finally:
        _CLOCK["deadline"] = None
        secs = time.time() - t0
        n_err = 0
        try:
            n_err = sum(1 for r in csv.DictReader(open(path(ERRORS_CSV))) if r["phase"] == name) if os.path.exists(path(ERRORS_CSV)) else 0
        except Exception:
            pass
        if status == "done" and n_err:
            status = "done with errors"
        if status == "done":
            mark_done(name, title)
        log_phase(name, secs, status)
        _set_status(name, status=status, title=title, seconds=round(secs, 1), errors=n_err, note=note, finished=now())
        print(f"  {name}: {status}, {secs / 60:.1f} min")
    if hard is not None:
        raise hard


def q_est(kind):
    """Expected seconds of one calculation of this kind: running mean of this session, else the prior (scaled by the measured speed factor)."""
    n = Q["counts"].get(kind, 0)
    if n:
        return Q["means"][kind]
    k = kind
    while True:
        if k in Q["prior"]:
            return Q["prior"][k] * Q["scale"]
        if ":" not in k:
            break
        k = k.rsplit(":", 1)[0]
    return PRIOR_DEFAULT_S * Q["scale"]


def q_record(kind, seconds):
    n = Q["counts"].get(kind, 0)
    Q["means"][kind] = (Q["means"].get(kind, 0.0) * n + float(seconds)) / (n + 1)
    Q["counts"][kind] = n + 1


def q_admit(kind):
    """(True, '') or (False, reason).  The deadline test and the cap test of the module header."""
    est = q_est(kind)
    now_ = time.time()
    if Q["t0"] and not Q["final"] and now_ + ADMIT_FACTOR * est > Q["t0"] + 60.0 * Q["deadline_min"]:
        return False, f"deadline: {q_now_min():.1f} min used, {ADMIT_FACTOR} x {est:.0f} s would pass minute {Q['deadline_min']:.0f}"
    if Q["cap_t"] and now_ + est > Q["cap_t"]:
        return False, f"phase cap of {Q['caps'].get(Q['phase'], '?')} min ({Q['phase']}): {est:.0f} s would pass it"
    return True, ""


@contextlib.contextmanager
def q_time_limit(sec):
    """Interrupt ONE calculation at the hard stop (where the operating system allows it: main thread, SIGALRM)."""
    if not sec or threading.current_thread() is not threading.main_thread() or not hasattr(signal, "setitimer"):
        yield
        return

    def _h(signum, frame):
        raise QTimeout("hard stop reached inside a calculation")
    old = signal.signal(signal.SIGALRM, _h)
    signal.setitimer(signal.ITIMER_REAL, float(sec))
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)


def _limited(fn):
    def g():
        remain = q_hard_stop_t() - time.time() - 20.0
        if "CALC_MAX_SEC" in globals():                     # GPU notebook: the per-attempt limit of b1v5_gpu follows the hard stop
            globals()["CALC_MAX_SEC"] = int(max(20, min(CALC_MAX_SEC_BASE, remain)))
        Q["calc_end"] = time.time() + CALC_TOTAL_SEC if "CALC_MAX_SEC" in globals() else None
        try:
            with q_time_limit(max(20.0, remain)):
                return fn()
        finally:
            Q["calc_end"] = None
    return g


def q_refuse(label, kind, reason):
    Q["refused"].append(dict(phase=Q["phase"], label=label, kind=kind, reason=reason, est_s=round(q_est(kind), 1), minute=round(q_now_min(), 1)))


def qrun(csv_name, system, pid, level, method, basis, fn, kind, coordinate="", retry_failed=False, label=None, quiet=True):
    """run_point with the admission test in front.  Returns the new row, or None (already done, or refused: see Q['refused'])."""
    if have(csv_name, system, pid, level, method, basis, retry_failed):
        return None
    ok, why = q_admit(kind)
    if not ok:
        q_refuse(label or f"{system}:{pid}:{method}:{basis}", kind, why)
        return None
    t0 = time.time()
    row = run_point(csv_name, system, pid, level, method, basis, _limited(fn), coordinate=coordinate, retry_failed=retry_failed, quiet=quiet)
    if row is not None and str(row["converged"]) in ("1", "True"):                    # run_point returns converged as an int, the CSV has text
        q_record(kind, time.time() - t0)
    return row


@contextlib.contextmanager
def qphase(name, title, cap_min, final=False, force=None):
    """A phase of the queue: `cap_min` is its cap.  A phase with refused calculations is left 'stopped (time limit)', never DONE, so that the next
    Run all continues.  final=True: figures / check / summary phases, exempt from the hard stop."""
    if force is None:
        force = bool(globals().get("FORCE_RERUN", False))               # the FORCE_RERUN switch of the settings cell
    Q["phase"], Q["cap_t"], Q["final"] = name, time.time() + 60.0 * cap_min, bool(final)
    Q["caps"][name] = cap_min
    Q["started"].append((name, round(q_now_min(), 1)))
    n0 = len(Q["refused"])
    try:
        with phase(name, title, max_min=cap_min + 2.0, force=force) as ph:
            yield ph
            k = len(Q["refused"]) - n0
            if k and not ph.skip:
                raise PhaseTimeout(f"{k} calculation(s) not run (time): a second Run all continues here")
    finally:
        Q["cap_t"], Q["final"] = None, False
        if len(Q["refused"]) > n0:
            print(f"  ** {len(Q['refused']) - n0} calculation(s) NOT RUN (time) in {name}; first: {Q['refused'][n0]['label']} ({Q['refused'][n0]['reason'][:80]})")


def q_save_log():
    """queue_log.json (means, counts, refused of this session) and queue_not_run.csv (appended)."""
    try:
        old = read_json("queue_log.json", {}) or {}
        sessions = old.get("sessions", [])
        sessions.append(dict(at=now(), elapsed_min=round(q_now_min(), 1), deadline_min=Q["deadline_min"], means={k: round(v, 2) for k, v in Q["means"].items()},
                             counts=dict(Q["counts"]), n_refused=len(Q["refused"]), phases=Q["started"], hardware=hw()))
        write_json("queue_log.json", dict(sessions=sessions))
        p = path("queue_not_run.csv")
        new = not os.path.exists(p)
        with open(p, "a", newline="") as fh:
            w = csv.writer(fh)
            if new:
                w.writerow(["time", "phase", "label", "kind", "estimate_s", "minute", "reason"])
            for r in Q["refused"]:
                w.writerow([now(), r["phase"], r["label"], r["kind"], r["est_s"], r["minute"], r["reason"]])
    except Exception as ex:
        print("could not write the queue log:", repr(ex)[:100])


def queue_lines():
    """Lines for RUN_SUMMARY.md: minutes against caps, what was not run (time), measured mean seconds per kind."""
    st = read_json("phase_status.json", {}) or {}
    L = ["", "## Deadline queue", "",
         f"Deadline: no calculation is started after minute {Q['deadline_min']:.0f}; {Q['final_min']:.0f} more minutes are reserved for figures, check and summary. "
         f"Minutes used at this point: {q_now_min():.1f}. Admission test: elapsed + {ADMIT_FACTOR} x running mean < deadline; a phase also stops at its cap.", "",
         "| phase | minutes used (this and earlier sessions) | cap (min) | status |", "|---|---|---|---|"]
    for k in sorted(st):
        if k in Q["caps"]:
            L.append(f"| {k} | {st[k].get('seconds', 0) / 60:.1f} | {Q['caps'][k]} | {st[k].get('status', '')} |")
    n = len(Q["refused"])
    L += ["", f"**Not run (time) in this session: {n} calculation(s).** Details in queue_not_run.csv; a second Run all continues with them."]
    if n:
        by = {}
        for r in Q["refused"]:
            by.setdefault(r["phase"], []).append(r)
        for ph, rs in by.items():
            L.append(f"- {ph}: {len(rs)} not run, first {rs[0]['label']} ({rs[0]['reason'][:90]}); last {rs[-1]['label']}")
    else:
        L.append("- none: every admitted calculation of the lists ran.")
    if Q["means"]:
        L += ["", "Measured mean seconds per calculation (this session): " + ", ".join(f"{k} {v:.1f} s (n={Q['counts'][k]})" for k, v in sorted(Q["means"].items()))]
    L += ["", f"Speed factor applied to the sandbox priors: {Q['scale']:.2f}"]
    return L


# ------------------------------------------------------------------ the geometry pack in a notebook
def pack_setup(b64_text, expected_sha256, systems=None, allow_mismatch=False, build_note=""):
    """Decode the embedded geometry pack, verify its SHA-256 and every structure's SHA-256, install it into geoms/ and the phase CSV files.
    Returns dict(ok, sha256, n_verified, problems, installed) and writes geompack_installed.json."""
    out = dict(ok=False, expected_sha256=expected_sha256, problems=[])
    try:
        pack, sha, match = load_pack_text(b64_text, expected_sha256)
        out["sha256"] = sha
        n_ok, bad = verify_hashes(pack)
        out.update(n_verified=n_ok, problems=bad, format=pack.get("format"), design=pack.get("design"), builder=pack.get("builder"))
        if not match:
            out["problems"].insert(0, f"pack SHA-256 {sha} differs from the expected {expected_sha256}")
        out["ok"] = bool(match and not bad)
        if out["ok"] or allow_mismatch:
            out["installed"] = install_pack(pack, systems, sha=sha)
        out["systems_in_pack"] = pack.get("systems")
        out["has"] = dict(ext78=all(s in pack.get("systems", []) for s in EXT_REACTIONS))
        out["_pack"] = pack
    except Exception as ex:
        out["problems"].append("could not read the pack: " + repr(ex)[:150])
    meta = {k: v for k, v in out.items() if k != "_pack"}
    write_json("geompack_installed.json", dict(meta, at=now()))
    return out


# ------------------------------------------------------------------ pip with a time cap (cell 1 of the notebooks)
def pip_install_capped(specs, timeout_min):
    """pip install with a wall-clock cap; returns (ok, seconds).  Never raises."""
    import subprocess
    t0 = time.time()
    try:
        r = subprocess.run([sys.executable, "-m", "pip", "install", "-q"] + list(specs), capture_output=True, text=True, timeout=60 * timeout_min)
        print(("ok   " if r.returncode == 0 else "FAIL ") + " ".join(specs) + f"  ({time.time() - t0:.0f} s)")
        if r.returncode:
            print(r.stderr[-400:])
        return r.returncode == 0, time.time() - t0
    except Exception as ex:
        print("pip step stopped:", repr(ex)[:120], f"({time.time() - t0:.0f} s)")
        return False, time.time() - t0
