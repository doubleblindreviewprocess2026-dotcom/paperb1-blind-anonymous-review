"""Tests of the 1-hour check and frozen criteria (b1v5_1h.py + b1v5_1h_check.py) on SYNTHETIC data (nothing here is a chemistry result).
    python test_1h_check.py        (needs numpy, pandas, matplotlib, openpyxl, pyscf not needed)
The modules are executed in ONE namespace in the order of the notebook cells (as Colab does), then the six-reaction design is set."""
import os, sys, tempfile, shutil, json, types
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault("MPLBACKEND", "Agg")
import check_test_synthetic as T          # the synthetic data generator of the full-size test (eight reactions)

ORDER = ["b1v5_common.py", "b1v5_cpu.py", "b1v5_v4_stored.py", "b1v5_lit_table.py", "b1v5_check.py", "b1v5_gpu_figs.py", "b1v5_1h.py", "b1v5_1h_check.py"]
OK = []


def chk(cond, text):
    OK.append(bool(cond))
    print(("  ok  " if cond else "  BAD ") + text)
    return bool(cond)


def namespace():
    ns = {"__name__": "__main__"}
    for f in ORDER:
        exec(compile(open(os.path.join(HERE, f)).read(), f, "exec"), ns)
    return ns


def items(rep, item, scope_start):
    return [x for x in rep["items"] if x["item"] == item and x["scope"].startswith(scope_start)]


def main():
    tmp = tempfile.mkdtemp(prefix="b1v5_1h_synth_")
    T.build(tmp, "A")                                   # writes into the module state of b1v5_common imported by T (separate namespace): write files only
    ns = namespace()
    ns["setup_root"](tmp, mount=False)
    ns["set_design"](ns["DESIGN_REACTIONS"])
    chk(ns["DESIGN"]["reactions"] == ["hcn", "nh3", "h2co", "ch3nc", "sn2_ff", "sn2_clcl"], "design = the six reactions")
    chk(ns["REACTIONS"] == ns["DESIGN"]["reactions"] and ns["SYSTEM_ORDER"][-1] == "ch4", "REACTIONS and SYSTEM_ORDER rebound")
    chk(all(l["system"] in ns["REACTIONS"] for l in ns["LIT"]) and len(ns["LIT"]) < len(ns["LIT_ALL"]), f"LIT restricted: {len(ns['LIT'])} of {len(ns['LIT_ALL'])} rows")
    rep = ns["run_check"](tmp)
    print("\n##### assertions")
    for lid in ("L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18"):
        x = [i for i in rep["items"] if i["item"] in ("R1", "R2") and i["scope"].startswith(lid + " ")]
        chk(len(x) == 2 and all(i["status"] == "NOT RUN" for i in x), f"{lid}: R1 and R2 NOT RUN (outside the design)")
    chk(not any(i["status"] == "FAIL" and i["scope"][:3] in ("L11", "L12", "L13", "L14", "L15", "L16", "L17", "L18") for i in rep["items"]), "no FAIL for any out-of-design literature row")
    for s in ("sn2_fcl", "sn2_ohf"):
        for it in ("I1", "I4", "I5", "I6", "I8", "I12"):
            x = items(rep, it, s)
            chk(x and all(i["status"] == "NOT RUN" for i in x), f"{it} {s}: NOT RUN")
    i11 = items(rep, "I11", "own vs")
    chk(i11 and i11[0]["status"] == "NOT RUN" and "D2" in i11[0]["message"], "I11 own vs database: NOT RUN (decision D2)")
    chk(items(rep, "I11", "Skala//Skala") and items(rep, "I11", "Skala//Skala")[0]["status"] == "NOT RUN", "I11 Skala//Skala: NOT RUN without G4")
    for ex_ in [("R1", "L1 ", "PASS"), ("R1", "L2 ", "PASS"), ("R2", "L1 ", "PASS"), ("R2", "L7 ", "FAIL"), ("I1", "hcn", "PASS"), ("I2", "ch4", "FAIL"), ("I5", "nh3", "FAIL"), ("I9", "v4 energies", "PASS")]:
        got = [i["status"] for i in items(rep, ex_[0], ex_[1])]
        chk(got == [ex_[2]], f"{ex_[0]} {ex_[1]}: expected {ex_[2]}, got {got}")
    p1 = items(rep, "P1", "geometry pack")
    chk(p1 and p1[0]["status"] == "NOT RUN", "P1: NOT RUN when no pack was installed (never FAIL)")
    p2 = items(rep, "P2", "gradient test")
    chk(p2 and p2[0]["status"] == "NOT RUN", "P2: NOT RUN when no gradient rows exist")
    dec = {k: v[0] for k, v in rep["decisions"].items()}
    for k, v in dict(hcn="MAIN", nh3="MAIN", h2co="MAIN", ch3nc="SUPPLEMENT", sn2_ff="SUPPLEMENT", sn2_clcl="MAIN").items():
        chk(dec[k] == v, f"decision {k}: expected {v}, got {dec[k]}")
    chk(dec["sn2_fcl"].startswith("NOT RUN") and dec["sn2_ohf"].startswith("NOT RUN"), f"extension reactions: {dec['sn2_fcl']}")
    chk(rep["paper"].startswith("4 of 6 design reactions"), "paper level counts n of 6: " + rep["paper"])
    chk(rep["counts"]["FAIL"] > 0 and os.path.exists(os.path.join(tmp, "check_report.md")), "report written")
    md = open(os.path.join(tmp, "check_report.md")).read()
    chk("1-hour design" in md and "NOT RUN" in md, "check_report.md says 1-hour design and lists NOT RUN")
    # six of six: pack and gradient rows present, everything verified -> FULL PAPER as defined
    for s in ("ch3nc", "sn2_ff"):
        pass
    d6 = {k: ("MAIN", "") for k in ns["DESIGN"]["reactions"]}
    d6.update({k: ("NOT RUN (extension)", "") for k in ("sn2_fcl", "sn2_ohf")})
    chk(ns["paper_level"](d6).startswith("6 of 6 design reactions") and "FULL PAPER" in ns["paper_level"](d6), "six MAIN reactions in at least 3 classes: " + ns["paper_level"](d6)[:80])
    d5 = dict(d6); d5["ch3nc"] = ("SUPPLEMENT", "")
    chk(ns["paper_level"](d5).startswith("5 of 6") and "short communication" in ns["paper_level"](d5), "five of six: add extras or short communication")
    d3 = dict(d6); d3.update(ch3nc=("SUPPLEMENT", ""), sn2_ff=("SUPPLEMENT", ""), sn2_clcl=("INCOMPLETE", ""))
    chk("DO NOT SUBMIT" in ns["paper_level"](d3) and "PROVISIONAL" in ns["paper_level"](d3), "three MAIN, one incomplete: provisional, do not submit")
    # frozen criteria
    F = ns["FROZEN"]
    chk("5 of 6" in F["F1"]["finding"] + F["F1"]["supports"] + F["F5"]["supports"] and "6 of 8" not in json.dumps(F), "FROZEN: F1 and F5 say 5 of 6, no '6 of 8' left")
    chk("5 of 6" in F["F4"]["supports"] and ns["F3_MIN_REACTIONS"] == 3, "FROZEN: F4 5 of 6, F3 needs 3 reactions")
    chk(ns["FROZEN_DATE"].startswith("2026-11-15"), "frozen date kept")
    # findings_outcome on the synthetic data: missing data -> NOT RUN (data missing), never an exception
    D = ns["Data"]()
    fo = ns["findings_outcome"](D)
    chk([r["finding"] for r in fo] == ["F1", "F2", "F3", "F4", "F5", "F6"] and all(r["verdict"] for r in fo), "findings_outcome: F1 to F6 with verdicts: " + ", ".join(f"{r['finding']} {r['verdict']}(n={r['n']})" for r in fo))
    chk(all(r["frozen_on"].startswith("2026-11-15") for r in fo), "frozen_on written")
    # empty folder
    tmp2 = tempfile.mkdtemp(prefix="b1v5_1h_empty_")
    ns["setup_root"](tmp2, mount=False)
    repC = ns["run_check"](tmp2)
    chk(repC["counts"]["FAIL"] == 0, f"empty folder: no FAIL, counts {repC['counts']}")
    fo2 = ns["findings_outcome"](ns["Data"]())
    chk(all(r["verdict"].startswith("NOT RUN") for r in fo2), "empty folder: every finding NOT RUN (data missing): " + str([r["verdict"] for r in fo2]))
    # gradient rows are not read as energies
    ns["setup_root"](tmp, mount=False)
    ns["append_row"]("phase00_verify.csv", ns["make_row"]("hcn", "HCN", ns["LV_OWN"], "DF-B3LYP-D3BJ-grad", "def2-TZVP", -93.0, True, "", "", "", 1.0, dict(rms=1e-5, tol=3e-4)))
    rows, _ = ns["load_rows"]()
    chk(not any(str(r["method"]).endswith("-grad") for r in rows), "load_rows drops the gradient-test rows")
    for t in (tmp, tmp2):
        shutil.rmtree(t, ignore_errors=True)
    print(f"\nRESULT: {sum(OK)} of {len(OK)} assertions passed")
    return 0 if all(OK) else 1


if __name__ == "__main__":
    sys.exit(main())
