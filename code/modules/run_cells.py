"""Helper for the tests: executes the code cells of a notebook one after the other in ONE namespace (as Colab does),
so that a test can change a function between two cells (for example to inject a crash)."""
import json, sys, time, traceback, io, contextlib


def load_cells(nb_path):
    nb = json.load(open(nb_path))
    return [("".join(c["source"]), i) for i, c in enumerate(nb["cells"]) if c["cell_type"] == "code"]


def run_cells(nb_path, hooks=None, only=None, ns=None, quiet=True, log=None):
    """hooks: {code-cell number: function(ns)} called just before that cell.  Returns (ns, timings)."""
    ns = ns if ns is not None else {"__name__": "__main__"}
    cells, timings = load_cells(nb_path), []
    for k, (src, idx) in enumerate(cells):
        if only is not None and k not in only:
            continue
        if hooks and k in hooks:
            hooks[k](ns)
        t0 = time.time()
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf) if quiet else contextlib.nullcontext():
                exec(compile(src, f"<cell {k}>", "exec"), ns)
            status = "ok"
        except BaseException as ex:
            status = "EXCEPTION " + repr(ex)[:200]
            buf.write("\n" + traceback.format_exc())
        out = buf.getvalue()
        if log:
            open(log, "a").write(f"\n######## cell {k} ({status}, {time.time() - t0:.1f} s)\n{out}")
        timings.append((k, status, round(time.time() - t0, 1)))
        print(f"cell {k}: {status} {time.time() - t0:.1f} s", flush=True)
        if status != "ok" and not (hooks and hooks.get("continue_on_error")):
            pass
    return ns, timings
