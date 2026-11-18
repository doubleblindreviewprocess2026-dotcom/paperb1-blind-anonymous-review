"""HTML pages of the site (used by build_site.py). Every number comes from the CSV or JSON files in data/."""
import csv, glob, hashlib, json, os, re
import data_dictionary as dd
from sitelib import esc, fmt, sgn, sign_class, human_size, png_size, table, page, colab_buttons, MINUS
from build_data import RX_ORDER, RX_LABEL, M_ORDER, M_LABEL, rd, fl

HERE = os.path.dirname(os.path.abspath(__file__))
SITE = os.path.dirname(HERE)
D = os.path.join(SITE, "data")

SW = {"Skala-1.1": "--c-skala", "Skala-1.1-noD3": "--c-skala", "wB97M-V": "--c-wb97", "B3LYP-D3BJ": "--c-b3lyp", "r2SCAN": "--c-r2scan", "PBE": "--c-pbe",
      "CCSD(T)-fc": "--c-cc", "CCSD(T)+dMP2": "--c-cc", "GFN2-xTB": "--c-xtb"}


def sw(m):
    style = "background:transparent;border:2px solid var(%s)" % SW[m] if m == "Skala-1.1-noD3" else "background:var(%s)" % SW[m]
    return '<span class="sw" style="%s"></span>' % style


def mhead(m, short=False):
    return sw(m) + esc(M_LABEL[m])


def errtd(v, d=2):
    if v is None:
        return '<td class="num">–</td>'
    return '<td class="num %s" title="%s">%s</td>' % (sign_class(v), esc("%.4f" % v), sgn(v, d))


def pagelink_cards(pre, items):
    out = ['<div class="grid">']
    for slug, title, text in items:
        out.append('<a class="card" href="%s%s/"><h3>%s</h3><p>%s</p></a>' % (pre, slug, esc(title), esc(text)))
    out.append("</div>")
    return "\n".join(out)


def data_link(rel, pre):
    return '<a href="%s%s">%s</a>' % (pre, esc(rel), esc(rel))


# ---------------------------------------------------------------- index
def index_page(R, anonymous):
    rx = {r["key"]: r for r in R["reactions"]}
    rows = {(r["r"], r["m"]): r for r in R["rows"] if r["b"] == "main"}
    pub = {p["m"]: p for p in R["published"] if p["b"] == "main"}
    head = [k for k in RX_ORDER if rx[k]["headline"]]
    sk = {k: rows[(k, "Skala-1.1")]["err"] for k in head}
    pos = [v for v in sk.values() if v > 0]
    neg = [k for k, v in sk.items() if v < 0]
    for m in ("B3LYP-D3BJ", "r2SCAN", "PBE"):
        assert all(rows[(k, m)]["err"] < 0 for k in head), m
    assert len(neg) == 1 and neg[0] == "sn2_clcl"
    words = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five"}
    s1 = ("This site holds the data, code and interactive views behind a study of the learned density functional Skala-1.1 on the barriers and reaction paths of "
          "five small reactions: HCN → HNC, NH₃ inversion, H₂CO → H₂ + CO, F⁻ + CH₃F and Cl⁻ + CH₃Cl.")
    s2 = ("Against literature barriers (n = %d reactions, B3LYP-D3(BJ) geometries, def2-TZVP or def2-TZVPD), Skala-1.1 has a mean absolute error of %s kcal/mol, "
          "against %s for ωB97M-V, %s for B3LYP-D3(BJ), %s for r²SCAN and %s for PBE." % (
              pub["Skala-1.1"]["n"], fmt(pub["Skala-1.1"]["mae"]), fmt(pub["wB97M-V"]["mae"]), fmt(pub["B3LYP-D3BJ"]["mae"]), fmt(pub["r2SCAN"]["mae"]), fmt(pub["PBE"]["mae"])))
    s3 = ("The signed error of Skala-1.1 is positive for %s reactions (%s to %s kcal/mol) and negative for Cl⁻ + CH₃Cl (%s), "
          "while PBE, r²SCAN and B3LYP-D3(BJ) put all five barriers too low." % (words[len(pos)], sgn(min(pos)), sgn(max(pos)), sgn(sk["sn2_clcl"])))
    s4 = "Five reactions are a small sample: the numbers describe these reactions only, and the limitations are listed on the About page."
    cards = ['<div class="grid g5">']
    for k in head:
        r = rx[k]
        cards.append('<div class="card"><h3>%s</h3><p><span class="tag">%s</span></p><p>Reference %s kcal/mol (%s)</p><p>Skala-1.1: <b class="stat">%s</b> &nbsp; ωB97M-V: <b class="stat">%s</b></p></div>' % (
            esc(r["label"]), esc(r["cls"]), fmt(r["ref"]), esc(r["ref_source"]), sgn(rows[(k, "Skala-1.1")]["err"]), sgn(rows[(k, "wB97M-V")]["err"])))
    cards.append("</div>")
    cols = ["Skala-1.1", "Skala-1.1-noD3", "wB97M-V", "B3LYP-D3BJ", "r2SCAN", "PBE", "CCSD(T)+dMP2"]
    trs = []
    for k in head:
        r = rx[k]
        ref = "%s <span class=\"small\">%s</span>" % (fmt(r["ref"]), esc(r["ref_source"]))
        if k == "h2co":
            ref += " †"
        cells = ['<td>%s</td>' % esc(r["label"]), '<td class="num">%s</td>' % ref.replace('<span class="small">', '<br><span class="small">')]
        for m in cols:
            row = rows[(k, m)]
            cells.append(errtd(None if row["wd"] else row["err"]))
        trs.append("<tr>%s</tr>" % "".join(cells))
    for lab, key, dd_ in (("Mean absolute error", "mae", 2), ("Mean signed error", "mse", 2), ("Largest absolute error", "max", 2)):
        cells = ['<td colspan="2"><b>%s</b></td>' % lab]
        for m in cols:
            p = pub.get(m)
            cells.append('<td class="num" title="%s">%s</td>' % (esc("%.3f" % p[key]), (fmt(p[key]) if key != "mse" else sgn(p[key]))) if p else '<td class="num">–</td>')
        trs.append('<tr class="total">%s</tr>' % "".join(cells))
    trs.append('<tr class="sub"><td colspan="2">n (reactions)</td>%s</tr>' % "".join('<td class="num">%d</td>' % pub[m]["n"] for m in cols))
    short = {"Skala-1.1-noD3": "Skala, no D3", "CCSD(T)+dMP2": "CCSD(T)+\u0394MP2"}
    hdr = "".join('<th scope="col" class="num">%s%s</th>' % (sw(m), esc(short.get(m, M_LABEL[m]))) for m in cols)
    tbl = ('<div class="tablewrap" tabindex="0" role="region" aria-label="Signed barrier errors per reaction and method"><table class="compact">'
           '<caption>Signed error of the forward barrier (computed minus literature reference), kcal/mol</caption>'
           '<thead><tr><th scope="col">Reaction</th><th scope="col" class="num">Reference</th>%s</tr></thead><tbody>%s</tbody></table></div>' % (hdr, "\n".join(trs)))
    zpe = rx["h2co"]["zpe_added"]
    foot = ('<p class="small">† The H₂CO reference includes zero-point energy; a B3LYP harmonic zero-point difference of %s kcal/mol was added to every computed H₂CO barrier. '
            'The CCSD(T)+ΔMP2 column has n = %d: there is no usable coupled-cluster barrier for H₂CO. The SN2 barriers are measured from the ion-molecule complex. '
            'Colour is a second cue only: orange fill = computed barrier above the reference, blue fill = below. Source: <a href="data/csv/study/master_errors.csv">data/csv/study/master_errors.csv</a> and '
            '<a href="data/csv/study/summary_stats.csv">summary_stats.csv</a>.</p>' % (sgn(zpe, 2), pub["CCSD(T)+dMP2"]["n"]))
    nlink = ('<div class="note"><p>CH₃NC → CH₃CN was also computed. Its literature value (38.25 kcal/mol) could not be checked in the source, so it enters no statistic; its numbers are in the '
             '<a href="results/">results table</a> and can be added to the statistics there by hand.</p></div>')
    body = """<h1>Skala-1.1 on reaction barriers and reaction paths</h1>
<p class="lead">%s %s</p>
<p class="lead">%s %s</p>
<div class="btnrow"><a class="btn" href="results/">Results table</a><a class="btn alt" href="demos/">Interactive demos</a><a class="btn alt" href="figures/">Figures</a><a class="btn alt" href="data/">Data files</a><a class="btn alt" href="code/">Notebooks and code</a></div>
<h2 id="reactions">The five reactions</h2>
%s
%s
<h2 id="errors">Signed barrier errors</h2>
%s
%s
<div class="panel"><h3 style="margin-top:0">The same numbers as a strip chart</h3>
<div class="legend" id="legend"></div>
<div class="chart" id="strip" data-src="data/json/results.json"></div>
<p class="small">Each row is one reaction. Marker position = signed error of the forward barrier (kcal/mol), computed minus reference; the grey band marks ±1 kcal/mol. Circle: Skala-1.1; ring: Skala-1.1 without D3; square: ωB97M-V; triangle: B3LYP-D3(BJ); diamond: r²SCAN; cross: PBE; plus: CCSD(T)+ΔMP2.</p></div>
<h2 id="where">Where to go</h2>
%s
<h2>How to read the numbers</h2>
<ul>
<li>Every error on this site is <b>computed minus reference</b>, in kcal/mol. A positive error means the computed barrier is too high.</li>
<li>All energies are single points at B3LYP-D3(BJ) geometries optimised for this study (def2-TZVP; def2-TZVPD for the two SN2 reactions). The exceptions are labelled.</li>
<li>Every number comes from a CSV file in <a href="data/">data</a>. The interactive pages read JSON files that are built from those CSV files by <code>tools/build_site.py</code>.</li>
</ul>
""" % (esc(s1), esc(s2), esc(s3), esc(s4), "\n".join(cards), nlink, tbl, foot,
       pagelink_cards("", [("results", "Results", "Sortable table of every reaction and method, reference choices, and a mean-error calculator that reproduces the published statistics."),
                           ("demos", "Demos", "Reaction paths, the SN2 symmetric stretch, the HCN T1 map, cost against accuracy, basis-set and D3 effects."),
                           ("figures", "Figures", "All figures with captions, as PNG, PDF and SVG, with the data behind each."),
                           ("data", "Data", "Every CSV with a data dictionary, geometries as xyz, checksums."),
                           ("code", "Code", "Colab notebooks, Python modules, how to run, how to re-check the numbers."),
                           ("checks", "Checks", "The automatic check of the run: what passed, what warned, what failed, and why."),
                           ("about", "About", "How to cite, licences, software versions, limitations.")]))
    return page("", "Skala-1.1 reaction barriers and paths: companion site", body, ["js/plot.js", "js/index.js"],
                "Data, code and interactive views for a study of Skala-1.1 barriers and reaction paths of five small molecules.", anonymous)


# ---------------------------------------------------------------- results
def results_page(R, anonymous):
    rx = R["reactions"]
    xtb = json.load(open(os.path.join(D, "records", "phase13_xtb_barriers.json")))["barrier_TS_minus_reactant_or_complex"]
    lit = {r["key"]: r["ref"] for r in rx}
    xrows = []
    for k in RX_ORDER:
        v = xtb.get(k)
        r = [r for r in rx if r["key"] == k][0]
        note = ""
        if k == "h2co":
            note = "GFN2-xTB found no peak on the scan"
        if k == "ch3nc":
            note = "reference unverified"
        xrows.append([esc(r["label"]), fmt(v) if v is not None else "–", fmt(lit[k]), sgn(v - lit[k]) if v is not None else "–", esc(note)])
    xt = table(["Reaction", "GFN2-xTB barrier at its own geometry", "Literature", "Difference", "Note"], xrows, num_cols=(1, 2, 3), caption="GFN2-xTB scout: barrier from GFN2-xTB energies at GFN2-xTB geometries (kcal/mol)",
               label="GFN2-xTB scout barriers")
    body = """<h1>Results</h1>
<p class="lead">One row for each reaction and method: the computed barrier, the reference, and the signed error (computed minus reference, kcal/mol). Change the reference, the basis set or the Skala-1.1 variant below; the table and the statistics follow.</p>
<noscript><p class="note warn">This page needs JavaScript to build its table. The same numbers are in <a href="../data/csv/study/master_errors.csv">master_errors.csv</a>.</p></noscript>
<div id="results-app" data-src="../data/json/results.json">
<div class="controls" role="group" aria-label="Settings of the table">
<fieldset><legend>Reference</legend><div class="opts">
<label><input type="radio" name="ref" value="lit" checked> Literature value</label>
<label><input type="radio" name="ref" value="own"> Our CCSD(T)+ΔMP2 (same geometry)</label></div></fieldset>
<fieldset><legend>Basis set</legend><div class="opts">
<label><input type="radio" name="basis" value="main" checked> def2-TZVP (TZVPD for SN2)</label>
<label><input type="radio" name="basis" value="large"> def2-QZVPP (HCN, NH₃, H₂CO only)</label></div></fieldset>
<fieldset><legend>Skala-1.1</legend><div class="opts">
<label><input type="radio" name="d3" value="both" checked> with and without D3</label>
<label><input type="radio" name="d3" value="with"> with D3</label>
<label><input type="radio" name="d3" value="without"> without D3</label></div></fieldset>
<fieldset><legend>Other methods</legend><div class="opts" id="method-boxes"></div></fieldset>
</div>
<p class="small" id="mode-note" role="status"></p>
<div class="tablewrap" tabindex="0" role="region" aria-label="Barriers and errors, sortable"><table id="res-table" class="compact">
<caption>Forward barrier per reaction and method</caption>
<thead><tr></tr></thead><tbody></tbody></table></div>
<p class="small" id="table-foot"></p>

<h2 id="calculator">Mean errors for a set of reactions</h2>
<p>Tick the reactions to include. The mean absolute error (MAE), mean signed error (MSE) and largest absolute error are computed in your browser from the rows above. For the published sets the result is compared with the statistics in <a href="../data/csv/study/summary_stats.csv">summary_stats.csv</a>
and <a href="../data/csv/study/derived_leave_one_out.csv">derived_leave_one_out.csv</a>. CH₃NC is unticked by default because its literature reference could not be verified; tick it to see its effect.</p>
<div class="controls" role="group" aria-label="Reactions in the statistics"><fieldset><legend>Reactions</legend><div class="opts" id="rx-boxes"></div></fieldset>
<div><button type="button" class="btn alt" id="rx-default">Reset to the published set</button></div></div>
<div class="tablewrap" tabindex="0" role="region" aria-label="Mean errors for the chosen reactions"><table id="stat-table">
<caption id="stat-cap">Mean errors</caption>
<thead><tr><th scope="col">Method</th><th scope="col" class="num">n</th><th scope="col" class="num">MAE</th><th scope="col" class="num">MSE</th><th scope="col" class="num">Largest |error|</th><th scope="col">Published value for this set</th></tr></thead>
<tbody></tbody></table></div>
<p class="small" id="stat-note" role="status"></p>

<h2 id="scout">GFN2-xTB scout</h2>
<p>GFN2-xTB is not a functional in the comparison above. It was run as a scout: the same barrier, from GFN2-xTB energies at GFN2-xTB geometries.</p>
%s
<p class="small">Source: <a href="../data/records/phase13_xtb_barriers.json">data/records/phase13_xtb_barriers.json</a> and the literature values of the table above. The effect of using GFN2-xTB geometries with the functionals is in <a href="../data/csv/results/scout_table.csv">scout_table.csv</a> and Figure 11.</p>
</div>
<h2>Where the numbers come from</h2>
<ul>
<li>Barriers and errors: <a href="../data/csv/study/master_errors.csv">data/csv/study/master_errors.csv</a> (one row per scored comparison; this page uses the rows with <code>is_primary_barrier = True</code>).</li>
<li>Errors against our own coupled-cluster value: <a href="../data/csv/study/derived_internal_reference.csv">derived_internal_reference.csv</a>. H₂CO is missing there: the CCSD(T) calculation at its transition state is unreliable (see <a href="../checks/">Checks</a>).</li>
<li>Reverse barriers, reaction energies and SN2 barriers from separated reactants are in <code>master_errors.csv</code> and <a href="../data/csv/tables/tableS1_all_scored_rows.csv">tableS1_all_scored_rows.csv</a>.</li>
</ul>
""" % xt
    return page("results", "Results", body, ["js/plot.js", "js/results.js"], "Barriers and signed errors of Skala-1.1 and four functionals for five small reactions, with a mean-error calculator.", anonymous)


# ---------------------------------------------------------------- demos
def demos_page(anonymous):
    body = """<h1>Interactive demos</h1>
<p class="lead">Five views of the data. Each one reads a JSON file that is built from the CSV files named under it, and says how many points it shows. Nothing is computed here that is not in a CSV file, except sums, means and the parabola fits that the CSV files also report.</p>
<ul class="toc"><li><a href="#demo-paths">Reaction paths</a></li><li><a href="#demo-sn2">SN2 symmetric stretch</a></li><li><a href="#demo-t1">HCN T1 map</a></li><li><a href="#demo-cost">Cost and accuracy</a></li><li><a href="#demo-basis">Basis set and D3</a></li></ul>
<noscript><p class="note warn">The demos need JavaScript.</p></noscript>

<section class="demo" id="demo-paths" aria-labelledby="h-paths">
<h2 id="h-paths">1. Reaction paths: HCN → HNC, NH₃ inversion, H₂CO → H₂ + CO</h2>
<p>Energy along a relaxed path (B3LYP-D3(BJ) geometries, def2-TZVP), relative to the first point. The black crosses are the coupled-cluster reference; coloured symbols are the functionals (their energy is the reference plus their signed error). Move the slider to step along the path and see the geometry.</p>
<div class="controls" id="paths-controls">
<fieldset><legend>Reaction</legend><div class="opts" id="paths-rx"></div></fieldset>
<fieldset><legend>Reference</legend><div class="opts" id="paths-ref"></div></fieldset>
<fieldset><legend>Show</legend><div class="opts"><label><input type="radio" name="paths-y" value="energy" checked> Relative energy</label><label><input type="radio" name="paths-y" value="error"> Signed error</label></div></fieldset>
<fieldset><legend>Methods</legend><div class="opts" id="paths-methods"></div></fieldset>
</div>
<div class="note warn hide" id="paths-warn"><p>H₂CO is a failed path: the scan is not smooth at the transition state, and the coupled-cluster reference is unstable there (points p05 to p06 and the transition state have no usable reference). It is shown for completeness and supports no conclusion.</p></div>
<div class="legend" id="paths-legend"></div>
<div class="two"><div class="chart" id="paths-chart"></div>
<div><label for="paths-slider">Path point: <b id="paths-pt"></b></label><br><input type="range" id="paths-slider" min="0" max="10" value="0" step="1" aria-describedby="paths-readout">
<svg class="mol" id="paths-mol" width="260" height="190" viewBox="0 0 260 190" role="img" aria-label="Geometry at the selected path point"></svg>
<div class="readout" id="paths-readout" role="status" aria-live="polite"></div></div></div>
<p class="small" id="paths-data"></p>
</section>

<section class="demo" id="demo-sn2" aria-labelledby="h-sn2">
<h2 id="h-sn2">2. SN2 symmetric stretch: where is the transition state of each method?</h2>
<p>Both C–X bonds are stretched together (the rest of the molecule fixed). For a D₃h SN2 transition state this symmetric stretch has the transition state as a <em>minimum</em>, so the lowest point of each curve is that method's barrier, measured from the ion-molecule complex. Open symbols are calculations whose SCF did not converge: they are shown but not used. The slider steps through the calculated distances.</p>
<div class="controls">
<fieldset><legend>Reaction</legend><div class="opts" id="sn2-rx"></div></fieldset>
<fieldset><legend>Curves</legend><div class="opts">
<label><input type="checkbox" id="sn2-skala" checked> Skala-1.1</label><label><input type="checkbox" id="sn2-wb97" checked> ωB97M-V</label>
<label><input type="checkbox" id="sn2-nod3"> Skala-1.1 without D3</label></div></fieldset>
<fieldset><legend>Annotations</legend><div class="opts">
<label><input type="checkbox" id="sn2-fit" checked> Parabola through three points, and its minimum</label>
<label><input type="checkbox" id="sn2-refs" checked> Reference barriers and B3LYP transition state</label></div></fieldset>
</div>
<div class="legend" id="sn2-legend"></div>
<div class="chart" id="sn2-chart"></div>
<label for="sn2-slider">C–X distance: <b id="sn2-d"></b></label><br><input type="range" id="sn2-slider" min="0" max="10" value="0" step="1" aria-describedby="sn2-readout">
<div class="readout" id="sn2-readout" role="status" aria-live="polite"></div>
<div class="tablewrap" tabindex="0" role="region" aria-label="Minima of the stretch curves"><table id="sn2-min" class="compact"></table></div>
<p class="small" id="sn2-data"></p>
</section>

<section class="demo" id="demo-t1" aria-labelledby="h-t1">
<h2 id="h-t1">3. HCN T1 map: where is the coupled-cluster reference reliable?</h2>
<p>T1 diagnostic of CCSD(T)/def2-SVP on a rigid grid: angle θ (0° = H on the carbon side) against distance r of H from the centre of mass of CN. Cells above the chosen T1 threshold are orange. Cells more than 120 kcal/mol above the lowest grid energy (outside the \"window\") are faded. The usual flag levels are 0.02 and 0.03.</p>
<div class="controls"><div><label for="t1-thr">T1 threshold: <b id="t1-val"></b></label><br><input type="range" id="t1-thr" min="0.010" max="0.060" step="0.001" value="0.020"></div>
<fieldset><legend>Preset</legend><div class="pillrow" id="t1-presets"></div></fieldset></div>
<div class="readout" id="t1-counts" role="status" aria-live="polite"></div>
<div class="chart"><canvas id="t1-canvas" class="map" tabindex="0" role="img" aria-label="HCN T1 map; use the arrow keys to read single cells"></canvas></div>
<div class="legend"><span><span class="sw" id="t1-sw-ok"></span>T1 below threshold</span><span><span class="sw" id="t1-sw-bad"></span>T1 at or above threshold</span><span><span class="sw" style="background:var(--grid)"></span>outside the 120 kcal/mol window (faded)</span></div>
<div class="readout" id="t1-cell" role="status" aria-live="polite">Hover over a cell, or focus the map and use the arrow keys.</div>
<p class="small" id="t1-data"></p>
</section>

<section class="demo" id="demo-cost" aria-labelledby="h-cost">
<h2 id="h-cost">4. Cost against accuracy</h2>
<p>Median wall time of one single-point calculation against the mean absolute error of the five-reaction headline set. Wall times are from one run of the GPU notebook (Tesla T4) and indicate order of magnitude, not a benchmark.</p>
<div class="controls"><fieldset><legend>Structures</legend><div class="opts" id="cost-sets"></div></fieldset></div>
<div class="legend" id="cost-legend"></div>
<div class="chart" id="cost-chart"></div>
<div class="tablewrap" tabindex="0" role="region" aria-label="Times and errors"><table id="cost-table" class="compact"></table></div>
<p class="readout" id="cost-ratio" role="status"></p>
<p class="small" id="cost-data"></p>
</section>

<section class="demo" id="demo-basis" aria-labelledby="h-basis">
<h2 id="h-basis">5. Basis-set and D3 effects</h2>
<p>Left: the change of the signed error from def2-TZVP to def2-QZVPP for the three reactions that were run at both (arrow from TZVP to QZVPP). Right: the D3 term of Skala-1.1, the barrier with D3 minus the barrier without D3 from the same SCF, for the six reactions.</p>
<div class="controls"><fieldset><legend>Left panel shows</legend><div class="opts"><label><input type="radio" name="basis-y" value="err" checked> Signed error</label><label><input type="radio" name="basis-y" value="bar"> Change of the barrier</label></div></fieldset></div>
<div class="legend" id="basis-legend"></div>
<div class="two"><div class="chart" id="basis-chart"></div><div class="chart" id="d3-chart"></div></div>
<div class="tablewrap" tabindex="0" role="region" aria-label="Basis set and D3 numbers"><table id="basis-table" class="compact"></table></div>
<p class="small" id="basis-data"></p>
</section>
"""
    return page("demos", "Demos", body, ["js/plot.js", "js/demos_a.js", "js/demos_b.js"], "Interactive views of the reaction paths, the SN2 stretch, the HCN T1 map, cost and accuracy, and basis-set and D3 effects.", anonymous)


# ---------------------------------------------------------------- figures
def figures_page(anonymous):
    figs = json.load(open(os.path.join(D, "json", "figure_captions.json"), encoding="utf-8"))
    out = ["<h1>Figures</h1>", '<p class="lead">The figures of the paper and the supplementary information, numbered as in those documents. Each figure is available as PNG, PDF and SVG, with the CSV file that holds the plotted values.</p>']
    out.append('<ul class="toc"><li><a href="#main">Main figures</a></li><li><a href="#si">Supplementary figures</a></li></ul>')
    for si in (False, True):
        out.append('<h2 id="%s">%s</h2><div class="figgrid">' % ("si" if si else "main", "Supplementary figures" if si else "Main figures"))
        for f in figs:
            if f["si"] != si:
                continue
            fid = f["id"]
            png = os.path.join(SITE, "figures", "png", fid + ".png")
            w, h = png_size(png)
            data = "data/csv/figures/%s_data.csv" % fid
            has_data = os.path.exists(os.path.join(SITE, data))
            links = ['<a href="../figures/png/%s.png">PNG</a>' % fid, '<a href="../figures/pdf/%s.pdf">PDF</a>' % fid, '<a href="../figures/svg/%s.svg">SVG</a>' % fid]
            if has_data:
                links.append('<a href="../%s">CSV data</a>' % data)
            anchor = f["label"].lower().replace(" ", "-")
            out.append('<figure class="fig" id="%s"><a href="../figures/png/%s.png"><img src="../figures/png/%s.png" width="%d" height="%d" loading="lazy" alt="%s"></a>'
                       '<figcaption><b>%s. %s.</b> %s</figcaption><div class="dl">%s<span class="small mut">%s</span></div></figure>' % (
                           anchor, fid, fid, w, h, esc("%s: %s" % (f["label"], f["title"])), esc(f["label"]), esc(f["title"]), esc(f["caption"]), " ".join(links), esc(fid)))
        out.append("</div>")
    out.append('<h2 id="summary-graphic">Summary graphic</h2><figure class="fig" id="graphical-abstract"><a href="../figures/png/graphical_abstract.png"><img src="../figures/png/graphical_abstract.png" width="%d" height="%d" loading="lazy" alt="%s"></a>'
               '<figcaption>Signed errors of Skala-1.1, \u03c9B97M-V and PBE for the five reactions (left; the grey band is \u00b11 kcal/mol) and the symmetric stretch of Cl\u207b + CH\u2083Cl (right). Values from <a href="../data/csv/study/master_errors.csv">master_errors.csv</a>. PNG only.</figcaption></figure>' % (
                   png_size(os.path.join(SITE, "figures", "png", "graphical_abstract.png")) + (esc("Summary graphic: signed errors of Skala-1.1, wB97M-V and PBE for five reactions, and the SN2 symmetric stretch"),)))
    out.append('<p class="small">The script that draws the figures is in <code>figures/code/</code>; see <a href="../code/#figbuild">Code: rebuild the figures</a>.</p>')
    out.append('<p class="small">Colour-vision checks of the figures (simulated deuteranopia and grey scale) were made for the paper; the symbols differ in shape as well as colour. Figure files are listed with checksums in <a href="../SHA256SUMS.txt">SHA256SUMS.txt</a>.</p>')
    return page("figures", "Figures", "\n".join(out), [], "All figures of the paper and the supplementary information with captions, as PNG, PDF and SVG.", anonymous)


# ---------------------------------------------------------------- data
def file_info(path):
    with open(path, newline="", encoding="utf-8") as f:
        rd_ = csv.reader(f)
        hdr = next(rd_)
        n = sum(1 for _ in rd_)
    return hdr, n


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def data_page(anonymous, strict=True):
    figs = {f["id"]: f for f in json.load(open(os.path.join(D, "json", "figure_captions.json"), encoding="utf-8"))}
    groups = {g[0]: g for g in dd.GROUPS}
    out = ["<h1>Data</h1>",
           '<p class="lead">Every file behind the numbers and figures: the study tables, the paper tables, the data of each figure, and the calculation records of the run. Column meanings are listed under each file. Energies are in kcal/mol unless the column name says hartree; distances in ångström; times in seconds.</p>',
           '<div class="btnrow"><a class="btn" href="../b1_geometries.zip" download>All geometries (zip)</a><a class="btn alt" href="../SHA256SUMS.txt">SHA-256 checksums</a><a class="btn alt" href="#geometries">Geometry files</a><a class="btn alt" href="#records">JSON records</a></div>',
           '<p><label for="filesearch">Find a file or column: </label><input type="search" id="filesearch" class="search" placeholder="e.g. signed_error, T1, phase11" aria-controls="filelist"></p>',
           '<p class="small" id="filecount" role="status"></p><div id="filelist">']
    undocumented = []
    toc = []
    for gk in ("study", "tables", "figures", "results"):
        g = groups[gk]
        files = sorted(glob.glob(os.path.join(D, "csv", gk, "*.csv")))
        out.append('<h2 id="%s">%s</h2><p>%s</p>' % (gk, esc(g[1]), esc(g[2])))
        toc.append('<li><a href="#%s">%s</a></li>' % (gk, esc(g[1])))
        for f in files:
            b = os.path.basename(f)
            hdr, n = file_info(f)
            size = os.path.getsize(f)
            if b in dd.FILES:
                title, source, notes = dd.FILES[b]
            elif b.endswith("_data.csv") and re.match(r"fig\w*?\d+_", b) and [i for i in figs if i.split("_")[0] == b.split("_")[0]]:
                fg = figs[[i for i in figs if i.split("_")[0] == b.split("_")[0]][0]]
                title, source, notes = "Data of %s" % fg["label"], "Figure script", "Values drawn in %s (%s)." % (fg["label"], fg["title"])
            else:
                undocumented.append(b); title, source, notes = b, "", ""
            crow = []
            for c in hdr:
                text, found = dd.describe(b, c)
                if not found:
                    undocumented.append("%s:%s" % (b, c))
                crow.append("<tr><td><code>%s</code></td><td>%s</td></tr>" % (esc(c), esc(text)))
            rel = "data/csv/%s/%s" % (gk, b)
            out.append('<details class="panel filerow" data-search="%s"><summary><code>%s</code> <span class="mut">– %s</span></summary>'
                       '<p>%s</p><p class="small">%d rows, %d columns, %s. Produced by: %s. <a href="../%s" download>Download</a> &middot; SHA-256 <code class="hash">%s</code></p>'
                       '<div class="tablewrap" tabindex="0" role="region" aria-label="Columns of %s"><table class="compact"><thead><tr><th scope="col">Column</th><th scope="col">Meaning</th></tr></thead><tbody>%s</tbody></table></div></details>' % (
                           esc((b + " " + title + " " + " ".join(hdr)).lower()), esc(b), esc(title), esc(notes) if notes else "", n, len(hdr), human_size(size), esc(source) if source else "the notebooks",
                           esc(rel), sha256(f), esc(b), "".join(crow)))
    if undocumented and strict:
        raise SystemExit("not in the data dictionary: %s" % undocumented)
    out.append("</div>")
    # geometries
    gdir = os.path.join(D, "geometries")
    gl = sorted(os.path.basename(p) for p in glob.glob(os.path.join(gdir, "*.xyz")))
    sysnames = ["hcn", "nh3", "h2co", "ch3nc", "sn2_ff", "sn2_clcl", "ch4"]
    labels = dict(RX_LABEL); labels["ch4"] = "CH₄ C–H stretch (reference test)"
    rows = []
    for s in sysnames:
        for g in gl:
            if g.startswith(s + "_"):
                rest = g[len(s) + 1:-4]
                m = re.match(r"(.+)_(own_b3lyp|xtb_gfn2|path_own_scan)$", rest)
                species, level = (m.group(1), m.group(2)) if m else (rest, "")
                if level == "path_own_scan":
                    species = "path (one frame per point)"
                lv = {"own_b3lyp": "B3LYP-D3(BJ) optimised", "xtb_gfn2": "GFN2-xTB", "path_own_scan": "relaxed scan, B3LYP-D3(BJ)"}.get(level, level)
                rows.append([esc(labels[s]), esc(species), esc(lv), '<a href="../data/geometries/%s">%s</a>' % (esc(g), esc(g)), human_size(os.path.getsize(os.path.join(gdir, g)))])
    out.append('<h2 id="geometries">Geometries</h2><p>Cartesian coordinates in ångström in the xyz format (first line: number of atoms). The path files hold one frame per path point in the order of the path tables. '
               'Species names: TS = transition state; RC = reaction complex; HCN, HNC, NH3 and so on are minima; Clm and Fm are the chloride and fluoride ions. One zip with all files: <a href="../b1_geometries.zip">b1_geometries.zip</a>.</p>')
    out.append(table(["System", "Species", "Level", "File", "Size"], rows, caption=None, label="Geometry files", cls="compact"))
    # records
    recs = [("config.json", "Settings of the run (design, reactions, deadline)."), ("tiers.json", "Reference tier of each system."), ("phase10_timing.json", "GPU notebook gate: routes, hardware, software versions, Skala weights check."),
            ("phase13_xtb_barriers.json", "GFN2-xTB barriers at GFN2-xTB geometries (used on the Results page)."), ("geompack_installed.json", "Check of the geometry pack at installation (hash and structure count)."),
            ("phase00_verify_summary.json", "Summary of the gradient test of the stored geometries."), ("phase05_hcn_grid.npz", "The rigid HCN grid as NumPy arrays (the same values as phase05_hcn_grid.csv)."),
            ("path_frames.json", "Atoms and coordinates of every path point (read by the path demo).")]
    out.append('<h2 id="records">JSON records</h2>')
    out.append(table(["File", "Content", "SHA-256 (first 16)"], [['<a href="../data/records/%s">%s</a>' % (r, r), esc(t), '<code class="hash">%s</code>' % sha256(os.path.join(D, "records", r))[:16]] for r, t in recs],
                     label="JSON records", cls="compact"))
    out.append('<h2 id="json">Files the pages read</h2><p>The interactive pages load these files, which <code>tools/build_site.py</code> writes from the CSV files above: <a href="../data/json/results.json">results.json</a>, <a href="../data/json/paths.json">paths.json</a>, '
               '<a href="../data/json/sn2.json">sn2.json</a>, <a href="../data/json/t1grid.json">t1grid.json</a>, <a href="../data/json/cost.json">cost.json</a>, <a href="../data/json/basis.json">basis.json</a>, '
               '<a href="../data/json/figure_captions.json">figure_captions.json</a>. Each has a <code>sources</code> entry naming its CSV files.</p>')
    out.append('<h2 id="hashes">Checksums</h2><p>SHA-256 of every file in <code>data/</code>, <code>code/</code> and <code>figures/</code>: <a href="../SHA256SUMS.txt">SHA256SUMS.txt</a>. Check a download with <code>sha256sum -c SHA256SUMS.txt</code> from the site folder. '
               'The geometry pack of the notebooks has its own hash on the <a href="../code/#geompack">Code page</a>.</p>')
    out.insert(2, '<ul class="toc">%s<li><a href="#geometries">Geometries</a></li><li><a href="#records">JSON records</a></li><li><a href="#hashes">Checksums</a></li></ul>' % "".join(toc))
    return page("data", "Data", "\n".join(out), ["js/data.js"], "Every CSV file of the study with a data dictionary, geometries in xyz format and checksums.", anonymous)


# ---------------------------------------------------------------- code
MODULE_INFO = [
    ("b1v5_common.py", "Shared helpers: result tables (one CSV per phase, written after every calculation), marker files, geometry files and the database of the six reactions."),
    ("b1v5_cpu.py", "CPU calculations: B3LYP-D3(BJ) optimisation, CCSD(T) with T1 and stability, DF-MP2, GFN2-xTB."),
    ("b1v5_gpu.py", "GPU calculations: the five functionals, the SCF rescue ladder, geometry optimiser, phases G-1 to G4 (full-size notebook)."),
    ("b1v5_gpu_figs.py", "Figures, tables and the frozen-criteria outcome table, regenerated from the CSV files."),
    ("b1v5_1h.py", "The deadline queue and the six-reaction design of the one-hour notebooks."),
    ("b1v5_1h_cpu.py", "The phases of the one-hour CPU notebook."),
    ("b1v5_1h_gpu.py", "The phases of the one-hour GPU notebook."),
    ("b1v5_1h_check.py", "The automatic check and the frozen criteria for the six-reaction design."),
    ("b1v5_check.py", "The automatic check: compares with the literature table, runs the internal sanity checks, writes check_report.csv."),
    ("b1v5_lit_table.py", "The literature table: every reference value with its source, zero-point convention and verification status."),
    ("b1v5_v4_stored.py", "Stored numbers of an earlier version of the study, used to check that the code reproduces them."),
    ("b1v5_geompack.py", "The geometry pack: builder, hash check and installation into the working folder."),
    ("b1v5_diag_sn2.py", "Diagnostic of the Cl⁻ + CH₃Cl transition state with Skala-1.1 (SCF state, density fitting, symmetric stretch)."),
    ("build_cpu_1h.py", "Builds the CPU notebook from the modules (needs nbformat)."),
    ("build_gpu_1h.py", "Builds the GPU notebook from the modules."),
    ("build_diag_sn2.py", "Builds the first SN2 diagnostic notebook."),
    ("build_diag_sym.py", "Builds the symmetric-stretch diagnostic notebook."),
    ("build_collect_zip.py", "Builds the notebook that packs the results folder into one zip."),
    ("build_geompack.py", "Builds the geometry pack (optimisations, Hessians, scans, xTB geometries)."),
    ("check_test_synthetic.py", "Test of the check on synthetic data (made-up energies; no chemistry)."),
    ("test_1h_check.py", "Test of the one-hour check and frozen criteria on synthetic data."),
    ("test_1h_notebooks.py", "Runs the notebooks cell by cell in a stand-in mode that needs no GPU and no Colab."),
    ("run_cells.py", "Runs the code cells of a notebook in one namespace, as Colab does (used by the tests)."),
]
NB_INFO = [
    ("B1v5_CPU_1h_v2.ipynb", "CPU notebook (reference side)", "Colab CPU runtime, 2 vCPU", "about 28 minutes for the core, up to 50 minutes with the tail",
     "Checks the geometry pack, then CCSD(T) (frozen core, T1, RHF stability) at the stationary points and along the paths, the MP2 basis correction, GFN2-xTB scout points, one all-electron CCSD(T) barrier, the HCN grid, and the automatic check.",
     "phase00 to phase08 CSV files, phase05_hcn_grid.csv, check_report.csv, RUN_SUMMARY.md"),
    ("B1v5_GPU_1h_v3.ipynb", "GPU notebook (functionals)", "Colab GPU runtime, T4 is enough", "about 50 minutes (deadline queue)",
     "Skala-1.1 (with D3, and without D3 from the same SCF), ωB97M-V, B3LYP-D3(BJ), r²SCAN and PBE at the stationary points, Skala-1.1, ωB97M-V and B3LYP-D3(BJ) along the paths, at GFN2-xTB geometries, def2-QZVPP for HCN, NH₃, H₂CO; then figures, tables, frozen criteria and the check.",
     "phase10 to phase13 CSV files, findings_outcome.csv, check_report.csv, figures, tables"),
    ("B1v5_GPU_diag_sn2.ipynb", "Diagnostic 1: Cl⁻ + CH₃Cl with Skala-1.1", "Colab GPU runtime, T4", "about 25 minutes",
     "Tests whether the Skala-1.1 energy at the transition state is a genuine SCF state (density fitting, other starting densities, xTB geometry).", "diag_sn2_clcl.csv, DIAG_SN2_SUMMARY.md"),
    ("B1v5_GPU_diag_sym.ipynb", "Diagnostic 2: symmetric stretch", "Colab GPU runtime, T4", "20 to 30 minutes (cap 40)",
     "Scans the symmetric C–X stretch of Cl⁻ + CH₃Cl and F⁻ + CH₃F until the minimum is inside the range, for Skala-1.1 and ωB97M-V (demo 2 on the Demos page).", "diag_sn2_clcl.csv (extended), DIAG_SYM_SUMMARY.md, figure"),
    ("B1v5_collect_zip.ipynb", "Collect the results", "any Colab runtime", "about 1 minute",
     "Packs the results folder in Google Drive into one zip and starts the download. Leaves out the cache folder (Skala weights). Changes nothing in Drive.", "one zip file"),
]


def code_page(anonymous):
    pre = "../"
    nbrows = []
    for fn, title, rt, tm, what, outs in NB_INFO:
        nbrows.append('<div class="panel"><h3>%s <span class="tag">%s</span></h3><p>%s</p><p class="small">Runtime: %s. Time: %s. Writes: %s.</p><div class="btnrow">%s</div></div>' % (
            esc(title), esc(fn), esc(what), esc(rt), esc(tm), esc(outs), colab_buttons("code/notebooks/" + fn, pre)))
    mods = [[ '<a href="../code/modules/%s">%s</a>' % (m, m), esc(t)] for m, t in MODULE_INFO]
    gp_sha = open(os.path.join(SITE, "code", "geompack", "b1v5_geompack_1h.json.sha256")).read().split()[0]
    gp_real = sha256(os.path.join(SITE, "code", "geompack", "b1v5_geompack_1h.json"))
    assert gp_sha == gp_real, (gp_sha, gp_real)
    body = """<h1>Code and notebooks</h1>
<p class="lead">The calculations ran in Google Colab notebooks. The notebooks are published without outputs; the Python modules they paste into their cells are published as plain files. You can re-run the whole study on Colab, or only re-check the published numbers on your own computer in a few seconds.</p>
<div data-cfg-warning class="note warn hide"><p>The Colab and GitHub buttons use the repository name set in <code>js/config.js</code>, which still holds a placeholder. Set <code>repo</code> there (see README.md) and the buttons work.</p></div>
<ul class="toc"><li><a href="#notebooks">Notebooks</a></li><li><a href="#run">How to run</a></li><li><a href="#ladder">SCF rescue ladder</a></li><li><a href="#quirks">Quirks</a></li><li><a href="#recheck">Re-check the numbers</a></li><li><a href="#figbuild">Rebuild the figures</a></li><li><a href="#modules">Modules</a></li><li><a href="#geompack">Geometry pack</a></li></ul>

<h2 id="notebooks">Notebooks</h2>
<p>Run the CPU notebook first, then the GPU notebook. The two diagnostics and the collector are optional.</p>
%s

<h2 id="run">How to run on Colab</h2>
<ol class="steps">
<li>Click <b>Open in Colab</b> for the CPU notebook. (If the button is dashed, <code>js/config.js</code> still holds the placeholder repository name; download the notebook and upload it to Colab instead.)</li>
<li>Use the default CPU runtime. For the GPU notebook choose <i>Runtime &gt; Change runtime type &gt; GPU</i>; a T4 is enough, an L4 or A100 reaches further down the list of calculations.</li>
<li>Choose <i>Runtime &gt; Run all</i>. The first cell mounts Google Drive and shows a popup: click <b>Allow</b> at once. The notebook retries the mount up to three times.</li>
<li>Wait. Every calculation is appended to a CSV file in the Drive folder <code>MyDrive/skala-b1-v5</code> as soon as it finishes, so a disconnect loses nothing.</li>
<li>If the runtime stops or a calculation was refused by the deadline queue, choose <i>Run all</i> again: finished calculations are skipped and the notebook continues where it stopped.</li>
<li>Run the GPU notebook after the CPU notebook: its last phases read the coupled-cluster results from the same Drive folder. At the end the folder holds the CSV files that this site publishes in <a href="../data/#results">data / results</a>.</li>
<li>To fetch the folder, run <code>B1v5_collect_zip</code> and download the zip.</li>
</ol>
<p><b>Expected times.</b> CPU notebook: about 28 minutes for the core list and up to 50 minutes with the tail. GPU notebook: about 50 minutes on a T4 (it completes the core list there). Diagnostic 1: about 25 minutes. Diagnostic 2: 20 to 30 minutes. Collector: about a minute. The times are estimates for a free Colab runtime and vary with the hardware you are given.</p>
<p><b>Outputs.</b> One CSV per phase (<code>phase03_ccsdt.csv</code>, <code>phase11_dft.csv</code>, ...), geometries, <code>check_report.csv</code> and <code>check_report.md</code>, <code>findings_outcome.csv</code>, tables, figures and <code>RUN_SUMMARY.md</code>. The published copies are in <a href="../data/">data</a>.</p>
<p><b>The deadline queue.</b> Both one-hour notebooks use one ordered list of calculations. Before each calculation the notebook asks whether it can finish before minute 50 (elapsed time plus 1.5 times the running mean time of that kind of calculation must stay below 50 minutes, and the phase must stay inside its own cap). A calculation that does not fit is not started; it is listed as \"not run (time)\" in <code>RUN_SUMMARY.md</code> and <code>queue_not_run.csv</code>. All of them were completed in a later session of the study.</p>

<h2 id="ladder">SCF rescue ladder</h2>
<p>When a self-consistent-field calculation does not converge, the notebooks try these steps in order and stop at the first that converges. The energy table records which attempt gave each energy (column <code>attempt</code>).</p>
<ol class="steps"><li>Default start (the usual initial guess).</li><li>Warm start from the density of a neighbouring point or a related calculation.</li><li>Fresh start with a different initial guess.</li><li>Level shift.</li><li>Start from the PBE density, with a level shift of 0.3 and up to 300 cycles.</li><li>Density fitting, with a level shift of 0.5 and up to 400 cycles.</li></ol>
<p>The Cl⁻ + CH₃Cl transition state with Skala-1.1 converged only on the last rung. Diagnostic 1 shows that density fitting changes that energy by about 0.05 kcal/mol, and the symmetric-stretch demo shows it with its neighbours. A barrier never mixes density-fitted and unfitted energies (check I14).</p>

<h2 id="quirks">Quirks to know</h2>
<ul>
<li><b>Drive first.</b> The mount popup is the only click the notebooks need. If you are late, the mount fails and is retried; if all three attempts fail, the notebook works in a local folder, and the last cell zips and downloads the results instead of releasing the runtime.</li>
<li><b>The last cell releases the runtime</b> to stop billing, unless Drive was not mounted. Set <code>AUTO_DISCONNECT = False</code> in the settings cell to keep it.</li>
<li><b>Failed points.</b> Set <code>RETRY_FAILED = True</code> to retry calculations that failed, and <code>FORCE_RERUN = True</code> to recompute finished ones.</li>
<li><b>Stand-in mode.</b> <code>STANDIN = True</code> runs the whole notebook on a plain CPU with cheap functionals in tiny basis sets, to test the machinery. Its numbers are not results.</li>
<li><b>Extension items.</b> <code>EXTENSION</code> flags run extra reactions (F⁻ + CH₃Cl and OH⁻ + CH₃F) in a session of their own. They are not part of the paper.</li>
<li><b>Skala without D3</b> is the same SCF with the stored D3 energy removed; it is not a second calculation.</li>
<li><b>Geometries are not recomputed.</b> They come from the geometry pack (below); <code>REBUILD_GEOMS = True</code> rebuilds them, which takes hours.</li>
</ul>

<h2 id="recheck">Re-check the numbers on your computer</h2>
<p>The automatic check that the notebooks run at the end is a plain Python module. <code>code/recheck_numbers.py</code> runs it again on the published calculation records and compares the result with the published <a href="../data/csv/results/check_report.csv">check_report.csv</a>. Needs Python 3 and numpy.</p>
<pre><code>git clone &lt;the repository&gt; &amp;&amp; cd &lt;the repository&gt;
python code/recheck_numbers.py</code></pre>
<p>Expected output:</p>
<pre><code>items: 139 (published report: 139)
counts: PASS 54, WARN 14, FAIL 5, NOT RUN 40, N/A 4, INFO 22
check_report.csv: identical to the published file</code></pre>
<p>The result of the check, with every warning and failure explained, is on the <a href="../checks/">Checks page</a>.</p>

<h2 id="figbuild">Rebuild the figures</h2>
<p>The figures, the tables and their captions are drawn by a Python script that reads only the CSV, JSON and xyz files of this repository (<code>data/</code>). It needs Python 3, numpy, pandas, matplotlib and Pillow (<code>requirements.txt</code>), takes less than a minute, and writes to a new folder, so the published files are not overwritten.</p>
<pre><code>pip install -r requirements.txt
python figures/code/build_figures.py
python figures/code/compare_with_published.py</code></pre>
<p>The second script compares the rebuilt data files, tables and PNG pixels with the published ones and prints <code>all agree</code>. Files and details: <a href="../figures/code/README.md">figures/code/README.md</a>.</p>

<h2 id="modules">Modules</h2>
<p>The notebooks contain these files pasted verbatim into their code cells (edit the module, then rebuild the notebook with the matching <code>build_*.py</code>).</p>
%s

<h2 id="geompack">Geometry pack</h2>
<p>The optimisations, transition states (one imaginary mode each), Hessians, zero-point energies, relaxed scans and GFN2-xTB geometries were computed once and are stored in one JSON file that each notebook carries (compressed) and verifies by its SHA-256 before use.</p>
<dl class="kv"><dt>File</dt><dd><a href="../code/geompack/b1v5_geompack_1h.json">code/geompack/b1v5_geompack_1h.json</a></dd><dt>SHA-256</dt><dd><code>%s</code></dd><dt>Read me</dt><dd><a href="../code/geompack/GEOMPACK_README.md">GEOMPACK_README.md</a> (format, contents, how it was built)</dd></dl>
<p>Check it yourself: <code>sha256sum code/geompack/b1v5_geompack_1h.json</code> should print the hash above. The same structures are available as xyz files on the <a href="../data/#geometries">Data page</a>.</p>
""" % ("\n".join(nbrows), table(["Module", "What it does"], mods, label="Python modules", cls="compact"), gp_sha)
    return page("code", "Code", body, [], "Colab notebooks, Python modules, run instructions, the SCF rescue ladder and a script that re-checks the published numbers.", anonymous)


# ---------------------------------------------------------------- checks
EXPLAIN = {
    ("I2", "hcn"): "T1 of the CCSD(T) calculation reaches 0.033 at path point p05 (the transition state: 0.026), above the 0.03 flag. The coupled-cluster reference near the HCN transition state is less certain than elsewhere. Point p05 is left out of the reference statistics and its neighbours are flagged.",
    ("I3", "hcn"): "The restricted Hartree-Fock solution is unstable at path points p03 and p04, so they are also left out of the reference. The HCN path reference has 8 of 11 points once the MP2 basis correction is applied.",
    ("I3", "h2co"): "The Hartree-Fock solution is unstable at the transition state and at path points p05 and p06. There is no usable CCSD(T) barrier for H₂CO: its reference is the literature value (81.9 kcal/mol) only, and the path reference has no values around the transition state.",
    ("I3", "ch4"): "The Hartree-Fock solution is unstable at the longest C–H stretch (p08, 1.70 Å), so the CH₄ reference is cut there. CH₄ is a bond-stretch test of the reference, not one of the five reactions.",
    ("I8", "h2co"): "Changing from def2-TZVP to def2-QZVPP lowers the Skala-1.1 barrier of H₂CO by 2.32 kcal/mol (ωB97M-V: 0.17), above the 1.0 limit. This is the largest basis-set sensitivity in the data. The headline table uses def2-TZVP (error +1.08); at def2-QZVPP the error is −1.24.",
    ("R2", "L1 hcn:barrier_fwd"): "ωB97M-V at def2-QZVPP gives 44.28 against the reference 48.07 (−3.79). The check suggests a wrong geometry level or functional string. We did not find a cause; the value is reported as it is and remains unexplained.",
    ("R2", "L3 hcn:reaction_energy"): "The same ωB97M-V anomaly appears in the HNC minus HCN energy (−2.89).",
    ("R1", "L6 ch3nc:barrier_fwd"): "The CH₃NC literature value is from a conference abstract and its zero-point convention is unknown, so it is not scored. Our CCSD(T)+ΔMP2 value is shown for information.",
    ("R2", "L6 ch3nc:barrier_fwd"): "As above: the value of ωB97M-V is shown but not scored; the large basis was not run for CH₃NC.",
    ("R3", "L6 ch3nc:barrier_fwd"): "The literature value of the CH₃NC barrier is unverified. It is excluded from every error statistic on this site until the source is confirmed.",
    ("R2", "L7 sn2_ff:barrier_fwd"): "The barrier of F⁻ + CH₃F from separated reactants is a small number (−0.34), so the relative tolerance rule applies; the difference is 0.61 kcal/mol. The large basis was not run. No action.",
    ("R3", "L19 ch3nc:reaction_energy"): "There is no literature value for the CH₃NC reaction energy; the table row is a placeholder and is in no statistic.",
    ("I2", "h2co"): "T1 is 0.023 at the transition state and at p05, above the 0.02 flag: the H₂CO coupled-cluster calculation is borderline even before the instability found by I3.",
    ("I2", "sn2_ff"): "T1 is 0.023 for F⁻ and 0.020 at the transition state, slightly above the 0.02 flag. Our CCSD(T)+ΔMP2 barrier agrees with the literature to +0.22 kcal/mol.",
    ("I5", "h2co"): "The path has one maximum, but there is no coupled-cluster energy at the transition state to compare it with, because that point is excluded (I3).",
    ("I10", "hcn"): "The all-electron CCSD(T) barrier of HCN is 0.35 kcal/mol above the frozen-core barrier, over the 0.3 limit. Both values are reported.",
    ("I12", "h2co"): "GFN2-xTB found no maximum along the H₂CO scan; the stored xTB transition state is the middle scan point, which has the B3LYP transition-state coordinates. The H₂CO geometry penalty therefore tests only the minimum.",
    ("I12", "ch3nc CCSD(T)"): "At the GFN2-xTB geometry the CCSD(T) barrier of CH₃NC is 42.87 against 39.97 at the B3LYP geometry (+2.90). The xTB geometry is poor for this molecule (C–N–C angle 60.2° against 69.0°), not specifically for Skala-1.1.",
    ("I12", "sn2_ff wB97M-V"): "At the GFN2-xTB geometry (asymmetric, shorter C–F bond) ωB97M-V gives 11.37 against 12.56 at the B3LYP geometry (−1.19); several functionals show a lower barrier there.",
}
ITEMS = [("R1", "Our CCSD(T) (with the MP2 basis correction) against the literature value."), ("R2", "ωB97M-V in the larger basis against the literature value (anchor functional)."),
         ("R3", "References that are unverified or missing: never scored."), ("R4", "Signed errors of all functionals against the literature (reported, not scored)."),
         ("P1, P2", "Geometry pack: hash check, and B3LYP gradient test of the stored structures."), ("I1", "SCF or calculation convergence of every stationary point and path point."),
         ("I2", "T1 diagnostic of each CCSD(T) point: flag above 0.02, fail above 0.03."), ("I3", "Stability of the Hartree-Fock solution of each CCSD(T) point."),
         ("I4", "Frequencies: one imaginary mode at the transition state, none at minima."), ("I5", "The maximum of the path against the transition-state energy."),
         ("I6", "The transition state lies above the minima."), ("I7", "Ordering of the functionals (PBE lowest)."),
         ("I8", "Change of the barrier from def2-TZVP to def2-QZVPP."), ("I9", "Reproduction of stored numbers of an earlier version."), ("I10", "All-electron against frozen-core CCSD(T) (HCN)."),
         ("I11", "Own geometry against database geometry; Skala-relaxed geometries."), ("I12", "Scout test: barrier at GFN2-xTB geometries against B3LYP geometries."),
         ("I13", "Spread of the standard functionals (information)."), ("I14", "No barrier mixes density-fitted and unfitted energies."), ("I15", "Size of the MP2 basis-set correction of the reference.")]


def checks_page(anonymous):
    rows = rd("results/check_report.csv")
    counts = {}
    for r in rows:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    cards = ['<div class="grid">']
    for st in ("PASS", "WARN", "FAIL", "INFO", "NOT RUN", "N/A"):
        cards.append('<div class="card"><p>%s</p><p class="big">%d</p></div>' % (st, counts.get(st, 0)))
    cards.append("</div>")
    ex_rows = []
    missing = []
    for r in rows:
        if r["status"] in ("WARN", "FAIL"):
            key = (r["item"], r["scope"])
            if key not in EXPLAIN:
                missing.append(key)
                continue
            val = " ".join(x for x in (r["value"], ("against " + r["reference"]) if r["reference"] else "") if x)
            ex_rows.append(['<span class="tag %s">%s</span>' % ("bad" if r["status"] == "FAIL" else "warn", r["status"]), esc(r["item"]), esc(r["scope"]), esc(val), esc(EXPLAIN[key])])
    if missing:
        raise SystemExit("no explanation for: %s" % missing)
    ex_rows.sort(key=lambda x: (0 if "FAIL" in x[0] else 1))
    allrows = []
    for r in rows:
        st = r["status"]
        allrows.append('<tr data-status="%s"><td class="status-%s">%s</td><td>%s</td><td>%s</td><td class="num">%s</td><td class="num">%s</td><td class="num">%s</td><td>%s</td><td>%s</td></tr>' % (
            esc(st), esc(st.replace(" ", "-")), esc(st), esc(r["item"]), esc(r["scope"]), esc(r["value"]), esc(r["reference"]), esc(r["difference"]), esc(r["tolerance"]), esc(r["message"])))
    fo = rd("results/findings_outcome.csv")
    t4 = {r["finding"]: r for r in rd("tables/table4_findings.csv")}
    fr = []
    for r in fo:
        t = t4[r["finding"]]
        vcls = {"supports": "ok", "negative": "bad", "inconclusive": "warn"}.get(r["verdict"], "")
        fr.append([esc(r["finding"]), esc(r["description"]), esc(r["criterion_supports"]), esc(r["criterion_negative"]), esc(r["criterion_inconclusive"]), esc(r["measured"]), esc(r["n"]),
                   '<span class="tag %s">%s</span>' % (vcls, esc(r["verdict"])), esc(t["sensitivity_note"])])
    frozen = fo[0]["frozen_on"]
    body = """<h1>Checks</h1>
<p class="lead">At the end of the run, an automatic check compared the results with the literature table and ran internal sanity tests. It produced %d lines. Nothing was hidden: every warning and failure is explained below, and the full report is searchable here and in <a href="../data/csv/results/check_report.csv">check_report.csv</a>.</p>
%s
<p>NOT RUN is not a failure: those items belong to reactions or tests outside the six-reaction design (extension reactions, tests that need data not collected), or to systems with too few points for the test. You can reproduce the report with <a href="../code/#recheck">code/recheck_numbers.py</a>.</p>
<h2 id="explained">Warnings and failures, explained</h2>
%s
<h2 id="items">What the items mean</h2>
%s
<h2 id="criteria">The six frozen criteria</h2>
<p>Six findings were stated with their criteria for \"supports\", \"negative\" and \"inconclusive\" before any result existed (%s). The table gives the criteria as written, the measured values and the verdict. The notes explain how sensitive the verdict is. Source: <a href="../data/csv/results/findings_outcome.csv">findings_outcome.csv</a> and <a href="../data/csv/tables/table4_findings.csv">table4_findings.csv</a>.</p>
%s
<h2 id="full">The full report</h2>
<div class="controls"><fieldset><legend>Status</legend><div class="pillrow" id="status-pills"></div></fieldset>
<div><label for="checksearch">Search: </label><input type="search" id="checksearch" class="search" placeholder="e.g. hcn, I3, wB97M-V"></div></div>
<p class="small" id="check-count" role="status"></p>
<div class="tablewrap" tabindex="0" role="region" aria-label="Full check report"><table id="check-table" class="compact"><thead><tr><th scope="col">Status</th><th scope="col">Item</th><th scope="col">Scope</th><th scope="col" class="num">Value</th><th scope="col" class="num">Reference</th><th scope="col" class="num">Difference</th><th scope="col">Tolerance</th><th scope="col">Message</th></tr></thead>
<tbody>%s</tbody></table></div>
""" % (len(rows), "".join(cards), table(["Status", "Item", "Scope", "Value", "What it means"], ex_rows, label="Warnings and failures", cls="compact"),
       table(["Item", "Meaning"], [[esc(a), esc(b)] for a, b in ITEMS], cls="compact", label="Check items"), esc(frozen),
       table(["", "Finding", "Supports if", "Negative if", "Inconclusive otherwise", "Measured", "n", "Verdict", "Sensitivity of the verdict"], fr, cls="compact", label="Frozen criteria"), "\n".join(allrows))
    return page("checks", "Checks", body, ["js/checks.js"], "The automatic check of the run with every warning and failure explained, and the six frozen criteria.", anonymous)


# ---------------------------------------------------------------- about
def about_page(anonymous):
    t = {r["item"]: r["setting"] for r in rd("tables/tableS6_software_hardware.csv")}
    ver = [["Skala", "skala 2026.9 (PyPI, August 2026); GPU build skala-cuda12x 2026.9"], ["Skala weights", "microsoft/skala-1.1, Hugging Face revision <code>99b5ed87e5f69d9216e1f9e30148b922eaea1241</code>"],
           ["Skala weights SHA-256", "<code>%s</code>" % esc(t["Skala weights SHA-256"])], ["PySCF", "2.12.1 (GPU notebook), 2.15.0 (CPU notebook and geometry pack)"],
           ["GPU4PySCF", "gpu4pyscf-cuda12x 1.8.1 with cupy-cuda12x 14.2.0"], ["Dispersion", "pyscf-dispersion 1.7.0 (D3(BJ))"], ["Geometry optimiser", "geomeTRIC 1.1.1"],
           ["GFN2-xTB", "tblite 0.8.0"], ["Python", "3.13"], ["Geometry pack SHA-256", "<code>%s</code>" % esc(open(os.path.join(SITE, "code", "geompack", "b1v5_geompack_1h.json.sha256")).read().split()[0])]]
    lim = ["<b>Five reactions.</b> The headline statistics rest on n = 5 reactions (n = 3 at the larger basis, n = 4 for CCSD(T)). The 95 %% bootstrap interval of the Skala-1.1 mean absolute error is 0.66 to 2.42 kcal/mol; differences between functionals of a few tenths of a kcal/mol are not resolved.",
           "<b>One geometry protocol.</b> All headline energies are single points at B3LYP-D3(BJ) geometries. GFN2-xTB geometries were used only in the scout test, and Skala-1.1-relaxed geometries were not computed.",
           "<b>An unexplained ωB97M-V value.</b> ωB97M-V underestimates the HCN → HNC barrier by 3.5 kcal/mol (and the HCN → HNC reaction energy by 2.7). We did not find the cause. It dominates the ωB97M-V mean absolute error: without HCN it is 0.46 instead of 1.08.",
           "<b>References.</b> The five literature values of the headline set were checked in their sources (the H₂CO value from the abstract of the paper); the H₂CO value includes zero-point energy, which was corrected with a B3LYP harmonic zero-point difference (−4.99 kcal/mol) for every method. The CH₃NC value could not be verified and is in no statistic. Our CCSD(T) is unusable for H₂CO (unstable reference at the transition state).",
           "<b>No stability analysis of the DFT solutions.</b> Convergence of each SCF was recorded, but no test whether a lower-energy solution exists. The Cl⁻ + CH₃Cl transition state with Skala-1.1 converged only with density fitting.",
           "<b>Skala-1.1 without D3</b> is the same SCF with the D3 energy removed, not a separate optimisation of the functional without D3.",
           "<b>Paths.</b> The HCN path reference has 8 of 11 points; the H₂CO path is a failed path and supports no conclusion.",
           "<b>Timings</b> come from one run on Colab-class hardware and are not a benchmark.",
           "<b>Not on this site:</b> the text of the manuscript and any correspondence. This site is the data, code and views behind it."]
    body = """<h1>About</h1>
<div data-cfg-warning class="note warn hide"><p>Some values in <code>js/config.js</code> (repository, DOI, contact, copyright holder) are still placeholders. Replace them before you publish; see README.md.</p></div>
<h2 id="cite">How to cite</h2>
<div class="panel cite"><pre id="citation" data-citation></pre><button type="button" class="btn alt" data-copy="#citation">Copy</button></div>
<p class="small">The archive entry (DOI) is a copy of the same files; if the website and the archive ever differ, the archive is the one to cite.</p>
<h2 id="licence">Licence</h2>
<ul><li><b>Code</b> (the Python modules, notebooks and the scripts in <code>tools/</code>, <code>code/</code> and <code>figures/code/</code>, and the site code in <code>js/</code> and <code>css/</code>): MIT licence, file <a href="../LICENSE">LICENSE</a>.</li>
<li><b>Data, figures and tables</b> (everything in <code>data/</code> and <code>figures/</code> except the scripts in <code>figures/code/</code>): Creative Commons Attribution 4.0 International (CC BY 4.0), file <a href="../LICENSE-DATA.txt">LICENSE-DATA.txt</a>. You may reuse and adapt them, also commercially, if you give credit.</li></ul>
<h2 id="versions">Software versions</h2>
%s
<p class="small">Hardware: CPU notebook Intel Xeon 2.20 GHz, 2 cores, 13 GB; GPU notebook Intel Xeon 2.00 GHz, 8 cores, 51 GB with a Tesla T4. Full list: <a href="../data/csv/tables/tableS6_software_hardware.csv">tableS6_software_hardware.csv</a>, <a href="../data/records/phase10_timing.json">phase10_timing.json</a>.</p>
<h2 id="limits">Limitations</h2>
<ul>%s</ul>
<h2 id="contact">Contact</h2>
<p data-hide-anon>Questions and corrections: <span data-cfg="contact"></span>. Please say which page and which file you mean.</p>
<p data-show-anon>Contact details are withheld in this anonymous copy.</p>
<h2>How this site is made</h2>
<p>The pages are plain HTML, CSS and JavaScript with no external libraries, fonts or services. <code>tools/build_site.py</code> rebuilds every page and JSON file from the CSV files in <code>data/</code>; <code>tools/sync_sources.py</code> copies the sources into this folder. The site stores nothing about you except your choice of light or dark theme, kept in your browser.</p>
""" % (table(["Item", "Setting"], [[esc(a), b] for a, b in ver], cls="compact", label="Software versions"), "".join("<li>%s</li>" % x.replace("%%", "%") for x in lim))
    return page("about", "About", body, [], "How to cite the study, licences, software versions and limitations.", anonymous)
