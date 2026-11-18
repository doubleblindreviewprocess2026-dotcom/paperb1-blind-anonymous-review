/* Demos 1 and 2: reaction paths and the SN2 symmetric stretch. Data: data/json/paths.json and data/json/sn2.json. */
(function () {
  "use strict";
  var $ = function (id) { return document.getElementById(id); };
  var fmt = B1.fmt, sgn = B1.sgn, esc = B1.esc;
  var NS = "http://www.w3.org/2000/svg";

  function swatch(m) {
    return '<span class="sw" style="' + (m === "Skala-1.1-noD3" ? "background:transparent;border:2px solid var(" + B1.METHODS[m].color + ")" : "background:var(" + B1.METHODS[m].color + ")") + '"></span>';
  }

  // =================================================================== paths
  (function () {
    if (!$("paths-chart")) return;
    var P, S = { rx: "hcn", ref: "CCSD(T)+dMP2", y: "energy", idx: 0, on: { "Skala-1.1": true, "Skala-1.1-noD3": false, "wB97M-V": true, "B3LYP-D3BJ": true } };
    var METHODS = ["Skala-1.1", "Skala-1.1-noD3", "wB97M-V", "B3LYP-D3BJ"];
    var REFLAB = { "CCSD(T)+dMP2": "CCSD(T)+ΔMP2 (MP2 basis correction)", "CCSD(T) raw": "CCSD(T) raw (def2-TZVP)" };
    var RADIUS = { H: 0.31, C: 0.76, N: 0.71, O: 0.66 }, COL = { H: "#cbd5e1", C: "#4b5563", N: "#2563eb", O: "#dc2626" };

    B1.loadJSON("../data/json/paths.json").then(function (j) { P = j; init(); }).catch(function (e) { B1.showLoadError($("paths-chart"), e); });

    function init() {
      var box = $("paths-rx");
      Object.keys(P).forEach(function (k) {
        var l = document.createElement("label");
        l.innerHTML = '<input type="radio" name="paths-rx" value="' + k + '"' + (k === S.rx ? " checked" : "") + "> " + esc(P[k].label);
        box.appendChild(l);
        l.firstChild.addEventListener("change", function () { S.rx = k; S.idx = 0; S.ref = P[k].series["CCSD(T)+dMP2"] ? "CCSD(T)+dMP2" : "CCSD(T) raw"; refs(); draw(); });
      });
      var mb = $("paths-methods");
      METHODS.forEach(function (m) {
        var l = document.createElement("label");
        l.innerHTML = '<input type="checkbox"' + (S.on[m] ? " checked" : "") + "> " + swatch(m) + esc(m === "Skala-1.1-noD3" ? "Skala-1.1 without D3" : B1.METHODS[m].label);
        mb.appendChild(l);
        l.firstChild.addEventListener("change", function (ev) { S.on[m] = ev.target.checked; draw(); });
      });
      Array.prototype.forEach.call(document.getElementsByName("paths-y"), function (r) { r.addEventListener("change", function () { S.y = r.value; draw(); }); });
      $("paths-slider").addEventListener("input", function (ev) { S.idx = +ev.target.value; draw(); });
      refs(); draw();
      B1.onResize($("paths-chart"), draw);
    }
    function refs() {
      var box = $("paths-ref"), ser = P[S.rx].series;
      box.innerHTML = "";
      Object.keys(ser).sort().reverse().forEach(function (k) {
        var l = document.createElement("label");
        l.innerHTML = '<input type="radio" name="paths-ref" value="' + esc(k) + '"' + (k === S.ref ? " checked" : "") + "> " + esc(REFLAB[k] || k);
        box.appendChild(l);
        l.firstChild.addEventListener("change", function () { S.ref = k; draw(); });
      });
    }

    function draw() {
      var d = P[S.rx], ser = d.series[S.ref], pts = d.points, n = pts.length;
      $("paths-warn").classList.toggle("hide", S.rx !== "h2co");
      var sl = $("paths-slider"); sl.max = n - 1; sl.value = S.idx;
      var cs = pts.map(function (p) { return p.coord; }), lo = Math.min.apply(null, cs), hi = Math.max.apply(null, cs), pad = (hi - lo) * 0.04;
      var xrev = cs[0] > cs[n - 1];
      // y range
      var ys = [];
      pts.forEach(function (p) {
        var r = ser.ref[p.id];
        if (S.y === "energy") { if (r !== undefined) ys.push(r); METHODS.forEach(function (m) { if (S.on[m] && ser.err[m] && ser.err[m][p.id] !== undefined && r !== undefined) ys.push(r + ser.err[m][p.id]); }); }
        else METHODS.forEach(function (m) { if (S.on[m] && ser.err[m] && ser.err[m][p.id] !== undefined) ys.push(ser.err[m][p.id]); });
      });
      var ylo = Math.min.apply(null, ys.concat(S.y === "error" ? [0] : [0])), yhi = Math.max.apply(null, ys.concat([0])), yp = (yhi - ylo) * 0.08 || 1;
      var legend = '<span><span class="sw" style="background:var(--c-cc)"></span>' + esc(REFLAB[S.ref]) + " (reference)</span>" + METHODS.filter(function (m) { return S.on[m]; }).map(function (m) { return "<span>" + swatch(m) + esc(m === "Skala-1.1-noD3" ? "Skala-1.1 without D3" : B1.METHODS[m].label) + "</span>"; }).join("");
      $("paths-legend").innerHTML = legend;
      var el = $("paths-chart");
      var p = new B1.Plot(el, { height: 300, x: [lo - pad, hi + pad], xrev: xrev, y: [ylo - yp, yhi + yp], xlabel: d.coord_label, ylabel: S.y === "energy" ? "Energy relative to the first point (kcal/mol)" : "Signed error (kcal/mol)",
        label: "Energy or signed error along the path of " + d.label, zero: S.y === "error" });
      if (S.y === "error") p.band(-1, 1);
      p.vline(pts[S.idx].coord, { color: "--focus", width: 2, dash: "4 3" });
      if (S.y === "energy") {
        var line = [];
        pts.forEach(function (q) { line.push(ser.ref[q.id] !== undefined ? [q.coord, ser.ref[q.id]] : null); });
        p.line(line.filter(Boolean), { color: "--c-cc", width: 1, dash: "2 3" });
        pts.forEach(function (q) { if (ser.ref[q.id] !== undefined) p.marker(q.coord, ser.ref[q.id], "x", { color: "--c-cc", size: 7, top: true, title: "reference " + q.id + ": " + fmt(ser.ref[q.id]) }); });
      }
      METHODS.forEach(function (m) {
        if (!S.on[m] || !ser.err[m]) return;
        var meta = B1.METHODS[m], line = [];
        pts.forEach(function (q) {
          var e = ser.err[m][q.id], r = ser.ref[q.id];
          if (e === undefined || (S.y === "energy" && r === undefined)) return;
          var y = S.y === "energy" ? r + e : e;
          line.push([q.coord, y]);
        });
        p.line(line, { color: meta.color, width: 1.5 });
        line.forEach(function (pt) { p.marker(pt[0], pt[1], meta.marker, { color: meta.color, size: 4.5, title: m + ": " + fmt(pt[1]) }); });
      });
      // flag points without a reference or excluded from the reference statistics
      var missing = pts.filter(function (q) { return ser.ref[q.id] === undefined; }).map(function (q) { return q.id; });
      readout(d, ser, missing);
      molecule(d);
      $("paths-data").textContent = "Data: " + d.sources + ". " + n + " path points; reference values at " + (n - missing.length) + " of them" + (missing.length ? " (missing: " + missing.join(", ") + ", excluded because T1 is above 0.03 or the Hartree-Fock solution is unstable)" : "") + ". Methods: " + METHODS.length + " functional variants. Basis " + d.basis + ".";
    }

    function readout(d, ser, missing) {
      var q = d.points[S.idx], r = ser.ref[q.id], html = "<b>" + q.id + "</b>: " + esc(d.coord_label) + " = <b>" + fmt(q.coord, Math.abs(q.coord) < 10 ? 3 : 1) + "</b><br>";
      if (r === undefined) html += "No reference value at this point (excluded from the reference).<br>";
      else html += "Reference: <b>" + fmt(r) + "</b> kcal/mol above the first point.<br>";
      METHODS.forEach(function (m) {
        if (!S.on[m]) return;
        var e = ser.err[m] && ser.err[m][q.id];
        html += swatch(m) + esc(m === "Skala-1.1-noD3" ? "Skala-1.1 without D3" : B1.METHODS[m].label) + ": " + (e === undefined || r === undefined ? "no value against this reference" : "<b>" + fmt(r + e) + "</b> (error " + sgn(e) + ")") + "<br>";
      });
      $("paths-readout").innerHTML = html;
      $("paths-pt").textContent = q.id + " (" + fmt(q.coord, Math.abs(q.coord) < 10 ? 3 : 1) + ")";
    }

    function molecule(d) {
      var svg = $("paths-mol"), pts = d.points, q = pts[S.idx], el = d.elements, mx = 0, my = 0;
      pts.forEach(function (p) { p.xyz.forEach(function (a) { mx = Math.max(mx, Math.abs(a[0])); my = Math.max(my, Math.abs(a[1])); }); });
      var sc = Math.min(110 / (mx || 1), 70 / (my || 1)), cx = 130, cy = 95, html = "", i, j;
      var order = q.xyz.map(function (a, k) { return k; }).sort(function (a, b) { return q.xyz[a][2] - q.xyz[b][2]; });
      for (i = 0; i < el.length; i++) for (j = i + 1; j < el.length; j++) {
        var a = q.xyz[i], b = q.xyz[j], dist = Math.sqrt(Math.pow(a[0] - b[0], 2) + Math.pow(a[1] - b[1], 2) + Math.pow(a[2] - b[2], 2));
        if (dist < 1.3 * (RADIUS[el[i]] + RADIUS[el[j]])) html += '<line x1="' + (cx + a[0] * sc) + '" y1="' + (cy - a[1] * sc) + '" x2="' + (cx + b[0] * sc) + '" y2="' + (cy - b[1] * sc) + '"/>';
      }
      order.forEach(function (k) {
        var a = q.xyz[k], r = el[k] === "H" ? 9 : 13;
        html += '<circle cx="' + (cx + a[0] * sc) + '" cy="' + (cy - a[1] * sc) + '" r="' + r + '" style="fill:' + COL[el[k]] + '"/><text x="' + (cx + a[0] * sc) + '" y="' + (cy - a[1] * sc) + '"' + (el[k] === "H" ? ' style="fill:#111"' : "") + ">" + el[k] + "</text>";
      });
      svg.innerHTML = html;
    }
  })();

  // =================================================================== SN2
  (function () {
    if (!$("sn2-chart")) return;
    var X, S = { rx: "sn2_clcl", idx: null };
    B1.loadJSON("../data/json/sn2.json").then(function (j) { X = j; init(); }).catch(function (e) { B1.showLoadError($("sn2-chart"), e); });
    function init() {
      var box = $("sn2-rx");
      ["sn2_clcl", "sn2_ff"].forEach(function (k) {
        var l = document.createElement("label");
        l.innerHTML = '<input type="radio" name="sn2-rx" value="' + k + '"' + (k === S.rx ? " checked" : "") + "> " + esc(X[k].label);
        box.appendChild(l);
        l.firstChild.addEventListener("change", function () { S.rx = k; S.idx = null; draw(); });
      });
      ["sn2-skala", "sn2-wb97", "sn2-nod3", "sn2-fit", "sn2-refs"].forEach(function (id) { $(id).addEventListener("change", draw); });
      $("sn2-slider").addEventListener("input", function (ev) { S.idx = +ev.target.value; draw(); });
      draw(); B1.onResize($("sn2-chart"), draw);
    }
    function dists(e) {
      var s = {};
      ["Skala-1.1", "wB97M-V"].forEach(function (m) { e.methods[m].points.forEach(function (p) { s[p.x.toFixed(4)] = p.x; }); });
      return Object.keys(s).map(function (k) { return s[k]; }).sort(function (a, b) { return a - b; });
    }
    function draw() {
      var e = X[S.rx], ds = dists(e), sl = $("sn2-slider");
      if (S.idx === null) { S.idx = 0; ds.forEach(function (d, k) { if (Math.abs(d - e.ts_x) < 5e-4) S.idx = k; }); }
      sl.max = ds.length - 1; if (S.idx > ds.length - 1) S.idx = 0; sl.value = S.idx;
      var showS = $("sn2-skala").checked, showW = $("sn2-wb97").checked, showN = $("sn2-nod3").checked, fit = $("sn2-fit").checked, refs = $("sn2-refs").checked;
      var series = [];
      if (showS) series.push({ m: "Skala-1.1", y: "y", col: "--c-skala", shape: "circle" });
      if (showN) series.push({ m: "Skala-1.1", y: "y_noD3", col: "--c-skala", shape: "plus", dash: "5 3", name: "Skala-1.1 without D3" });
      if (showW) series.push({ m: "wB97M-V", y: "y", col: "--c-wb97", shape: "square" });
      var ys = [], xs = ds.slice();
      series.forEach(function (s) { e.methods[s.m].points.forEach(function (p) { if (p.conv && p[s.y] !== null) ys.push(p[s.y]); }); });
      if (refs) { ys.push(e.ref_ccsdt, e.ref_lit); }
      if (!ys.length) ys = [0, 10];
      var ylo = Math.floor(Math.min.apply(null, ys) - 1), yhi = Math.ceil(Math.max.apply(null, ys) + 1);
      var lo = Math.min.apply(null, xs) - 0.03, hi = Math.max.apply(null, xs) + 0.03;
      var legend = "";
      series.forEach(function (s) { legend += "<span>" + swatch(s.m === "Skala-1.1" && s.y === "y_noD3" ? "Skala-1.1-noD3" : s.m) + esc(s.name || B1.METHODS[s.m].label) + "</span>"; });
      legend += '<span><svg width="14" height="14" aria-hidden="true"><circle cx="7" cy="7" r="4.5" fill="var(--card)" stroke="var(--fg)" stroke-width="2"/></svg>SCF not converged (not used)</span>';
      if (fit) legend += '<span><svg width="16" height="14" aria-hidden="true"><path d="M8 1 L14 7 L8 13 L2 7Z" fill="var(--fg)"/></svg>minimum of the parabola</span>';
      $("sn2-legend").innerHTML = legend;
      var el = $("sn2-chart");
      var p = new B1.Plot(el, { height: 330, x: [lo, hi], y: [ylo, yhi], xlabel: "C–X distance, both bonds (Å)", ylabel: "Energy relative to the complex (kcal/mol)", zero: false,
        label: "Energy of the symmetric stretch of " + e.label + " for Skala-1.1 and ωB97M-V" });
      if (refs) {
        p.hline(e.ref_ccsdt, { color: "--c-cc", dash: "6 3" }); p.hline(e.ref_lit, { color: "--c-cc", dash: "1 3", width: 1.5 });
        p.text(p.x1 - 4, p.sy(e.ref_ccsdt) - 4, "CCSD(T)+ΔMP2 barrier " + fmt(e.ref_ccsdt), { anchor: "end", cls: "ax" });
        p.text(p.x1 - 4, p.sy(e.ref_lit) + 12, "literature " + fmt(e.ref_lit), { anchor: "end", cls: "ax" });
        p.vline(e.ts_x, { color: "--mut", dash: "2 3" });
        p.text(p.sx(e.ts_x) + 3, p.y1 + 10, "B3LYP TS", { cls: "ax" });
      }
      p.vline(ds[S.idx], { color: "--focus", width: 2, dash: "4 3" });
      series.forEach(function (s) {
        var md = e.methods[s.m], conv = md.points.filter(function (q) { return q.conv && q[s.y] !== null; }).sort(function (a, b) { return a.x - b.x; });
        var ln = p.line(conv.map(function (q) { return [q.x, q[s.y]]; }), { color: s.col, width: 1.3, dash: s.dash });
        conv.forEach(function (q) { p.marker(q.x, q[s.y], s.shape, { color: s.col, size: 4.6, title: (s.name || s.m) + ", " + q.x.toFixed(3) + " Å: " + fmt(q[s.y]) }); });
        if (s.y === "y") md.points.filter(function (q) { return !q.conv; }).forEach(function (q) { p.marker(q.x, q[s.y], "ring", { color: s.col, size: 5, title: s.m + ", " + q.x.toFixed(3) + " Å: SCF not converged (" + q.attempt + ")" }); });
        if (fit && s.y === "y") {
          var c = md.coef, a = Math.min.apply(null, md.fit_x) - 0.12, b = Math.max.apply(null, md.fit_x) + 0.12;
          p.curve(function (x) { return c[0] * x * x + c[1] * x + c[2]; }, a, b, { color: s.col, width: 1, dash: "2 2" });
          p.marker(md.vertex.x, md.vertex.y, "diamond", { color: s.col, size: 5.5, top: true, title: s.m + " minimum: " + md.vertex.x.toFixed(3) + " Å, " + fmt(md.vertex.y) + " kcal/mol" });
        }
      });
      readout(e, ds[S.idx], showS, showW, showN);
      minTable(e);
      $("sn2-d").textContent = ds[S.idx].toFixed(3) + " Å (point " + (S.idx + 1) + " of " + ds.length + ")";
      var npts = 0; ["Skala-1.1", "wB97M-V"].forEach(function (m) { npts += e.methods[m].points.length; });
      $("sn2-data").textContent = "Data: " + e.sources + ". " + npts + " calculated points for this reaction (" + e.methods["Skala-1.1"].points.length + " Skala-1.1, " + e.methods["wB97M-V"].points.length + " ωB97M-V), all with density fitting, def2-TZVPD. The parabola goes through the three converged points named in the table; its minimum is computed from them.";
    }
    function at(md, d) { return md.points.filter(function (q) { return Math.abs(q.x - d) < 5e-4; }); }
    function readout(e, d, showS, showW, showN) {
      var html = "C–X = <b>" + d.toFixed(3) + " Å</b><br>", vs = {}, one = function (name, m, key, col) {
        var pts = at(e.methods[m], d), cv = pts.filter(function (q) { return q.conv && q[key] !== null; });
        if (cv.length) { vs[name] = cv[0][key]; return swatch(m) + name + ": <b>" + fmt(cv[0][key]) + "</b> kcal/mol" + (cv[0].ts ? " (the B3LYP transition-state distance)" : "") + " <span class=\"small mut\">(" + esc(cv[0].attempt) + ")</span><br>"; }
        if (pts.length) return swatch(m) + name + ": SCF did not converge (" + pts.length + " attempt" + (pts.length > 1 ? "s" : "") + ": " + pts.map(function (q) { return esc(q.attempt.replace(" (not converged)", "")); }).join(", ") + "); the last iterate, " + fmt(pts[0][key]) + ", is not used.<br>";
        return swatch(m) + name + ": no calculation at this distance.<br>";
      };
      if (showS) html += one("Skala-1.1", "Skala-1.1", "y");
      if (showN) html += one("Skala-1.1 without D3", "Skala-1.1", "y_noD3");
      if (showW) html += one("ωB97M-V", "wB97M-V", "y");
      if (vs["Skala-1.1"] !== undefined && vs["ωB97M-V"] !== undefined) html += "Skala-1.1 minus ωB97M-V: <b>" + sgn(vs["Skala-1.1"] - vs["ωB97M-V"]) + "</b> kcal/mol";
      $("sn2-readout").innerHTML = html;
    }
    function minTable(e) {
      var rows = "";
      ["Skala-1.1", "wB97M-V"].forEach(function (m) {
        var d = e.methods[m];
        rows += "<tr><td>" + swatch(m) + esc(B1.METHODS[m].label) + '</td><td class="num">' + d.vertex.x.toFixed(3) + '</td><td class="num">' + fmt(d.vertex.y) + '</td><td class="num">' + fmt(d.energy_at_ts) + '</td><td class="num">' + sgn(d.err_vs_ccsdt) + '</td><td class="num">' + sgn(d.err_vs_lit) + "</td><td>" + d.fit_x.map(function (x) { return x.toFixed(3); }).join(", ") + "</td></tr>";
      });
      $("sn2-min").innerHTML = "<caption>Minimum of the symmetric stretch, " + esc(e.label) + " (kcal/mol, Å)</caption><thead><tr><th scope=\"col\">Method</th><th scope=\"col\" class=\"num\">Minimum at C–X</th><th scope=\"col\" class=\"num\">Barrier at the minimum</th><th scope=\"col\" class=\"num\">At the B3LYP TS distance (" + e.ts_x.toFixed(3) + ")</th><th scope=\"col\" class=\"num\">Minimum minus CCSD(T)+ΔMP2</th><th scope=\"col\" class=\"num\">Minimum minus literature</th><th scope=\"col\">Three points used (Å)</th></tr></thead><tbody>" + rows + "</tbody>";
    }
  })();
})();
