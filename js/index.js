/* Home page: strip chart of the headline signed errors (data/json/results.json). */
(function () {
  "use strict";
  var el = document.getElementById("strip");
  if (!el) return;
  var R;
  var SHOW = ["Skala-1.1", "Skala-1.1-noD3", "wB97M-V", "B3LYP-D3BJ", "r2SCAN", "PBE", "CCSD(T)+dMP2"];
  document.getElementById("legend").innerHTML = B1.legendHTML(SHOW);

  function draw() {
    var rx = R.reactions.filter(function (r) { return r.headline; }), n = rx.length;
    var rows = {};
    R.rows.forEach(function (r) { if (r.b === "main") rows[r.r + "|" + r.m] = r; });
    var lo = 0, hi = 0;
    rx.forEach(function (r) { SHOW.forEach(function (m) { var row = rows[r.key + "|" + m]; if (row && !row.wd) { lo = Math.min(lo, row.err); hi = Math.max(hi, row.err); } }); });
    lo = Math.floor(lo - 0.5); hi = Math.ceil(hi + 0.5);
    var p = new B1.Plot(el, { height: n * 58 + (el.clientWidth < 480 ? 64 : 48), x: [lo, hi], y: [0, n], yticks: [], zero: false, xlabel: "Signed error of the barrier, computed minus reference (kcal/mol)",
      margin: { l: 12, r: 14, t: 8, b: el.clientWidth < 480 ? 60 : 44 }, label: "Strip chart of signed barrier errors for five reactions and seven methods. Values are in the table above." });
    p.xband(-1, 1);
    for (var q = 1; q < n; q++) p.hline(q, { color: "--line", width: 1 });
    p.vline(0, { color: "--fg", width: 1.2 });
    rx.forEach(function (r, i) {
      var yc = n - i - 0.55;
      p.text(p.x0 + 4, p.sy(n - i - 0.06) + 2, r.label, { cls: "lbl" });
      SHOW.forEach(function (m, j) {
        var row = rows[r.key + "|" + m];
        if (!row || row.wd || row.err === null) return;
        var meta = B1.METHODS[m] || { color: "--fg", marker: "circle" };
        var shape = m === "CCSD(T)+dMP2" ? "plus" : meta.marker;
        p.marker(row.err, yc + (3 - j) * 0.075, shape, { color: meta.color, size: 4.6, title: B1.METHODS[m].label + ", " + r.label + ": " + B1.sgn(row.err) + " kcal/mol" });
      });
    });
  }
  B1.loadJSON(el.getAttribute("data-src")).then(function (j) { R = j; draw(); B1.onResize(el, draw); }).catch(function (e) { B1.showLoadError(el, e); });
})();
