/* Demos 3 to 5: HCN T1 map, cost against accuracy, basis-set and D3 effects.
   Data: data/json/t1grid.json, cost.json, basis.json (built from the CSV files named under each demo). */
(function () {
  "use strict";
  var $ = function (id) { return document.getElementById(id); };
  var fmt = B1.fmt, sgn = B1.sgn, esc = B1.esc;
  function swatch(m) {
    return '<span class="sw" style="' + (m === "Skala-1.1-noD3" ? "background:transparent;border:2px solid var(" + B1.METHODS[m].color + ")" : "background:var(" + B1.METHODS[m].color + ")") + '"></span>';
  }

  // =================================================================== T1 map
  (function () {
    var cv = $("t1-canvas");
    if (!cv) return;
    var T, thr = 0.02, cur = null, geom = null;
    var nth, nr;
    B1.loadJSON("../data/json/t1grid.json").then(function (j) { T = j; nth = T.theta.length; nr = T.r.length; init(); }).catch(function (e) { B1.showLoadError(cv.parentNode, e); });

    function init() {
      var sl = $("t1-thr");
      sl.addEventListener("input", function () { thr = +sl.value; draw(); });
      [0.02, 0.03].forEach(function (v) {
        var b = document.createElement("button"); b.type = "button"; b.textContent = "T1 = " + v.toFixed(2); b.setAttribute("aria-pressed", "false");
        b.addEventListener("click", function () { sl.value = v; thr = v; draw(); });
        $("t1-presets").appendChild(b);
      });
      cv.addEventListener("mousemove", function (ev) { var c = cellAt(ev); if (c) { cur = c; draw(); } });
      cv.addEventListener("click", function (ev) { var c = cellAt(ev); if (c) { cur = c; draw(); } });
      cv.addEventListener("keydown", function (ev) {
        var k = ev.key, c = cur || [0, 0];
        if (k === "ArrowRight") c = [Math.min(nth - 1, c[0] + 1), c[1]];
        else if (k === "ArrowLeft") c = [Math.max(0, c[0] - 1), c[1]];
        else if (k === "ArrowUp") c = [c[0], Math.min(nr - 1, c[1] + 1)];
        else if (k === "ArrowDown") c = [c[0], Math.max(0, c[1] - 1)];
        else return;
        ev.preventDefault(); cur = c; draw();
      });
      draw(); B1.onResize(cv.parentNode, draw);
    }
    function cellAt(ev) {
      if (!geom) return null;
      var r = cv.getBoundingClientRect(), x = (ev.clientX - r.left) / r.width * geom.W, y = (ev.clientY - r.top) / r.height * geom.H;
      var i = Math.floor((x - geom.l) / geom.cw), j = nr - 1 - Math.floor((y - geom.t) / geom.ch);
      return i >= 0 && i < nth && j >= 0 && j < nr ? [i, j] : null;
    }
    function draw() {
      var W = Math.max(300, Math.floor(cv.parentNode.clientWidth)), H = Math.round(Math.max(240, Math.min(440, W * 0.5))), dpr = window.devicePixelRatio || 1;
      cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); cv.style.width = W + "px"; cv.style.height = H + "px";
      var g = cv.getContext("2d"); g.setTransform(dpr, 0, 0, dpr, 0, 0);
      var l = 48, r = 8, t = 8, b = 38, cw = (W - l - r) / nth, ch = (H - t - b) / nr;
      geom = { W: W, H: H, l: l, t: t, cw: cw, ch: ch };
      var ok = B1.css("--c-skala"), bad = B1.css("--c-wb97"), fg = B1.css("--fg"), mut = B1.css("--mut"), grid = B1.css("--grid"), focus = B1.css("--focus");
      g.clearRect(0, 0, W, H);
      var nWin = 0, okWin = 0, okAll = 0, nAll = 0, i, j;
      for (i = 0; i < nth; i++) for (j = 0; j < nr; j++) {
        var v = T.t1[i][j], inw = T.win[i][j] === 1, below = v < thr;
        nAll++; if (below) okAll++;
        if (inw) { nWin++; if (below) okWin++; }
        g.globalAlpha = inw ? 1 : 0.3;
        g.fillStyle = below ? ok : bad;
        g.fillRect(l + i * cw, t + (nr - 1 - j) * ch, Math.ceil(cw), Math.ceil(ch));
      }
      g.globalAlpha = 1;
      g.strokeStyle = grid; g.lineWidth = 0.5;
      for (i = 0; i <= nth; i++) { g.beginPath(); g.moveTo(l + i * cw, t); g.lineTo(l + i * cw, t + nr * ch); g.stroke(); }
      for (j = 0; j <= nr; j++) { g.beginPath(); g.moveTo(l, t + j * ch); g.lineTo(l + nth * cw, t + j * ch); g.stroke(); }
      g.fillStyle = mut; g.font = "12px system-ui, sans-serif"; g.textAlign = "center";
      var th = T.theta, step = W < 480 ? 60 : 30, k;
      for (k = 0; k < nth; k++) if (th[k] % step === 0) g.fillText(th[k].toFixed(0) + "°", l + (k + 0.5) * cw, t + nr * ch + 14);
      g.fillStyle = fg; g.fillText("θ (0° = H on the C side)", l + nth * cw / 2, H - 6);
      g.textAlign = "right"; g.fillStyle = mut;
      for (k = 0; k < nr; k += 2) g.fillText(T.r[k].toFixed(2), l - 5, t + (nr - 1 - k + 0.5) * ch + 4);
      g.save(); g.translate(11, t + nr * ch / 2); g.rotate(-Math.PI / 2); g.textAlign = "center"; g.fillStyle = fg; g.fillText("r (Å)", 0, 0); g.restore();
      if (cur) { g.strokeStyle = focus; g.lineWidth = 3; g.strokeRect(l + cur[0] * cw + 1, t + (nr - 1 - cur[1]) * ch + 1, cw - 2, ch - 2); }
      $("t1-val").textContent = thr.toFixed(3);
      $("t1-sw-ok").style.background = ok; $("t1-sw-bad").style.background = bad;
      $("t1-counts").innerHTML = "In the " + T.window_kcal + " kcal/mol window: <b>" + okWin + "</b> of <b>" + nWin + "</b> cells have T1 below " + thr.toFixed(3) + " (<b>" + (100 * okWin / nWin).toFixed(1) + " %</b>). All " + nAll + " grid cells: " + okAll + " below.";
      Array.prototype.forEach.call($("t1-presets").children, function (b2) { b2.setAttribute("aria-pressed", Math.abs(parseFloat(b2.textContent.replace("T1 = ", "")) - thr) < 5e-4 ? "true" : "false"); });
      if (cur) {
        var tv = T.t1[cur[0]][cur[1]];
        $("t1-cell").innerHTML = "θ = <b>" + T.theta[cur[0]].toFixed(1) + "°</b>, r = <b>" + T.r[cur[1]].toFixed(3) + " Å</b>: T1 = <b>" + tv.toFixed(4) + "</b> (" + (tv < thr ? "below" : "at or above") + " the threshold); " + fmt(T.erel[cur[0]][cur[1]], 1) + " kcal/mol above the lowest grid energy, " + (T.win[cur[0]][cur[1]] ? "inside" : "outside") + " the window.";
      }
      cv.setAttribute("aria-label", "HCN T1 map: " + okWin + " of " + nWin + " cells in the window have T1 below " + thr.toFixed(3) + ". Use the arrow keys to read single cells.");
      $("t1-data").textContent = "Data: " + T.sources + ". " + nAll + " grid cells (" + nth + " angles x " + nr + " distances); " + nWin + " inside the window. CCSD(T)/def2-SVP, rigid grid, one calculation per cell.";
    }
  })();

  // =================================================================== cost
  (function () {
    if (!$("cost-chart")) return;
    var C, set = "all";
    var METHODS = ["PBE", "r2SCAN", "B3LYP-D3BJ", "Skala-1.1", "wB97M-V"];
    var SETS = { all: "Structures with all five methods (n = 19; the medians of Figure 12)", tz: "All def2-TZVP and def2-TZVPD structures (any method present)", qz: "def2-QZVPP structures" };
    B1.loadJSON("../data/json/cost.json").then(function (j) { C = j; init(); }).catch(function (e) { B1.showLoadError($("cost-chart"), e); });
    function init() {
      var box = $("cost-sets");
      Object.keys(SETS).forEach(function (k) {
        var l = document.createElement("label"); l.innerHTML = '<input type="radio" name="cost-set" value="' + k + '"' + (k === set ? " checked" : "") + "> " + esc(SETS[k]);
        box.appendChild(l); l.firstChild.addEventListener("change", function () { set = k; draw(); });
      });
      $("cost-legend").innerHTML = B1.legendHTML(METHODS);
      draw(); B1.onResize($("cost-chart"), draw);
    }
    function pick() {
      return C.rows.filter(function (r) {
        if (set === "all") return METHODS.every(function (m) { return r[m] !== null; });
        if (set === "tz") return r.basis !== "def2-QZVPP";
        return r.basis === "def2-QZVPP";
      });
    }
    function draw() {
      var rows = pick(), med = {}, n = {};
      METHODS.forEach(function (m) { var v = rows.filter(function (r) { return r[m] !== null; }).map(function (r) { return r[m]; }); n[m] = v.length; med[m] = v.length ? B1.median(v) : null; });
      var vals = METHODS.map(function (m) { return med[m]; }).filter(function (v) { return v !== null; });
      var lo = Math.floor(Math.log(Math.min.apply(null, vals)) / Math.LN10 * 4) / 4 - 0.05, hi = Math.ceil(Math.log(Math.max.apply(null, vals)) / Math.LN10 * 4) / 4 + 0.15;
      var ticksAll = [1, 2, 5, 10, 20, 50, 100, 200], tk = ticksAll.filter(function (t) { var v = Math.log(t) / Math.LN10; return v >= lo && v <= hi; }).map(function (t) { return Math.log(t) / Math.LN10; });
      var p = new B1.Plot($("cost-chart"), { height: 320, x: [lo, hi], y: [0, 6], xticks: tk, xfmt: function (v) { return String(Math.round(Math.pow(10, v))); }, zero: false,
        xlabel: "Median wall time of one single point (s, log scale)", ylabel: "Mean absolute error, five reactions (kcal/mol)", label: "Cost against accuracy for five methods" });
      METHODS.forEach(function (m) {
        if (med[m] === null) return;
        var x = Math.log(med[m]) / Math.LN10, meta = B1.METHODS[m], pub = C.published[m];
        p.line([[x, pub.loo_min], [x, pub.loo_max]], { color: meta.color, width: 1.2 });
        p.marker(x, C.mae[m].mae, meta.marker === "ring" ? "circle" : meta.marker, { color: meta.color, size: 6.5, title: meta.label + ": median " + fmt(med[m]) + " s, MAE " + fmt(C.mae[m].mae) + " kcal/mol" });
        p.text(p.sx(x) + 9, p.sy(C.mae[m].mae) + 4, meta.label, { cls: "ax" });
      });
      var base = med["Skala-1.1"], html = "<caption>Median wall time and mean absolute error (kcal/mol); bar = range of the MAE when one reaction is left out (published)</caption><thead><tr><th scope=\"col\">Method</th><th scope=\"col\" class=\"num\">Structures</th><th scope=\"col\" class=\"num\">Median time (s)</th><th scope=\"col\" class=\"num\">Time / Skala-1.1</th><th scope=\"col\" class=\"num\">MAE</th><th scope=\"col\" class=\"num\">MAE, one reaction left out</th></tr></thead><tbody>";
      METHODS.forEach(function (m) {
        var pub = C.published[m];
        html += "<tr><td>" + swatch(m) + esc(B1.METHODS[m].label) + '</td><td class="num">' + n[m] + '</td><td class="num">' + (med[m] === null ? "–" : fmt(med[m])) + '</td><td class="num">' + (med[m] === null || !base ? "–" : fmt(med[m] / base)) + '</td><td class="num">' + fmt(C.mae[m].mae) + '</td><td class="num">' + fmt(pub.loo_min) + " to " + fmt(pub.loo_max) + "</td></tr>";
      });
      $("cost-table").innerHTML = html + "</tbody>";
      function ratio(a, b) {
        var v = rows.filter(function (r) { return r[a] !== null && r[b] !== null; }).map(function (r) { return r[a] / r[b]; });
        return v.length ? { m: B1.median(v), n: v.length } : null;
      }
      var r1 = ratio("wB97M-V", "Skala-1.1"), r2 = ratio("B3LYP-D3BJ", "Skala-1.1"), r3 = ratio("Skala-1.1", "PBE");
      $("cost-ratio").innerHTML = "Median per-structure ratios in this set: ωB97M-V / Skala-1.1 = <b>" + (r1 ? fmt(r1.m) + "</b> (n = " + r1.n + ")" : "–</b>") + "; B3LYP-D3(BJ) / Skala-1.1 = <b>" + (r2 ? fmt(r2.m) + "</b> (n = " + r2.n + ")" : "–</b>") + "; Skala-1.1 / PBE = <b>" + (r3 ? fmt(r3.m) + "</b> (n = " + r3.n + ")" : "–</b>") + ".";
      $("cost-data").textContent = "Data: " + C.sources + ". " + C.rows.length + " structures in the timing table; " + rows.length + " in the chosen set. The MAE is the headline mean absolute error of each method (n = 5 reactions), whatever set is chosen.";
    }
  })();

  // =================================================================== basis and D3
  (function () {
    if (!$("basis-chart")) return;
    var Bd, mode = "err";
    var METHODS = ["Skala-1.1", "Skala-1.1-noD3", "wB97M-V"];
    var RXL = { hcn: "HCN → HNC", nh3: "NH₃ inversion", h2co: "H₂CO → H₂ + CO", sn2_ff: "F⁻ + CH₃F", sn2_clcl: "Cl⁻ + CH₃Cl", ch3nc: "CH₃NC → CH₃CN" };
    B1.loadJSON("../data/json/basis.json").then(function (j) { Bd = j; init(); }).catch(function (e) { B1.showLoadError($("basis-chart"), e); });
    function init() {
      Array.prototype.forEach.call(document.getElementsByName("basis-y"), function (r) { r.addEventListener("change", function () { mode = r.value; draw(); }); });
      $("basis-legend").innerHTML = B1.legendHTML(METHODS).replace("Skala-1.1 without D3", "Skala-1.1 without D3");
      draw(); B1.onResize($("basis-chart"), draw);
    }
    function get(r, m) { return Bd.basis.filter(function (x) { return x.r === r && x.m === m; })[0]; }
    function draw() {
      var rxs = ["hcn", "nh3", "h2co"], rowsN = rxs.length * 4;
      var vals = [];
      rxs.forEach(function (r) { METHODS.forEach(function (m) { var b = get(r, m); if (mode === "err") vals.push(b.err_main, b.err_large); else vals.push(0, b.shift); }); });
      var lo = Math.floor(Math.min.apply(null, vals.concat([0])) - 0.5), hi = Math.ceil(Math.max.apply(null, vals.concat([0])) + 0.5);
      var short = { "Skala-1.1": "Skala-1.1", "Skala-1.1-noD3": "no D3", "wB97M-V": "\u03c9B97M-V" };
      var padL = 1.7;
      var p = new B1.Plot($("basis-chart"), { height: rowsN * 24 + ($("basis-chart").clientWidth < 480 ? 68 : 52), x: [lo - padL, hi], y: [0, rowsN], yticks: [], zero: false, margin: { l: 12, r: 14, t: 8, b: $("basis-chart").clientWidth < 480 ? 60 : 44 },
        xticks: B1.plotTicks(lo, hi, 6), xlabel: mode === "err" ? "Signed error (kcal/mol); filled = def2-TZVP, open = def2-QZVPP" : "Change of the barrier from def2-TZVP to def2-QZVPP (kcal/mol)", label: "Basis-set dependence for three reactions" });
      if (mode === "err") p.xband(-1, 1);
      p.vline(0, { color: "--fg", width: 1.2 });
      rxs.forEach(function (r, i) {
        var top = rowsN - i * 4;
        if (i > 0) p.hline(top, { color: "--line", width: 1 });
        p.text(p.x0 + 4, p.sy(top - 0.1) + 12, RXL[r], { cls: "lbl" });
        METHODS.forEach(function (m, j) {
          var b = get(r, m), y = top - 1.3 - j * 0.9, meta = B1.METHODS[m];
          p.text(p.x0 + 12, p.sy(y) + 4, short[m], { cls: "ax" });
          if (mode === "err") {
            p.line([[b.err_main, y], [b.err_large, y]], { color: meta.color, width: 1.6 });
            p.marker(b.err_main, y, "circle", { color: meta.color, size: 5, title: m + " " + RXL[r] + ", def2-TZVP: " + sgn(b.err_main) });
            p.marker(b.err_large, y, "ring", { color: meta.color, size: 5, title: m + " " + RXL[r] + ", def2-QZVPP: " + sgn(b.err_large) });
          } else {
            p.line([[0, y], [b.shift, y]], { color: meta.color, width: 5 });
            p.marker(b.shift, y, m === "Skala-1.1-noD3" ? "ring" : "circle", { color: meta.color, size: 4.5, title: m + " " + RXL[r] + ": " + sgn(b.shift) });
          }
        });
      });
      // D3 panel
      var d3 = Bd.d3, n = d3.length, dv = d3.map(function (x) { return x.d3; });
      var dl = Math.floor(Math.min.apply(null, dv.concat([0])) * 10 - 1) / 10, dh = Math.ceil(Math.max.apply(null, dv.concat([0])) * 10 + 1) / 10;
      var q = new B1.Plot($("d3-chart"), { height: n * 40 + 52, x: [dl, dh], y: [0, n], yticks: [], zero: false, margin: { l: 12, r: 14, t: 8, b: 44 },
        xlabel: "D3 contribution to the barrier (kcal/mol)", label: "D3 contribution of Skala-1.1 for six reactions" });
      q.vline(0, { color: "--fg", width: 1.2 });
      d3.forEach(function (x, i) {
        var yc = n - i - 0.68;
        if (i > 0) q.hline(n - i, { color: "--line", width: 1 });
        q.text(q.x0 + 4, q.sy(n - i - 0.1) + 12, RXL[x.r] + (x.headline ? "" : " (not in the headline set)"), { cls: "lbl" });
        q.line([[0, yc], [x.d3, yc]], { color: "--c-skala", width: 6 });
        q.marker(x.d3, yc, "circle", { color: "--c-skala", size: 4.5, title: RXL[x.r] + ": D3 contribution " + sgn(x.d3, 3) });
      });
      // table
      var html = "<caption>Basis set and D3 numbers (kcal/mol)</caption><thead><tr><th scope=\"col\">Reaction</th><th scope=\"col\">Method</th><th scope=\"col\" class=\"num\">Barrier, TZVP</th><th scope=\"col\" class=\"num\">Barrier, QZVPP</th><th scope=\"col\" class=\"num\">Change</th><th scope=\"col\" class=\"num\">Error, TZVP</th><th scope=\"col\" class=\"num\">Error, QZVPP</th></tr></thead><tbody>";
      rxs.forEach(function (r) { METHODS.forEach(function (m) { var b = get(r, m); html += "<tr><td>" + RXL[r] + "</td><td>" + swatch(m) + esc(m === "Skala-1.1-noD3" ? "Skala-1.1 without D3" : B1.METHODS[m].label) + '</td><td class="num">' + fmt(b.main) + '</td><td class="num">' + fmt(b.large) + '</td><td class="num">' + sgn(b.shift) + '</td><td class="num">' + sgn(b.err_main) + '</td><td class="num">' + sgn(b.err_large) + "</td></tr>"; }); });
      html += '</tbody><tbody><tr><th scope="colgroup" colspan="7">D3 contribution of Skala-1.1: barrier with D3 minus barrier without D3 (electronic barrier, same SCF)</th></tr><tr><th scope="col">Reaction</th><th scope="col" class="num">With D3</th><th scope="col" class="num">Without D3</th><th scope="col" class="num">D3 contribution</th><th scope="col" class="num">Error with D3</th><th scope="col" class="num">Error without D3</th><th></th></tr>';
      d3.forEach(function (x) { html += "<tr><td>" + RXL[x.r] + '</td><td class="num">' + fmt(x.with) + '</td><td class="num">' + fmt(x.without) + '</td><td class="num">' + sgn(x.d3, 3) + '</td><td class="num">' + sgn(x.err_with) + '</td><td class="num">' + sgn(x.err_without) + "</td><td></td></tr>"; });
      $("basis-table").innerHTML = html + "</tbody>";
      $("basis-data").textContent = "Data: " + Bd.sources + ". Basis sets: " + Bd.basis.length / 3 + " reactions x 3 methods at both basis sets; D3: " + d3.length + " reactions. The H₂CO barriers with D3 are electronic barriers; the zero-point correction cancels in the D3 difference.";
    }
  })();
})();
