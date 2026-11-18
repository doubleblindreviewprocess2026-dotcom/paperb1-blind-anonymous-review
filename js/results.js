/* Results page: sortable table and mean-error calculator. Reads data/json/results.json (built from data/csv/study/*.csv). */
(function () {
  "use strict";
  var app = document.getElementById("results-app");
  if (!app) return;
  var R, rxBy = {}, rowIdx = {}, mOrder = {}, rxOrder = {};
  var CC = ["CCSD(T)-fc", "CCSD(T)+dMP2"];
  var S = { ref: "lit", basis: "main", d3: "both", on: {}, sort: { key: "rx", dir: 1 }, ticked: {} };
  var $ = function (id) { return document.getElementById(id); };
  var fmt = B1.fmt, sgn = B1.sgn, esc = B1.esc;

  B1.loadJSON(app.getAttribute("data-src")).then(function (j) { R = j; init(); }).catch(function (e) { B1.showLoadError($("mode-note"), e); });

  function label(m) { for (var i = 0; i < R.methods.length; i++) if (R.methods[i].key === m) return R.methods[i].label; return m; }

  function init() {
    R.reactions.forEach(function (r, i) { rxBy[r.key] = r; rxOrder[r.key] = i; });
    R.methods.forEach(function (m, i) { mOrder[m.key] = i; S.on[m.key] = true; });
    R.rows.forEach(function (r) { rowIdx[r.r + "|" + r.m + "|" + r.b] = r; });
    // method boxes (Skala variants are controlled by the D3 radio)
    var box = $("method-boxes");
    R.methods.forEach(function (m) {
      if (m.key.indexOf("Skala") === 0) return;
      var l = document.createElement("label");
      l.innerHTML = '<input type="checkbox" checked value="' + esc(m.key) + '"> <span class="sw" style="background:var(' + B1.METHODS[m.key].color + ')"></span>' + esc(m.label);
      box.appendChild(l);
      l.firstChild.addEventListener("change", function (ev) { S.on[m.key] = ev.target.checked; buildTicks(); draw(); });
    });
    Array.prototype.forEach.call(app.querySelectorAll('input[name="ref"],input[name="basis"],input[name="d3"]'), function (inp) {
      inp.addEventListener("change", function () { S[inp.name] = inp.value; modeChanged(); });
    });
    $("rx-default").addEventListener("click", function () { resetTicks(); draw(); });
    buildHead();
    modeChanged();
    B1.onResize && document.addEventListener("themechange", function () { /* table colours are CSS only */ });
  }

  // ---- table head with sort buttons
  var COLS = [["rx", "Reaction"], ["m", "Method"], ["basis", "Basis"], ["comp", "Computed", 1], ["ref", "Reference", 1], ["err", "Signed error", 1], ["abs", "|Error|", 1], ["note", "Note"]];
  function buildHead() {
    var tr = $("res-table").querySelector("thead tr");
    COLS.forEach(function (c) {
      var th = document.createElement("th");
      th.scope = "col"; if (c[2]) th.className = "num";
      var b = document.createElement("button"); b.type = "button"; b.className = "sort"; b.textContent = c[1];
      b.addEventListener("click", function () {
        if (S.sort.key === c[0]) S.sort.dir = -S.sort.dir; else { S.sort.key = c[0]; S.sort.dir = 1; }
        draw();
      });
      th.appendChild(b); th.setAttribute("data-key", c[0]);
      tr.appendChild(th);
    });
  }

  function modeChanged() {
    // the own-CCSD(T) reference exists only at the main basis
    var own = app.querySelector('input[name="ref"][value="own"]'), large = app.querySelector('input[name="basis"][value="large"]');
    own.disabled = S.basis === "large"; large.disabled = S.ref === "own";
    var d3 = app.querySelectorAll('input[name="d3"]');
    var note = "";
    if (S.ref === "own") note = "Errors against our own CCSD(T)+ΔMP2 barrier at the same geometry (electronic energies, no zero-point correction). H₂CO is not available: its coupled-cluster reference is unusable. The CCSD(T) rows are the reference itself and are hidden.";
    else if (S.basis === "large") note = "def2-QZVPP was run for HCN, NH₃ and H₂CO and for Skala-1.1, Skala-1.1 without D3 and ωB97M-V only.";
    $("mode-note").textContent = note;
    resetTicks();
    draw();
  }

  // ---- rows of the current view
  function viewRows() {
    var out = [];
    R.rows.forEach(function (r) {
      if (r.b !== S.basis) return;
      if (r.m === "Skala-1.1" && S.d3 === "without") return;
      if (r.m === "Skala-1.1-noD3" && S.d3 === "with") return;
      if (r.m.indexOf("Skala") !== 0 && !S.on[r.m]) return;
      var rxn = rxBy[r.r], notes = [], err = r.err, ref = rxn.ref, comp = r.comp;
      if (S.ref === "own") {
        if (CC.indexOf(r.m) >= 0) return;
        err = r.own; ref = rxn.own_ref; comp = r.elec;
        if (!rxn.own_usable) { notes.push("our CCSD(T) is unusable here"); err = null; }
      } else {
        if (r.zpe) notes.push("zero-point " + sgn(r.zpe) + " included");
        if (!rxn.verified) notes.push("reference unverified, in no published statistic");
      }
      if (r.wd) notes.push("withdrawn by the check, in no statistic");
      var usable = !r.wd && err !== null && err !== undefined;
      out.push({ rx: r.r, m: r.m, basis: r.basis, comp: comp, ref: ref, err: err, abs: err === null || err === undefined ? null : Math.abs(err), note: notes.join("; "), usable: usable, wd: r.wd });
    });
    return out;
  }

  function sortRows(rows) {
    var k = S.sort.key, d = S.sort.dir;
    rows.sort(function (a, b) {
      var x, y;
      if (k === "rx") { x = rxOrder[a.rx]; y = rxOrder[b.rx]; if (x === y) { x = mOrder[a.m]; y = mOrder[b.m]; } }
      else if (k === "m") { x = mOrder[a.m]; y = mOrder[b.m]; if (x === y) { x = rxOrder[a.rx]; y = rxOrder[b.rx]; } }
      else if (k === "basis" || k === "note") { x = a[k]; y = b[k]; }
      else { x = a[k]; y = b[k]; if (x === null || x === undefined) return 1; if (y === null || y === undefined) return -1; }
      return (x < y ? -1 : x > y ? 1 : 0) * d;
    });
    return rows;
  }

  function td(txt, cls, title) { return "<td" + (cls ? ' class="' + cls + '"' : "") + (title ? ' title="' + esc(title) + '"' : "") + ">" + txt + "</td>"; }
  function draw() {
    var rows = sortRows(viewRows()), html = "", n = 0;
    rows.forEach(function (r) {
      var sc = r.err === null || r.err === undefined ? "" : (r.err > 0 ? " pos" : r.err < 0 ? " neg" : "");
      html += "<tr" + (r.wd || !rxBy[r.rx].verified ? ' class="flag"' : "") + ">" +
        td(esc(rxBy[r.rx].label), "nowrap") + td('<span class="sw" style="' + (r.m === "Skala-1.1-noD3" ? "background:transparent;border:2px solid var(" + B1.METHODS[r.m].color + ")" : "background:var(" + B1.METHODS[r.m].color + ")") + '"></span>' + esc(label(r.m))) +
        td(esc(r.basis), "nowrap") + td(fmt(r.comp), "num", r.comp === null ? "" : r.comp.toFixed(4)) + td(fmt(r.ref), "num", r.ref === null ? "" : String(r.ref)) +
        td(sgn(r.err), "num" + sc, r.err === null || r.err === undefined ? "" : r.err.toFixed(4)) + td(fmt(r.abs), "num") + td(esc(r.note), "small") + "</tr>";
      n++;
    });
    $("res-table").querySelector("tbody").innerHTML = html;
    Array.prototype.forEach.call($("res-table").querySelectorAll("th"), function (th) {
      var k = th.getAttribute("data-key");
      if (k === S.sort.key) th.setAttribute("aria-sort", S.sort.dir > 0 ? "ascending" : "descending"); else th.removeAttribute("aria-sort");
    });
    $("table-foot").textContent = n + " rows. Energies in kcal/mol. Hover over a number for four decimals. Signed error = computed minus reference. SN2 barriers are measured from the ion-molecule complex. Rows in grey italics are in no published statistic.";
    drawCalc(viewRows());
  }

  // ---- calculator
  function defaultSet() {
    if (S.ref === "own") return ["hcn", "nh3", "ch3nc", "sn2_ff", "sn2_clcl"];
    if (S.basis === "large") return ["hcn", "nh3", "h2co"];
    return R.reactions.filter(function (r) { return r.headline; }).map(function (r) { return r.key; });
  }
  function resetTicks() {
    S.ticked = {};
    defaultSet().forEach(function (k) { S.ticked[k] = true; });
    buildTicks();
  }
  function available() {
    var av = {};
    viewRows().forEach(function (r) { if (r.usable) av[r.rx] = true; });
    return av;
  }
  function buildTicks() {
    var box = $("rx-boxes"), av = available();
    box.innerHTML = "";
    R.reactions.forEach(function (r) {
      var l = document.createElement("label"), dis = !av[r.key];
      l.innerHTML = '<input type="checkbox" value="' + r.key + '"' + (S.ticked[r.key] && !dis ? " checked" : "") + (dis ? " disabled" : "") + "> " + esc(r.label) +
        (!r.verified && S.ref === "lit" ? ' <span class="tag warn" title="literature value could not be verified">reference unverified</span>' : "") + (dis ? ' <span class="small mut">(no data in this view)</span>' : "");
      box.appendChild(l);
      l.firstChild.addEventListener("change", function (ev) { S.ticked[r.key] = ev.target.checked; drawCalc(viewRows()); });
    });
  }

  function setEq(a, b) { if (a.length !== b.length) return false; return a.every(function (x) { return b.indexOf(x) >= 0; }); }
  function tickedList() { return Object.keys(S.ticked).filter(function (k) { return S.ticked[k]; }); }

  function drawCalc(rows) {
    var sel = tickedList(), by = {}, methods = [];
    rows.forEach(function (r) {
      if (!r.usable || sel.indexOf(r.rx) < 0) return;
      (by[r.m] = by[r.m] || []).push(r.err);
    });
    R.methods.forEach(function (m) { if (by[m.key]) methods.push(m.key); });
    var pubSet = S.ref === "own" ? ["hcn", "nh3", "ch3nc", "sn2_ff", "sn2_clcl"] : (S.basis === "large" ? ["hcn", "nh3", "h2co"] : R.reactions.filter(function (r) { return r.headline; }).map(function (r) { return r.key; }));
    var isFull = setEq(sel, pubSet), removed = null;
    if (!isFull && S.ref === "lit" && S.basis === "main" && sel.length === pubSet.length - 1 && sel.every(function (k) { return pubSet.indexOf(k) >= 0; }))
      removed = pubSet.filter(function (k) { return sel.indexOf(k) < 0; })[0];
    var html = "", allMatch = true, anyPub = false;
    methods.forEach(function (m) {
      var e = by[m], n = e.length, sa = 0, ss = 0, mx = 0;
      e.forEach(function (x) { sa += Math.abs(x); ss += x; mx = Math.max(mx, Math.abs(x)); });
      var mae = sa / n, mse = ss / n, pubTxt = "", pub = null;
      if (isFull) {
        if (S.ref === "own") pub = (R.published_internal.filter(function (p) { return p.m === m; })[0]);
        else pub = R.published.filter(function (p) { return p.m === m && p.b === S.basis; })[0];
      } else if (removed) {
        var l = R.leave_one_out.filter(function (p) { return p.m === m && p.removed === removed; })[0];
        if (l) pub = { mae: l.mae, mse: l.mse, max: null, n: l.n };
      }
      if (pub && pub.n === n) {
        anyPub = true;
        var okm = Math.abs(pub.mae - mae) <= 0.0015 && Math.abs(pub.mse - mse) <= 0.0015 && (pub.max === null || Math.abs(pub.max - mx) <= 0.0015);
        if (!okm) allMatch = false;
        pubTxt = '<span class="tag ' + (okm ? "ok" : "bad") + '">' + (okm ? "matches" : "differs") + "</span> MAE " + fmt(pub.mae, 3) + ", MSE " + sgn(pub.mse, 3) + (pub.max === null ? "" : ", max " + fmt(pub.max, 3)) +
          (pub.lo !== undefined && pub.lo !== null && isFull ? " (95 % interval of the MAE " + fmt(pub.lo) + " to " + fmt(pub.hi) + ")" : "");
      } else if (pub) pubTxt = '<span class="small mut">published n = ' + pub.n + ", here n = " + n + "</span>";
      else pubTxt = '<span class="small mut">no published value for this set</span>';
      html += "<tr><td>" + '<span class="sw" style="' + (m === "Skala-1.1-noD3" ? "background:transparent;border:2px solid var(" + B1.METHODS[m].color + ")" : "background:var(" + B1.METHODS[m].color + ")") + '"></span>' + esc(label(m)) + '</td><td class="num">' + n + '</td><td class="num" title="' + mae.toFixed(4) + '">' + fmt(mae) +
        '</td><td class="num ' + (mse > 0 ? "pos" : mse < 0 ? "neg" : "") + '" title="' + mse.toFixed(4) + '">' + sgn(mse) + '</td><td class="num">' + fmt(mx) + "</td><td>" + pubTxt + "</td></tr>";
    });
    $("stat-table").querySelector("tbody").innerHTML = html || '<tr><td colspan="6">Tick at least one reaction.</td></tr>';
    var cap = sel.length + " reaction" + (sel.length === 1 ? "" : "s") + ": " + sel.map(function (k) { return rxBy[k].label; }).join(", ");
    $("stat-cap").textContent = "Mean errors (kcal/mol) for " + cap;
    var msg = "";
    if (isFull) msg = "This is the published set" + (S.ref === "own" ? " of Table 3 (ii)" : (S.basis === "large" ? " at the larger basis" : "")) + ". ";
    else if (removed) msg = "This is the published leave-one-out set without " + rxBy[removed].label + ". ";
    if (anyPub) msg += allMatch ? "All computed values match the published statistics (to 0.0015 kcal/mol)." : "Some values differ from the published statistics: check the ticked reactions.";
    if (sel.indexOf("ch3nc") >= 0 && S.ref === "lit") msg += " CH₃NC is included although its reference is unverified.";
    if (sel.indexOf("h2co") >= 0 && S.ref === "lit" && S.basis === "main") msg += " The CCSD(T) rows of H₂CO were withdrawn by the check and are skipped, so their n is one lower.";
    $("stat-note").textContent = msg.trim();
  }
})();
