/* Checks page: filter the full report by status and text. */
(function () {
  "use strict";
  var tb = document.querySelector("#check-table tbody");
  if (!tb) return;
  var rows = Array.prototype.slice.call(tb.rows), pills = document.getElementById("status-pills"), box = document.getElementById("checksearch"), count = document.getElementById("check-count");
  var statuses = ["All"], tally = {};
  rows.forEach(function (r) { var s = r.getAttribute("data-status"); if (!tally[s]) { tally[s] = 0; statuses.push(s); } tally[s]++; });
  var cur = "All";
  statuses.forEach(function (s) {
    var b = document.createElement("button"); b.type = "button"; b.textContent = s + (s === "All" ? " (" + rows.length + ")" : " (" + tally[s] + ")");
    b.setAttribute("aria-pressed", s === "All" ? "true" : "false"); b.setAttribute("data-s", s);
    b.addEventListener("click", function () { cur = s; Array.prototype.forEach.call(pills.children, function (c) { c.setAttribute("aria-pressed", c.getAttribute("data-s") === s ? "true" : "false"); }); apply(); });
    pills.appendChild(b);
  });
  function apply() {
    var q = box.value.trim().toLowerCase(), n = 0;
    rows.forEach(function (r) {
      var ok = (cur === "All" || r.getAttribute("data-status") === cur) && (!q || r.textContent.toLowerCase().indexOf(q) >= 0);
      r.hidden = !ok; if (ok) n++;
    });
    count.textContent = n + " of " + rows.length + " lines shown.";
  }
  box.addEventListener("input", apply);
  apply();
})();
