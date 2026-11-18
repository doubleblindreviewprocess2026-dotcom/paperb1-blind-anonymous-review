/* Data page: filter the file list by file name, title or column name. */
(function () {
  "use strict";
  var box = document.getElementById("filesearch");
  if (!box) return;
  var items = Array.prototype.slice.call(document.querySelectorAll("#filelist details.filerow")), count = document.getElementById("filecount");
  function apply() {
    var q = box.value.trim().toLowerCase(), n = 0;
    items.forEach(function (d) {
      var ok = !q || d.getAttribute("data-search").indexOf(q) >= 0;
      d.hidden = !ok; if (ok) n++;
      if (q && ok && items.length - n > 0 && n <= 3) d.open = true;
      if (!q) d.open = false;
    });
    count.textContent = n + " of " + items.length + " files shown.";
  }
  box.addEventListener("input", apply);
  apply();
})();
