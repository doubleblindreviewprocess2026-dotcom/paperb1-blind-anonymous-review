/* Shared helpers: theme switch, configuration (js/config.js), Colab buttons, citation box, JSON loading, number formats. No libraries. */
(function () {
  "use strict";
  var CFG = window.SITE_CONFIG || {};
  var B1 = (window.B1 = window.B1 || {});
  B1.cfg = CFG;

  // site root URL, from the address of this script (works under any host and sub-folder)
  var me = document.currentScript && document.currentScript.src;
  B1.rootURL = me ? me.replace(/js\/common\.js(\?.*)?$/, "") : location.href.replace(/[^/]*$/, "");

  function isPlaceholder(v) { return !v || /PLACEHOLDER/.test(v); }
  B1.isPlaceholder = isPlaceholder;

  // robots: only when anonymous: true
  if (CFG.anonymous === true) {
    var m = document.createElement("meta"); m.name = "robots"; m.content = "noindex"; document.head.appendChild(m);
  }

  // ---------- theme
  function currentTheme() {
    var t = document.documentElement.getAttribute("data-theme");
    if (t) return t;
    return window.matchMedia && matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
  }
  function setThemeLabel(btn) { btn.textContent = currentTheme() === "dark" ? "Light theme" : "Dark theme"; btn.setAttribute("aria-label", "Switch to the " + (currentTheme() === "dark" ? "light" : "dark") + " theme"); }
  function initTheme() {
    var btn = document.querySelector("button.theme");
    if (!btn) return;
    setThemeLabel(btn);
    btn.addEventListener("click", function () {
      var next = currentTheme() === "dark" ? "light" : "dark";
      document.documentElement.setAttribute("data-theme", next);
      try { localStorage.setItem("b1-theme", next); } catch (e) { /* storage blocked: fine */ }
      setThemeLabel(btn);
      document.dispatchEvent(new CustomEvent("themechange"));
    });
    if (window.matchMedia) {
      var mq = matchMedia("(prefers-color-scheme: dark)");
      var h = function () { setThemeLabel(btn); document.dispatchEvent(new CustomEvent("themechange")); };
      if (mq.addEventListener) mq.addEventListener("change", h);
    }
  }

  // ---------- configuration
  B1.colabURL = function (path) {
    var repo = CFG.repo || "REPO_PLACEHOLDER";
    return "https://colab.research.google.com/github/" + repo + "/blob/" + (CFG.branch || "main") + "/" + (CFG.folder || "") + path;
  };
  B1.githubURL = function (path) {
    return "https://github.com/" + (CFG.repo || "REPO_PLACEHOLDER") + "/blob/" + (CFG.branch || "main") + "/" + (CFG.folder || "") + (path || "");
  };
  B1.citation = function () {
    var anon = CFG.anonymous === true;
    var title = CFG.paperTitle || "";
    var lines = [];
    lines.push("Companion website and data for: " + title + ".");
    if (anon) {
      lines.push("Anonymous copy for review. Website: " + B1.rootURL);
    } else {
      lines.push((CFG.holder && !isPlaceholder(CFG.holder) ? CFG.holder + ". " : "") + (CFG.year || "") + ". Website: " + B1.rootURL);
      lines.push("Archive (same files): " + (isPlaceholder(CFG.doi) ? "https://doi.org/DOI_PLACEHOLDER" : "https://doi.org/" + CFG.doi));
      lines.push("Code and data repository: " + (isPlaceholder(CFG.repo) ? "https://github.com/REPO_PLACEHOLDER" : "https://github.com/" + CFG.repo));
    }
    return lines.join("\n");
  };
  function applyConfig() {
    var anon = CFG.anonymous === true;
    var set = function (sel, fn) { Array.prototype.forEach.call(document.querySelectorAll(sel), fn); };
    set("[data-cfg]", function (el) {
      var k = el.getAttribute("data-cfg"), v = CFG[k];
      if (k === "doi" && v && !isPlaceholder(v)) { el.textContent = v; if (el.tagName === "A") el.href = "https://doi.org/" + v; return; }
      if (k === "repo" && v && !isPlaceholder(v) && el.tagName === "A") { el.href = "https://github.com/" + v; }
      el.textContent = v == null ? "" : v;
      if (isPlaceholder(v)) el.classList.add("placeholder");
    });
    set("[data-hide-anon]", function (el) { if (anon) el.classList.add("hide"); });
    set("[data-show-anon]", function (el) { if (!anon) el.classList.add("hide"); });
    set("a[data-colab]", function (a) {
      a.href = B1.colabURL(a.getAttribute("data-colab"));
      a.target = "_blank"; a.rel = "noopener";
      if (isPlaceholder(CFG.repo)) { a.classList.add("pending"); a.title = "Set repo in js/config.js to make this button work"; }
    });
    set("a[data-github]", function (a) {
      a.href = B1.githubURL(a.getAttribute("data-github")); a.target = "_blank"; a.rel = "noopener";
    });
    set("[data-cfg-warning]", function (el) {
      var bad = !anon && (isPlaceholder(CFG.repo) || isPlaceholder(CFG.doi) || isPlaceholder(CFG.contact) || isPlaceholder(CFG.holder));
      el.classList.toggle("hide", !bad);
    });
    set("[data-citation]", function (el) { el.textContent = B1.citation(); });
    set("button[data-copy]", function (b) {
      b.addEventListener("click", function () {
        var t = document.querySelector(b.getAttribute("data-copy"));
        if (!t) return;
        var txt = t.textContent;
        var done = function () { var o = b.textContent; b.textContent = "Copied"; setTimeout(function () { b.textContent = o; }, 1500); };
        if (navigator.clipboard && navigator.clipboard.writeText) { navigator.clipboard.writeText(txt).then(done, function () { fallbackCopy(txt, done); }); } else fallbackCopy(txt, done);
      });
    });
  }
  function fallbackCopy(txt, done) {
    var ta = document.createElement("textarea"); ta.value = txt; ta.setAttribute("readonly", ""); ta.style.position = "fixed"; ta.style.opacity = "0";
    document.body.appendChild(ta); ta.select(); try { document.execCommand("copy"); done(); } catch (e) { /* ignore */ } ta.remove();
  }

  // ---------- data loading
  B1.loadJSON = function (url) {
    return fetch(url).then(function (r) {
      if (!r.ok) throw new Error(url + ": HTTP " + r.status);
      return r.json();
    });
  };
  B1.showLoadError = function (el, err) {
    if (!el) return;
    var local = location.protocol === "file:";
    el.innerHTML = "";
    var d = document.createElement("div"); d.className = "note bad"; d.setAttribute("role", "alert");
    d.textContent = local
      ? "This page loads its data files with fetch(), which browsers block for files opened from disk. Serve the folder (python -m http.server in the site folder, then open http://localhost:8000/) or open the hosted copy."
      : "The data file could not be loaded (" + (err && err.message ? err.message : err) + ").";
    el.appendChild(d);
  };

  // ---------- number formats (U+2212 for the minus sign)
  B1.fmt = function (x, d) { if (x === null || x === undefined || x === "" || isNaN(x)) return "–"; var s = Number(x).toFixed(d === undefined ? 2 : d); return s.replace("-", "−"); };
  B1.sgn = function (x, d) { if (x === null || x === undefined || x === "" || isNaN(x)) return "–"; var v = Number(x); var s = v.toFixed(d === undefined ? 2 : d); if (Number(s) === 0) return (s.replace("-", "")); return (v > 0 ? "+" : "−") + s.replace("-", ""); };
  B1.esc = function (s) { return String(s).replace(/[&<>"]/g, function (c) { return ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", "\"": "&quot;" })[c]; }); };
  B1.css = function (v) { return getComputedStyle(document.documentElement).getPropertyValue(v).trim(); };
  B1.median = function (a) { if (!a.length) return NaN; var s = a.slice().sort(function (x, y) { return x - y; }), m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; };

  // method labels and colours shared by the pages
  B1.METHODS = {
    "Skala-1.1": { label: "Skala-1.1", color: "--c-skala", marker: "circle" },
    "Skala-1.1-noD3": { label: "Skala-1.1 without D3", color: "--c-skala", marker: "ring" },
    "wB97M-V": { label: "ωB97M-V", color: "--c-wb97", marker: "square" },
    "B3LYP-D3BJ": { label: "B3LYP-D3(BJ)", color: "--c-b3lyp", marker: "triangle" },
    "r2SCAN": { label: "r²SCAN", color: "--c-r2scan", marker: "diamond" },
    "PBE": { label: "PBE", color: "--c-pbe", marker: "cross" },
    "CCSD(T)-fc": { label: "CCSD(T), frozen core", color: "--c-cc", marker: "x" },
    "CCSD(T)+dMP2": { label: "CCSD(T)+ΔMP2", color: "--c-cc", marker: "plus" },
    "CCSD(T) raw": { label: "CCSD(T), frozen core", color: "--c-cc", marker: "x" },
    "GFN2-xTB": { label: "GFN2-xTB", color: "--c-xtb", marker: "star" }
  };
  B1.REACTION_LABEL = { hcn: "HCN → HNC", nh3: "NH₃ inversion", h2co: "H₂CO → H₂ + CO", ch3nc: "CH₃NC → CH₃CN", sn2_ff: "F⁻ + CH₃F", sn2_clcl: "Cl⁻ + CH₃Cl" };

  document.addEventListener("DOMContentLoaded", function () { initTheme(); applyConfig(); });
})();
