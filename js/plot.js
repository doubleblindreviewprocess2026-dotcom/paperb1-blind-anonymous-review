/* Minimal SVG plot helper (no libraries). Colours come from the CSS custom properties of css/style.css, so the plots follow the theme.
   Usage: var p = new B1.Plot(container, {height: 300, x: [0, 10], y: [-5, 5], xlabel: "...", ylabel: "...", label: "text for screen readers"});
          p.line(points, {color: "--c-skala"}); p.marker(x, y, "circle", {color: "--c-skala", filled: true, title: "text"}); p.finish(); */
(function () {
  "use strict";
  var B1 = (window.B1 = window.B1 || {});
  var NS = "http://www.w3.org/2000/svg";

  function el(name, attrs, parent) {
    var e = document.createElementNS(NS, name);
    for (var k in attrs) if (Object.prototype.hasOwnProperty.call(attrs, k)) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }
  function niceStep(span, n) {
    var raw = span / n, p = Math.pow(10, Math.floor(Math.log(raw) / Math.LN10)), f = raw / p;
    return (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * p;
  }
  function ticks(lo, hi, n) {
    var st = niceStep(hi - lo, n || 6), t = [], v = Math.ceil(lo / st - 1e-9) * st;
    for (; v <= hi + 1e-9; v += st) t.push(Math.abs(v) < st * 1e-9 ? 0 : +v.toFixed(10));
    return t;
  }
  function tickText(v, step) {
    var d = step >= 1 ? 0 : Math.min(4, Math.ceil(-Math.log(step) / Math.LN10 - 1e-9));
    return B1.fmt(v, d);
  }

  function Plot(container, opt) {
    this.opt = opt || {};
    var o = this.opt;
    this.container = container;
    this.w = Math.max(280, Math.floor(container.clientWidth || container.getBoundingClientRect().width || 320));
    this.h = o.height || Math.round(Math.min(420, Math.max(240, this.w * 0.62)));
    var narrow = this.w < 480;
    this.m = o.margin || { l: narrow ? 46 : 58, r: narrow ? 10 : 18, t: o.title ? 26 : 12, b: o.xlabel ? (narrow && o.xlabel.length > 38 ? 60 : 46) : 30 };
    container.textContent = "";
    this.svg = el("svg", { viewBox: "0 0 " + this.w + " " + this.h, width: this.w, height: this.h, role: "img", "aria-label": o.label || "plot" }, container);
    this.svg.style.maxWidth = "100%";
    this.gGrid = el("g", { "class": "gridl" }, this.svg);
    this.gBack = el("g", {}, this.svg);
    this.gData = el("g", {}, this.svg);
    this.gAxis = el("g", { "class": "axis" }, this.svg);
    this.gTop = el("g", {}, this.svg);
    var xr = o.x, yr = o.y;
    this.x0 = this.m.l; this.x1 = this.w - this.m.r; this.y0 = this.h - this.m.b; this.y1 = this.m.t;
    this.xd = o.xrev ? [xr[1], xr[0]] : xr; this.yd = yr;
    this.drawAxes();
  }
  Plot.prototype.sx = function (v) { return this.x0 + (v - this.xd[0]) / (this.xd[1] - this.xd[0]) * (this.x1 - this.x0); };
  Plot.prototype.sy = function (v) { return this.y0 + (v - this.yd[0]) / (this.yd[1] - this.yd[0]) * (this.y1 - this.y0); };
  Plot.prototype.drawAxes = function () {
    var o = this.opt, i, t, g;
    var xl = Math.min(o.x[0], o.x[1]), xh = Math.max(o.x[0], o.x[1]);
    var xt = o.xticks || ticks(xl, xh, this.w < 480 ? 4 : 7), yt = o.yticks || ticks(o.y[0], o.y[1], this.h < 300 ? 5 : 7);
    var xs = xt.length > 1 ? Math.abs(xt[1] - xt[0]) : 1, ys = yt.length > 1 ? Math.abs(yt[1] - yt[0]) : 1;
    for (i = 0; i < yt.length; i++) {
      var y = this.sy(yt[i]);
      el("line", { x1: this.x0, x2: this.x1, y1: y, y2: y }, this.gGrid);
      t = el("text", { x: this.x0 - 6, y: y + 4, "text-anchor": "end", "class": "ax" }, this.gAxis); t.textContent = tickText(yt[i], ys);
    }
    for (i = 0; i < xt.length; i++) {
      var x = this.sx(xt[i]);
      el("line", { x1: x, x2: x, y1: this.y0, y2: this.y1 }, this.gGrid);
      t = el("text", { x: x, y: this.y0 + 15, "text-anchor": "middle", "class": "ax" }, this.gAxis); t.textContent = o.xfmt ? o.xfmt(xt[i]) : tickText(xt[i], xs);
    }
    el("line", { x1: this.x0, x2: this.x1, y1: this.y0, y2: this.y0 }, this.gAxis);
    el("line", { x1: this.x0, x2: this.x0, y1: this.y0, y2: this.y1 }, this.gAxis);
    if (o.xlabel) {
      t = el("text", { x: (this.x0 + this.x1) / 2, y: this.h - 6, "text-anchor": "middle", "class": "lbl" }, this.gAxis); t.textContent = o.xlabel;
      var avail = this.w - 16, bw = 0;
      try { bw = t.getComputedTextLength(); } catch (e) { bw = 0; }
      if (bw > avail) {
        // too wide for a narrow screen: break at the first "(" or ";" into two lines
        var cut = o.xlabel.search(/ \(|; /);
        if (cut > 0) {
          t.textContent = o.xlabel.slice(0, cut).trim(); t.setAttribute("y", this.h - 22);
          var t2 = el("text", { x: (this.x0 + this.x1) / 2, y: this.h - 6, "text-anchor": "middle", "class": "lbl" }, this.gAxis);
          t2.textContent = o.xlabel.slice(cut).replace(/^[;\s]+/, "");
        }
      }
    }
    if (o.ylabel) { t = el("text", { transform: "translate(13," + (this.y0 + this.y1) / 2 + ") rotate(-90)", "text-anchor": "middle", "class": "lbl" }, this.gAxis); t.textContent = o.ylabel; }
    if (o.title) { t = el("text", { x: this.x0, y: 16, "class": "lbl", "font-weight": "600" }, this.gAxis); t.textContent = o.title; }
    if (o.zero !== false && o.y[0] < 0 && o.y[1] > 0) el("line", { x1: this.x0, x2: this.x1, y1: this.sy(0), y2: this.sy(0), "class": "zero" }, this.gBack);
  };
  Plot.prototype.band = function (y0, y1, cls) {
    el("rect", { x: this.x0, width: this.x1 - this.x0, y: Math.min(this.sy(y0), this.sy(y1)), height: Math.abs(this.sy(y1) - this.sy(y0)), "class": cls || "band" }, this.gBack);
  };
  Plot.prototype.xband = function (a, b, cls) {
    var p = this.sx(a), q = this.sx(b);
    el("rect", { x: Math.min(p, q), width: Math.abs(q - p), y: this.y1, height: this.y0 - this.y1, "class": cls || "band" }, this.gBack);
  };
  Plot.prototype.hline = function (y, o) {
    o = o || {};
    var e = el("line", { x1: this.x0, x2: this.x1, y1: this.sy(y), y2: this.sy(y) }, this.gData);
    e.style.stroke = "var(" + (o.color || "--mut") + ")"; e.style.strokeWidth = o.width || 1; if (o.dash) e.setAttribute("stroke-dasharray", o.dash);
    return e;
  };
  Plot.prototype.vline = function (x, o) {
    o = o || {};
    var e = el("line", { x1: this.sx(x), x2: this.sx(x), y1: this.y0, y2: this.y1 }, o.top ? this.gTop : this.gData);
    e.style.stroke = "var(" + (o.color || "--mut") + ")"; e.style.strokeWidth = o.width || 1; if (o.dash) e.setAttribute("stroke-dasharray", o.dash);
    return e;
  };
  Plot.prototype.line = function (pts, o) {
    o = o || {};
    var d = "", i, started = false;
    for (i = 0; i < pts.length; i++) {
      if (pts[i] === null || pts[i][1] === null || isNaN(pts[i][1])) { started = false; continue; }
      d += (started ? "L" : "M") + this.sx(pts[i][0]).toFixed(1) + " " + this.sy(pts[i][1]).toFixed(1);
      started = true;
    }
    var e = el("path", { d: d, fill: "none" }, this.gData);
    e.style.stroke = "var(" + (o.color || "--fg") + ")"; e.style.strokeWidth = o.width || 1.8;
    if (o.dash) e.setAttribute("stroke-dasharray", o.dash);
    e.setAttribute("stroke-linejoin", "round");
    return e;
  };
  Plot.prototype.curve = function (fn, a, b, o) {
    var pts = [], n = 80, i;
    for (i = 0; i <= n; i++) { var x = a + (b - a) * i / n; pts.push([x, fn(x)]); }
    return this.line(pts, o);
  };
  Plot.prototype.text = function (x, y, s, o) {
    o = o || {};
    var t = el("text", { x: x, y: y, "text-anchor": o.anchor || "start", "class": o.cls || "lbl" }, o.top ? this.gTop : this.gData);
    if (o.color) t.style.fill = "var(" + o.color + ")";
    t.textContent = s;
    return t;
  };
  // marker shapes: circle, ring, square, triangle, diamond, cross (x), plus, star
  Plot.prototype.marker = function (x, y, shape, o) {
    o = o || {};
    if (isNaN(x) || isNaN(y) || y === null) return null;
    var cx = this.sx(x), cy = this.sy(y), r = o.size || 5, col = "var(" + (o.color || "--fg") + ")", e, d;
    var filled = o.filled !== false && shape !== "ring" && shape !== "cross" && shape !== "plus" && shape !== "x";
    var g = el("g", { "class": "pt" }, o.top ? this.gTop : this.gData);
    switch (shape) {
      case "square": e = el("rect", { x: cx - r, y: cy - r, width: 2 * r, height: 2 * r }, g); break;
      case "triangle": e = el("path", { d: "M" + cx + " " + (cy - r * 1.15) + "L" + (cx + r * 1.1) + " " + (cy + r * 0.85) + "L" + (cx - r * 1.1) + " " + (cy + r * 0.85) + "Z" }, g); break;
      case "diamond": e = el("path", { d: "M" + cx + " " + (cy - r * 1.3) + "L" + (cx + r * 1.3) + " " + cy + "L" + cx + " " + (cy + r * 1.3) + "L" + (cx - r * 1.3) + " " + cy + "Z" }, g); break;
      case "cross": case "x": e = el("path", { d: "M" + (cx - r) + " " + (cy - r) + "L" + (cx + r) + " " + (cy + r) + "M" + (cx - r) + " " + (cy + r) + "L" + (cx + r) + " " + (cy - r) }, g); break;
      case "plus": e = el("path", { d: "M" + (cx - r * 1.2) + " " + cy + "L" + (cx + r * 1.2) + " " + cy + "M" + cx + " " + (cy - r * 1.2) + "L" + cx + " " + (cy + r * 1.2) }, g); break;
      case "star": d = ""; for (var k = 0; k < 10; k++) { var rr = k % 2 ? r * 0.5 : r * 1.3, a = -Math.PI / 2 + k * Math.PI / 5; d += (k ? "L" : "M") + (cx + rr * Math.cos(a)).toFixed(1) + " " + (cy + rr * Math.sin(a)).toFixed(1); } e = el("path", { d: d + "Z" }, g); break;
      default: e = el("circle", { cx: cx, cy: cy, r: r }, g);
    }
    if (filled) { e.style.fill = col; e.style.stroke = o.edge ? "var(--bg)" : col; e.style.strokeWidth = o.edge ? 1 : 1; }
    else { e.style.fill = shape === "ring" || o.filled === false ? "var(--card)" : "none"; e.style.stroke = col; e.style.strokeWidth = 2; }
    if (o.title) { var t = el("title", {}, g); t.textContent = o.title; }
    return g;
  };
  Plot.prototype.finish = function () { return this.svg; };

  B1.Plot = Plot;
  B1.plotTicks = ticks;
  // redraw a chart when its container changes width (rotation, window resize)
  B1.onResize = function (container, fn) {
    var last = container.clientWidth, timer = null;
    var handler = function () {
      if (timer) clearTimeout(timer);
      timer = setTimeout(function () { if (container.clientWidth !== last) { last = container.clientWidth; fn(); } }, 120);
    };
    window.addEventListener("resize", handler);
    document.addEventListener("themechange", fn);
  };
  B1.legendHTML = function (items) {
    return items.map(function (it) {
      var m = B1.METHODS[it] || { label: it, color: "--fg" };
      return '<span><span class="sw" style="background:var(' + m.color + ')"></span>' + B1.esc(m.label) + "</span>";
    }).join("");
  };
})();
