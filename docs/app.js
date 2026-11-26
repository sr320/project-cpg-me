/* CpG O/E Atlas: renders docs/data/dashboard.json (built by analysis/build_dashboard.py). */
(() => {
  "use strict";
  const SVG = "http://www.w3.org/2000/svg";
  const FEATURE_LABEL = { cds: "CDS", gene_body: "Gene body", upstream: "Upstream 1 kb" };
  const state = { data: null, species: null, feature: "gene_body", goslimFeature: "gene_body",
                  sortKey: null, sortDir: 1 };
  const $ = (id) => document.getElementById(id);
  const tip = $("tooltip");

  // ---------- helpers ----------
  const el = (tag, attrs = {}, text) => {
    const e = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
    if (text !== undefined) e.textContent = text;
    return e;
  };
  const sv = (tag, attrs = {}, text) => {
    const e = document.createElementNS(SVG, tag);
    for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
    if (text !== undefined) e.textContent = text;
    return e;
  };
  const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  const fmt = (x, d = 2) => (x === null || x === undefined || Number.isNaN(x) ? "–" : Number(x).toFixed(d));
  const fmtInt = (x) => (x === null || x === undefined ? "–" : Number(x).toLocaleString("en-US"));
  const fmtP = (p) => {
    if (p === null || p === undefined) return "–";
    if (p === 0) return "< 1e-300";
    if (p < 0.001) { const e = Math.floor(Math.log10(p)); return `${(p / 10 ** e).toFixed(1)}e${e}`; }
    return p.toFixed(p < 0.01 ? 3 : 2);
  };
  const pdf = (x, mu, sd, w) => (w * Math.exp(-0.5 * ((x - mu) / sd) ** 2)) / (sd * Math.sqrt(2 * Math.PI));
  const bestFeature = (sp) => (sp.fits.gene_body ? "gene_body" : "cds");

  function showTip(evt, rows) {
    tip.replaceChildren();
    for (const r of rows) {
      const d = el("div", { class: "row" });
      if (r.key) { const k = el("span", { class: "key-line" }); k.style.borderTopColor = r.key; d.append(k); }
      if (r.v !== undefined) d.append(el("span", { class: "v" }, r.v));
      if (r.l !== undefined) d.append(el("span", { class: "l" }, r.l));
      tip.append(d);
    }
    tip.hidden = false;
    moveTip(evt);
  }
  function moveTip(evt) {
    const pad = 14, w = tip.offsetWidth, h = tip.offsetHeight;
    let x = evt.clientX + pad, y = evt.clientY + pad;
    if (x + w > window.innerWidth - 8) x = evt.clientX - w - pad;
    if (y + h > window.innerHeight - 8) y = evt.clientY - h - pad;
    tip.style.left = `${Math.max(8, x)}px`;
    tip.style.top = `${Math.max(8, y)}px`;
  }
  const hideTip = () => { tip.hidden = true; };
  function hoverable(node, rows) {
    node.setAttribute("tabindex", "0");
    node.addEventListener("pointerenter", (e) => showTip(e, rows()));
    node.addEventListener("pointermove", moveTip);
    node.addEventListener("pointerleave", hideTip);
    node.addEventListener("focus", () => {
      const b = node.getBoundingClientRect();
      showTip({ clientX: b.right, clientY: b.top }, rows());
    });
    node.addEventListener("blur", hideTip);
  }

  function segmented(container, options, current, onChange) {
    container.replaceChildren();
    for (const [value, label, enabled] of options) {
      const b = el("button", { type: "button", role: "radio", "aria-checked": String(value === current) }, label);
      if (!enabled) b.disabled = true;
      b.addEventListener("click", () => onChange(value));
      container.append(b);
    }
  }

  // ---------- header, tiles ----------
  function renderMeta(d) {
    $("meta").textContent = `${d.species.length} species · data built ${d.generated} · commit ${d.git_commit}`;
    $("foot").textContent = `Built from results/ by analysis/build_dashboard.py. Numbers are summaries; per-gene tables are in the repository's results directories (not published here).`;
  }

  function renderTiles(d) {
    const focal = d.species.filter((s) => s.role === "focal");
    const withGB = focal.filter((s) => s.fits.gene_body);
    const bimodalGB = withGB.filter((s) => s.fits.gene_body.bimodal);
    const tested = focal.filter((s) => s.goslim_test);
    const supports = tested.filter((s) => s.goslim_test[bestFeature(s)]?.supports_2015);
    const nVal = d.species.reduce((a, s) => a + (s.validation?.length || 0), 0);
    const ctl = d.species.find((s) => s.role === "control_unmethylated");
    const tiles = [
      ["Species analysed", `${d.species.length}`, `${focal.length} focal, ${d.species.filter((s) => s.role.startsWith("control")).length} controls` +
        (d.species.some((s) => s.role === "tier_check") ? `, ${d.species.filter((s) => s.role === "tier_check").length} tier check` : "")],
      ["Bimodal in gene bodies", `${bimodalGB.length}`, `of ${withGB.length} focal species`, `/ ${withGB.length}`],
      ["Housekeeping genes lower", `${supports.length}`, ctl && ctl.goslim_test ? `of ${tested.length} focal; ${ctl.organism} control: ${ctl.goslim_test[bestFeature(ctl)]?.supports_2015 ? "yes" : "no gap"}` : `of ${tested.length} focal species`, `/ ${tested.length}`],
      ["WGBS validation datasets", `${nVal}`, nVal ? "gene methylation vs CpG O/E" : "alignments running"],
    ];
    const box = $("tiles");
    box.replaceChildren();
    for (const [label, value, sub, of] of tiles) {
      const t = el("div", { class: "tile" });
      t.append(el("div", { class: "label" }, label));
      const v = el("div", { class: "value" }, value);
      if (of) v.append(" ", el("small", {}, of));
      t.append(v, el("div", { class: "sub" }, sub));
      box.append(t);
    }
  }

  // ---------- species table ----------
  const ROLE_LABEL = { control_methylated: "methylated control", control_unmethylated: "unmethylated control",
    tier_check: "transcriptome check" };
  const COLS = [
    { key: "group", label: "Group", get: (s) => s.group, cls: "grp" },
    { key: "organism", label: "Species", get: (s) => s.organism, cls: "sp" },
    { key: "genes", label: "Genes", get: (s) => s.fits.cds?.n, num: true, show: fmtInt },
    { key: "median", label: "Median CpG O/E (CDS)", get: (s) => s.fits.cds?.median, num: true, show: (x) => fmt(x) },
    { key: "d_cds", label: "Ashman's D, CDS", get: (s) => s.fits.cds?.ashman_d, num: true, show: (x) => fmt(x) },
    { key: "d_gb", label: "Ashman's D, gene body", get: (s) => s.fits.gene_body?.ashman_d, num: true, show: (x) => fmt(x) },
    { key: "bimodal", label: "Bimodal (gene body)", get: (s) => (s.fits.gene_body ? (s.fits.gene_body.bimodal ? 1 : 0) : -1),
      chip: (s) => (s.fits.gene_body ? (s.fits.gene_body.bimodal ? ["yes", "yes"] : ["no", ""]) : ["–", ""]) },
    { key: "low", label: "Low class share (gene body)", get: (s) => s.fits.gene_body?.low_weight ?? s.fits.cds?.low_weight, num: true, show: (x) => (x == null ? "–" : `${Math.round(x * 100)} %`) },
    { key: "hk", label: "Median, HK vs IND", get: (s) => s.goslim_test?.[bestFeature(s)]?.median_housekeeping, num: true,
      show: (_, s) => { const t = s.goslim_test?.[bestFeature(s)]; return t ? `${fmt(t.median_housekeeping)} vs ${fmt(t.median_inducible)}` : "–"; } },
    { key: "rb", label: "Rank-biserial", get: (s) => s.goslim_test?.[bestFeature(s)]?.rank_biserial, num: true, show: (x) => (x == null ? "–" : (x > 0 ? "+" : "") + fmt(x)) },
    { key: "p", label: "p (HK lower)", get: (s) => s.goslim_test?.[bestFeature(s)]?.p_housekeeping_lower, num: true, show: fmtP },
    { key: "go", label: "GO coverage", get: (s) => (s.annotation_stats ? s.annotation_stats.n_with_go / s.annotation_stats.n_genes : null), num: true, show: (x) => (x == null ? "–" : `${Math.round(x * 100)} %`) },
  ];

  function renderTable() {
    const d = state.data, table = $("species-table");
    table.replaceChildren();
    const thead = el("thead"), hr = el("tr");
    for (const c of COLS) {
      const th = el("th", { scope: "col", class: c.num ? "num" : "", tabindex: "0" }, c.label);
      if (state.sortKey === c.key) th.setAttribute("aria-sort", state.sortDir > 0 ? "ascending" : "descending");
      const sort = () => {
        state.sortDir = state.sortKey === c.key ? -state.sortDir : 1;
        state.sortKey = c.key;
        renderTable();
      };
      th.addEventListener("click", sort);
      th.addEventListener("keydown", (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); sort(); } });
      hr.append(th);
    }
    thead.append(hr);
    let rows = [...d.species];
    if (state.sortKey) {
      const c = COLS.find((x) => x.key === state.sortKey);
      rows.sort((a, b) => {
        const va = c.get(a), vb = c.get(b);
        if (va == null && vb == null) return 0;
        if (va == null) return 1;
        if (vb == null) return -1;
        return (va < vb ? -1 : va > vb ? 1 : 0) * state.sortDir;
      });
    }
    const tbody = el("tbody");
    for (const s of rows) {
      const tr = el("tr", { tabindex: "0" });
      if (s.accession === state.species) tr.classList.add("selected");
      for (const c of COLS) {
        const td = el("td", { class: [c.num ? "num" : "", c.cls || ""].join(" ").trim() });
        if (c.chip) {
          const [txt, k] = c.chip(s);
          td.append(el("span", { class: `chip ${k}` }, txt));
        } else {
          const v = c.get(s);
          td.textContent = c.show ? c.show(v, s) : v ?? "–";
        }
        if (c.key === "organism" && s.role !== "focal") {
          td.append(" ", el("span", { class: "chip ctl" }, ROLE_LABEL[s.role] || s.role));
        }
        if (c.key === "organism" && s.tier === "tsa") {
          td.append(" ", el("span", { class: "chip ctl" }, "transcriptome, CDS only"));
        }
        tr.append(td);
      }
      const pick = () => selectSpecies(s.accession, true);
      tr.addEventListener("click", pick);
      tr.addEventListener("keydown", (e) => { if (e.key === "Enter") pick(); });
      tbody.append(tr);
    }
    table.append(thead, tbody);
  }

  // ---------- distribution chart ----------
  function selectSpecies(acc, scroll) {
    state.species = acc;
    const sp = state.data.species.find((s) => s.accession === acc);
    if (!sp.hist[state.feature]) state.feature = bestFeature(sp);
    $("species-select").value = acc;
    renderTable();
    renderDist();
    if (scroll) $("h-dist").scrollIntoView({ behavior: "smooth", block: "start" });
  }

  function renderDist() {
    const d = state.data, sp = d.species.find((s) => s.accession === state.species);
    segmented($("feature-seg"), ["cds", "gene_body", "upstream"].map((f) => [f, FEATURE_LABEL[f], !!sp.hist[f]]),
      state.feature, (f) => { state.feature = f; renderDist(); });
    const fit = sp.fits[state.feature], counts = sp.hist[state.feature], bins = d.bins;
    const n = counts.reduce((a, b) => a + b, 0), bw = bins[1] - bins[0];
    $("dist-title").textContent = `${sp.organism} · ${FEATURE_LABEL[state.feature]} · n = ${fmtInt(fit.n)} genes`;
    $("dist-title").classList.add("sp-title");

    const box = $("dist-chart"), W = Math.max(300, box.clientWidth), H = 260;
    const m = { l: 12, r: 12, t: 12, b: 30 }, iw = W - m.l - m.r, ih = H - m.t - m.b;
    const dens = counts.map((c) => c / (n * bw));
    const grid = Array.from({ length: 241 }, (_, i) => i / 120);
    const curves = fit.low_mean != null ? [
      [grid.map((x) => pdf(x, fit.low_mean, fit.low_sd, fit.low_weight)), css("--series-2"), "low"],
      [grid.map((x) => pdf(x, fit.high_mean, fit.high_sd, fit.high_weight)), css("--series-1"), "high"],
    ] : [];
    const ymax = Math.max(...dens, ...curves.flatMap((c) => c[0])) * 1.05;
    const X = (x) => m.l + (x / 2) * iw, Y = (y) => m.t + ih - (y / ymax) * ih;

    const svg = sv("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
      "aria-label": `Histogram of ${FEATURE_LABEL[state.feature]} CpG O/E for ${sp.organism}` });
    svg.append(sv("line", { x1: m.l, x2: W - m.r, y1: m.t + ih, y2: m.t + ih, stroke: css("--axis") }));
    for (const t of [0, 0.5, 1, 1.5, 2]) {
      svg.append(sv("line", { x1: X(t), x2: X(t), y1: m.t + ih, y2: m.t + ih + 4, stroke: css("--axis") }));
      svg.append(sv("text", { x: X(t), y: H - 12, "text-anchor": t === 0 ? "start" : t === 2 ? "end" : "middle", class: "tick-label" }, t.toFixed(1)));
    }
    svg.append(sv("text", { x: W - m.r, y: H - 0, "text-anchor": "end" }, "CpG O/E"));
    const barW = iw / counts.length;
    counts.forEach((c, i) => {
      const y = Y(dens[i]);
      const g = sv("g");
      const h = m.t + ih - y;
      const bar = sv("rect", { x: X(bins[i]) + 0.5, y, width: Math.max(0, barW - 1), height: Math.max(0, h), fill: css("--hist"), rx: Math.min(2, barW / 3) });
      const hit = sv("rect", { x: X(bins[i]), y: m.t, width: barW, height: ih, fill: "transparent" });
      g.append(bar, hit);
      hoverable(hit, () => [{ v: fmtInt(c), l: "genes" }, { l: `CpG O/E ${bins[i].toFixed(3)}–${bins[i + 1].toFixed(3)}` }]);
      hit.addEventListener("pointerenter", () => bar.setAttribute("fill", css("--hist-hover")));
      hit.addEventListener("pointerleave", () => bar.setAttribute("fill", css("--hist")));
      svg.append(g);
    });
    for (const [vals, color] of curves) {
      const path = vals.map((y, i) => `${i ? "L" : "M"}${X(grid[i]).toFixed(1)},${Y(y).toFixed(1)}`).join("");
      svg.append(sv("path", { d: path, fill: "none", stroke: color, "stroke-width": 2,
        "stroke-dasharray": fit.bimodal ? "" : "2 4", "stroke-linecap": "round", "pointer-events": "none" }));
    }
    if (fit.bimodal && fit.threshold != null) {
      svg.append(sv("line", { x1: X(fit.threshold), x2: X(fit.threshold), y1: m.t, y2: m.t + ih,
        stroke: css("--ink-2"), "stroke-dasharray": "4 3", "pointer-events": "none" }));
    }
    box.replaceChildren(svg);

    const lg = $("dist-legend");
    lg.replaceChildren();
    const key = (cls, color, text) => {
      const s = el("span"), k = el("span", { class: cls });
      if (cls === "key-line") k.style.borderTopColor = color;
      if (!fit.bimodal && cls === "key-line") k.style.borderTopStyle = "dotted";
      s.append(k, text);
      lg.append(s);
    };
    if (curves.length) {
      key("key-line", css("--series-2"), `low CpG O/E component (${Math.round(fit.low_weight * 100)} % at ${fmt(fit.low_mean)})`);
      key("key-line", css("--series-1"), `high component (${fmt(fit.high_mean)})`);
      if (fit.bimodal && fit.threshold != null) key("key-dash", null, `threshold ${fmt(fit.threshold)}`);
    }
    const verdict = fit.bimodal ? "Bimodal" : "Not bimodal";
    $("dist-note").textContent = `${verdict}: ΔBIC ${fmt(fit.delta_bic, 0)}, Ashman's D ${fmt(fit.ashman_d)} ` +
      `(rule: ΔBIC > 10 and D > 2). Median ${fmt(fit.median)}; GpC O/E median ${fmt(fit.gpc_oe_median)}.` +
      (fit.bimodal ? "" : " Dotted curves: the two-component fit, shown for reference.");
  }

  // ---------- effect-size dot chart ----------
  function renderEffect() {
    const sp = state.data.species.filter((s) => s.goslim_test);
    const box = $("effect-chart"), W = Math.max(300, box.clientWidth);
    const rowH = 24, labelW = Math.min(190, W * 0.42);
    const m = { l: labelW, r: 16, t: 8, b: 30 }, H = m.t + m.b + rowH * sp.length, iw = W - m.l - m.r;
    const vals = sp.flatMap((s) => ["cds", "gene_body"].map((f) => s.goslim_test[f]?.rank_biserial)).filter((v) => v != null);
    const lo = Math.min(-0.5, Math.floor(Math.min(...vals) * 10) / 10), hi = Math.max(0.2, Math.ceil(Math.max(...vals) * 10) / 10);
    const X = (x) => m.l + ((x - lo) / (hi - lo)) * iw;
    const svg = sv("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Rank-biserial effect per species, CDS and gene body" });
    for (let t = Math.ceil(lo * 10) / 10; t <= hi + 1e-9; t += 0.1) {
      const tt = Math.round(t * 10) / 10;
      svg.append(sv("line", { x1: X(tt), x2: X(tt), y1: m.t, y2: H - m.b, stroke: tt === 0 ? css("--axis") : css("--grid"), "stroke-width": tt === 0 ? 1.5 : 1 }));
      if (Math.abs(tt * 10) % 2 === 0) svg.append(sv("text", { x: X(tt), y: H - 12, "text-anchor": "middle", class: "tick-label" }, (tt > 0 ? "+" : "") + tt.toFixed(1)));
    }
    svg.append(sv("text", { x: X(lo), y: H, "text-anchor": "start" }, "← housekeeping lower"));
    svg.append(sv("text", { x: X(hi), y: H, "text-anchor": "end" }, "higher →"));
    const series = [["gene_body", css("--series-1"), "Gene body"], ["cds", css("--series-2"), "CDS"]];
    sp.forEach((s, i) => {
      const y = m.t + rowH * i + rowH / 2;
      const lab = sv("text", { x: m.l - 8, y: y + 4, "text-anchor": "end", class: "ink2" });
      lab.append(sv("tspan", { "font-style": "italic" }, s.organism));
      svg.append(lab);
      svg.append(sv("line", { x1: m.l, x2: W - m.r, y1: y, y2: y, stroke: css("--grid"), "stroke-dasharray": "1 3" }));
      const pts = series.map(([f]) => s.goslim_test[f]?.rank_biserial).filter((v) => v != null);
      if (pts.length === 2) svg.append(sv("line", { x1: X(pts[0]), x2: X(pts[1]), y1: y, y2: y, stroke: css("--axis"), "stroke-width": 2 }));
      for (const [f, color, label] of series) {
        const t = s.goslim_test[f];
        if (!t) continue;
        const c = sv("circle", { cx: X(t.rank_biserial), cy: y, r: 5, fill: color, stroke: css("--surface"), "stroke-width": 2 });
        const hit = sv("circle", { cx: X(t.rank_biserial), cy: y, r: 12, fill: "transparent" });
        hoverable(hit, () => [
          { l: `${s.organism} · ${label}` },
          { key: color, v: (t.rank_biserial > 0 ? "+" : "") + fmt(t.rank_biserial), l: "rank-biserial" },
          { v: `${fmt(t.median_housekeeping)} vs ${fmt(t.median_inducible)}`, l: "median HK vs IND" },
          { v: fmtP(t.p_housekeeping_lower), l: `p · n ${fmtInt(t.n_housekeeping)} / ${fmtInt(t.n_inducible)}` },
        ]);
        svg.append(c, hit);
      }
    });
    box.replaceChildren(svg);
    const lg = $("effect-legend");
    lg.replaceChildren();
    for (const [, color, label] of series) {
      const s = el("span"), k = el("span", { class: "key-dot" });
      k.style.background = color;
      s.append(k, label);
      lg.append(s);
    }
  }

  // ---------- GO Slim heatmap ----------
  function hexToRgb(h) { const n = parseInt(h.slice(1), 16); return [(n >> 16) & 255, (n >> 8) & 255, n & 255]; }
  function mix(a, b, t) { return a.map((v, i) => Math.round(v + (b[i] - v) * t)); }
  function divergingScale() {
    const dark = css("--surface").toLowerCase() === "#1a1a19";
    const stops = (dark ? ["#86b6ef", "#2a5f9e", css("--mid"), "#9e3f3e", "#f08a88"]
                        : ["#1c5cab", "#86b6ef", css("--mid"), "#f0a3a0", "#c23b3a"]).map(hexToRgb);
    return (v) => {
      const t = Math.max(-1, Math.min(1, v)) * 2 + 2, i = Math.min(3, Math.floor(t));
      const c = mix(stops[i], stops[i + 1], t - i);
      return `rgb(${c[0]},${c[1]},${c[2]})`;
    };
  }

  function renderHeat() {
    const d = state.data, g = d.goslim[state.goslimFeature];
    segmented($("goslim-seg"), [["gene_body", "Gene body", !!d.goslim.gene_body], ["cds", "CDS", !!d.goslim.cds]],
      state.goslimFeature, (f) => { state.goslimFeature = f; renderHeat(); });
    const spIdx = [...new Set(g.cells.map((c) => c[1]))].sort((a, b) => a - b);
    const species = spIdx.map((i) => d.species[i]);
    const lim = 0.3, color = divergingScale();
    const cellW = 44, cellH = 15, labelW = 300, headH = 110;
    const W = labelW + cellW * species.length + 8, H = headH + cellH * g.terms.length + 4;
    const svg = sv("svg", { viewBox: `0 0 ${W} ${H}`, width: W, height: H, role: "img",
      "aria-label": "GO Slim term by species heatmap of median CpG O/E shift" });
    svg.style.width = `${W}px`;
    species.forEach((s, j) => {
      const x = labelW + cellW * j + cellW / 2;
      svg.append(sv("text", { x, y: headH - 8, transform: `rotate(-45 ${x} ${headH - 8})`, "text-anchor": "start", class: "ink2", "font-style": "italic" }, s.organism));
    });
    g.terms.forEach((t, i) => {
      const y = headH + cellH * i;
      const tag = t.class === "housekeeping" ? "  [HK]" : t.class === "inducible" ? "  [IND]" : "";
      const lab = sv("text", { x: labelW - 8, y: y + cellH - 4, "text-anchor": "end", class: tag ? "ink" : "ink2" });
      lab.append(sv("tspan", {}, t.name));
      if (tag) lab.append(sv("tspan", { "font-weight": "600" }, tag));
      svg.append(lab);
    });
    const colOf = new Map(spIdx.map((si, j) => [si, j]));
    for (const [ti, si, shift, q, n, med] of g.cells) {
      const x = labelW + cellW * colOf.get(si), y = headH + cellH * ti;
      const r = sv("rect", { x: x + 1, y: y + 1, width: cellW - 2, height: cellH - 2, rx: 2,
        fill: color(shift / lim), opacity: q < 0.05 ? 1 : 0.35 });
      const sp = d.species[si], term = g.terms[ti];
      hoverable(r, () => [
        { l: `${term.name}${term.class ? ` (${term.class})` : ""}` },
        { l: sp.organism },
        { v: (shift > 0 ? "+" : "") + fmt(shift, 3), l: "vs other annotated genes" },
        { v: fmt(med, 3), l: `median CpG O/E · ${fmtInt(n)} genes` },
        { v: q < 0.001 ? fmtP(q) : fmt(q, 3), l: "q (BH)" },
      ]);
      r.addEventListener("pointerenter", () => r.setAttribute("stroke", css("--ink")));
      r.addEventListener("pointerleave", () => r.removeAttribute("stroke"));
      svg.append(r);
    }
    $("heat-chart").replaceChildren(svg);
    const lg = $("heat-legend");
    lg.replaceChildren();
    const bar = el("span", { class: "bar" });
    bar.style.background = `linear-gradient(90deg, ${[-1, -0.5, 0, 0.5, 1].map((v) => color(v)).join(",")})`;
    lg.append(el("span", {}, `−${lim}`), bar, el("span", {}, `+${lim}`), el("span", {}, "median CpG O/E, term − other genes (clipped)"));
  }

  // ---------- validation ----------
  function renderValidation() {
    const d = state.data, box = $("validation");
    box.replaceChildren();
    const rows = d.species.flatMap((s) => (s.validation || []).flatMap((v) =>
      Object.entries(v.features).map(([f, st]) => ({ s, v, f, st }))));
    if (!rows.length) {
      box.append(el("p", { class: "empty" }, "Whole-genome bisulfite sequencing (WGBS) alignments are running for Magallana gigas (gonad; eight female gonads) and Apis mellifera (sperm). Results will appear here: correlation between CpG O/E and measured gene methylation, and how often the low-CpG O/E class matches methylated genes."));
      return;
    }
    const wrap = el("div", { class: "table-scroll" }), t = el("table", { class: "data" });
    const head = ["Species", "Dataset", "Feature", "Genes", "Spearman ρ", "AUC", "Misclassified", "Precision", "Recall", "Median meth, low / high class"];
    const thead = el("thead"), hr = el("tr");
    head.forEach((h, i) => hr.append(el("th", { scope: "col", class: i >= 3 ? "num" : "" }, h)));
    thead.append(hr);
    const tb = el("tbody");
    for (const { s, v, f, st } of rows) {
      const tr = el("tr");
      const lc = st.low_class || {};
      const cells = [[s.organism, "sp"], [v.label, ""], [FEATURE_LABEL[f], ""], [fmtInt(st.n_genes), "num"],
        [fmt(st.spearman_rho), "num"], [fmt(st.auc_low_oe_predicts_methylated), "num"],
        [lc.misclassification == null ? "–" : `${Math.round(lc.misclassification * 100)} %`, "num"],
        [fmt(lc.precision), "num"], [fmt(lc.recall), "num"],
        [`${fmt(lc.median_meth_low)} / ${fmt(lc.median_meth_high)}`, "num"]];
      for (const [txt, cls] of cells) tr.append(el("td", { class: cls }, txt));
      const what = (v.source || "").split(";")[0].replace(/\s*WGBS\b/, "");
      if (what) tr.children[1].append(el("div", { class: "sub" }, what));
      tr.title = v.source;
      tb.append(tr);
    }
    t.append(thead, tb);
    wrap.append(t);
    box.append(wrap, el("p", { class: "note" }, "Low class: posterior probability > 0.5 for the low CpG O/E mixture component. Methylated gene: weighted methylation ≥ 0.5."));
  }

  // ---------- boot ----------
  function renderAll() {
    renderTiles(state.data);
    renderTable();
    renderDist();
    renderEffect();
    renderHeat();
    renderValidation();
  }

  fetch("data/dashboard.json")
    .then((r) => { if (!r.ok) throw new Error(r.status); return r.json(); })
    .then((d) => {
      state.data = d;
      renderMeta(d);
      const sel = $("species-select");
      for (const s of d.species) sel.append(el("option", { value: s.accession }, `${s.organism} (${s.group}${s.tier === "tsa" ? ", transcriptome" : ""})`));
      sel.addEventListener("change", () => selectSpecies(sel.value, false));
      state.species = d.species[0].accession;
      sel.value = state.species;
      renderAll();
      let t;
      window.addEventListener("resize", () => { clearTimeout(t); t = setTimeout(() => { renderDist(); renderEffect(); }, 150); });
      window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", renderAll);
    })
    .catch((e) => {
      $("tiles").replaceChildren(el("p", { class: "empty" }, `Could not load data/dashboard.json (${e.message}).`));
    });
})();
