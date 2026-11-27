/* CpG O/E Atlas: renders docs/data/dashboard.json (built by analysis/build_dashboard.py). */
(() => {
  "use strict";
  const SVG = "http://www.w3.org/2000/svg";
  const FEATURE_LABEL = { cds: "CDS", gene_body: "Gene body", upstream: "Upstream 1 kb" };
  const state = { data: null, species: null, feature: "gene_body", goslimFeature: "gene_body", covFeature: "gene_body",
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
  // Narrow charts: "Drosophila melanogaster" -> "D. melanogaster".
  const spName = (name, W) => (W < 520 ? name.replace(/^([A-Z])[a-z]+ /, "$1. ") : name);
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
    renderPCA();
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
    svg.append(sv("text", { x: X(lo), y: H, "text-anchor": "start" }, W < 520 ? "← HK lower" : "← housekeeping lower"));
    svg.append(sv("text", { x: X(hi), y: H, "text-anchor": "end" }, "higher →"));
    const series = [["gene_body", css("--series-1"), "Gene body"], ["cds", css("--series-2"), "CDS"]];
    sp.forEach((s, i) => {
      const y = m.t + rowH * i + rowH / 2;
      const lab = sv("text", { x: m.l - 8, y: y + 4, "text-anchor": "end", class: "ink2" });
      lab.append(sv("tspan", { "font-style": "italic" }, spName(s.organism, W)));
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

  // ---------- covariate check ----------
  const signed = (x) => (x == null ? "–" : (x > 0 ? "+" : "") + fmt(x));
  function renderCov() {
    const d = state.data, card = $("cov-card");
    const all = d.species.filter((s) => s.covariate);
    if (!all.length) { card.hidden = true; return; }
    card.hidden = false;
    const has = (f) => all.some((s) => s.covariate[f]);
    if (!has(state.covFeature)) state.covFeature = has("gene_body") ? "gene_body" : "cds";
    const f = state.covFeature;
    segmented($("cov-seg"), [["gene_body", "Gene body", has("gene_body")], ["cds", "CDS", has("cds")]],
      f, (k) => { state.covFeature = k; renderCov(); });
    const sp = all.filter((s) => s.covariate[f]);

    // summary: focal species with a significant raw gap, and how many keep it after adjustment
    const withGap = sp.filter((s) => s.covariate[f].rb_raw < 0 && s.covariate[f].p_raw < 0.05);
    const kept = withGap.filter((s) => s.covariate[f].rb_residual < 0 && s.covariate[f].p_residual < 0.05);
    const lost = withGap.filter((s) => !kept.includes(s)).map((s) => s.organism);
    const sum = $("cov-summary");
    sum.replaceChildren(el("strong", {}, "Result: "),
      `in ${kept.length} of ${withGap.length} species with a housekeeping-lower gap (${f === "cds" ? "CDS" : "gene body"}), the gap stays significant after adjustment` +
      (lost.length ? `; it does not in ${lost.length > 1 ? lost.slice(0, -1).join(", ") + " and " + lost.at(-1) : lost[0]}. ` : ". ") +
      "Where WGBS is available, housekeeping genes are directly and far more methylated, before and after the same adjustment.");

    // dumbbell: before (hollow) -> after (filled)
    const box = $("cov-chart"), W = Math.max(300, box.clientWidth);
    const rowH = 24, labelW = Math.min(190, W * 0.42);
    const m = { l: labelW, r: 16, t: 8, b: 30 }, H = m.t + m.b + rowH * sp.length, iw = W - m.l - m.r;
    const vals = sp.flatMap((s) => [s.covariate[f].rb_raw, s.covariate[f].rb_residual]).filter((v) => v != null);
    const lo = Math.min(-0.5, Math.floor(Math.min(...vals) * 10) / 10), hi = Math.max(0.2, Math.ceil(Math.max(...vals) * 10) / 10);
    const X = (x) => m.l + ((x - lo) / (hi - lo)) * iw;
    const svg = sv("svg", { viewBox: `0 0 ${W} ${H}`, role: "img",
      "aria-label": `Housekeeping vs inducible rank-biserial before and after GC and length adjustment, ${FEATURE_LABEL[f]}` });
    for (let t = Math.ceil(lo * 10) / 10; t <= hi + 1e-9; t += 0.1) {
      const tt = Math.round(t * 10) / 10;
      svg.append(sv("line", { x1: X(tt), x2: X(tt), y1: m.t, y2: H - m.b, stroke: tt === 0 ? css("--axis") : css("--grid"), "stroke-width": tt === 0 ? 1.5 : 1 }));
      if (Math.abs(tt * 10) % 2 === 0) svg.append(sv("text", { x: X(tt), y: H - 12, "text-anchor": "middle", class: "tick-label" }, (tt > 0 ? "+" : "") + tt.toFixed(1)));
    }
    svg.append(sv("text", { x: X(lo), y: H, "text-anchor": "start" }, W < 520 ? "← HK lower" : "← housekeeping lower"));
    svg.append(sv("text", { x: X(hi), y: H, "text-anchor": "end" }, "higher →"));
    const cBefore = css("--ink-2"), cAfter = css("--series-1");
    sp.forEach((s, i) => {
      const c = s.covariate[f], y = m.t + rowH * i + rowH / 2;
      const lab = sv("text", { x: m.l - 8, y: y + 4, "text-anchor": "end", class: "ink2" });
      lab.append(sv("tspan", { "font-style": "italic" }, spName(s.organism, W)));
      if (s.role !== "focal") lab.append(sv("tspan", {}, " (control)"));
      svg.append(lab);
      svg.append(sv("line", { x1: m.l, x2: W - m.r, y1: y, y2: y, stroke: css("--grid"), "stroke-dasharray": "1 3" }));
      svg.append(sv("line", { x1: X(c.rb_raw), x2: X(c.rb_residual), y1: y, y2: y, stroke: css("--axis"), "stroke-width": 2 }));
      svg.append(sv("circle", { cx: X(c.rb_raw), cy: y, r: 4.5, fill: css("--surface"), stroke: cBefore, "stroke-width": 1.5 }));
      svg.append(sv("circle", { cx: X(c.rb_residual), cy: y, r: 5, fill: cAfter, stroke: css("--surface"), "stroke-width": 2 }));
      const x0 = Math.min(X(c.rb_raw), X(c.rb_residual)) - 10, x1 = Math.max(X(c.rb_raw), X(c.rb_residual)) + 10;
      const hit = sv("rect", { x: x0, y: y - rowH / 2, width: x1 - x0, height: rowH, fill: "transparent" });
      hoverable(hit, () => [
        { l: `${s.organism} · ${FEATURE_LABEL[f]}` + (s.role === "control_unmethylated" ? " · unmethylated control"
          : s.role === "control_methylated" ? " · methylated control" : "") },
        { v: signed(c.rb_raw), l: `before adjustment · p ${fmtP(c.p_raw)}` },
        { key: cAfter, v: signed(c.rb_residual), l: `after GC and length · p ${fmtP(c.p_residual)}` },
        { v: `${fmt(c.median_gc_hk)} vs ${fmt(c.median_gc_ind)}`, l: "median GC, HK vs IND" },
        { v: `${Math.round(c.covariate_r2 * 100)} %`, l: "of CpG O/E variance explained by GC and length" },
        { v: signed(c.rb_gpc), l: "GpC O/E control (composition only)" },
      ]);
      svg.append(hit);
    });
    box.replaceChildren(svg);
    const lg = $("cov-legend");
    lg.replaceChildren();
    const a = el("span"), ka = el("span", { class: "key-dot" });
    ka.style.background = "transparent"; ka.style.border = `1.5px solid ${cBefore}`;
    a.append(ka, "Before adjustment");
    const b = el("span"), kb = el("span", { class: "key-dot" });
    kb.style.background = cAfter;
    b.append(kb, "After GC and length adjustment");
    lg.append(a, b);

    // WGBS direct check
    const wbox = $("cov-wgbs");
    wbox.replaceChildren();
    // One row per study: samples from the same project (e.g. eight female gonads) are merged,
    // with ranges for their values.
    const groups = new Map();
    for (const s of d.species) {
      const src = Object.fromEntries((s.validation || []).map((v) => [v.label, v.source || ""]));
      for (const w of s.covariate_wgbs || []) {
        const desc = (src[w.dataset] || w.dataset).split(";")[0].replace(/\s*WGBS\b/, "").replace(/, sample .*$/, "");
        const key = `${s.accession}|${desc}`;
        if (!groups.has(key)) groups.set(key, { s, desc, ws: [] });
        groups.get(key).ws.push(w);
      }
    }
    if (!groups.size) { wbox.append(el("p", { class: "empty" }, "No WGBS datasets.")); return; }
    const wrap = el("div", { class: "table-scroll" }), t = el("table", { class: "data compact" });
    const head = ["Dataset", "Median meth HK / IND", "Methylated HK / IND", "Effect", "Adjusted"];
    const thead = el("thead"), hr = el("tr");
    head.forEach((h, i) => hr.append(el("th", { scope: "col", class: i >= 1 ? "num" : "" }, h)));
    thead.append(hr);
    const tb = el("tbody");
    const range = (ws, k, f) => {
      const v = ws.map((w) => w[k]).filter((x) => x != null);
      const lo = Math.min(...v), hi = Math.max(...v);
      return f(lo) === f(hi) ? f(lo) : lo < 0 ? `${f(lo)} to ${f(hi)}` : `${f(lo)}–${f(hi)}`;
    };
    const pct = (x) => `${Math.round(x * 100)}`;
    for (const { s, desc, ws } of groups.values()) {
      const tr = el("tr");
      const ds = el("td");
      const what = ws.length > 1 ? `${desc}, ${ws.length} samples`
        : desc.includes(ws[0].dataset.split("_")[0]) ? desc : `${desc} (${ws[0].dataset})`;
      ds.append(el("span", { class: "sp" }, s.organism), el("div", { class: "sub" }, what));
      tr.append(ds);
      [[`${range(ws, "median_meth_hk", fmt)} / ${range(ws, "median_meth_ind", fmt)}`, "num"],
        [`${range(ws, "frac_methylated_hk", pct)} / ${range(ws, "frac_methylated_ind", pct)} %`, "num"],
        [range(ws, "rb_meth", signed), "num"], [range(ws, "rb_meth_adjusted", signed), "num"]]
        .forEach(([txt, cls]) => tr.append(el("td", { class: cls }, txt)));
      tr.title = ws.map((w) => `${w.dataset}: p ${fmtP(w.p_meth)} before, ${fmtP(w.p_meth_adjusted)} after adjustment; n ${w.n_housekeeping} / ${w.n_inducible}`).join("\n");
      tb.append(tr);
    }
    t.append(thead, tb);
    wrap.append(t);
    wbox.append(wrap, el("p", { class: "note" }, "HK: housekeeping, IND: inducible genes. Methylated: share of genes with methylation ≥ 0.5. Effect: rank-biserial correlation of measured methylation; below zero, housekeeping genes are more methylated. Adjusted: the same after GC and length adjustment. Hover a row for p-values."));
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
    // Size the label column to the longest term name (with its bold [HK]/[IND] tag).
    const probe = sv("svg", { width: 0, height: 0 });
    probe.style.position = "absolute";
    $("heat-chart").append(probe);
    const longest = Math.max(...g.terms.map((t) => {
      const tag = t.class === "housekeeping" ? "  [HK]" : t.class === "inducible" ? "  [IND]" : "";
      const lab = sv("text", { class: tag ? "ink" : "ink2" });
      lab.append(sv("tspan", {}, t.name));
      if (tag) lab.append(sv("tspan", { "font-weight": "600" }, tag));
      probe.append(lab);
      return lab.getComputedTextLength();
    }));
    probe.remove();
    const cellW = 44, cellH = 15, labelW = Math.ceil(longest) + 16, headH = 110;
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

  // ---------- trait space (PCA) ----------
  const CLADES = ["coral", "bivalve", "hemiptera", "tunicate", "seastar", "control"];
  const INV = ["invasive", "introduced", "native", "uncertain"];
  const INV_COLOR = { invasive: "--cat-2", introduced: "--cat-4", native: "--cat-1", uncertain: "--muted" };
  const pca = { set: "summary", feat: "gene_body", clade: "all", center: false, color: "group", view: "2d",
                x: 0, y: 1, z: 2, cmp: "native", centroids: true, weight: true, hidden: new Set(), rot: null, res: null };
  const SUMMARY = [
    ["median", "median", (f) => f.median], ["mean", "mean", (f) => f.mean],
    ["low_mean", "low mode, mean", (f) => f.low_mean], ["high_mean", "high mode, mean", (f) => f.high_mean],
    ["low_sd", "low mode, SD", (f) => f.low_sd], ["high_sd", "high mode, SD", (f) => f.high_sd],
    ["low_weight", "low mode, share", (f) => f.low_weight], ["frac_low", "low posterior share", (f) => f.frac_low_posterior],
    ["ashman_d", "Ashman's D", (f) => f.ashman_d],
    ["dbic", "log ΔBIC", (f) => (f.delta_bic == null ? null : Math.sign(f.delta_bic) * Math.log10(1 + Math.abs(f.delta_bic)))],
    ["gpc", "GpC O/E median", (f) => f.gpc_oe_median],
  ];

  // Jacobi eigen-decomposition of a symmetric matrix; returns eigenpairs sorted by value, descending.
  function eigSym(M) {
    const n = M.length, A = M.map((r) => r.slice()), V = A.map((_, i) => A.map((__, j) => (i === j ? 1 : 0)));
    for (let sweep = 0; sweep < 60; sweep++) {
      let off = 0;
      for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += A[i][j] ** 2;
      if (off < 1e-20) break;
      for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
        if (Math.abs(A[p][q]) < 1e-18) continue;
        const th = (A[q][q] - A[p][p]) / (2 * A[p][q]);
        const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1)), c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < n; k++) { const a = A[k][p], b = A[k][q]; A[k][p] = c * a - s * b; A[k][q] = s * a + c * b; }
        for (let k = 0; k < n; k++) { const a = A[p][k], b = A[q][k]; A[p][k] = c * a - s * b; A[q][k] = s * a + c * b; }
        for (let k = 0; k < n; k++) { const a = V[k][p], b = V[k][q]; V[k][p] = c * a - s * b; V[k][q] = s * a + c * b; }
      }
    }
    return A.map((_, i) => ({ val: Math.max(0, A[i][i]), vec: V.map((r) => r[i]) })).sort((a, b) => b.val - a.val);
  }

  // Seeded RNG so the permutation p-value does not jitter between renders.
  function rng(seed) { let s = seed >>> 0; return () => ((s = (s * 1664525 + 1013904223) >>> 0) / 4294967296); }

  function profile(sp, f) {
    const d = state.data;
    if (pca.set === "summary") {
      const fit = sp.fits[f];
      if (!fit) return null;
      const v = SUMMARY.map(([, , g]) => g(fit));
      v.push(sp.goslim_test?.[f]?.rank_biserial ?? null);
      return v.some((x) => x == null || Number.isNaN(x)) ? null : v;
    }
    if (pca.set === "shape") {
      const h = sp.hist[f];
      if (!h) return null;
      const n = h.reduce((a, b) => a + b, 0);
      return n ? h.map((c) => c / n) : null;
    }
    const g = d.goslim[f];
    if (!g || !pca._goslimBySp?.[f]) return null;
    const e = pca._goslimBySp[f].get(d.species.indexOf(sp));
    return e ? (pca.weight ? e.w : e.raw) : null;
  }
  function featureNames(f) {
    if (pca.set === "summary") return [...SUMMARY.map((x) => x[1]), "HK vs IND effect"];
    if (pca.set === "shape") return state.data.bins.slice(0, -1).map((b, i) => `CpG O/E ${b.toFixed(3)}–${state.data.bins[i + 1].toFixed(3)}`);
    return state.data.goslim[f].terms.map((t) => t.name);
  }

  function computePCA() {
    const d = state.data, f = pca.feat;
    if (pca.set === "goslim" && !pca._goslimBySp) {
      // Per species, two term vectors: raw shifts, and shifts shrunk toward 0 by their reliability
      // (empirical Bayes): w = tau2 / (tau2 + noise), noise = (pi/2) sigma^2 / n for a median of
      // n genes, sigma = typical per-gene SD of CpG O/E, tau2 = between-species variance of the
      // term's shift beyond noise. Terms resting on few genes are pulled toward "no shift".
      pca._goslimBySp = {};
      pca._goslimInfo = {};
      const med = (a) => { const b = [...a].sort((x, y) => x - y), h = b.length >> 1; return b.length % 2 ? b[h] : (b[h - 1] + b[h]) / 2; };
      for (const [ff, g] of Object.entries(d.goslim)) {
        const vars = d.species.map((sp) => sp.fits[ff]).filter((ft) => ft && ft.low_sd != null)
          .map((ft) => ft.low_weight * ft.low_sd ** 2 + ft.high_weight * ft.high_sd ** 2 +
            ft.low_weight * ft.high_weight * (ft.low_mean - ft.high_mean) ** 2);
        const sigma2 = med(vars); // median across species, so a degenerate fit cannot set it
        const noise = (n) => (Math.PI / 2) * sigma2 / n;
        const byTerm = g.terms.map(() => []);
        for (const c of g.cells) if (d.species[c[1]].role === "focal") byTerm[c[0]].push(c);
        const tau2 = byTerm.map((cs) => {
          if (cs.length < 3) return 0;
          const m = cs.reduce((a, c) => a + c[2], 0) / cs.length;
          const v = cs.reduce((a, c) => a + (c[2] - m) ** 2, 0) / (cs.length - 1);
          return Math.max(0, v - cs.reduce((a, c) => a + noise(c[4]), 0) / cs.length);
        });
        const m = new Map();
        const ws = [];
        for (const [ti, si, shift, , n] of g.cells) {
          if (!m.has(si)) m.set(si, { raw: new Array(g.terms.length).fill(0), w: new Array(g.terms.length).fill(0) });
          const w = tau2[ti] > 0 ? tau2[ti] / (tau2[ti] + noise(n)) : 0;
          m.get(si).raw[ti] = shift;
          m.get(si).w[ti] = w * shift;
          ws.push(w);
        }
        pca._goslimBySp[ff] = m;
        pca._goslimInfo[ff] = { sigma: Math.sqrt(sigma2), medW: med(ws), nZero: tau2.filter((t) => t === 0).length };
      }
    }
    const inClade = (s) => pca.clade === "all" || s.group === pca.clade;
    const pts = [], dropped = [];
    for (const s of d.species) {
      if (!inClade(s)) continue;
      const v = profile(s, f);
      const fit = s.role === "focal";
      if (!v) { if (fit) dropped.push(s); continue; }
      if (!fit && pca.center && pca.clade === "all" && s.group === "control") continue;
      pts.push({ s, v, fit });
    }
    const train = pts.filter((p) => p.fit);
    const p = train[0]?.v.length || 0;
    if (train.length < 4 || !p) return { empty: true, dropped };
    // clade means (from fitted species), optionally removed
    const cmean = {};
    for (const c of CLADES) {
      const rows = train.filter((t) => t.s.group === c);
      if (rows.length) cmean[c] = rows[0].v.map((_, j) => rows.reduce((a, r) => a + r.v[j], 0) / rows.length);
    }
    const base = (pt) => (pca.center && pca.clade === "all" && cmean[pt.s.group]
      ? pt.v.map((x, j) => x - cmean[pt.s.group][j]) : pt.v.slice());
    for (const pt of pts) pt.b = base(pt);
    const mu = new Array(p).fill(0), sd = new Array(p).fill(1);
    for (const t of train) t.b.forEach((x, j) => { mu[j] += x / train.length; });
    let clip = Infinity;
    if (pca.set === "summary") {
      // Robust scaling: median / MAD, clipped at ±3, so one degenerate mixture fit
      // (e.g. a "mode" made of a single gene) cannot dominate a component.
      const med = (a) => { const b = [...a].sort((x, y) => x - y), h = b.length >> 1; return b.length % 2 ? b[h] : (b[h - 1] + b[h]) / 2; };
      for (let j = 0; j < p; j++) {
        const col = train.map((t) => t.b[j]);
        mu[j] = med(col);
        const mad = 1.4826 * med(col.map((x) => Math.abs(x - mu[j])));
        sd[j] = mad > 1e-12 ? mad : 1;
      }
      clip = 3;
    }
    for (const pt of pts) pt.z = pt.b.map((x, j) => Math.max(-clip, Math.min(clip, (x - mu[j]) / sd[j])));
    if (clip < Infinity) { // re-centre clipped scores so the PCA is on a centred matrix
      const m2 = new Array(p).fill(0);
      for (const t of train) t.z.forEach((x, j) => { m2[j] += x / train.length; });
      for (const pt of pts) pt.z = pt.z.map((x, j) => x - m2[j]);
    }
    const C = Array.from({ length: p }, () => new Array(p).fill(0));
    for (const t of train) for (let i = 0; i < p; i++) { const zi = t.z[i]; if (!zi) continue; for (let j = i; j < p; j++) C[i][j] += zi * t.z[j]; }
    for (let i = 0; i < p; i++) for (let j = i; j < p; j++) { C[i][j] /= train.length - 1; C[j][i] = C[i][j]; }
    const eig = eigSym(C).slice(0, Math.min(10, p, train.length - 1));
    for (const e of eig) { // orient: largest |loading| positive
      let k = 0; e.vec.forEach((x, j) => { if (Math.abs(x) > Math.abs(e.vec[k])) k = j; });
      if (e.vec[k] < 0) e.vec = e.vec.map((x) => -x);
    }
    const total = C.reduce((a, r, i) => a + r[i], 0);
    for (const pt of pts) pt.pc = eig.map((e) => e.vec.reduce((a, w, j) => a + w * pt.z[j], 0));
    return { pts, train, eig, total, names: featureNames(f), dropped };
  }

  // Within-clade permutation test of invasive vs comparison centroids in the full (scaled) profile space.
  function invasionTest(res) {
    const lab = (s) => {
      const st = s.invasion;
      if (pca.cmp === "native") return st === "invasive" ? 1 : st === "native" ? 0 : null;
      if (pca.cmp === "nonnative") return st === "invasive" ? 1 : st === "native" || st === "introduced" ? 0 : null;
      return st === "invasive" || st === "introduced" ? 1 : st === "native" ? 0 : null; // any non-native vs native
    };
    const strata = [];
    for (const c of CLADES) {
      const rows = res.train.filter((t) => t.s.group === c && lab(t.s) !== null);
      const n1 = rows.filter((r) => lab(r.s) === 1).length, n0 = rows.length - n1;
      if (n1 && n0) strata.push({ c, rows: rows.map((r) => r.z), y: rows.map((r) => lab(r.s)), n1, n0 });
    }
    if (!strata.length) return null;
    const dist = (rows, y) => {
      const p = rows[0].length, a = new Array(p).fill(0), b = new Array(p).fill(0);
      let na = 0, nb = 0;
      rows.forEach((r, i) => { if (y[i]) { na++; r.forEach((x, j) => (a[j] += x)); } else { nb++; r.forEach((x, j) => (b[j] += x)); } });
      return Math.sqrt(a.reduce((s, x, j) => s + (x / na - b[j] / nb) ** 2, 0));
    };
    const w = strata.map((st) => (st.n1 * st.n0) / (st.n1 + st.n0)), W = w.reduce((a, b) => a + b, 0);
    const obs = strata.map((st) => dist(st.rows, st.y));
    const T = obs.reduce((a, x, i) => a + w[i] * x, 0) / W;
    const B = 4999, r = rng(20261127), ge = strata.map(() => 0);
    let geT = 0;
    for (let b = 0; b < B; b++) {
      let t = 0;
      strata.forEach((st, i) => {
        const y = st.y.slice();
        for (let k = y.length - 1; k > 0; k--) { const j = Math.floor(r() * (k + 1)); [y[k], y[j]] = [y[j], y[k]]; }
        const dd = dist(st.rows, y);
        if (dd >= obs[i] - 1e-12) ge[i]++;
        t += w[i] * dd;
      });
      if (t / W >= T - 1e-12) geT++;
    }
    return { strata: strata.map((st, i) => ({ ...st, d: obs[i], p: (ge[i] + 1) / (B + 1) })), T, p: (geT + 1) / (B + 1) };
  }

  const shapePath = (st, x, y, r) => {
    if (st === "invasive") { const h = r * 1.25; return `M${x},${y - h}L${x + h * 0.95},${y + h * 0.65}L${x - h * 0.95},${y + h * 0.65}Z`; }
    if (st === "introduced") { const h = r * 1.2; return `M${x},${y - h}L${x + h},${y}L${x},${y + h}L${x - h},${y}Z`; }
    if (st === "uncertain") { const h = r * 0.8; return `M${x - h},${y - h}h${2 * h}v${2 * h}h${-2 * h}Z`; }
    return `M${x - r},${y}a${r},${r} 0 1,0 ${2 * r},0a${r},${r} 0 1,0 ${-2 * r},0`;
  };
  const markColor = (s) => css(pca.color === "group" ? `--cat-${CLADES.indexOf(s.group) + 1}` : INV_COLOR[s.invasion] || "--muted");

  function renderPCA() {
    const d = state.data;
    // controls
    const hasF = (f) => d.species.some((s) => (pca.set === "goslim" ? d.goslim[f] : s.fits[f]));
    if (!hasF(pca.feat)) pca.feat = "cds";
    const feats = pca.set === "goslim" ? ["gene_body", "cds"] : ["gene_body", "cds", "upstream"];
    if (!feats.includes(pca.feat)) pca.feat = "gene_body";
    segmented($("pca-feat-seg"), feats.map((f) => [f, FEATURE_LABEL[f], hasF(f)]), pca.feat, (f) => { pca.feat = f; renderPCA(); });
    segmented($("pca-color-seg"), [["group", "Clade", true], ["invasion", "Invasion", true]], pca.color, (c) => { pca.color = c; pca.hidden.clear(); renderPCA(); });
    segmented($("pca-view-seg"), [["2d", "2D", true], ["3d", "3D", true]], pca.view, (v) => { pca.view = v; renderPCA(); });
    segmented($("pca-cmp-seg"), [["native", "vs native", true], ["nonnative", "vs all others", true], ["any", "non-native vs native", true]],
      pca.cmp, (c) => { pca.cmp = c; renderPCA(); });
    $("pca-z-wrap").hidden = pca.view !== "3d";
    $("pca-weight-wrap").hidden = pca.set !== "goslim";
    $("pca-center").disabled = pca.clade !== "all";

    const res = computePCA();
    pca.res = res;
    if (res.empty) {
      $("pca-chart").replaceChildren(el("p", { class: "empty" }, "Too few species with this profile for a PCA."));
      for (const id of ["pca-legend", "pca-test", "pca-load", "pca-scree"]) $(id).replaceChildren();
      $("pca-note").textContent = "";
      return;
    }
    const k = res.eig.length;
    for (const [id, key] of [["pca-x", "x"], ["pca-y", "y"], ["pca-z", "z"]]) {
      const sel = $(id);
      if (pca[key] >= k) pca[key] = Math.min(k - 1, { x: 0, y: 1, z: 2 }[key]);
      sel.replaceChildren(...res.eig.map((e, i) => el("option", { value: i }, `PC${i + 1} (${Math.round((100 * e.val) / res.total)} %)`)));
      sel.value = pca[key];
    }
    if (pca.view === "3d") draw3D(res); else draw2D(res);
    renderPcaLegend(res);
    renderPcaSide(res);
    const nFit = res.train.length, nProj = res.pts.length - nFit;
    $("pca-note").textContent = `${nFit} focal species fitted` + (nProj ? `, ${nProj} projected (hollow)` : "") +
      (res.dropped.length ? `. Not in this view (no ${FEATURE_LABEL[pca.feat].toLowerCase()} profile): ${res.dropped.map((s) => s.organism).join(", ")}.` : ".") +
      (pca.view === "3d" ? " Drag to rotate." : "") +
      (pca.set === "goslim" && pca.weight && pca._goslimInfo?.[pca.feat]
        ? ` Gene-count weighting: median weight ${pca._goslimInfo[pca.feat].medW.toFixed(2)}; ` +
          `${pca._goslimInfo[pca.feat].nZero} terms with no between-species signal beyond noise get weight 0.` : "");
  }

  function pointTip(pt) {
    const s = pt.s, ax = pca.view === "3d" ? [pca.x, pca.y, pca.z] : [pca.x, pca.y];
    return [
      { l: `${s.organism} · ${s.group}${s.tier === "tsa" ? " · transcriptome" : ""}${pt.fit ? "" : ` · ${ROLE_LABEL[s.role] || s.role} (projected)`}` },
      { key: INV_COLOR[s.invasion] ? css(INV_COLOR[s.invasion]) : undefined, v: s.invasion, l: s.invasion_source ? `· ${s.invasion_source}` : "" },
      ...(s.invasion_basis ? [{ l: s.invasion_basis }] : []),
      { l: ax.map((i) => `PC${i + 1} ${pt.pc[i] >= 0 ? "+" : ""}${pt.pc[i].toFixed(2)}`).join("   ") },
    ];
  }
  const visible = (s) => !pca.hidden.has(pca.color === "group" ? s.group : s.invasion);

  function addMark(svg, pt, x, y, r, opacity = 1) {
    const s = pt.s, color = markColor(s), sel = s.accession === state.species;
    const g = sv("g", { opacity: visible(s) ? opacity : 0.08 });
    g.append(sv("path", { d: shapePath(s.invasion, x, y, r), fill: pt.fit ? color : css("--surface"),
      stroke: pt.fit ? css("--surface") : color, "stroke-width": pt.fit ? 1.5 : 2 }));
    if (sel) g.append(sv("circle", { cx: x, cy: y, r: r + 5, fill: "none", stroke: css("--ink"), "stroke-width": 1.5 }));
    if (visible(s)) {
      const hit = sv("circle", { cx: x, cy: y, r: Math.max(10, r + 4), fill: "transparent", style: "cursor:pointer" });
      hoverable(hit, () => pointTip(pt));
      hit.addEventListener("click", () => selectSpecies(s.accession, false));
      g.append(hit);
    }
    svg.append(g);
    return g;
  }

  function draw2D(res) {
    const box = $("pca-chart"), W = Math.max(300, box.clientWidth), H = Math.min(560, Math.max(320, W * 0.72));
    box.classList.remove("rot");
    const m = { l: 40, r: 14, t: 12, b: 36 }, iw = W - m.l - m.r, ih = H - m.t - m.b;
    const xs = res.pts.map((p) => p.pc[pca.x]), ys = res.pts.map((p) => p.pc[pca.y]);
    const pad = (lo, hi) => { const r = hi - lo || 1; return [lo - r * 0.06, hi + r * 0.06]; };
    const [x0, x1] = pad(Math.min(...xs), Math.max(...xs)), [y0, y1] = pad(Math.min(...ys), Math.max(...ys));
    const X = (v) => m.l + ((v - x0) / (x1 - x0)) * iw, Y = (v) => m.t + ih - ((v - y0) / (y1 - y0)) * ih;
    const svg = sv("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": `PCA scatter, PC${pca.x + 1} vs PC${pca.y + 1}` });
    const ticks = (lo, hi) => { const st = 10 ** Math.floor(Math.log10((hi - lo) / 4)); const step = [1, 2, 5, 10].map((f) => f * st).find((s) => (hi - lo) / s <= 6); const out = []; for (let t = Math.ceil(lo / step) * step; t <= hi; t += step) out.push(+t.toFixed(10)); return out; };
    for (const t of ticks(x0, x1)) {
      svg.append(sv("line", { x1: X(t), x2: X(t), y1: m.t, y2: m.t + ih, stroke: t === 0 ? css("--axis") : css("--grid") }));
      svg.append(sv("text", { x: X(t), y: m.t + ih + 14, "text-anchor": "middle", class: "tick-label" }, String(t)));
    }
    for (const t of ticks(y0, y1)) {
      svg.append(sv("line", { x1: m.l, x2: m.l + iw, y1: Y(t), y2: Y(t), stroke: t === 0 ? css("--axis") : css("--grid") }));
      svg.append(sv("text", { x: m.l - 6, y: Y(t) + 4, "text-anchor": "end", class: "tick-label" }, String(t)));
    }
    const pct = (i) => Math.round((100 * res.eig[i].val) / res.total);
    svg.append(sv("text", { x: m.l + iw, y: H - 2, "text-anchor": "end", class: "ink2" }, `PC${pca.x + 1} (${pct(pca.x)} % of variance) →`));
    svg.append(sv("text", { x: 12, y: m.t + 4, transform: `rotate(-90 12 ${m.t + 4})`, "text-anchor": "end", class: "ink2" }, `PC${pca.y + 1} (${pct(pca.y)} %) →`));
    // clade centroids: native/comparison (ring) -> invasive (triangle), joined
    if (pca.centroids) {
      for (const c of CLADES) {
        const rows = res.train.filter((t) => t.s.group === c);
        const inv = rows.filter((t) => t.s.invasion === "invasive"), nat = rows.filter((t) => t.s.invasion === "native");
        if (!inv.length || !nat.length || pca.hidden.has(c)) continue;
        const cen = (a) => [a.reduce((s, t) => s + t.pc[pca.x], 0) / a.length, a.reduce((s, t) => s + t.pc[pca.y], 0) / a.length];
        const [ix, iy] = cen(inv), [nx, ny] = cen(nat), col = css(`--cat-${CLADES.indexOf(c) + 1}`);
        svg.append(sv("line", { x1: X(nx), y1: Y(ny), x2: X(ix), y2: Y(iy), stroke: col, "stroke-width": 2, "stroke-dasharray": "5 4", opacity: 0.8 }));
        svg.append(sv("path", { d: shapePath("native", X(nx), Y(ny), 9), fill: "none", stroke: col, "stroke-width": 2.5 }));
        svg.append(sv("path", { d: shapePath("invasive", X(ix), Y(iy), 10), fill: "none", stroke: col, "stroke-width": 2.5 }));
      }
    }
    const order = [...res.pts].sort((a, b) => (a.s.accession === state.species) - (b.s.accession === state.species));
    for (const pt of order) addMark(svg, pt, X(pt.pc[pca.x]), Y(pt.pc[pca.y]), 5);
    const sel = res.pts.find((p) => p.s.accession === state.species);
    if (sel) svg.append(sv("text", { x: X(sel.pc[pca.x]) + 11, y: Y(sel.pc[pca.y]) + 4, class: "ink", "font-style": "italic", "pointer-events": "none" }, sel.s.organism));
    box.replaceChildren(svg);
  }

  function draw3D(res) {
    const box = $("pca-chart"), W = Math.max(300, box.clientWidth), H = Math.min(560, Math.max(340, W * 0.72));
    box.classList.add("rot");
    if (!pca.rot) pca.rot = { a: -0.6, b: 0.35 };
    const ax = [pca.x, pca.y, pca.z];
    const sc = ax.map((i) => Math.max(...res.pts.map((p) => Math.abs(p.pc[i]))) || 1);
    const S = Math.max(...sc); // one scale for all three axes: distances stay comparable
    const R = Math.min(W, H) * 0.42, cx = W / 2, cy = H / 2;
    const svg = sv("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": `3D PCA, PC${ax[0] + 1}, PC${ax[1] + 1}, PC${ax[2] + 1}` });
    const draw = () => {
      svg.replaceChildren();
      const { a, b } = pca.rot, ca = Math.cos(a), sa = Math.sin(a), cb = Math.cos(b), sb = Math.sin(b);
      const proj = ([x, y, z]) => { // rotate about vertical (a) then tilt (b); y-up
        const x1 = ca * x + sa * z, z1 = -sa * x + ca * z;
        const y2 = cb * y - sb * z1, z2 = sb * y + cb * z1;
        return [cx + (x1 / S) * R, cy - (y2 / S) * R, z2 / S];
      };
      const pct = (i) => Math.round((100 * res.eig[i].val) / res.total);
      [[1, 0, 0], [0, 1, 0], [0, 0, 1]].forEach((u, i) => {
        const p0 = proj(u.map((v) => -v * S)), p1 = proj(u.map((v) => v * S));
        svg.append(sv("line", { x1: p0[0], y1: p0[1], x2: p1[0], y2: p1[1], stroke: css("--axis"), "stroke-width": 1.5 }));
        svg.append(sv("text", { x: p1[0], y: p1[1] - 4, "text-anchor": "middle", class: "ink2", "font-weight": "600" }, `PC${ax[i] + 1} (${pct(ax[i])} %)`));
      });
      const items = res.pts.map((pt) => ({ pt, q: proj(ax.map((i) => pt.pc[i])) })).sort((u, v) => u.q[2] - v.q[2]);
      for (const { pt, q } of items) addMark(svg, pt, q[0], q[1], 4 + 1.6 * (q[2] + 1), 0.55 + 0.45 * ((q[2] + 1) / 2));
    };
    draw();
    let drag = null;
    svg.addEventListener("pointerdown", (e) => {
      drag = { x: e.clientX, y: e.clientY, id: e.pointerId, captured: false };
    });
    svg.addEventListener("pointermove", (e) => {
      if (!drag) return;
      const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
      if (!drag.captured) {
        if (Math.abs(dx) + Math.abs(dy) < 3) return;
        drag.captured = true; svg.setPointerCapture(drag.id); box.classList.add("dragging");
      }
      hideTip();
      pca.rot.a += dx * 0.01; pca.rot.b = Math.max(-1.5, Math.min(1.5, pca.rot.b + dy * 0.01));
      drag.x = e.clientX; drag.y = e.clientY;
      draw();
    });
    const end = () => { drag = null; box.classList.remove("dragging"); };
    svg.addEventListener("pointerup", end);
    svg.addEventListener("pointercancel", end);
    box.replaceChildren(svg);
  }

  function renderPcaLegend(res) {
    const lg = $("pca-legend");
    lg.replaceChildren();
    const present = new Set(res.pts.map((p) => (pca.color === "group" ? p.s.group : p.s.invasion)));
    const keys = (pca.color === "group" ? CLADES : INV).filter((k) => present.has(k));
    for (const k of keys) {
      const b = el("button", { type: "button", class: "lg" + (pca.hidden.has(k) ? " off" : ""), "aria-pressed": String(!pca.hidden.has(k)) });
      const dot = el("span", { class: "key-dot" });
      dot.style.background = css(pca.color === "group" ? `--cat-${CLADES.indexOf(k) + 1}` : INV_COLOR[k]);
      const n = res.pts.filter((p) => (pca.color === "group" ? p.s.group : p.s.invasion) === k).length;
      b.append(dot, `${k} (${n})`);
      b.addEventListener("click", () => { pca.hidden.has(k) ? pca.hidden.delete(k) : pca.hidden.add(k); renderPCA(); });
      lg.append(b);
    }
    const shapes = el("span", { class: "ink2" });
    const mini = (st) => { const s = sv("svg", { width: 14, height: 14, viewBox: "0 0 14 14" }); s.append(sv("path", { d: shapePath(st, 7, 7.5, 4.5), fill: css("--ink-2") })); return s; };
    lg.append(el("span", {}, "  shape:"));
    for (const st of INV) { const s = el("span"); s.append(mini(st), st); lg.append(s); }
    if (pca.centroids && pca.view === "2d") {
      const s = el("span"), k = el("span", { class: "key-dash" });
      s.append(k, "clade centroid, native (ring) → invasive (triangle)");
      lg.append(s);
    }
    lg.append(el("span", { class: "hint" }, "Click a colour key to fade it."));
  }

  function renderPcaSide(res) {
    // test
    const box = $("pca-test");
    box.replaceChildren();
    const t = invasionTest(res);
    if (!t) { box.append(el("p", { class: "test-sum" }, "No clade has both groups in this view.")); }
    else {
      const tb = el("table", { class: "test-table" }), hr = el("tr");
      for (const h of ["Clade", "n inv / other", "Distance", "p"]) hr.append(el("th", {}, h));
      tb.append(hr);
      for (const st of t.strata) {
        const tr = el("tr");
        for (const v of [st.c, `${st.n1} / ${st.n0}`, fmt(st.d), fmtP(st.p)]) tr.append(el("td", {}, v));
        tb.append(tr);
      }
      const tr = el("tr", { class: "total" });
      for (const v of ["Combined", `${t.strata.reduce((a, s) => a + s.n1, 0)} / ${t.strata.reduce((a, s) => a + s.n0, 0)}`, fmt(t.T), fmtP(t.p)]) tr.append(el("td", {}, v));
      tb.append(tr);
      box.append(tb);
      const verdict = t.p < 0.05 ? "Invasive species sit apart from the comparison species more than shuffled labels would (p < 0.05)."
        : "No clear separation: invasive species are no further from the comparison species than shuffled labels.";
      box.append(el("p", { class: "test-sum" }, `${verdict} Distance between group centroids in the full profile space` +
        `${pca.set === "summary" ? " (in robust SD units)" : ""}; p from 4,999 within-clade shuffles. Labels: GISD and WRiMS where they cover the species.`));
    }
    // loadings
    const lb = $("pca-load");
    lb.replaceChildren();
    const list = el("div", { class: "load-list" });
    const axes = pca.view === "3d" ? [pca.x, pca.y, pca.z] : [pca.x, pca.y];
    for (const i of axes) {
      const vec = res.eig[i].vec, top = vec.map((v, j) => [j, v]).sort((a, b) => Math.abs(b[1]) - Math.abs(a[1])).slice(0, 5);
      const mx = Math.max(...vec.map(Math.abs));
      list.append(el("div", { class: "pc" }, `PC${i + 1} · ${Math.round((100 * res.eig[i].val) / res.total)} % of variance`));
      for (const [j, v] of top) {
        const r = el("div", { class: "load-row" }), bar = el("div", { class: "bar" }), fill = el("span");
        const w = (Math.abs(v) / mx) * 50;
        fill.style.width = `${w}%`; fill.style.left = v >= 0 ? "50%" : `${50 - w}%`;
        bar.append(fill);
        r.append(el("span", { class: "nm", title: res.names[j] }, res.names[j]), bar, el("span", { class: "v" }, (v > 0 ? "+" : "") + v.toFixed(2)));
        list.append(r);
      }
    }
    lb.append(list);
    // scree
    const sb = $("pca-scree"), W = 300, H = 70, n = res.eig.length, bw = (W - 10) / n;
    const svg = sv("svg", { viewBox: `0 0 ${W} ${H + 16}`, role: "img", "aria-label": "Variance explained per component" });
    res.eig.forEach((e, i) => {
      const f = e.val / res.total, h = f * H / Math.max(...res.eig.map((x) => x.val / res.total));
      const on = axes.includes(i);
      const r = sv("rect", { x: 5 + i * bw + 1, y: H - h, width: bw - 2, height: h, rx: 2, fill: on ? css("--series-1") : css("--hist") });
      hoverable(r, () => [{ v: `${(100 * f).toFixed(1)} %`, l: `of variance, PC${i + 1}` }]);
      svg.append(r, sv("text", { x: 5 + i * bw + bw / 2, y: H + 12, "text-anchor": "middle", class: "tick-label" }, `${i + 1}`));
    });
    const wrap = el("div", { class: "chart" });
    wrap.append(svg);
    sb.replaceChildren(el("p", { class: "chart-title" }, "Variance explained by each component (plotted ones in blue)"), wrap);
  }

  function initPCA() {
    const cl = $("pca-clade");
    cl.append(el("option", { value: "all" }, "All clades"));
    for (const c of CLADES.filter((c) => c !== "control")) cl.append(el("option", { value: c }, `${c} only`));
    cl.addEventListener("change", () => { pca.clade = cl.value; pca.hidden.clear(); renderPCA(); });
    $("pca-set").addEventListener("change", (e) => { pca.set = e.target.value; renderPCA(); });
    $("pca-center").addEventListener("change", (e) => { pca.center = e.target.checked; renderPCA(); });
    $("pca-weight").addEventListener("change", (e) => { pca.weight = e.target.checked; renderPCA(); });
    $("pca-centroids").addEventListener("change", (e) => { pca.centroids = e.target.checked; renderPCA(); });
    for (const [id, key] of [["pca-x", "x"], ["pca-y", "y"], ["pca-z", "z"]]) $(id).addEventListener("change", (e) => { pca[key] = +e.target.value; renderPCA(); });
  }

  // ---------- boot ----------
  function renderAll() {
    renderTiles(state.data);
    renderTable();
    renderDist();
    renderEffect();
    renderPCA();
    renderCov();
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
      initPCA();
      renderAll();
      let t;
      window.addEventListener("resize", () => { clearTimeout(t); t = setTimeout(() => { renderDist(); renderEffect(); renderPCA(); renderCov(); }, 150); });
      window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", renderAll);
    })
    .catch((e) => {
      $("tiles").replaceChildren(el("p", { class: "empty" }, `Could not load data/dashboard.json (${e.message}).`));
    });
})();
