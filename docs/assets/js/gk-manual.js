/* GPUMDkit website — manual page: sidebar scroll-spy, tabs, annotated screenshots
   and the searchable command reference (from command_reference.tsv). */
(function () {
  "use strict";
  const GK = window.GK, $ = GK.$, $$ = GK.$$, esc = GK.esc;
  const zh = () => GK.lang === "zh";

  /* ------------------------------------------------------------ sidebar scroll-spy */
  const links = $$(".man-side a[href^='#']");
  const chaps = links.map(a => $(a.getAttribute("href"))).filter(Boolean);
  function setOn(id) {
    links.forEach(a => a.classList.toggle("on", a.getAttribute("href") === "#" + id));
    const cur = $(".man-side a.on");
    if (cur && getComputedStyle($(".man-side")).overflowX === "auto") cur.scrollIntoView({ inline: "center", block: "nearest" });
  }
  function spy() {
    let best = chaps[0];
    chaps.forEach(c => { if (c.getBoundingClientRect().top < innerHeight * 0.35) best = c; });
    if (best) setOn(best.id);
  }
  addEventListener("scroll", spy, { passive: true });
  spy();

  /* ------------------------------------------------------------ tabs: [data-tabs] buttons toggle sibling panes */
  $$("[data-tabs]").forEach(group => {
    const scope = group.closest("[data-tabscope]") || group.parentElement;
    const buttons = $$("[data-pane]", group);
    function show(id) {
      buttons.forEach(b => b.classList.toggle("active", b.dataset.pane === id));
      $$(".pane", scope).forEach(p => { if (p.dataset.group === group.dataset.tabs) p.classList.toggle("on", p.id === id); });
      group.dispatchEvent(new CustomEvent("pane", { detail: id }));
    }
    buttons.forEach(b => b.addEventListener("click", () => show(b.dataset.pane)));
  });

  /* ------------------------------------------------------------ annotated screenshots: legend <-> region */
  $$("[data-shot]").forEach(box => {
    const rgs = $$(".rg", box), btns = $$(".shot-legend button", box);
    function on(n) {
      rgs.forEach(r => r.classList.toggle("on", r.dataset.n === n));
      btns.forEach(b => b.classList.toggle("on", b.dataset.n === n));
    }
    btns.forEach(b => {
      b.addEventListener("mouseenter", () => on(b.dataset.n));
      b.addEventListener("focus", () => on(b.dataset.n));
      b.addEventListener("click", () => on(b.dataset.n));
    });
    rgs.forEach(r => $("b", r).addEventListener("mouseenter", () => on(r.dataset.n)));
    box.addEventListener("mouseleave", () => on(""));
  });

  /* ------------------------------------------------------------ anatomy highlights */
  $$(".anat").forEach(box => {
    const hs = $$(".hl", box), lis = $$("ol li", box);
    function on(n) { hs.forEach(h => h.classList.toggle("on", h.dataset.n === n)); lis.forEach((li, i) => li.classList.toggle("on", String(i + 1) === n)); }
    hs.forEach(h => h.addEventListener("mouseenter", () => on(h.dataset.n)));
    lis.forEach((li, i) => li.addEventListener("mouseenter", () => on(String(i + 1))));
    box.addEventListener("mouseleave", () => on(""));
  });

  /* screenshots open full size */
  $$(".shot img").forEach(img => { img.addEventListener("click", () => GK.lightbox([{ src: img.getAttribute("src"), title: img.alt }], 0)); });

  /* ------------------------------------------------------------ figure gallery */
  const figs = $$(".fig");
  const figItems = figs.map(f => ({ src: $("img", f).getAttribute("src"), title: $("img", f).alt, command: $("code", f).textContent }));
  figs.forEach(f => f.addEventListener("click", () => {
    const vis = figs.filter(x => !x.classList.contains("hide"));
    GK.lightbox(vis.map(x => figItems[figs.indexOf(x)]), vis.indexOf(f));
  }));
  const figTabs = $("[data-tabs='figs']");
  if (figTabs) figTabs.addEventListener("pane", e => {
    const cat = e.detail.replace("fg-", "");
    figs.forEach(f => f.classList.toggle("hide", cat !== "all" && f.dataset.cat !== cat));
  });

  /* ------------------------------------------------------------ command reference */
  const refBody = $("[data-ref-body]"), refFilters = $("[data-ref-filters]"), refSearch = $("[data-ref-search]"), refCount = $("[data-ref-count]");
  const CAT = { main: ["Main", "主要"], format: ["Format conversion", "格式转换"], sampling: ["Sampling", "采样"], calc: ["Calculators", "计算器"],
    analyze: ["Analyzers", "分析器"], utility: ["Utilities", "工具"], plots: ["Plots", "绘图"] };
  let rows = [], cat = "all";

  function loadRows() {
    const plots = Object.entries((window.GK_DATA || {}).plots || {}).flatMap(([p, desc]) => p.split("/").map((n, i) => {
      const name = i === 0 ? n : p.split("/")[0].replace(/\d*$/, "") + n;
      return { cat: "plots", syntax: "gpumdkit.sh -plt " + name, en: desc, zh: desc, status: "" };
    }));
    return fetch("command_reference.tsv").then(r => { if (!r.ok) throw new Error(r.status); return r.text(); }).then(txt => {
      const lines = txt.replace(/\r/g, "").split("\n").filter(Boolean);
      const head = lines.shift().split("\t");
      const ix = k => head.indexOf(k);
      rows = lines.map(l => { const c = l.split("\t"); return { cat: c[ix("category")], syntax: c[ix("syntax")], en: c[ix("description_en")], zh: c[ix("description_zh")], status: c[ix("status")] }; })
        .map(r => (r.cat === "plot" ? Object.assign(r, { cat: "plots" }) : r));
      // the TSV lists only a few plot commands with arguments; add the remaining types from the site data
      rows = rows.concat(plots.filter(p => !rows.some(r => r.syntax === p.syntax || r.syntax.startsWith(p.syntax + " "))));
    });
  }
  function renderFilters() {
    if (!refFilters) return;
    const cats = ["all"].concat(Object.keys(CAT).filter(k => rows.some(r => r.cat === k)));
    refFilters.innerHTML = cats.map(k => `<button class="chip${k === cat ? " active" : ""}" data-cat="${k}">${k === "all" ? (zh() ? "全部" : "All") : CAT[k][zh() ? 1 : 0]}</button>`).join("");
  }
  function renderRows() {
    if (!refBody) return;
    const q = ((refSearch && refSearch.value) || "").trim().toLowerCase();
    const list = rows.filter(r => (cat === "all" || r.cat === cat) && (!q || (r.syntax + " " + r.en + " " + r.zh).toLowerCase().includes(q)));
    refBody.innerHTML = list.map(r => `<tr><td class="cmd"><button data-copy="${esc(r.syntax)}" title="${zh() ? "点击复制" : "Click to copy"}"><code>${esc(r.syntax)}</code></button></td>
      <td>${esc(zh() ? r.zh : r.en)}</td><td><span class="st">${esc(CAT[r.cat] ? CAT[r.cat][zh() ? 1 : 0] : r.cat)}</span></td></tr>`).join("") ||
      `<tr><td colspan="3" class="ref-empty">${zh() ? "没有匹配的命令" : "No matching commands"}</td></tr>`;
    if (refCount) refCount.textContent = list.length + (zh() ? " 条命令" : list.length === 1 ? " command" : " commands");
  }
  if (refBody) {
    loadRows().then(() => { renderFilters(); renderRows(); }).catch(() => {
      refBody.innerHTML = `<tr><td colspan="3" class="ref-empty">${zh() ? "无法加载命令表（本地 file:// 预览不支持），请通过网站访问，或直接查看" : "The command table could not be loaded (not available from a file:// preview). Open the website, or read"} <a class="blue" href="command_reference.tsv">command_reference.tsv</a>.</td></tr>`;
    });
    if (refFilters) refFilters.addEventListener("click", e => { const b = e.target.closest("[data-cat]"); if (b) { cat = b.dataset.cat; renderFilters(); renderRows(); } });
    if (refSearch) refSearch.addEventListener("input", renderRows);
    GK.on("lang", () => { renderFilters(); renderRows(); });
  }

  /* placeholders follow the language */
  function placeholders() {
    if (refSearch) refSearch.placeholder = zh() ? "筛选命令…" : "Filter commands…";
  }
  placeholders();
  GK.on("lang", placeholders);

  /* open a hash target even when it sits inside a collapsed FAQ item */
  if (location.hash) { const t = $(location.hash); if (t && t.closest("details")) t.closest("details").open = true; }
})();
