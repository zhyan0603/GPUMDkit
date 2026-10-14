/* GPUMDkit website — publications page: search, year filter, grouped timeline */
(function () {
  "use strict";
  const GK = window.GK, $ = GK.$, $$ = GK.$$;
  const list = $("[data-pub-list]"), filters = $("[data-years]"), search = $("[data-search]"), count = $("[data-result-count]");
  if (!list) return;
  let pubs = [], year = "all";

  const mark = (text, q) => {
    const s = GK.esc(text);
    if (!q) return s;
    const i = text.toLowerCase().indexOf(q);
    return i < 0 ? s : GK.esc(text.slice(0, i)) + `<mark class="hit">${GK.esc(text.slice(i, i + q.length))}</mark>` + GK.esc(text.slice(i + q.length));
  };
  function renderFilters() {
    const years = Array.from(new Set(pubs.map(p => p.year))).sort((a, b) => b - a);
    filters.innerHTML = `<button class="chip${year === "all" ? " active" : ""}" data-y="all">${GK.lang === "zh" ? GK.t("pp.all") : "All years"}</button>` +
      years.map(y => `<button class="chip${year === y ? " active" : ""}" data-y="${y}">${y}</button>`).join("");
  }
  function render() {
    const q = (search.value || "").trim().toLowerCase();
    const hit = pubs.filter(p => (year === "all" || p.year === year) && (!q || [p.title, p.authors, p.journal].join(" ").toLowerCase().includes(q)));
    const zh = GK.lang === "zh";
    count.textContent = zh ? `${hit.length} ${GK.t("pp.count")}` : `${hit.length} ${hit.length === 1 ? "paper" : "papers"}`;
    if (!hit.length) { list.innerHTML = `<div class="empty">${zh ? "没有匹配的论文" : "No matching papers"}</div>`; return; }
    const byYear = {};
    hit.forEach(p => { (byYear[p.year] = byYear[p.year] || []).push(p); });
    list.innerHTML = Object.keys(byYear).sort((a, b) => b - a).map(y => `<div class="year-block">
      <div class="y"><b>${y}</b><span>${byYear[y].length} ${zh ? GK.t("pp.papers") : byYear[y].length === 1 ? "paper" : "papers"}</span></div>
      <div class="pub-list">${byYear[y].map((p, i) => `<a class="card pub-row reveal" style="--d:${i * 0.04}s" href="${p.link}" target="_blank" rel="noopener">
        <div><h4>${mark(p.title, q)}</h4><div class="au">${mark(p.authors, q)}</div><div class="jr">${mark(p.journal, q)}</div></div>
        <span class="go"><svg><use href="#i-arrow"/></svg></span></a>`).join("")}</div></div>`).join("");
    $$(".reveal", list).forEach(el => GK.observe(el));
  }
  GK.json("publications.json").then(data => {
    pubs = data;
    $$("[data-pub-count]").forEach(el => { el.dataset.count = data.length; GK.countUp(el); });
    renderFilters(); render();
  }).catch(() => { list.innerHTML = `<div class="empty">Publications could not be loaded.</div>`; });
  filters.addEventListener("click", e => { const b = e.target.closest("[data-y]"); if (!b) return; year = b.dataset.y; renderFilters(); render(); });
  search.addEventListener("input", render);
  GK.on("lang", () => { renderFilters(); render(); });
})();
