/* GPUMDkit website — gallery page: search, category filters, masonry, lightbox */
(function () {
  "use strict";
  const GK = window.GK, $ = GK.$, $$ = GK.$$;
  const grid = $("[data-masonry]"), filters = $("[data-filters]"), search = $("[data-search]"), count = $("[data-result-count]");
  if (!grid) return;

  // category from the plot command
  const CATS = [
    ["train", "NEP training", "NEP 训练", /-plt (train|prediction|train_test|parity_density|des|dimer|force_errors|restart|net_force)\b/],
    ["md", "MD & structure", "MD 与结构", /-plt (thermo|rdf|xrd|vac|cohesive)\b/],
    ["diff", "Diffusion", "扩散", /-plt (msd|msd_all|msd_conv|sdc|D|sigma|doas)\b|D\/sigma/],
    ["heat", "Heat transport", "热输运", /-plt (emd|nemd|hnemd|viscosity)\b/],
    ["phon", "Phonons & charge", "声子与电荷", /-plt (pdos|phonon|charge|born_charge)\b/],
  ];
  const catOf = cmd => (CATS.find(c => c[3].test(cmd)) || ["other"])[0];
  let items = [], active = "all";

  function renderFilters() {
    const zh = GK.lang === "zh";
    filters.innerHTML = `<button class="chip${active === "all" ? " active" : ""}" data-cat="all">${GK.t("gp.all") === "gp.all" ? "All" : GK.t("gp.all")}</button>` +
      CATS.map(c => `<button class="chip${active === c[0] ? " active" : ""}" data-cat="${c[0]}">${zh ? c[2] : c[1]}</button>`).join("");
  }
  function apply() {
    const q = (search.value || "").trim().toLowerCase();
    let n = 0;
    $$(".m-item", grid).forEach((el, i) => {
      const it = items[i];
      const ok = (active === "all" || it.cat === active) && (!q || (it.title + " " + it.command).toLowerCase().includes(q));
      el.classList.toggle("hide", !ok); if (ok) n++;
    });
    count.textContent = `${n} ${GK.lang === "zh" ? GK.t("gp.count") : n === 1 ? "figure" : "figures"}`;
  }
  GK.json("gallery.json").then(data => {
    items = data.map(d => ({ src: d.image, title: d.title, command: d.command, cat: catOf(d.command) }));
    grid.innerHTML = items.map((it, i) => {
      const c = CATS.find(x => x[0] === it.cat);
      return `<figure class="m-item reveal" style="--d:${(i % 6) * 0.05}s" data-i="${i}"><div class="im"><img src="${it.src}" alt="${GK.esc(it.title)}" loading="lazy"></div>
        <figcaption><div class="meta"><b>${GK.esc(it.title)}</b>${c ? `<span class="tag" data-cat-label="${c[0]}">${GK.lang === "zh" ? c[2] : c[1]}</span>` : ""}</div>
        <div class="cmdline"><code>${GK.esc(it.command)}</code><button class="copy-btn" data-copy="${GK.esc(it.command)}" aria-label="Copy command"><svg class="cp"><use href="#i-copy"/></svg><svg class="ok"><use href="#i-check"/></svg></button></div></figcaption></figure>`;
    }).join("");
    $$(".reveal", grid).forEach(el => GK.observe(el));
    renderFilters(); apply();
    grid.addEventListener("click", e => {
      if (e.target.closest(".copy-btn")) return;
      const f = e.target.closest(".m-item");
      if (!f) return;
      const visible = $$(".m-item:not(.hide)", grid).map(el => +el.dataset.i);
      GK.lightbox(visible.map(i => items[i]), visible.indexOf(+f.dataset.i));
    });
  }).catch(() => { grid.innerHTML = `<div class="empty">Gallery could not be loaded.</div>`; });

  filters.addEventListener("click", e => { const b = e.target.closest("[data-cat]"); if (!b) return; active = b.dataset.cat; renderFilters(); apply(); });
  search.addEventListener("input", apply);
  GK.on("lang", () => {
    renderFilters(); apply();
    $$("[data-cat-label]", grid).forEach(t => { const c = CATS.find(x => x[0] === t.dataset.catLabel); t.textContent = GK.lang === "zh" ? c[2] : c[1]; });
  });
})();
