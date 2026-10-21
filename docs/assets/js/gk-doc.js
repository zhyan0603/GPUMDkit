/* GPUMDkit website — tutorial pages: language switch, table-of-contents highlight, image zoom. */
(function () {
  "use strict";
  const GK = window.GK, $ = GK.$, $$ = GK.$$;
  const page = document.body.dataset.lang, alt = document.body.dataset.alt;

  // each language has its own page: switching language opens the matching page
  GK.on("lang", () => { if (alt && GK.lang !== page) location.href = alt + location.hash; });

  // highlight the heading in view in the right-hand contents
  const links = $$(".doc-toc a[href^='#']");
  const heads = links.map(a => document.getElementById(decodeURIComponent(a.getAttribute("href").slice(1)))).filter(Boolean);
  function spy() {
    let cur = heads[0];
    heads.forEach(h => { if (h.getBoundingClientRect().top < innerHeight * 0.3) cur = h; });
    links.forEach(a => a.classList.toggle("on", cur && a.getAttribute("href") === "#" + cur.id));
  }
  if (heads.length) { addEventListener("scroll", spy, { passive: true }); spy(); }

  // annotated screenshots: hovering an entry highlights its region
  $$("[data-shot]").forEach(box => {
    const rgs = $$(".rg", box), btns = $$(".shot-legend button", box);
    const on = n => { rgs.forEach(r => r.classList.toggle("on", r.dataset.n === n)); btns.forEach(b => b.classList.toggle("on", b.dataset.n === n)); };
    btns.forEach(b => ["mouseenter", "focus", "click"].forEach(ev => b.addEventListener(ev, () => on(b.dataset.n))));
    rgs.forEach(r => $("b", r).addEventListener("mouseenter", () => on(r.dataset.n)));
    box.addEventListener("mouseleave", () => on(""));
  });

  // figures open full size
  $$(".doc-body img").forEach(img => img.addEventListener("click", () => GK.lightbox([{ src: img.currentSrc || img.src, title: img.alt }], 0)));

  // open a hash target that sits inside a collapsed block
  if (location.hash) {
    try { const t = document.getElementById(decodeURIComponent(location.hash.slice(1))); if (t && t.closest("details")) t.closest("details").open = true; } catch (e) { /* ignore */ }
  }
})();
