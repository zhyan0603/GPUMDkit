/* GPUMDkit website — shared behaviour for every page:
   theme, language, navigation, reveal animations, counters, spotlight cards,
   copy buttons, toast, lightbox, command palette and GitHub stats. */
(function () {
  "use strict";

  const $ = (s, r = document) => r.querySelector(s);
  const $$ = (s, r = document) => Array.from(r.querySelectorAll(s));
  const store = {
    get(k) { try { return localStorage.getItem(k); } catch (e) { return null; } },
    set(k, v) { try { localStorage.setItem(k, v); } catch (e) { /* private mode */ } },
  };
  const reduced = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  // English strings for text created by scripts (Chinese ones live in gk-i18n.js)
  const EN = {
    "pal.ph": "Search commands, plot types or docs…",
    "pal.sections": "Sections",
    "pal.docs": "Documentation",
    "pal.cli": "Command-line options",
    "pal.plots": "Plot types",
    "pal.calc": "Calculators",
    "pal.empty": "No matches",
    "pal.nav": "navigate",
    "pal.open": "open / copy",
    "pal.close": "close",
    "toast.copied": "Copied to clipboard",
    "com.authors": "Show all {n} authors",
    "com.authors.less": "Show less",
  };

  const GK = (window.GK = {
    $, $$, store, reduced,
    root: document.body.dataset.root || "",
    lang: store.get("gk-lang") || "en",
    t(key) {
      const zh = (window.GK_I18N && window.GK_I18N.zh) || {};
      if (GK.lang === "zh" && zh[key] != null) return zh[key];
      return EN[key] != null ? EN[key] : (GK._en[key] != null ? GK._en[key] : key);
    },
    _en: {},
    on(evt, fn) { document.addEventListener("gk:" + evt, fn); },
    emit(evt, detail) { document.dispatchEvent(new CustomEvent("gk:" + evt, { detail })); },
    docHref(key) {
      const d = window.GK_DOCS && window.GK_DOCS[key];
      return GK.root + (d ? (GK.lang === "zh" ? d[1] : d[0]) : "manual.html");
    },
    esc(s) { return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;"); },
  });

  /* ------------------------------------------------------------ theme */
  function setTheme(t) {
    document.documentElement.setAttribute("data-theme", t);
    store.set("gk-theme", t);
    GK.emit("theme", t);
  }
  $$("[data-theme-toggle]").forEach(b => b.addEventListener("click", () => {
    setTheme(document.documentElement.getAttribute("data-theme") === "dark" ? "light" : "dark");
  }));

  /* ------------------------------------------------------------ split-word headings */
  // Latin text animates word by word; CJK text phrase by phrase, cut after punctuation,
  // so a Chinese heading only wraps between phrases (e.g. after "，")
  const CJK = /[　-鿿＀-￯]/, CUT = "，。、；：！？）」』";
  function tokens(text) {
    const out = [];
    (text.match(/\s+|\S+/g) || []).forEach(t => {
      if (!CJK.test(t)) { out.push(t); return; }
      let cur = "";
      for (const ch of t) { cur += ch; if (CUT.includes(ch)) { out.push(cur); cur = ""; } }
      if (cur) out.push(cur);
    });
    return out;
  }
  function splitWords(el) {
    const walk = node => {
      Array.from(node.childNodes).forEach(c => {
        if (c.nodeType === 3) {
          const frag = document.createDocumentFragment();
          tokens(c.textContent).forEach(tok => {
            if (/^\s+$/.test(tok)) { frag.appendChild(document.createTextNode(tok)); return; }
            const s = document.createElement("span"); s.className = "w"; s.textContent = tok; frag.appendChild(s);
          });
          node.replaceChild(frag, c);
        } else if (c.nodeType === 1 && !c.classList.contains("w") && c.tagName !== "CODE" && c.tagName !== "BR") {
          walk(c);
        }
      });
    };
    walk(el);
    $$(".w", el).forEach((w, i) => w.style.setProperty("--i", i));
  }

  /* ------------------------------------------------------------ language */
  $$("[data-i18n]").forEach(el => { GK._en[el.dataset.i18n] = el.innerHTML; });
  $$("[data-i18n-ph]").forEach(el => { GK._en[el.dataset.i18nPh] = el.getAttribute("placeholder"); });
  function applyLang(lang) {
    GK.lang = lang;
    store.set("gk-lang", lang);
    document.documentElement.lang = lang === "zh" ? "zh-CN" : "en";
    const zh = (window.GK_I18N && window.GK_I18N.zh) || {};
    $$("[data-i18n]").forEach(el => {
      const k = el.dataset.i18n;
      el.innerHTML = lang === "zh" && zh[k] != null ? zh[k] : GK._en[k];
      if (el.classList.contains("split-words")) splitWords(el);
    });
    $$("[data-i18n-ph]").forEach(el => { const k = el.dataset.i18nPh; el.setAttribute("placeholder", lang === "zh" && zh[k] ? zh[k] : GK._en[k]); });
    $$("[data-doc]").forEach(a => { a.href = GK.docHref(a.dataset.doc); });
    $$("[data-lang-toggle]").forEach(b => { b.textContent = lang === "zh" ? "EN" : "中文"; });
    GK.emit("lang", lang);
  }
  $$("[data-lang-toggle]").forEach(b => b.addEventListener("click", () => applyLang(GK.lang === "zh" ? "en" : "zh")));
  $$(".split-words").forEach(splitWords);
  if (GK.lang !== "en") applyLang(GK.lang); else $$("[data-doc]").forEach(a => { a.href = GK.docHref(a.dataset.doc); });

  /* ------------------------------------------------------------ navigation */
  const nav = $("#nav"), progress = $(".progress");
  function onScroll() {
    const y = window.scrollY;
    if (nav) nav.classList.toggle("scrolled", y > 8);
    if (progress) {
      const h = document.documentElement.scrollHeight - innerHeight;
      progress.style.transform = `scaleX(${h > 0 ? y / h : 0})`;
    }
  }
  addEventListener("scroll", onScroll, { passive: true });
  onScroll();
  const menuBtn = $("#menu-btn"), links = $("#nav-links");
  if (menuBtn && links) {
    menuBtn.addEventListener("click", e => { e.stopPropagation(); links.classList.toggle("open"); });
    document.addEventListener("click", e => { if (!links.contains(e.target)) links.classList.remove("open"); });
    $$("a", links).forEach(a => a.addEventListener("click", () => links.classList.remove("open")));
  }
  // highlight the nav link of the section in view
  const spy = $$('.nav-links a[href^="#"]').map(a => [a, $(a.getAttribute("href"))]).filter(p => p[1]);
  if (spy.length) {
    const io = new IntersectionObserver(es => es.forEach(e => {
      spy.forEach(([a, s]) => { if (s === e.target) a.classList.toggle("active", e.isIntersecting); });
    }), { rootMargin: "-45% 0px -50% 0px" });
    spy.forEach(([, s]) => io.observe(s));
  }

  /* ------------------------------------------------------------ reveal + counters */
  function countUp(el) {
    const to = +el.dataset.count || 0, dur = reduced ? 0 : 1400, t0 = performance.now();
    const step = now => {
      const p = dur ? Math.min(1, (now - t0) / dur) : 1;
      el.textContent = Math.round(to * (1 - Math.pow(1 - p, 4)));
      if (p < 1) requestAnimationFrame(step);
    };
    requestAnimationFrame(step);
  }
  GK.countUp = countUp;
  const rio = new IntersectionObserver(entries => entries.forEach(e => {
    if (!e.isIntersecting) return;
    e.target.classList.add("in");
    $$("[data-count]", e.target).forEach(c => { c._seen = true; countUp(c); });
    rio.unobserve(e.target);
  }), { threshold: 0.12, rootMargin: "0px 0px -6% 0px" });
  GK.observe = el => rio.observe(el);
  $$(".reveal, .split-words").forEach(el => rio.observe(el));

  /* ------------------------------------------------------------ spotlight cards + magnetic buttons */
  document.addEventListener("pointermove", e => {
    const card = e.target.closest && e.target.closest(".card");
    if (card) {
      const r = card.getBoundingClientRect();
      card.style.setProperty("--mx", e.clientX - r.left + "px");
      card.style.setProperty("--my", e.clientY - r.top + "px");
    }
  }, { passive: true });
  if (!reduced && matchMedia("(pointer: fine)").matches) {
    $$(".btn-primary").forEach(b => {
      b.addEventListener("pointermove", e => {
        const r = b.getBoundingClientRect();
        b.style.transform = `translate(${(e.clientX - r.left - r.width / 2) * 0.12}px, ${(e.clientY - r.top - r.height / 2) * 0.22}px)`;
      });
      b.addEventListener("pointerleave", () => { b.style.transform = ""; });
    });
  }

  /* ------------------------------------------------------------ toast + copy */
  let toastEl, toastTimer;
  GK.toast = (msg) => {
    if (!toastEl) {
      toastEl = document.createElement("div"); toastEl.className = "toast"; toastEl.setAttribute("role", "status");
      document.body.appendChild(toastEl);
    }
    toastEl.innerHTML = `<svg><use href="#i-check"/></svg><span>${GK.esc(msg)}</span>`;
    toastEl.classList.add("show");
    clearTimeout(toastTimer); toastTimer = setTimeout(() => toastEl.classList.remove("show"), 1800);
  };
  GK.copy = async (text, btn) => {
    try { await navigator.clipboard.writeText(text); }
    catch (e) {
      const ta = document.createElement("textarea"); ta.value = text; ta.style.position = "fixed"; ta.style.opacity = "0";
      document.body.appendChild(ta); ta.select(); try { document.execCommand("copy"); } catch (_) { /* ignore */ } ta.remove();
    }
    GK.toast(GK.t("toast.copied"));
    if (btn) { btn.classList.add("done"); setTimeout(() => btn.classList.remove("done"), 1400); }
  };
  document.addEventListener("click", e => {
    const b = e.target.closest("[data-copy], [data-copy-from], [data-copy-pre]");
    if (!b) return;
    let text = b.dataset.copy;
    if (b.hasAttribute("data-copy-from")) text = ($(b.dataset.copyFrom) || {}).textContent;
    if (b.hasAttribute("data-copy-pre")) text = b.parentElement.querySelector("pre").textContent;
    if (text) GK.copy(text.trim(), b);
  });

  /* ------------------------------------------------------------ hero install tabs */
  $$("[data-install]").forEach(box => {
    const code = $("[data-install-code]", box);
    $$("[data-cmd]", box).forEach(tab => tab.addEventListener("click", () => {
      $$("[data-cmd]", box).forEach(x => x.classList.toggle("active", x === tab));
      code.textContent = tab.dataset.cmd;
    }));
  });

  /* ------------------------------------------------------------ lightbox */
  let lb, lbItems = [], lbIdx = 0;
  function lbRender() {
    const it = lbItems[lbIdx];
    if (!it) return;
    $(".lb-img img", lb).src = it.src;
    $(".lb-img img", lb).alt = it.title;
    $(".lb-bar b", lb).textContent = it.title;
    const cmd = $(".lb-bar .cmd", lb);
    cmd.style.display = it.command ? "" : "none";
    $("code", cmd).textContent = it.command || "";
    $(".copy-btn", cmd).dataset.copy = it.command || "";
  }
  GK.lightbox = (items, idx = 0) => {
    if (!lb) {
      lb = document.createElement("div"); lb.className = "lightbox"; lb.setAttribute("role", "dialog"); lb.setAttribute("aria-modal", "true");
      lb.innerHTML = `<div class="lb-inner"><button class="lb-close" aria-label="Close"><svg><use href="#i-close"/></svg></button>
        <div class="lb-img"><img alt=""></div>
        <div class="lb-bar"><b></b><div class="cmd"><code></code><button class="copy-btn" aria-label="Copy"><svg class="cp"><use href="#i-copy"/></svg><svg class="ok"><use href="#i-check"/></svg></button></div></div>
        <button class="lb-nav lb-prev" aria-label="Previous"><svg><use href="#i-left"/></svg></button><button class="lb-nav lb-next" aria-label="Next"><svg><use href="#i-right"/></svg></button></div>`;
      document.body.appendChild(lb);
      lb.addEventListener("click", e => { if (e.target === lb || e.target.closest(".lb-close")) GK.lightbox.close(); });
      $(".lb-prev", lb).addEventListener("click", () => { lbIdx = (lbIdx - 1 + lbItems.length) % lbItems.length; lbRender(); });
      $(".lb-next", lb).addEventListener("click", () => { lbIdx = (lbIdx + 1) % lbItems.length; lbRender(); });
    }
    lbItems = items; lbIdx = idx; lbRender();
    lb.classList.add("open"); document.body.style.overflow = "hidden";
  };
  GK.lightbox.close = () => { if (lb) { lb.classList.remove("open"); document.body.style.overflow = ""; } };
  document.addEventListener("keydown", e => {
    if (!lb || !lb.classList.contains("open")) return;
    if (e.key === "Escape") GK.lightbox.close();
    if (e.key === "ArrowLeft") $(".lb-prev", lb).click();
    if (e.key === "ArrowRight") $(".lb-next", lb).click();
  });

  /* ------------------------------------------------------------ data helpers */
  const jsonCache = {};
  GK.json = url => (jsonCache[url] = jsonCache[url] || fetch(url).then(r => { if (!r.ok) throw new Error(r.status); return r.json(); }));

  /* ------------------------------------------------------------ site numbers (stats.json) */
  // every [data-stat="key"] shows stats.json[key]; the number written in the HTML is the fallback
  GK.stats = {};
  function fillStats() {
    $$("[data-stat]").forEach(el => {
      const v = GK.stats[el.dataset.stat];
      if (v == null) return;
      if (el.hasAttribute("data-count")) { el.dataset.count = v; if (el._seen) countUp(el); }
      else (el.querySelector(".w") || el).textContent = v;   // .w: word span inside an animated heading
    });
  }
  GK.json(GK.root + "stats.json").then(s => { GK.stats = s || {}; fillStats(); GK.emit("stats", GK.stats); }).catch(() => {});
  GK.on("lang", fillStats);

  /* ------------------------------------------------------------ GitHub stats (cached for 6 h) */
  (function gh() {
    const fmt = n => n >= 1000 ? (n / 1000).toFixed(1).replace(/\.0$/, "") + "k" : String(n);
    const fill = d => {
      $$("[data-gh-stars]").forEach(e => { e.textContent = "★ " + fmt(d.stars); if (e.closest(".gh-stats")) e.textContent = fmt(d.stars); });
      $$("[data-gh-forks]").forEach(e => { e.textContent = fmt(d.forks); });
    };
    try {
      const c = JSON.parse(store.get("gk-gh") || "null");
      if (c && Date.now() - c.t < 6 * 3600e3) { fill(c); return; }
    } catch (e) { /* ignore */ }
    fetch("https://api.github.com/repos/zhyan0603/GPUMDkit").then(r => r.ok ? r.json() : Promise.reject())
      .then(j => { const d = { stars: j.stargazers_count, forks: j.forks_count, t: Date.now() }; store.set("gk-gh", JSON.stringify(d)); fill(d); })
      .catch(() => { $$("[data-gh-stars]").forEach(e => { if (!e.closest(".gh-stats")) e.textContent = "★"; }); });
  })();

  /* ------------------------------------------------------------ command palette (⌘K) */
  const SECTIONS = [["features", "Features", "功能"], ["terminal", "Live terminal", "在线终端"], ["workflow", "Active-learning loop", "主动学习闭环"],
    ["console", "Web console", "网页控制台"], ["agent", "Agent skills", "Agent 技能"], ["install", "Installation", "安装"],
    ["gallery", "Gallery", "图库"], ["publications", "Publications", "出版物"], ["community", "Community & citation", "社区与引用"]];
  const ICON = { sec: "i-arrow", doc: "i-doc", cli: "i-term", plt: "i-plot", calc: "i-calc" };
  function paletteItems() {
    const D = window.GK_DATA || {}, onHome = document.body.dataset.page === "home", zh = GK.lang === "zh";
    const out = [];
    SECTIONS.forEach(([id, en, z]) => out.push({ g: "pal.sections", k: "sec", title: zh ? z : en, sub: "#" + id, href: (onHome ? "" : GK.root + "index.html") + "#" + id }));
    out.push({ g: "pal.sections", k: "sec", title: zh ? "使用手册" : "User manual", sub: "manual.html", href: GK.root + "manual.html" });
    out.push({ g: "pal.sections", k: "sec", title: zh ? "完整图库" : "Full gallery", sub: "gallery.html", href: GK.root + "gallery.html" });
    out.push({ g: "pal.sections", k: "sec", title: zh ? "全部出版物" : "All publications", sub: "publications.html", href: GK.root + "publications.html" });
    Object.entries(window.GK_DOCS || {}).forEach(([key, d]) => out.push({ g: "pal.docs", k: "doc", title: zh ? d[3] : d[2], sub: zh ? d[1] : d[0], href: GK.docHref(key) }));
    Object.entries(D.flags || {}).forEach(([f, desc]) => out.push({ g: "pal.cli", k: "cli", title: "gpumdkit.sh " + f, code: true, sub: desc, copy: "gpumdkit.sh " + f }));
    Object.entries(D.plots || {}).forEach(([p, desc]) => {
      p.split("/").forEach((name, i) => {
        const full = i === 0 ? name : p.split("/")[0].replace(/\d*$/, "") + name;
        out.push({ g: "pal.plots", k: "plt", title: "gpumdkit.sh -plt " + full, code: true, sub: desc, copy: "gpumdkit.sh -plt " + full });
      });
    });
    Object.entries(D.calcs || {}).forEach(([c, desc]) => out.push({ g: "pal.calc", k: "calc", title: "gpumdkit.sh -calc " + c, code: true, sub: desc, copy: "gpumdkit.sh -calc " + c }));
    return out;
  }
  function score(q, s, from = 0) {
    // substring hits rank first (earlier and at word starts is better); otherwise a subsequence match
    const low = s.toLowerCase(), at = low.indexOf(q, from);
    if (at >= 0) {
      const start = at === 0 || /[\s\-_/.]/.test(low[at - 1]);
      return { sc: 1000 - (at - from) * 2 + (start ? 300 : 0) + (low.length - from === q.length ? 200 : 0) - low.length * 0.1,
        idx: Array.from({ length: q.length }, (_, k) => at + k) };
    }
    let i = 0, sc = 0, run = 0; const idx = [];
    for (let j = from; j < low.length && i < q.length; j++) {
      if (low[j] === q[i]) { run++; sc += 1 + run * 2 + (j === from || /[\s\-_/.]/.test(low[j - 1]) ? 4 : 0); idx.push(j); i++; } else run = 0;
    }
    return i === q.length ? { sc: sc - low.length * 0.05, idx } : null;
  }
  // command titles share the "gpumdkit.sh " prefix; match after it so it cannot dominate
  const skip = t => (t.startsWith("gpumdkit.sh ") ? 12 : 0);
  function hl(text, idx) {
    if (!idx || !idx.length) return GK.esc(text);
    const set = new Set(idx); let out = "";
    for (let i = 0; i < text.length; i++) out += set.has(i) ? `<mark>${GK.esc(text[i])}</mark>` : GK.esc(text[i]);
    return out;
  }
  let pal, palItems = [], palSel = 0, palView = [];
  function palRender() {
    const q = $("input", pal).value.trim().toLowerCase();
    let view;
    if (!q) view = palItems.filter(it => it.k !== "cli" && it.k !== "plt" && it.k !== "calc").concat(palItems.filter(it => it.k === "cli").slice(0, 8));
    else {
      view = palItems.map(it => {
        const mt = score(q, it.title, skip(it.title));
        const ms = mt ? null : score(q, it.sub || "");
        if (!mt && !(ms && ms.sc >= 1000)) return null;   // descriptions only count on substring hits
        return Object.assign({}, it, { m: mt || ms, inTitle: !!mt });
      }).filter(Boolean).sort((a, b) => (b.inTitle - a.inTitle) || (b.m.sc - a.m.sc)).slice(0, 60);
    }
    palView = view; palSel = Math.min(palSel, Math.max(0, view.length - 1));
    const list = $(".pal-list", pal);
    if (!view.length) { list.innerHTML = `<div class="pal-empty">${GK.t("pal.empty")}</div>`; return; }
    let html = "", last = "";
    view.forEach((it, i) => {
      if (it.g !== last && !q) { html += `<div class="pal-group">${GK.t(it.g)}</div>`; last = it.g; }
      const t = it.inTitle ? hl(it.title, it.m && it.m.idx) : GK.esc(it.title);
      html += `<div class="pal-item${i === palSel ? " sel" : ""}" data-i="${i}"><div class="pi"><svg><use href="#${ICON[it.k]}"/></svg></div>
        <div class="pt"><b>${it.code ? `<code>${t}</code>` : t}</b><small>${GK.esc(it.sub || "")}</small></div><span class="pk">${it.copy ? "⧉" : "↵"}</span></div>`;
    });
    list.innerHTML = html;
    const sel = $(".pal-item.sel", list); if (sel) sel.scrollIntoView({ block: "nearest" });
  }
  function palRun(it) {
    if (!it) return;
    closePalette();
    if (it.copy) GK.copy(it.copy);
    else if (it.href) { if (it.href.startsWith("#")) { const t = $(it.href); if (t) t.scrollIntoView({ behavior: reduced ? "auto" : "smooth" }); } else location.href = it.href; }
  }
  function openPalette() {
    if (!pal) {
      pal = document.createElement("div"); pal.className = "palette"; pal.setAttribute("role", "dialog"); pal.setAttribute("aria-modal", "true");
      pal.innerHTML = `<div class="pal-box"><div class="pal-in"><svg><use href="#i-search"/></svg><input type="text" autocomplete="off" spellcheck="false"><kbd>Esc</kbd></div>
        <div class="pal-list"></div><div class="pal-foot"><span><kbd>↑</kbd><kbd>↓</kbd><i data-k="pal.nav"></i></span><span><kbd>↵</kbd><i data-k="pal.open"></i></span><span><kbd>Esc</kbd><i data-k="pal.close"></i></span></div></div>`;
      document.body.appendChild(pal);
      pal.addEventListener("click", e => { if (e.target === pal) closePalette(); const it = e.target.closest(".pal-item"); if (it) palRun(palView[+it.dataset.i]); });
      pal.addEventListener("mousemove", e => { const it = e.target.closest(".pal-item"); if (it && +it.dataset.i !== palSel) { palSel = +it.dataset.i; $$(".pal-item", pal).forEach(x => x.classList.toggle("sel", x === it)); } });
      const inp = $("input", pal);
      inp.addEventListener("input", () => { palSel = 0; palRender(); });
      inp.addEventListener("keydown", e => {
        if (e.key === "ArrowDown") { e.preventDefault(); palSel = Math.min(palView.length - 1, palSel + 1); palRender(); }
        else if (e.key === "ArrowUp") { e.preventDefault(); palSel = Math.max(0, palSel - 1); palRender(); }
        else if (e.key === "Enter") { e.preventDefault(); palRun(palView[palSel]); }
        else if (e.key === "Escape") closePalette();
      });
    }
    palItems = paletteItems();
    $("input", pal).placeholder = GK.t("pal.ph");
    $$("[data-k]", pal).forEach(i => { i.textContent = " " + GK.t(i.dataset.k); });
    $("input", pal).value = ""; palSel = 0; palRender();
    pal.classList.add("open");
    setTimeout(() => $("input", pal).focus(), 30);
  }
  function closePalette() { if (pal) pal.classList.remove("open"); }
  GK.palette = openPalette;
  $$("[data-palette]").forEach(b => b.addEventListener("click", openPalette));
  document.addEventListener("keydown", e => {
    const typing = /INPUT|TEXTAREA/.test(document.activeElement.tagName) || document.activeElement.isContentEditable;
    if ((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === "k") { e.preventDefault(); pal && pal.classList.contains("open") ? closePalette() : openPalette(); }
    else if (e.key === "/" && !typing) { e.preventDefault(); openPalette(); }
  });
})();
