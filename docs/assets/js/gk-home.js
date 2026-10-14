/* GPUMDkit website — home page: particle logo, feature visuals, live terminal,
   active-learning loop, web console, gallery rows, publications and citation. */
(function () {
  "use strict";
  const GK = window.GK, $ = GK.$, $$ = GK.$$, D = window.GK_DATA || {};
  const isDark = () => document.documentElement.getAttribute("data-theme") === "dark";
  const palette = () => isDark() ? { blue: [76, 155, 232], copper: [226, 164, 104] } : { blue: [24, 103, 174], copper: [193, 136, 82] };
  function rng(seed) { let s = seed >>> 0; return () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; }; }
  const sleep = ms => new Promise(r => setTimeout(r, ms));
  // run `start` while the element is on screen and `stop` when it leaves
  function whileVisible(el, start, stop) {
    if (!el) return;
    new IntersectionObserver(es => es.forEach(e => (e.isIntersecting ? start() : stop && stop())), { threshold: 0.05 }).observe(el);
  }

  /* ================================================================ hero particle logo */
  (function heroParticles() {
    const cv = $("#hero-canvas");
    if (!cv) return;
    const ctx = cv.getContext("2d");
    let W = 0, H = 0, dpr = 1, P = [], bonds = [], step = 8, raf = 0, running = false, t0 = performance.now();
    const mouse = { x: -1e4, y: -1e4, on: false };

    function logoPoints(S) {
      // draw the official logo geometry (measured from the PNG) and sample its pixels
      const oc = document.createElement("canvas"); oc.width = oc.height = S;
      const g = oc.getContext("2d");
      const k = (S * 0.82) / 540, ox = S / 2 - 400 * k, oy = S / 2 - 293 * k;
      g.setTransform(k, 0, 0, k, ox, oy);
      g.lineWidth = 38; g.strokeStyle = "#0000ff";
      g.beginPath(); g.arc(400, 293, 246, -28 * Math.PI / 180, (-28 - 342) * Math.PI / 180, true); g.stroke();
      g.font = '800 1000px "Helvetica Neue", Helvetica, Arial, sans-serif';
      const m = g.measureText("G"), gw = m.actualBoundingBoxLeft + m.actualBoundingBoxRight, gh = m.actualBoundingBoxAscent + m.actualBoundingBoxDescent;
      const s = (424 - 167) / gh, sx = (511 - 278) / (gw * s);
      g.save(); g.translate(278 + m.actualBoundingBoxLeft * s * sx, 167 + m.actualBoundingBoxAscent * s); g.scale(s * sx, s);
      g.fillStyle = "#0000ff"; g.fillText("G", 0, 0); g.restore();
      g.fillStyle = "#ff0000"; g.beginPath(); g.arc(400 + 246 * Math.cos(21 * Math.PI / 180), 293 - 246 * Math.sin(21 * Math.PI / 180), 30, 0, 2 * Math.PI); g.fill();
      const data = g.getImageData(0, 0, S, S).data, pts = [];
      const step = Math.max(6, Math.round(S / 74)), J = rng(3);
      for (let y = 0; y < S; y += step) for (let x = (y / step) % 2 ? step / 2 : 0; x < S; x += step) {
        const i = (Math.round(y) * S + Math.round(x)) * 4;
        // hexagonal grid plus a little jitter so it reads as atoms, not pixels
        if (data[i + 3] > 140) pts.push([x + (J() - 0.5) * step * 0.45, y + (J() - 0.5) * step * 0.45, data[i] > data[i + 2] ? 1 : 0]);
      }
      pts.step = step;
      return pts;
    }
    function build() {
      const r = cv.getBoundingClientRect();
      dpr = Math.min(2, window.devicePixelRatio || 1);
      W = r.width; H = r.height;
      cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
      const S = Math.round(Math.min(W, H));
      const pts = logoPoints(S), R = rng(7);
      const offX = (W - S) / 2, offY = (H - S) / 2;
      P = pts.map(([x, y, c]) => {
        const a = R() * Math.PI * 2, d = S * (0.55 + R() * 0.6);
        return { hx: x + offX, hy: y + offY, x: W / 2 + Math.cos(a) * d, y: H / 2 + Math.sin(a) * d, vx: 0, vy: 0, c,
          s: 1.7 + R() * 1.3 + (c ? 0.4 : 0), ph: R() * 6.28, delay: R() * 900 };
      });
      // bonds between neighbouring atoms (by their home positions)
      step = pts.step;
      const lim = (step * 1.32) ** 2;
      bonds = [];
      for (let i = 0; i < P.length; i++) for (let j = i + 1; j < P.length; j++) {
        const dx = P[i].hx - P[j].hx, dy = P[i].hy - P[j].hy;
        if (dx * dx + dy * dy < lim) bonds.push(i, j);
      }
      t0 = performance.now();
    }
    function frame(now) {
      const el = now - t0, pal = palette(), t = now / 1000;
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      ctx.clearRect(0, 0, W, H);
      const R2 = 110 * 110;
      // physics
      for (const p of P) {
        const ready = el > p.delay;
        const tx = p.hx + Math.sin(t * 1.6 + p.ph) * 0.9, ty = p.hy + Math.cos(t * 1.3 + p.ph) * 0.9;
        const k = ready ? 0.045 : 0.004;
        p.vx += (tx - p.x) * k; p.vy += (ty - p.y) * k;
        if (mouse.on) {
          const dx = p.x - mouse.x, dy = p.y - mouse.y, d2 = dx * dx + dy * dy;
          if (d2 < R2) { const d = Math.sqrt(d2) || 1, f = Math.pow(1 - d / 110, 2) * 5.2; p.vx += (dx / d) * f; p.vy += (dy / d) * f; }
        }
        p.vx *= 0.82; p.vy *= 0.82; p.x += p.vx; p.y += p.vy;
      }
      // threads from the cursor to nearby particles
      if (mouse.on) {
        ctx.lineWidth = 1;
        let n = 0;
        for (const p of P) {
          const dx = p.x - mouse.x, dy = p.y - mouse.y, d2 = dx * dx + dy * dy;
          if (d2 < 150 * 150 && n++ < 70) {
            ctx.strokeStyle = `rgba(${pal.blue.join(",")},${0.22 * (1 - Math.sqrt(d2) / 150)})`;
            ctx.beginPath(); ctx.moveTo(mouse.x, mouse.y); ctx.lineTo(p.x, p.y); ctx.stroke();
          }
        }
      }
      // bonds fade in once the logo has formed and snap when stretched
      const form = Math.min(1, Math.max(0, (el - 900) / 900));
      if (form > 0) {
        const brk = (step * 2.1) ** 2;
        ctx.lineWidth = 1;
        ctx.strokeStyle = `rgba(${pal.blue.join(",")},${0.2 * form})`;
        ctx.beginPath();
        for (let b = 0; b < bonds.length; b += 2) {
          const p = P[bonds[b]], q = P[bonds[b + 1]], dx = p.x - q.x, dy = p.y - q.y;
          if (dx * dx + dy * dy < brk) { ctx.moveTo(p.x, p.y); ctx.lineTo(q.x, q.y); }
        }
        ctx.stroke();
      }
      // particles, batched by colour; fast ones glow a little
      for (const c of [0, 1]) {
        const col = c ? pal.copper : pal.blue;
        ctx.fillStyle = `rgba(${col.join(",")},0.92)`;
        ctx.beginPath();
        for (const p of P) {
          if (p.c !== c) continue;
          const sp = Math.min(1.8, Math.hypot(p.vx, p.vy) * 0.15);
          ctx.moveTo(p.x + p.s + sp, p.y); ctx.arc(p.x, p.y, p.s + sp, 0, Math.PI * 2);
        }
        ctx.fill();
      }
      if (running) raf = requestAnimationFrame(frame);
    }
    function start() { if (running || GK.reduced) return; running = true; raf = requestAnimationFrame(frame); }
    function stop() { running = false; cancelAnimationFrame(raf); }
    cv.addEventListener("pointermove", e => { const r = cv.getBoundingClientRect(); mouse.x = e.clientX - r.left; mouse.y = e.clientY - r.top; mouse.on = true; });
    cv.addEventListener("pointerleave", () => { mouse.on = false; });
    cv.addEventListener("click", e => {
      const r = cv.getBoundingClientRect(), cx = e.clientX - r.left, cy = e.clientY - r.top;
      for (const p of P) { const dx = p.x - cx, dy = p.y - cy, d = Math.hypot(dx, dy) || 1, f = 26 + Math.random() * 22; p.vx += dx / d * f; p.vy += dy / d * f; }
    });
    let rt; addEventListener("resize", () => { clearTimeout(rt); rt = setTimeout(() => { build(); if (GK.reduced) staticDraw(); }, 150); });
    function staticDraw() { for (const p of P) { p.x = p.hx; p.y = p.hy; } running = true; frame(performance.now()); running = false; }
    build();
    if (GK.reduced) staticDraw();
    document.addEventListener("visibilitychange", () => (document.hidden ? stop() : start()));
    whileVisible(cv, start, stop);
    GK.on("theme", () => { if (GK.reduced) staticDraw(); });
  })();

  /* ================================================================ formats marquee */
  (function marquees() {
    const IN = ["VASP OUTCAR", "XDATCAR", "POSCAR", "CIF", "CP2K log", "ABACUS SCF / MD", "LAMMPS dump", "LAMMPS data", "DeepMD npy", "ASE traj", "MTP cfg", "extxyz"];
    const OUT = ["train.xyz", "test.xyz", "nep.txt", "loss.out", "thermo.out", "msd.out", "kappa.out", "neighbor.out", "POSCAR", "LAMMPS data", "DeepMD npy", "active.xyz"];
    const fill = (el, list, cls) => { if (!el) return; const html = list.map(x => `<span class="fmt ${cls}"><i></i>${x}</span>`).join(""); el.innerHTML = html + html; };
    fill($('[data-marquee="in"]'), IN, "");
    fill($('[data-marquee="out"]'), OUT, "to");
  })();

  /* ================================================================ publications count + latest */
  GK.json("publications.json").then(pubs => {
    $$("[data-pub-count]").forEach(el => {
      el.dataset.count = pubs.length;
      if (el._seen || el.closest(".in")) GK.countUp(el);
      else el.textContent = "0";
    });
    const grid = $("[data-pub-grid]");
    if (grid) {
      grid.innerHTML = pubs.slice(0, 4).map(p => `<a class="card pub reveal" href="${p.link}" target="_blank" rel="noopener">
        <div class="yr"><b>${p.year}</b><svg><use href="#i-out"/></svg></div><h4>${GK.esc(p.title)}</h4><div class="au">${GK.esc(p.authors)}</div>
        <div class="jr">${GK.esc(p.journal)}</div></a>`).join("");
      $$(".reveal", grid).forEach((el, i) => { el.style.setProperty("--d", i * 0.06 + "s"); GK.observe(el); });
    }
  }).catch(() => { $$("[data-pub-count]").forEach(el => { el.textContent = "—"; }); });

  /* ================================================================ gallery data: feature strip + rows */
  GK.json("gallery.json").then(items => {
    const all = items.map(it => ({ src: it.image, title: it.title, command: it.command }));
    const strip = $("[data-plot-strip]");
    if (strip) {
      const html = items.map((it, i) => `<div class="thumb" data-gi="${i}"><img src="${it.image}" alt="${GK.esc(it.title)}" loading="lazy"><span>${GK.esc(it.command)}</span></div>`).join("");
      strip.innerHTML = html + html;
    }
    const rows = [items.filter((_, i) => i % 2 === 0), items.filter((_, i) => i % 2 === 1)];
    rows.forEach((row, r) => {
      const el = $(`[data-g-row="${r}"]`);
      if (!el) return;
      const html = row.map(it => `<div class="g-card" data-gi="${items.indexOf(it)}"><div class="im"><img src="${it.image}" alt="${GK.esc(it.title)}" loading="lazy"></div>
        <div class="meta"><b>${GK.esc(it.title)}</b><code>${GK.esc(it.command)}</code></div></div>`).join("");
      el.innerHTML = html + html;
    });
    document.addEventListener("click", e => {
      const c = e.target.closest("[data-gi]");
      if (c) GK.lightbox(all, +c.dataset.gi);
    });
    GK._gallery = all;
  }).catch(() => {});

  /* ================================================================ feature visuals */
  (function convViz() {
    const track = $("[data-conv-files]"), cnt = $("[data-conv-count]"), body = $("[data-xyz-body]"), box = $("[data-xyz]");
    if (!track) return;
    const files = ["OUTCAR", "XDATCAR", "cp2k.log", "dump.lammpstrj", "POSCAR", "model.cif", "set.000/*.npy", "md.traj", "train.cfg", "running_scf.log"];
    const html = files.map(f => `<div class="f"><i></i>${f}</div>`).join("");
    track.innerHTML = html + html;
    // live preview of the frames being written, in the header layout used by out2xyz.sh
    const R = rng(29), f4 = x => (x < 0 ? "" : " ") + x.toFixed(4);
    let n = 2600, cfg = 37, timer = 0;
    function frame() {
      const a = (10.86 + (R() - 0.5) * 0.04).toFixed(2), e = -(1046.2 + R() * 3.5);
      const atom = (el, x) => `<div><span class="el">${el}</span>${f4(x + R() * 0.06)} ${f4(x + R() * 0.06)} ${f4(x + R() * 0.06)}  ${f4((R() - 0.5) * 0.6)} ${f4((R() - 0.5) * 0.6)} ${f4((R() - 0.5) * 0.6)}</div>`;
      body.innerHTML = `<div class="v">192</div><div class="cm"><span class="k">Config_type</span>=iter01_${cfg} <span class="k">Weight</span>=1.0 ` +
        `<span class="k">Lattice</span>="${a} 0 0 0 ${a} 0 0 0 ${a}" <span class="k">Energy</span>=<span class="v">${e.toFixed(4)}</span> pbc="T T T" ` +
        `<span class="k">Properties</span>=species:S:1:pos:R:3:forces:R:3</div>` + atom("Si", 0) + atom("O", 1.3575);
    }
    function tick() {
      n += 1 + Math.floor(R() * 4); if (R() < 0.3) cfg++;
      if (n > 9999) { n = 2600; cfg = 37; }
      cnt.textContent = n.toLocaleString("en-US");
      frame();
      box.classList.remove("flash"); void box.offsetWidth; box.classList.add("flash");
    }
    cnt.textContent = n.toLocaleString("en-US"); frame();
    whileVisible(track, () => { if (!timer && !GK.reduced) timer = setInterval(tick, 900); }, () => { clearInterval(timer); timer = 0; });
  })();

  (function fpsViz() {
    const cv = $("[data-fps]"), lbl = $("[data-fps-n]");
    if (!cv) return;
    const ctx = cv.getContext("2d"), R = rng(11), pts = [];
    const centers = [[0.22, 0.3, 0.1], [0.62, 0.26, 0.08], [0.44, 0.7, 0.12], [0.8, 0.66, 0.07], [0.15, 0.78, 0.05]];
    const gauss = () => { let u = 0; for (let k = 0; k < 6; k++) u += R(); return u / 3 - 1; };
    for (let i = 0; i < 260; i++) { const c = centers[Math.floor(R() * centers.length)]; pts.push([Math.min(0.96, Math.max(0.04, c[0] + gauss() * c[2] * 1.5)), Math.min(0.92, Math.max(0.08, c[1] + gauss() * c[2] * 1.5))]); }
    const order = [0], dmin = pts.map(() => Infinity);
    for (let k = 0; k < 32; k++) {
      const last = pts[order[order.length - 1]]; let best = -1, bi = 0;
      pts.forEach((p, i) => { const d = (p[0] - last[0]) ** 2 + (p[1] - last[1]) ** 2; dmin[i] = Math.min(dmin[i], d); if (dmin[i] > best) { best = dmin[i]; bi = i; } });
      order.push(bi);
    }
    let raf = 0, run = false, t0 = 0, w = 0, h = 0, dpr = 1;
    const size = () => { const r = cv.getBoundingClientRect(); dpr = Math.min(2, devicePixelRatio || 1); w = r.width; h = r.height; cv.width = w * dpr; cv.height = h * dpr; };
    const draw = now => {
      const el = (now - t0) / 1000, cyc = el % 9, n = Math.min(order.length, Math.floor(Math.max(0, cyc - 0.6) / 0.18));
      const pal = palette();
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, w, h);
      ctx.fillStyle = isDark() ? "rgba(160,170,185,0.35)" : "rgba(120,130,145,0.38)";
      ctx.beginPath(); for (const [x, y] of pts) { ctx.moveTo(x * w + 2.4, y * h); ctx.arc(x * w, y * h, 2.4, 0, 6.283); } ctx.fill();
      for (let k = 0; k < n; k++) {
        const [x, y] = pts[order[k]], age = cyc - 0.6 - k * 0.18;
        if (age < 0.5) { ctx.strokeStyle = `rgba(${pal.blue.join(",")},${0.6 * (1 - age / 0.5)})`; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.arc(x * w, y * h, 5 + age * 40, 0, 6.283); ctx.stroke(); }
        ctx.fillStyle = `rgb(${(k === n - 1 ? pal.copper : pal.blue).join(",")})`;
        ctx.beginPath(); ctx.arc(x * w, y * h, 4.4, 0, 6.283); ctx.fill();
      }
      lbl.textContent = n;
      if (run) raf = requestAnimationFrame(draw);
    };
    size(); addEventListener("resize", size);
    whileVisible(cv, () => { if (run) return; run = true; t0 = performance.now(); raf = requestAnimationFrame(draw); }, () => { run = false; cancelAnimationFrame(raf); });
  })();

  (function folders() {
    const el = $("[data-folders]"); if (!el) return;
    el.innerHTML = Array.from({ length: 24 }, (_, i) => `<i class="fd${i >= 20 ? " cu" : ""}" style="animation-delay:${(i % 8 + Math.floor(i / 8)) * 0.09}s"></i>`).join("");
  })();
  (function bars() {
    const el = $("[data-bars]"); if (!el) return;
    const R = rng(5);
    el.innerHTML = Array.from({ length: 16 }, (_, i) => { const hgt = 30 + 70 * Math.exp(-Math.pow((i - 7) / 4.2, 2)) * (0.8 + 0.4 * R()); return `<i style="height:${hgt}%;animation-delay:${-R() * 3}s;animation-duration:${2.4 + R() * 1.6}s"></i>`; }).join("");
  })();

  /* ================================================================ live terminal */
  (function terminal() {
    const body = $("#term-body");
    if (!body) return;
    const esc = GK.esc;
    const PROMPT = '<span class="tp-u">user@gpumdkit</span> <span class="tp-d">~/project</span> <span class="tp-s">$</span> ';
    const out = document.createElement("div"); body.appendChild(out);
    const inLine = document.createElement("div"); inLine.className = "term-in";
    inLine.innerHTML = `<span class="pfx">${PROMPT}</span><input type="text" spellcheck="false" autocomplete="off" autocapitalize="off" aria-label="Terminal input">`;
    body.appendChild(inLine);
    const input = $("input", inLine), pfx = $(".pfx", inLine);
    // abort: set while a side-panel command interrupts the running one; runId: latest side-panel click
    let mode = "shell", busy = false, hist = [], hIdx = 0, touched = false, abort = false, runId = 0;

    const hlCmd = s => esc(s).replace(/^(\S+)/, '<span class="t-cmd">$1</span>').replace(/(\s)(-[\w-]+)/g, '$1<span class="t-flag">$2</span>');
    function hl(line) {
      let s = esc(line);
      if (/^\s*[+|>]?[-=+]{6,}[+<]?\s*$/.test(line) || /^\s*-{10,}\s*$/.test(line)) return `<span class="t-box">${s}</span>`;
      if (/Code path:/.test(line)) return `<span class="t-dim">${s}</span>`;
      if (/^\s*[_|\\/ ()<>]*[_\\/|][_|\\/ ()<>]*$/.test(line) && line.trim().length > 20) return `<span class="t-art">${s}</span>`;
      s = s.replace(/^(\s*)\|/, '$1<span class="t-box">|</span>').replace(/\|(\s*)$/, '<span class="t-box">|</span>$1');
      s = s.replace(/(\s)(\d{1,3}\))/g, '$1<span class="t-num">$2</span>').replace(/^(\s*)(\d\))/, '$1<span class="t-num">$2</span>');
      s = s.replace(/(\s)(-[a-z][\w-]*)/g, '$1<span class="t-flag">$2</span>');
      s = s.replace(/(#{3,})/g, '<span class="t-num">$1</span>');
      s = s.replace(/(Conversion complete\.|Converged OUTCAR files: \d+)/g, '<span class="t-ok">$1</span>');
      return s;
    }
    const scroll = () => { body.scrollTop = body.scrollHeight; };
    function line(html) { const d = document.createElement("div"); d.className = "term-line"; d.innerHTML = html || " "; out.appendChild(d); return d; }
    async function print(lines, opt = {}) {
      const delay = opt.delay != null ? opt.delay : (lines.length > 30 ? 4 : 12);
      for (const l of lines) { if (abort) return; line(opt.raw ? l : hl(l)); scroll(); if (delay) await sleep(delay); }
    }
    function setMode(m) { mode = m; pfx.innerHTML = m === "shell" ? PROMPT : " "; }
    function echo(cmd) { line((mode === "shell" ? PROMPT : " ") + (mode === "shell" ? hlCmd(cmd) : esc(cmd))); }
    const note = s => [`<span class="t-note"> ${esc(s)}</span>`];

    const galleryFor = type => {
      const items = GK._gallery || [];
      const alias = { test: "prediction", sigma: "D/sigma", D: "D/sigma" };
      const want = "-plt " + (alias[type] || type);
      return items.findIndex(it => it.command && it.command.split(" ").slice(1).join(" ") === want);
    };
    async function timeDemo(kind) {
      const pad = (s, n) => (s + " ".repeat(n)).slice(0, n);
      const ctr = (s, n) => { const l = Math.floor((n - s.length) / 2); return " ".repeat(Math.max(0, l)) + s + " ".repeat(Math.max(0, n - s.length - l)); };
      if (kind === "nep" || kind === "gnep") {
        await print(["+-----------------+-----------+-----------------+---------------------+", "|       Step      | Time Diff |    Time Left    |    Finish Time      |", "+-----------------+-----------+-----------------+---------------------+"], { delay: 20 });
        const t0 = Date.now();
        for (let i = 0; i < 6; i++) {
          await sleep(420); if (abort) return;
          const step = 61000 + i * 100, left = Math.round((100000 - step) / 100 * 37.5);
          const fin = new Date(t0 + left * 1000).toISOString().replace("T", " ").slice(0, 19);
          await print([`| ${pad(String(step), 15)} | ${pad((37 + (i % 3 === 1 ? 1 : 0)) + " s", 9)} | ${pad(`${Math.floor(left / 3600)} h ${Math.floor(left % 3600 / 60)} m ${left % 60} s`, 15)} | ${pad(fin, 19)} |`], { delay: 0 });
        }
      } else {
        await print([" ----------------- System Information ----------------", " total frames: 2000000", " -----------------------------------------------------",
          [ctr("Current Frame", 15), ctr("Speed (steps/s)", 15), ctr("Total Time", 15), ctr("Time Left", 15), ctr("Estimated End", 20)].join(" ")], { delay: 20 });
        const t0 = Date.now(), f = x => `${Math.floor(x / 3600)}h ${Math.floor(x % 3600 / 60)}m ${x % 60}s`;
        for (let i = 0; i < 6; i++) {
          await sleep(420); if (abort) return;
          const fr = 520000 + i * 5300, sp = 5310 + ((i * 37) % 23) - 11, left = Math.round((2000000 - fr) / sp);
          const end = new Date(t0 + left * 1000).toISOString().replace("T", " ").slice(0, 19);
          await print([[ctr(String(fr), 15), ctr(sp.toFixed(2), 15), ctr(f(Math.round(2000000 / sp)), 15), ctr(f(left), 15), ctr(end, 20)].join(" ")], { delay: 0 });
        }
      }
      await print(["^C", ...note("(demo values · press Ctrl+C to stop monitoring in a real terminal)")], { raw: true, delay: 30 });
    }
    async function out2xyz() {
      const L = D.out2xyz || [];
      for (const l of L) {
        if (abort) return;
        if (/Progress:/.test(l)) {
          const d = line("");
          for (let k = 0; k <= 50 && !abort; k += 5) { d.innerHTML = hl(` Progress: [${"#".repeat(k)}${" ".repeat(50 - k)}] ${String(k * 2).padStart(3)}% (${k === 50 ? 1 : 0}/1)`); scroll(); await sleep(45); }
        } else { line(hl(l)); scroll(); await sleep(60); }
      }
    }
    async function run(cmd) {
      const c = cmd.trim();
      if (mode === "menu") return menuInput(c);
      if (mode === "sub") return subInput(c);
      if (!c) return;
      const tok = c.split(/\s+/), t0 = tok[0];
      if (t0 === "clear") { out.innerHTML = ""; return; }
      if (t0 === "help") return print(["Demo commands:", "  gpumdkit.sh              interactive menu (then type 1-8, 0 to exit)", "  gpumdkit.sh -h           all command-line options",
        "  gpumdkit.sh -plt -h      list plot types", "  gpumdkit.sh -plt <type>  e.g. train, test, thermo, msd, rdf, emd", "  gpumdkit.sh -calc        calculator tools",
        "  gpumdkit.sh -out2xyz .   OUTCAR -> extxyz", "  gpumdkit.sh -time nep    training ETA", "  gpumdkit.sh -skill       agent skill info", "  ls, pwd, clear"], { raw: false, delay: 10 });
      if (t0 === "ls") return print(["OUTCAR  POSCAR  loss.out  nep.in  nep.txt  run.in  thermo.out  train.xyz"], { delay: 0 });
      if (t0 === "pwd") return print(["/home/user/project"], { delay: 0 });
      if (t0 === "whoami") return print(["user"], { delay: 0 });
      if (t0 === "cd") return;
      if (t0 === "echo") return print([tok.slice(1).join(" ")], { delay: 0 });
      if (t0 !== "gpumdkit.sh" && t0 !== "./gpumdkit.sh") return print([`zsh: command not found: ${t0}`], { delay: 0 });

      const f = tok[1];
      if (!f) { await print(D.main || [], { delay: 10 }); setMode("menu"); return; }
      if (f === "-h" || f === "-help") return print(D.help || []);
      if (f === "-plt") {
        const ty = tok[2];
        if (!ty || ty === "-h") return print(D.plt || []);
        const gi = galleryFor(ty);
        if (gi >= 0) {
          const it = GK._gallery[gi];
          await print([`<span class="t-dim"> Running plt_scripts for "${esc(ty)}" …</span>`], { raw: true, delay: 160 });
          if (abort) return;
          const d = line(""); d.innerHTML = `<span class="t-fig" data-gi="${gi}"><img src="${it.src}" alt="${esc(it.title)}"></span>`;
          const img = $("img", d); img.addEventListener("load", scroll); scroll();
          return print([`<span class="t-ok"> ✓</span> ${esc(it.title)} <span class="t-dim">(click the figure to enlarge)</span>`], { raw: true, delay: 0 });
        }
        if (D.plots && Object.keys(D.plots).some(k => k.split("/").includes(ty) || k === ty || (k.startsWith("thermo2") && /^thermo[23]$/.test(ty))))
          return print(note(`(web demo) "${ty}" plots your own output files; run it in a folder with the data. See gpumdkit.sh -plt -h.`), { raw: true });
        await print([` Unknown plot type: ${ty}`, " Available types are listed below."], { delay: 20 });
        return print(D.plt || []);
      }
      if (f === "-calc") {
        if (!tok[2] || tok[2] === "-h") return print(D.calc || []);
        if (D.calcs && D.calcs[tok[2]]) return print(note(`(web demo) -calc ${tok[2]}: ${D.calcs[tok[2]]}. It needs your input files; run it locally.`), { raw: true });
        return print([" See the codes in calculators folder for more details", " Code path: ~/GPUMDkit/Scripts/calculators"]);
      }
      if (f === "-out2xyz") {
        if (!tok[2] || tok[2] === "-h") return print([" Usage: -out2xyz dir_name (eg. gpumdkit.sh -out2xyz .)", " See the source code of out2xyz.sh for more details", " Code path: ~/GPUMDkit/Scripts/format_conversion/out2xyz.sh"]);
        return out2xyz();
      }
      if (f === "-time") {
        const m = tok[2];
        if (["nep", "gnep", "gpumd", "wpe"].includes(m)) return timeDemo(m === "wpe" ? "gpumd" : m);
        return print([" Supported modes: gpumd, wpe, nep or gnep.", " See the codes in analyzer folder for more details", " Code path: ~/GPUMDkit/Scripts/analyzer/time_consuming_*.sh"]);
      }
      if (f === "-skill") return print((D.skill || []).slice(0, 34));
      if (f === "-server") {
        await print(note("(web demo) -server starts the GPUMDkit web console for this directory. Scroll down to see it."), { raw: true });
        return;
      }
      if (D.flags && D.flags[f]) return print(note(`(web demo) gpumdkit.sh ${f}: ${D.flags[f].replace(/\s{2,}/g, " ")}. Add -h for its usage when running locally.`), { raw: true });
      return print([` Unknown: ${f}`, " Use -h for help, or create ~/.gpumdkit.in to define custom commands"]);
    }
    async function menuInput(c) {
      if (D.menus && D.menus[c]) {
        const m = D.menus[c];
        await print(m.lines, { delay: 8 });
        if (m.waits) { setMode("sub"); mode = "sub"; return; }
        await print(D.citation || [], { delay: 8 }); setMode("shell"); return;
      }
      if (c === "8") { await print([" ------------>>"].concat(D.help || []), { delay: 4 }); await print(D.citation || [], { delay: 8 }); setMode("shell"); return; }
      if (c === "0") { await print(D.citation || [], { delay: 8 }); setMode("shell"); return; }
      return print([" ------------>>", " Please reinput function number:"], { delay: 10 });
    }
    async function subInput(c) {
      if (c === "000") { await print(D.main || [], { delay: 8 }); setMode("menu"); return; }
      if (/^\d{3}$/.test(c) || /^[a-z][\w]+$/.test(c)) {
        await print(note(`(web demo) Function ${c} would now ask for its input files. Run gpumdkit.sh locally to continue.`), { raw: true });
        await print(D.citation || [], { delay: 8 }); setMode("shell"); return;
      }
      return print([" ------------>>", " Please reinput function number..."], { delay: 10 });
    }

    function complete() {
      const v = input.value, words = v.split(" ");
      let cands = [], prefix = "";
      if (mode !== "shell") return;
      if (words.length === 1) cands = ["gpumdkit.sh", "clear", "help", "ls", "pwd"].filter(x => x.startsWith(words[0]));
      else if (words[0] === "gpumdkit.sh") {
        const cur = words[words.length - 1], prev = words[words.length - 2];
        if (words.length === 2) cands = Object.keys(D.flags || {}).concat(["-help"]).filter(x => x.startsWith(cur));
        else if (prev === "-plt") cands = Object.keys(D.plots || {}).flatMap(k => k === "thermo2/3" ? ["thermo2", "thermo3"] : [k]).concat(["test"]).filter(x => x.startsWith(cur));
        else if (prev === "-calc") cands = Object.keys(D.calcs || {}).filter(x => x.startsWith(cur));
        else if (prev === "-time") cands = ["gpumd", "wpe", "nep", "gnep"].filter(x => x.startsWith(cur));
        prefix = words.slice(0, -1).join(" ") + " ";
      }
      cands = Array.from(new Set(cands)).sort();
      if (!cands.length) return;
      if (cands.length === 1) { input.value = (prefix || "") + cands[0] + " "; return; }
      let common = cands[0];
      for (const x of cands) while (!x.startsWith(common)) common = common.slice(0, -1);
      const cur = words[words.length - 1];
      if (common.length > cur.length) { input.value = (prefix || "") + common; return; }
      line(PROMPT + hlCmd(v)); line(`<span class="t-dim">${cands.map(esc).join("  ")}</span>`); scroll();
    }
    async function submit(cmd) {
      if (busy) return;
      busy = true; inLine.style.visibility = "hidden";
      echo(cmd);
      if (cmd.trim() && mode === "shell") { hist.push(cmd); hIdx = hist.length; }
      try { await run(cmd); } catch (e) { line(`<span class="t-err">${esc(String(e))}</span>`); }
      busy = false; inLine.style.visibility = ""; input.value = ""; scroll();
    }
    input.addEventListener("keydown", e => {
      touched = true;
      if (e.key === "Enter") { e.preventDefault(); submit(input.value); }
      else if (e.key === "Tab") { e.preventDefault(); complete(); }
      else if (e.key === "ArrowUp") { e.preventDefault(); if (hIdx > 0) { hIdx--; input.value = hist[hIdx]; } }
      else if (e.key === "ArrowDown") { e.preventDefault(); if (hIdx < hist.length - 1) { hIdx++; input.value = hist[hIdx]; } else { hIdx = hist.length; input.value = ""; } }
      else if (e.key === "l" && e.ctrlKey) { e.preventDefault(); out.innerHTML = ""; }
      else if (e.key === "c" && e.ctrlKey) { e.preventDefault(); line((mode === "shell" ? PROMPT : " ") + esc(input.value) + "^C"); input.value = ""; setMode("shell"); scroll(); }
    });
    body.addEventListener("click", e => { if (!e.target.closest(".t-fig") && !getSelection().toString()) input.focus({ preventScroll: true }); });

    // side-panel commands always start from a clean shell prompt: a running demo or an open
    // gpumdkit.sh menu is left with Ctrl+C first, as you would in a real terminal
    async function typeAndRun(cmd, speed = 24) {
      const my = ++runId;
      touched = true;
      let cut = false;
      if (busy) { abort = true; while (busy) await sleep(16); abort = false; cut = true; }
      if (my !== runId) return;
      if (cut || mode !== "shell") { line((mode === "shell" ? "" : " ") + esc(input.value) + "^C"); setMode("shell"); }
      input.value = "";
      for (const ch of cmd) { if (my !== runId) return; input.value += ch; scroll(); await sleep(speed); }
      await sleep(120);
      if (my !== runId || busy) return;
      submit(cmd);
    }
    // suggestion list
    const SUG = [["gpumdkit.sh", "interactive menu", "交互菜单"], ["gpumdkit.sh -h", "every command-line option", "全部命令行选项"],
      ["gpumdkit.sh -plt train", "plot NEP training", "绘制 NEP 训练结果"], ["gpumdkit.sh -plt -h", "list the {plots} plot types", "列出 {plots} 种绘图"],
      ["gpumdkit.sh -calc", "calculator tools", "计算器工具"], ["gpumdkit.sh -out2xyz .", "OUTCAR → extxyz", "OUTCAR → extxyz"],
      ["gpumdkit.sh -time nep", "live training ETA", "实时训练剩余时间"], ["gpumdkit.sh -skill", "agent skill info", "Agent 技能信息"]];
    const side = $("[data-term-cmds]"), term = $("#term");
    let picked = null;
    const renderSide = () => {
      side.innerHTML = SUG.map(([c, en, zh]) => `<button data-run="${esc(c)}"${c === picked ? ' class="on"' : ""}><div><code>${esc(c)}</code><small>${(GK.lang === "zh" ? zh : en).replace("{plots}", GK.stats.plot_types || 42)}</small></div><span class="ar">↵</span></button>`).join("");
    };
    const pick = c => { picked = c; $$("[data-run]", side).forEach(b => b.classList.toggle("on", b.dataset.run === c)); };
    renderSide(); GK.on("lang", renderSide); GK.on("stats", renderSide);
    side.addEventListener("click", e => {
      const b = e.target.closest("[data-run]");
      if (!b) return;
      pick(b.dataset.run);
      // on narrow screens the list sits below the terminal: bring the terminal back into view
      const r = term.getBoundingClientRect();
      if (r.bottom < 140 || r.top > innerHeight - 140) term.scrollIntoView({ behavior: GK.reduced ? "auto" : "smooth", block: "start" });
      typeAndRun(b.dataset.run);
    });
    input.addEventListener("input", () => pick(null));

    // greet, then auto-play the interactive menu once when the terminal first scrolls into view
    line('<span class="t-dim">Welcome to the GPUMDkit web demo. Type <span class="t-flag">help</span> for a list of commands.</span>');
    let played = false;
    whileVisible($("#term"), async () => {
      if (played) return; played = true;
      await sleep(500);
      if (touched) return;
      for (const ch of "gpumdkit.sh") { if (touched) return; input.value += ch; await sleep(70); }
      await sleep(200); if (touched) { input.value = ""; return; }
      await submit("gpumdkit.sh");
      await sleep(900); if (touched) return;
      input.value = "4"; await sleep(350); if (touched) return;
      await submit("4");
    });
  })();

  /* ================================================================ active-learning loop */
  (function loop() {
    const svg = $("[data-loop-svg]");
    if (!svg) return;
    const STEPS = [
      { ic: "i-db", en: ["Label with DFT", "Batch-prepare VASP or CP2K single-point jobs from sampled structures, one folder per frame, with shared inputs linked."],
        zh: ["DFT 标注", "将采样得到的结构批量准备为 VASP 或 CP2K 单点计算，每帧一个文件夹，并自动链接共享输入文件。"], cmds: ["gpumdkit.sh → 3 → 301"] },
      { ic: "i-conv", en: ["Collect the data", "Turn finished OUTCARs or CP2K logs into an extxyz training set, then check its composition and energy, force and virial ranges."],
        zh: ["整理数据", "把计算完成的 OUTCAR 或 CP2K 日志转换为 extxyz 训练集，并检查成分以及能量、力和维里的分布范围。"], cmds: ["-out2xyz .", "-cp2k2xyz", "-analyze_comp", "-range"] },
      { ic: "i-brain", en: ["Train NEP", "Run nep, follow the loss and the ETA live, and inspect parity plots once it converges."],
        zh: ["训练 NEP", "运行 nep，实时查看 loss 与剩余时间，收敛后检查 parity 图。"], cmds: ["-time nep", "-plt train", "-plt test"] },
      { ic: "i-atom", en: ["Explore with MD", "Prepare GPUMD or LAMMPS sampling runs with the new potential across conditions, and keep an eye on their progress."],
        zh: ["MD 探索", "用新势函数准备 GPUMD 或 LAMMPS 在不同条件下的采样模拟，并随时掌握运行进度。"], cmds: ["gpumdkit.sh → 3 → 302", "-time gpumd", "-plt thermo"] },
      { ic: "i-filter", en: ["Select new structures", "Pick the most informative frames by farthest-point sampling or by maximum force deviation, then send them back to DFT."],
        zh: ["筛选新结构", "通过最远点采样或最大力偏差挑选最有价值的结构，再送回 DFT 标注，开始下一轮。"], cmds: ["gpumdkit.sh → 2 → 203", "gpumdkit.sh → 2 → 205", "-filter_dist_pbc"] },
    ];
    const C = 280, R = 205, N = STEPS.length, ns = "http://www.w3.org/2000/svg";
    const ang = i => (-90 + i * 360 / N) * Math.PI / 180;
    const pos = i => [C + R * Math.cos(ang(i)), C + R * Math.sin(ang(i))];
    let html = `<circle class="loop-ring" cx="${C}" cy="${C}" r="${R}"/><path class="loop-arc" data-arc/>
      <circle r="7" fill="var(--copper)" data-pulse style="filter:drop-shadow(0 0 8px var(--copper))"/>`;
    STEPS.forEach((s, i) => {
      const [x, y] = pos(i);
      html += `<g class="loop-node" data-node="${i}" transform="translate(${x} ${y})" tabindex="0" role="button">
        <circle class="bg" r="44"/><svg class="ic" x="-15" y="-15" width="30" height="30"><use href="#${s.ic}"/></svg>
        <text y="70" data-node-label>${s.en[0]}</text><text class="n" y="-56">0${i + 1}</text></g>`;
    });
    svg.innerHTML = html;
    const arc = $("[data-arc]", svg), pulse = $("[data-pulse]", svg);
    // all five steps sit in one grid cell, so the card always has the height of the longest step
    const panes = $("[data-loop-panes]"), bars = $("[data-loop-steps]");
    bars.innerHTML = STEPS.map((_, i) => `<button aria-label="Step ${i + 1}" data-bar="${i}"><i></i></button>`).join("");
    let cur = 0, pausedUntil = 0, startT = performance.now();
    function arcPath(a0, a1) {
      const [x0, y0] = [C + R * Math.cos(a0), C + R * Math.sin(a0)], [x1, y1] = [C + R * Math.cos(a1), C + R * Math.sin(a1)];
      return `M ${x0} ${y0} A ${R} ${R} 0 ${a1 - a0 > Math.PI ? 1 : 0} 1 ${x1} ${y1}`;
    }
    function renderPanes() {
      const zh = GK.lang === "zh", esc = GK.esc;
      panes.innerHTML = STEPS.map((s, i) => {
        const L = zh ? s.zh : s.en;
        return `<div class="loop-pane" data-pane="${i}"><div class="step">${zh ? `第 ${i + 1} 步 / 共 5 步` : `Step ${i + 1} / 5`}</div>
          <h3>${esc(L[0])}</h3><p>${esc(L[1])}</p>
          <div class="cmdrow">${s.cmds.map(c => `<span class="chip">${esc(c)}</span>`).join("")}<a class="chip" href="${GK.docHref("active_learning_workflow")}">${zh ? "工作流文档 →" : "Workflow docs →"}</a></div></div>`;
      }).join("");
    }
    function select(i) {
      cur = i;
      $$("[data-node]", svg).forEach((n, k) => n.classList.toggle("on", k === i));
      $$("[data-pane]", panes).forEach((p, k) => { p.classList.toggle("on", k === i); p.setAttribute("aria-hidden", String(k !== i)); });
      $$("[data-bar]", bars).forEach((b, k) => {
        const fill = $("i", b);
        fill.style.transition = "none"; fill.style.transform = `scaleX(${k < i ? 1 : 0})`;
        if (k === i) { void fill.offsetWidth; fill.style.transition = "transform 3.4s linear"; fill.style.transform = "scaleX(1)"; }
      });
      startT = performance.now();
    }
    function relabel() { $$("[data-node-label]", svg).forEach((t, i) => { t.textContent = (GK.lang === "zh" ? STEPS[i].zh : STEPS[i].en)[0]; }); renderPanes(); select(cur); }
    GK.on("lang", relabel);
    svg.addEventListener("click", e => { const n = e.target.closest("[data-node]"); if (n) { pausedUntil = performance.now() + 9000; select(+n.dataset.node); } });
    svg.addEventListener("keydown", e => { if (e.key === "Enter") { const n = e.target.closest("[data-node]"); if (n) select(+n.dataset.node); } });
    bars.addEventListener("click", e => { const b = e.target.closest("[data-bar]"); if (b) { pausedUntil = performance.now() + 9000; select(+b.dataset.bar); } });
    let raf = 0, running = false;
    function frame(now) {
      if (now > pausedUntil && now - startT > 3400) select((cur + 1) % N);
      // arc grows from the previous node to the current one; a pulse travels along it
      const p = Math.min(1, (now - startT) / 900), a0 = ang((cur + N - 1) % N), a1base = ang(cur);
      const a1 = a1base < a0 ? a1base + Math.PI * 2 : a1base;
      const e = 1 - Math.pow(1 - p, 3), aa = a0 + (a1 - a0) * e;
      arc.setAttribute("d", arcPath(a0, aa));
      pulse.setAttribute("cx", C + R * Math.cos(aa)); pulse.setAttribute("cy", C + R * Math.sin(aa));
      if (running) raf = requestAnimationFrame(frame);
    }
    select(0); relabel();
    whileVisible(svg, () => { if (running) return; running = true; startT = performance.now(); raf = requestAnimationFrame(frame); }, () => { running = false; cancelAnimationFrame(raf); });
  })();

  /* ================================================================ web console showcase */
  (function consoleShow() {
    const stage = $(".console-stage"), browser = $("[data-tilt]");
    if (!stage) return;
    const imgs = $$("[data-console-img]"), tabs = $$("[data-console-tab]"), url = $("[data-console-url]");
    const URLS = ["127.0.0.1:8888 / neptrain2", "127.0.0.1:8888 / neptrain2 / diffusion / 300K"];
    let cur = 0, timer = 0;
    function show(i) {
      cur = i;
      imgs.forEach((im, k) => { im.classList.toggle("off", k !== i); im.classList.remove("pan"); });
      tabs.forEach((t, k) => t.classList.toggle("active", k === i));
      url.textContent = URLS[i];
      requestAnimationFrame(() => requestAnimationFrame(() => imgs[i].classList.add("pan")));
    }
    tabs.forEach(t => t.addEventListener("click", () => { show(+t.dataset.consoleTab); clearInterval(timer); timer = setInterval(() => show((cur + 1) % 2), 7000); }));
    whileVisible(stage, () => { show(cur); if (!timer) timer = setInterval(() => show((cur + 1) % 2), 7000); }, () => { clearInterval(timer); timer = 0; });
    if (!GK.reduced && matchMedia("(pointer: fine)").matches) {
      stage.addEventListener("pointermove", e => {
        const r = stage.getBoundingClientRect(), x = (e.clientX - r.left) / r.width - 0.5, y = (e.clientY - r.top) / r.height - 0.5;
        browser.style.transform = `rotateY(${x * 8}deg) rotateX(${-y * 6}deg)`;
      });
      stage.addEventListener("pointerleave", () => { browser.style.transform = ""; });
    }
  })();

  /* ================================================================ citation */
  (function cite() {
    const tog = $("[data-authors-toggle]"), au = $("[data-authors]"), view = $("[data-bib]");
    if (!au) return;
    // "Zihan Yan" -> "Yan, Zihan"
    const names = au.textContent.split(",").map(s => s.trim()).filter(Boolean).map(n => { const p = n.split(" "); return p.pop() + ", " + p.join(" "); });
    if (tog) {
      // the label counts the authors in the list; the toggle only shows while the list is cut off
      const label = () => { tog.textContent = GK.t(au.classList.contains("open") ? "com.authors.less" : "com.authors").replace("{n}", names.length); };
      const fit = () => { tog.hidden = !au.classList.contains("open") && au.scrollHeight <= au.clientHeight + 2; };
      tog.addEventListener("click", () => { au.classList.toggle("open"); label(); });
      label(); fit();
      addEventListener("resize", fit); GK.on("lang", () => { label(); fit(); });
      if (document.fonts) document.fonts.ready.then(fit);
    }
    const fields = [["title", "GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP"], ["author", names.join(" and\n             ")], ["journal", "MGE Advances"],
      ["volume", "4"], ["pages", "e70074"], ["year", "2026"], ["doi", "10.1002/mgea.70074"]];
    const bib = "@article{yan2026gpumdkit,\n" + fields.map(([k, v]) => `  ${k.padEnd(7)} = {${v}}`).join(",\n") + "\n}";
    // on-page preview: one line per field, the author list shortened to first two … last
    if (view) {
      const short = names.length > 3 ? `${names[0]} and ${names[1]} and … and ${names[names.length - 1]}` : names.join(" and ");
      view.innerHTML = `<div><span class="bk">@article</span>{yan2026gpumdkit,</div>` +
        fields.map(([k, v], i) => `<div>  <span class="bf">${k.padEnd(7)}</span> = {<span class="bv">${GK.esc(k === "author" ? short : v)}</span>}${i < fields.length - 1 ? "," : ""}</div>`).join("") + "<div>}</div>";
    }
    $$("[data-copy-bib]").forEach(b => b.addEventListener("click", () => GK.copy(bib, b.classList.contains("copy-btn") ? b : null)));
  })();
})();
