function el(id) { return document.getElementById(id); }
function esc(s) {
  return String(s)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}
function joinRel(base, name) { return base ? base + "/" + name : name; }
function fmtSize(n) {
  if (n == null || n < 0) return "-";
  if (n < 1024) return n + " B";
  if (n < 1048576) return (n / 1024).toFixed(1) + " KB";
  return (n / 1048576).toFixed(1) + " MB";
}
function fmtAgo(sec) {
  if (sec == null) return "";
  if (sec < 5) return "now";
  if (sec < 60) return Math.floor(sec) + "s";
  if (sec < 3600) return Math.floor(sec / 60) + "m";
  if (sec < 86400) return Math.floor(sec / 3600) + "h";
  return Math.floor(sec / 86400) + "d";
}
function fmtDur(sec) {
  if (sec == null || !isFinite(sec)) return "--";
  if (sec < 60) return Math.floor(sec) + "s";
  if (sec < 3600) return Math.floor(sec / 60) + "m " + Math.floor(sec % 60) + "s";
  return Math.floor(sec / 3600) + "h " + Math.floor((sec % 3600) / 60) + "m";
}
function fmtRate(bps) {
  if (bps == null || !isFinite(bps) || bps <= 0) return "";
  var perMin = bps * 60;
  if (perMin < 1024) return perMin.toFixed(0) + " B/min";
  if (perMin < 1048576) return (perMin / 1024).toFixed(1) + " KB/min";
  return (perMin / 1048576).toFixed(2) + " MB/min";
}
function fmtGenRate(gps) {
  if (gps == null || !isFinite(gps) || gps <= 0) return "--";
  if (gps >= 1) return gps.toFixed(1) + " gen/s";
  if (gps * 60 >= 1) return (gps * 60).toFixed(1) + " gen/min";
  return (gps * 3600).toFixed(1) + " gen/h";
}

var curPath = "";
var rootDir = "";
var running = false;
var cmdHistory = [];
var lastFiles = null;
var lastSubdirs = null;
var lastNep = null;
var lastPollTs = 0;
var viewer = { name: "", rel: "", text: null, mode: "text" };
var lossText = null;
var toastTimer = null;
var pollTimer = null;
var pollActive = false;
var scanSeq = 0;
var termBusy = false;

function toast(msg) {
  var t = el("toast");
  t.textContent = msg;
  t.classList.add("show");
  if (toastTimer) clearTimeout(toastTimer);
  toastTimer = setTimeout(function() { t.classList.remove("show"); }, 2500);
}

function showLogin(msg) {
  el("loginView").classList.remove("hidden");
  el("appView").classList.add("hidden");
  if (msg) el("loginMsg").textContent = msg;
}
function hideLogin() {
  el("loginView").classList.add("hidden");
  el("appView").classList.remove("hidden");
  el("loginMsg").textContent = "";
}

async function api(url, opts) {
  var resp = await fetch(url, opts);
  if (resp.status === 401) showLogin("Please sign in to continue.");
  return resp;
}

async function doLogin() {
  var pw = el("pwInput").value;
  if (!pw) return;
  el("loginBtn").disabled = true;
  try {
    var resp = await fetch("/api/login", {
      method: "POST",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify({password: pw})
    });
    var data = {};
    try { data = await resp.json(); } catch (e) {}
    if (resp.ok) {
      el("pwInput").value = "";
      hideLogin();
      loadScan(curPath);
    } else {
      el("loginMsg").textContent = data.error || "Login failed.";
    }
  } catch (e) {
    el("loginMsg").textContent = "Network error: " + e;
  }
  el("loginBtn").disabled = false;
}

async function loadScan(path, silent) {
  var seq = ++scanSeq;
  var resp;
  try {
    resp = await api("/api/scan?path=" + encodeURIComponent(path));
  } catch (e) {
    if (!silent) el("sbErr").textContent = "network error";
    return false;
  }
  if (seq !== scanSeq) return false;
  if (resp.status === 401) return false;
  var data = null;
  try { data = await resp.json(); } catch (e) {}
  if (seq !== scanSeq) return false;
  if (!resp.ok || !data) {
    if (!silent) el("sbErr").textContent = (data && data.error) || "scan failed";
    return false;
  }
  hideLogin();
  curPath = data.path || "";
  rootDir = data.root || "";
  var nowS = Math.floor(Date.now() / 1000);
  var pollGap = lastPollTs ? (Date.now() / 1000 - lastPollTs) : 0;
  var prevF = {}, prevS = {};
  if (lastFiles) lastFiles.forEach(function(f) { prevF[f.name] = f; });
  if (lastSubdirs) lastSubdirs.forEach(function(s) { prevS[s.name] = s; });
  var prevNepDone = lastNep ? lastNep.done : null;
  lastFiles = data.files || [];
  lastSubdirs = data.subdirs || [];
  lastNep = data.training ? {done: data.training.done} : null;
  lastPollTs = Date.now() / 1000;
  var liveFiles = (data.files || []).filter(function(f) {
    return prevF[f.name] && f.size > prevF[f.name].size;
  }).map(function(f) {
    var grew = f.size - prevF[f.name].size;
    return {name: f.name, grew: grew, rate: pollGap > 0 ? grew / pollGap : null};
  });
  renderCrumb(curPath);
  renderDirList(data.subdirs || [], prevS, nowS);
  renderFileList(data.files || [], prevF);
  var liveCount = liveFiles.length;
  el("sbErr").textContent = "";
  el("sbRoot").textContent = rootDir;
  el("sbRoot").title = rootDir;
  el("sbCwd").textContent = "/" + curPath;
  el("sbLive").textContent = liveCount ? liveCount + " live" : "";
  el("sbEntries").textContent = ((data.dirs || []).length + (data.files || []).length) + " entries";
  var srvNow = data.now || nowS;
  var fresh = false;
  var filesArr = data.files || [];
  for (var fi = 0; fi < filesArr.length; fi++) {
    if (filesArr[fi].mtime != null && srvNow - filesArr[fi].mtime < 180) { fresh = true; break; }
  }
  if (!fresh) {
    var subsArr = data.subdirs || [];
    for (var si = 0; si < subsArr.length; si++) {
      if (subsArr[si].newest != null && srvNow - subsArr[si].newest < 180) { fresh = true; break; }
    }
  }
  pollActive = fresh || (data.training && !data.training.finished);
  renderMonitor(data, {liveFiles: liveFiles, pollGap: pollGap, prevNepDone: prevNepDone});
  if (!running) renderRecs(data.recommendations || []);
  renderFigures(data.files || []);
  el("termCwdLabel").textContent = "/" + curPath;
  el("termPathLabel").textContent = "/" + curPath;
  return true;
}

function renderCrumb(path) {
  var c = el("crumb");
  c.innerHTML = "";
  c.appendChild(crumbLink("root", ""));
  if (!path) return;
  var parts = path.split("/");
  var acc = "";
  for (var i = 0; i < parts.length; i++) {
    acc = acc ? acc + "/" + parts[i] : parts[i];
    var sep = document.createElement("span");
    sep.className = "sep";
    sep.textContent = " / ";
    c.appendChild(sep);
    if (i === parts.length - 1) {
      var here = document.createElement("span");
      here.className = "here";
      here.textContent = parts[i];
      c.appendChild(here);
    } else {
      c.appendChild(crumbLink(parts[i], acc));
    }
  }
}

function crumbLink(label, rel) {
  var a = document.createElement("a");
  a.href = "#";
  a.textContent = label;
  a.onclick = function(ev) { ev.preventDefault(); loadScan(rel); };
  return a;
}

function renderDirList(subs, prevS, nowS) {
  var dl = el("dirList");
  dl.innerHTML = "";
  if (!subs.length) {
    dl.innerHTML = '<p class="muted" style="padding:4px 8px;margin:0">(none)</p>';
    return;
  }
  subs.forEach(function(s) {
    var row = document.createElement("div");
    row.className = "drow";
    var dot = document.createElement("span");
    dot.className = "dot";
    var prev = prevS[s.name];
    if (prev && s.newest != null && prev.newest != null && s.newest > prev.newest) dot.classList.add("live");
    row.appendChild(dot);
    var nm = document.createElement("span");
    nm.className = "dname";
    nm.textContent = s.name + "/";
    row.appendChild(nm);
    var meta = document.createElement("span");
    meta.className = "dmeta";
    var bits = [s.count + " items"];
    if (s.newest != null) bits.push(fmtAgo(nowS - s.newest));
    if (s.loss) bits.push("NEP");
    meta.textContent = bits.join(" · ");
    row.appendChild(meta);
    row.onclick = function() { loadScan(joinRel(curPath, s.name)); };
    dl.appendChild(row);
  });
}

function renderFileList(files, prevF) {
  var fl = el("fileList");
  fl.innerHTML = "";
  var filter = el("fileFilter").value.trim().toLowerCase();
  var shown = 0;
  files.forEach(function(f) {
    if (filter && f.name.toLowerCase().indexOf(filter) < 0) return;
    shown++;
    var row = document.createElement("div");
    row.className = "frow";
    var dot = document.createElement("span");
    dot.className = "dot";
    var prev = prevF[f.name];
    if (prev && f.size > prev.size) dot.classList.add("live");
    row.appendChild(dot);
    var a = document.createElement("a");
    a.href = "#";
    a.textContent = f.name;
    a.onclick = function(ev) { ev.preventDefault(); fileClicked(f.name); };
    row.appendChild(a);
    var sz = document.createElement("span");
    sz.className = "fsize";
    sz.textContent = fmtSize(f.size);
    row.appendChild(sz);
    fl.appendChild(row);
  });
  if (!shown) fl.innerHTML = '<p class="muted" style="padding:4px 8px;margin:0">(no files)</p>';
}

function renderMonitor(data, ctx) {
  var body = el("monBody");
  var card = el("monCard");
  body.innerHTML = "";
  var training = data.training;
  if (!training) lossText = null;
  var parts = [];
  if (training) {
    var pct = training.total > 0 ? (100 * training.done / training.total) : 0;
    var rate = null;
    if (ctx.prevNepDone != null && ctx.pollGap > 2) {
      var dr = training.done - ctx.prevNepDone;
      if (dr > 0) rate = dr / ctx.pollGap;
    }
    var eta = rate ? (training.total - training.done) / rate : null;
    var div = document.createElement("div");
    div.innerHTML =
      '<div class="stat-row">' +
      '<div class="stat' + (training.finished ? " done" : "") + '"><div class="stat-v">' + pct.toFixed(1) + '%</div><div class="stat-l">progress</div></div>' +
      '<div class="stat"><div class="stat-v">' + training.done.toLocaleString() + ' / ' + training.total.toLocaleString() + '</div><div class="stat-l">generation</div></div>' +
      '<div class="stat"><div class="stat-v">' + fmtGenRate(rate) + '</div><div class="stat-l">rate</div></div>' +
      '<div class="stat"><div class="stat-v">' + (eta != null ? fmtDur(eta) : "--") + '</div><div class="stat-l">eta</div></div>' +
      '</div>' +
      '<div class="nep-bar-track"><div class="nep-bar-fill" style="width:' + Math.min(100, pct) + '%"></div></div>' +
      '<div class="chart-wrap"><canvas class="chart" id="lossCanvas" height="220"></canvas><span class="pulse-dot hidden" id="lossDot"></span></div>' +
      '<div class="lg-row" id="lossLegend"></div>' +
      '<p class="muted" style="margin:6px 0 0">loss.out loss functions, log-log scale; the pulsing dot marks the latest record</p>';
    parts.push(div);
  }
  if (ctx.liveFiles.length) {
    var lf = document.createElement("div");
    var rows = ctx.liveFiles.map(function(f) {
      return '<div class="row" style="justify-content:space-between;padding:3px 0">' +
        '<span class="mono" style="font-size:12.5px">' + esc(f.name) + "</span>" +
        '<span class="muted mono">+' + fmtSize(f.grew) + (f.rate ? " · " + fmtRate(f.rate) : "") + "</span></div>";
    }).join("");
    lf.innerHTML = '<div style="font-weight:600;font-size:13px;margin:6px 0 4px">Growing files</div>' + rows;
    parts.push(lf);
  }
  if (!parts.length) {
    card.classList.add("hidden");
    return;
  }
  card.classList.remove("hidden");
  el("monMeta").textContent = training ? (training.finished ? "training done" : "training") : ctx.liveFiles.length + " growing";
  parts.forEach(function(p) { body.appendChild(p); });
  if (training) {
    fetchLoss().then(drawLoss).catch(function() {});
  }
}

function renderRecs(recs) {
  var card = el("recCard");
  var box = el("recBox");
  box.innerHTML = "";
  if (!recs.length) { card.classList.add("hidden"); return; }
  card.classList.remove("hidden");
  el("recCount").textContent = recs.length + " actions";
  var list = document.createElement("div");
  list.className = "rec-list";
  recs.forEach(function(rec) {
    var row = document.createElement("div");
    row.className = "rec";
    var icon = document.createElement("span");
    icon.className = "ric";
    icon.innerHTML = ICON_SVG;
    var c = document.createElement("code");
    c.className = "rcmd";
    c.textContent = rec.command.replace(/^gpumdkit\.sh\s+/, "");
    c.title = (rec.description || "") + "  |  click to copy: " + rec.command;
    c.onclick = function() { copyText(rec.command); };
    var ev = document.createElement("span");
    ev.className = "rev";
    (rec.evidence || []).forEach(function(e) {
      var b = document.createElement("span");
      b.className = "badge ok";
      b.textContent = e;
      ev.appendChild(b);
    });
    var btn = document.createElement("button");
    btn.textContent = "Run";
    btn.onclick = function() { runAction(rec, btn); };
    row.appendChild(icon);
    row.appendChild(c);
    row.appendChild(ev);
    row.appendChild(btn);
    list.appendChild(row);
  });
  box.appendChild(list);
}

var ICON_SVG = '<svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polyline points="3 17 9 11 13 15 21 7"/><polyline points="15 7 21 7 21 13"/></svg>';
var ICON_SUN = '<svg width="17" height="17" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="5"/><line x1="12" y1="1" x2="12" y2="3"/><line x1="12" y1="21" x2="12" y2="23"/><line x1="4.22" y1="4.22" x2="5.64" y2="5.64"/><line x1="18.36" y1="18.36" x2="19.78" y2="19.78"/><line x1="1" y1="12" x2="3" y2="12"/><line x1="21" y1="12" x2="23" y2="12"/><line x1="4.22" y1="19.78" x2="5.64" y2="18.36"/><line x1="18.36" y1="5.64" x2="19.78" y2="4.22"/></svg>';
var ICON_MOON = '<svg width="17" height="17" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"/></svg>';

var themeManual = null;
try { themeManual = localStorage.getItem("gk_theme"); } catch (e) {}
function applyTheme() {
  var dark = themeManual === "dark" ||
    (themeManual !== "light" && window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches);
  document.body.classList.toggle("dark", dark);
  var btn = el("themeBtn");
  if (btn) btn.innerHTML = dark ? ICON_SUN : ICON_MOON;
  drawLoss();
}
function toggleTheme() {
  var dark = document.body.classList.contains("dark");
  themeManual = dark ? "light" : "dark";
  try { localStorage.setItem("gk_theme", themeManual); } catch (e) {}
  applyTheme();
}

function setRailActive(btn) {
  var nodes = document.querySelectorAll(".rail-btn");
  for (var i = 0; i < nodes.length; i++) nodes[i].classList.remove("active");
  btn.classList.add("active");
}
function railScroll(targetId, fallbackMsg) {
  var target = el(targetId);
  if (target && !target.classList.contains("hidden")) {
    target.scrollIntoView({behavior: "smooth", block: "start"});
  } else {
    toast(fallbackMsg);
  }
}

function copyText(text) {
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(text).then(function() {
      toast("Copied to clipboard");
    }, function() {
      toast(fallbackCopy(text) ? "Copied to clipboard" : "Copy failed");
    });
  } else {
    toast(fallbackCopy(text) ? "Copied to clipboard" : "Copy failed");
  }
}
function fallbackCopy(text) {
  var ta = document.createElement("textarea");
  ta.value = text;
  ta.style.position = "fixed";
  ta.style.opacity = "0";
  document.body.appendChild(ta);
  ta.select();
  var ok = false;
  try { ok = document.execCommand("copy"); } catch (e) {}
  document.body.removeChild(ta);
  return ok;
}

function fileClicked(name) {
  var dot = name.lastIndexOf(".");
  var ext = dot >= 0 ? name.substring(dot + 1).toLowerCase() : "";
  if (ext === "png" || ext === "jpg" || ext === "jpeg") openLightbox(name);
  else viewFile(name);
}

function pollTick() {
  if (document.hidden) { schedule(); return; }
  loadScan(curPath, true).then(schedule);
}
function schedule() {
  if (pollTimer) clearTimeout(pollTimer);
  pollTimer = null;
  var sel = parseInt(el("refreshSel").value, 10);
  if (sel <= 0) { el("sbRefresh").textContent = "auto off"; return; }
  var delay = pollActive ? sel * 1000 : 60000;
  el("sbRefresh").textContent = pollActive ? "auto " + sel + "s" : "idle · 60s check";
  pollTimer = setTimeout(pollTick, delay);
}
function startPolling() {
  schedule();
}
async function runAction(rec, btn) {
  if (running) return;
  running = true;
  btn.disabled = true;
  btn.textContent = "Running...";
  el("outCard").classList.remove("hidden");
  el("output").textContent = "$ " + rec.command + "\n\n(running, please wait...)";
  el("outCard").scrollIntoView({behavior: "smooth", block: "nearest"});
  var t0 = Date.now();
  try {
    var resp = await api("/api/run", {
      method: "POST",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify({action: rec.action, path: curPath})
    });
    var data = null;
    try { data = await resp.json(); } catch (e) {}
    if (!resp.ok || !data) {
      el("output").textContent = (data && data.error) || "The command could not be started.";
      addHistory(rec.command, -1, 0);
    } else {
      var text = "$ " + data.command + "\n\n" + (data.output || "");
      if (data.returncode !== 0) text += "\n[exit code " + data.returncode + "]";
      el("output").textContent = text;
      addHistory(data.command, data.returncode, (Date.now() - t0) / 1000);
      if (data.image) openLightbox(data.image);
      if (data.returncode === 0) toast(data.image ? "Done: " + data.image : "Done");
    }
  } catch (e) {
    el("output").textContent = "Network error: " + e;
  }
  btn.disabled = false;
  btn.textContent = "Run";
  running = false;
  loadScan(curPath, true);
}

function addHistory(cmd, rc, dur) {
  cmdHistory.unshift({cmd: cmd, rc: rc, dur: dur, at: new Date()});
  if (cmdHistory.length > 30) cmdHistory.pop();
  renderHistory();
}
function renderHistory() {
  var card = el("histCard");
  var list = el("histList");
  if (!cmdHistory.length) { card.classList.add("hidden"); return; }
  card.classList.remove("hidden");
  list.innerHTML = "";
  cmdHistory.forEach(function(h) {
    var d = document.createElement("div");
    d.className = "hist-item";
    var rc = h.rc === 0
      ? '<span class="badge ok">ok</span>'
      : '<span class="badge" style="background:#fdeaea;color:#b03030">' + (h.rc < 0 ? "err" : "exit " + h.rc) + "</span>";
    d.innerHTML = rc + '<span class="hcmd">' + esc(h.cmd) + "</span>" +
      '<span class="hmeta">' + fmtDur(h.dur) + " · " + h.at.toLocaleTimeString() + "</span>";
    list.appendChild(d);
  });
}

function arrMinMax(a) {
  var mn = a[0], mx = a[0];
  for (var i = 1; i < a.length; i++) {
    if (a[i] < mn) mn = a[i];
    if (a[i] > mx) mx = a[i];
  }
  return [mn, mx];
}
function fmtNum(v) {
  if (Math.abs(v) >= 1e5 || (Math.abs(v) < 1e-3 && v !== 0)) return v.toExponential(1);
  return String(Math.round(v * 1000) / 1000);
}
function fmtAxis(v) {
  if (v === 0) return "0";
  var e = Math.round(Math.log10(Math.abs(v)));
  if (Math.abs(v / Math.pow(10, e) - 1) < 1e-9) return "1e" + e;
  if (Math.abs(v) >= 1e4 || Math.abs(v) < 1e-2) return v.toExponential(1).replace("e+", "e");
  return String(Math.round(v * 100) / 100);
}
function hexToRgba(hex, alpha) {
  var m = /^#?([0-9a-f]{6})$/i.exec(hex);
  if (!m) return "rgba(31,119,180," + alpha + ")";
  var v = parseInt(m[1], 16);
  return "rgba(" + ((v >> 16) & 255) + "," + ((v >> 8) & 255) + "," + (v & 255) + "," + alpha + ")";
}
function axisTicks(lo, hi, isLog, count) {
  var ticks = [];
  if (isLog) {
    var k0 = Math.ceil(lo - 1e-9), k1 = Math.floor(hi + 1e-9);
    if (k1 - k0 >= 0 && k1 - k0 <= 8) {
      for (var k = k0; k <= k1; k++) ticks.push({pos: k, label: "1e" + k});
      return ticks;
    }
  }
  var n = count || 4;
  for (var i = 0; i <= n; i++) {
    var pos = lo + (hi - lo) * i / n;
    var orig = isLog ? Math.pow(10, pos) : pos;
    ticks.push({pos: pos, label: fmtAxis(orig)});
  }
  return ticks;
}
function themeColors() {
  var grid = "#e6ebef", axis = "#8a97a8";
  try {
    var cs = getComputedStyle(document.body);
    var g = cs.getPropertyValue("--grid").trim();
    var a = cs.getPropertyValue("--axis").trim();
    if (g) grid = g;
    if (a) axis = a;
  } catch (e) {}
  return {grid: grid, axis: axis};
}
function drawSeries(canvas, series, opts) {
  opts = opts || {};
  var theme = themeColors();
  var dpr = window.devicePixelRatio || 1;
  var w = canvas.clientWidth || 600;
  var h = canvas.clientHeight || 220;
  canvas.width = w * dpr;
  canvas.height = h * dpr;
  var ctx = canvas.getContext("2d");
  ctx.scale(dpr, dpr);
  ctx.clearRect(0, 0, w, h);
  var flat = [];
  series.forEach(function(s) {
    var pts = [];
    for (var i = 0; i < s.xs.length; i++) {
      var xv = opts.logX ? Math.log10(s.xs[i]) : s.xs[i];
      var yv = opts.logY ? Math.log10(s.ys[i]) : s.ys[i];
      if (!isFinite(xv) || !isFinite(yv)) continue;
      pts.push([xv, yv]);
    }
    s._pts = pts;
    for (var j = 0; j < pts.length; j++) flat.push(pts[j]);
  });
  if (!flat.length) {
    ctx.fillStyle = theme.axis;
    ctx.font = "12px monospace";
    ctx.fillText("(no data)", 12, 24);
    return null;
  }
  var xsAll = flat.map(function(p) { return p[0]; });
  var ysAll = flat.map(function(p) { return p[1]; });
  var xr = arrMinMax(xsAll), yr = arrMinMax(ysAll);
  var xmin = xr[0], xmax = xr[1], ymin = yr[0], ymax = yr[1];
  if (xmax === xmin) xmax = xmin + 1;
  if (ymax === ymin) ymax = ymin + Math.abs(ymin) * 0.1 + 1;
  var pad = {l: 62, r: 14, t: 14, b: 42};
  var pw = w - pad.l - pad.r, ph = h - pad.t - pad.b;
  function X(x) { return pad.l + pw * (x - xmin) / (xmax - xmin); }
  function Y(y) { return pad.t + ph * (1 - (y - ymin) / (ymax - ymin)); }
  ctx.strokeStyle = theme.grid;
  ctx.lineWidth = 1;
  ctx.font = "10px monospace";
  ctx.fillStyle = theme.axis;
  var xTicks = axisTicks(xmin, xmax, opts.logX, 4);
  var yTicks = axisTicks(ymin, ymax, opts.logY, 4);
  yTicks.forEach(function(t) {
    var yy = Y(t.pos);
    ctx.beginPath();
    ctx.moveTo(pad.l, yy);
    ctx.lineTo(w - pad.r, yy);
    ctx.stroke();
    ctx.fillText(t.label, 6, yy + 3);
  });
  xTicks.forEach(function(t) {
    var xx = X(t.pos);
    ctx.beginPath();
    ctx.moveTo(xx, pad.t);
    ctx.lineTo(xx, h - pad.b);
    ctx.stroke();
    ctx.fillText(t.label, xx - 14, h - pad.b + 16);
  });
  if (opts.xLabel) {
    ctx.fillText(opts.xLabel, pad.l + pw / 2 - ctx.measureText(opts.xLabel).width / 2, h - 8);
  }
  if (opts.yLabel) {
    ctx.save();
    ctx.translate(12, pad.t + ph / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText(opts.yLabel, -ctx.measureText(opts.yLabel).width / 2, 0);
    ctx.restore();
  }
  var lastPoints = [];
  series.forEach(function(s) {
    if (!s._pts.length) return;
    ctx.strokeStyle = s.color || "rgb(24, 103, 174)";
    ctx.lineWidth = 2;
    ctx.beginPath();
    for (var i = 0; i < s._pts.length; i++) {
      var px = X(s._pts[i][0]), py = Y(s._pts[i][1]);
      if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
    }
    ctx.stroke();
    if (s._pts.length === 1) {
      ctx.fillStyle = s.color || "rgb(24, 103, 174)";
      ctx.beginPath();
      ctx.arc(X(s._pts[0][0]), Y(s._pts[0][1]), 3, 0, 2 * Math.PI);
      ctx.fill();
    }
    var lp = s._pts[s._pts.length - 1];
    lastPoints.push({px: X(lp[0]), py: Y(lp[1]), color: s.color || "rgb(24, 103, 174)"});
  });
  return lastPoints;
}

async function fetchLoss() {
  var rel = joinRel(curPath, "loss.out");
  var resp = await api("/api/file?path=" + encodeURIComponent(rel));
  if (resp.ok) {
    lossText = await resp.text();
  } else {
    lossText = null;
  }
}
var LOSS_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"];

function drawLoss() {
  var canvas = el("lossCanvas");
  if (!canvas || !lossText) return;
  var parsed = parseNumeric(lossText);
  if (!parsed) return;
  var rows = parsed.rows;
  var nSeries, labels;
  if (parsed.cols >= 7) {
    nSeries = 6;
    labels = ["Total", "L1-Reg", "L2-Reg", "Energy-train", "Force-train", "Virial-train"];
  } else if (parsed.cols === 6) {
    nSeries = 4;
    labels = ["Loss", "Energy-train", "Force-train", "Virial-train"];
  } else {
    nSeries = parsed.cols - 1;
    labels = [];
    for (var ci = 1; ci <= nSeries; ci++) labels.push("c" + ci);
  }
  var series = [];
  for (var c = 1; c <= nSeries; c++) {
    var xs = [], ys = [];
    for (var r = 0; r < rows.length; r++) {
      xs.push(r + 1);
      ys.push(rows[r][c]);
    }
    series.push({xs: xs, ys: ys, color: LOSS_COLORS[(c - 1) % LOSS_COLORS.length], label: labels[c - 1]});
  }
  var step = rows.length > 1 ? (rows[1][0] - rows[0][0]) : null;
  var xlabel = step === 100 ? "Generation/100" : (step === 1 ? "Epoch" : "record #");
  var last = drawSeries(canvas, series, {logX: true, logY: true, xLabel: xlabel, yLabel: "Loss functions"});
  renderLossLegend(labels);
  placePulseDot(last);
}

function renderLossLegend(labels) {
  var legend = el("lossLegend");
  if (!legend) return;
  legend.innerHTML = "";
  labels.forEach(function(lb, i) {
    var item = document.createElement("span");
    item.className = "lg-item";
    var dot = document.createElement("span");
    dot.className = "lg-dot";
    dot.style.background = LOSS_COLORS[i % LOSS_COLORS.length];
    var tx = document.createElement("span");
    tx.textContent = lb;
    item.appendChild(dot);
    item.appendChild(tx);
    legend.appendChild(item);
  });
}

function placePulseDot(lastPoints) {
  var dot = el("lossDot");
  if (!dot) return;
  if (!lastPoints || !lastPoints.length) {
    dot.classList.add("hidden");
    return;
  }
  var p = lastPoints[0];
  dot.style.left = p.px + "px";
  dot.style.top = p.py + "px";
  dot.style.background = p.color;
  dot.style.setProperty("--pc", hexToRgba(p.color, 0.55));
  dot.classList.remove("hidden");
}

function parseNumeric(text) {
  var rows = [];
  var cols = 0;
  var lines = text.split(/\r?\n/);
  for (var i = 0; i < lines.length; i++) {
    var line = lines[i].trim();
    if (!line || line.charAt(0) === "#") continue;
    var parts = line.split(/\s+/);
    var nums = [];
    var ok = true;
    for (var j = 0; j < parts.length; j++) {
      var v = parseFloat(parts[j]);
      if (isNaN(v)) { ok = false; break; }
      nums.push(v);
    }
    if (ok && nums.length) {
      rows.push(nums);
      if (nums.length > cols) cols = nums.length;
    }
  }
  if (rows.length < 2 || cols < 2) return null;
  for (var r = 0; r < rows.length; r++) {
    if (rows[r].length !== cols) return null;
  }
  return {rows: rows, cols: cols};
}

function parseKV(text) {
  var items = [];
  text.split(/\r?\n/).forEach(function(line) {
    var t = line.trim();
    if (t === "") { items.push({type: "blank", raw: line}); return; }
    if (t.charAt(0) === "#" || t.charAt(0) === "!") { items.push({type: "comment", raw: line}); return; }
    var parts = t.split(/\s+/);
    if (parts.length >= 2 && !/^[0-9.-]/.test(parts[0])) {
      items.push({type: "kv", kw: parts[0], val: parts.slice(1).join(" ")});
    } else {
      items.push({type: "other", raw: line});
    }
  });
  return items;
}

function viewFile(name) {
  var rel = joinRel(curPath, name);
  el("viewerOverlay").classList.remove("hidden");
  el("viewerTitle").textContent = name;
  el("viewMsg").textContent = "";
  el("viewTabs").innerHTML = "";
  el("viewBody").innerHTML = '<p class="muted">(loading...)</p>';
  viewer = { name: name, rel: rel, text: null, mode: "text" };
  modalOpened("viewer");
  api("/api/file?path=" + encodeURIComponent(rel)).then(function(resp) {
    if (resp.status === 401) { closeViewer(); return; }
    if (resp.ok) {
      return resp.text().then(function(text) {
        if (viewer.rel !== rel) return;
        renderView(text);
      });
    }
    var ct = resp.headers.get("Content-Type") || "";
    if (ct.indexOf("application/json") >= 0) {
      return resp.json().then(function(j) {
        el("viewBody").innerHTML = "";
        el("viewMsg").textContent = " " + (j.error || "preview failed");
      });
    }
    el("viewBody").innerHTML = "";
    el("viewMsg").textContent = " preview failed";
  }).catch(function() {
    el("viewMsg").textContent = " network error";
  });
}

function closeViewer() {
  el("viewerOverlay").classList.add("hidden");
  modalClosed("viewer");
}

var openModals = {};
function modalOpened(id) {
  openModals[id] = true;
  document.body.style.overflow = "hidden";
}
function modalClosed(id) {
  delete openModals[id];
  if (!Object.keys(openModals).length) document.body.style.overflow = "";
}

function renderView(text) {
  viewer.text = text;
  var numeric = parseNumeric(text);
  var kv = parseKV(text);
  var kvCount = kv.filter(function(it) { return it.type === "kv"; }).length;
  var nonBlank = kv.filter(function(it) { return it.type !== "blank"; }).length;
  var isForm = kvCount >= 2 && kvCount * 2 >= nonBlank;
  var tabs = [];
  if (numeric) tabs.push(["data", "Data"]);
  if (isForm) tabs.push(["form", "Form"]);
  tabs.push(["text", "Text"]);
  tabs.push(["edit", "Edit"]);
  renderTabs(tabs);
  setMode(numeric ? "data" : (isForm ? "form" : "text"));
}

function renderTabs(tabs) {
  var bar = el("viewTabs");
  bar.innerHTML = "";
  tabs.forEach(function(t) {
    var b = document.createElement("button");
    b.className = "tab";
    b.textContent = t[1];
    b.setAttribute("data-mode", t[0]);
    b.onclick = function() { setMode(t[0]); };
    bar.appendChild(b);
  });
}

function setMode(mode) {
  viewer.mode = mode;
  var tabs = el("viewTabs").children;
  for (var i = 0; i < tabs.length; i++) {
    tabs[i].classList.toggle("active", tabs[i].getAttribute("data-mode") === mode);
  }
  var body = el("viewBody");
  body.innerHTML = "";
  el("viewMsg").textContent = "";
  if (mode === "data") renderDataMode(body);
  else if (mode === "form") renderFormMode(body);
  else if (mode === "text") renderTextMode(body);
  else if (mode === "edit") renderEditMode(body);
}

function renderDataMode(body) {
  var numeric = parseNumeric(viewer.text);
  if (!numeric) { body.innerHTML = '<p class="muted">not a numeric table</p>'; return; }
  var wrap = document.createElement("div");
  wrap.className = "table-scroll";
  var table = document.createElement("table");
  table.className = "data";
  var thead = document.createElement("thead");
  var htr = document.createElement("tr");
  for (var c = 0; c < numeric.cols; c++) {
    var th = document.createElement("th");
    th.textContent = "c" + c;
    (function(k) { th.onclick = function() { drawCol(k); }; })(c);
    htr.appendChild(th);
  }
  thead.appendChild(htr);
  table.appendChild(thead);
  var tbody = document.createElement("tbody");
  var maxRows = Math.min(numeric.rows.length, 200);
  for (var r = 0; r < maxRows; r++) {
    var tr = document.createElement("tr");
    for (var c2 = 0; c2 < numeric.cols; c2++) {
      var td = document.createElement("td");
      td.textContent = fmtNum(numeric.rows[r][c2]);
      tr.appendChild(td);
    }
    tbody.appendChild(tr);
  }
  table.appendChild(tbody);
  wrap.appendChild(table);
  body.appendChild(wrap);
  if (numeric.rows.length > maxRows) {
    var note = document.createElement("p");
    note.className = "muted";
    note.textContent = "showing first 200 of " + numeric.rows.length + " rows";
    body.appendChild(note);
  }
  var canvas = document.createElement("canvas");
  canvas.className = "chart";
  canvas.id = "colCanvas";
  canvas.height = 220;
  body.appendChild(canvas);
  var hint = document.createElement("p");
  hint.className = "muted";
  hint.textContent = "click a column header to plot it against the first column";
  body.appendChild(hint);
  drawCol(1);
}

function drawCol(k) {
  var numeric = parseNumeric(viewer.text);
  var canvas = el("colCanvas");
  if (!numeric || !canvas) return;
  drawSeries(canvas,
    [{xs: numeric.rows.map(function(r) { return r[0]; }),
      ys: numeric.rows.map(function(r) { return r[k]; }),
      color: "rgb(24, 103, 174)"}],
    {logY: false, xLabel: "c0", yLabel: "c" + k});
}

function renderFormMode(body) {
  var items = parseKV(viewer.text);
  var form = document.createElement("div");
  items.forEach(function(it, idx) {
    if (it.type === "kv") {
      var row = document.createElement("div");
      row.className = "form-row";
      var lab = document.createElement("label");
      lab.textContent = it.kw;
      lab.title = it.kw;
      var inp = document.createElement("input");
      inp.type = "text";
      inp.value = it.val;
      inp.setAttribute("data-idx", String(idx));
      row.appendChild(lab);
      row.appendChild(inp);
      form.appendChild(row);
    } else if (it.type === "comment" || it.type === "other") {
      var c = document.createElement("div");
      c.className = "form-comment";
      c.textContent = it.raw;
      form.appendChild(c);
    }
  });
  body.appendChild(form);
  var saveRow = document.createElement("div");
  saveRow.className = "row";
  saveRow.style.marginTop = "12px";
  var btn = document.createElement("button");
  btn.textContent = "Save";
  btn.onclick = function() {
    var inputs = el("viewBody").querySelectorAll("input[data-idx]");
    var byIdx = {};
    for (var i = 0; i < inputs.length; i++) byIdx[inputs[i].getAttribute("data-idx")] = inputs[i].value;
    var lines = items.map(function(it, idx) {
      if (it.type === "kv" && byIdx[String(idx)] !== undefined) return it.kw + " " + byIdx[String(idx)];
      return it.raw;
    });
    var eol = viewer.text && viewer.text.indexOf("\r\n") >= 0 ? "\r\n" : "\n";
    saveText(lines.join(eol));
  };
  saveRow.appendChild(btn);
  body.appendChild(saveRow);
}

function renderTextMode(body) {
  var pre = document.createElement("pre");
  pre.className = "preview";
  pre.textContent = viewer.text || "";
  body.appendChild(pre);
}

function renderEditMode(body) {
  var ta = document.createElement("textarea");
  ta.className = "editor";
  ta.id = "editArea";
  ta.value = viewer.text || "";
  body.appendChild(ta);
  var row = document.createElement("div");
  row.className = "row";
  row.style.marginTop = "10px";
  var save = document.createElement("button");
  save.textContent = "Save";
  save.onclick = function() { saveText(el("editArea").value); };
  var revert = document.createElement("button");
  revert.className = "ghost";
  revert.textContent = "Revert";
  revert.onclick = function() { el("editArea").value = viewer.text || ""; };
  row.appendChild(save);
  row.appendChild(revert);
  body.appendChild(row);
}

async function saveText(text) {
  var bytes = new Blob([text]).size;
  if (bytes > 200 * 1024) {
    el("viewMsg").textContent = " content exceeds the 200 KB edit limit";
    return;
  }
  try {
    var resp = await api("/api/save", {
      method: "POST",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify({path: viewer.rel, content: text})
    });
    var data = {};
    try { data = await resp.json(); } catch (e) {}
    if (resp.ok) {
      viewer.text = text;
      toast("Saved " + viewer.name);
      loadScan(curPath, true);
    } else {
      el("viewMsg").textContent = " " + (data.error || "save failed");
    }
  } catch (e) {
    el("viewMsg").textContent = " network error: " + e;
  }
}

function openLightbox(name) {
  var rel = joinRel(curPath, name);
  el("lightboxTitle").textContent = name;
  el("lightboxImg").src = "/api/image?path=" + encodeURIComponent(rel);
  el("lightboxOverlay").classList.remove("hidden");
  modalOpened("lightbox");
}
function closeLightbox() {
  el("lightboxOverlay").classList.add("hidden");
  el("lightboxImg").src = "";
  modalClosed("lightbox");
}

function renderFigures(files) {
  var card = el("figCard");
  var grid = el("figGrid");
  grid.innerHTML = "";
  var imgs = files.filter(function(f) { return /\.(png|jpe?g)$/i.test(f.name); });
  el("figCount").textContent = imgs.length ? imgs.length + " figures" : "";
  if (!imgs.length) { card.classList.add("hidden"); return; }
  card.classList.remove("hidden");
  imgs.forEach(function(f) {
    var cell = document.createElement("div");
    cell.className = "fig-thumb";
    var img = document.createElement("img");
    img.src = "/api/image?path=" + encodeURIComponent(joinRel(curPath, f.name));
    img.alt = f.name;
    img.loading = "lazy";
    var cap = document.createElement("div");
    cap.className = "fig-cap";
    cap.textContent = f.name;
    cell.appendChild(img);
    cell.appendChild(cap);
    cell.onclick = function() { openLightbox(f.name); };
    grid.appendChild(cell);
  });
}

function termAppend(text) {
  var out = el("termOut");
  out.textContent += text;
  out.scrollTop = out.scrollHeight;
}

async function termExec(raw) {
  var cmd = raw.trim();
  if (!cmd) return;
  termAppend("$ " + cmd + "\n");
  if (cmd === "clear") { el("termOut").textContent = ""; return; }
  if (cmd === "pwd") { termAppend((rootDir || "") + "/" + curPath + "\n"); return; }
  if (cmd === "cd" || cmd === "cd ~") {
    var ok0 = await loadScan("");
    termAppend(ok0 ? "-> " + (rootDir || "") + "\n" : "cannot go there\n");
    return;
  }
  if (cmd.indexOf("cd ") === 0 && cmd.indexOf("&&") < 0 && cmd.indexOf(";") < 0) {
    var arg = cmd.substring(3).trim().replace(/^["']|["']$/g, "");
    var dest;
    if (arg === "..") {
      dest = curPath ? curPath.split("/").slice(0, -1).join("/") : "";
    } else if (arg === "/") {
      dest = "";
    } else if (arg.charAt(0) === "/") {
      termAppend(" Error: stay inside the server root; use relative paths or plain 'cd'.\n");
      return;
    } else {
      dest = joinRel(curPath, arg);
    }
    var ok = await loadScan(dest);
    termAppend(ok ? "-> " + dest + "\n" : "no such directory: " + arg + "\n");
    return;
  }
  if (termBusy) { termAppend("(busy: the previous command is still running)\n"); return; }
  termBusy = true;
  var input = el("termInput");
  input.disabled = true;
  try {
    var resp = await api("/api/exec", {
      method: "POST",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify({path: curPath, command: cmd})
    });
    var data = {};
    try { data = await resp.json(); } catch (e) {}
    if (resp.ok) {
      if (data.output) termAppend(data.output + "\n");
      if (data.returncode !== 0) termAppend("[exit code " + data.returncode + "]\n");
      loadScan(curPath, true);
    } else {
      termAppend(" Error: " + (data.error || "command failed") + "\n");
    }
  } catch (e) {
    termAppend(" Network error: " + e + "\n");
  } finally {
    termBusy = false;
    input.disabled = false;
    input.focus();
  }
}

document.addEventListener("DOMContentLoaded", function() {
  var logoImg = el("logoImg");
  if (logoImg) {
    logoImg.onerror = function() {
      logoImg.style.display = "none";
      el("brandText").classList.remove("hidden");
    };
  }
  el("refreshBtn").onclick = function() { loadScan(curPath); };
  el("loginBtn").onclick = doLogin;
  el("pwInput").addEventListener("keydown", function(ev) {
    if (ev.key === "Enter") doLogin();
  });
  el("fileFilter").addEventListener("input", function() {
    renderFileList(lastFiles || [], {});
  });
  el("viewerClose").onclick = closeViewer;
  el("viewerOverlay").onclick = function(ev) {
    if (ev.target === el("viewerOverlay")) closeViewer();
  };
  el("lightboxClose").onclick = closeLightbox;
  el("lightboxOverlay").onclick = function(ev) {
    if (ev.target === el("lightboxOverlay")) closeLightbox();
  };
  document.addEventListener("keydown", function(ev) {
    if (ev.key !== "Escape") return;
    var t = ev.target;
    if (t && t.closest && t.closest("input, textarea")) return;
    if (!el("lightboxOverlay").classList.contains("hidden")) closeLightbox();
    else if (!el("viewerOverlay").classList.contains("hidden")) closeViewer();
  });
  el("termInput").addEventListener("keydown", function(ev) {
    if (ev.key === "Enter") {
      var value = el("termInput").value;
      el("termInput").value = "";
      termExec(value);
    }
  });
  el("termBody").addEventListener("click", function() {
    var sel = window.getSelection ? String(window.getSelection()) : "";
    if (sel) return;
    if (!termBusy) el("termInput").focus();
  });
  el("themeBtn").onclick = toggleTheme;
  el("refreshSel").onchange = startPolling;
  el("railFiles").onclick = function() {
    setRailActive(el("railFiles"));
    window.scrollTo({top: 0, behavior: "smooth"});
  };
  el("railMon").onclick = function() {
    setRailActive(el("railMon"));
    var mon = el("monCard");
    if (mon && !mon.classList.contains("hidden")) mon.scrollIntoView({behavior: "smooth", block: "start"});
    else railScroll("recCard", "Nothing to monitor in this directory");
  };
  el("railTerm").onclick = function() {
    setRailActive(el("railTerm"));
    railScroll("termCard", "Terminal unavailable");
    setTimeout(function() {
      if (!termBusy) el("termInput").focus();
    }, 400);
  };
  el("railHist").onclick = function() {
    setRailActive(el("railHist"));
    railScroll("histCard", "No command history yet");
  };
  var mq = window.matchMedia ? window.matchMedia("(prefers-color-scheme: dark)") : null;
  if (mq && mq.addEventListener) {
    mq.addEventListener("change", function() {
      if (!themeManual) applyTheme();
    });
  }
  applyTheme();
  loadScan("");
  startPolling();
});
