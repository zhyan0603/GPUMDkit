function el(id) { return document.getElementById(id); }
function esc(s) {
  return String(s)
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}
function joinRel(base, name) { return base ? base + "/" + name : name; }
function parentRel(path) { return path ? path.split("/").slice(0, -1).join("/") : ""; }
function fmtSize(n) {
  if (n == null || n < 0) return "-";
  if (n < 1024) return n + " B";
  if (n < 1048576) return (n / 1024).toFixed(1) + " KB";
  return (n / 1048576).toFixed(1) + " MB";
}
function fmtAgo(sec) {
  if (sec == null) return "";
  if (sec < 5) return "now";
  if (sec < 60) return Math.floor(sec) + "s ago";
  if (sec < 3600) return Math.floor(sec / 60) + "m ago";
  if (sec < 86400) return Math.floor(sec / 3600) + "h ago";
  return Math.floor(sec / 86400) + "d ago";
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
function arrMinMax(a) {
  var mn = a[0], mx = a[0];
  for (var i = 1; i < a.length; i++) {
    if (a[i] < mn) mn = a[i];
    if (a[i] > mx) mx = a[i];
  }
  return [mn, mx];
}
function median(a) {
  if (!a.length) return null;
  var s = a.slice().sort(function(x, y) { return x - y; });
  var mid = Math.floor(s.length / 2);
  return s.length % 2 ? s[mid] : (s[mid - 1] + s[mid]) / 2;
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
function toast(msg) {
  var t = el("toast");
  t.textContent = msg;
  t.classList.add("show");
  if (toastTimer) clearTimeout(toastTimer);
  toastTimer = setTimeout(function() { t.classList.remove("show"); }, 2400);
}

var ICON_SUN = '<svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round"><circle cx="12" cy="12" r="5"/><line x1="12" y1="1" x2="12" y2="3"/><line x1="12" y1="21" x2="12" y2="23"/><line x1="4.22" y1="4.22" x2="5.64" y2="5.64"/><line x1="18.36" y1="18.36" x2="19.78" y2="19.78"/><line x1="1" y1="12" x2="3" y2="12"/><line x1="21" y1="12" x2="23" y2="12"/><line x1="4.22" y1="19.78" x2="5.64" y2="18.36"/><line x1="18.36" y1="5.64" x2="19.78" y2="4.22"/></svg>';
var ICON_MOON = '<svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"/></svg>';
var ICON_CHEV = '<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polyline points="6 9 12 15 18 9"/></svg>';

var PREVIEW_MAX_CLIENT = 512 * 1024;
var LOSS_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"];
var KNOWN_COLS = {
  "msd.out": {4: ["t", "msd_x", "msd_y", "msd_z"], 7: ["t", "msd_x", "msd_y", "msd_z", "sdc_x", "sdc_y", "sdc_z"]},
  "sdc.out": {4: ["t", "vac_x", "vac_y", "vac_z"]},
  "loss.out": {6: ["gen", "Loss", "E_tr", "F_tr", "V_tr"], 10: ["gen", "L_total", "L1", "L2", "E_tr", "F_tr", "V_tr", "E_te", "F_te", "V_te"]},
};
var REC_CATEGORIES = {
  plt_msd: "Diffusion and transport",
  plt_sdc: "Diffusion and transport",
  plt_msd_sdc: "Diffusion and transport",
  plt_vac: "Diffusion and transport",
  plt_msd_conv: "Diffusion and transport",
  plt_thermo: "Thermodynamics",
  plt_train: "NEP training",
  plt_train_density: "NEP training",
  plt_train_test: "NEP training",
  plt_prediction: "NEP training",
  plt_sigma: "Arrhenius analysis",
  plt_D: "Arrhenius analysis",
};
var REC_CAT_ORDER = [
  "Diffusion and transport",
  "NEP training",
  "Thermodynamics",
  "Arrhenius analysis",
];

var curPath = "";
var rootDir = "";
var viewMode = "overview";
var running = false;
var termBusy = false;
var busyFlag = false;
var scanSeq = 0;
var pollActive = false;
var pollTimer = null;
var toastTimer = null;
var ageTimer = null;
var themeManual = null;
var connOk = null;
var lastGoodAt = null;
var lastScanServerNow = null;
var lastFiles = null;
var lastSubdirs = null;
var lastPath = null;
var nepSamples = [];
var lossData = null;
var cmdHistory = [];
var figRendered = {};
var panelTab = "output";
var viewer = null;
var pendingOpen = null;
var drawerOpen = false;

try { themeManual = localStorage.getItem("gk_theme"); } catch (e) {}

function setConn(ok) {
  if (connOk !== ok) {
    connOk = ok;
    el("connWrap").classList.toggle("bad", !ok);
  }
  if (ok) lastGoodAt = Date.now();
  updateConnText();
}
function updateConnText() {
  if (connOk === null) { el("connText").textContent = "connecting..."; return; }
  if (!connOk) {
    var age = lastGoodAt ? fmtAgo((Date.now() - lastGoodAt) / 1000) : "never";
    el("connText").textContent = "connection problem · last update " + age;
    return;
  }
  el("connText").textContent = "updated " + fmtAgo((Date.now() - lastGoodAt) / 1000);
}

async function api(url, opts) {
  var resp;
  try {
    resp = await fetch(url, opts);
  } catch (e) {
    setConn(false);
    throw e;
  }
  if (resp.status === 401) showLogin("Please sign in to continue.");
  return resp;
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
    setConn(false);
    if (!silent) el("sbErr").textContent = (data && data.error) || "scan failed";
    return false;
  }
  setConn(true);
  hideLogin();
  var pathChanged = data.path !== lastPath;
  curPath = data.path || "";
  rootDir = data.root || "";
  lastPath = curPath;
  lastScanServerNow = data.now;
  busyFlag = !!data.busy;
  el("sbErr").textContent = "";
  var nowS = data.now || Math.floor(Date.now() / 1000);
  var pollGap = lastPollTs ? (Date.now() / 1000 - lastPollTs) : 0;
  var prevF = {}, prevS = {};
  if (!pathChanged) {
    if (lastFiles) lastFiles.forEach(function(f) { prevF[f.name] = f; });
    if (lastSubdirs) lastSubdirs.forEach(function(s) { prevS[s.name] = s; });
  } else {
    nepSamples = [];
    lossData = null;
  }
  var prevNepDone = lastNepDone;
  var training = data.training;
  lastNepDone = training && !training.finished && training.read_ok !== false ? training.done : null;
  if (training && training.done != null && prevNepDone != null && training.done < prevNepDone) {
    nepSamples = [];
  }
  if (training && !training.finished) {
    nepSamples.push({done: training.done, ts: Date.now() / 1000, count: training.count});
    if (nepSamples.length > 8) nepSamples.shift();
  }
  lastPollTs = Date.now() / 1000;
  lastFiles = data.files || [];
  lastSubdirs = data.subdirs || [];
  var liveFiles = (data.files || []).filter(function(f) {
    return prevF[f.name] && f.size > prevF[f.name].size;
  }).map(function(f) {
    var grew = f.size - prevF[f.name].size;
    return {name: f.name, grew: grew, rate: pollGap > 0 ? grew / pollGap : null, age: f.mtime != null ? nowS - f.mtime : null};
  });
  liveFiles = liveFiles.filter(function(f) { return f.age != null && f.age < 180; });
  var fresh = false;
  var filesArr = data.files || [];
  for (var fi = 0; fi < filesArr.length; fi++) {
    if (filesArr[fi].mtime != null && nowS - filesArr[fi].mtime < 180) { fresh = true; break; }
  }
  if (!fresh) {
    var subsArr = data.subdirs || [];
    for (var si = 0; si < subsArr.length; si++) {
      if (subsArr[si].newest != null && nowS - subsArr[si].newest < 180) { fresh = true; break; }
    }
  }
  pollActive = fresh || (data.training && !data.training.finished);
  renderCrumb(curPath);
  renderBrowser(data, nowS, prevS, prevF);
  el("sbRoot").textContent = rootDir;
  el("sbRoot").title = rootDir;
  el("sbCwd").textContent = "/" + curPath;
  el("sbLive").textContent = liveFiles.length ? liveFiles.length + " live" : "";
  el("sbEntries").textContent = ((data.dirs || []).length + (data.files || []).length) + " entries";
  el("panelBusy").textContent = busyFlag ? "server busy: a command is running" : "";
  if (viewMode === "overview") {
    renderMonitor(data, liveFiles, nowS);
    renderRecs(data.recommendations || [], data.files || [], nowS);
  } else if (viewMode === "figures") {
    renderFigures(data.files || []);
  }
  schedule();
  return true;
}
var lastNepDone = null;
var lastPollTs = 0;

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
    sep.textContent = "/";
    sep.className = "crumbsep";
    c.appendChild(sep);
    if (i === parts.length - 1) {
      var here = document.createElement("span");
      here.textContent = parts[i];
      here.className = "crumbhere";
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

function renderBrowser(data, nowS, prevS, prevF) {
  el("curPathLabel").textContent = (rootDir || "") + "/" + curPath;
  el("curPathLabel").title = (rootDir || "") + "/" + curPath;
  var up = el("upBtn");
  up.disabled = !curPath;
  var fl = el("fileList");
  fl.innerHTML = "";
  var filter = el("fileFilter").value.trim().toLowerCase();
  var sortMode = el("sortSel").value;
  var dirs = (data.subdirs || []).slice();
  var files = (data.files || []).slice();
  files.sort(function(a, b) {
    if (sortMode === "mtime") return (b.mtime || 0) - (a.mtime || 0) || a.name.localeCompare(b.name);
    if (sortMode === "size") return (b.size || 0) - (a.size || 0) || a.name.localeCompare(b.name);
    return a.name.localeCompare(b.name);
  });
  var shown = 0;
  dirs.forEach(function(s) {
    if (filter && s.name.toLowerCase().indexOf(filter) < 0) return;
    shown++;
    var row = document.createElement("div");
    row.className = "brow dirrow";
    var dot = document.createElement("span");
    dot.className = "dot";
    var prev = prevS[s.name];
    if (prev && s.newest != null && prev.newest != null && s.newest > prev.newest) {
      dot.classList.add("live");
      row.title = "recently written";
    }
    row.appendChild(dot);
    var nm = document.createElement("span");
    nm.className = "bname";
    nm.textContent = s.name + "/";
    row.appendChild(nm);
    var meta = document.createElement("span");
    meta.className = "bmeta";
    var bits = [s.count + " items"];
    if (s.newest != null) bits.push(fmtAgo(nowS - s.newest));
    if (s.loss) bits.push("NEP");
    meta.textContent = bits.join(" · ");
    row.appendChild(meta);
    row.onclick = function() { loadScan(joinRel(curPath, s.name)); };
    fl.appendChild(row);
  });
  files.forEach(function(f) {
    if (filter && f.name.toLowerCase().indexOf(filter) < 0) return;
    shown++;
    var row = document.createElement("div");
    row.className = "brow";
    var dot = document.createElement("span");
    dot.className = "dot";
    var prev = prevF[f.name];
    if (prev && f.size > prev.size && f.mtime != null && nowS - f.mtime < 180) {
      dot.classList.add("live");
      row.title = "recently written";
    }
    row.appendChild(dot);
    var a = document.createElement("a");
    a.href = "#";
    a.textContent = f.name;
    a.onclick = function(ev) { ev.preventDefault(); fileClicked(f); };
    row.appendChild(a);
    var meta = document.createElement("span");
    meta.className = "bmeta";
    var bits = [];
    if (f.mtime != null) bits.push(fmtAgo(nowS - f.mtime));
    bits.push(fmtSize(f.size));
    meta.textContent = bits.join(" · ");
    row.appendChild(meta);
    fl.appendChild(row);
  });
  el("browserCount").textContent = shown ? shown + " shown" : "";
  if (!shown) {
    var empty = document.createElement("div");
    empty.className = "empty";
    empty.textContent = filter ? "no files match the filter" : "empty directory";
    fl.appendChild(empty);
  }
}

function fileClicked(f) {
  var ext = f.name.substring(f.name.lastIndexOf(".") + 1).toLowerCase();
  if (ext === "png" || ext === "jpg" || ext === "jpeg") {
    setView("figures");
    openLightbox(joinRel(curPath, f.name), f.name);
  } else {
    setView("data");
    viewFile(f);
  }
}

function setView(mode) {
  viewMode = mode;
  var tabs = document.querySelectorAll(".viewtab");
  for (var i = 0; i < tabs.length; i++) {
    var active = tabs[i].getAttribute("data-view") === mode;
    tabs[i].classList.toggle("active", active);
    tabs[i].setAttribute("aria-pressed", active ? "true" : "false");
  }
  el("viewOverview").classList.toggle("hidden", mode !== "overview");
  el("viewData").classList.toggle("hidden", mode !== "data");
  el("viewFigures").classList.toggle("hidden", mode !== "figures");
  if (drawerOpen) toggleDrawer(false);
}
function toggleDrawer(open) {
  drawerOpen = open == null ? !drawerOpen : open;
  el("browserPanel").classList.toggle("open", drawerOpen);
}

function trainingRate() {
  var deltas = [];
  for (var i = 1; i < nepSamples.length; i++) {
    var d = nepSamples[i].done - nepSamples[i - 1].done;
    var t = nepSamples[i].ts - nepSamples[i - 1].ts;
    if (d > 0 && t > 1) deltas.push(d / t);
  }
  if (deltas.length < 2) return {rate: null, estimating: true};
  return {rate: median(deltas), estimating: false};
}

function statTile(value, label, done) {
  return '<div class="stat' + (done ? " done" : "") + '"><div class="stat-v">' + value + '</div><div class="stat-l">' + label + '</div></div>';
}

function renderMonitor(data, liveFiles, nowS) {
  var body = el("monBody");
  body.innerHTML = "";
  var monMeta = el("monMeta");
  var t = data.training;
  if (!t) {
    monMeta.textContent = "";
    body.innerHTML = '<div class="empty">No training records in this directory (nep.in and loss.out not both found).</div>';
    return;
  }
  var parts = [];
  var evidence = [];
  var notes = [];
  if (t.loss_mtime == null) {
    monMeta.textContent = "read failure";
    body.innerHTML = '<div class="empty">loss.out could not be read.</div>';
    return;
  }
  if (t.loss_empty) {
    monMeta.textContent = "no records";
    body.innerHTML = '<div class="empty">loss.out is empty or contains no valid records yet.</div>';
    return;
  }
  var writeAge = nowS - t.loss_mtime;
  var stateText;
  if (t.finished) {
    stateText = "finished · reached the target generation count";
    monMeta.textContent = "finished";
  } else if (writeAge < 180) {
    stateText = "updating · last write " + fmtAgo(writeAge);
    monMeta.textContent = "updating";
  } else {
    stateText = "no update for " + fmtAgo(writeAge) + " · the job may have ended, or it writes sparsely";
    monMeta.textContent = "no recent update";
  }
  var rateInfo = trainingRate();
  var eta = null;
  if (!t.finished && t.has_target && rateInfo.rate) {
    eta = (t.total - t.done) / rateInfo.rate;
  }
  var useBar = !t.multi_run && t.has_target;
  var tiles = '<div class="stat-row">';
  if (useBar) {
    var pct = t.total > 0 ? 100 * t.done / t.total : 0;
    tiles += statTile(pct.toFixed(1) + "%", "progress", t.finished);
  } else {
    tiles += statTile(t.count.toLocaleString(), "records", t.finished);
  }
  tiles += statTile(t.done.toLocaleString() + " / " + t.total.toLocaleString(), "generation");
  tiles += statTile(t.finished ? "--" : fmtGenRate(rateInfo.rate), "rate");
  if (!t.has_target) {
    tiles += statTile("--", "eta");
    notes.push("nep.in has no generation line; the target defaults to " + t.total.toLocaleString() + ", so no exact ETA is shown");
  } else if (t.finished) {
    tiles += statTile("--", "eta");
  } else {
    tiles += statTile(rateInfo.estimating ? "estimating" : (eta != null ? fmtDur(eta) : "--"), "eta");
    if (rateInfo.estimating) notes.push("ETA needs a few more refresh samples");
  }
  tiles += "</div>";
  evidence.push("loss.out · " + t.count.toLocaleString() + " records · last write " + fmtAgo(writeAge));
  evidence.push("target: " + (t.has_target ? "nep.in generation " + t.total.toLocaleString() : "default (no generation line in nep.in)"));
  if (t.multi_run) {
    notes.push("loss.out contains records from more than one run (generations restart at " + t.seg_start_gen + "); the chart shows the whole history and progress uses the latest run");
  }
  if (t.skipped > 0) {
    notes.push(t.skipped + " invalid or incomplete lines were skipped while parsing");
  }
  parts.push(
    '<div class="stat-row">' + tiles.replace('<div class="stat-row">', "").replace("</div>", "") + "</div>"
  );
  var barHtml = "";
  if (useBar) {
    var pct2 = Math.min(100, 100 * t.done / t.total);
    barHtml = '<div class="bar-track"><div class="bar-fill' + (t.finished ? " done" : "") + '" style="width:' + pct2 + '%"></div></div>';
  }
  parts.push(barHtml);
  parts.push('<div class="evidence">' + evidence.join(" · ") + "</div>");
  parts.push('<div class="muted" style="margin-top:4px">' + stateText + "</div>");
  parts.push('<div class="chart-wrap"><canvas class="chart" id="lossCanvas" height="230"></canvas><span class="pulse-dot hidden" id="lossDot"></span></div>');
  parts.push('<div class="lg-row" id="lossLegend"></div>');
  parts.push('<div id="lossNotes"></div>');
  if (liveFiles.length) {
    var rows = liveFiles.map(function(f) {
      return '<div class="row" style="justify-content:space-between;padding:2px 0">' +
        '<span class="mono" style="font-size:12px">' + esc(f.name) + "</span>" +
        '<span class="muted mono">+' + fmtSize(f.grew) + (f.rate ? " · " + fmtRate(f.rate) : "") + "</span></div>";
    }).join("");
    parts.push('<div style="font-weight:600;font-size:12.5px;margin:10px 0 3px">Growing files</div>' + rows);
  }
  body.innerHTML = parts.join("");
  var notesEl = el("lossNotes");
  notes.forEach(function(n) {
    var d = document.createElement("div");
    d.className = "chart-note";
    d.textContent = n;
    notesEl.appendChild(d);
  });
  fetchLoss().then(function() {
    if (curPath !== lossFetchPath) return;
    drawLossChart();
  }).catch(function() {});
}
var lossFetchPath = null;
async function fetchLoss() {
  lossFetchPath = curPath;
  var resp = await api("/api/loss?path=" + encodeURIComponent(curPath));
  if (resp.status === 401) return;
  if (curPath !== lossFetchPath) return;
  if (resp.ok) {
    lossData = await resp.json();
    if (curPath !== lossFetchPath) { lossData = null; return; }
  } else {
    lossData = null;
  }
}

function drawLossChart() {
  var canvas = el("lossCanvas");
  if (!canvas || !lossData || lossData.empty) {
    var dot = el("lossDot");
    if (dot) dot.classList.add("hidden");
    return;
  }
  var xs = lossData.xs || [];
  var series = (lossData.series || []).map(function(s, i) {
    return {xs: xs, ys: s.values, color: LOSS_COLORS[i % LOSS_COLORS.length], label: s.label};
  });
  var xLabel = lossData.x_mode === "generation" ? "generation" : "record #";
  var last = drawSeries(canvas, series, {logX: true, logY: true, xLabel: xLabel, yLabel: "Loss functions"});
  var legend = el("lossLegend");
  if (legend) {
    legend.innerHTML = "";
    series.forEach(function(s) {
      var item = document.createElement("span");
      item.className = "lg-item";
      var d = document.createElement("span");
      d.className = "lg-dot";
      d.style.background = s.color;
      item.appendChild(d);
      item.appendChild(document.createTextNode(s.label));
      legend.appendChild(item);
    });
  }
  var dot = el("lossDot");
  if (dot) {
    if (last && last.length) {
      var p = last[0];
      dot.style.left = p.px + "px";
      dot.style.top = p.py + "px";
      dot.style.background = p.color;
      dot.style.setProperty("--pc", hexToRgba(p.color, 0.55));
      dot.classList.remove("hidden");
    } else {
      dot.classList.add("hidden");
    }
  }
  var notesEl = el("lossNotes");
  if (notesEl && lossData.sampled) {
    var d = document.createElement("div");
    d.className = "chart-note";
    d.textContent = "chart downsampled to " + lossData.points + " of " + lossData.count + " records";
    notesEl.appendChild(d);
  }
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
function drawSeries(canvas, series, opts) {
  opts = opts || {};
  var theme = themeColors();
  var dpr = window.devicePixelRatio || 1;
  var w = canvas.clientWidth || 600;
  var h = canvas.clientHeight || 230;
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
  var pad = {l: 60, r: 14, t: 12, b: 40};
  var pw = w - pad.l - pad.r, ph = h - pad.t - pad.b;
  function X(x) { return pad.l + pw * (x - xmin) / (xmax - xmin); }
  function Y(y) { return pad.t + ph * (1 - (y - ymin) / (ymax - ymin)); }
  ctx.strokeStyle = theme.grid;
  ctx.lineWidth = 1;
  ctx.font = "10px monospace";
  ctx.fillStyle = theme.axis;
  axisTicks(ymin, ymax, opts.logY, 4).forEach(function(t) {
    var yy = Y(t.pos);
    ctx.beginPath();
    ctx.moveTo(pad.l, yy);
    ctx.lineTo(w - pad.r, yy);
    ctx.stroke();
    ctx.fillText(t.label, 6, yy + 3);
  });
  axisTicks(xmin, xmax, opts.logX, 4).forEach(function(t) {
    var xx = X(t.pos);
    ctx.beginPath();
    ctx.moveTo(xx, pad.t);
    ctx.lineTo(xx, h - pad.b);
    ctx.stroke();
    ctx.fillText(t.label, xx - 14, h - pad.b + 16);
  });
  if (opts.xLabel) ctx.fillText(opts.xLabel, pad.l + pw / 2 - ctx.measureText(opts.xLabel).width / 2, h - 8);
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
  canvas._plot = {X: X, xmin: xmin, xmax: xmax, series: series, opts: opts};
  return lastPoints;
}

function attachHover(canvas, readoutEl, labelsFn) {
  canvas.onmousemove = function(ev) {
    var rect = canvas.getBoundingClientRect();
    var px = ev.clientX - rect.left;
    var plot = canvas._plot;
    if (!plot || !plot.series.length || !readoutEl) return;
    var xs = plot.series[0].xs;
    if (!xs || !xs.length) return;
    var frac = (px - 60) / (canvas.clientWidth - 74);
    if (frac < 0) frac = 0;
    if (frac > 1) frac = 1;
    var idx = Math.round(frac * (xs.length - 1));
    readoutEl.textContent = labelsFn(idx, plot.series);
  };
  canvas.onmouseleave = function() {
    if (readoutEl) readoutEl.textContent = "";
  };
}

function renderRecs(recs, files, nowS) {
  var box = el("recBox");
  box.innerHTML = "";
  el("recCount").textContent = recs.length ? recs.length + " actions" : "";
  if (!recs.length) {
    box.innerHTML = '<div class="empty">No recognized GPUMDkit output files. Run a simulation, or open a directory that contains results.</div>';
    return;
  }
  var fileMap = {};
  files.forEach(function(f) { fileMap[f.name] = f; });
  var groups = {};
  recs.forEach(function(rec) {
    var cat = REC_CATEGORIES[rec.action] || "Other";
    if (!groups[cat]) groups[cat] = [];
    groups[cat].push(rec);
  });
  REC_CAT_ORDER.concat(Object.keys(groups).filter(function(c) { return REC_CAT_ORDER.indexOf(c) < 0; })).forEach(function(cat) {
    var list = groups[cat];
    if (!list) return;
    var grp = document.createElement("div");
    grp.className = "recgroup";
    var head = document.createElement("div");
    head.className = "recgroup-head";
    head.innerHTML = '<span class="chev">' + ICON_CHEV + "</span>" + esc(cat) + ' <span class="muted">' + list.length + "</span>";
    head.onclick = function() { grp.classList.toggle("closed"); };
    grp.appendChild(head);
    var body = document.createElement("div");
    body.className = "recgroup-body";
    list.forEach(function(rec) {
      body.appendChild(buildRecCard(rec, fileMap, nowS));
    });
    grp.appendChild(body);
    box.appendChild(grp);
  });
}

function buildRecCard(rec, fileMap, nowS) {
  var div = document.createElement("div");
  div.className = "rec";
  var head = document.createElement("div");
  head.className = "rec-title";
  head.textContent = rec.description;
  if (rec.stale) {
    var badge = document.createElement("span");
    badge.className = "badge warn";
    badge.textContent = "result older than inputs";
    badge.style.marginLeft = "8px";
    head.appendChild(badge);
  }
  var io = document.createElement("div");
  io.className = "rec-io";
  io.textContent = "input: " + rec.evidence.join(", ") + "  ->  output: " + rec.produces;
  var result = document.createElement("div");
  result.className = "rec-result";
  var outFile = fileMap[rec.produces];
  if (outFile && outFile.mtime != null) {
    result.textContent = "result: " + rec.produces + " · " + fmtAgo(nowS - outFile.mtime);
  } else {
    result.textContent = "result: not generated yet";
    result.className += " muted";
  }
  var actions = document.createElement("div");
  actions.className = "rec-actions";
  var run = document.createElement("button");
  run.textContent = "Generate figure";
  run.disabled = busyFlag || running;
  run.title = busyFlag ? "server busy: another command is running" : "";
  run.onclick = function() { runAction(rec, run); };
  actions.appendChild(run);
  if (outFile) {
    var open = document.createElement("button");
    open.className = "ghost small";
    open.textContent = "Open result";
    open.onclick = function() {
      setView("figures");
      openLightbox(joinRel(curPath, rec.produces), rec.produces);
    };
    actions.appendChild(open);
  }
  var cmdToggle = document.createElement("button");
  cmdToggle.className = "ghost small";
  cmdToggle.textContent = "View command";
  var cmdrow = document.createElement("div");
  cmdrow.className = "rec-cmdrow hidden";
  cmdrow.innerHTML = '<code class="rec-cmd">' + esc(rec.command) + "</code>";
  var copyBtn = document.createElement("button");
  copyBtn.className = "ghost small";
  copyBtn.textContent = "Copy";
  copyBtn.onclick = function() { copyText(rec.command); };
  cmdrow.appendChild(copyBtn);
  cmdToggle.onclick = function() {
    var hidden = cmdrow.classList.contains("hidden");
    cmdrow.classList.toggle("hidden", !hidden);
    cmdToggle.textContent = hidden ? "Hide command" : "View command";
  };
  actions.appendChild(cmdToggle);
  div.appendChild(head);
  div.appendChild(io);
  div.appendChild(result);
  div.appendChild(actions);
  div.appendChild(cmdrow);
  return div;
}

function viewerDirty() {
  return !!(viewer && viewer.draft != null && viewer.draft !== viewer.text);
}

function viewFile(f) {
  if (viewerDirty()) {
    pendingOpen = f;
    el("viewerConfirm").classList.remove("hidden");
    return;
  }
  openViewer(f);
}

function openViewer(f) {
  var rel = joinRel(curPath, f.name);
  var card = el("viewerBody");
  card.innerHTML = '<div class="empty">(loading...)</div>';
  viewer = {
    rel: rel, name: f.name, size: f.size, baseMtime: null,
    truncated: false, text: null, draft: null, mode: "text"
  };
  el("viewerConfirm").classList.add("hidden");
  el("viewerConflict").classList.add("hidden");
  var oversized = f.size != null && f.size > PREVIEW_MAX_CLIENT;
  var url = "/api/file?path=" + encodeURIComponent(rel) + (oversized ? "&partial=1" : "");
  api(url).then(function(resp) {
    if (resp.status === 401) { openViewerReset(); return; }
    if (viewer === null || viewer.rel !== rel) return;
    if (resp.ok) {
      var ct = resp.headers.get("Content-Type") || "";
      if (ct.indexOf("application/json") >= 0) {
        return resp.json().then(function(j) {
          if (viewer === null || viewer.rel !== rel) return;
          viewer.truncated = true;
          viewer.baseMtime = j.mtime;
          viewer.text = j.head + "\n[... " + (j.size - j.head_bytes - j.tail_bytes).toLocaleString() + " bytes omitted ...]\n" + j.tail;
          viewer.size = j.size;
          renderViewer();
        });
      }
      return resp.text().then(function(text) {
        if (viewer === null || viewer.rel !== rel) return;
        viewer.baseMtime = parseInt(resp.headers.get("X-File-Mtime") || "0", 10) || null;
        viewer.text = text;
        renderViewer();
      });
    }
    if (resp.status === 413) {
      return fetchPartial(rel).then(function(j) {
        if (viewer === null || viewer.rel !== rel) return;
        if (!j) {
          openViewerReset();
          el("viewerMsg").textContent = " file too large to preview";
          return;
        }
        viewer.truncated = true;
        viewer.baseMtime = j.mtime;
        viewer.text = j.head + "\n[... " + (j.size - j.head_bytes - j.tail_bytes).toLocaleString() + " bytes omitted ...]\n" + j.tail;
        viewer.size = j.size;
        renderViewer();
      });
    }
    var ct2 = resp.headers.get("Content-Type") || "";
    if (ct2.indexOf("application/json") >= 0) {
      return resp.json().then(function(j) {
        if (viewer === null || viewer.rel !== rel) return;
        openViewerReset();
        el("viewerMsg").textContent = " " + (j.error || "preview failed");
      });
    }
    openViewerReset();
    el("viewerMsg").textContent = " preview failed";
  }).catch(function() {
    if (viewer !== null && viewer.rel === rel) {
      openViewerReset();
      el("viewerMsg").textContent = " network error";
    }
  });
}

function fetchPartial(rel) {
  return api("/api/file?path=" + encodeURIComponent(rel) + "&partial=1").then(function(resp) {
    if (resp.ok) return resp.json();
    return null;
  }).catch(function() { return null; });
}

function openViewerReset() {
  el("viewerBody").innerHTML = '<div class="empty">Select a file from the browser to preview, plot, or edit.</div><p class="err" id="viewerMsg"></p>';
}

function renderViewer() {
  if (!viewer) return;
  var card = el("viewerBody");
  card.innerHTML = "";
  el("viewerMsg").textContent = "";
  var head = document.createElement("div");
  head.className = "viewer-head";
  var title = document.createElement("strong");
  title.textContent = viewer.name + " · " + fmtSize(viewer.size);
  head.appendChild(title);
  if (viewerDirty()) {
    var dirty = document.createElement("span");
    dirty.className = "badge warn";
    dirty.textContent = "unsaved changes";
    head.appendChild(dirty);
  }
  var tabs = document.createElement("div");
  tabs.className = "tabs";
  tabs.id = "viewTabs";
  head.appendChild(tabs);
  card.appendChild(head);
  var banner = document.createElement("div");
  if (viewer.truncated) {
    banner.className = "chart-note";
    banner.textContent = "bounded preview: first 128 KB and last 64 KB shown; editing is disabled because the full content is not loaded";
    card.appendChild(banner);
  }
  var msg = document.createElement("p");
  msg.className = "err";
  msg.id = "viewerMsg";
  card.appendChild(msg);
  var confirmBar = document.createElement("div");
  confirmBar.className = "confbar warn hidden";
  confirmBar.id = "viewerConfirm";
  confirmBar.innerHTML = "This file has unsaved changes.";
  var keepBtn = document.createElement("button");
  keepBtn.textContent = "Keep editing";
  var discardBtn = document.createElement("button");
  discardBtn.className = "danger";
  discardBtn.textContent = "Discard and open the new file";
  confirmBar.appendChild(keepBtn);
  confirmBar.appendChild(discardBtn);
  card.appendChild(confirmBar);
  keepBtn.onclick = function() {
    confirmBar.classList.add("hidden");
    pendingOpen = null;
  };
  discardBtn.onclick = function() {
    confirmBar.classList.add("hidden");
    var target = pendingOpen;
    pendingOpen = null;
    viewer.draft = null;
    if (target) openViewer(target);
  };
  var conflictBar = document.createElement("div");
  conflictBar.className = "confbar conflict hidden";
  conflictBar.id = "viewerConflict";
  conflictBar.textContent = "The file changed on disk while you were editing.";
  var rereadBtn = document.createElement("button");
  rereadBtn.textContent = "Reload from disk";
  var keepDraftBtn = document.createElement("button");
  keepDraftBtn.textContent = "Keep my draft";
  var overBtn = document.createElement("button");
  overBtn.className = "danger";
  overBtn.textContent = "Overwrite disk anyway";
  conflictBar.appendChild(rereadBtn);
  conflictBar.appendChild(keepDraftBtn);
  conflictBar.appendChild(overBtn);
  card.appendChild(conflictBar);
  rereadBtn.onclick = function() {
    conflictBar.classList.add("hidden");
    var f = {name: viewer.name, size: null};
    viewer.draft = null;
    openViewer(f);
  };
  keepDraftBtn.onclick = function() { conflictBar.classList.add("hidden"); };
  overBtn.onclick = function() { conflictBar.classList.add("hidden"); saveDraft(true); };
  var body = document.createElement("div");
  body.id = "viewBody";
  card.appendChild(body);
  var text = viewer.draft != null ? viewer.draft : viewer.text;
  var numeric = viewer.truncated ? null : parseNumeric(text || "");
  var kv = viewer.truncated ? [] : parseKV(text || "");
  var kvCount = kv.filter(function(it) { return it.type === "kv"; }).length;
  var nonBlank = kv.filter(function(it) { return it.type !== "blank"; }).length;
  var isForm = !viewer.truncated && kvCount >= 2 && kvCount * 2 >= nonBlank;
  var tabList = [];
  if (numeric) tabList.push(["data", "Data"]);
  if (isForm) tabList.push(["form", "Form"]);
  tabList.push(["text", "Text"]);
  if (!viewer.truncated) tabList.push(["edit", "Edit"]);
  if (viewer.mode !== "text" && viewer.mode !== "edit" && viewer.mode !== "form" && viewer.mode !== "data") viewer.mode = "text";
  if (viewer.mode === "data" && !numeric) viewer.mode = isForm ? "form" : "text";
  if (viewer.mode === "form" && !isForm) viewer.mode = "text";
  if (viewer.mode === "edit" && viewer.truncated) viewer.mode = "text";
  tabList.forEach(function(t) {
    var b = document.createElement("button");
    b.className = "tab";
    b.textContent = t[1];
    b.setAttribute("data-mode", t[0]);
    b.onclick = function() { setViewerMode(t[0]); };
    tabs.appendChild(b);
  });
  setViewerMode(viewer.mode);
}

function setViewerMode(mode) {
  if (!viewer) return;
  viewer.mode = mode;
  var tabs = el("viewTabs");
  if (tabs) {
    for (var i = 0; i < tabs.children.length; i++) {
      tabs.children[i].classList.toggle("active", tabs.children[i].getAttribute("data-mode") === mode);
    }
  }
  var body = el("viewBody");
  if (!body) return;
  body.innerHTML = "";
  if (mode === "data") renderDataMode(body);
  else if (mode === "form") renderFormMode(body);
  else if (mode === "text") renderTextMode(body);
  else if (mode === "edit") renderEditMode(body);
}

function colNames(name, cols) {
  var known = KNOWN_COLS[name];
  if (known && known[cols]) return known[cols];
  var names = [];
  for (var i = 0; i < cols; i++) names.push("c" + i);
  return names;
}

function renderDataMode(body) {
  var text = viewer.draft != null ? viewer.draft : viewer.text;
  var numeric = parseNumeric(text || "");
  if (!numeric) { body.innerHTML = '<div class="empty">not a numeric table</div>'; return; }
  var names = colNames(viewer.name, numeric.cols);
  var state = viewer._plotState || {x: 0, ys: [1], logX: false, logY: false};
  if (state.ys.length > 1 || numeric.cols <= state.ys[0]) state.ys = state.ys.filter(function(c) { return c < numeric.cols; });
  if (!state.ys.length) state.ys = [Math.min(1, numeric.cols - 1)];
  viewer._plotState = state;
  var ctl = document.createElement("div");
  ctl.className = "chartctl";
  var xSel = document.createElement("select");
  xSel.setAttribute("aria-label", "X column");
  for (var c = 0; c < numeric.cols; c++) {
    var o = document.createElement("option");
    o.value = c;
    o.textContent = names[c];
    if (c === state.x) o.selected = true;
    xSel.appendChild(o);
  }
  xSel.onchange = function() { state.x = parseInt(xSel.value, 10); renderDataMode(body); };
  ctl.appendChild(document.createTextNode("X "));
  ctl.appendChild(xSel);
  ctl.appendChild(document.createTextNode(" Y "));
  for (var c2 = 0; c2 < numeric.cols; c2++) {
    (function(ci) {
      var lab = document.createElement("label");
      lab.className = "chk";
      var cb = document.createElement("input");
      cb.type = "checkbox";
      cb.checked = state.ys.indexOf(ci) >= 0;
      cb.onchange = function() {
        if (cb.checked) { if (state.ys.indexOf(ci) < 0) state.ys.push(ci); }
        else state.ys = state.ys.filter(function(v) { return v !== ci; });
        renderDataMode(body);
      };
      lab.appendChild(cb);
      lab.appendChild(document.createTextNode(names[ci]));
      ctl.appendChild(lab);
    })(c2);
  }
  var scaleX = document.createElement("label");
  scaleX.className = "chk";
  var cbX = document.createElement("input");
  cbX.type = "checkbox";
  cbX.checked = state.logX;
  cbX.onchange = function() { state.logX = cbX.checked; renderDataMode(body); };
  scaleX.appendChild(cbX);
  scaleX.appendChild(document.createTextNode("log X"));
  ctl.appendChild(scaleX);
  var scaleY = document.createElement("label");
  scaleY.className = "chk";
  var cbY = document.createElement("input");
  cbY.type = "checkbox";
  cbY.checked = state.logY;
  cbY.onchange = function() { state.logY = cbY.checked; renderDataMode(body); };
  scaleY.appendChild(cbY);
  scaleY.appendChild(document.createTextNode("log Y"));
  ctl.appendChild(scaleY);
  body.appendChild(ctl);
  var canvas = document.createElement("canvas");
  canvas.className = "chart";
  canvas.id = "colCanvas";
  canvas.height = 230;
  body.appendChild(canvas);
  var readout = document.createElement("div");
  readout.className = "hoverread";
  body.appendChild(readout);
  var xs = numeric.rows.map(function(r) { return r[state.x]; });
  var series = state.ys.map(function(ci, i) {
    return {xs: xs, ys: numeric.rows.map(function(r) { return r[ci]; }), color: LOSS_COLORS[i % LOSS_COLORS.length], label: names[ci]};
  });
  drawSeries(canvas, series, {logX: state.logX, logY: state.logY, xLabel: names[state.x], yLabel: "value"});
  attachHover(canvas, readout, function(idx, sers) {
    var parts = [names[state.x] + "=" + fmtNum(numeric.rows[idx][state.x])];
    state.ys.forEach(function(ci) {
      parts.push(names[ci] + "=" + fmtNum(numeric.rows[idx][ci]));
    });
    return "row " + (idx + 1) + " · " + parts.join(" · ");
  });
  var wrap = document.createElement("div");
  wrap.className = "table-scroll";
  var table = document.createElement("table");
  table.className = "data";
  var thead = document.createElement("thead");
  var htr = document.createElement("tr");
  for (var c3 = 0; c3 < numeric.cols; c3++) {
    var th = document.createElement("th");
    th.textContent = names[c3];
    htr.appendChild(th);
  }
  thead.appendChild(htr);
  table.appendChild(thead);
  var tbody = document.createElement("tbody");
  var maxRows = Math.min(numeric.rows.length, 200);
  for (var r = 0; r < maxRows; r++) {
    var tr = document.createElement("tr");
    for (var c4 = 0; c4 < numeric.cols; c4++) {
      var td = document.createElement("td");
      td.textContent = fmtNum(numeric.rows[r][c4]);
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

function renderFormMode(body) {
  var text = viewer.draft != null ? viewer.draft : viewer.text;
  var items = parseKV(text || "");
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
  var inputs = form.querySelectorAll("input[data-idx]");
  for (var i = 0; i < inputs.length; i++) {
    inputs[i].addEventListener("input", function() {
      var byIdx = {};
      for (var k = 0; k < inputs.length; k++) byIdx[inputs[k].getAttribute("data-idx")] = inputs[k].value;
      var lines = items.map(function(it, idx2) {
        if (it.type === "kv" && byIdx[String(idx2)] !== undefined) return it.kw + " " + byIdx[String(idx2)];
        return it.raw;
      });
      setDraft(joinEol(lines));
    });
  }
  var saveRow = document.createElement("div");
  saveRow.className = "row";
  saveRow.style.marginTop = "10px";
  var btn = document.createElement("button");
  btn.textContent = "Save";
  btn.onclick = function() { saveDraft(false); };
  saveRow.appendChild(btn);
  body.appendChild(saveRow);
}

function joinEol(lines) {
  var base = viewer.text || "";
  var eol = base.indexOf("\r\n") >= 0 ? "\r\n" : "\n";
  return lines.join(eol);
}

function renderTextMode(body) {
  var pre = document.createElement("pre");
  pre.className = "preview";
  pre.textContent = viewer.draft != null ? viewer.draft : (viewer.text || "");
  body.appendChild(pre);
}

function renderEditMode(body) {
  var ta = document.createElement("textarea");
  ta.className = "editor";
  ta.id = "editArea";
  ta.value = viewer.draft != null ? viewer.draft : (viewer.text || "");
  ta.addEventListener("input", function() { setDraft(ta.value); });
  body.appendChild(ta);
  var row = document.createElement("div");
  row.className = "row";
  row.style.marginTop = "8px";
  var save = document.createElement("button");
  save.textContent = "Save";
  save.onclick = function() { saveDraft(false); };
  var revert = document.createElement("button");
  revert.className = "ghost";
  revert.textContent = "Revert to saved";
  revert.onclick = function() {
    viewer.draft = null;
    renderViewer();
  };
  row.appendChild(save);
  row.appendChild(revert);
  body.appendChild(row);
}

function setDraft(text) {
  if (!viewer) return;
  viewer.draft = text;
  if (viewer.draft === viewer.text) viewer.draft = null;
  var head = el("viewerBody").querySelector(".viewer-head");
  if (head) {
    var existing = head.querySelector(".badge.warn");
    if (viewerDirty() && !existing) {
      var dirty = document.createElement("span");
      dirty.className = "badge warn";
      dirty.textContent = "unsaved changes";
      head.insertBefore(dirty, head.querySelector(".tabs"));
    } else if (!viewerDirty() && existing) {
      existing.remove();
    }
  }
}

async function saveDraft(force) {
  if (!viewer || viewer.draft == null) return;
  var content = viewer.draft;
  var saveRel = viewer.rel;
  var bytes = new Blob([content]).size;
  if (bytes > 200 * 1024) {
    el("viewerMsg").textContent = " content exceeds the 200 KB edit limit";
    return;
  }
  try {
    var resp = await api("/api/save", {
      method: "POST",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify({path: saveRel, content: content, base_mtime: force ? undefined : viewer.baseMtime})
    });
    var data = {};
    try { data = await resp.json(); } catch (e) {}
    if (resp.status === 409 && data.error === "conflict") {
      viewer.baseMtime = data.mtime;
      el("viewerConflict").classList.remove("hidden");
      return;
    }
    if (resp.ok) {
      if (viewer && viewer.rel === saveRel) {
        viewer.text = content;
        viewer.draft = null;
        viewer.baseMtime = data.mtime || viewer.baseMtime;
        renderViewer();
      }
      toast("Saved " + viewer.name);
      loadScan(curPath, true);
    } else {
      el("viewerMsg").textContent = " " + (data.error || "save failed");
    }
  } catch (e) {
    el("viewerMsg").textContent = " network error: " + e;
  }
}

function renderFigures(files) {
  var box = el("figBox");
  var imgs = files.filter(function(f) { return /\.(png|jpe?g)$/i.test(f.name); });
  el("figCount").textContent = imgs.length ? imgs.length + " figures" : "";
  if (!imgs.length) {
    box.innerHTML = '<div class="empty">No PNG or JPEG images in this directory.</div>';
    figRendered = {};
    return;
  }
  var keep = {};
  imgs.forEach(function(f) { keep[f.name + ":" + f.size] = true; });
  var existing = box.querySelectorAll(".fig-thumb");
  for (var i = existing.length - 1; i >= 0; i--) {
    var node = existing[i];
    var key = node.getAttribute("data-key");
    if (!keep[key]) {
      node.remove();
      delete figRendered[key];
    }
  }
  if (box.querySelector(".empty")) box.querySelector(".empty").remove();
  imgs.forEach(function(f) {
    var key = f.name + ":" + f.size;
    if (figRendered[key]) return;
    var cell = document.createElement("div");
    cell.className = "fig-thumb";
    cell.setAttribute("data-key", key);
    var img = document.createElement("img");
    img.src = "/api/image?path=" + encodeURIComponent(joinRel(curPath, f.name));
    img.alt = f.name;
    img.loading = "lazy";
    var cap = document.createElement("div");
    cap.className = "fig-cap";
    cap.textContent = f.name;
    cell.appendChild(img);
    cell.appendChild(cap);
    cell.onclick = function() { openLightbox(joinRel(curPath, f.name), f.name); };
    box.appendChild(cell);
    figRendered[key] = true;
  });
}

function openLightbox(rel, name) {
  el("lightboxTitle").textContent = name;
  el("lightboxImg").src = "/api/image?path=" + encodeURIComponent(rel);
  el("lightboxOverlay").classList.remove("hidden");
  el("lightboxClose").focus();
}
function closeLightbox() {
  el("lightboxOverlay").classList.add("hidden");
  el("lightboxImg").src = "";
}

function setPanelTab(tab) {
  panelTab = tab;
  var tabs = document.querySelectorAll(".ptab");
  for (var i = 0; i < tabs.length; i++) {
    var active = tabs[i].getAttribute("data-tab") === tab;
    tabs[i].classList.toggle("active", active);
    tabs[i].setAttribute("aria-pressed", active ? "true" : "false");
  }
  el("pOutput").classList.toggle("hidden", tab !== "output");
  el("pHistory").classList.toggle("hidden", tab !== "history");
  el("pTerminal").classList.toggle("hidden", tab !== "terminal");
  expandPanel(true);
}
function expandPanel(open) {
  el("bottomPanel").classList.toggle("closed", open === false);
}

function showOutput(text, truncated) {
  el("output").textContent = text;
  el("outTruncNote").classList.toggle("hidden", !truncated);
  setPanelTab("output");
}

async function runAction(rec, btn) {
  if (running) return;
  var runDir = curPath;
  var runDirLabel = "/" + runDir;
  running = true;
  btn.disabled = true;
  btn.textContent = "Running...";
  showOutput("$ " + rec.command + "\n\ndirectory: " + (runDir ? runDir : "(root)") + "\n(running, please wait...)", false);
  var t0 = Date.now();
  var record = {
    cmd: rec.command, dir: runDir, start: new Date(), dur: null,
    status: "running", rc: null, output: "", truncated: false, image: null
  };
  try {
    var resp = await api("/api/run", {
      method: "POST",
      headers: {"Content-Type": "application/json"},
      body: JSON.stringify({action: rec.action, path: runDir})
    });
    var data = null;
    try { data = await resp.json(); } catch (e) {}
    if (!resp.ok || !data) {
      record.status = data && data.error && data.error.indexOf("timed out") >= 0 ? "timeout" : "failed";
      record.output = (data && data.error) || "The command could not be started.";
      showOutput("$ " + rec.command + "\n\ndirectory: " + (runDir ? runDir : "(root)") + "\n\n" + record.output, false);
    } else {
      record.status = data.returncode === 0 ? "ok" : "exit " + data.returncode;
      record.rc = data.returncode;
      record.output = data.output || "";
      record.truncated = !!data.truncated;
      record.image = data.image || null;
      var text = "$ " + data.command + "\n\ndirectory: " + (runDir ? runDir : "(root)") + "\n\n" + (data.output || "");
      if (data.returncode !== 0) text += "\n[exit code " + data.returncode + "]";
      if (data.truncated) text += "\n[output truncated at 64 KB]";
      if (data.image) text += "\n[result: " + data.image + "]";
      showOutput(text, data.truncated);
      if (data.returncode === 0) toast(data.image ? "Done: " + data.image : "Done");
    }
  } catch (e) {
    record.status = "network";
    record.output = "";
    showOutput("$ " + rec.command + "\n\nNetwork error: " + e + "\nThe server may still be running the command; its state is unknown.", false);
  }
  record.dur = (Date.now() - t0) / 1000;
  addRecord(record);
  btn.disabled = false;
  btn.textContent = "Generate figure";
  running = false;
  loadScan(curPath, true);
}

function addRecord(record) {
  cmdHistory.unshift(record);
  if (cmdHistory.length > 40) cmdHistory.pop();
  renderHistory();
}
function renderHistory() {
  var list = el("histList");
  list.innerHTML = "";
  el("histCount").textContent = cmdHistory.length ? "(" + cmdHistory.length + ")" : "";
  if (!cmdHistory.length) {
    list.innerHTML = '<div class="empty">No commands have been run in this session. History is kept for this browser session only.</div>';
    return;
  }
  cmdHistory.forEach(function(h) {
    var d = document.createElement("div");
    d.className = "hist-item";
    var badge = h.status === "ok"
      ? '<span class="badge ok">ok</span>'
      : '<span class="badge">' + esc(h.status) + "</span>";
    var meta = (h.dir ? "in /" + esc(h.dir) : "in (root)") + " · " + h.start.toLocaleTimeString() + " · " + fmtDur(h.dur);
    d.innerHTML = badge + '<span class="hcmd">' + esc(h.cmd) + "</span>" +
      '<span class="hmeta">' + meta + "</span>";
    var acts = document.createElement("span");
    acts.className = "hacts";
    var outBtn = document.createElement("button");
    outBtn.className = "ghost small";
    outBtn.textContent = "output";
    outBtn.onclick = function() { showOutput("$ " + h.cmd + "\n\n" + h.output + (h.rc ? "\n[exit code " + h.rc + "]" : ""), h.truncated); };
    acts.appendChild(outBtn);
    var copyBtn = document.createElement("button");
    copyBtn.className = "ghost small";
    copyBtn.textContent = "copy";
    copyBtn.onclick = function() { copyText(h.cmd); };
    acts.appendChild(copyBtn);
    if (h.image) {
      var resBtn = document.createElement("button");
      resBtn.className = "ghost small";
      resBtn.textContent = "result";
      resBtn.title = "open " + h.image + " generated in " + (h.dir ? h.dir : "the root directory");
      resBtn.onclick = function() {
        setView("figures");
        openLightbox(joinRel(h.dir, h.image), h.image);
      };
      acts.appendChild(resBtn);
    }
    d.appendChild(acts);
    list.appendChild(d);
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
      dest = parentRel(curPath);
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
  if (busyFlag) termAppend("(note: the server reports another command is already running)\n");
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
      if (data.truncated) termAppend("[output truncated at 64 KB]\n");
      if (data.returncode !== 0) termAppend("[exit code " + data.returncode + "]\n");
      loadScan(curPath, true);
    } else {
      termAppend(" Error: " + (data.error || "command failed") + "\n");
    }
  } catch (e) {
    termAppend(" Network error: " + e + "\n The server may still be running the command; its state is unknown.\n");
  } finally {
    termBusy = false;
    input.disabled = false;
    input.focus();
  }
}

function pollTick() {
  if (document.hidden) { schedule(); return; }
  loadScan(curPath, true).then(schedule, schedule);
}
function schedule() {
  if (pollTimer) clearTimeout(pollTimer);
  pollTimer = null;
  var sel = parseInt(el("refreshSel").value, 10);
  if (sel <= 0) { el("sbRefresh").textContent = "auto off"; return; }
  var delay = pollActive ? sel * 1000 : 60000;
  el("sbRefresh").textContent = pollActive ? "active · auto " + sel + "s" : "idle · 60s check";
  pollTimer = setTimeout(pollTick, delay);
}

function applyTheme() {
  var dark = themeManual === "dark" ||
    (themeManual !== "light" && window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches);
  document.body.classList.toggle("dark", dark);
  var btn = el("themeBtn");
  if (btn) btn.innerHTML = dark ? ICON_SUN : ICON_MOON;
  if (viewMode === "overview" && lossData && !lossData.empty) drawLossChart();
  if (viewer && viewer.mode === "data") setViewerMode("data");
}
function toggleTheme() {
  var dark = document.body.classList.contains("dark");
  themeManual = dark ? "light" : "dark";
  try { localStorage.setItem("gk_theme", themeManual); } catch (e) {}
  applyTheme();
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
    renderBrowser({subdirs: lastSubdirs || [], files: lastFiles || []}, lastScanServerNow || Math.floor(Date.now() / 1000), {}, {});
  });
  el("sortSel").onchange = function() {
    renderBrowser({subdirs: lastSubdirs || [], files: lastFiles || []}, lastScanServerNow || Math.floor(Date.now() / 1000), {}, {});
  };
  el("upBtn").onclick = function() { loadScan(parentRel(curPath)); };
  el("copyPathBtn").onclick = function() {
    copyText((rootDir || "") + "/" + curPath);
  };
  el("curPathLabel").onclick = function() { copyText((rootDir || "") + "/" + curPath); };
  el("drawerBtn").onclick = function() { toggleDrawer(); };
  var viewtabs = document.querySelectorAll(".viewtab");
  for (var i = 0; i < viewtabs.length; i++) {
    (function(btn) {
      btn.onclick = function() { setView(btn.getAttribute("data-view")); };
    })(viewtabs[i]);
  }
  var ptabs = document.querySelectorAll(".ptab");
  for (var j = 0; j < ptabs.length; j++) {
    (function(btn) {
      btn.onclick = function() { setPanelTab(btn.getAttribute("data-tab")); };
    })(ptabs[j]);
  }
  el("panelToggle").onclick = function() {
    expandPanel(el("bottomPanel").classList.contains("closed"));
  };
  expandPanel(false);
  el("lightboxClose").onclick = closeLightbox;
  el("lightboxOverlay").onclick = function(ev) {
    if (ev.target === el("lightboxOverlay")) closeLightbox();
  };
  document.addEventListener("keydown", function(ev) {
    if (ev.key !== "Escape") return;
    var t = ev.target;
    if (t && t.closest && t.closest("input, textarea, select")) return;
    if (!el("lightboxOverlay").classList.contains("hidden")) closeLightbox();
  });
  el("termInput").addEventListener("keydown", function(ev) {
    if (ev.key === "Enter") {
      var value = el("termInput").value;
      el("termInput").value = "";
      termExec(value);
    }
  });
  el("termBody") && el("termBody").addEventListener("click", function() {
    var sel = window.getSelection ? String(window.getSelection()) : "";
    if (sel) return;
    if (!termBusy) el("termInput").focus();
  });
  el("themeBtn").onclick = toggleTheme;
  el("refreshSel").onchange = schedule;
  var mq = window.matchMedia ? window.matchMedia("(prefers-color-scheme: dark)") : null;
  if (mq && mq.addEventListener) {
    mq.addEventListener("change", function() {
      if (!themeManual) applyTheme();
    });
  }
  ageTimer = setInterval(updateConnText, 5000);
  applyTheme();
  loadScan("");
  schedule();
});
