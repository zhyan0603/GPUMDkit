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
  if (sec < 5) return "刚刚";
  if (sec < 60) return Math.floor(sec) + " 秒前";
  if (sec < 3600) return Math.floor(sec / 60) + " 分钟前";
  if (sec < 86400) return Math.floor(sec / 3600) + " 小时前";
  return Math.floor(sec / 86400) + " 天前";
}
function fmtDur(sec) {
  if (sec == null || !isFinite(sec)) return "--";
  if (sec < 60) return Math.floor(sec) + "s";
  if (sec < 3600) return Math.floor(sec / 60) + "m " + Math.floor(sec % 60) + "s";
  return Math.floor(sec / 3600) + "h " + Math.floor((sec % 3600) / 60) + "m";
}
function fmtGenRate(gps) {
  if (gps == null || !isFinite(gps) || gps <= 0) return "--";
  if (gps >= 1) return gps.toFixed(1) + " gen/s";
  if (gps * 60 >= 1) return (gps * 60).toFixed(1) + " gen/min";
  return (gps * 3600).toFixed(1) + " gen/h";
}
function fmtNum(v) {
  if (Math.abs(v) >= 1e5 || (Math.abs(v) < 1e-3 && v !== 0)) return v.toExponential(2);
  return String(Number(v.toPrecision(4)));
}
function fmtAxis(v) {
  if (v === 0) return "0";
  var e = Math.round(Math.log10(Math.abs(v)));
  if (Math.abs(v / Math.pow(10, e) - 1) < 1e-9) return "1e" + e;
  if (Math.abs(v) >= 1e4 || Math.abs(v) < 1e-2) return v.toExponential(1).replace("e+", "e");
  return String(Number(v.toPrecision(3)));
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
      toast("已复制");
    }, function() {
      toast(fallbackCopy(text) ? "已复制" : "复制失败");
    });
  } else {
    toast(fallbackCopy(text) ? "已复制" : "复制失败");
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

var ICON_SUN = '<svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round"><circle cx="12" cy="12" r="4.5"/><line x1="12" y1="2" x2="12" y2="4.5"/><line x1="12" y1="19.5" x2="12" y2="22"/><line x1="2" y1="12" x2="4.5" y2="12"/><line x1="19.5" y1="12" x2="22" y2="12"/><line x1="4.6" y1="4.6" x2="6.4" y2="6.4"/><line x1="17.6" y1="17.6" x2="19.4" y2="19.4"/><line x1="4.6" y1="19.4" x2="6.4" y2="17.6"/><line x1="17.6" y1="6.4" x2="19.4" y2="4.6"/></svg>';
var ICON_MOON = '<svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"><path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"/></svg>';
var ICON_CHEV = '<svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polyline points="6 9 12 15 18 9"/></svg>';

var PREVIEW_MAX_CLIENT = 512 * 1024;
var LOSS_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"];
var KNOWN_COLS = {
  "msd.out": {4: ["t", "msd_x", "msd_y", "msd_z"], 7: ["t", "msd_x", "msd_y", "msd_z", "sdc_x", "sdc_y", "sdc_z"]},
  "sdc.out": {4: ["t", "vac_x", "vac_y", "vac_z"]},
  "loss.out": {6: ["gen", "Loss", "E_tr", "F_tr", "V_tr"], 10: ["gen", "L_total", "L1", "L2", "E_tr", "F_tr", "V_tr", "E_te", "F_te", "V_te"]},
};
var REC_CATEGORIES = {
  plt_msd: "扩散与输运",
  plt_sdc: "扩散与输运",
  plt_msd_sdc: "扩散与输运",
  plt_vac: "扩散与输运",
  plt_msd_conv: "扩散与输运",
  plt_thermo: "热力学",
  plt_train: "NEP 训练",
  plt_train_density: "NEP 训练",
  plt_train_test: "NEP 训练",
  plt_prediction: "NEP 训练",
  plt_sigma: "Arrhenius 分析",
  plt_D: "Arrhenius 分析",
};
var REC_CAT_ORDER = ["扩散与输运", "NEP 训练", "热力学", "Arrhenius 分析"];
var FIG_PAGE = 12;

var curPath = "";
var rootDir = "";
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
var lastPollTs = 0;
var lastNepDone = null;
var nepSamples = [];
var cmdHistory = [];
var drawerTab = "output";
var sidebarOpen = true;
var pendingOpen = null;
var figShown = FIG_PAGE;
var lbList = [];
var lbIndex = -1;
var lbReturnFocus = null;
var viewer = null;
var lossView = { data: null, pending: false, reason: null, logX: true, logY: true, hidden: {}, fetchedAt: 0 };
var plotState = null;
var lastRecs = null;

try {
  themeManual = localStorage.getItem("gk_theme");
  var qp = new URLSearchParams(location.search).get("theme");
  if (qp === "dark" || qp === "light") themeManual = qp;
} catch (e) {}

function setConn(ok) {
  if (connOk !== ok) {
    connOk = ok;
    el("connWrap").classList.toggle("bad", !ok);
  }
  if (ok) lastGoodAt = Date.now();
  updateConnText();
}
function updateConnText() {
  if (connOk === null) { el("connText").textContent = "连接中..."; return; }
  if (!connOk) {
    var age = lastGoodAt ? fmtAgo((Date.now() - lastGoodAt) / 1000) : "从未";
    el("connText").textContent = "连接异常 · 上次更新 " + age;
    return;
  }
  el("connText").textContent = "更新于 " + fmtAgo((Date.now() - lastGoodAt) / 1000);
}

async function api(url, opts) {
  var resp;
  try {
    resp = await fetch(url, opts);
  } catch (e) {
    setConn(false);
    throw e;
  }
  if (resp.status === 401) showLogin("请先登录。");
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
      el("loginMsg").textContent = data.error || "登录失败。";
    }
  } catch (e) {
    el("loginMsg").textContent = "网络错误: " + e;
  }
  el("loginBtn").disabled = false;
}

async function loadScan(path, silent) {
  var seq = ++scanSeq;
  var resp;
  try {
    resp = await api("/api/scan?path=" + encodeURIComponent(path));
  } catch (e) {
    if (!silent) showBanner("网络错误");
    return false;
  }
  if (seq !== scanSeq) return false;
  if (resp.status === 401) return false;
  var data = null;
  try { data = await resp.json(); } catch (e) {}
  if (seq !== scanSeq) return false;
  if (!resp.ok || !data) {
    setConn(false);
    showBanner((data && data.error) || "读取目录失败");
    return false;
  }
  setConn(true);
  el("errorBanner").classList.add("hidden");
  hideLogin();
  var pathChanged = data.path !== lastPath;
  curPath = data.path || "";
  rootDir = data.root || "";
  lastPath = curPath;
  lastScanServerNow = data.now;
  busyFlag = !!data.busy;
  var nowS = data.now || Math.floor(Date.now() / 1000);
  var pollGap = lastPollTs ? (Date.now() / 1000 - lastPollTs) : 0;
  var prevF = {}, prevS = {};
  if (!pathChanged) {
    if (lastFiles) lastFiles.forEach(function(f) { prevF[f.name] = f; });
    if (lastSubdirs) lastSubdirs.forEach(function(s) { prevS[s.name] = s; });
  } else {
    nepSamples = [];
    lossView = { data: null, pending: false, reason: null, logX: true, logY: true, hidden: {}, fetchedAt: 0 };
    figShown = FIG_PAGE;
    plotState = null;
  }
  var prevNepDone = lastNepDone;
  var training = data.training;
  lastNepDone = training && !training.finished ? training.done : null;
  if (training && training.done != null && prevNepDone != null && training.done < prevNepDone) {
    nepSamples = [];
  }
  if (training && !training.finished) {
    nepSamples.push({done: training.done, ts: Date.now() / 1000});
    if (nepSamples.length > 8) nepSamples.shift();
  }
  lastPollTs = Date.now() / 1000;
  lastFiles = data.files || [];
  lastSubdirs = data.subdirs || [];
  lastRecs = data.recommendations || [];
  lastTraining = data.training || null;
  var liveFiles = (data.files || []).filter(function(f) {
    return prevF[f.name] && f.size > prevF[f.name].size;
  }).map(function(f) {
    var grew = f.size - prevF[f.name].size;
    return {name: f.name, grew: grew, rate: pollGap > 0 ? grew / pollGap : null, age: f.mtime != null ? nowS - f.mtime : null};
  });
  liveFiles = liveFiles.filter(function(f) { return f.age != null && f.age < 180; });
  var fresh = liveFiles.length > 0;
  if (!fresh) {
    var filesArr = data.files || [];
    for (var fi = 0; fi < filesArr.length; fi++) {
      if (filesArr[fi].mtime != null && nowS - filesArr[fi].mtime < 180) { fresh = true; break; }
    }
  }
  if (!fresh) {
    var subsArr = data.subdirs || [];
    for (var si = 0; si < subsArr.length; si++) {
      if (subsArr[si].newest != null && nowS - subsArr[si].newest < 180) { fresh = true; break; }
    }
  }
  pollActive = fresh || (data.training && !data.training.finished);
  renderCrumb(curPath);
  renderSidebar(data, nowS, prevS, prevF);
  renderCanvas(data, liveFiles, nowS);
  el("runningChip").classList.toggle("hidden", !(running || busyFlag));
  schedule();
  if (hasLossFile(data) && !lossView.data) fetchLoss(seq);
  return true;
}

function hasLossFile(data) {
  return (data.files || []).some(function(f) { return f.name === "loss.out"; });
}

async function fetchLoss(seq) {
  var path = curPath;
  lossView.pending = true;
  lossView.reason = null;
  try {
    var resp = await api("/api/loss?path=" + encodeURIComponent(path));
    var data = null;
    try { data = await resp.json(); } catch (e) {}
    if (seq !== scanSeq || path !== curPath) return;
    lossView.pending = false;
    if (resp.ok && data) {
      if (data.empty) {
        lossView.data = null;
        lossView.reason = data.reason || "empty";
      } else {
        lossView.data = data;
        lossView.fetchedAt = Date.now() / 1000;
      }
    } else {
      lossView.data = null;
      lossView.reason = "error";
    }
  } catch (e) {
    if (seq === scanSeq && path === curPath) {
      lossView.pending = false;
      lossView.reason = "error";
    }
  }
  if (seq === scanSeq && path === curPath) renderCanvasRoot();
}

function showBanner(msg) {
  el("errorText").textContent = msg + " · 显示的可能是旧数据";
  el("errorBanner").classList.remove("hidden");
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

function renderSidebar(data, nowS, prevS, prevF) {
  var fl = el("fileList");
  fl.innerHTML = "";
  var filter = el("fileFilter").value.trim().toLowerCase();
  var dirs = (data.subdirs || []).slice();
  var files = (data.files || []).slice().sort(function(a, b) {
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
      row.title = "最近有写入";
    }
    row.appendChild(dot);
    var nm = document.createElement("span");
    nm.className = "bname";
    nm.textContent = s.name + "/";
    row.appendChild(nm);
    var meta = document.createElement("span");
    meta.className = "bmeta";
    var bits = [];
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
      row.title = "最近有写入";
    }
    row.appendChild(dot);
    var a = document.createElement("a");
    a.href = "#";
    a.textContent = f.name;
    a.onclick = function(ev) { ev.preventDefault(); fileClicked(f, a); };
    row.appendChild(a);
    var meta = document.createElement("span");
    meta.className = "bmeta";
    var bits = [];
    if (f.mtime != null) bits.push(fmtAgo(nowS - f.mtime));
    meta.textContent = bits.join(" · ");
    row.appendChild(meta);
    fl.appendChild(row);
  });
  if (!shown) {
    var empty = document.createElement("div");
    empty.className = "muted";
    empty.style.padding = "8px";
    empty.textContent = filter ? "无匹配文件" : "空目录";
    fl.appendChild(empty);
  }
}

function fileClicked(f, source) {
  var ext = f.name.substring(f.name.lastIndexOf(".") + 1).toLowerCase();
  if (ext === "png" || ext === "jpg" || ext === "jpeg") {
    openLightbox(f, source);
  } else {
    openPreview(f);
  }
}

function renderCanvas(data, liveFiles, nowS) {
  var title = curPath ? curPath.split("/").pop() : (rootDir ? rootDir.split("/").pop() : "root");
  el("dirTitle").textContent = title;
  el("dirTitle").title = (rootDir || "") + "/" + curPath;
  var bits = [((data.dirs || []).length + (data.files || []).length) + " 项"];
  bits.push("更新于 " + fmtAgo(0));
  if (liveFiles.length) bits.push(liveFiles.length + " 个文件在增长");
  el("statusLine").textContent = bits.join(" · ");
  renderCanvasRoot();
}

function buildLossSection() {
  var sec = document.createElement("section");
  sec.className = "sec";
  var head = document.createElement("div");
  head.className = "sec-head";
  head.innerHTML = "<h2>训练曲线</h2><span class='sec-meta'>loss.out</span>";
  sec.appendChild(head);
  if (lossView.pending) {
    var loading = document.createElement("div");
    loading.className = "muted";
    loading.style.marginTop = "10px";
    loading.textContent = "正在读取 loss.out...";
    sec.appendChild(loading);
    return sec;
  }
  if (!lossView.data) {
    var reason = lossView.reason === "unreadable" ? "无法读取" : (lossView.reason === "empty" ? "文件为空" : "无有效记录");
    var note = document.createElement("div");
    note.className = "chart-note";
    note.style.marginTop = "10px";
    note.innerHTML = "loss.out " + reason + " · ";
    var openBtn = document.createElement("button");
    openBtn.className = "note-link";
    openBtn.textContent = "打开文件";
    openBtn.onclick = function() {
      var f = (lastFiles || []).find(function(x) { return x.name === "loss.out"; });
      if (f) openPreview(f);
    };
    note.appendChild(openBtn);
    sec.appendChild(note);
    return sec;
  }
  var wrap = document.createElement("div");
  wrap.className = "losswrap";
  wrap.id = "lossWrap";
  var chartCol = document.createElement("div");
  chartCol.className = "losschart";
  var canvas = document.createElement("canvas");
  canvas.className = "chart";
  canvas.id = "lossCanvas";
  chartCol.appendChild(canvas);
  var bar = document.createElement("div");
  bar.className = "chartbar";
  var summary = document.createElement("span");
  summary.className = "scale-summary";
  summary.id = "scaleSummary";
  bar.appendChild(summary);
  var scaleBtn = document.createElement("button");
  scaleBtn.className = "linkbtn";
  scaleBtn.textContent = "坐标";
  scaleBtn.setAttribute("aria-label", "展开坐标设置");
  scaleBtn.onclick = function() {
    panel.classList.toggle("hidden");
  };
  bar.appendChild(scaleBtn);
  var spacer = document.createElement("span");
  spacer.className = "grow";
  bar.appendChild(spacer);
  var zoomBtn = document.createElement("button");
  zoomBtn.className = "linkbtn";
  zoomBtn.textContent = "放大查看";
  zoomBtn.onclick = openChartZoom;
  bar.appendChild(zoomBtn);
  chartCol.appendChild(bar);
  var panel = document.createElement("div");
  panel.className = "scalepanel hidden";
  panel.id = "scalePanel";
  panel.appendChild(scaleControl("X", "logX"));
  panel.appendChild(scaleControl("Y", "logY"));
  chartCol.appendChild(panel);
  var legend = document.createElement("div");
  legend.className = "lg-row";
  legend.id = "lossLegend";
  chartCol.appendChild(legend);
  var hover = document.createElement("div");
  hover.className = "hoverread";
  hover.id = "lossHover";
  chartCol.appendChild(hover);
  var notes = document.createElement("div");
  notes.id = "lossNotes";
  chartCol.appendChild(notes);
  wrap.appendChild(chartCol);
  var sumBox = document.createElement("div");
  sumBox.className = "losssummary";
  sumBox.id = "lossSummary";
  wrap.appendChild(sumBox);
  sec.appendChild(wrap);
  return sec;
}

function scaleControl(axis, key) {
  var wrap = document.createElement("span");
  wrap.className = "row";
  wrap.style.gap = "6px";
  var lab = document.createElement("span");
  lab.textContent = axis;
  wrap.appendChild(lab);
  var seg = document.createElement("span");
  seg.className = "seg";
  ["Linear", "Log"].forEach(function(mode) {
    var b = document.createElement("button");
    b.className = "segbtn" + ((mode === "Log") === !!lossView[key] ? " on" : "");
    b.textContent = mode;
    b.setAttribute("data-mode", mode);
    b.onclick = function() {
      lossView[key] = mode === "Log";
      updateScaleSummary();
      seg.querySelectorAll(".segbtn").forEach(function(x) {
        x.classList.toggle("on", (x.getAttribute("data-mode") === "Log") === !!lossView[key]);
      });
      sizeAndDrawLoss();
    };
    seg.appendChild(b);
  });
  wrap.appendChild(seg);
  return wrap;
}

function updateScaleSummary() {
  var s = el("scaleSummary");
  if (s) s.textContent = "X " + (lossView.logX ? "Log" : "Linear") + " · Y " + (lossView.logY ? "Log" : "Linear");
}

function lossChartSize(availW) {
  if (availW >= 860) {
    var w = Math.min(availW - 230, 760);
    var h = Math.max(300, Math.min(420, Math.round(w / 1.85)));
    return {w: w, h: h, stacked: false};
  }
  var w2 = Math.min(availW, 760);
  var h2 = Math.max(280, Math.min(460, Math.round(w2 / 1.5)));
  return {w: w2, h: h2, stacked: true};
}

function sizeAndDrawLoss() {
  var canvas = el("lossCanvas");
  var wrap = el("lossWrap");
  if (!canvas || !wrap) return;
  var avail = (el("contentCol") || document.body).clientWidth - 64;
  if (avail < 80) return;
  var geo = lossChartSize(avail);
  wrap.classList.toggle("stacked", geo.stacked);
  canvas.style.width = geo.w + "px";
  canvas.style.height = geo.h + "px";
  drawLossTo(canvas, geo.w, geo.h);
  renderLossSummary();
  renderLossLegend("lossLegend");
  renderLossNotes();
  updateScaleSummary();
}

function visibleLossSeries() {
  var d = lossView.data;
  if (!d) return [];
  var out = [];
  (d.series || []).forEach(function(s, i) {
    if (lossView.hidden[s.label]) return;
    out.push({xs: d.xs, ys: s.values, color: LOSS_COLORS[i % LOSS_COLORS.length], label: s.label});
  });
  return out;
}

function drawLossTo(canvas, w, h) {
  var series = visibleLossSeries();
  var d = lossView.data;
  drawSeries(canvas, series, {
    logX: lossView.logX, logY: lossView.logY,
    xLabel: d && d.x_mode === "generation" ? "generation" : "record #",
    yLabel: "Loss functions"
  });
  attachHover(canvas, el("lossHover"), function(idx) {
    var parts = ["x=" + fmtNum(d.xs[idx])];
    series.forEach(function(s) {
      parts.push(s.label + "=" + fmtNum(s.ys[idx]));
    });
    return "第 " + (idx + 1) + " 个采样点 · " + parts.join(" · ");
  });
}

function renderLossSummary() {
  var box = el("lossSummary");
  if (!box || !lossView.data) return;
  var d = lossView.data;
  var t = null;
  if (lastTraining) t = lastTraining;
  box.innerHTML = "";
  function line(label, value) {
    var l = document.createElement("div");
    l.className = "sline";
    l.innerHTML = "<span class='slabel'>" + label + "</span><span class='sval'>" + value + "</span>";
    box.appendChild(l);
  }
  line("记录", d.count.toLocaleString() + " 条");
  line("最新", "代数 " + (d.last_gen != null ? d.last_gen.toLocaleString() : "--"));
  if (d.multi_run) {
    line("范围", "多次训练记录（最新段自 " + d.seg_start_gen + " 起）");
  }
  if (t && t.has_target) {
    line("目标", t.total.toLocaleString());
    if (!t.multi_run) {
      var pct = Math.min(100, 100 * t.done / t.total);
      var track = document.createElement("div");
      track.className = "bar-track";
      var fill = document.createElement("div");
      fill.className = "bar-fill" + (t.finished ? " done" : "");
      fill.style.width = pct + "%";
      track.appendChild(fill);
      box.appendChild(track);
      var rate = trainingRate();
      if (t.finished) {
        line("状态", "已完成 · 达到目标代数");
      } else {
        line("速率", fmtGenRate(rate.rate));
        if (rate.estimating) {
          line("预计", "估算中（需更多刷新样本）");
        } else if (rate.rate) {
          line("预计剩余", fmtDur((t.total - t.done) / rate.rate));
        }
        line("状态", t.loss_mtime != null ? "最后写入 " + fmtAgo(lastScanServerNow - t.loss_mtime) : "--");
      }
    } else {
      line("状态", "多次训练记录，进度按最新段");
    }
  } else if (t) {
    line("目标", "未指定（nep.in 无 generation 行）");
    line("状态", t.loss_mtime != null ? "最后写入 " + fmtAgo(lastScanServerNow - t.loss_mtime) : "--");
  } else {
    line("目标", "未发现 nep.in，不显示进度");
    line("状态", d.mtime != null ? "最后写入 " + fmtAgo(lastScanServerNow - d.mtime) : "--");
  }
}
var lastTraining = null;

function renderLossLegend(legendId) {
  var legend = el(legendId);
  if (!legend || !lossView.data) return;
  legend.innerHTML = "";
  (lossView.data.series || []).forEach(function(s, i) {
    var item = document.createElement("button");
    item.className = "lg-item" + (lossView.hidden[s.label] ? " off" : "");
    var dot = document.createElement("span");
    dot.className = "lg-dot";
    dot.style.background = LOSS_COLORS[i % LOSS_COLORS.length];
    item.appendChild(dot);
    item.appendChild(document.createTextNode(s.label));
    item.title = lossView.hidden[s.label] ? "点击显示该系列" : "点击隐藏该系列";
    item.onclick = function() {
      lossView.hidden[s.label] = !lossView.hidden[s.label];
      renderLossLegend(legendId);
      if (legendId === "lossLegend") sizeAndDrawLoss();
      else drawZoomChart();
    };
    legend.appendChild(item);
  });
}

function renderLossNotes() {
  var notes = el("lossNotes");
  if (!notes || !lossView.data) return;
  notes.innerHTML = "";
  var d = lossView.data;
  if (d.sampled) {
    notes.appendChild(noteLine("曲线降采样至 " + d.points + " / " + d.count + " 条记录"));
  }
  if (d.skipped > 0) {
    notes.appendChild(noteLine(d.skipped + " 行无效或不完整记录被跳过"));
  }
  if (d.multi_run) {
    notes.appendChild(noteLine("loss.out 包含多轮训练（代数回退），进度按最新一段"));
  }
  if (lossView.logY) {
    var nonpos = 0;
    (d.series || []).forEach(function(s) {
      if (lossView.hidden[s.label]) return;
      s.values.forEach(function(v) { if (v <= 0) nonpos++; });
    });
    if (nonpos > 0) {
      var line = noteLine(nonpos + " 个非正值未在 Log 纵轴下显示");
      var sw = document.createElement("button");
      sw.className = "note-link";
      sw.textContent = "切换 Y 为 Linear";
      sw.onclick = function() { lossView.logY = false; sizeAndDrawLoss(); };
      line.appendChild(sw);
      notes.appendChild(line);
    }
  }
}
function noteLine(text) {
  var d = document.createElement("div");
  d.className = "chart-note";
  d.textContent = text;
  return d;
}

function openChartZoom() {
  el("chartZoom").classList.remove("hidden");
  drawZoomChart();
}
function drawZoomChart() {
  var canvas = el("zoomCanvas");
  if (!canvas || !lossView.data) return;
  var w = Math.min(window.innerWidth * 0.92, 1150);
  var h = Math.min(window.innerHeight * 0.72, 620);
  canvas.style.width = w + "px";
  canvas.style.height = h + "px";
  drawLossTo(canvas, w, h);
  renderLossLegend("zoomLegend");
  el("zoomScale").textContent = "X " + (lossView.logX ? "Log" : "Linear") + " · Y " + (lossView.logY ? "Log" : "Linear");
}

function buildFiguresSection(imgs) {
  var sec = document.createElement("section");
  sec.className = "sec";
  var head = document.createElement("div");
  head.className = "sec-head";
  head.innerHTML = "<h2>图片</h2><span class='sec-meta'>" + imgs.length + " 张</span>";
  sec.appendChild(head);
  var grid = document.createElement("div");
  grid.className = "figgrid" + (imgs.length === 1 ? " single" : "");
  grid.id = "figGrid";
  var show = imgs.slice(0, figShown);
  show.forEach(function(f) {
    grid.appendChild(figThumb(f));
  });
  sec.appendChild(grid);
  if (imgs.length > figShown) {
    var more = document.createElement("div");
    more.className = "figmores";
    var btn = document.createElement("button");
    btn.className = "ghost";
    btn.textContent = "显示更多（" + (imgs.length - figShown) + "）";
    btn.onclick = function() {
      figShown += FIG_PAGE;
      renderCanvasRoot();
    };
    more.appendChild(btn);
    sec.appendChild(more);
  }
  return sec;
}
function figThumb(f) {
  var cell = document.createElement("button");
  cell.className = "fig-thumb";
  var img = document.createElement("img");
  img.src = "/api/image?path=" + encodeURIComponent(joinRel(curPath, f.name));
  img.alt = f.name;
  img.loading = "lazy";
  var cap = document.createElement("div");
  cap.className = "fig-cap";
  cap.textContent = f.name;
  cap.title = f.name;
  cell.appendChild(img);
  cell.appendChild(cap);
  cell.onclick = function() { openLightbox(f, cell); };
  return cell;
}

function openLightbox(f, source) {
  lbList = (lastFiles || []).filter(function(x) { return /\.(png|jpe?g)$/i.test(x.name); });
  lbIndex = lbList.findIndex(function(x) { return x.name === f.name; });
  lbReturnFocus = source || null;
  showLightboxAt(lbIndex);
  el("lightbox").classList.remove("hidden");
  el("lbClose").focus();
}
function showLightboxAt(idx) {
  if (idx < 0 || idx >= lbList.length) return;
  lbIndex = idx;
  var f = lbList[idx];
  el("lbImg").src = "/api/image?path=" + encodeURIComponent(joinRel(curPath, f.name));
  el("lbTitle").textContent = f.name;
  el("lbCount").textContent = (idx + 1) + " / " + lbList.length;
}
function lightboxStep(delta) {
  if (lbIndex < 0) return;
  var next = lbIndex + delta;
  if (next < 0 || next >= lbList.length) return;
  showLightboxAt(next);
}
function closeLightbox() {
  el("lightbox").classList.add("hidden");
  el("lbImg").src = "";
  if (lbReturnFocus && lbReturnFocus.focus) lbReturnFocus.focus();
  lbReturnFocus = null;
}

function buildAnalysesSection(recs, files, nowS) {
  var fileMap = {};
  files.forEach(function(f) { fileMap[f.name] = f; });
  var sec = document.createElement("section");
  sec.className = "sec";
  var head = document.createElement("div");
  head.className = "sec-head";
  head.innerHTML = "<h2>可用分析</h2><span class='sec-meta'>" + recs.length + " 项</span>";
  sec.appendChild(head);
  var groups = {};
  recs.forEach(function(rec) {
    var cat = REC_CATEGORIES[rec.action] || "其他";
    if (!groups[cat]) groups[cat] = [];
    groups[cat].push(rec);
  });
  var cats = REC_CAT_ORDER.concat(Object.keys(groups).filter(function(c) { return REC_CAT_ORDER.indexOf(c) < 0 && c !== "其他"; }));
  if (groups["其他"]) cats.push("其他");
  cats.forEach(function(cat) {
    var list = groups[cat];
    if (!list) return;
    var box = document.createElement("div");
    box.className = "actrows";
    if (cat !== "扩散与输运" || list.length > 2) {
      var lab = document.createElement("div");
      lab.className = "chart-note";
      lab.style.margin = "10px 2px 4px";
      lab.textContent = cat;
      box.appendChild(lab);
    }
    list.forEach(function(rec) {
      box.appendChild(actRow(rec, fileMap, nowS));
    });
    sec.appendChild(box);
  });
  return sec;
}

function actRow(rec, fileMap, nowS) {
  var row = document.createElement("div");
  row.className = "actrow";
  var head = document.createElement("div");
  head.className = "actrow-head";
  var desc = document.createElement("span");
  desc.className = "adesc";
  desc.textContent = descOf(rec.action);
  head.appendChild(desc);
  var inp = document.createElement("span");
  inp.className = "ain";
  inp.textContent = rec.evidence.join(", ");
  head.appendChild(inp);
  var res = document.createElement("span");
  res.className = "ares";
  var outFile = fileMap[rec.produces];
  if (outFile && outFile.mtime != null) {
    res.textContent = "结果 " + fmtAgo(nowS - outFile.mtime);
    if (rec.stale) res.classList.add("stale");
    res.title = rec.stale ? "结果早于输入文件，可能需要重新生成" : rec.produces;
  } else {
    res.textContent = "未生成";
  }
  head.appendChild(res);
  var chev = document.createElement("span");
  chev.className = "achev";
  chev.innerHTML = ICON_CHEV;
  head.appendChild(chev);
  head.onclick = function() { row.classList.toggle("open"); };
  row.appendChild(head);
  var body = document.createElement("div");
  body.className = "actrow-body";
  var cmdrow = document.createElement("div");
  cmdrow.className = "cmdrow";
  var code = document.createElement("code");
  code.className = "rec-cmd";
  code.textContent = rec.command;
  cmdrow.appendChild(code);
  var copyBtn = document.createElement("button");
  copyBtn.className = "ghost small";
  copyBtn.textContent = "复制";
  copyBtn.onclick = function(ev) { ev.stopPropagation(); copyText(rec.command); };
  cmdrow.appendChild(copyBtn);
  body.appendChild(cmdrow);
  var actbar = document.createElement("div");
  actbar.className = "actbar";
  var run = document.createElement("button");
  run.className = "ghost small";
  run.textContent = "生成图像";
  run.disabled = busyFlag || running;
  run.title = busyFlag ? "服务器正在执行其他命令" : "在当前目录运行 " + rec.command;
  run.onclick = function(ev) { ev.stopPropagation(); runAction(rec, run); };
  actbar.appendChild(run);
  if (outFile) {
    var open = document.createElement("button");
    open.className = "ghost small";
    open.textContent = "查看结果";
    open.onclick = function(ev) {
      ev.stopPropagation();
      openLightbox(outFile, open);
    };
    actbar.appendChild(open);
  }
  body.appendChild(actbar);
  row.appendChild(body);
  return row;
}

function descOf(action) {
  var map = {
    plt_msd: "MSD 曲线",
    plt_sdc: "SDC 曲线",
    plt_msd_sdc: "MSD 与 SDC",
    plt_vac: "VAC 曲线",
    plt_thermo: "热力学量曲线",
    plt_train: "训练损失与 parity",
    plt_train_density: "训练 parity 密度图",
    plt_train_test: "训练/测试对比",
    plt_prediction: "预测 parity 图",
    plt_msd_conv: "MSD 收敛性检查",
    plt_sigma: "电导率 Arrhenius 图",
    plt_D: "扩散系数 Arrhenius 图",
  };
  return map[action] || action;
}

function buildFilesSection(files) {
  var sec = document.createElement("section");
  sec.className = "sec";
  var head = document.createElement("div");
  head.className = "sec-head";
  head.innerHTML = "<h2>文件</h2><span class='sec-meta'>" + files.length + " 个</span>";
  sec.appendChild(head);
  var box = document.createElement("div");
  var nowS = lastScanServerNow || 0;
  var show = files.slice(0, 12);
  show.forEach(function(f) {
    var row = document.createElement("div");
    row.className = "filerow";
    var a = document.createElement("a");
    a.href = "#";
    a.className = "fname";
    a.textContent = f.name;
    a.onclick = function(ev) { ev.preventDefault(); openPreview(f); };
    row.appendChild(a);
    var meta = document.createElement("span");
    meta.className = "fmeta";
    meta.textContent = [fmtSize(f.size), f.mtime != null ? fmtAgo(nowS - f.mtime) : ""].filter(Boolean).join(" · ");
    row.appendChild(meta);
    box.appendChild(row);
  });
  if (files.length > 12) {
    var more = document.createElement("div");
    more.className = "muted";
    more.style.padding = "6px 10px";
    more.textContent = "以及另外 " + (files.length - 12) + " 个文件，见左侧文件栏";
    box.appendChild(more);
  }
  sec.appendChild(box);
  return sec;
}

function buildSubdirsSection(subdirs) {
  var sec = document.createElement("section");
  sec.className = "sec";
  var head = document.createElement("div");
  head.className = "sec-head";
  head.innerHTML = "<h2>子目录</h2><span class='sec-meta'>" + subdirs.length + " 个</span>";
  sec.appendChild(head);
  var nowS = lastScanServerNow || 0;
  subdirs.forEach(function(s) {
    var row = document.createElement("div");
    row.className = "subrow";
    var nm = document.createElement("span");
    nm.className = "sname";
    nm.textContent = s.name + "/";
    row.appendChild(nm);
    var meta = document.createElement("span");
    meta.className = "smeta";
    var bits = [s.count + " 项"];
    if (s.newest != null) bits.push(fmtAgo(nowS - s.newest));
    meta.textContent = bits.join(" · ");
    row.appendChild(meta);
    row.onclick = function() { loadScan(joinRel(curPath, s.name)); };
    sec.appendChild(row);
  });
  return sec;
}

function renderCanvasRoot() {
  var root = el("canvasRoot");
  root.innerHTML = "";
  var data = { files: lastFiles || [], subdirs: lastSubdirs || [] };
  var imgs = (data.files || []).filter(function(f) { return /\.(png|jpe?g)$/i.test(f.name); });
  var hasLoss = hasLossFile(data);
  var recs = lastRecs || [];
  var hasAny = (data.files || []).length + (data.subdirs || []).length;

  if (!hasAny) {
    var empty = document.createElement("div");
    empty.className = "bigempty";
    empty.innerHTML = "<p>此目录为空</p>";
    if (curPath) {
      var up = document.createElement("button");
      up.className = "ghost";
      up.textContent = "返回上级";
      up.onclick = function() { loadScan(parentRel(curPath)); };
      empty.appendChild(up);
    }
    root.appendChild(empty);
    return;
  }
  if (hasLoss) root.appendChild(buildLossSection());
  if (imgs.length) root.appendChild(buildFiguresSection(imgs));
  if (recs.length) root.appendChild(buildAnalysesSection(recs, data.files || [], lastScanServerNow || 0));
  var plainFiles = !hasLoss && !imgs.length;
  if (plainFiles && (data.files || []).length) root.appendChild(buildFilesSection(data.files || []));
  if (!hasLoss && !imgs.length && (data.subdirs || []).length && !(data.files || []).length) {
    root.appendChild(buildSubdirsSection(data.subdirs || []));
  }
  if (hasLoss) sizeAndDrawLoss();
}

function viewerDirty() {
  return !!(viewer && viewer.draft != null && viewer.draft !== viewer.text);
}

function openPreview(f) {
  if (viewerDirty()) {
    pendingOpen = f;
    el("pvConfirm").classList.remove("hidden");
    return;
  }
  loadPreview(f);
}

function loadPreview(f) {
  var rel = joinRel(curPath, f.name);
  viewer = {
    rel: rel, name: f.name, size: f.size, baseMtime: null,
    truncated: false, text: null, draft: null, mode: "preview", editing: false
  };
  el("pvConfirm").classList.add("hidden");
  el("pvConflict").classList.add("hidden");
  el("pvMsg").textContent = "";
  el("pvBody").innerHTML = '<div class="muted">读取中...</div>';
  el("pvTitle").textContent = f.name;
  el("pvMeta").textContent = fmtSize(f.size);
  el("previewLayer").classList.remove("hidden");
  el("pvClose").focus();
  var oversized = f.size != null && f.size > PREVIEW_MAX_CLIENT;
  var url = "/api/file?path=" + encodeURIComponent(rel) + (oversized ? "&partial=1" : "");
  api(url).then(function(resp) {
    if (resp.status === 401) { closePreviewReset(); return; }
    if (viewer === null || viewer.rel !== rel) return;
    if (resp.ok) {
      var ct = resp.headers.get("Content-Type") || "";
      if (ct.indexOf("application/json") >= 0) {
        return resp.json().then(function(j) { applyPartial(rel, j); });
      }
      return resp.text().then(function(text) {
        if (viewer === null || viewer.rel !== rel) return;
        viewer.baseMtime = parseInt(resp.headers.get("X-File-Mtime") || "0", 10) || null;
        viewer.text = text;
        renderPreview();
      });
    }
    if (resp.status === 413) {
      return fetchPartial(rel).then(function(j) {
        if (viewer === null || viewer.rel !== rel) return;
        if (!j) { previewError(rel, "文件过大，无法预览"); return; }
        applyPartial(rel, j);
      });
    }
    var ct2 = resp.headers.get("Content-Type") || "";
    if (ct2.indexOf("application/json") >= 0) {
      return resp.json().then(function(j) { previewError(rel, (j && j.error) || "预览失败"); });
    }
    previewError(rel, "预览失败");
  }).catch(function() {
    if (viewer !== null && viewer.rel === rel) previewError(rel, "网络错误");
  });
}

function applyPartial(rel, j) {
  if (viewer === null || viewer.rel !== rel) return;
  viewer.truncated = true;
  viewer.baseMtime = j.mtime;
  viewer.text = j.head + "\n[... 已省略 " + (j.size - j.head_bytes - j.tail_bytes).toLocaleString() + " 字节 ...]\n" + j.tail;
  viewer.size = j.size;
  el("pvMeta").textContent = fmtSize(j.size);
  renderPreview();
}
function previewError(rel, msg) {
  el("pvBody").innerHTML = "";
  el("pvMsg").textContent = msg;
}
function fetchPartial(rel) {
  return api("/api/file?path=" + encodeURIComponent(rel) + "&partial=1").then(function(resp) {
    if (resp.ok) return resp.json();
    return null;
  }).catch(function() { return null; });
}
function closePreviewReset() {
  el("pvBody").innerHTML = '<div class="muted">选择左侧文件以预览。</div>';
}

function closePreview() {
  if (viewerDirty()) {
    el("pvConfirm").classList.remove("hidden");
    return;
  }
  el("previewLayer").classList.add("hidden");
  viewer = null;
  lbReturnFocus && lbReturnFocus.focus && lbReturnFocus.focus();
}

function renderPreview() {
  if (!viewer) return;
  var body = el("pvBody");
  body.innerHTML = "";
  el("pvMsg").textContent = "";
  el("pvDirty").classList.toggle("hidden", !viewerDirty());
  var text = viewer.draft != null ? viewer.draft : viewer.text;
  var numeric = viewer.truncated ? null : parseNumeric(text || "");
  var kv = viewer.truncated ? [] : parseKV(text || "");
  var kvCount = kv.filter(function(it) { return it.type === "kv"; }).length;
  var nonBlank = kv.filter(function(it) { return it.type !== "blank"; }).length;
  var isForm = !viewer.truncated && kvCount >= 2 && kvCount * 2 >= nonBlank;
  var modes = [];
  if (numeric) modes.push(["preview", numeric ? "预览" : "预览"]);
  modes.push(["text", "文本"]);
  if (isForm) modes.push(["form", "表单"]);
  if (!viewer.truncated) el("pvEditBtn").classList.remove("hidden");
  else el("pvEditBtn").classList.add("hidden");
  var modeBar = el("pvModes");
  modeBar.innerHTML = "";
  modeBar.classList.toggle("hidden", modes.length <= 1 && !numeric);
  modes.forEach(function(m) {
    var b = document.createElement("button");
    b.className = "segbtn";
    b.textContent = m[1];
    b.setAttribute("data-mode", m[0]);
    b.onclick = function() { viewer.mode = m[0]; renderPreview(); };
    modeBar.appendChild(b);
  });
  if (viewer.mode === "plot") viewer.mode = numeric ? "plot" : "preview";
  if (viewer.mode === "form" && !isForm) viewer.mode = "preview";
  modeBar.querySelectorAll(".segbtn").forEach(function(b) {
    b.classList.toggle("on", b.getAttribute("data-mode") === viewer.mode);
  });
  if (viewer.mode === "form") renderFormMode(body);
  else if (viewer.editing) renderEditMode(body);
  else if (viewer.mode === "preview" && numeric) renderTableMode(body, numeric);
  else renderTextMode(body);
}

function renderTableMode(body, numeric) {
  var names = colNames(viewer.name, numeric.cols);
  var wrap = document.createElement("div");
  wrap.className = "table-scroll";
  var table = document.createElement("table");
  table.className = "data";
  var thead = document.createElement("thead");
  var htr = document.createElement("tr");
  for (var c = 0; c < numeric.cols; c++) {
    var th = document.createElement("th");
    th.textContent = names[c];
    htr.appendChild(th);
  }
  thead.appendChild(htr);
  table.appendChild(thead);
  var tbody = document.createElement("tbody");
  var maxRows = Math.min(numeric.rows.length, 100);
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
    note.textContent = "显示前 100 行（共 " + numeric.rows.length + " 行）";
    body.appendChild(note);
  }
  var opts = document.createElement("div");
  opts.className = "plotopts";
  opts.id = "plotOpts";
  var head = document.createElement("div");
  head.className = "plotopts-head";
  head.innerHTML = "<span class='achev'>" + ICON_CHEV + "</span> 绘图";
  head.onclick = function() { opts.classList.toggle("open"); };
  opts.appendChild(head);
  var optsBody = document.createElement("div");
  optsBody.className = "plotopts-body";
  opts.appendChild(optsBody);
  body.appendChild(opts);
  optsBody.appendChild(buildFilePlot(numeric, names));
}

function buildFilePlot(numeric, names) {
  var holder = document.createElement("div");
  var state = plotState && plotState.rel === viewer.rel ? plotState : {rel: viewer.rel, x: 0, ys: [1], logX: false, logY: false, open: false};
  plotState = state;
  var ctl = document.createElement("div");
  ctl.className = "chartctl";
  var xSel = document.createElement("select");
  xSel.setAttribute("aria-label", "X 列");
  for (var c = 0; c < numeric.cols; c++) {
    var o = document.createElement("option");
    o.value = c;
    o.textContent = names[c];
    if (c === state.x) o.selected = true;
    xSel.appendChild(o);
  }
  xSel.onchange = function() { state.x = parseInt(xSel.value, 10); redraw(); };
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
        redraw();
      };
      lab.appendChild(cb);
      lab.appendChild(document.createTextNode(names[ci]));
      ctl.appendChild(lab);
    })(c2);
  }
  ["logX", "logY"].forEach(function(key) {
    var lab = document.createElement("label");
    lab.className = "chk";
    var cb = document.createElement("input");
    cb.type = "checkbox";
    cb.checked = state[key];
    cb.onchange = function() { state[key] = cb.checked; redraw(); };
    lab.appendChild(cb);
    lab.appendChild(document.createTextNode(key === "logX" ? "log X" : "log Y"));
    ctl.appendChild(lab);
  });
  holder.appendChild(ctl);
  var canvas = document.createElement("canvas");
  canvas.className = "chart";
  canvas.style.width = "100%";
  canvas.style.height = "300px";
  holder.appendChild(canvas);
  var readout = document.createElement("div");
  readout.className = "hoverread";
  holder.appendChild(readout);
  function redraw() {
    var xs = numeric.rows.map(function(r) { return r[state.x]; });
    var series = state.ys.map(function(ci, i) {
      return {xs: xs, ys: numeric.rows.map(function(r) { return r[ci]; }), color: LOSS_COLORS[i % LOSS_COLORS.length], label: names[ci]};
    });
    drawSeries(canvas, series, {logX: state.logX, logY: state.logY, xLabel: names[state.x], yLabel: "value"});
    attachHover(canvas, readout, function(idx) {
      var parts = [names[state.x] + "=" + fmtNum(numeric.rows[idx][state.x])];
      state.ys.forEach(function(ci) { parts.push(names[ci] + "=" + fmtNum(numeric.rows[idx][ci])); });
      return "第 " + (idx + 1) + " 行 · " + parts.join(" · ");
    });
  }
  redraw();
  return holder;
}

function renderTextMode(body) {
  var pre = document.createElement("pre");
  pre.className = "preview";
  pre.textContent = viewer.draft != null ? viewer.draft : (viewer.text || "");
  body.appendChild(pre);
  if (viewer.truncated) {
    var note = document.createElement("p");
    note.className = "chart-note";
    note.textContent = "有界预览：仅显示开头 128 KB 与结尾 64 KB；编辑已禁用";
    body.appendChild(note);
  }
}

function renderEditMode(body) {
  var ta = document.createElement("textarea");
  ta.className = "editor";
  ta.value = viewer.draft != null ? viewer.draft : (viewer.text || "");
  ta.addEventListener("input", function() { setDraft(ta.value); });
  body.appendChild(ta);
  var row = document.createElement("div");
  row.className = "row";
  row.style.marginTop = "8px";
  var save = document.createElement("button");
  save.className = "primary small";
  save.textContent = "保存";
  save.onclick = function() { saveDraft(false); };
  var revert = document.createElement("button");
  revert.className = "ghost small";
  revert.textContent = "放弃修改";
  revert.onclick = function() { viewer.draft = null; renderPreview(); };
  row.appendChild(save);
  row.appendChild(revert);
  body.appendChild(row);
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
  var row = document.createElement("div");
  row.className = "row";
  row.style.marginTop = "10px";
  var btn = document.createElement("button");
  btn.className = "primary small";
  btn.textContent = "保存";
  btn.onclick = function() { saveDraft(false); };
  row.appendChild(btn);
  body.appendChild(row);
}

function colNames(name, cols) {
  var known = KNOWN_COLS[name];
  if (known && known[cols]) return known[cols];
  var names = [];
  for (var i = 0; i < cols; i++) names.push("c" + i);
  return names;
}

function setDraft(text) {
  if (!viewer) return;
  viewer.draft = text;
  if (viewer.draft === viewer.text) viewer.draft = null;
  el("pvDirty").classList.toggle("hidden", !viewerDirty());
}

function joinEol(lines) {
  var base = viewer.text || "";
  var eol = base.indexOf("\r\n") >= 0 ? "\r\n" : "\n";
  return lines.join(eol);
}

async function saveDraft(force) {
  if (!viewer || viewer.draft == null) return;
  var content = viewer.draft;
  var saveRel = viewer.rel;
  var bytes = new Blob([content]).size;
  if (bytes > 200 * 1024) {
    el("pvMsg").textContent = " 内容超过 200 KB 编辑上限";
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
      el("pvConflict").classList.remove("hidden");
      return;
    }
    if (resp.ok) {
      if (viewer && viewer.rel === saveRel) {
        viewer.text = content;
        viewer.draft = null;
        viewer.baseMtime = data.mtime || viewer.baseMtime;
        viewer.editing = false;
        el("pvEditBtn").textContent = "编辑";
        renderPreview();
      }
      toast("已保存 " + viewer.name);
      loadScan(curPath, true);
    } else {
      el("pvMsg").textContent = " " + (data.error || "保存失败");
    }
  } catch (e) {
    el("pvMsg").textContent = " 网络错误: " + e;
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

function themeColors() {
  var grid = "#e8e8ed", axis = "#8a8a8f";
  try {
    var cs = getComputedStyle(document.body);
    var g = cs.getPropertyValue("--grid").trim();
    var a = cs.getPropertyValue("--axis").trim();
    if (g) grid = g;
    if (a) axis = a;
  } catch (e) {}
  return {grid: grid, axis: axis};
}
function axisTicks(lo, hi, isLog) {
  var ticks = [];
  if (isLog) {
    var k0 = Math.ceil(lo - 1e-9), k1 = Math.floor(hi + 1e-9);
    if (k1 - k0 >= 0 && k1 - k0 <= 7) {
      for (var k = k0; k <= k1; k++) ticks.push({pos: k, label: "1e" + k});
      return ticks;
    }
  }
  var n = 4;
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
  var h = canvas.clientHeight || 300;
  if (w < 40 || h < 40) return null;
  canvas.width = Math.round(w * dpr);
  canvas.height = Math.round(h * dpr);
  var ctx = canvas.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  var flat = [];
  series.forEach(function(s) {
    var pts = [];
    for (var i = 0; i < s.xs.length; i++) {
      var xv = opts.logX && s.xs[i] > 0 ? Math.log10(s.xs[i]) : (opts.logX ? NaN : s.xs[i]);
      var yv = opts.logY && s.ys[i] > 0 ? Math.log10(s.ys[i]) : (opts.logY ? NaN : s.ys[i]);
      pts.push([xv, yv]);
      if (isFinite(xv) && isFinite(yv)) flat.push([xv, yv]);
    }
    s._pts = pts;
  });
  if (!flat.length) {
    ctx.fillStyle = theme.axis;
    ctx.font = "12px " + "ui-monospace, monospace";
    ctx.fillText("无数据", 14, 26);
    return null;
  }
  var xsAll = flat.map(function(p) { return p[0]; });
  var ysAll = flat.map(function(p) { return p[1]; });
  var xr = arrMinMax(xsAll), yr = arrMinMax(ysAll);
  var xmin = xr[0], xmax = xr[1], ymin = yr[0], ymax = yr[1];
  if (xmax === xmin) xmax = xmin + 1;
  if (ymax === ymin) ymax = ymin + Math.abs(ymin) * 0.1 + 1;
  var pad = {l: 58, r: 14, t: 12, b: 42};
  var pw = w - pad.l - pad.r, ph = h - pad.t - pad.b;
  function X(x) { return pad.l + pw * (x - xmin) / (xmax - xmin); }
  function Y(y) { return pad.t + ph * (1 - (y - ymin) / (ymax - ymin)); }
  ctx.strokeStyle = theme.grid;
  ctx.lineWidth = 1;
  ctx.font = "10px ui-monospace, monospace";
  ctx.fillStyle = theme.axis;
  axisTicks(ymin, ymax, opts.logY).forEach(function(t) {
    var yy = Y(t.pos);
    ctx.beginPath();
    ctx.moveTo(pad.l, yy);
    ctx.lineTo(w - pad.r, yy);
    ctx.stroke();
    ctx.fillText(t.label, 6, yy + 3);
  });
  axisTicks(xmin, xmax, opts.logX).forEach(function(t) {
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
    ctx.lineWidth = 1.8;
    ctx.beginPath();
    var started = false;
    var lastValid = null;
    for (var i = 0; i < s._pts.length; i++) {
      var xv = s._pts[i][0], yv = s._pts[i][1];
      if (!isFinite(xv) || !isFinite(yv)) {
        if (started) { ctx.stroke(); ctx.beginPath(); started = false; }
        continue;
      }
      var px = X(xv), py = Y(yv);
      if (!started) { ctx.moveTo(px, py); started = true; }
      else ctx.lineTo(px, py);
      lastValid = [px, py];
    }
    if (started) ctx.stroke();
    if (lastValid && s._pts.length < 4) {
      ctx.fillStyle = s.color || "rgb(24, 103, 174)";
      ctx.beginPath();
      ctx.arc(lastValid[0], lastValid[1], 2.6, 0, 2 * Math.PI);
      ctx.fill();
    }
    if (lastValid) lastPoints.push({px: lastValid[0], py: lastValid[1], color: s.color || "rgb(24, 103, 174)"});
  });
  canvas._xs = series.length ? series[0].xs : [];
  return lastPoints;
}

function attachHover(canvas, readoutEl, labelsFn) {
  canvas.onmousemove = function(ev) {
    var plot = {l: 58, r: 14};
    var xs = canvas._xs;
    if (!xs || !xs.length || !readoutEl) return;
    var rect = canvas.getBoundingClientRect();
    var px = ev.clientX - rect.left;
    var frac = (px - plot.l) / (canvas.clientWidth - plot.l - plot.r);
    if (frac < 0) frac = 0;
    if (frac > 1) frac = 1;
    var idx = Math.round(frac * (xs.length - 1));
    readoutEl.textContent = labelsFn(idx);
  };
  canvas.onmouseleave = function() {
    if (readoutEl) readoutEl.textContent = "";
  };
}

function setDrawerTab(tab) {
  drawerTab = tab;
  document.querySelectorAll("#drawerTabs .segbtn").forEach(function(b) {
    b.classList.toggle("on", b.getAttribute("data-tab") === tab);
  });
  el("dOutput").classList.toggle("hidden", tab !== "output");
  el("dHistory").classList.toggle("hidden", tab !== "history");
  el("dTerminal").classList.toggle("hidden", tab !== "terminal");
}
function openDrawer(tab) {
  if (tab) setDrawerTab(tab);
  el("drawer").classList.remove("hidden");
  el("drawerMask").classList.remove("hidden");
}
function closeDrawer() {
  el("drawer").classList.add("hidden");
  el("drawerMask").classList.add("hidden");
}

function showOutput(text, truncated) {
  el("output").textContent = text;
  el("outTruncNote").classList.toggle("hidden", !truncated);
  openDrawer("output");
}

async function runAction(rec, btn) {
  if (running) return;
  var runDir = curPath;
  running = true;
  btn.disabled = true;
  btn.textContent = "运行中...";
  el("runningChip").classList.remove("hidden");
  showOutput("$ " + rec.command + "\n\n目录: " + (runDir ? runDir : "（根目录）") + "\n（运行中，请稍候...）", false);
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
      record.status = data && data.error && data.error.indexOf("timed out") >= 0 ? "超时" : "失败";
      record.output = (data && data.error) || "命令未能启动。";
      showOutput("$ " + rec.command + "\n\n" + record.output, false);
    } else {
      record.status = data.returncode === 0 ? "成功" : "退出码 " + data.returncode;
      record.rc = data.returncode;
      record.output = data.output || "";
      record.truncated = !!data.truncated;
      record.image = data.image || null;
      var text = "$ " + data.command + "\n\n目录: " + (runDir ? runDir : "（根目录）") + "\n\n" + (data.output || "");
      if (data.returncode !== 0) text += "\n[退出码 " + data.returncode + "]";
      if (data.truncated) text += "\n[输出已截断至 64 KB]";
      showOutput(text, data.truncated);
      if (curPath === runDir) {
        if (data.returncode === 0) toast(data.image ? "已生成 " + data.image : "完成");
      } else {
        toast("原目录已完成: " + (data.image || rec.command));
      }
    }
  } catch (e) {
    record.status = "网络中断";
    showOutput("$ " + rec.command + "\n\n网络错误: " + e + "\n命令可能仍在服务器上运行，状态未知。", false);
  }
  record.dur = (Date.now() - t0) / 1000;
  cmdHistory.unshift(record);
  if (cmdHistory.length > 40) cmdHistory.pop();
  renderHistory();
  btn.disabled = false;
  btn.textContent = "生成图像";
  running = false;
  loadScan(curPath, true);
}

function renderHistory() {
  var list = el("histList");
  list.innerHTML = "";
  if (!cmdHistory.length) {
    list.innerHTML = '<div class="muted">本次会话尚未运行命令。</div>';
    return;
  }
  cmdHistory.forEach(function(h) {
    var d = document.createElement("div");
    d.className = "hist-item";
    var badge = h.status === "成功"
      ? '<span class="badge">' + h.status + "</span>"
      : '<span class="badge warn">' + esc(h.status) + "</span>";
    var meta = (h.dir ? "/" + esc(h.dir) : "根目录") + " · " + h.start.toLocaleTimeString() + " · " + fmtDur(h.dur);
    d.innerHTML = badge + '<span class="hcmd">' + esc(h.cmd) + "</span>" +
      '<span class="hmeta">' + meta + "</span>";
    var acts = document.createElement("span");
    acts.className = "hacts";
    var outBtn = document.createElement("button");
    outBtn.className = "ghost small";
    outBtn.textContent = "输出";
    outBtn.onclick = function() { showOutput("$ " + h.cmd + "\n\n" + h.output + (h.rc ? "\n[退出码 " + h.rc + "]" : ""), h.truncated); };
    acts.appendChild(outBtn);
    var copyBtn = document.createElement("button");
    copyBtn.className = "ghost small";
    copyBtn.textContent = "复制";
    copyBtn.onclick = function() { copyText(h.cmd); };
    acts.appendChild(copyBtn);
    if (h.image) {
      var resBtn = document.createElement("button");
      resBtn.className = "ghost small";
      resBtn.textContent = "结果";
      resBtn.title = "查看 " + h.image + "（来自 " + (h.dir || "根目录") + "）";
      resBtn.onclick = function() {
        var f = {name: h.image};
        lbList = [f];
        lbIndex = 0;
        el("lbImg").src = "/api/image?path=" + encodeURIComponent(joinRel(h.dir, h.image));
        el("lbTitle").textContent = h.image;
        el("lbCount").textContent = h.dir || "根目录";
        el("lightbox").classList.remove("hidden");
        el("lbClose").focus();
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
    termAppend(ok0 ? "-> " + (rootDir || "") + "\n" : "无法前往\n");
    return;
  }
  if (cmd.indexOf("cd ") === 0 && cmd.indexOf("&&") < 0 && cmd.indexOf(";") < 0) {
    var arg = cmd.substring(3).trim().replace(/^["']|["']$/g, "");
    var dest;
    if (arg === "..") dest = parentRel(curPath);
    else if (arg === "/") dest = "";
    else if (arg.charAt(0) === "/") { termAppend(" 错误: 请使用相对路径或直接 cd 返回根目录。\n"); return; }
    else dest = joinRel(curPath, arg);
    var ok = await loadScan(dest);
    termAppend(ok ? "-> " + dest + "\n" : "没有这个目录: " + arg + "\n");
    return;
  }
  if (termBusy) { termAppend("（上一条命令仍在运行）\n"); return; }
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
      if (data.truncated) termAppend("[输出已截断至 64 KB]\n");
      if (data.returncode !== 0) termAppend("[退出码 " + data.returncode + "]\n");
      loadScan(curPath, true);
    } else {
      termAppend(" 错误: " + (data.error || "命令失败") + "\n");
    }
  } catch (e) {
    termAppend(" 网络错误: " + e + "\n 命令可能仍在运行，状态未知。\n");
  } finally {
    termBusy = false;
    input.disabled = false;
    input.focus();
  }
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

function pollTick() {
  if (document.hidden) { schedule(); return; }
  loadScan(curPath, true).then(schedule, schedule);
}
function schedule() {
  if (pollTimer) clearTimeout(pollTimer);
  pollTimer = null;
  var sel = parseInt(el("refreshSel").value, 10);
  if (sel <= 0) return;
  var delay = pollActive ? sel * 1000 : 60000;
  pollTimer = setTimeout(pollTick, delay);
}

function applyTheme() {
  var dark = themeManual === "dark" ||
    (themeManual !== "light" && window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches);
  document.body.classList.toggle("dark", dark);
  var btn = el("themeBtn");
  if (btn) btn.innerHTML = dark ? ICON_SUN : ICON_MOON;
  if (lossView.data && el("lossCanvas")) sizeAndDrawLoss();
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
  el("loginBtn").onclick = doLogin;
  el("pwInput").addEventListener("keydown", function(ev) {
    if (ev.key === "Enter") doLogin();
  });
  el("retryBtn").onclick = function() { loadScan(curPath); };
  el("fileFilter").addEventListener("input", function() {
    renderSidebar({subdirs: lastSubdirs || [], files: lastFiles || []}, lastScanServerNow || 0, {}, {});
  });
  el("sidebarBtn").onclick = function() {
    if (window.innerWidth <= 900) {
      el("sidebarPanel").classList.toggle("open");
      return;
    }
    sidebarOpen = !sidebarOpen;
    el("layout").classList.toggle("nosidebar", !sidebarOpen);
  };
  el("themeBtn").onclick = toggleTheme;
  el("refreshSel").onchange = schedule;
  el("moreBtn").onclick = function(ev) {
    ev.stopPropagation();
    el("moreMenu").classList.toggle("hidden");
  };
  document.addEventListener("click", function(ev) {
    if (!el("moreMenu").classList.contains("hidden") && !el("moreMenu").contains(ev.target) && ev.target !== el("moreBtn")) {
      el("moreMenu").classList.add("hidden");
    }
  });
  el("menuCommands").onclick = function() {
    el("moreMenu").classList.add("hidden");
    openDrawer("output");
  };
  el("menuRefresh").onclick = function() {
    el("moreMenu").classList.add("hidden");
    loadScan(curPath);
  };
  el("runningChip").onclick = function() { openDrawer("output"); };
  document.querySelectorAll("#drawerTabs .segbtn").forEach(function(b) {
    b.onclick = function() { setDrawerTab(b.getAttribute("data-tab")); };
  });
  el("drawerClose").onclick = closeDrawer;
  el("drawerMask").onclick = closeDrawer;
  el("lbClose").onclick = closeLightbox;
  el("lightbox").onclick = function(ev) {
    if (ev.target === el("lightbox")) closeLightbox();
  };
  el("pvClose").onclick = closePreview;
  el("previewLayer").onclick = function(ev) {
    if (ev.target === el("previewLayer")) closePreview();
  };
  el("pvKeepBtn").onclick = function() {
    el("pvConfirm").classList.add("hidden");
    pendingOpen = null;
  };
  el("pvDiscardBtn").onclick = function() {
    el("pvConfirm").classList.add("hidden");
    var target = pendingOpen;
    pendingOpen = null;
    if (viewer) viewer.draft = null;
    if (target) loadPreview(target);
  };
  el("pvRereadBtn").onclick = function() {
    el("pvConflict").classList.add("hidden");
    var f = {name: viewer.name, size: null};
    if (viewer) viewer.draft = null;
    loadPreview(f);
  };
  el("pvKeepDraftBtn").onclick = function() { el("pvConflict").classList.add("hidden"); };
  el("pvOverwriteBtn").onclick = function() { el("pvConflict").classList.add("hidden"); saveDraft(true); };
  el("pvEditBtn").onclick = function() {
    if (!viewer || viewer.truncated) return;
    viewer.editing = !viewer.editing;
    el("pvEditBtn").textContent = viewer.editing ? "完成编辑" : "编辑";
    if (!viewer.editing) viewer.mode = "preview";
    renderPreview();
  };
  el("zoomClose").onclick = function() { el("chartZoom").classList.add("hidden"); };
  el("chartZoom").onclick = function(ev) {
    if (ev.target === el("chartZoom")) el("chartZoom").classList.add("hidden");
  };
  document.addEventListener("keydown", function(ev) {
    if (ev.key !== "Escape") return;
    var t = ev.target;
    if (t && t.closest && t.closest("input, textarea, select")) return;
    if (!el("chartZoom").classList.contains("hidden")) { el("chartZoom").classList.add("hidden"); return; }
    if (!el("lightbox").classList.contains("hidden")) { closeLightbox(); return; }
    if (!el("previewLayer").classList.contains("hidden")) { closePreview(); return; }
    if (!el("drawer").classList.contains("hidden")) { closeDrawer(); return; }
    if (el("sidebarPanel").classList.contains("open")) { el("sidebarPanel").classList.remove("open"); return; }
    if (!el("moreMenu").classList.contains("hidden")) el("moreMenu").classList.add("hidden");
  });
  document.addEventListener("keydown", function(ev) {
    if (el("lightbox").classList.contains("hidden")) return;
    if (ev.key === "ArrowLeft") lightboxStep(-1);
    if (ev.key === "ArrowRight") lightboxStep(1);
  });
  el("termInput").addEventListener("keydown", function(ev) {
    if (ev.key === "Enter") {
      var value = el("termInput").value;
      el("termInput").value = "";
      termExec(value);
    }
  });
  el("termOut").addEventListener("click", function() {
    var sel = window.getSelection ? String(window.getSelection()) : "";
    if (sel) return;
    if (!termBusy) el("termInput").focus();
  });
  var mq = window.matchMedia ? window.matchMedia("(prefers-color-scheme: dark)") : null;
  if (mq && mq.addEventListener) {
    mq.addEventListener("change", function() {
      if (!themeManual) applyTheme();
    });
  }
  if (window.ResizeObserver) {
    var ro = new ResizeObserver(function() {
      if (lossView.data && el("lossCanvas")) sizeAndDrawLoss();
    });
    ro.observe(el("contentCol") || document.body);
  }
  ageTimer = setInterval(updateConnText, 5000);
  applyTheme();
  loadScan("");
  schedule();
});
