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
  if (sec < 5) return T("justNow");
  if (sec < 60) return Math.floor(sec) + T("sAgo");
  if (sec < 3600) return Math.floor(sec / 60) + T("mAgo");
  if (sec < 86400) return Math.floor(sec / 3600) + T("hAgo");
  return Math.floor(sec / 86400) + T("dAgo");
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
      toast(T("copied"));
    }, function() {
      toast(fallbackCopy(text) ? T("copied") : T("copyFail"));
    });
  } else {
    toast(fallbackCopy(text) ? T("copied") : T("copyFail"));
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

/* ---- i18n ---- */
var LANG = "en";
var STRINGS = {
  en: {
    loading: "Loading working directory…",
    login: "Sign in", loginSub: "This server requires a password.", loginBtn: "Sign in",
    retry: "Retry", close: "Close", edit: "Edit", doneEdit: "Done",
    unsaved: "unsaved", dirtyWarn: "This file has unsaved changes.",
    keepEditing: "Keep editing", discardOpen: "Discard and open",
    conflictWarn: "The file changed on disk while you were editing.",
    reload: "Reload", keepDraft: "Keep my draft", overwrite: "Overwrite",
    lossCurve: "Loss curve", commands: "Commands & output",
    output: "Output", history: "History", terminal: "Terminal",
    truncated: "output truncated at 64 KB on the server",
    termNote: "Each command runs independently in the current directory; session state like export is not kept; 120 s limit per command",
    termPlaceholder: "command + Enter (cd / clear / pwd built in)",
    lbHint: "← → navigate · Esc close",
    copied: "Copied", copyFail: "Copy failed",
    justNow: "just now", sAgo: "s ago", mAgo: "m ago", hAgo: "h ago", dAgo: "d ago",
    files: "files", file: "file", dirs: "dirs", dir: "dir",
    training: "training", finished: "finished", updated: "updated",
    lossTitle: "Loss curve", figTitle: "Figures", actTitle: "Actions",
    filesTitle: "Files", subdirsTitle: "Subdirectories", emptyTitle: "Empty",
    genFigure: "Figure", viewResult: "View", copy: "Copy", run: "Run", running: "Running…",
    showMore: "Show more", filter: "Filter…",
    emptyDir: "This directory is empty.", goUp: "Go up",
    noRecords: "no valid records", unreadable: "unreadable", emptyFile: "empty",
    readingLoss: "Reading loss.out…", openFile: "open file",
    records: "records", latestGen: "latest gen", target: "target",
    rate: "rate", eta: "ETA", estimating: "estimating", lastWrite: "last write",
    status: "status", multiRun: "multi-run log", latestSeg: "latest seg",
    progress: "progress", done: "done", noTarget: "no target in nep.in",
    noNepIn: "no nep.in found", notGenerated: "not generated",
    resultAt: "result", stale: "stale", showNMore: "Show N more",
    trainDone: "Training finished (target reached)",
    trainActive: "Training active · written",
    trainIdle: "No recent write",
    trainNoNep: "loss.out found · no nep.in",
    trainNone: "No training records",
    noCmdYet: "No commands run in this session.",
    saved: "Saved", netErr: "network error", saveFail: "save failed",
    overLimit: "content exceeds the 200 KB edit limit",
    preview: "Table", textMode: "Text", formMode: "Form",
    plotOpts: "Plot options", showFirst: "showing first 100 of",
    rowsLabel: "rows", truncatedNote: "Bounded preview: first 128 KB and last 64 KB shown; editing disabled",
    errDirRead: "Directory read failed. Content below may be stale.",
    connBad: "connection problem · updated",
    updatedAt: "updated", itemsLabel: "items",
    langBtn: "中文",
  },
  zh: {
    loading: "正在读取工作目录…",
    login: "登录", loginSub: "此服务器需要密码。", loginBtn: "登录",
    retry: "重试", close: "关闭", edit: "编辑", doneEdit: "完成",
    unsaved: "未保存", dirtyWarn: "当前文件有未保存的修改。",
    keepEditing: "继续编辑", discardOpen: "放弃并打开",
    conflictWarn: "文件在磁盘上已被外部修改。",
    reload: "重新读取", keepDraft: "保留我的修改", overwrite: "覆盖",
    lossCurve: "训练曲线", commands: "命令与输出",
    output: "输出", history: "记录", terminal: "命令",
    truncated: "输出在服务端截断至 64 KB",
    termNote: "每条命令都在当前目录独立执行；export 等会话状态不会保留；单条限时 120 秒",
    termPlaceholder: "输入命令后回车（内建 cd / clear / pwd）",
    lbHint: "← → 切换 · Esc 关闭",
    copied: "已复制", copyFail: "复制失败",
    justNow: "刚刚", sAgo: " 秒前", mAgo: " 分钟前", hAgo: " 小时前", dAgo: " 天前",
    files: "个文件", file: "个文件", dirs: "个子目录", dir: "个子目录",
    training: "训练", finished: "已完成", updated: "更新于",
    lossTitle: "训练曲线", figTitle: "图片", actTitle: "可用操作",
    filesTitle: "文件", subdirsTitle: "子目录", emptyTitle: "空目录",
    genFigure: "生成图像", viewResult: "查看", copy: "复制", run: "运行", running: "运行中…",
    showMore: "显示更多", filter: "筛选…",
    emptyDir: "此目录为空。", goUp: "返回上级",
    noRecords: "无有效记录", unreadable: "无法读取", emptyFile: "文件为空",
    readingLoss: "正在读取 loss.out…", openFile: "打开文件",
    records: "记录", latestGen: "最新代数", target: "目标",
    rate: "速率", eta: "预计剩余", estimating: "估算中", lastWrite: "最后写入",
    status: "状态", multiRun: "多轮训练记录", latestSeg: "最新段",
    progress: "进度", done: "已完成", noTarget: "nep.in 未指定目标",
    noNepIn: "未发现 nep.in", notGenerated: "未生成",
    resultAt: "结果", stale: "已过期", showNMore: "再显示 N 张",
    trainDone: "训练已完成（达到目标代数）",
    trainActive: "训练进行中 · 有写入",
    trainIdle: "近期无写入",
    trainNoNep: "发现 loss.out · 无 nep.in",
    trainNone: "未发现训练记录",
    noCmdYet: "本次会话尚未运行命令。",
    saved: "已保存", netErr: "网络错误", saveFail: "保存失败",
    overLimit: "内容超过 200 KB 编辑上限",
    preview: "表格", textMode: "文本", formMode: "表单",
    plotOpts: "绘图选项", showFirst: "显示前 100 行（共",
    rowsLabel: "行）", truncatedNote: "有界预览：仅显示开头 128 KB 与结尾 64 KB；编辑已禁用",
    errDirRead: "目录读取失败，下方内容可能是旧数据。",
    connBad: "连接异常 · 更新于",
    updatedAt: "更新于", itemsLabel: "项",
    langBtn: "EN",
  },
};
function T(key) {
  var dict = STRINGS[LANG] || STRINGS.en;
  return dict[key] != null ? dict[key] : (STRINGS.en[key] != null ? STRINGS.en[key] : key);
}
function applyLang() {
  document.querySelectorAll("[data-i18n]").forEach(function(node) {
    node.textContent = T(node.getAttribute("data-i18n"));
  });
  document.querySelectorAll("[data-i18n-opt]").forEach(function(node) {
    node.textContent = T(node.getAttribute("data-i18n-opt"));
  });
  var ph = el("termInput");
  if (ph) ph.placeholder = T("termPlaceholder");
  var fl = el("fileFilter");
  if (fl) fl.placeholder = T("filter");
  el("langBtn").textContent = T("langBtn");
  document.documentElement.lang = LANG;
  renderAll();
}
function setLang(lang) {
  LANG = lang;
  try { localStorage.setItem("gk_lang", lang); } catch (e) {}
  applyLang();
}

var ICON_SUN = '<svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round"><circle cx="12" cy="12" r="4.5"/><line x1="12" y1="2" x2="12" y2="4.5"/><line x1="12" y1="19.5" x2="12" y2="22"/><line x1="2" y1="12" x2="4.5" y2="12"/><line x1="19.5" y1="12" x2="22" y2="12"/><line x1="4.6" y1="4.6" x2="6.4" y2="6.4"/><line x1="17.6" y1="17.6" x2="19.4" y2="19.4"/><line x1="4.6" y1="19.4" x2="6.4" y2="17.6"/><line x1="17.6" y1="6.4" x2="19.4" y2="4.6"/></svg>';
var ICON_MOON = '<svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"><path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"/></svg>';
var ICON_FOLDER = '<svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor"><path d="M10 4H4a2 2 0 0 0-2 2v12a2 2 0 0 0 2 2h16a2 2 0 0 0 2-2V8a2 2 0 0 0-2-2h-8l-2-2z"/></svg>';
var ICON_CHART = '<svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polyline points="22 12 18 12 15 21 9 3 6 12 2 12"/></svg>';
var ICON_IMG = '<svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><rect x="3" y="3" width="18" height="18" rx="2"/><circle cx="8.5" cy="8.5" r="1.5"/><polyline points="21 15 16 10 5 21"/></svg>';
var ICON_GEAR = '<svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><circle cx="12" cy="12" r="3"/><path d="M19.4 15a1.65 1.65 0 0 0 .33 1.82l.06.06a2 2 0 1 1-2.83 2.83l-.06-.06a1.65 1.65 0 0 0-1.82-.33 1.65 1.65 0 0 0-1 1.51V21a2 2 0 1 1-4 0v-.09A1.65 1.65 0 0 0 9 19.4a1.65 1.65 0 0 0-1.82.33l-.06.06a2 2 0 1 1-2.83-2.83l.06-.06a1.65 1.65 0 0 0 .33-1.82 1.65 1.65 0 0 0-1.51-1H3a2 2 0 1 1 0-4h.09A1.65 1.65 0 0 0 4.6 9a1.65 1.65 0 0 0-.33-1.82l-.06-.06a2 2 0 1 1 2.83-2.83l.06.06a1.65 1.65 0 0 0 1.82.33H9a1.65 1.65 0 0 0 1-1.51V3a2 2 0 1 1 4 0v.09a1.65 1.65 0 0 0 1 1.51 1.65 1.65 0 0 0 1.82-.33l.06-.06a2 2 0 1 1 2.83 2.83l-.06.06a1.65 1.65 0 0 0-.33 1.82V9a1.65 1.65 0 0 0 1.51 1H21a2 2 0 1 1 0 4h-.09a1.65 1.65 0 0 0-1.51 1z"/></svg>';
var ICON_DOC = '<svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"/><polyline points="14 2 14 8 20 8"/></svg>';
var ICON_CHEV = '<svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polyline points="6 9 12 15 18 9"/></svg>';

var PREVIEW_MAX_CLIENT = 512 * 1024;
var LOSS_COLORS = ["#1867ae", "#c18852", "#3fa06a", "#a64747", "#7b5ea7", "#8c564b"];
var KNOWN_COLS = {
  "msd.out": {4: ["t", "msd_x", "msd_y", "msd_z"], 7: ["t", "msd_x", "msd_y", "msd_z", "sdc_x", "sdc_y", "sdc_z"]},
  "sdc.out": {4: ["t", "vac_x", "vac_y", "vac_z"]},
  "loss.out": {6: ["gen", "Loss", "E_tr", "F_tr", "V_tr"], 10: ["gen", "L_total", "L1", "L2", "E_tr", "F_tr", "V_tr", "E_te", "F_te", "V_te"]},
};
var REC_NAMES = {
  plt_msd: "MSD curve", plt_sdc: "SDC curve", plt_msd_sdc: "MSD & SDC",
  plt_vac: "VAC curve", plt_thermo: "Thermo curves",
  plt_train: "Training loss & parity", plt_train_density: "Parity density",
  plt_train_test: "Train/test compare", plt_prediction: "Prediction parity",
  plt_msd_conv: "MSD convergence", plt_sigma: "Arrhenius sigma", plt_D: "Arrhenius D",
};
var FIG_PAGE = 10;

var curPath = "";
var rootDir = "";
var running = false;
var termBusy = false;
var busyFlag = false;
var scanSeq = 0;
var pollTimer = null;
var toastTimer = null;
var themeManual = null;
var connOk = null;
var lastGoodAt = null;
var firstScanDone = false;
var lastScanServerNow = null;
var lastFiles = null;
var lastSubdirs = null;
var lastPath = null;
var lastRecs = null;
var lastTraining = null;
var lastPollTs = 0;
var lastNepDone = null;
var nepSamples = [];
var cmdHistory = [];
var drawerTab = "output";
var pendingOpen = null;
var figShown = FIG_PAGE;
var lbList = [];
var lbIndex = -1;
var lbReturnFocus = null;
var viewer = null;
var lossView = { data: null, pending: false, reason: null, logX: true, logY: true, hidden: {}, fetchedAt: 0 };
var plotState = null;

try {
  themeManual = localStorage.getItem("gk_theme");
  var savedLang = localStorage.getItem("gk_lang");
  if (savedLang === "zh" || savedLang === "en") LANG = savedLang;
  var qp = new URLSearchParams(location.search).get("theme");
  if (qp === "dark" || qp === "light") themeManual = qp;
  var ql = new URLSearchParams(location.search).get("lang");
  if (ql === "zh" || ql === "en") LANG = ql;
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
  if (connOk === null) { el("connText").textContent = "…"; return; }
  if (!connOk) {
    el("connText").textContent = T("connBad") + " " + (lastGoodAt ? fmtAgo((Date.now() - lastGoodAt) / 1000) : "—");
    return;
  }
  el("connText").textContent = T("updatedAt") + " " + fmtAgo((Date.now() - lastGoodAt) / 1000);
}

async function api(url, opts) {
  var resp;
  try {
    resp = await fetch(url, opts);
  } catch (e) {
    setConn(false);
    throw e;
  }
  if (resp.status === 401) showLogin(T("loginSub"));
  return resp;
}

function showLogin(msg) {
  el("loginView").classList.remove("hidden");
  el("appView").classList.add("hidden");
  el("loading-overlay").classList.add("hidden");
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
      el("loginMsg").textContent = data.error || T("saveFail");
    }
  } catch (e) {
    el("loginMsg").textContent = T("netErr") + ": " + e;
  }
  el("loginBtn").disabled = false;
}

async function loadScan(path, silent) {
  var seq = ++scanSeq;
  var resp;
  try {
    resp = await api("/api/scan?path=" + encodeURIComponent(path));
  } catch (e) {
    if (!silent) showBanner(T("netErr"));
    return false;
  }
  if (seq !== scanSeq) return false;
  if (resp.status === 401) return false;
  var data = null;
  try { data = await resp.json(); } catch (e) {}
  if (seq !== scanSeq) return false;
  if (!resp.ok || !data) {
    setConn(false);
    showBanner((data && data.error) || T("errDirRead"));
    return false;
  }
  setConn(true);
  if (!firstScanDone) {
    firstScanDone = true;
    el("loading-overlay").classList.add("hidden");
  }
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
  _prevF = prevF; _prevS = prevS; _nowS = nowS;
  renderAll();
  if (hasLossFile(data) && !lossView.data) fetchLoss(seq);
  return true;
}
var _prevF = {};
var _prevS = {};
var _nowS = 0;

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
  if (seq === scanSeq && path === curPath) renderAll();
}

function showBanner(msg) {
  el("errorText").textContent = msg;
  el("errorBanner").classList.remove("hidden");
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

function fileIcon(name, isDir) {
  if (isDir) return ICON_FOLDER;
  if (/\.(png|jpe?g)$/i.test(name)) return ICON_IMG;
  if (/\.out$/.test(name)) return ICON_CHART;
  if (/\.(in|txt|xyz|sh|py)$/.test(name)) return ICON_GEAR;
  return ICON_DOC;
}

function renderAll() {
  if (!firstScanDone) return;
  renderStatusLine();
  renderCrumb();
  renderRail();
  renderMain();
}

function renderStatusLine() {
  var files = lastFiles || [];
  var subdirs = lastSubdirs || [];
  var title = curPath ? curPath.split("/").pop() : (rootDir ? rootDir.replace(/\/+$/, "").split("/").pop() : "root");
  el("dirTitle").textContent = title;
  el("dirTitle").title = (rootDir || "") + "/" + curPath;
  el("stFiles").textContent = files.length + " " + T(files.length === 1 ? "file" : "files") + (subdirs.length ? " · " + subdirs.length + " " + T(subdirs.length === 1 ? "dir" : "dirs") : "");
  var t = lastTraining;
  var trainTxt;
  if (t && !t.loss_empty) {
    if (t.finished) trainTxt = T("training") + " " + T("finished");
    else if (t.has_target) trainTxt = T("training") + " " + Math.round(100 * t.done / t.total) + "%";
    else trainTxt = T("training");
  } else if (hasLossFile({files: files})) {
    trainTxt = T("trainNoNep");
  } else {
    trainTxt = T("trainNone");
  }
  el("stTrain").textContent = trainTxt;
  el("stUpdated").textContent = T("updated") + " " + fmtAgo(0);
}

function renderCrumb() {
  var c = el("crumb");
  c.innerHTML = "";
  c.appendChild(crumbLink("root", ""));
  if (curPath) {
    var parts = curPath.split("/");
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
}
function crumbLink(label, rel) {
  var a = document.createElement("a");
  a.href = "#";
  a.textContent = label;
  a.onclick = function(ev) { ev.preventDefault(); loadScan(rel); };
  return a;
}

function renderRail() {
  var fl = el("fileList");
  fl.innerHTML = "";
  var filterText = el("fileFilter").value.trim().toLowerCase();
  var nowS = _nowS;
  var shown = 0;
  (lastSubdirs || []).forEach(function(s) {
    if (filterText && s.name.toLowerCase().indexOf(filterText) < 0) return;
    shown++;
    var row = document.createElement("div");
    row.className = "frow dir";
    var prev = _prevS[s.name];
    var live = prev && s.newest != null && prev.newest != null && s.newest > prev.newest;
    row.innerHTML = '<span class="fic">' + ICON_FOLDER + "</span>";
    if (live) row.innerHTML += '<span class="live" title="' + T("updatedAt") + '"></span>';
    var nm = document.createElement("span");
    nm.className = "fname";
    nm.textContent = s.name + "/";
    row.appendChild(nm);
    var meta = document.createElement("span");
    meta.className = "fmeta";
    var bits = [];
    if (s.newest != null) bits.push(fmtAgo(nowS - s.newest));
    meta.textContent = bits.join(" ");
    row.appendChild(meta);
    row.onclick = function() { loadScan(joinRel(curPath, s.name)); };
    fl.appendChild(row);
  });
  (lastFiles || []).forEach(function(f) {
    if (filterText && f.name.toLowerCase().indexOf(filterText) < 0) return;
    shown++;
    var row = document.createElement("div");
    row.className = "frow";
    var prev = _prevF[f.name];
    var live = prev && f.size > prev.size && f.mtime != null && nowS - f.mtime < 180;
    row.innerHTML = '<span class="fic">' + fileIcon(f.name, false) + "</span>";
    if (live) row.innerHTML += '<span class="live"></span>';
    var a = document.createElement("a");
    a.href = "#";
    a.className = "fname";
    a.textContent = f.name;
    a.onclick = function(ev) { ev.preventDefault(); fileClicked(f, a); };
    row.appendChild(a);
    var meta = document.createElement("span");
    meta.className = "fmeta";
    var bits = [fmtSize(f.size)];
    if (f.mtime != null) bits.push(fmtAgo(nowS - f.mtime));
    meta.textContent = bits.join(" ");
    row.appendChild(meta);
    fl.appendChild(row);
  });
  if (!shown) {
    var empty = document.createElement("div");
    empty.className = "rail-empty";
    empty.textContent = filterText ? "—" : T("emptyDir");
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

function makePanel(title, metaText) {
  var panel = document.createElement("section");
  panel.className = "panel";
  var head = document.createElement("div");
  head.className = "panel-head";
  var h = document.createElement("h2");
  h.textContent = title;
  head.appendChild(h);
  if (metaText) {
    var m = document.createElement("span");
    m.className = "pmeta";
    m.textContent = metaText;
    head.appendChild(m);
  }
  var acts = document.createElement("span");
  acts.className = "pacts";
  head.appendChild(acts);
  panel.appendChild(head);
  panel._acts = acts;
  return panel;
}

function renderMain() {
  var root = el("canvasRoot");
  root.innerHTML = "";
  var files = lastFiles || [];
  var subdirs = lastSubdirs || [];
  var imgs = files.filter(function(f) { return /\.(png|jpe?g)$/i.test(f.name); });
  var hasLoss = files.some(function(f) { return f.name === "loss.out"; });
  var recs = lastRecs || [];
  var hasAny = files.length + subdirs.length;
  var nowS = _nowS || lastScanServerNow || 0;

  if (!hasAny) {
    var empty = document.createElement("div");
    empty.className = "emptybox";
    empty.textContent = T("emptyDir");
    if (curPath) {
      var up = document.createElement("br");
      empty.appendChild(up);
      var upBtn = document.createElement("button");
      upBtn.type = "button";
      upBtn.className = "btn soft";
      upBtn.textContent = T("goUp");
      upBtn.onclick = function() { loadScan(parentRel(curPath)); };
      empty.appendChild(upBtn);
    }
    root.appendChild(empty);
    return;
  }
  if (hasLoss) root.appendChild(buildLossPanel(nowS));
  if (imgs.length) root.appendChild(buildFigPanel(imgs));
  if (recs.length) root.appendChild(buildActPanel(recs, files, nowS));
  if (!hasLoss && !imgs.length && files.length) root.appendChild(buildFilesPanel(files, nowS));
  if (!hasLoss && !imgs.length && subdirs.length && !files.length) root.appendChild(buildSubdirsPanel(subdirs, nowS));
  if (hasLoss && lossView.data) sizeAndDrawLoss();
}

/* ---- panels ---- */
function buildLossPanel(nowS) {
  var panel = makePanel(T("lossTitle"), "loss.out");
  var body = document.createElement("div");
  body.className = "panel-body";
  if (lossView.pending) {
    body.innerHTML = '<div class="muted">' + esc(T("readingLoss")) + "</div>";
    panel.appendChild(body);
    return panel;
  }
  if (!lossView.data) {
    var reason = lossView.reason === "unreadable" ? T("unreadable") : (lossView.reason === "empty" ? T("emptyFile") : T("noRecords"));
    var note = document.createElement("p");
    note.className = "muted";
    note.style.margin = "0";
    note.textContent = "loss.out " + reason + " · ";
    var openBtn = document.createElement("button");
    openBtn.type = "button";
    openBtn.className = "note-link";
    openBtn.textContent = T("openFile");
    openBtn.onclick = function() {
      var f = (lastFiles || []).find(function(x) { return x.name === "loss.out"; });
      if (f) openPreview(f);
    };
    note.appendChild(openBtn);
    body.appendChild(note);
    panel.appendChild(body);
    return panel;
  }
  var scaleSeg = document.createElement("span");
  scaleSeg.className = "seg";
  ["logX", "logY"].forEach(function(axis) {
    var b = document.createElement("button");
    b.type = "button";
    b.className = "segbtn" + (lossView[axis] ? " on" : "");
    b.setAttribute("data-axis", axis);
    b.textContent = (axis === "logX" ? "X " : "Y ") + (lossView[axis] ? "Log" : "Lin");
    b.onclick = function() {
      lossView[axis] = !lossView[axis];
      b.classList.toggle("on", lossView[axis]);
      b.textContent = (axis === "logX" ? "X " : "Y ") + (lossView[axis] ? "Log" : "Lin");
      sizeAndDrawLoss();
      if (!el("chartZoom").classList.contains("hidden")) drawZoomChart();
    };
    scaleSeg.appendChild(b);
  });
  panel._acts.appendChild(scaleSeg);
  var zoomBtn = document.createElement("button");
  zoomBtn.type = "button";
  zoomBtn.className = "btn ghost small";
  zoomBtn.textContent = "⤢";
  zoomBtn.title = "zoom";
  zoomBtn.onclick = openChartZoom;
  panel._acts.appendChild(zoomBtn);

  var layout = document.createElement("div");
  layout.className = "loss-layout";
  layout.id = "lossLayout";
  var chartCol = document.createElement("div");
  var canvas = document.createElement("canvas");
  canvas.className = "chart";
  canvas.id = "lossCanvas";
  canvas.style.height = "270px";
  chartCol.appendChild(canvas);
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
  layout.appendChild(chartCol);
  var sumBox = document.createElement("div");
  sumBox.className = "loss-summary";
  sumBox.id = "lossSummary";
  layout.appendChild(sumBox);
  body.appendChild(layout);
  panel.appendChild(body);
  return panel;
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

function sizeAndDrawLoss() {
  var canvas = el("lossCanvas");
  if (!canvas || !lossView.data) return;
  var w = canvas.clientWidth;
  if (w < 60) return;
  drawLossTo(canvas);
  renderLossSummary();
  renderLossLegend("lossLegend");
  renderLossNotes();
}

function drawLossTo(canvas) {
  var series = visibleLossSeries();
  var d = lossView.data;
  drawSeries(canvas, series, {
    logX: lossView.logX, logY: lossView.logY,
    xLabel: d && d.x_mode === "generation" ? "generation" : "record #",
    yLabel: "Loss functions"
  });
  attachHover(canvas, canvas.id === "zoomCanvas" ? null : el("lossHover"), function(idx) {
    var parts = ["x=" + fmtNum(d.xs[idx])];
    series.forEach(function(s) {
      parts.push(s.label + "=" + fmtNum(s.ys[idx]));
    });
    return "# " + (idx + 1) + " · " + parts.join(" · ");
  });
}

function renderLossSummary() {
  var box = el("lossSummary");
  if (!box || !lossView.data) return;
  var d = lossView.data;
  var t = lastTraining;
  box.innerHTML = "";
  function line(label, value) {
    var l = document.createElement("div");
    l.className = "ls-line";
    l.innerHTML = "<span>" + esc(label) + "</span>";
    var v = document.createElement("strong");
    v.textContent = value;
    l.appendChild(v);
    box.appendChild(l);
  }
  line(T("records"), d.count.toLocaleString());
  line(T("latestGen"), d.last_gen != null ? d.last_gen.toLocaleString() : "--");
  if (d.multi_run) {
    line(T("status"), T("multiRun"));
    line(T("latestSeg"), "≥ " + d.seg_start_gen);
  }
  if (t && t.has_target) {
    line(T("target"), t.total.toLocaleString());
    if (!t.multi_run) {
      var track = document.createElement("div");
      track.className = "bar-track";
      var fill = document.createElement("div");
      fill.className = "bar-fill" + (t.finished ? " done" : "");
      fill.style.width = Math.min(100, 100 * t.done / t.total) + "%";
      track.appendChild(fill);
      box.appendChild(track);
      var rate = trainingRate();
      if (t.finished) {
        line(T("status"), T("done"));
      } else {
        line(T("rate"), fmtGenRate(rate.rate));
        line(T("eta"), rate.estimating ? T("estimating") : (rate.rate ? fmtDur((t.total - t.done) / rate.rate) : "--"));
        line(T("lastWrite"), t.loss_mtime != null ? fmtAgo(lastScanServerNow - t.loss_mtime) : "--");
      }
    }
  } else {
    line(T("target"), t ? T("noTarget") : T("noNepIn"));
    line(T("lastWrite"), d.mtime != null ? fmtAgo(lastScanServerNow - d.mtime) : "--");
  }
}

function renderLossLegend(legendId) {
  var legend = el(legendId);
  if (!legend || !lossView.data) return;
  legend.innerHTML = "";
  (lossView.data.series || []).forEach(function(s, i) {
    var item = document.createElement("button");
    item.type = "button";
    item.className = "lg-item" + (lossView.hidden[s.label] ? " off" : "");
    var dot = document.createElement("span");
    dot.className = "lg-dot";
    dot.style.background = LOSS_COLORS[i % LOSS_COLORS.length];
    item.appendChild(dot);
    item.appendChild(document.createTextNode(s.label));
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
  if (d.sampled) notes.appendChild(noteLine(d.points + " / " + d.count + " pts"));
  if (d.skipped > 0) notes.appendChild(noteLine(d.skipped + " skipped"));
  if (d.multi_run) notes.appendChild(noteLine("multi-run log"));
  if (lossView.logY) {
    var nonpos = 0;
    (d.series || []).forEach(function(s) {
      if (lossView.hidden[s.label]) return;
      s.values.forEach(function(v) { if (v <= 0) nonpos++; });
    });
    if (nonpos > 0) {
      var line = noteLine(nonpos + " non-positive values hidden (log Y)");
      var sw = document.createElement("button");
      sw.type = "button";
      sw.className = "note-link";
      sw.textContent = "Y → Linear";
      sw.onclick = function() {
        lossView.logY = false;
        sizeAndDrawLoss();
        var b = document.querySelector('#lossLayout .seg [data-axis="logY"]');
        if (b) { b.classList.remove("on"); b.textContent = "Y Lin"; }
      };
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
  var w = Math.min(window.innerWidth * 0.9, 1120);
  var h = Math.min(window.innerHeight * 0.7, 600);
  canvas.style.width = w + "px";
  canvas.style.height = h + "px";
  drawLossTo(canvas);
  renderLossLegend("zoomLegend");
  el("zoomScale").textContent = "X " + (lossView.logX ? "Log" : "Linear") + " · Y " + (lossView.logY ? "Log" : "Linear");
}

function buildFigPanel(imgs) {
  var panel = makePanel(T("figTitle"), imgs.length + "");
  var body = document.createElement("div");
  body.className = "panel-body";
  var strip = document.createElement("div");
  strip.className = "fig-strip";
  imgs.slice(0, figShown).forEach(function(f) {
    var cell = document.createElement("button");
    cell.type = "button";
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
    strip.appendChild(cell);
  });
  if (imgs.length > figShown) {
    var more = document.createElement("button");
    more.type = "button";
    more.className = "fig-more";
    more.textContent = "+ " + (imgs.length - figShown);
    more.onclick = function() {
      figShown += FIG_PAGE;
      renderMain();
    };
    strip.appendChild(more);
  }
  body.appendChild(strip);
  panel.appendChild(body);
  return panel;
}

function buildActPanel(recs, files, nowS) {
  var panel = makePanel(T("actTitle"), recs.length + "");
  var body = document.createElement("div");
  body.className = "panel-body";
  var fileMap = {};
  files.forEach(function(f) { fileMap[f.name] = f; });
  recs.forEach(function(rec) {
    var row = document.createElement("div");
    row.className = "act-line";
    var nm = document.createElement("span");
    nm.className = "aname";
    nm.textContent = REC_NAMES[rec.action] || rec.action;
    nm.title = rec.command;
    row.appendChild(nm);
    var inp = document.createElement("span");
    inp.className = "ain";
    inp.textContent = rec.evidence.join(", ");
    inp.title = rec.command;
    row.appendChild(inp);
    var res = document.createElement("span");
    res.className = "ares";
    var outFile = fileMap[rec.produces];
    if (outFile && outFile.mtime != null) {
      res.textContent = T("resultAt") + " " + fmtAgo(nowS - outFile.mtime) + (rec.stale ? " · " + T("stale") : "");
      if (rec.stale) res.classList.add("stale");
      res.title = rec.produces;
    } else {
      res.textContent = T("notGenerated");
    }
    row.appendChild(res);
    var acts = document.createElement("span");
    acts.className = "acts";
    var run = document.createElement("button");
    run.type = "button";
    run.className = "btn primary small";
    run.textContent = T("genFigure");
    run.disabled = busyFlag || running;
    run.title = busyFlag ? "busy" : rec.command;
    run.onclick = function() { runAction(rec, run); };
    acts.appendChild(run);
    if (outFile) {
      var open = document.createElement("button");
      open.type = "button";
      open.className = "btn ghost small";
      open.textContent = T("viewResult");
      open.onclick = function() { openLightbox(outFile, open); };
      acts.appendChild(open);
    }
    var copyBtn = document.createElement("button");
    copyBtn.type = "button";
    copyBtn.className = "btn ghost small";
    copyBtn.textContent = T("copy");
    copyBtn.onclick = function() { copyText(rec.command); };
    acts.appendChild(copyBtn);
    row.appendChild(acts);
    body.appendChild(row);
  });
  panel.appendChild(body);
  return panel;
}

function buildFilesPanel(files, nowS) {
  var panel = makePanel(T("filesTitle"), files.length + "");
  var body = document.createElement("div");
  body.className = "panel-body";
  var table = document.createElement("table");
  table.className = "filetable";
  files.slice(0, 14).forEach(function(f) {
    var tr = document.createElement("tr");
    var td1 = document.createElement("td");
    td1.innerHTML = '<span class="fic" style="display:inline-flex;color:var(--ink-faint)">' + fileIcon(f.name, false) + "</span> ";
    var a = document.createElement("span");
    a.className = "fname";
    a.textContent = f.name;
    a.style.cursor = "pointer";
    a.onclick = function() { openPreview(f); };
    td1.appendChild(a);
    var td2 = document.createElement("td");
    td2.className = "fmeta";
    td2.textContent = [fmtSize(f.size), f.mtime != null ? fmtAgo(nowS - f.mtime) : ""].filter(Boolean).join(" · ");
    tr.appendChild(td1);
    tr.appendChild(td2);
    table.appendChild(tr);
  });
  body.appendChild(table);
  if (files.length > 14) {
    var more = document.createElement("p");
    more.className = "muted";
    more.style.margin = "8px 0 0";
    more.textContent = "+" + (files.length - 14);
    body.appendChild(more);
  }
  panel.appendChild(body);
  return panel;
}

function buildSubdirsPanel(subdirs, nowS) {
  var panel = makePanel(T("subdirsTitle"), subdirs.length + "");
  var body = document.createElement("div");
  body.className = "panel-body";
  subdirs.forEach(function(s) {
    var row = document.createElement("div");
    row.className = "subline";
    row.innerHTML = '<span class="fic">' + ICON_FOLDER + "</span>";
    var nm = document.createElement("span");
    nm.className = "sname";
    nm.textContent = s.name + "/";
    row.appendChild(nm);
    var meta = document.createElement("span");
    meta.className = "smeta";
    var bits = [s.count + " " + T("itemsLabel")];
    if (s.newest != null) bits.push(fmtAgo(nowS - s.newest));
    meta.textContent = bits.join(" · ");
    row.appendChild(meta);
    row.onclick = function() { loadScan(joinRel(curPath, s.name)); };
    body.appendChild(row);
  });
  panel.appendChild(body);
  return panel;
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

/* ---- preview ---- */
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
  el("pvBody").innerHTML = '<div class="muted">…</div>';
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
        if (!j) { previewError(rel, "too large"); return; }
        applyPartial(rel, j);
      });
    }
    previewError(rel, "preview failed");
  }).catch(function() {
    if (viewer !== null && viewer.rel === rel) previewError(rel, T("netErr"));
  });
}
function applyPartial(rel, j) {
  if (viewer === null || viewer.rel !== rel) return;
  viewer.truncated = true;
  viewer.baseMtime = j.mtime;
  viewer.text = j.head + "\n[... " + (j.size - j.head_bytes - j.tail_bytes).toLocaleString() + " bytes omitted ...]\n" + j.tail;
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
  el("pvBody").innerHTML = "";
}
function closePreview() {
  if (viewerDirty()) {
    el("pvConfirm").classList.remove("hidden");
    return;
  }
  el("previewLayer").classList.add("hidden");
  viewer = null;
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
  if (numeric) modes.push(["preview", T("preview")]);
  modes.push(["text", T("textMode")]);
  if (isForm) modes.push(["form", T("formMode")]);
  el("pvEditBtn").classList.toggle("hidden", viewer.truncated);
  var modeBar = el("pvModes");
  modeBar.innerHTML = "";
  modeBar.classList.toggle("hidden", modes.length <= 1);
  modes.forEach(function(m) {
    var b = document.createElement("button");
    b.type = "button";
    b.className = "segbtn";
    b.textContent = m[1];
    b.setAttribute("data-mode", m[0]);
    b.onclick = function() { viewer.mode = m[0]; viewer.editing = false; el("pvEditBtn").textContent = T("edit"); renderPreview(); };
    modeBar.appendChild(b);
  });
  if (viewer.mode === "form" && !isForm) viewer.mode = "preview";
  modeBar.querySelectorAll(".segbtn").forEach(function(b) {
    b.classList.toggle("on", b.getAttribute("data-mode") === viewer.mode && !viewer.editing);
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
    note.textContent = T("showFirst") + " " + numeric.rows.length + " " + T("rowsLabel");
    body.appendChild(note);
  }
  var opts = document.createElement("div");
  opts.className = "plotopts";
  var head = document.createElement("div");
  head.className = "plotopts-head";
  head.innerHTML = "<span class='achev'>" + ICON_CHEV + "</span> " + T("plotOpts");
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
  var state = plotState && plotState.rel === viewer.rel ? plotState : {rel: viewer.rel, x: 0, ys: [1], logX: false, logY: false};
  plotState = state;
  var ctl = document.createElement("div");
  ctl.className = "chartctl";
  var xSel = document.createElement("select");
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
  canvas.style.height = "280px";
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
      return "# " + (idx + 1) + " · " + parts.join(" · ");
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
    note.textContent = T("truncatedNote");
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
  save.type = "button";
  save.className = "btn primary small";
  save.textContent = T("saved") === "已保存" ? "保存" : "Save";
  save.onclick = function() { saveDraft(false); };
  var revert = document.createElement("button");
  revert.type = "button";
  revert.className = "btn ghost small";
  revert.textContent = LANG === "zh" ? "放弃修改" : "Revert";
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
  row.style.marginTop = "8px";
  var btn = document.createElement("button");
  btn.type = "button";
  btn.className = "btn primary small";
  btn.textContent = LANG === "zh" ? "保存" : "Save";
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
    el("pvMsg").textContent = T("overLimit");
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
        el("pvEditBtn").textContent = T("edit");
        renderPreview();
      }
      toast(T("saved") + " " + viewer.name);
      loadScan(curPath, true);
    } else {
      el("pvMsg").textContent = " " + (data.error || T("saveFail"));
    }
  } catch (e) {
    el("pvMsg").textContent = " " + T("netErr") + ": " + e;
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

/* ---- chart engine ---- */
function themeColors() {
  var grid = "#e8ebef", axis = "#8994a0";
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
  var h = canvas.clientHeight || 280;
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
    ctx.font = "11px ui-monospace, monospace";
    ctx.fillText("no data", 14, 24);
    return null;
  }
  var xr = arrMinMax(flat.map(function(p) { return p[0]; }));
  var yr = arrMinMax(flat.map(function(p) { return p[1]; }));
  var xmin = xr[0], xmax = xr[1], ymin = yr[0], ymax = yr[1];
  if (xmax === xmin) xmax = xmin + 1;
  if (ymax === ymin) ymax = ymin + Math.abs(ymin) * 0.1 + 1;
  var pad = {l: 54, r: 12, t: 10, b: 38};
  var pw = w - pad.l - pad.r, ph = h - pad.t - pad.b;
  function X(x) { return pad.l + pw * (x - xmin) / (xmax - xmin); }
  function Y(y) { return pad.t + ph * (1 - (y - ymin) / (ymax - ymin)); }
  ctx.strokeStyle = theme.grid;
  ctx.lineWidth = 1;
  ctx.font = "9.5px ui-monospace, monospace";
  ctx.fillStyle = theme.axis;
  axisTicks(ymin, ymax, opts.logY).forEach(function(t) {
    var yy = Y(t.pos);
    ctx.beginPath();
    ctx.moveTo(pad.l, yy);
    ctx.lineTo(w - pad.r, yy);
    ctx.stroke();
    ctx.fillText(t.label, 5, yy + 3);
  });
  axisTicks(xmin, xmax, opts.logX).forEach(function(t) {
    var xx = X(t.pos);
    ctx.beginPath();
    ctx.moveTo(xx, pad.t);
    ctx.lineTo(xx, h - pad.b);
    ctx.stroke();
    ctx.fillText(t.label, xx - 13, h - pad.b + 15);
  });
  if (opts.xLabel) ctx.fillText(opts.xLabel, pad.l + pw / 2 - ctx.measureText(opts.xLabel).width / 2, h - 6);
  if (opts.yLabel) {
    ctx.save();
    ctx.translate(11, pad.t + ph / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText(opts.yLabel, -ctx.measureText(opts.yLabel).width / 2, 0);
    ctx.restore();
  }
  var lastPoints = [];
  series.forEach(function(s) {
    if (!s._pts.length) return;
    ctx.strokeStyle = s.color || "#1867ae";
    ctx.lineWidth = 1.7;
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
      ctx.fillStyle = s.color || "#1867ae";
      ctx.beginPath();
      ctx.arc(lastValid[0], lastValid[1], 2.4, 0, 2 * Math.PI);
      ctx.fill();
    }
    if (lastValid) lastPoints.push({px: lastValid[0], py: lastValid[1], color: s.color || "#1867ae"});
  });
  canvas._xs = series.length ? series[0].xs : [];
  return lastPoints;
}
function attachHover(canvas, readoutEl, labelsFn) {
  canvas.onmousemove = function(ev) {
    var xs = canvas._xs;
    if (!xs || !xs.length) return;
    var rect = canvas.getBoundingClientRect();
    var px = ev.clientX - rect.left;
    var frac = (px - 54) / (canvas.clientWidth - 66);
    if (frac < 0) frac = 0;
    if (frac > 1) frac = 1;
    var idx = Math.round(frac * (xs.length - 1));
    if (readoutEl) readoutEl.textContent = labelsFn(idx);
  };
  canvas.onmouseleave = function() {
    if (readoutEl) readoutEl.textContent = "";
  };
}

/* ---- drawer & exec ---- */
function setDrawerTab(tab) {
  drawerTab = tab;
  document.querySelectorAll(".drawer-head .segbtn[data-tab]").forEach(function(b) {
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
  btn.textContent = T("running");
  showOutput("$ " + rec.command + "\n\ndir: " + (runDir || "(root)") + "\n…", false);
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
      record.status = "failed";
      record.output = (data && data.error) || "failed to start";
      showOutput("$ " + rec.command + "\n\n" + record.output, false);
    } else {
      record.status = data.returncode === 0 ? "ok" : "exit " + data.returncode;
      record.rc = data.returncode;
      record.output = data.output || "";
      record.truncated = !!data.truncated;
      record.image = data.image || null;
      var text = "$ " + data.command + "\n\ndir: " + (runDir || "(root)") + "\n\n" + (data.output || "");
      if (data.returncode !== 0) text += "\n[exit " + data.returncode + "]";
      if (data.truncated) text += "\n[truncated]";
      showOutput(text, data.truncated);
      if (data.returncode === 0) toast(data.image ? data.image : "done");
    }
  } catch (e) {
    record.status = "network";
    showOutput("$ " + rec.command + "\n\n" + T("netErr") + ": " + e, false);
  }
  record.dur = (Date.now() - t0) / 1000;
  cmdHistory.unshift(record);
  if (cmdHistory.length > 40) cmdHistory.pop();
  renderHistory();
  btn.disabled = false;
  btn.textContent = T("genFigure");
  running = false;
  loadScan(runDir, true);
}
function renderHistory() {
  var list = el("histList");
  list.innerHTML = "";
  if (!cmdHistory.length) {
    list.innerHTML = '<div class="muted">' + esc(T("noCmdYet")) + "</div>";
    return;
  }
  cmdHistory.forEach(function(h) {
    var d = document.createElement("div");
    d.className = "hist-item";
    var badge = h.status === "ok" ? '<span class="badge">' + h.status + "</span>" : '<span class="badge dirty">' + esc(h.status) + "</span>";
    var meta = (h.dir ? "/" + esc(h.dir) : "root") + " · " + h.start.toLocaleTimeString() + " · " + fmtDur(h.dur);
    d.innerHTML = badge + '<span class="hcmd">' + esc(h.cmd) + "</span>" + '<span class="hmeta">' + meta + "</span>";
    var acts = document.createElement("span");
    acts.className = "hacts";
    var outBtn = document.createElement("button");
    outBtn.type = "button";
    outBtn.className = "btn ghost small";
    outBtn.textContent = T("output");
    outBtn.onclick = function() { showOutput("$ " + h.cmd + "\n\n" + h.output + (h.rc ? "\n[exit " + h.rc + "]" : ""), h.truncated); };
    acts.appendChild(outBtn);
    if (h.image) {
      var resBtn = document.createElement("button");
      resBtn.type = "button";
      resBtn.className = "btn ghost small";
      resBtn.textContent = T("viewResult");
      resBtn.onclick = function() {
        el("lbImg").src = "/api/image?path=" + encodeURIComponent(joinRel(h.dir, h.image));
        el("lbTitle").textContent = h.image;
        el("lbCount").textContent = h.dir || "root";
        el("lightbox").classList.remove("hidden");
        el("lbClose").focus();
      };
      acts.appendChild(resBtn);
    }
    var copyBtn = document.createElement("button");
    copyBtn.type = "button";
    copyBtn.className = "btn ghost small";
    copyBtn.textContent = T("copy");
    copyBtn.onclick = function() { copyText(h.cmd); };
    acts.appendChild(copyBtn);
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
    termAppend(ok0 ? "-> " + (rootDir || "") + "\n" : "cannot go\n");
    return;
  }
  if (cmd.indexOf("cd ") === 0 && cmd.indexOf("&&") < 0 && cmd.indexOf(";") < 0) {
    var arg = cmd.substring(3).trim().replace(/^["']|["']$/g, "");
    var dest;
    if (arg === "..") dest = parentRel(curPath);
    else if (arg === "/") dest = "";
    else if (arg.charAt(0) === "/") { termAppend(" use relative paths or plain cd\n"); return; }
    else dest = joinRel(curPath, arg);
    var ok = await loadScan(dest);
    termAppend(ok ? "-> " + dest + "\n" : "no such directory: " + arg + "\n");
    return;
  }
  if (termBusy) { termAppend("(busy)\n"); return; }
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
      if (data.truncated) termAppend("[truncated]\n");
      if (data.returncode !== 0) termAppend("[exit " + data.returncode + "]\n");
      loadScan(curPath, true);
    } else {
      termAppend(" error: " + (data.error || "failed") + "\n");
    }
  } catch (e) {
    termAppend(" network error: " + e + "\n");
  } finally {
    termBusy = false;
    input.disabled = false;
    input.focus();
  }
}

/* ---- polling: manual by default ---- */
function pollTick() {
  if (document.hidden) { schedulePoll(); return; }
  loadScan(curPath, true).then(schedulePoll, schedulePoll);
}
function schedulePoll() {
  if (pollTimer) clearTimeout(pollTimer);
  pollTimer = null;
  var secs = parseInt(el("refreshSel").value, 10);
  if (secs <= 0) return;
  pollTimer = setTimeout(pollTick, secs * 1000);
}

/* ---- theme ---- */
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
  el("langBtn").onclick = function() { setLang(LANG === "en" ? "zh" : "en"); };
  el("loginBtn").onclick = doLogin;
  el("pwInput").addEventListener("keydown", function(ev) {
    if (ev.key === "Enter") doLogin();
  });
  el("retryBtn").onclick = function() { loadScan(curPath); };
  el("refreshBtn").onclick = function() { loadScan(curPath); };
  el("refreshSel").onchange = schedulePoll;
  el("fileFilter").addEventListener("input", renderRail);
  el("themeBtn").onclick = toggleTheme;
  el("drawerBtn").onclick = function() { openDrawer("output"); };
  document.querySelectorAll(".drawer-head .segbtn[data-tab]").forEach(function(b) {
    b.onclick = function() { setDrawerTab(b.getAttribute("data-tab")); };
  });
  el("drawerClose").onclick = closeDrawer;
  el("drawerMask").onclick = closeDrawer;
  el("zoomClose").onclick = function() { el("chartZoom").classList.add("hidden"); };
  el("chartZoom").onclick = function(ev) {
    if (ev.target === el("chartZoom")) el("chartZoom").classList.add("hidden");
  };
  el("lbClose").onclick = closeLightbox;
  el("lightbox").onclick = function(ev) {
    if (ev.target === el("lightbox")) closeLightbox();
  };
  el("pvClose").onclick = closePreview;
  el("previewLayer").onclick = function(ev) {
    if (ev.target === el("previewLayer")) closePreview();
  };
  el("pvKeepBtn").onclick = function() { el("pvConfirm").classList.add("hidden"); pendingOpen = null; };
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
    el("pvEditBtn").textContent = viewer.editing ? T("doneEdit") : T("edit");
    renderPreview();
  };
  document.addEventListener("keydown", function(ev) {
    if (ev.key !== "Escape") return;
    var t = ev.target;
    if (t && t.closest && t.closest("input, textarea, select")) return;
    if (!el("chartZoom").classList.contains("hidden")) { el("chartZoom").classList.add("hidden"); return; }
    if (!el("lightbox").classList.contains("hidden")) { closeLightbox(); return; }
    if (!el("previewLayer").classList.contains("hidden")) { closePreview(); return; }
    if (!el("drawer").classList.contains("hidden")) { closeDrawer(); return; }
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
    ro.observe(document.body);
  }
  applyLang();
  applyTheme();
  loadScan("");
});
