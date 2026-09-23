"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     server.py
Category:   Utilities
Purpose:    Start a local web server for the current working directory.
            It scans the directory, recommends allowlisted GPUMDkit plot
            commands based on detected output files, runs them on click,
            and renders the generated figures in the browser. It also
            monitors growing output files, reports NEP training progress,
            previews data files as tables or charts in a modal viewer,
            edits small keyword files (run.in / nep.in) through forms,
            shows figure thumbnails with a lightbox, and provides a web
            terminal for arbitrary quick commands (login required).
            Supports optional password login. The web UI assets live
            in the web/ folder next to this script and are reloaded
            automatically when they change on disk.
Usage:      gpumdkit.sh -server [port] [-b <address>] [-pw <password>]
            python3 server.py [port] [-b <address>] [-pw <password>]
Arguments:
  port              Port to listen on (default: 8888)
  -b <address>      Bind address (default: 127.0.0.1; 0.0.0.0 for LAN)
  -pw <password>    Password for web login. Required when binding to a
                    non-loopback address; prompted interactively if
                    omitted there. Optional on loopback.
Output:
  (no files; serves the web UI and a JSON API over HTTP)
Author:     Zihan YAN (yanzihan@westlake.edu.cn)
Last-modified: 2026-09-22
=============================================================================
"""

import array
import getpass
import glob
import hashlib
import hmac
import ipaddress
import json
import os
import secrets
import socket
import stat
import subprocess
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

DEFAULT_PORT = 8888
DEFAULT_BIND = "127.0.0.1"
SESSION_TTL = 12 * 3600
PBKDF2_ROUNDS = 200000
POST_MAX_BYTES = 512 * 1024
PREVIEW_MAX_BYTES = 512 * 1024
PARTIAL_HEAD_BYTES = 128 * 1024
PARTIAL_TAIL_BYTES = 64 * 1024
EDIT_MAX_BYTES = 200 * 1024
EXEC_TIMEOUT = 120
EXEC_MAX_CHARS = 4096
EXEC_OUTPUT_CHARS = 65536
LOSS_POINT_CAP = 600
LOSS_SERIES_MAX = 6
IMAGE_MAX_BYTES = 20 * 1024 * 1024
LISTING_LIMIT = 2000
SUBDIR_LIMIT = 100
SUBDIR_STAT_LIMIT = 100
COOKIE_NAME = "gpumdkit_token"
IMAGE_TYPES = {".png": "image/png", ".jpg": "image/jpeg", ".jpeg": "image/jpeg"}

RUN_ACTIONS = {
    "plt_msd": {
        "argv": ["-plt", "msd", "save"],
        "requires": ["msd.out"],
        "produces": "msd.png",
        "timeout": 900,
    },
    "plt_sdc": {
        "argv": ["-plt", "sdc", "save"],
        "requires": ["msd.out"],
        "produces": "sdc.png",
        "timeout": 900,
    },
    "plt_msd_sdc": {
        "argv": ["-plt", "msd_sdc", "save"],
        "requires": ["msd.out"],
        "produces": "msd_sdc.png",
        "timeout": 900,
    },
    "plt_vac": {
        "argv": ["-plt", "vac", "save"],
        "requires": ["sdc.out"],
        "produces": "vac.png",
        "timeout": 900,
    },
    "plt_thermo": {
        "argv": ["-plt", "thermo", "save"],
        "requires": ["thermo.out"],
        "produces": "thermo.png",
        "timeout": 900,
    },
    "plt_train": {
        "argv": ["-plt", "train", "save"],
        "requires": ["loss.out", "energy_train.out", "force_train.out"],
        "produces": "train.png",
        "timeout": 900,
    },
    "plt_train_density": {
        "argv": ["-plt", "train_density", "save"],
        "requires": ["loss.out", "energy_train.out", "force_train.out"],
        "produces": "train_density.png",
        "timeout": 900,
    },
    "plt_train_test": {
        "argv": ["-plt", "train_test", "save"],
        "requires": [
            "energy_train.out",
            "force_train.out",
            "stress_train.out",
            "energy_test.out",
            "force_test.out",
            "stress_test.out",
        ],
        "produces": "train_test.png",
        "timeout": 900,
    },
    "plt_prediction": {
        "argv": ["-plt", "prediction", "save"],
        "requires": ["energy_train.out", "force_train.out"],
        "produces": "prediction.png",
        "timeout": 900,
    },
    "plt_msd_conv": {
        "argv": ["-plt", "msd_conv", "save"],
        "requires_glob": ["msd_step*.out"],
        "produces": "msd_convergence.png",
        "timeout": 900,
    },
    "plt_sigma": {
        "argv": ["-plt", "sigma", "save"],
        "batch": {"suffix": "K", "requires": ["msd.out", "thermo.out"], "min": 2},
        "produces": "Arrhenius_sigma.png",
        "timeout": 900,
    },
    "plt_D": {
        "argv": ["-plt", "D", "save"],
        "batch": {"suffix": "K", "requires": ["msd.out"], "min": 2},
        "produces": "Arrhenius_D.png",
        "timeout": 900,
    },
}

DESCRIPTIONS = {
    "plt_msd": "Plot MSD (x, y, z) with slope annotations",
    "plt_sdc": "Plot SDC (x, y, z) from msd.out",
    "plt_msd_sdc": "Plot MSD and SDC together",
    "plt_vac": "Plot VAC from sdc.out",
    "plt_thermo": "Plot thermodynamic properties from thermo.out",
    "plt_train": "Plot NEP training loss and parity panels",
    "plt_train_density": "Plot training parity density panels",
    "plt_train_test": "Compare train vs test energy/force/stress",
    "plt_prediction": "Plot prediction parity for train/test",
    "plt_msd_conv": "MSD convergence check across msd_step*.out",
    "plt_sigma": "Arrhenius plot of conductivity across temperature dirs",
    "plt_D": "Arrhenius plot of diffusivity across temperature dirs",
}

USAGE_LINES = [
    " Usage: gpumdkit.sh -server [port] [-b <address>] [-pw <password>]",
    "    or: python3 server.py [port] [-b <address>] [-pw <password>]",
    "",
    " Arguments:",
    "   port             Port to listen on (default: 8888)",
    "   -b <address>     Bind address (default: 127.0.0.1; 0.0.0.0 for LAN)",
    "   -pw <password>   Web login password; required for non-loopback binds",
    "",
    " Examples:",
    "   gpumdkit.sh -server",
    "   gpumdkit.sh -server 9000 -b 0.0.0.0 -pw mysecret",
    "",
    " Notes:",
    "   The server sandbox is the working directory where it is started.",
    "   Remote access without opening ports: ssh -L <port>:127.0.0.1:<port>",
    "   The web terminal executes real shell commands after login.",
    "",
]

ROOT = os.path.realpath(os.getcwd())
KIT_ROOT = os.environ.get("GPUMDkit_path") or os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir, os.pardir)
)
GPUMDKIT_SH = os.path.join(KIT_ROOT, "gpumdkit.sh")
LOGO_SQUARE = os.path.join(KIT_ROOT, "docs", "Gallery", "gpumdkit_logo.png")
LOGO_LATERAL = os.path.join(KIT_ROOT, "docs", "Gallery", "gpumdkit_logo_lateral.png")

SALT = None
PASSWORD_HASH = None
SESSIONS = {}
SESSIONS_LOCK = threading.Lock()
AUTH_LOCK = threading.Lock()
LOGIN_FAILURES = 0
LOCKOUT_UNTIL = 0.0
EXEC_LOCK = threading.Lock()
LOSS_CACHE = {}
LOSS_LOCK = threading.Lock()


def parse_args(argv):
    """Parse command-line arguments; return (port, bind, password)."""
    port = DEFAULT_PORT
    bind = DEFAULT_BIND
    password = None
    i = 0
    n = len(argv)
    while i < n:
        arg = argv[i]
        if arg in ("-h", "--help"):
            for line in USAGE_LINES:
                print(line)
            sys.exit(0)
        if arg == "-b":
            if i + 1 >= n:
                print(" Error: option '-b' requires an address argument.")
                sys.exit(1)
            bind = argv[i + 1]
            i += 2
            continue
        if arg == "-pw":
            if i + 1 >= n:
                print(" Error: option '-pw' requires a password argument.")
                sys.exit(1)
            password = argv[i + 1]
            i += 2
            continue
        if arg.startswith("-"):
            print(f" Error: unknown option '{arg}'. See -h for usage.")
            sys.exit(1)
        try:
            port = int(arg)
        except ValueError:
            print(f" Error: '{arg}' is not a valid port number.")
            sys.exit(1)
        if not 1 <= port <= 65535:
            print(f" Error: port {port} is out of range 1-65535.")
            sys.exit(1)
        i += 1
    if not bind:
        print(" Error: the bind address must not be empty.")
        sys.exit(1)
    if password == "":
        password = None
    return port, bind, password


def is_loopback(address):
    """Return True if the address only accepts local connections."""
    if address == "localhost":
        return True
    try:
        return ipaddress.ip_address(address).is_loopback
    except ValueError:
        return False


def primary_lan_ip():
    """Return the primary non-loopback IPv4 address of this host, or None."""
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    try:
        sock.connect(("10.255.255.255", 1))
        return sock.getsockname()[0]
    except OSError:
        return None
    finally:
        sock.close()


def resolve_in_root(rel):
    """Resolve a relative path inside the sandbox root; None on escape."""
    rel = (rel or "").strip()
    if "\x00" in rel:
        return None
    rel = rel.replace("\\", "/")
    while rel.startswith("/"):
        rel = rel[1:]
    candidate = os.path.join(ROOT, rel) if rel else ROOT
    target = os.path.realpath(candidate)
    if target == ROOT or target.startswith(ROOT + os.sep):
        return target
    return None


def count_columns(path):
    """Return the number of columns in the first data line of a file."""
    try:
        with open(path, "r", errors="replace") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                return len(line.split())
    except OSError:
        pass
    return 0


def is_temp_dir(name):
    """Return True for temperature-style directory names such as 500K."""
    return name.endswith("K") and name[:-1].isdigit()


def _loss_parse_append(path, entry):
    """Parse complete lines beyond entry["offset"] into the loss cache."""
    try:
        with open(path, "rb") as handle:
            handle.seek(entry["offset"])
            chunk = handle.read()
    except OSError:
        return False
    if not chunk:
        return True
    ends_newline = chunk.endswith(b"\n")
    lines = chunk.split(b"\n")
    complete = lines[:-1]
    consumed = entry["offset"]
    cols = entry["cols"]
    arrays = entry["arrays"]
    for raw in complete:
        consumed += len(raw) + 1
        text = raw.decode("utf-8", errors="replace").strip()
        if not text or text.startswith("#"):
            continue
        parts = text.split()
        if cols == 0:
            cols = len(parts)
            arrays = [array.array("d") for _ in range(cols)]
        if arrays is None or len(parts) != cols:
            entry["skipped"] += 1
            continue
        values = []
        valid = True
        for token in parts:
            try:
                value = float(token)
            except ValueError:
                valid = False
                break
            if value != value or value in (float("inf"), float("-inf")):
                valid = False
                break
            values.append(value)
        if not valid:
            entry["skipped"] += 1
            continue
        gen = values[0]
        if entry["last_gen"] is not None and gen < entry["last_gen"]:
            entry["multi_run"] = True
            entry["seg_start_gen"] = gen
        if entry["first_gen"] is None:
            entry["first_gen"] = gen
        entry["last_gen"] = gen
        for i in range(cols):
            arrays[i].append(values[i])
        entry["count"] += 1
    if arrays:
        entry["arrays"] = arrays
    entry["cols"] = cols
    if ends_newline:
        entry["offset"] = entry["offset"] + len(chunk)
    else:
        entry["offset"] = consumed
    return True


def loss_entry(path):
    """Return the cached loss record for a file, reparsing on change.

    Appends are parsed incrementally from the stored offset; a shrinking
    file (truncation or rewrite) triggers a full reparse. Rewrites that
    keep the file at least as large are treated as appends; the README
    documents this boundary.
    """
    try:
        info = os.stat(path)
        size = info.st_size
        mtime = info.st_mtime
    except OSError:
        with LOSS_LOCK:
            LOSS_CACHE.pop(path, None)
        return None
    with LOSS_LOCK:
        entry = LOSS_CACHE.get(path)
        if entry and entry["size"] == size and entry["mtime"] == mtime:
            return entry
        if entry is None or size < entry["offset"]:
            entry = {
                "size": 0,
                "mtime": 0.0,
                "offset": 0,
                "count": 0,
                "skipped": 0,
                "cols": 0,
                "arrays": None,
                "first_gen": None,
                "last_gen": None,
                "seg_start_gen": None,
                "multi_run": False,
            }
        if not _loss_parse_append(path, entry):
            return None
        entry["size"] = size
        entry["mtime"] = mtime
        LOSS_CACHE[path] = entry
        return entry


def loss_series_labels(cols, count):
    """Return loss series labels matching plt_train.py conventions."""
    if cols >= 7:
        return ["Total", "L1-Reg", "L2-Reg", "Energy-train", "Force-train", "Virial-train"][:count]
    if cols == 6:
        return ["Loss", "Energy-train", "Force-train", "Virial-train"][:count]
    return ["c" + str(i) for i in range(1, count + 1)]


def loss_payload(target):
    """Build the bounded loss chart payload for a directory."""
    path = os.path.join(target, "loss.out")
    entry = loss_entry(path)
    payload = {"path": "loss.out", "now": int(time.time())}
    if entry is None:
        payload["empty"] = True
        payload["reason"] = "unreadable"
        return payload
    payload["mtime"] = int(entry["mtime"])
    payload["size"] = entry["size"]
    payload["skipped"] = entry["skipped"]
    if entry["count"] == 0 or not entry["arrays"]:
        payload["empty"] = True
        payload["reason"] = "no valid records" if entry["cols"] else "empty"
        payload["count"] = entry["count"]
        return payload
    n = entry["count"]
    cols = entry["cols"]
    if n <= LOSS_POINT_CAP:
        idxs = list(range(n))
        payload["sampled"] = False
    else:
        step = n / LOSS_POINT_CAP
        idxs = []
        last = -1
        for i in range(LOSS_POINT_CAP):
            k = int(i * step)
            if k != last:
                idxs.append(k)
                last = k
        if idxs[-1] != n - 1:
            idxs.append(n - 1)
        payload["sampled"] = True
    payload["empty"] = False
    payload["count"] = n
    payload["points"] = len(idxs)
    payload["cols"] = cols
    payload["first_gen"] = entry["first_gen"]
    payload["last_gen"] = entry["last_gen"]
    payload["multi_run"] = entry["multi_run"]
    payload["seg_start_gen"] = entry["seg_start_gen"]
    arrays = entry["arrays"]
    if entry["multi_run"] or entry["first_gen"] is None:
        payload["x_mode"] = "index"
        payload["xs"] = [i + 1 for i in idxs]
    else:
        payload["x_mode"] = "generation"
        payload["xs"] = [arrays[0][i] for i in idxs]
    n_series = 6 if cols >= 7 else (4 if cols == 6 else cols - 1)
    n_series = min(n_series, cols - 1, LOSS_SERIES_MAX)
    labels = loss_series_labels(cols, n_series)
    series = []
    for c in range(1, n_series + 1):
        series.append(
            {"label": labels[c - 1], "values": [arrays[c][i] for i in idxs]}
        )
    payload["series"] = series
    return payload


def training_status(target):
    """Summarize NEP training state from nep.in and the loss cache."""
    nep_in = os.path.join(target, "nep.in")
    loss_file = os.path.join(target, "loss.out")
    if not os.path.isfile(nep_in) or not os.path.isfile(loss_file):
        return None
    total = None
    try:
        with open(nep_in, "r", errors="replace") as handle:
            for line in handle:
                parts = line.split()
                if parts and parts[0] == "generation":
                    try:
                        total = int(float(parts[1]))
                    except (ValueError, IndexError):
                        total = None
                    break
    except OSError:
        return None
    has_target = total is not None
    if not has_target:
        total = 100000
    entry = loss_entry(loss_file)
    if entry is None:
        return {
            "total": total,
            "has_target": has_target,
            "done": 0,
            "count": 0,
            "skipped": 0,
            "multi_run": False,
            "first_gen": None,
            "last_gen": None,
            "seg_start_gen": None,
            "loss_mtime": None,
            "loss_size": None,
            "loss_empty": True,
            "finished": False,
        }
    last_gen = entry["last_gen"]
    done = last_gen if last_gen is not None else 0
    return {
        "total": total,
        "has_target": has_target,
        "done": min(done, total),
        "count": entry["count"],
        "skipped": entry["skipped"],
        "multi_run": entry["multi_run"],
        "first_gen": entry["first_gen"],
        "last_gen": last_gen,
        "seg_start_gen": entry["seg_start_gen"],
        "loss_mtime": int(entry["mtime"]) if entry["mtime"] else None,
        "loss_size": entry["size"],
        "loss_empty": entry["count"] == 0,
        "finished": has_target and last_gen is not None and last_gen >= total,
    }


def newest_mtime(paths):
    """Return the newest existing mtime among paths, or None."""
    newest = None
    for path in paths[:100]:
        try:
            mtime = os.path.getmtime(path)
        except OSError:
            continue
        if newest is None or mtime > newest:
            newest = mtime
    return newest


def recommendation_stale(target, spec):
    """True when an existing result file is older than its newest input."""
    output = os.path.join(target, spec["produces"])
    if not os.path.isfile(output):
        return False
    try:
        output_mtime = os.path.getmtime(output)
    except OSError:
        return False
    inputs = [os.path.join(target, name) for name in spec.get("requires", [])]
    for pattern in spec.get("requires_glob", []):
        inputs.extend(glob.glob(os.path.join(target, pattern))[:100])
    batch = spec.get("batch")
    if batch:
        try:
            names = [
                name
                for name in os.listdir(target)
                if is_temp_dir(name) and os.path.isdir(os.path.join(target, name))
            ]
        except OSError:
            names = []
        for name in names[:100]:
            inputs.extend(os.path.join(target, name, req) for req in batch["requires"])
    newest = newest_mtime(inputs)
    return newest is not None and output_mtime < newest


def build_recommendations(target):
    """Detect output files and return recommended allowlisted commands."""
    recs = []

    def present(name):
        return os.path.isfile(os.path.join(target, name))

    def add(action, evidence):
        spec = RUN_ACTIONS[action]
        recs.append(
            {
                "action": action,
                "command": "gpumdkit.sh " + " ".join(spec["argv"]),
                "description": DESCRIPTIONS[action],
                "evidence": evidence,
                "produces": spec["produces"],
                "stale": recommendation_stale(target, spec),
            }
        )

    if present("msd.out"):
        columns = count_columns(os.path.join(target, "msd.out"))
        if columns >= 7:
            add("plt_msd_sdc", ["msd.out (7 cols)"])
            add("plt_sdc", ["msd.out (7 cols)"])
        add("plt_msd", ["msd.out"])
    if present("sdc.out"):
        add("plt_vac", ["sdc.out"])
    if present("thermo.out"):
        add("plt_thermo", ["thermo.out"])
    train_core = ["loss.out", "energy_train.out", "force_train.out"]
    if all(present(name) for name in train_core):
        add("plt_train", train_core)
        add("plt_train_density", train_core)
    train_test_files = [
        "energy_train.out",
        "force_train.out",
        "stress_train.out",
        "energy_test.out",
        "force_test.out",
        "stress_test.out",
    ]
    if all(present(name) for name in train_test_files):
        add("plt_train_test", train_test_files)
    if present("energy_train.out") and present("force_train.out"):
        add("plt_prediction", ["energy_train.out", "force_train.out"])
    step_files = glob.glob(os.path.join(target, "msd_step*.out"))
    if step_files:
        add("plt_msd_conv", [f"{len(step_files)} msd_step*.out"])
    try:
        entries = [
            name
            for name in os.listdir(target)
            if is_temp_dir(name) and os.path.isdir(os.path.join(target, name))
        ]
    except OSError:
        entries = []
    sigma_dirs = []
    d_dirs = []
    for name in entries:
        full = os.path.join(target, name)
        has_msd = os.path.isfile(os.path.join(full, "msd.out"))
        has_thermo = os.path.isfile(os.path.join(full, "thermo.out"))
        if has_msd and has_thermo:
            sigma_dirs.append(name)
        if has_msd:
            d_dirs.append(name)
    if len(sigma_dirs) >= 2:
        add("plt_sigma", [f"{len(sigma_dirs)} temperature dirs"])
    if len(d_dirs) >= 2:
        add("plt_D", [f"{len(d_dirs)} temperature dirs"])
    return recs


def subdir_summaries(target, dir_names):
    """Return compact per-subdirectory summaries for the sidebar view."""
    summaries = []
    for name in dir_names[:SUBDIR_LIMIT]:
        full = os.path.join(target, name)
        count = 0
        newest = None
        has_msd = has_thermo = has_loss = False
        try:
            with os.scandir(full) as entries:
                stat_count = 0
                for entry in entries:
                    count += 1
                    if stat_count >= SUBDIR_STAT_LIMIT or not entry.is_file(
                        follow_symlinks=False
                    ):
                        continue
                    stat_count += 1
                    try:
                        mtime = entry.stat(follow_symlinks=False).st_mtime
                    except OSError:
                        continue
                    if newest is None or mtime > newest:
                        newest = mtime
                    if entry.name == "msd.out":
                        has_msd = True
                    elif entry.name == "thermo.out":
                        has_thermo = True
                    elif entry.name == "loss.out":
                        has_loss = True
        except OSError:
            continue
        summaries.append(
            {
                "name": name,
                "count": count,
                "newest": int(newest) if newest is not None else None,
                "msd": has_msd,
                "thermo": has_thermo,
                "loss": has_loss,
            }
        )
    return summaries


def scan_directory(target):
    """List directory entries, recommendations, and training status."""
    try:
        names = sorted(os.listdir(target))
    except OSError:
        return None
    truncated = len(names) > LISTING_LIMIT
    dirs = []
    files = []
    for name in names[:LISTING_LIMIT]:
        full = os.path.join(target, name)
        if os.path.isdir(full):
            dirs.append({"name": name, "type": "dir"})
            continue
        try:
            info = os.stat(full)
            size = info.st_size
            mtime = int(info.st_mtime)
        except OSError:
            size = -1
            mtime = None
        files.append({"name": name, "type": "file", "size": size, "mtime": mtime})
    return {
        "dirs": dirs,
        "files": files,
        "subdirs": subdir_summaries(target, [d["name"] for d in dirs]),
        "recommendations": build_recommendations(target),
        "training": training_status(target),
        "truncated": truncated,
    }


def create_session():
    """Create a login session and return its token."""
    token = secrets.token_urlsafe(32)
    with SESSIONS_LOCK:
        SESSIONS[token] = time.time() + SESSION_TTL
    return token


def check_session(token):
    """Return True if the token belongs to a valid, unexpired session."""
    if not token:
        return False
    now = time.time()
    with SESSIONS_LOCK:
        for key in [k for k, exp in SESSIONS.items() if exp < now]:
            del SESSIONS[key]
        expiry = SESSIONS.get(token)
    return expiry is not None and now < expiry


def execute_action(action, target_dir):
    """Run an allowlisted GPUMDkit command; return (payload, HTTP status)."""
    spec = RUN_ACTIONS[action]
    env = dict(os.environ)
    env["MPLBACKEND"] = "Agg"
    argv = ["bash", GPUMDKIT_SH] + spec["argv"]
    try:
        proc = subprocess.run(
            argv,
            cwd=target_dir,
            env=env,
            capture_output=True,
            text=True,
            timeout=spec["timeout"],
        )
    except subprocess.TimeoutExpired:
        return {"error": f"command timed out after {spec['timeout']} seconds"}, 504
    except OSError as exc:
        return {"error": f"failed to start the command: {exc}"}, 500
    output = proc.stdout or ""
    if proc.stderr:
        output = output + ("\n" if output else "") + proc.stderr
    image = None
    if os.path.isfile(os.path.join(target_dir, spec["produces"])):
        image = spec["produces"]
    payload = {
        "ok": proc.returncode == 0,
        "returncode": proc.returncode,
        "command": "gpumdkit.sh " + " ".join(spec["argv"]),
        "output": output[:EXEC_OUTPUT_CHARS],
        "truncated": len(output) > EXEC_OUTPUT_CHARS,
        "image": image,
    }
    return payload, 200


def print_banner(bind, port, loopback, has_password):
    """Print the startup summary for the terminal."""
    if loopback:
        display_host = "127.0.0.1" if bind != "localhost" else "localhost"
    elif bind == "0.0.0.0":
        display_host = primary_lan_ip() or "<this-host-IP>"
    else:
        display_host = bind
    print(" GPUMDkit web server starting (Scripts/utils/server/server.py)")
    print(f" Root directory : {ROOT}")
    print(f" gpumdkit.sh    : {GPUMDKIT_SH}")
    print(f" Web UI         : http://{display_host}:{port}")
    if loopback:
        print(f" Remote access  : ssh -L {port}:127.0.0.1:{port} <user>@<this-host>")
    else:
        print(" Warning        : the server is exposed to the network; keep the")
        print("                  password private and stop the server after use.")
    if has_password:
        print(" Password       : enabled")
    else:
        print(" Password       : disabled (loopback bind, no -pw given)")
    print(" Press Ctrl+C to stop the server.")
    sys.stdout.flush()


WEB_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "web")
WEB_CACHE = {}
WEB_LOCK = threading.Lock()


def web_asset(name):
    """Return the current content of a web asset, reloading on change."""
    path = os.path.join(WEB_DIR, name)
    try:
        mtime = os.path.getmtime(path)
    except OSError:
        return None
    with WEB_LOCK:
        cached = WEB_CACHE.get(name)
        if cached is not None and cached[0] == mtime:
            return cached[1]
        try:
            with open(path, "r", encoding="utf-8") as handle:
                content = handle.read()
        except OSError:
            return None
        WEB_CACHE[name] = (mtime, content)
        return content


class KitHandler(BaseHTTPRequestHandler):
    """HTTP handler serving the GPUMDkit web UI and its JSON API."""

    protocol_version = "HTTP/1.1"
    server_version = "GPUMDkitWeb/1.0"
    timeout = 120

    def do_GET(self):
        """Route GET requests."""
        self._route("GET")

    def do_POST(self):
        """Route POST requests."""
        self._route("POST")

    def log_message(self, fmt, *log_args):
        """Write one sanitized access-log line to stdout."""
        msg = (fmt % log_args).replace("\r", " ").replace("\n", " ")
        sys.stdout.write(" %s %s\n" % (self.log_date_time_string(), msg))
        sys.stdout.flush()

    def _route(self, method):
        parsed = urlparse(self.path)
        path = parsed.path
        query = parse_qs(parsed.query)
        try:
            if method == "GET":
                if path == "/":
                    self._serve_web("index.html", "text/html; charset=utf-8")
                elif path == "/style.css":
                    self._serve_web("style.css", "text/css; charset=utf-8")
                elif path == "/app.js":
                    self._serve_web("app.js", "application/javascript; charset=utf-8")
                elif path in ("/logo.png", "/favicon.ico"):
                    self._serve_logo(LOGO_SQUARE)
                elif path == "/logo_lateral.png":
                    self._serve_logo(LOGO_LATERAL)
                elif path == "/api/scan":
                    self._api_scan(query)
                elif path == "/api/loss":
                    self._api_loss(query)
                elif path == "/api/file":
                    self._api_file(query)
                elif path == "/api/image":
                    self._api_image(query)
                else:
                    self._send_json({"error": "not found"}, 404)
            else:
                if path == "/api/login":
                    self._api_login()
                elif path == "/api/run":
                    self._api_run()
                elif path == "/api/exec":
                    self._api_exec()
                elif path == "/api/save":
                    self._api_save()
                else:
                    self._send_json({"error": "not found"}, 404)
        except (BrokenPipeError, ConnectionResetError):
            self.close_connection = True
        except Exception as exc:
            try:
                self._send_json({"error": f"internal error: {exc}"}, 500)
            except Exception:
                self.close_connection = True

    def _send_body(self, status, body, content_type, extra_headers=None):
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Cache-Control", "no-store")
        for key, value in extra_headers or []:
            self.send_header(key, value)
        self.end_headers()
        self.wfile.write(body)

    def _send_json(self, obj, status=200, extra_headers=None):
        body = json.dumps(obj).encode("utf-8")
        self._send_body(status, body, "application/json; charset=utf-8", extra_headers)

    def _session_token(self):
        cookie = self.headers.get("Cookie", "")
        for part in cookie.split(";"):
            key, _, value = part.strip().partition("=")
            if key == COOKIE_NAME:
                return value
        return self.headers.get("X-Auth-Token")

    def _authorized(self):
        if PASSWORD_HASH is None:
            return True
        return check_session(self._session_token())

    def _serve_web(self, name, content_type):
        body = web_asset(name)
        if body is None:
            self._send_json({"error": f"web asset '{name}' not found"}, 404)
            return
        self._send_body(200, body.encode("utf-8"), content_type)

    def _serve_logo(self, path):
        if not os.path.isfile(path):
            self._send_json({"error": "logo not found"}, 404)
            return
        try:
            with open(path, "rb") as handle:
                data = handle.read(IMAGE_MAX_BYTES)
        except OSError:
            self._send_json({"error": "cannot read the logo"}, 404)
            return
        self._send_body(200, data, "image/png")

    def _read_json(self):
        try:
            length = int(self.headers.get("Content-Length", "0") or 0)
        except ValueError:
            self._send_json({"error": "invalid Content-Length header"}, 400)
            return None
        if length <= 0 or length > POST_MAX_BYTES:
            self._send_json({"error": "missing or oversized request body"}, 400)
            return None
        raw = self.rfile.read(length)
        try:
            body = json.loads(raw.decode("utf-8"))
        except (ValueError, UnicodeDecodeError):
            self._send_json({"error": "invalid JSON body"}, 400)
            return None
        if not isinstance(body, dict):
            self._send_json({"error": "invalid JSON body"}, 400)
            return None
        return body

    def _api_login(self):
        global LOGIN_FAILURES, LOCKOUT_UNTIL
        body = self._read_json()
        if body is None:
            return
        if PASSWORD_HASH is None:
            self._send_json({"ok": True, "message": "no password required"})
            return
        password = body.get("password")
        if not isinstance(password, str):
            self._send_json({"error": "password is required"}, 400)
            return
        with AUTH_LOCK:
            now = time.time()
            if now < LOCKOUT_UNTIL:
                self._send_json({"error": "too many failed attempts, retry later"}, 429)
                return
            digest = hashlib.pbkdf2_hmac(
                "sha256", password.encode("utf-8"), SALT, PBKDF2_ROUNDS
            )
            if not hmac.compare_digest(digest, PASSWORD_HASH):
                LOGIN_FAILURES += 1
                if LOGIN_FAILURES >= 5:
                    LOCKOUT_UNTIL = now + 60
                    LOGIN_FAILURES = 0
                self._send_json({"error": "wrong password"}, 401)
                return
            LOGIN_FAILURES = 0
        token = create_session()
        cookie = (
            f"{COOKIE_NAME}={token}; Path=/; HttpOnly; SameSite=Strict; "
            f"Max-Age={SESSION_TTL}"
        )
        self._send_json({"ok": True}, 200, extra_headers=[("Set-Cookie", cookie)])

    def _api_scan(self, query):
        if not self._authorized():
            self._send_json({"error": "unauthorized"}, 401)
            return
        rel = (query.get("path") or [""])[0]
        target = resolve_in_root(rel)
        if target is None or not os.path.isdir(target):
            self._send_json({"error": "invalid directory"}, 400)
            return
        data = scan_directory(target)
        if data is None:
            self._send_json({"error": "cannot read the directory"}, 400)
            return
        rel_norm = os.path.relpath(target, ROOT)
        data["path"] = "" if rel_norm == "." else rel_norm
        data["root"] = ROOT
        data["now"] = int(time.time())
        data["busy"] = EXEC_LOCK.locked()
        self._send_json(data)

    def _api_loss(self, query):
        if not self._authorized():
            self._send_json({"error": "unauthorized"}, 401)
            return
        rel = (query.get("path") or [""])[0]
        target = resolve_in_root(rel)
        if target is None or not os.path.isdir(target):
            self._send_json({"error": "invalid directory"}, 400)
            return
        payload = loss_payload(target)
        rel_norm = os.path.relpath(target, ROOT)
        payload["dir"] = "" if rel_norm == "." else rel_norm
        self._send_json(payload)

    def _api_file(self, query):
        if not self._authorized():
            self._send_json({"error": "unauthorized"}, 401)
            return
        rel = (query.get("path") or [""])[0]
        target = resolve_in_root(rel)
        if target is None or not os.path.isfile(target):
            self._send_json({"error": "file not found"}, 400)
            return
        try:
            info = os.stat(target)
            size = info.st_size
            mtime = int(info.st_mtime)
        except OSError:
            self._send_json({"error": "cannot stat the file"}, 400)
            return
        if size > PREVIEW_MAX_BYTES:
            if (query.get("partial") or [""])[0] == "1":
                self._api_file_partial(target, size)
                return
            self._send_json(
                {
                    "error": "file too large to preview (limit "
                    f"{PREVIEW_MAX_BYTES // 1024} KB, size {size} bytes); "
                    "a bounded head/tail preview is available"
                },
                413,
            )
            return
        try:
            with open(target, "rb") as handle:
                data = handle.read(PREVIEW_MAX_BYTES + 1)
        except OSError:
            self._send_json({"error": "cannot read the file"}, 400)
            return
        if b"\x00" in data:
            self._send_json({"error": "binary file, no text preview"}, 400)
            return
        text = data.decode("utf-8", errors="replace")
        self._send_body(
            200,
            text.encode("utf-8"),
            "text/plain; charset=utf-8",
            extra_headers=[("X-File-Mtime", str(mtime))],
        )

    def _api_file_partial(self, target, size):
        try:
            with open(target, "rb") as handle:
                head = handle.read(PARTIAL_HEAD_BYTES)
        except OSError:
            self._send_json({"error": "cannot read the file"}, 400)
            return
        if b"\x00" in head:
            self._send_json({"error": "binary file, no text preview"}, 400)
            return
        tail = b""
        if size > PARTIAL_HEAD_BYTES:
            try:
                with open(target, "rb") as handle:
                    handle.seek(max(0, size - PARTIAL_TAIL_BYTES))
                    tail = handle.read()
            except OSError:
                tail = b""
        payload = {
            "truncated": True,
            "size": size,
            "head_bytes": len(head),
            "tail_bytes": len(tail),
            "head": head.decode("utf-8", errors="replace"),
            "tail": tail.decode("utf-8", errors="replace"),
            "mtime": int(os.stat(target).st_mtime),
        }
        self._send_json(payload)

    def _api_image(self, query):
        if not self._authorized():
            self._send_json({"error": "unauthorized"}, 401)
            return
        rel = (query.get("path") or [""])[0]
        target = resolve_in_root(rel)
        if target is None or not os.path.isfile(target):
            self._send_json({"error": "image not found"}, 400)
            return
        ext = os.path.splitext(target)[1].lower()
        if ext not in IMAGE_TYPES:
            self._send_json({"error": "unsupported image type"}, 400)
            return
        try:
            size = os.path.getsize(target)
        except OSError:
            self._send_json({"error": "cannot stat the file"}, 400)
            return
        if size > IMAGE_MAX_BYTES:
            self._send_json({"error": "image too large"}, 413)
            return
        try:
            with open(target, "rb") as handle:
                data = handle.read(IMAGE_MAX_BYTES + 1)
        except OSError:
            self._send_json({"error": "cannot read the image"}, 400)
            return
        self._send_body(200, data, IMAGE_TYPES[ext])

    def _api_run(self):
        if not self._authorized():
            self._send_json({"error": "unauthorized"}, 401)
            return
        body = self._read_json()
        if body is None:
            return
        action = body.get("action")
        rel = body.get("path", "")
        if action not in RUN_ACTIONS:
            self._send_json({"error": "unknown action"}, 400)
            return
        if not isinstance(rel, str):
            self._send_json({"error": "invalid path"}, 400)
            return
        target = resolve_in_root(rel)
        if target is None or not os.path.isdir(target):
            self._send_json({"error": "invalid directory"}, 400)
            return
        spec = RUN_ACTIONS[action]
        for required in spec.get("requires", []):
            if not os.path.isfile(os.path.join(target, required)):
                self._send_json(
                    {"error": f"required file '{required}' not found in the selected directory"},
                    400,
                )
                return
        for pattern in spec.get("requires_glob", []):
            if not glob.glob(os.path.join(target, pattern)):
                self._send_json(
                    {"error": f"no file matches '{pattern}' in the selected directory"},
                    400,
                )
                return
        batch = spec.get("batch")
        if batch:
            try:
                entries = [
                    name
                    for name in os.listdir(target)
                    if is_temp_dir(name) and os.path.isdir(os.path.join(target, name))
                ]
            except OSError:
                entries = []
            ok_dirs = []
            for name in entries:
                full = os.path.join(target, name)
                if all(
                    os.path.isfile(os.path.join(full, req)) for req in batch["requires"]
                ):
                    ok_dirs.append(name)
            if len(ok_dirs) < batch["min"]:
                self._send_json(
                    {"error": f"needs at least {batch['min']} temperature directories with " + " and ".join(batch["requires"])},
                    400,
                )
                return
        if not EXEC_LOCK.acquire(blocking=False):
            self._send_json({"error": "another command is already running"}, 409)
            return
        try:
            payload, status = execute_action(action, target)
        finally:
            EXEC_LOCK.release()
        self._send_json(payload, status)

    def _api_exec(self):
        if not self._authorized():
            self._send_json({"error": "unauthorized"}, 401)
            return
        body = self._read_json()
        if body is None:
            return
        command = body.get("command")
        rel = body.get("path", "")
        if not isinstance(command, str) or not command.strip():
            self._send_json({"error": "command is required"}, 400)
            return
        if len(command) > EXEC_MAX_CHARS:
            self._send_json(
                {"error": f"command longer than {EXEC_MAX_CHARS} characters"}, 400
            )
            return
        if not isinstance(rel, str):
            self._send_json({"error": "invalid path"}, 400)
            return
        target = resolve_in_root(rel)
        if target is None or not os.path.isdir(target):
            self._send_json({"error": "invalid directory"}, 400)
            return
        if not EXEC_LOCK.acquire(blocking=False):
            self._send_json({"error": "another command is already running"}, 409)
            return
        try:
            env = dict(os.environ)
            env["TERM"] = "dumb"
            env["PAGER"] = "cat"
            try:
                proc = subprocess.run(
                    ["bash", "-c", command],
                    cwd=target,
                    env=env,
                    capture_output=True,
                    text=True,
                    timeout=EXEC_TIMEOUT,
                )
            except subprocess.TimeoutExpired:
                self._send_json(
                    {"error": f"command timed out after {EXEC_TIMEOUT} seconds"}, 504
                )
                return
            except OSError as exc:
                self._send_json({"error": f"failed to start: {exc}"}, 500)
                return
        finally:
            EXEC_LOCK.release()
        output = proc.stdout or ""
        if proc.stderr:
            output = output + ("\n" if output else "") + proc.stderr
        self._send_json(
            {
                "ok": proc.returncode == 0,
                "returncode": proc.returncode,
                "output": output[:EXEC_OUTPUT_CHARS],
                "truncated": len(output) > EXEC_OUTPUT_CHARS,
            }
        )

    def _api_save(self):
        if not self._authorized():
            self._send_json({"error": "unauthorized"}, 401)
            return
        body = self._read_json()
        if body is None:
            return
        rel = body.get("path")
        content = body.get("content")
        if not isinstance(rel, str) or not isinstance(content, str):
            self._send_json({"error": "path and content are required"}, 400)
            return
        target = resolve_in_root(rel)
        if target is None or not os.path.isfile(target):
            self._send_json({"error": "file not found"}, 400)
            return
        base_mtime = body.get("base_mtime")
        if base_mtime is not None:
            if not isinstance(base_mtime, int):
                self._send_json({"error": "invalid base_mtime"}, 400)
                return
            try:
                current_mtime = int(os.stat(target).st_mtime)
            except OSError:
                self._send_json({"error": "cannot stat the file"}, 400)
                return
            if current_mtime != base_mtime:
                self._send_json(
                    {"error": "conflict", "mtime": current_mtime},
                    409,
                )
                return
        data = content.encode("utf-8")
        if len(data) > EDIT_MAX_BYTES:
            self._send_json(
                {"error": f"edited content exceeds the {EDIT_MAX_BYTES // 1024} KB limit"},
                413,
            )
            return
        tmp = target + ".kit_tmp_" + secrets.token_hex(6)
        try:
            with open(tmp, "wb") as handle:
                handle.write(data)
            os.chmod(tmp, stat.S_IMODE(os.stat(target).st_mode))
            os.replace(tmp, target)
        except OSError as exc:
            try:
                os.unlink(tmp)
            except OSError:
                pass
            self._send_json({"error": f"failed to save the file: {exc}"}, 500)
            return
        try:
            saved_mtime = int(os.stat(target).st_mtime)
        except OSError:
            saved_mtime = None
        self._send_json({"ok": True, "size": len(data), "mtime": saved_mtime})


def main():
    """Parse options, configure auth, and serve until Ctrl+C."""
    global SALT, PASSWORD_HASH
    if not os.path.isfile(GPUMDKIT_SH):
        print(f" Error: gpumdkit.sh not found at '{GPUMDKIT_SH}'.")
        print("        Set GPUMDkit_path or start the server via 'gpumdkit.sh -server'.")
        sys.exit(1)
    port, bind, password = parse_args(sys.argv[1:])
    loopback = is_loopback(bind)
    if not loopback and not password:
        try:
            password = getpass.getpass(" Set a password for web login: ")
        except Exception:
            password = ""
        if not password:
            print(" Error: a password is required when binding to a non-loopback address.")
            sys.exit(1)
    if password:
        SALT = secrets.token_bytes(16)
        PASSWORD_HASH = hashlib.pbkdf2_hmac(
            "sha256", password.encode("utf-8"), SALT, PBKDF2_ROUNDS
        )
    try:
        httpd = ThreadingHTTPServer((bind, port), KitHandler)
    except OSError as exc:
        print(f" Error: cannot listen on {bind}:{port} ({exc}).")
        sys.exit(1)
    httpd.daemon_threads = True
    print_banner(bind, port, loopback, bool(password))
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print(" Server stopped.")
    finally:
        httpd.server_close()


if __name__ == "__main__":
    main()
