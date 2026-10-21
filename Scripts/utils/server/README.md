# GPUMDkit Web Console (`Scripts/utils/server/`)

A zero-dependency web console for GPUMDkit: it scans a working directory,
recommends allowlisted GPUMDkit plot commands, runs them, monitors NEP
training and MD jobs, renders figures, previews data files, edits small
keyword files with conflict detection, and provides a web terminal.

This README is for maintainers and AI agents modifying or debugging the
server. Read it before changing anything.

## Quick facts

- Entry: `gpumdkit.sh -server [port] [-b <address>] [-pw <password>]`
  (routes to `Scripts/utils/server/server.py`)
- Python 3.7+ (3.9 tested), **standard library only**
- Frontend: vanilla HTML/CSS/JS, **no frameworks, no CDN, no build step**
  (must keep working offline and on HPC login nodes)
- The server sandbox is the working directory where it is started

## Start and access the Web Console

Choose the directory that should be visible in the browser, then start the
server from that directory. Activate the environment that provides GPUMDkit
and its plotting dependencies before starting it:

```bash
conda activate gpumdkit
cd /path/to/working-directory
gpumdkit.sh -server
```

The default address is `127.0.0.1:8888`. For local use, open
`http://127.0.0.1:8888` in a browser on the same machine. The process stays
in the foreground; press `Ctrl+C` in that terminal to stop it.

### Remote access through SSH (recommended)

On the remote host, start the server from the desired root directory and keep
that terminal open:

```bash
conda activate gpumdkit
cd /path/to/working-directory
gpumdkit.sh -server
```

On the local workstation, open a second terminal and forward the remote
loopback port:

```bash
ssh -N -L 8888:127.0.0.1:8888 <user>@<host>
```

Then open `http://127.0.0.1:8888` locally. The browser is connected to the
remote service through SSH; files remain on the remote host unless explicitly
downloaded. Stop the tunnel with `Ctrl+C` in the local terminal and stop the
server with `Ctrl+C` in the remote terminal.

If local port `8888` is already occupied, use another local port while keeping
the remote port unchanged, for example:

```bash
ssh -N -L 8890:127.0.0.1:8888 <user>@<host>
```

Open `http://127.0.0.1:8890`. If the remote port is occupied, start the
server on another port (for example `gpumdkit.sh -server 9000`) and forward
that remote port (`ssh -N -L 8888:127.0.0.1:9000 <user>@<host>`).

### Command options

```text
gpumdkit.sh -server [port] [-b <address>] [-pw <password>]
gpumdkit.sh -server -h
```

| Option | Default | Purpose |
|---|---|---|
| `port` | `8888` | TCP port for the web service |
| `-b <address>` | `127.0.0.1` | Address to bind; loopback is recommended with an SSH tunnel |
| `-pw <password>` | No password on loopback | Require a password for web login; non-loopback binds require one |

Avoid binding to `0.0.0.0` for routine remote use. That exposes the service
on the host's network interfaces and requires a password; the service does
not provide HTTPS. Prefer an SSH tunnel. Do not place a password in shared
scripts or shell history.

### What the console does (and does not do)

- Browse and filter files under the server's starting directory; preview,
  create, upload, download, and edit supported files.
- Run the configured GPUMDkit plotting actions and view generated figures.
- Read NEP and MD output files for progress estimates. It does not launch a
  simulation, submit scheduler jobs, or prove that a process is still alive.
- The terminal runs real shell commands on the remote host. It is not confined
  to the file-browser root; only use it with trusted users and keep the service
  private.

The full bilingual user guide is in
[`docs/tutorials/en/remote_web_console.md`](../../../docs/tutorials/en/remote_web_console.md)
and [`docs/tutorials/zh/远程网页控制台.md`](../../../docs/tutorials/zh/远程网页控制台.md).

## File layout

```
Scripts/utils/server/
├── server.py        # HTTP server, security, API, recommendation rules
└── web/             # frontend assets, served at /, /style.css, /app.js
    ├── index.html   # topbar, file toolbar, sidebar, content canvas, overlays
    ├── style.css    # design tokens (light + dark) and components
    └── app.js       # all frontend logic
```

`web_asset(name)` caches file contents keyed by mtime and reloads them
when the file changes on disk, so editing `app.js` or `style.css` takes
effect on the next browser refresh without a server restart. The frontend
files are real files; do not inline them back into Python strings.

## HTTP API

All `/api/*` endpoints except `/api/login` require a session. Auth is a
cookie (`gpumdkit_token`, HttpOnly, SameSite=Strict) or an `X-Auth-Token`
header. Without a configured password (loopback bind, no `-pw`), all
endpoints are open.

| Endpoint | Method | Purpose |
|---|---|---|
| `/` , `/style.css`, `/app.js` | GET | Web assets (no auth) |
| `/logo.png`, `/logo_lateral.png`, `/favicon.ico` | GET | Brand logos from `docs/Gallery/` (no auth) |
| `/api/scan?path=<rel>` | GET | Listing, recommendations (+ staleness), training and MD simulation status, subdir summaries, `busy`, server `now` |
| `/api/loss?path=<rel>` | GET | Bounded loss chart data (see below) |
| `/api/file?path=<rel>[&partial=1]` | GET | Text preview, limit 512 KB, `X-File-Mtime` header; `partial=1` returns head/tail JSON for oversized files |
| `/api/image?path=<rel>` | GET | PNG/JPEG only, limit 20 MB |
| `/api/download?path=<rel>` | GET | Download one file from the sandbox |
| `/api/login` | POST | `{password}` -> session cookie (5 fails -> 60 s lockout) |
| `/api/run` | POST | `{action, path}` runs an allowlisted command |
| `/api/exec` | POST | `{path, command}` runs a shell command (120 s timeout) |
| `/api/save` | POST | `{path, content, base_mtime?}` atomic overwrite; 409 conflict when the file changed externally |
| `/api/create` | POST | `{path, name, kind}` creates an empty file or directory without overwriting |
| `/api/upload?path=<rel>&name=<name>` | POST | Upload a streamed `application/octet-stream` body; existing names are refused |

All timestamp arithmetic must use the server-provided `now` field, never
the client clock (client and server may differ through SSH tunnels).

Create and upload operations are limited to the resolved working-directory
sandbox and refuse to overwrite an existing entry. Uploads are streamed to a
temporary sibling file before an exclusive link makes the completed file
visible.

## Loss data path (`/api/loss`)

`loss.out` can grow large, so it is never sent through the preview API.
The server keeps a cache per file (keyed by size and mtime) and parses
only appended complete lines; a shrinking file triggers a full reparse.
A half-written trailing line is excluded and re-read on the next request.
The response carries at most 600 sampled points (bucket sampling that
keeps the first and last point of each bucket), the record count,
`first_gen`/`last_gen`, `multi_run` detection (generation numbers that
restart), and `x_mode`: `generation` when the file is a single
monotonic run, otherwise `record #`. Invalid or incomplete lines are
counted in `skipped`. Known boundary: a rewrite that keeps the file at
least as large as the cached offset is treated as an append.

`training_status()` uses the same cache: target generations come from
`nep.in` (`has_target` is false when the file has no `generation` line
and the default 100000 is assumed), progress uses the last complete
valid record, and `loss_mtime` is the evidence time shown to the user.
File mtime shows when data was written, not whether a process is alive.

## MD simulation monitoring (`/api/scan`)

When `neighbor.out` exists, the monitor parses its complete records
incrementally and uses the reported step plus the strict sum of `run <integer>`
commands in `run.in` for progress. The chart is capped at 1200 sampled points
and plots radial/angular `actual` counts with the corresponding `max` values
reported in the same file. Dashed max lines are references, not a standalone
physical-stability criterion. Recent observed step deltas provide a speed and
ETA estimate; the first sample has no ETA yet.

Only when `neighbor.out` is absent does the monitor use `thermo.out`, matching
`gpumdkit.sh -time gpumd`: it multiplies the count of non-empty, non-comment
rows by the first positive `dump_thermo` interval. This fallback is approximate.
The API does not launch or inspect a GPUMD process, and file mtime alone is not
treated as evidence that a simulation is running.

## MD 模拟监测（`/api/scan`）

存在 `neighbor.out` 时，监测器会增量解析完整记录，并使用其中报告的步数
以及 `run.in` 中严格匹配的 `run <integer>` 总和计算进度。图表最多返回
1200 个采样点，绘制径向/角向 `actual` 数量及同一文件报告的对应 `max` 值。
`max` 虚线仅作参考，不能单独作为物理稳定性判据。速度和 ETA 根据最近观测到的
步数变化估算；首次采样时尚无 ETA。

只有在 `neighbor.out` 不存在时，监测器才使用 `thermo.out`，与
`gpumdkit.sh -time gpumd` 的规则一致：非空且非注释行数乘以第一个正数
`dump_thermo` 间隔得到近似进度。API 不会启动或检查 GPUMD 进程，也不会仅凭文件
修改时间判断模拟正在运行。

## Security model (do not weaken)

- Default bind `127.0.0.1`; non-loopback binds require a password
  (PBKDF2-hashed in memory, never written to disk).
- `resolve_in_root()` sandboxes every path: `realpath` + prefix check,
  so `..` and symlink escapes are rejected. This protects the file
  APIs and the run/exec working directory; the terminal itself runs
  arbitrary shell commands and is **not** a filesystem sandbox.
- `/api/run` executes only the fixed `RUN_ACTIONS` allowlist with
  server-side validation.
- `/api/exec` is arbitrary command execution by design (the terminal);
  it sits behind login, the exec lock, a 120 s timeout, and a
  4096-char limit. Keep the default loopback posture.
- `EXEC_LOCK` serializes `/api/run` and `/api/exec`; the scan response
  exposes its state as `busy` so the UI can explain refusals.
- All responses set `X-Content-Type-Options: nosniff` and
  `Cache-Control: no-store`.

## Recommendation engine

`RUN_ACTIONS` maps an action id to `argv`, required input files
(`requires`), glob requirements (`requires_glob`), batch requirements
(`batch`: temperature subdirectories ending in `K` with per-dir files),
the produced PNG, and a timeout. `build_recommendations()` turns
detected files into cards and marks `stale` when an existing result is
older than the newest input (mtime comparison; staleness is a hint, not
a scientific validity claim).

The required/produced file names were verified against
`Scripts/plt_scripts/*.py` by reading their `np.loadtxt`/`savefig`
calls:

| action | requires | produces |
|---|---|---|
| `plt_msd` | `msd.out` | `msd.png` |
| `plt_sdc` | `msd.out` | `sdc.png` |
| `plt_msd_sdc` | `msd.out` (7 columns) | `msd_sdc.png` |
| `plt_vac` | `sdc.out` | `vac.png` |
| `plt_thermo` | `thermo.out` | `thermo.png` |
| `plt_train` | `loss.out energy_train.out force_train.out` | `train.png` |
| `plt_train_density` | same as `plt_train` | `train_density.png` |
| `plt_train_test` | train/test energy/force/stress `.out` (6 files) | `train_test.png` |
| `plt_prediction` | `energy_train.out force_train.out` | `prediction.png` |
| `plt_msd_conv` | `msd_step*.out` (glob) | `msd_convergence.png` |
| `plt_sigma` | >= 2 `NNNK` dirs each with `msd.out thermo.out` | `Arrhenius_sigma.png` |
| `plt_D` | >= 2 `NNNK` dirs each with `msd.out` | `Arrhenius_D.png` |

Commands that need extra user arguments (e.g. `plt_msd_all`) are
excluded on purpose. When adding an action, verify the real script
behavior first and mirror it here.

## Frontend architecture (web/app.js)

Layout: one content canvas per directory. Topbar (logo, connection status,
refresh interval, theme, more menu) / location and file-action toolbar
(home, parent, breadcrumbs, create, upload, selected-file download) / left
file browser (current-directory name and item counts, filter, collapsible) / one
scrolling content column whose sections appear only when
their data exists: training curve (only when loss.out is present and
parseable), MD monitor (when neighbor.out or thermo.out is present), figures
(only when images exist), available analyses (only when recommendations
exist), plain file listing (only for directories without loss, MD status, or
images), subdirectory list, or a single empty state.
There are no fixed view tabs, no bottom panel, no persistent status
bar; secondary surfaces are on-demand overlays: a lightbox for images,
a preview sheet for files, a chart zoom layer, and a command drawer
opened from the more menu or when a command runs.

Polling is off by default and follows the interval selected in the topbar.
When the page is hidden, scheduled scans are skipped and then resume when it
is visible again.

State globals: `curPath`, `rootDir`, `connOk`, `lastGoodAt`,
`lastScanServerNow`, `lastFiles`, `lastSubdirs`, `lastPath`,
`lastRecs`, `lastTraining`, `lastSimulation`, `nepSamples`, `lossView` (chart data plus
per-directory scale state `logX`/`logY`/`hidden`), `plotState`,
`cmdHistory`, `viewer`, `drawerTab`, `figShown`, `lbList`/`lbIndex`.

Rendering: `loadScan()` discards out-of-order responses via `scanSeq`,
resets cross-directory baselines when the path changes, and re-renders
browser plus canvas. `fetchLoss(seq)` captures its own path and
sequence; late responses never draw into another directory. Polling is
adaptive (selected interval while jobs are active, 60 s idle check).

Loss chart: geometry is computed from the measured content width
(`lossChartSize`): side-by-side curve plus summary at >= 860 px
(curve 600-760 px wide, ratio ~1.85, capped on wide screens), stacked
with full container width below that; a ResizeObserver redraws on
layout changes. Scales default to log-log but are switchable per axis
with a visible summary (`X Log · Y Log`); the legend toggles series;
hover shows real sampled values; missing values break the line; a zoom
layer re-renders at viewport size. Non-positive values in log mode are
counted and reported with a linear-switch link. Training summary
(target, progress, ETA) appears only when nep.in exists; loss curves
render without nep.in.

Viewer (preview sheet): read-only by default with an explicit edit
action; one draft state shared across the form and text modes;
switching files with unsaved changes asks first; saves carry
`base_mtime` from open time and a 409 conflict offers reload, keep
draft, or explicit overwrite; files above 512 KB get a bounded
head/tail preview with editing disabled. Numeric preview is a table
with an expandable plot options block (X/Y column selection, log
axes); known column names only for reliably identified formats.

Execution records: command, directory, start time, duration, status
(ok / exit code / failed / timeout / network), output, truncation
flag, and result file bound to the directory where the command ran.
Network errors are reported as unknown state. History is session-only.

Terminal: each command is an independent `bash -c`; session state such
as `export` is not kept, `cd`/`clear`/`pwd` are intercepted and sync
the UI. There is no cancel button because the backend cannot stop a
running child reliably.

## Validation checklist (run every time you change code)

```bash
bash -n gpumdkit.sh
python3 -B -c 'from pathlib import Path; p=Path("Scripts/utils/server/server.py"); compile(p.read_text(), str(p), "exec")'
node --check Scripts/utils/server/web/app.js
git diff --check
```

For logic changes, smoke-test the endpoints with curl in a temp
directory (4-column `msd.out`, `nep.in` plus multi-column `loss.out`,
a > 512 KB text file), then start the server, log in, and walk through
scan, loss, file preview (normal and partial), run, exec, and save
conflict. Remove every temporary file afterward. Node-based stub tests
work well for pure functions: stub `document`/`window`, `eval` app.js,
and call `parseNumeric`, `drawSeries`, `lossChartSize`, `trainingRate`.

Browser evidence can be produced with headless Chrome (adjust the path
per machine):

```bash
"/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" \
  --headless=new --disable-gpu --hide-scrollbars \
  --window-size=1440,900 --virtual-time-budget=5000 \
  --screenshot=shot.png "http://127.0.0.1:8888/"
```

Notes: `--window-size` must use a comma (`1440,900`); an `x` separator
is silently ignored and renders at an unrelated viewport. Use
`--dump-dom` plus section checks for structural acceptance, and append
`?theme=dark` to capture the dark theme. The canvas `style` attribute
in the dumped DOM proves the JS geometry code ran. Interactive flows
(lightbox keyboard navigation, preview editing) still need a real
browser session; mark them unverified otherwise.

## Pitfalls already fixed (do not regress)

- **Browser-global collisions**: `var history` silently loses to
  `window.history` (read-only). The variable is `cmdHistory`. Before
  adding a top-level `var`/`function` in `app.js`, check it is not a
  `window` property (`history`, `name`, `status`, `event`, `origin`,
  `close`, `open`, `top`, `self`, `length`, `screen`, ...).
- **`.hidden` vs later `display` rules**: the utility class must stay
  `.hidden { display: none !important; }`.
- **Clock skew**: compare mtimes against the server `now` only.
- **Stale async responses**: guard every fetch continuation against
  current state (`viewer.rel`, `scanSeq`); capture the working
  directory when a run starts (`runAction`, `saveDraft`) and bind
  results to it. The loss fetch must capture its own path and sequence
  per call — an earlier version relied on a shared mutable
  `lossFetchPath` global, which let a late response from directory A
  render into directory B.
- **Cross-directory baselines**: file growth and ETA samples must be
  reset when the directory changes.
- **ASCII only, no emojis** in code and terminal output; user-facing
  `print`/`echo` lines start with a leading space (GPUMDkit
  convention).

## Integration points

- `gpumdkit.sh`: the `-server` case routes to
  `Scripts/utils/server/server.py`; update the help table row and
  `Scripts/utils/completion.sh` if the CLI surface changes.
- Logos: `docs/Gallery/gpumdkit_logo.png` (square, login + favicon)
  and `gpumdkit_logo_lateral.png` (header). Missing files degrade to
  text.
- `docs/updates.info`: record user-visible feature changes and bug
  fixes here (not documentation-only changes); keep the five newest
  entries.
- Repository-wide agent conventions: see `AGENTS.md` and
  `skills/gpumdkit-skill/references/contributing.md`.
