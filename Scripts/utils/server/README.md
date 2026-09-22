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

## File layout

```
Scripts/utils/server/
├── server.py        # HTTP server, security, API, recommendation rules
└── web/             # frontend assets, served at /, /style.css, /app.js
    ├── index.html   # topbar, browser, view tabs, bottom panel, status bar
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
| `/api/scan?path=<rel>` | GET | Listing, recommendations (+ staleness), training status, subdir summaries, `busy`, server `now` |
| `/api/loss?path=<rel>` | GET | Bounded loss chart data (see below) |
| `/api/file?path=<rel>[&partial=1]` | GET | Text preview, limit 512 KB, `X-File-Mtime` header; `partial=1` returns head/tail JSON for oversized files |
| `/api/image?path=<rel>` | GET | PNG/JPEG only, limit 20 MB |
| `/api/login` | POST | `{password}` -> session cookie (5 fails -> 60 s lockout) |
| `/api/run` | POST | `{action, path}` runs an allowlisted command |
| `/api/exec` | POST | `{path, command}` runs a shell command (120 s timeout) |
| `/api/save` | POST | `{path, content, base_mtime?}` atomic overwrite; 409 conflict when the file changed externally |

All timestamp arithmetic must use the server-provided `now` field, never
the client clock (client and server may differ through SSH tunnels).

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

Layout: topbar (logo, breadcrumb, connection status, refresh interval,
theme toggle) / left file browser (unified dirs+files, filter, name/
modified/size sort, path copy) / main view tabs (Overview, Data,
Figures) / collapsible bottom panel (Output, History, Terminal) /
status bar. View state is independent of the current directory and
survives auto refresh. Narrow screens collapse the browser into a
drawer.

State globals: `curPath`, `rootDir`, `viewMode`, `running`, `termBusy`,
`busyFlag`, `scanSeq`, `pollActive`, `connOk`, `lastGoodAt`,
`lastScanServerNow`, `lastFiles`, `lastSubdirs`, `lastPath`,
`nepSamples`, `lossData`, `cmdHistory`, `figRendered`, `panelTab`,
`viewer`, `pendingOpen`, `drawerOpen`.

Rendering: `loadScan()` fetches `/api/scan`, discards out-of-order
responses via `scanSeq`, resets cross-directory baselines (growth
deltas, ETA samples, loss data) when the path changes, and re-renders
only the current view. Polling is adaptive: the selected interval while
jobs are active (any file or subdir written within 180 s, or unfinished
training), a 60 s idle check otherwise. Connection failures keep the
old data, mark it stale in the topbar, and recover silently.

Training monitor: distinguishes updating, no-recent-update, finished,
read failure, and empty states, with evidence lines from server time.
ETA is the median of recent positive rate samples; samples reset on
directory change or generation regression; without an explicit target
or enough samples no precise ETA is shown.

Viewer (Data view): one draft state shared across the Form and Edit
tabs (switching tabs never restores stale disk content); switching
files with unsaved changes asks first. Saves carry `base_mtime` from
open time; a 409 conflict offers reload, keep-draft, or explicit
overwrite. Files above 512 KB get a bounded head/tail preview with
editing disabled. The Data tab has X/Y column selection, multi-series
plotting, linear/log axes, hover value readout, and known column names
only for reliably identified formats (`msd.out`, `sdc.out`,
`loss.out`); unknown formats show `c0..cN`.

Execution records: command, directory, start time, duration, status
(ok / exit code / failed / timeout / network), output, truncation flag,
and result file bound to the directory where the command ran. Network
errors are reported as unknown state, never as command failure.
History is session-only by design; do not persist terminal commands.

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
and call `parseNumeric`, `drawSeries`, `trainingRate`, `schedule`.

Browser-only behavior (layout at 1440/1024/390 px, keyboard focus,
theme switching, reduced motion) should be verified manually when a
browser is available; mark it as unverified otherwise.

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
  current state (`viewer.rel`, `scanSeq`, `lossFetchPath`); capture
  the working directory when a run starts (`runAction`, `saveDraft`)
  and bind results to it.
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
