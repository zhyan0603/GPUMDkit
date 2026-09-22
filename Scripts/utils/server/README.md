# GPUMDkit Web Console (`Scripts/utils/server/`)

A zero-dependency web console for GPUMDkit: it scans a working directory,
recommends allowlisted GPUMDkit plot commands, runs them, monitors running
NEP trainings and MD jobs, renders figures, previews data files, edits small
keyword files, and provides a web terminal.

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
    ├── index.html   # page skeleton (topbar, rail, sidebar, status bar)
    ├── style.css    # design tokens (light + dark) and components
    └── app.js       # all frontend logic
```

`web_asset(name)` in `server.py` caches file contents keyed by mtime and
reloads them when the file changes on disk. Editing `app.js` or `style.css`
therefore takes effect on the next browser refresh without restarting the
server. Use this when iterating on the UI.

`index.html` references assets via absolute paths (`/style.css`, `/app.js`).
Do not inline them back into Python strings; the single-string design was
deliberately removed because it was fragile to edit.

## HTTP API

All `/api/*` endpoints except `/api/login` require a session. Auth is a
cookie (`gpumdkit_token`, HttpOnly, SameSite=Strict) or an `X-Auth-Token`
header. Without a configured password (loopback bind, no `-pw`), all
endpoints are open.

| Endpoint | Method | Purpose |
|---|---|---|
| `/` , `/style.css`, `/app.js` | GET | Web assets (no auth) |
| `/logo.png`, `/logo_lateral.png`, `/favicon.ico` | GET | Brand logos from `docs/Gallery/` (no auth) |
| `/api/scan?path=<rel>` | GET | Directory listing, recommendations, training status, subdir summaries, server `now` |
| `/api/file?path=<rel>` | GET | Text preview, limit 512 KB, rejects binary |
| `/api/image?path=<rel>` | GET | PNG/JPEG only, limit 20 MB |
| `/api/login` | POST | `{password}` -> sets session cookie (5 fails -> 60 s lockout) |
| `/api/run` | POST | `{action, path}` runs an allowlisted command |
| `/api/exec` | POST | `{path, command}` runs an arbitrary shell command (120 s timeout) |
| `/api/save` | POST | `{path, content}` atomically overwrites an existing file, limit 200 KB |

Response of `/api/scan` (fields the frontend relies on): `dirs`, `files`
(`name/size/mtime`), `subdirs` (`name/count/newest/msd/thermo/loss`),
`recommendations`, `training` (`total/done/loss_file/finished` or `null`),
`truncated`, `path`, `root`, and `now` (server clock; always use it for
mtime comparisons, never the client clock).

## Security model (do not weaken)

- Default bind `127.0.0.1`; non-loopback binds require a password
  (PBKDF2-hashed in memory, never written to disk).
- `resolve_in_root()` sandboxes every path: `realpath` + prefix check, so
  `..` and symlink escapes are rejected.
- `/api/run` executes only the fixed `RUN_ACTIONS` allowlist; parameters
  are validated server-side. Never add an endpoint that takes a raw
  command from the recommendation UI.
- `/api/exec` is arbitrary command execution **by design** (the terminal);
  it sits behind login, the exec lock, a 120 s timeout, and a 4096-char
  limit. Keep the default loopback posture.
- `EXEC_LOCK` serializes `/api/run` and `/api/exec` (plot scripts write
  fixed filenames such as `msd.png`; concurrent runs would clobber).
- All responses set `X-Content-Type-Options: nosniff` and `Cache-Control:
  no-store`.

## Recommendation engine

`RUN_ACTIONS` maps an action id to `argv` (gpumdkit.sh arguments), required
input files (`requires`), glob requirements (`requires_glob`), batch
requirements (`batch`: temperature subdirectories ending in `K` with
per-dir files), the produced PNG, and a timeout.

`build_recommendations()` turns detected files into recommendation cards.
The required/produced file names were verified against the plot scripts
(`Scripts/plt_scripts/*.py`) by reading their `np.loadtxt`/`savefig` calls:

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

Column counting (`count_columns`) distinguishes 4-column (GPUMDkit
`-calc msd`) from 7-column (GPUMD `compute_msd`) `msd.out`. When adding an
action, verify the real script behavior first and mirror it here; commands
that need extra user arguments (e.g. `plt_msd_all`) are excluded on
purpose.

`training_status()` parses `generation` from `nep.in` (default 100000) and
the last generation from `loss.out`. Note `loss.out` is append-mode across
runs, so progress is approximate if the file archives multiple runs.

## Frontend architecture (web/app.js)

State globals: `curPath`, `rootDir`, `running`, `cmdHistory`, `lastFiles`,
`lastSubdirs`, `lastNep`, `lastPollTs`, `viewer`, `lossText`, `scanSeq`,
`termBusy`, `pollActive`, `themeManual`, `openModals`.

Flow: `loadScan(path)` fetches `/api/scan`, then renders breadcrumb,
sidebar (dirs + files), status bar, monitor card, recommendations, and the
figure grid. A `scanSeq` counter discards out-of-order responses (poll vs
navigation race). `pollActive` is computed from file freshness (any file
or subdir touched within 180 s, or unfinished training) using the server
`now` field; `schedule()` then polls at the user interval when active and
at 60 s when idle.

Viewer: clicking a file opens a modal with tabs Data (numeric table, click
a column header to plot), Form (keyword editor, preserves comments and
line endings), Text, Edit (raw, 200 KB limit). `viewer.rel` is re-checked
when a fetch resolves to prevent stale responses from cross-wiring file
names and content (this was a real data-corruption bug).

Loss chart: mirrors `plt_train.py` (`loglog(loss[:, 1:7])`): x axis is the
row index (1-based), not column 0, because `loss.out` is append-mode;
6 series for 10-column NEP files with matplotlib tab10 colors; log-log
axes; a pulsing dot marks the latest Total point.

Theme: light first, follows `prefers-color-scheme`, manual toggle persists
in `localStorage` (`gk_theme`). All colors are CSS variables in
`style.css` (`:root` and `body.dark`). Canvas charts read `--grid`/`--axis`
via `getComputedStyle` at draw time and redraw on theme changes.

Terminal: `cd`/`clear`/`pwd` are intercepted client-side and sync the UI
(`cd` calls `loadScan`); other commands POST `/api/exec`. Input is
disabled while a command is in flight.

## Validation checklist (run every time you change code)

```bash
bash -n gpumdkit.sh
python3 -B -c 'from pathlib import Path; p=Path("Scripts/utils/server/server.py"); compile(p.read_text(), str(p), "exec")'
node --check Scripts/utils/server/web/app.js
git diff --check
```

For logic changes, smoke-test the endpoints with curl in a temp directory
(fake `msd.out` with 4 numeric columns, `nep.in` + multi-column
`loss.out`), then start the server, log in, and walk through
`/api/scan`, `/api/run`, `/api/exec`, `/api/save`. Remove every temporary
file afterward.

Node-based stub tests also work well: stub `document`/`window`, then
`eval` `app.js` and call functions directly (see session history pattern:
`parseNumeric`, `drawSeries`, `applyTheme`, `schedule`).

## Pitfalls already fixed (do not regress)

- **Browser-global collisions**: `var history` silently loses to
  `window.history` (read-only), which broke command history and killed
  `runAction` midway. The variable is now `cmdHistory`. Before adding a
  top-level `var`/`function` in `app.js`, check it is not a `window`
  property (`history`, `name`, `status`, `event`, `origin`, `close`,
  `open`, `top`, `self`, `length`, `screen`, ...).
- **`.hidden` vs later `display` rules**: the utility class must stay
  `.hidden { display: none !important; }`. A later `.overlay { display:
  flex }` once overrode it and made both modals permanently visible with
  a dead Close button.
- **Clock skew**: always compare mtimes against the server-provided `now`,
  not `Date.now()` (client and server may be different machines via SSH
  tunnel).
- **Stale async responses**: guard every fetch continuation against
  current state (`viewer.rel`, `scanSeq`); browsers fire clicks faster
  than networks resolve.
- **`String.prototype` availability**: `endsWith` is fine, but keep ES5
  style elsewhere; no transpiler exists here.
- **ASCII only, no emojis** in code and terminal output; user-facing
  `print`/`echo` lines start with a leading space (GPUMDkit convention).

## Integration points

- `gpumdkit.sh`: the `-server` case routes to
  `Scripts/utils/server/server.py`; update the help table row and
  `Scripts/utils/completion.sh` if the CLI surface changes.
- Logos: `docs/Gallery/gpumdkit_logo.png` (square, login + favicon) and
  `gpumdkit_logo_lateral.png` (header). Missing files degrade to text.
- `docs/updates.info`: record user-visible feature changes and bug fixes
  here (not documentation-only changes); keep the five newest entries.
- Repository-wide agent conventions: see `AGENTS.md` and
  `skills/gpumdkit-skill/references/contributing.md`.
