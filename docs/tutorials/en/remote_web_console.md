<div align="center">
  <h1>Remote Web Console</h1>
  <p style="text-align: justify;">Browse your simulation and training folders, draw figures, follow NEP and MD progress, and edit small input files from a local browser, through an SSH tunnel.</p>
</div>

## What it is

`gpumdkit.sh -server` starts a browser-based console rooted at the directory from which it is launched. The browser displays files that **stay on that machine**; nothing is copied to your workstation unless you download it. Use it when you want to:

- check whether a NEP training or an MD run is going well, without logging in and reading `loss.out` by hand;
- draw the usual figures (training, MSD, thermo, Arrhenius) with one click and look at them immediately;
- change a keyword in `nep.in` or `run.in` and save it safely;
- run a quick shell command on the remote host.

It is **not** a job launcher: the console never starts GPUMD or NEP, never submits scheduler jobs, and cannot tell whether a process is still alive. It reads the output files that your runs write.

## Start the server on the remote host

Connect with SSH, activate the environment that contains GPUMDkit and its plotting packages, change to the directory you want to expose, and start the server:

```bash
ssh <user>@<host>
conda activate gpumdkit
cd /path/to/working-directory
gpumdkit.sh -server
```

Keep this terminal open. By default the server listens on `127.0.0.1:8888` on the remote host. The starting directory is the root of the file browser; the interface can open its sub-folders, but the file APIs cannot browse above it.

## Connect from your workstation

Open a second terminal **on your own computer** and create the tunnel:

```bash
ssh -N -L 8888:127.0.0.1:8888 <user>@<host>
```

Leave it running and open this address in your local browser:

```text
http://127.0.0.1:8888
```

```text
 your workstation                                  remote host
┌─────────────────────────┐   encrypted SSH tunnel   ┌──────────────────────────────┐
│ browser                 │ ◄──────────────────────► │ gpumdkit.sh -server          │
│ http://127.0.0.1:8888   │  ssh -N -L 8888:…:8888   │ 127.0.0.1:8888, files stay   │
└─────────────────────────┘                          └──────────────────────────────┘
```

The first `8888` is the **local** port, the last one is the **remote** server port. If a port is busy:

| Busy port | What to do |
|---|---|
| Local | Choose another local port: `ssh -N -L 8890:127.0.0.1:8888 <user>@<host>`, then open `http://127.0.0.1:8890`. |
| Remote | Start the server on another port, `gpumdkit.sh -server 9000`, and forward it: `ssh -N -L 8888:127.0.0.1:9000 <user>@<host>`. |

Stop the tunnel with `Ctrl+C` in the local terminal and the server with `Ctrl+C` in its remote terminal.

## A tour of the interface

The page has one content area per folder. **Sections appear only when their data exists**: a folder without `loss.out` shows no loss curve, a folder without images shows no figure strip. Hover a number below (or its entry on the right) to highlight the region, and click an image to enlarge it.

### A NEP training folder

<!-- annotated-shot: training -->

### A GPUMD run folder

<!-- annotated-shot: md -->

The screenshots show example runs; your folders will show your own data.

## What appears for which files

| Files in the folder | What the console shows |
|---|---|
| `loss.out` | The loss curve. With `nep.in` also a progress bar (generation against the `generation` target) and a parameter panel. |
| `neighbor.out` and `run.in` | The MD monitor: neighbor counts against their maxima, and progress against the `run` steps of `run.in`. |
| `thermo.out` (no `neighbor.out`) | An approximate MD progress taken from the number of rows and `dump_thermo`, the same rule as `gpumdkit.sh -time gpumd`. |
| `.png` / `.jpg` images | A strip of thumbnails; click one for the full-size viewer (`←` `→` to move, `Esc` to close). |
| Files needed by a plot (see below) | An **Actions** card for that plot. |
| Sub-folders | A list with item counts, so you can jump into `300K/`, `400K/`, … |

## Plot actions

Each card in **Actions** runs a fixed GPUMDkit plot script in the current folder. The card tells you whether the figure has been generated and when; **Run** (re)generates it and **View** opens it. A figure older than its input files is marked as out of date. **Run all** runs every card.

| Action | Needs in the folder | Produces |
|---|---|---|
| MSD | `msd.out` | `msd.png` |
| SDC | `msd.out` | `sdc.png` |
| MSD & SDC | `msd.out` with seven columns | `msd_sdc.png` |
| MSD convergence | `msd_step*.out` | `msd_convergence.png` |
| VAC | `sdc.out` | `vac.png` |
| Thermo | `thermo.out` | `thermo.png` |
| Training | `loss.out`, `energy_train.out`, `force_train.out` | `train.png` |
| Training density | the same files as Training | `train_density.png` |
| Train and test | the six train/test energy, force and stress `.out` files | `train_test.png` |
| Prediction | `energy_train.out`, `force_train.out` | `prediction.png` |
| Arrhenius conductivity | at least two `NNNK` folders, each with `msd.out` and `thermo.out` | `Arrhenius_sigma.png` |
| Arrhenius diffusivity | at least two `NNNK` folders, each with `msd.out` | `Arrhenius_D.png` |

Plots that need extra arguments, for example `-plt msd_all` with a list of elements, are not offered as cards; run them from the terminal drawer instead.

### Example: an Arrhenius plot from several temperatures

Keep one folder per temperature, named like `300K`, `400K`, `500K`, each containing the `msd.out` (and `thermo.out` for conductivity) of that run:

```text
diffusion/
├── 300K/   msd.out  thermo.out
├── 400K/   msd.out  thermo.out
├── 500K/   msd.out  thermo.out
└── 600K/   msd.out  thermo.out
```

Start the server in `diffusion/` (or open it in the file browser). With at least two temperature folders the **Arrhenius** cards appear. Click **Run**, then **View**. The figures are written next to the data (`Arrhenius_D.png`, `Arrhenius_sigma.png`). Whether the fit is meaningful, for example that every MSD is converged, is for you to check: use the MSD and MSD convergence cards inside each temperature folder first.

## Follow a NEP training

Open the training folder (the one holding `nep.in` and `loss.out`).

- **Loss curve.** Both axes are logarithmic by default; the **X Log** and **Y Log** buttons switch each axis to linear. Click a legend entry to hide or show a loss term, hover to read sampled values, and use the zoom button for a larger view. When the generation numbers restart (a continued run), the x axis becomes the record number.
- **Progress.** The target is the `generation` value in `nep.in`; if `nep.in` has none, the NEP default of 100000 is assumed and the card says so. Speed and ETA appear after the console has seen at least two samples of `loss.out`.
- **Parameters.** The panel on the right lists the keywords of `nep.in`; **open file** shows the whole file.

Set **auto refresh** in the top bar (5 s to 60 s) to keep the view current. It is off by default and pauses while the browser tab is hidden. A rising progress bar means `loss.out` is receiving new records; it does not prove that the process is alive.

## Follow an MD run

In a folder with `neighbor.out`, the monitor plots the radial and angular neighbor counts together with the maxima that `neighbor.out` reports. The dashed maximum lines are references, not a stability test. Progress compares the MD step with the sum of the `run` commands in `run.in`; ETA needs two samples.

## Edit files safely

Click a file name to open the preview sheet. Text files up to 512 KB are shown in full; larger files show their head and tail and cannot be edited. Numeric files are shown as a table with an optional plot.

<div class="img-pair">
  <img src="../../assets/img/console-edit-text.png" alt="Editing nep.in as text in the web console">
  <img src="../../assets/img/console-edit-form.png" alt="Editing nep.in in form mode in the web console">
</div>

Press **Edit** to make the file editable. For keyword files such as `nep.in` and `run.in` you can switch between **Text** and **Form** views; both edit the same draft. **Save** writes the file atomically.

If the file changed on disk after you opened it, for example because another process rewrote it, the save is refused and you choose between **Reload**, **Keep my draft** and **Overwrite**. Switching to another file with unsaved changes asks first.

**New** creates an empty file or folder, **Upload** sends files from your workstation, and **Download** fetches the selected remote file. Creating and uploading never overwrite an existing name.

> Editing `nep.in` does not restart anything. Start or restart `nep` and GPUMD yourself, from your own terminal or scheduler, and let the console show the new output.

## The command drawer

The terminal icon in the top bar opens the **Commands & output** drawer.

- **Output** lists what the Run buttons did: the command, the folder it ran in, the duration, the status and the output.
- **Terminal** runs shell commands on the remote host.

<div align="center">
  <img src="../../assets/img/console-terminal.png" alt="The terminal tab of the command drawer" width="46%">
</div>

Each command runs on its own with `bash -c` in the current folder and is stopped after 120 s. State such as `export` is not kept between commands, and `cd`, `pwd` and `clear` are built in. There is no cancel button. Because the terminal runs real commands, it is **not** limited to the folder shown in the file browser.

## Server options and safety

```text
gpumdkit.sh -server [port] [-b <address>] [-pw <password>]
gpumdkit.sh -server -h
```

| Option | Default | Description |
|---|---|---|
| `port` | `8888` | Port used by the web server |
| `-b <address>` | `127.0.0.1` | Network address to listen on |
| `-pw <password>` | Optional on loopback | Login password; required for any other address |

- **Prefer the SSH tunnel with the default loopback address.** Binding to another address such as `0.0.0.0` exposes the console on the host's network, requires a password (you are asked for one if `-pw` is omitted), and does not use HTTPS.
- After five wrong passwords the login is locked for 60 s.
- Only the file APIs are confined to the starting directory. The terminal is not, so run the server only for users you trust.
- Keep passwords out of shared scripts and shell history.

## Troubleshooting

| Symptom | What to check |
|---|---|
| The browser cannot connect | The remote server **and** the local tunnel must both be running, and the forwarded ports must match. |
| The wrong files appear | The root is the directory where `gpumdkit.sh -server` started. Stop it and start it from the intended directory. |
| An action card is missing | The folder lacks one of the required files in the table above. |
| A plot action fails | Start the server from the environment that has GPUMDkit and its plotting packages; read the error in the **Output** tab. |
| Progress does not change | It changes only when `loss.out` or the MD output receives new records. The console is not a process or scheduler monitor. |
| A save is refused | The file changed on disk. Choose Reload, Keep my draft or Overwrite. |
| A command seems to hang | Commands are stopped after 120 s; long jobs belong in your own terminal or scheduler. |
