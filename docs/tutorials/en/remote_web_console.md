<div align="center">
  <h1>Remote Web Console</h1>
  <p style="text-align: justify;">Browse and work with files on a remote GPUMDkit host from your local browser through an SSH tunnel.</p>
</div>

## What it is

`gpumdkit.sh -server` starts a browser-based console rooted at the directory
from which it is launched. The browser displays files that remain on that
machine; it does not copy the working directory to your workstation. Use this
when you want to inspect remote simulation/training outputs, run supported
plot actions, or perform small file operations without transferring the whole
project.

## Start the server on the remote host

Connect to the host with SSH, activate the environment containing GPUMDkit,
change to the directory you want to expose in the console, and start the
server:

```bash
ssh <user>@<host>
conda activate gpumdkit
cd /path/to/working-directory
gpumdkit.sh -server
```

Replace the placeholders with your SSH account, host, and actual working
directory. If you do not use Conda, activate or load the environment that
provides `gpumdkit.sh` and the Python packages required by the plot actions.
Keep this remote terminal open while using the console.

By default, the server listens on `127.0.0.1:8888` on the remote host. The
starting directory is the file-browser root; the UI can navigate its
subdirectories, but the file APIs cannot browse above it.

## Create the SSH tunnel from your workstation

Open another terminal on your local workstation and run:

```bash
ssh -N -L 8888:127.0.0.1:8888 <user>@<host>
```

Use the same account and host as for the remote login. Leave this tunnel
terminal running, then open this address in your local browser:

```text
http://127.0.0.1:8888
```

The first `8888` is the local port; the last `8888` is the remote server port.
If local port `8888` is busy, choose another local port, such as `8890`:

```bash
ssh -N -L 8890:127.0.0.1:8888 <user>@<host>
```

Then open `http://127.0.0.1:8890`. If the server uses a different remote port,
for example `9000`, start it with `gpumdkit.sh -server 9000` and forward that
port with `ssh -N -L 8888:127.0.0.1:9000 <user>@<host>`.

Stop the tunnel with `Ctrl+C` in the local terminal. Stop the web server with
`Ctrl+C` in its remote terminal.

## Server options

```text
gpumdkit.sh -server [port] [-b <address>] [-pw <password>]
gpumdkit.sh -server -h
```

| Option | Default | Description |
|---|---|---|
| `port` | `8888` | Port used by the web server |
| `-b <address>` | `127.0.0.1` | Network address to listen on |
| `-pw <password>` | Disabled on loopback | Require a password to sign in |

For an SSH tunnel, keep the default loopback bind. Do not bind to `0.0.0.0`
unless direct network access is intentional and permitted by the host's
security policy. A non-loopback bind requires a password; the web server does
not provide HTTPS. When binding beyond loopback, omitting `-pw` prompts for a
password interactively. Do not put passwords in shared scripts or shell
history.

## What you can do in the browser

- Navigate the current directory, filter its files, and open supported text
  and image previews.
- Create files/folders, upload files from the workstation, and download
  selected remote files.
- Run the available GPUMDkit plot actions and view their output figures.
- Inspect NEP training or MD progress when the corresponding output files are
  present. The console reads output files; it does not launch GPUMD, submit
  scheduler jobs, or determine process liveness from file timestamps.
- Use the terminal drawer for quick commands. Commands run on the remote host;
  this terminal is not restricted to the file-browser root.

Start the server from the Conda environment you intend to use. This ensures
that allowlisted plot actions inherit the expected GPUMDkit installation and
Python dependencies. The browser's upload/download controls are the explicit
ways to transfer files between your workstation and the remote host.

## Troubleshooting

| Symptom | Check |
|---|---|
| Browser cannot connect | Confirm the remote server and local SSH tunnel are both still running, and that the forwarded ports match. |
| The wrong files appear | The server root is the directory where `gpumdkit.sh -server` was started. Stop it and relaunch from the intended directory. |
| Plot actions fail | Start the server from the environment containing GPUMDkit and the plot dependencies; inspect the action output in the console. |
| Progress does not change | The page rereads output files, but progress changes only when those files receive new records. The monitor is not a process or scheduler monitor. |

