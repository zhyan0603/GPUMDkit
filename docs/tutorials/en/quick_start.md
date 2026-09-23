<div align="center">
  <h1>🚀 Quick Start</h1>
  <p style="text-align: justify;">This page takes you through installation, the first checks, and a small isolated exercise before pointing you to the next task for your existing data.</p>
</div>

## What it does {#what-it-does}

GPUMDkit brings common computational-materials tasks such as format conversion, structure analysis, property calculation, and visualization into one command-line tool, without requiring custom scripts.

## Before you start {#before-you-start}

Install GPUMDkit first, then check the command-line entry points. Conda is the recommended route:

### 1. Install with Conda (Recommended) {#install-with-conda-recommended}

```bash
conda create -n gpumdkit -c gpumdkit -c conda-forge gpumdkit
conda activate gpumdkit
```

This installs GPUMDkit and its standard dependencies. You do not need to clone the repository or configure environment variables first.

### 2. Check the installation {#verify-the-installation}

```bash
gpumdkit.sh -h
gpumdkit.sh -doctor
```

`-h` lists the available options. `-doctor` checks the Python environment and reports the status of the environment, common packages, and optional packages. `MISS` for an optional package is not an installation failure; install that package only when you need the corresponding feature.

## First practice: convert a POSCAR {#hello-world-example}

To avoid overwriting existing files, do the exercise in a new directory. The exercise runs in a subshell; the `mkdir` command intentionally omits `-p`, so an existing directory makes it fail and `exit 1` leaves the subshell before any existing file can be overwritten.

```bash
(
set -e
mkdir gpumdkit-first-check && cd gpumdkit-first-check || exit 1
cat > POSCAR << 'EOF'
Si
1.0
0 2.715 2.715
2.715 0 2.715
2.715 2.715 0
Si
1
direct
0 0 0
EOF
gpumdkit.sh -pos2exyz POSCAR model.xyz
head model.xyz
)
```

If `model.xyz` shows one Si atom and lattice information, this POSCAR-to-extxyz conversion and output check succeeded. It does not verify every calculator, plotter, or simulation feature.

`extxyz` is extended XYZ: the first line gives the atom count, the second stores lattice and structure-level properties, and later lines store each atom's element, coordinates, and per-atom properties. It is the native training-data format used by NEP.

## Existing data: choose the next task {#common-examples}

If you already have structures or trajectories, you can skip the practice above. Return to the [tutorial index](index.md) and use the Common Tasks table to choose format conversion, structure sampling, analysis, calculators, plotting, or workflow preparation; each page documents its inputs, outputs, and scope.

## Optional packages and agent skills {#optional-packages}

Some features require optional packages:

```bash
pip install neptrain calorine
```

The package also includes the English and Chinese agent skills. Run the following command to find their installed paths and the instructions for configuring them in your agent client:

```bash
gpumdkit.sh -skill
```

## Install from Source (Optional) {#install-from-source-optional}

If you need to use or modify the source code, clone the repository and run the installer:

```bash
git clone https://github.com/zhyan0603/GPUMDkit.git
cd GPUMDkit
source ./install.sh
```

The installer chooses a shell configuration file from `$SHELL`, writes the `GPUMDkit_path` and `PATH` settings there, and loads them for the current shell. It uses `~/.bashrc` by default; when `$SHELL` points to zsh, it uses `~/.zshrc`. If an existing GPUMDkit path is found, the installer prints the old path and asks whether to replace it. Before changing the rc file, it creates a backup.

## Read a command before running it {#read-a-command-before-running-it}

CLI examples use this notation:

| Notation | Meaning | Example |
|---|---|---|
| `<required>` | Replace with a value you must provide. Do not type the angle brackets. | `<output.xyz>` → `model.xyz` |
| `[optional]` | May be omitted; the script then uses its documented behavior. | `[max_corr_steps]` |
| `...` | One or more values may follow. | `<element...>` |

For an unfamiliar Python-backed command, run its `-h` option first to read the argument description. Choose time parameters such as `dt_fs` from your own trajectory-writing settings; they are not inferred from an MD input file. See [Plot Scripts](plot_scripts.md) for plot types and script-specific argument positions.

## Interactive mode {#interactive-mode}

Run `gpumdkit.sh` and select a module by number. Each module provides its own sub-menu.

## CLI mode {#cli-mode}

Direct commands use fixed positional arguments:

```bash
gpumdkit.sh -<option> [args...]
```

To see the complete option list, run:

```bash
gpumdkit.sh -h
```

## Notes {#notes}

- For Python-backed commands that expose option help, use `gpumdkit.sh -<option> -h`; use `gpumdkit.sh -plt -h` to list plot types. Plot-specific argument positions vary by script.
- For detailed module usage, use the task navigation on the [tutorial index](index.md).
