<div align="center">
  <h1>🚀 Quick Start</h1>
  <p style="text-align: justify;">This page helps you install GPUMDkit and run your first commands. You can use GPUMDkit through an interactive menu or direct command-line options.</p>
</div>

## What it does

GPUMDkit provides a single entry point for common tasks in computational materials science — format conversion, structure analysis, property calculation, and visualization — without writing custom scripts.

## Before you start

### 1. Prepare a Python environment

```bash
conda create -n gpumdkit python=3.12
conda activate gpumdkit
```

Some optional functions require additional packages:

```bash
pip install neptrain dpdata calorine
```

Other Python dependencies are loaded by the corresponding scripts. If a package is missing, Python will report it when that function is used.

### 2. Install GPUMDkit

Clone the repository and run the installer:

```bash
git clone https://github.com/zhyan0603/GPUMDkit.git
cd GPUMDkit
source ./install.sh
```

The installer writes `GPUMDkit_path` and `PATH` settings to your shell configuration and loads them for the current shell. It writes to `~/.bashrc` by default; if the current shell is zsh, it writes to `~/.zshrc`. If an existing GPUMDkit path is found, the installer prints the old path and asks whether to replace it. Before changing the rc file, it creates a backup.

Typical installation output:

```text
======================================================
  GPUMDkit Installation
======================================================
 [1/4] Detecting GPUMDkit directory...
       /path/to/GPUMDkit
 [2/4] Detecting shell configuration...
       Target: /Users/you/.bashrc
       Adding environment variables to /Users/you/.bashrc
       Success: Environment variables added.
 [3/4] Setting executable permissions...
       Added executable permission to gpumdkit.sh
 [4/4] Loading environment...

======================================================
  Installation Complete!  GPUMDkit is ready to use.
======================================================
```

If GPUMDkit was installed before, you may see:

```text
Existing GPUMDkit configuration found.
Existing path(s):
  - /old/path/to/GPUMDkit
New path:
  - /new/path/to/GPUMDkit

Replace the existing GPUMDkit configuration with the new path? [y/N]:
```

### 3. Verify the installation

```bash
gpumdkit.sh -h
```

This prints a help table listing all available options.

## Interactive mode

```bash
gpumdkit.sh
```

This opens the main menu:

```text
           ____ ____  _   _ __  __ ____  _    _ _
          / ___|  _ \| | | |  \/  |  _ \| | _(_) |_
         | |  _| |_) | | | | |\/| | | | | |/ / | __|
         | |_| |  __/| |_| | |  | | |_| |   <| | |_
          \____|_|    \___/|_|  |_|____/|_|\_\_|\__|

          GPUMDkit Version 1.5.6 (dev) (2026-06-17)
    Core Developer: Zihan YAN (yanzihan@westlake.edu.cn)
 Main Contributors: Denan LI, Xin WU, Zhoulin LIU & Chen HUA

 ---------------------- GPUMD ------------------------
 1) Format Conversion          2) Sample Structures
 3) Workflow                   4) Calculators
 5) Analyzer                   6) Visualization
 7) Utilities                  8) Developing...
 0) Exit
 ------------>>
 Input the function number:
```

Select a module by number. Each module provides sub-menus with specific functions.

## CLI mode

Direct commands use fixed positional arguments:

```bash
gpumdkit.sh -<option> [args...]
```

Examples:

```bash
gpumdkit.sh -pos2exyz POSCAR model.xyz
gpumdkit.sh -plt train
gpumdkit.sh -calc msd trajectory.xyz Li 10
```

The first example reads `POSCAR` and writes `model.xyz`.

## Common examples

### Convert a POSCAR to extxyz

```bash
gpumdkit.sh -pos2exyz POSCAR model.xyz
```

### Add GPUMD group labels

```bash
gpumdkit.sh -addgroup POSCAR Li Y Cl
```

Group labels are used by some GPUMD-related workflows that need atom grouping, such as species-specific MSD or diffusion calculations.

### Plot NEP training results

```bash
gpumdkit.sh -plt train
```

<div align="center">
  <img src="../../Gallery/train.png" alt="NEP training results" width="72%" />
</div>

### Plot NEP test results

```bash
gpumdkit.sh -plt test
```

<div align="center">
  <img src="../../Gallery/prediction.png" alt="NEP test results" width="72%" />
</div>

### Plot thermodynamic data

```bash
gpumdkit.sh -plt thermo
```

<div align="center">
  <img src="../../Gallery/thermo.png" alt="Thermo plot" width="72%" />
</div>

### Plot MSD and self-diffusion coefficient

```bash
gpumdkit.sh -plt msd
gpumdkit.sh -plt sdc
```

<div align="center">
  <img src="../../Gallery/msd.png" alt="MSD plot" width="45%" />
  <img src="../../Gallery/sdc.png" alt="SDC plot" width="45%" />
</div>

## Notes

- Use `gpumdkit.sh -h` to see all available options.
- Use `gpumdkit.sh -<option> -h` to get help for a specific option (e.g., `gpumdkit.sh -plt train -h`).
- For detailed usage of each module, see the corresponding tutorial pages linked in the [index](index.md).
