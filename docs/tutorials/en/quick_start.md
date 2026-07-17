<div align="center">
  <h1>🚀 Quick Start</h1>
  <p style="text-align: justify;">This page helps you install GPUMDkit and run your first commands. You can use GPUMDkit through an interactive menu or direct command-line options.</p>
</div>

## What it does

GPUMDkit provides a single entry point for common tasks in computational materials science — format conversion, structure analysis, property calculation, and visualization — without writing custom scripts.

## Before you start

### 1. Install with Conda (Recommended)

```bash
conda create -n gpumdkit -c gpumdkit -c conda-forge gpumdkit
conda activate gpumdkit
```

This installs GPUMDkit and its standard dependencies. You do not need to clone the repository or configure environment variables.

Some features require optional packages:

```bash
pip install neptrain calorine
```

### 2. Install from Source (Optional)

To work with the source code, clone the repository and run the installer:

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

For a more useful first check, run:

```bash
gpumdkit.sh -doctor
```

The terminal groups the result into an **Environment** section and **Common** and
**Optional packages** sections. `MISS` for an optional package is not an
installation failure: install it only when you need the feature named beside it.
For example, `dpdata` is needed for DeepMD conversion, while `calorine` is used
by NEP-assisted calculations. This makes it easier to diagnose a missing module
before starting a workflow.

## Read a command before running it

CLI examples use this notation:

| Notation | Meaning | Example |
|---|---|---|
| `<required>` | Replace with a value you must provide. Do not type the angle brackets. | `<output.xyz>` → `model.xyz` |
| `[optional]` | May be omitted; the script then uses its documented behavior. | `[max_corr_steps]` |
| `...` | One or more values may follow. | `<element...>` |

Before using an unfamiliar Python-backed command, ask the command itself for
its argument description. For example:

```text
$ gpumdkit.sh -calc msd -h
 Usage: gpumdkit.sh -calc msd <extxyz_file> <element_symbol> <dt_fs> [max_corr_steps]

 Arguments:
   extxyz_file       Path to the input extxyz trajectory file
   element_symbol    Chemical symbol (e.g., Li, O, Na)
   dt_fs             Time step between consecutive frames (fs)
   max_corr_steps    Max correlation lag steps (optional)
```

Here `dt_fs` is the interval **between stored trajectory frames**, in fs; it is
not automatically inferred from an MD input file. Choose it from your own
trajectory-writing settings. For the plot dispatcher, use `gpumdkit.sh -plt -h`
to list plot types; plot-specific argument positions are documented on the
[Plot Scripts](plot_scripts.md) page.

Interactive menus use `------------>>` to mark where input is expected. After
many CLI commands, the terminal also prints `Code path:`; that path identifies
the implementation responsible for the result or error.

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

          GPUMDkit Version 1.5.6 (dev) (2026-07-10)
    Core Developer: Zihan YAN (yanzihan@westlake.edu.cn)
 Main Contributors: Denan LI, Xin WU, Zhoulin LIU & Chen HUA

 ---------------------- GPUMD ------------------------
 1) Format Conversion          2) Sample Structures
 3) Workflow                   4) Calculators
 5) Analyzer                   6) Visualization
 7) Utilities                  8) Help                
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

## Hello World Example

A minimal end-to-end check that your install works — create a tiny silicon POSCAR, convert it to extxyz, and inspect the output:

```bash
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
```

If `model.xyz` shows a single Si atom with lattice info, your install works.

> **What is extxyz?** extxyz is the extended XYZ format: line 1 is the atom count, line 2 contains the lattice and per-structure properties (energy/forces/virial as needed), and lines 3+ list each atom with its per-atom properties. It is the native training-data format for NEP.

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
- Use `gpumdkit.sh -<option> -h` for Python-backed command help. Use `gpumdkit.sh -plt -h` to list plot types; plot-specific argument positions vary by script.
- For detailed usage of each module, see the corresponding tutorial pages linked in the [index](index.md).
