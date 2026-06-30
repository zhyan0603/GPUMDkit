<div align="center">
  <h1>🚀 Quick Start</h1>
  <p style="text-align: justify;">This page introduces the basic installation and several common commands. GPUMDkit can be used through the interactive menu or by direct command-line options.</p>
</div>

## 1. Prepare a Python Environment

You can use a clean conda environment:

```bash
conda create -n gpumdkit python=3.12
conda activate gpumdkit
```

Some optional functions use additional packages:

```bash
pip install neptrain dpdata calorine
```

Other Python dependencies are loaded by the corresponding scripts. If a package is missing, Python will report it when that function is used.

## 2. Install GPUMDkit

Clone the repository:

```bash
git clone https://github.com/zhyan0603/GPUMDkit.git
cd GPUMDkit
```

Run the installer:

```bash
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

Check the command:

```bash
gpumdkit.sh -h
```

The help table:

```text
+-------------------------------------------------------------------------------------------------------+
|                          GPUMDkit 1.5.6 (dev) (2026-06-17) Command Help                               |
+-------------------------------------------------------------------------------------------------------+
|                                          MAIN FUNCTIONS                                               |
+-------------------------------------------------------------------------------------------------------+
| -h            Show this help table            | -plt <type>        Plot and visualization tools       |
| -calc <type>  Calculator tools                | -time <gpumd|nep>  Time-consuming analyzer            |
| -update       Update GPUMDkit                 | -clean             Clean extra files in current dir   |
+-------------------------------------------------------------------------------------------------------+
|                                         FORMAT CONVERSION                                             |
+-------------------------------------------------------------------------------------------------------+
| -out2xyz      OUTCAR -> extxyz (shell)        | -out2exyz          OUTCAR -> extxyz (python)          |
| -cp2k2xyz     CP2K log -> xyz                 | -xdat2exyz         XDATCAR -> extxyz                  |
| -cif2pos      cif -> POSCAR                   | -cif2exyz          cif -> extxyz                      |
| -pos2exyz     POSCAR -> extxyz                | -exyz2pos          extxyz -> POSCAR                   |
| -pos2lmp      POSCAR -> LAMMPS data           | -lmp2exyz          LAMMPS dump -> extxyz              |
| -traj2exyz    ASE traj -> extxyz              | -replicate         Replicate structure                |
| -addgroup     Add group labels                | -addweight         Add structure weight in extxyz     |
| -clean_xyz    Clean extra info in extxyz      | -get_frame         Extract specific frame             |
| -frame_range  Extract frames by range         |                                                       |
+-------------------------------------------------------------------------------------------------------+
|                                            ANALYSIS                                                   |
+-------------------------------------------------------------------------------------------------------+
| -range        Energy/force/virial statistics  | -analyze_comp      Analyze composition                |
| -chem_species Analyze chemical species        | -cbc               Charge balance check               |
| -min_dist     Min distance (no PBC)           | -min_dist_pbc      Min distance with PBC              |
| -filter_dist  Filter by min_dist (no PBC)     | -filter_dist_pbc   Filter by min_dist (PBC)           |
| -pda          Probability density analysis    | -filter_box        Filter by box-edge length          |
| -pynep        Deprecated PyNEP sampling       | -nep_modifier      Modify NEP model interactively     |
+-------------------------------------------------------------------------------------------------------+
| Detailed usage: gpumdkit.sh -<option> -h    Plot details: gpumdkit.sh -plt <type> -h                  |
+-------------------------------------------------------------------------------------------------------+
```

## 3. Usage Modes

### Interactive Mode

```bash
gpumdkit.sh
```

Main menu:

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

### Command-Line Mode

Direct commands use fixed positional arguments:

```bash
gpumdkit.sh -pos2exyz POSCAR model.xyz
gpumdkit.sh -plt train
gpumdkit.sh -plt thermo
```

For example:

```bash
gpumdkit.sh -pos2exyz POSCAR model.xyz
```

means reading `POSCAR` and writing `model.xyz`.

## 4. Practical Examples

### Convert POSCAR to extxyz

```bash
gpumdkit.sh -pos2exyz POSCAR model.xyz
```

### Add GPUMD Group Labels

```bash
gpumdkit.sh -addgroup POSCAR Li Y Cl
```

Group labels are used by some GPUMD-related workflows and analyses that need atom grouping, such as species-specific MSD or diffusion calculations.

### Plot NEP Training Results

```bash
gpumdkit.sh -plt train
```

<div align="center">
  <img src="../../Gallery/train.png" alt="NEP training results" width="72%" />
</div>

### Plot NEP Test Results

```bash
gpumdkit.sh -plt test
```

<div align="center">
  <img src="../../Gallery/prediction.png" alt="NEP test results" width="72%" />
</div>

### Plot Thermodynamic Data

```bash
gpumdkit.sh -plt thermo
```

<div align="center">
  <img src="../../Gallery/thermo.png" alt="Thermo plot" width="72%" />
</div>

### Plot MSD and Self-Diffusion Coefficient

```bash
gpumdkit.sh -plt msd
gpumdkit.sh -plt sdc
```

<div align="center">
  <img src="../../Gallery/msd.png" alt="MSD plot" width="45%" />
  <img src="../../Gallery/sdc.png" alt="SDC plot" width="45%" />
</div>

## 5. Help

```bash
gpumdkit.sh -h
gpumdkit.sh -plt -h
gpumdkit.sh -calc -h
gpumdkit.sh -plt train -h
gpumdkit.sh -calc msd -h
```
