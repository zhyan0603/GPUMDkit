# GPUMDkit Tutorials

Welcome to **GPUMDkit**!

## What is GPUMDkit?

GPUMDkit is a powerful command-line toolkit that streamlines your tasks for GPUMD and NEP. It provides:

- **🎯 Two Usage Modes**: Interactive mode and command-line mode
- **🔄 Format Conversion**: Convert between VASP, LAMMPS, CP2K, CASTEP, ABACUS, and extxyz formats.
- **📊 Visualization**: Comprehensive plotting tools for NEP training, MD simulations, and analysis.
- **🧮 Calculators**: Compute ionic conductivity, descriptors, etc.
- **🔍 Analysis Tools**: Structure validation, filtering, composition analysis, and quality checks, etc.
- **⚙️ Workflow Automation**: Batch processing for DFT calculations and MD simulations.
- **🤖 Active Learning**: Semi-Automated NEP model training workflow.

## Quick Start

### Installation

```bash
# 1. Clone the repository
git clone https://github.com/zhyan0603/GPUMDkit.git

# 2. Add to your ~/.bashrc
export GPUMDkit_path=/path/to/GPUMDkit
export PATH=${GPUMDkit_path}:${PATH}
source ${GPUMDkit_path}/Scripts/utils/completion.sh

# 3. Reload your shell
source ~/.bashrc

# 4. Make executable
cd ${GPUMDkit_path}
chmod +x gpumdkit.sh
```

## Two Ways to Use GPUMDkit

### 🖱️ Interactive Mode (Beginner-Friendly)

Launch the interactive menu:
```bash
gpumdkit.sh
```

You'll see:
```
           ____ ____  _   _ __  __ ____  _    _ _
          / ___|  _ \| | | |  \/  |  _ \| | _(_) |_
         | |  _| |_) | | | | |\/| | | | | |/ / | __|
         | |_| |  __/| |_| | |  | | |_| |   <| | |_
          \____|_|    \___/|_|  |_|____/|_|\_\_|\__|

          GPUMDkit Version 1.5.5 (dev) (2026-05-10)
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

**Best for:**

- Learning GPUMDkit features
- Complex multi-step tasks
- When you need guided prompts

### ⚡ Command-Line Mode

Direct command execution:
```bash
gpumdkit.sh -h
```

You will see:

```
+-------------------------------------------------------------------------------------------------------+
|                          GPUMDkit 1.5.5 (dev) (2026-05-10) Command Help                               |
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
| -pda          Probability density analysis    | -hbond             Hydrogen-bond analysis             |
| -pynep        FPS sampling by PyNEP           |                                                       |
+-------------------------------------------------------------------------------------------------------+
| Detailed usage: gpumdkit.sh -<option> -h    Plot details: gpumdkit.sh -plt <type> -h                  |
+-------------------------------------------------------------------------------------------------------+
```

**Best for:**

- Quick one-off operations
- Scripting and automation
- Batch processing

## Some Functionalities

### 1. Format Conversion

Convert structure files between different formats:

#### Interactive mode (option 1)

```shell
 +-------------------------------------------------------------+
 |                   FORMAT CONVERSION TOOLS                   |
 +-------------------------------------------------------------+
 | 101) VASP to extxyz            106) Add group labels        |
 | 102) MTP to extxyz             107) Add weight to extxyz    |
 | 103) CP2K to extxyz            108) Extract frame extxyz    |
 | 104) ABACUS to extxyz          109) Clean XYZ info          |
 | 105) extxyz to POSCAR          110) Replicate structure     |
 +-------------------------------------------------------------+
 | out2exyz) OUTCAR to extxyz     xdat2exyz) XDATCAR to extxyz |
 | pos2exyz) POSCAR to extxyz     pos2lmp)   POSCAR to LAMMPS  |
 | cif2pos)  CIF to POSCAR        lmp2exyz)  LAMMPS to extxyz  |
 | cif2exyz) CIF to extxyz        traj2exyz) ASE traj to extxyz|
 +-------------------------------------------------------------+
 | 000) Return to main menu                                    |
 +-------------------------------------------------------------+
 Input the function number or converter keyword:
```

#### Command-line mode

```shell
# VASP → extxyz
gpumdkit.sh -out2xyz <dir>

# POSCAR → extxyz
gpumdkit.sh -pos2exyz POSCAR model.xyz

# CIF → extxyz
gpumdkit.sh -cif2exyz input.cif model.xyz

# LAMMPS → extxyz
gpumdkit.sh -lmp2exyz dump.lammpstrj Li Y Cl

# Add group labels
gpumdkit.sh -addgroup POSCAR Li Y Cl

# Add weights to structures
gpumdkit.sh -addweight train.xyz weighted.xyz 2.0
```

**📖 Learn more:** [Format Conversion Guide](format_conversion.md)

### 2. Structure Sampling

Select diverse structures for training:

#### Interactive mode (option 2)

```shell
 +------------------------------------------------------+
 |                 SAMPLE STRUCTURE TOOLS               |
 +------------------------------------------------------+
 | 201) Sample structures from extxyz                   |
 | 202) FPS sampling by PyNEP [deprecated]              |
 | 203) FPS sampling by NepTrain [preferred]            |
 | 204) Perturb structure                               |
 | 205) Select max force deviation structs              |
 +------------------------------------------------------+
 | 000) Return to the main menu                         |
 +------------------------------------------------------+
 Input the function number:
```

**📖 Learn more:** [Structure Sampling Guide](sample_structures.md)

### 3. Workflow

Batch processing for high-throughput calculations:

#### Interactive mode (option 3)

```
 +---------------------------------------------------------+
 |                      WORKFLOW TOOLS                     |
 +---------------------------------------------------------+
 | 301) SCF batch pretreatment                             |
 | 302) MD sample batch pretreatment (gpumd)               |
 | 303) MD sample batch pretreatment (lmp)                 |
 +---------------------------------------------------------+
 | 000) Return to the main menu                            |
 +---------------------------------------------------------+
 Input the function number:
```

**📖 Learn more:** [Workflow Guide](workflow.md) | [Active Learning](workflow_active_learning.md)

### 4. Calculators

Compute material properties:

#### Interactive mode (option 4)

```shell
 +----------------------------------------------------------+
 |                     CALCULATOR TOOLS                     |
 +----------------------------------------------------------+
 | 401) Calc ionic conductivity                             |
 | 402) Calc properties by nep                              |
 | 403) Calc descriptors of specific elements               |
 | 404) Calc density of atomistic states (DOAS)             |
 | 405) Calc nudged elastic band (NEB) by nep               |
 | 406) Build neighbor list                                 |
 | 407) Calc displacement from trajectory                   |
 | 408) Calc averaged structure                             |
 | 409) Calc octahedral tilt                                |
 | 410) Calc polarization for ABO3                          |
 | 411) Minimize structure by nep                           |
 | 412) Calc mean square displacement (MSD) from trajectory |
 +----------------------------------------------------------+
 | 000) Return to the main menu                             |
 +----------------------------------------------------------+
 Input the function number:
```

#### Command-line mode

```bash
# Ionic conductivity from MSD
gpumdkit.sh -calc ionic-cond Li 1

# Calculate descriptors
gpumdkit.sh -calc des train.xyz desc.npy nep.txt Li

# Density of atomistic states
gpumdkit.sh -calc doas structures.xyz nep.txt doas.txt
```

**📖 Learn more:** [Calculators Guide](calculators.md)

### 5. Structure Analysis

Validate and filter structure files:

#### Interactive mode (option 5)

```
 +------------------------------------------------------+
 |                    ANALYZER TOOLS                    |
 +------------------------------------------------------+
 | 501) Analyze composition of extxyz                   |
 | 502) Find outliers of extxyz                         |
 | 503) Analyze chemical species of extxyz              |
 | 504) Check charge balance of extxyz                  |
 | 505) Analyze energy/force/virial range               |
 | 506) Filter structures by minimum distance           |
 | 507) Get minimum interatomic distance                |
 | 508) Probability density analysis                    |
 +------------------------------------------------------+
 | 000) Return to the main menu                         |
 +------------------------------------------------------+
 Input the function number:
```

#### Command-line mode

```bash
# Check data quality
gpumdkit.sh -range train.xyz force
gpumdkit.sh -min_dist_pbc dump.xyz
gpumdkit.sh -analyze_comp train.xyz
```

**📖 Learn more:** [Analyzer Guide](analyzer.md)

### 6. Visualization & Plotting

Visualize simulation results and training data:

```bash
# NEP Training
gpumdkit.sh -plt train           # Training plot
gpumdkit.sh -plt prediction      # Prediction plot
gpumdkit.sh -plt force_errors    # Force error analysis

# MD Simulations
gpumdkit.sh -plt thermo          # Thermo info
gpumdkit.sh -plt msd             # Mean square displacement
gpumdkit.sh -plt rdf             # Radial distribution function

# Transport Properties
gpumdkit.sh -plt sdc             # Self-diffusion coefficient
gpumdkit.sh -plt emd             # Thermal conductivity
```

**📖 Learn more:** [Plot Scripts Guide](plot_scripts.md)

## Detailed Topic Guides

For in-depth information on specific topics:

- 📄 [Format Conversion](format_conversion.md) - Converting between file formats
- 📊 [Plot Scripts](plot_scripts.md) - Visualization and plotting
- 🧮 [Calculators](calculators.md) - Property calculations
- 🔍 [Analyzer](analyzer.md) - Data analysis and quality control
- 📦 [Structure Sampling](sample_structures.md) - Sampling methods and strategies
- ⚙️ [Workflow Automation](workflow.md) - Batch processing
- 🤖 [Active Learning](workflow_active_learning.md) - NEP model improvement
- 🔧 [Custom Commands](custom_commands.md) - Extending GPUMDkit

## Advanced Features

### Custom Commands

Extend GPUMDkit with your own shortcuts by editing `~/.gpumdkit.in`:

```bash
custom_quick_analyze() {
    gpumdkit.sh -range "$1" force
}
```

Then use: `gpumdkit.sh -quick_analyze train.xyz`

**Learn more:** [Custom Commands Guide](custom_commands.md)

### Scripting with GPUMDkit

GPUMDkit commands work great in shell scripts:

```bash
#!/bin/bash
for temp in 600 900 1200; do
    cd temp_${temp}/
    gpumdkit.sh -plt thermo save
    gpumdkit.sh -plt msd save
    (echo 201; echo "dump.xyz random 50") | gpumdkit.sh > /dev/null
    cd ..
done
```

---

**🎉 Ready to get started?** Try the interactive mode: `gpumdkit.sh`

Thank you for using GPUMDkit! For questions or feedback, please contact Zihan YAN (yanzihan@westlake.edu.cn) or open an issue on GitHub.
