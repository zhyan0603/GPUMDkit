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

        GPUMDkit Version 1.4.2 (dev) (2025-12-17)
  Core Developer: Zihan YAN (yanzihan@westlake.edu.cn)

 ----------------------- GPUMD -----------------------
 1) Format Conversion          2) Sample Structures
 3) Workflow                   4) Calculators
 5) Analyzer                   6) Developing ...
 0) Quit!
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
+==================================================================================================+
|                              GPUMDkit 1.4.2 (dev) (2025-12-17) Usage                             |
+======================================== Conversions =============================================+
| -out2xyz       Convert OUTCAR to extxyz       | -pos2exyz     Convert POSCAR to extxyz           |
| -cif2pos       Convert cif to POSCAR          | -pos2lmp      Convert POSCAR to LAMMPS           |
| -cif2exyz      Convert cif to extxyz          | -lmp2exyz     Convert LAMMPS-dump to extxyz      |
| -addgroup      Add group label                | -addweight    Add weight to the struct in extxyz |
| Developing...                                 | Developing...                                    |
+========================================= Analysis ===============================================+
| -range         Print range of energy etc.     | -max_rmse     Get max RMSE from extxyz           |
| -min_dist      Get min_dist between atoms     | -min_dist_pbc Get min_dist considering PBC       |
| -filter_box    Filter struct by box limits    | -filter_value Filter struct by value (efs)       |
| -filter_dist   Filter struct by min_dist      | -analyze_comp Analyze composition of extxyz      |
| -pynep         Sample struct by pynep         | Developing...                                    |
+====================================== Misc Utilities ============================================+
| -plt           Plot scripts                   | -get_frame     Extract the specified frame       |
| -calc          Calculators                    | -clean_xyz     Clean extra info in XYZ file      |
| -clean         Clear files for work_dir       | -time          Time consuming Analyzer           |
| -update        Update GPUMDkit                | Developing...                                    |
+==================================================================================================+
| For detailed usage and examples, use: gpumdkit.sh -<option> -h                                   |
+==================================================================================================+
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
 ------------>>
 101) Convert VASP to extxyz
 102) Convert mtp to extxyz
 103) Convert CP2K to extxyz
 104) Convert ABACUS to extxyz
 105) Convert extxyz to POSCAR
 000) Return to the main menu
 ------------>>
 Input the function number:
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
 ------------>>
 201) Sample structures from extxyz
 202) Sample structures by pynep
 203) Sample structures by neptrain
 204) Perturb structure
 205) Select max force deviation structs
 000) Return to the main menu
 ------------>>
 Input the function number:
```

**📖 Learn more:** [Structure Sampling Guide](sample_structures.md)

### 3. Workflow

Batch processing for high-throughput calculations:

#### Interactive mode (option 3)

```
 ------------>>
 301) SCF batch pretreatment
 302) MD sample batch pretreatment (gpumd)
 303) MD sample batch pretreatment (lmp)
 000) Return to the main menu
 ------------>>
 Input the function number:
```

**📖 Learn more:** [Workflow Guide](workflow.md) | [Active Learning](workflow_active_learning.md)

### 4. Calculators

Compute material properties:

#### Interactive mode (option 4)

```shell
 ------------>>
 401) Calc ionic conductivity
 402) Calc properties by nep
 403) Calc descriptors of specific elements
 404) Calc density of atomistic states (DOAS)
 405) Calc nudged elastic band (NEB) by nep
 000) Return to the main menu
 ------------>>
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
 ------------>>
 501) Analyze composition of extxyz
 502) Find outliers of extxyz
 000) Return to the main menu
 ------------>>
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
