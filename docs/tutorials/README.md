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

## Core Functionalities

### 1. Format Conversion

Convert structure files between different formats:

```bash
# VASP → extxyz
gpumdkit.sh -out2xyz ./vasp_calculations/

# POSCAR → extxyz
gpumdkit.sh -pos2exyz POSCAR model.xyz

# CIF → extxyz
gpumdkit.sh -cif2exyz structure.cif

# LAMMPS → extxyz
gpumdkit.sh -lmp2exyz dump.lammpstrj Li Y Cl

# Add group labels (required for NEP training)
gpumdkit.sh -addgroup POSCAR Li Y Cl

# Add weights to structures
gpumdkit.sh -addweight train.xyz weighted.xyz 2.0
```

**📖 Learn more:** [Format Conversion Guide](format_conversion.md)

### 2. Visualization & Plotting

Visualize simulation results and training data:

```bash
# NEP Training
gpumdkit.sh -plt train           # Training progress
gpumdkit.sh -plt prediction      # Prediction accuracy
gpumdkit.sh -plt force_errors    # Force error analysis

# MD Simulations
gpumdkit.sh -plt thermo          # Temperature, pressure, energy
gpumdkit.sh -plt msd             # Mean square displacement
gpumdkit.sh -plt rdf             # Radial distribution function

# Transport Properties
gpumdkit.sh -plt sdc             # Self-diffusion coefficient
gpumdkit.sh -plt emd             # Thermal conductivity

# Advanced Analysis
gpumdkit.sh -plt des umap desc.npy  # Descriptor visualization
gpumdkit.sh -plt charge          # Charge distribution
```

**📖 Learn more:** [Plot Scripts Guide](plot_scripts.md)

### 3. Calculators

Compute material properties:

```bash
# Ionic conductivity from MSD
gpumdkit.sh -calc ionic-cond Li 1

# NEP predictions
gpumdkit.sh -calc nep input.xyz output.xyz nep.txt

# Calculate descriptors
gpumdkit.sh -calc des umap train.xyz desc.npy nep.txt Li

# Density of atomistic states
gpumdkit.sh -calc doas structures.xyz nep.txt output.txt
```

**📖 Learn more:** [Calculators Guide](calculators.md)

### 4. Structure Analysis

Validate and filter structure files:

```bash
# Check data quality
gpumdkit.sh -range train.xyz force
gpumdkit.sh -min_dist_pbc model.xyz
gpumdkit.sh -analyze_comp train.xyz

# Filter outliers
gpumdkit.sh -filter_value train.xyz force 30
gpumdkit.sh -filter_dist train.xyz 1.5
gpumdkit.sh -filter_box train.xyz 20
```

**📖 Learn more:** [Analyzer Tools Guide](analyzer.md)

### 5. Structure Sampling

Select diverse structures for training:

```bash
# PyNEP farthest point sampling
gpumdkit.sh -pynep candidates.xyz selected.xyz nep.txt

# Generate perturbed structures
# (via interactive mode option 204)
```

**📖 Learn more:** [Structure Sampling Guide](sample_structures.md)

### 6. Workflow Automation

Batch processing for high-throughput calculations:

```bash
# Interactive mode
gpumdkit.sh
# Select: 3) Workflow
# - 301) SCF batch pretreatment (VASP/CP2K)
# - 302) MD batch pretreatment (GPUMD)
# - 303) MD batch pretreatment (LAMMPS)
```

**📖 Learn more:** [Workflow Guide](workflow_dev.md) | [Active Learning](workflow_active_learning.md)

## Complete Feature Reference

### Format Conversion Commands

| Command | Description | Example |
|---------|-------------|---------|
| `-out2xyz` | OUTCAR → extxyz | `gpumdkit.sh -out2xyz ./` |
| `-pos2exyz` | POSCAR → extxyz | `gpumdkit.sh -pos2exyz POSCAR model.xyz` |
| `-exyz2pos` | extxyz → POSCAR | `gpumdkit.sh -exyz2pos train.xyz` |
| `-cp2k2exyz` | CP2K → extxyz | `gpumdkit.sh -cp2k2exyz ./cp2k_out/` |
| `-abacus2exyz` | ABACUS → extxyz | `gpumdkit.sh -abacus2exyz ./` |
| `-cif2exyz` | CIF → extxyz | `gpumdkit.sh -cif2exyz structure.cif` |
| `-cif2pos` | CIF → POSCAR | `gpumdkit.sh -cif2pos structure.cif` |
| `-lmp2exyz` | LAMMPS → extxyz | `gpumdkit.sh -lmp2exyz dump.lmp Li Y` |
| `-pos2lmp` | POSCAR → LAMMPS | `gpumdkit.sh -pos2lmp POSCAR lmp.data Li Y` |
| `-addgroup` | Add group labels | `gpumdkit.sh -addgroup POSCAR Li Y Cl` |
| `-addweight` | Add structure weight | `gpumdkit.sh -addweight in.xyz out.xyz 2.0` |
| `-get_frame` | Extract frame | `gpumdkit.sh -get_frame dump.xyz 1000` |

### Analysis Commands

| Command | Description | Example |
|---------|-------------|---------|
| `-range` | Show property range | `gpumdkit.sh -range train.xyz force` |
| `-min_dist` | Min atomic distance | `gpumdkit.sh -min_dist model.xyz` |
| `-min_dist_pbc` | Min distance with PBC | `gpumdkit.sh -min_dist_pbc model.xyz` |
| `-filter_value` | Filter by threshold | `gpumdkit.sh -filter_value train.xyz force 30` |
| `-filter_dist` | Filter by distance | `gpumdkit.sh -filter_dist train.xyz 1.5` |
| `-filter_box` | Filter by box size | `gpumdkit.sh -filter_box train.xyz 20` |
| `-analyze_comp` | Composition analysis | `gpumdkit.sh -analyze_comp train.xyz` |
| `-time` | Estimate runtime | `gpumdkit.sh -time gpumd` or `-time nep` |

### Plotting Commands

| Command | Description | Input File |
|---------|-------------|------------|
| `-plt thermo` | Thermodynamic properties | `thermo.out` |
| `-plt train` | NEP training progress | `loss.out`, `*_train.out` |
| `-plt prediction` | NEP predictions | `*_test.out` |
| `-plt train_test` | Train & test comparison | `*_train.out`, `*_test.out` |
| `-plt msd` | Mean square displacement | `msd.out` |
| `-plt msd_all` | MSD per species | `msd.out` |
| `-plt sdc` | Self-diffusion coefficient | `msd.out` |
| `-plt rdf` | Radial distribution | `rdf.out` |
| `-plt vac` | Velocity autocorrelation | `sdc.out` |
| `-plt force_errors` | Force error metrics | `force_train.out` |
| `-plt charge` | Charge distribution | `charge_train.out` |
| `-plt lr` | Learning rate | `loss.out` |
| `-plt restart` | Restart parameters | `nep.restart` |
| `-plt des` | Descriptors | `descriptors.npy` |
| `-plt dimer` | Dimer interaction | NEP model |
| `-plt doas` | Density of states | DOAS files |
| `-plt emd` | EMD thermal conductivity | EMD outputs |
| `-plt nemd` | NEMD thermal transport | NEMD outputs |

### Calculator Commands

| Command | Description |
|---------|-------------|
| `-calc ionic-cond` | Ionic conductivity |
| `-calc nep` | NEP predictions |
| `-calc des` | Calculate descriptors |
| `-calc doas` | Density of atomistic states |

### Utility Commands

| Command | Description |
|---------|-------------|
| `-clean` | Clean working directory |
| `-update` | Update GPUMDkit |
| `-h` | Show help |

## Real-World Examples

### Example 1: Prepare NEP Training Data from VASP

```bash
# 1. Convert VASP calculations to extxyz
cd vasp_calculations/
gpumdkit.sh -out2xyz ./

# 2. Check data quality
gpumdkit.sh -range train.xyz force
gpumdkit.sh -min_dist_pbc train.xyz

# 3. Filter problematic structures
gpumdkit.sh -filter_value train.xyz force 30

# 4. Add group labels for NEP
gpumdkit.sh -addgroup POSCAR Li Y Cl

# 5. Now train.xyz and model.xyz are ready for NEP!
```

### Example 2: Monitor and Visualize NEP Training

```bash
# During training, in another terminal
gpumdkit.sh -plt train

# After training completes
gpumdkit.sh -plt train
gpumdkit.sh -plt prediction
gpumdkit.sh -plt force_errors
gpumdkit.sh -plt lr

# Check descriptor coverage
gpumdkit.sh -calc des umap train.xyz desc.npy nep.txt Li
gpumdkit.sh -plt des umap desc.npy
```

### Example 3: Analyze MD Simulation

```bash
# After GPUMD simulation
cd simulation_results/

# Thermodynamic properties
gpumdkit.sh -plt thermo

# Transport properties
gpumdkit.sh -plt msd
gpumdkit.sh -calc ionic-cond Li 1

# Structure
gpumdkit.sh -plt rdf
```

### Example 4: Active Learning Workflow

```bash
# Iteration 1: Initial training
cd iteration_01/
nep  # Train initial model

# Iteration 2: Active learning
cd ../iteration_02/

# Run MD with active mode, then select uncertain structures
gpumdkit.sh -pynep md_trajectory.xyz selected.xyz ../iteration_01/nep.txt

# Set up DFT calculations
gpumdkit.sh  # Interactive mode → 3) Workflow → 301

# After DFT, add to training set and retrain
cat ../iteration_01/train.xyz dft_results.xyz > train.xyz
nep  # Retrain
```

## Tips & Best Practices

💡 **Use Tab Completion**: Bash completion makes command entry faster
```bash
gpumdkit.sh -plt <Tab>  # Shows available plot types
```

💡 **Save Plots**: Add `save` to any plot command to export PNG
```bash
gpumdkit.sh -plt thermo save
```

💡 **Check Help**: Most commands have help options
```bash
gpumdkit.sh -h
gpumdkit.sh -plt -h
gpumdkit.sh -calc -h
```

💡 **Organize Your Data**: Keep training data in well-structured directories
```
project/
├── 00_raw_data/
├── 01_converted/
├── 02_filtered/
├── 03_training/
└── 04_models/
```

💡 **Quality Control**: Always validate data before training
```bash
gpumdkit.sh -range train.xyz force
gpumdkit.sh -min_dist_pbc train.xyz
gpumdkit.sh -analyze_comp train.xyz
```

## Getting Help

- **Documentation**: This guide and linked tutorials
- **Script Documentation**: See [Scripts directory README files](../../Scripts/)
- **GitHub Issues**: [Report bugs or request features](https://github.com/zhyan0603/GPUMDkit/issues)
- **Contact**: Zihan YAN (yanzihan@westlake.edu.cn)

## Detailed Topic Guides

For in-depth information on specific topics:

- 📄 [Format Conversion](format_conversion.md) - Converting between file formats
- 📊 [Plot Scripts](plot_scripts.md) - Visualization and plotting
- 🧮 [Calculators](calculators.md) - Property calculations
- 🔍 [Analyzer Tools](analyzer.md) - Data analysis and quality control
- 📦 [Structure Sampling](sample_structures.md) - Sampling methods and strategies
- ⚙️ [Workflow Automation](workflow_dev.md) - Batch processing
- 🤖 [Active Learning](workflow_active_learning.md) - NEP model improvement
- 🔧 [Custom Commands](custom_commands.md) - Extending GPUMDkit

## Advanced Features

### Custom Commands

Extend GPUMDkit with your own shortcuts by editing `~/.gpumdkit.in`:

```bash
custom_quick_analyze() {
    gpumdkit.sh -range "$1" force
    gpumdkit.sh -min_dist_pbc "$1"
    gpumdkit.sh -analyze_comp "$1"
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
    gpumdkit.sh -plt thermo
    gpumdkit.sh -plt msd
    gpumdkit.sh -calc ionic-cond Li 1
    cd ..
done
```

---

**🎉 Ready to get started?** Try the interactive mode: `gpumdkit.sh`

Thank you for using GPUMDkit! For questions or feedback, please contact Zihan YAN (yanzihan@westlake.edu.cn) or open an issue on GitHub.
