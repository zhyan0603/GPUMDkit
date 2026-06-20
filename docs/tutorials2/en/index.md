# GPUMDkit Tutorials

<div align="center">
  <p>
    <a href="../zh/index.md">中文</a> | <strong>English</strong>
  </p>
</div>

Welcome to **GPUMDkit** - a powerful command-line toolkit for GPUMD and NEP!

## What is GPUMDkit?

GPUMDkit streamlines your molecular dynamics workflow with:

- **Format Conversion**: Convert between VASP, LAMMPS, CP2K, ABACUS, CIF, and extxyz formats
- **Calculators**: Compute ionic conductivity, descriptors, MSD, NEB, and more
- **Analyzers**: Structure validation, filtering, composition analysis
- **Visualization**: 31+ plot types for training, transport, and structural analysis
- **Workflows**: Batch processing for DFT and MD simulations
- **Sampling**: Structure selection using uniform, random, or FPS methods

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/zhyan0603/GPUMDkit.git

# Install
cd GPUMDkit
source ./install.sh
```

### Dependencies

```bash
# Create conda environment
conda create -n gpumdkit python=3.12
conda activate gpumdkit

# Install required packages
pip install neptrain ase pymatgen dpdata numpy scipy matplotlib
```

### Basic Usage

```bash
# Interactive mode (menu-driven)
gpumdkit.sh

# Command-line mode
gpumdkit.sh -<option> [args...]

# Get help
gpumdkit.sh -h                    # General help
gpumdkit.sh -<option> -h          # Specific command help
gpumdkit.sh -plt -h               # Plot help
gpumdkit.sh -calc -h              # Calculator help
```

## Tutorials

| Tutorial | Description |
|----------|-------------|
| [Quick Start](quickstart.md) | Installation, configuration, and first steps |
| [Format Conversion](format_conversion.md) | Convert between structure file formats |
| [Calculators](calculators.md) | Compute material properties |
| [Analyzers](analyzers.md) | Structure analysis and validation |
| [Visualization](visualization.md) | Plotting and data visualization |
| [Workflows](workflows.md) | Batch processing and automation |
| [Structure Sampling](sampling.md) | Structure selection methods |
| [NEP Training Guide](nep_training.md) | Complete NEP model training workflow |

## Common Workflows

### Workflow 1: NEP Model Training

```bash
# 1. Convert DFT data
gpumdkit.sh -out2xyz ./vasp_results/

# 2. Add group labels
gpumdkit.sh -addgroup POSCAR Li Y Cl

# 3. Sample structures
gpumdkit.sh  # Select: 2) Sample Structures -> 203

# 4. Train NEP (external)
# ...

# 5. Validate
gpumdkit.sh -plt train
gpumdkit.sh -plt prediction
```

### Workflow 2: Ionic Conductivity

```bash
# 1. Calculate MSD
gpumdkit.sh -calc msd trajectory.xyz Li 10

# 2. Calculate conductivity
gpumdkit.sh -calc ionic-cond Li 1

# 3. Visualize
gpumdkit.sh -plt msd
gpumdkit.sh -plt sdc
```

### Workflow 3: Structure Analysis

```bash
# 1. Check composition
gpumdkit.sh -analyze_comp train.xyz

# 2. Check distances
gpumdkit.sh -min_dist_pbc train.xyz

# 3. Check force range
gpumdkit.sh -range train.xyz force
```

## Module Overview

| Module | Menu | CLI Flag | Description |
|--------|------|----------|-------------|
| Format Conversion | 1 | `-out2xyz`, `-pos2exyz`, etc. | Structure file conversion |
| Sample Structures | 2 | `-pynep` | Structure sampling |
| Workflows | 3 | - | Batch processing |
| Calculators | 4 | `-calc <type>` | Property calculations |
| Analyzers | 5 | `-range`, `-min_dist`, etc. | Structure analysis |
| Visualization | 6 | `-plt <type>` | Plotting tools |
| Utilities | 7 | `-time`, `-clean` | Utility functions |

## Output File Conventions

| File | Content |
|------|---------|
| `model.xyz` | Structure file in extxyz format |
| `train.xyz` | Training dataset |
| `nep.txt` | NEP model file |
| `msd.out` | Mean square displacement data |
| `thermo.out` | Thermodynamic properties |
| `loss.out` | Training loss history |
| `rdf.out` | Radial distribution function |
| `sdc.out` | Self-diffusion coefficient data |

## Getting Help

- **General help**: `gpumdkit.sh -h`
- **Plot help**: `gpumdkit.sh -plt -h`
- **Calculator help**: `gpumdkit.sh -calc -h`
- **Specific command**: `gpumdkit.sh -<option> -h`

## Contact

- **Developer**: Zihan YAN (yanzihan@westlake.edu.cn)
- **GitHub**: https://github.com/zhyan0603/GPUMDkit
- **Documentation**: https://zhyan0603.github.io/GPUMDkit/
