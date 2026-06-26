---
name: gpumdkit
description: >
  Use when working with GPUMD molecular dynamics or NEP potentials for computational materials science.
  Provides format conversion, structure analysis, visualization, calculation tools, and workflow automation.
  Supports VASP, LAMMPS, CP2K, ABACUS, CIF, and extxyz formats.
  Use when user asks about: molecular dynamics, NEP training, structure conversion, ionic conductivity,
  radial distribution function, mean square displacement, thermal conductivity, phonon DOS,
  structure sampling, farthest point sampling, or computational materials analysis.
allowed-tools: Bash(gpumdkit *) Bash(python3 *) Bash(python *)
---

# GPUMDkit - Molecular Dynamics Toolkit

## Overview

GPUMDkit is a command-line toolkit for GPUMD (Graphics Processing Units Molecular Dynamics) and
NEP (Neuroevolution Potential) programs. It streamlines common tasks in computational materials science.

## Quick Start

### Installation
```bash
git clone https://github.com/zhyan0603/GPUMDkit.git
cd GPUMDkit && source ./install.sh
```

### Dependencies
```bash
conda create -n gpumdkit python=3.12
conda activate gpumdkit
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

## Module Overview

| Module | Description | Skill Reference |
|--------|-------------|-----------------|
| Format Conversion | Convert between VASP, LAMMPS, CP2K, ABACUS, CIF, extxyz | See [gpumdkit-format-conversion](../gpumdkit-format-conversion/SKILL.md) |
| Calculators | Compute ionic conductivity, descriptors, MSD, NEB, etc. | See [gpumdkit-calculators](../gpumdkit-calculators/SKILL.md) |
| Analyzers | Structure validation, filtering, composition analysis | See [gpumdkit-analyzers](../gpumdkit-analyzers/SKILL.md) |
| Visualization | Plot training results, transport properties, structural data | See [gpumdkit-visualization](../gpumdkit-visualization/SKILL.md) |
| Workflows | Batch processing for DFT and MD simulations | See [gpumdkit-workflows](../gpumdkit-workflows/SKILL.md) |
| Sampling | Structure selection using uniform, random, or FPS methods | See [gpumdkit-sampling](../gpumdkit-sampling/SKILL.md) |

## Common Workflows

### Workflow 1: NEP Model Training Pipeline
```bash
# 1. Convert DFT data to extxyz
gpumdkit.sh -out2xyz ./vasp_results/

# 2. Add group labels
gpumdkit.sh -addgroup POSCAR Li Y Cl

# 3. Sample diverse structures
gpumdkit.sh  # Select: 2) Sample Structures -> 203) FPS by NepTrain

# 4. Train NEP model (external)
# ... train nep.txt ...

# 5. Validate training
gpumdkit.sh -plt train
gpumdkit.sh -plt prediction

# 6. Check data quality
gpumdkit.sh -range train.xyz force
gpumdkit.sh -min_dist_pbc train.xyz
```

### Workflow 2: Ionic Conductivity Calculation
```bash
# 1. Run MD simulation with GPUMD
# ... generate msd.out ...

# 2. Calculate MSD
gpumdkit.sh -calc msd trajectory.xyz Li 10

# 3. Calculate ionic conductivity
gpumdkit.sh -calc ionic-cond Li 1

# 4. Plot results
gpumdkit.sh -plt msd
gpumdkit.sh -plt sdc

# 5. Arrhenius analysis (multiple temperatures)
gpumdkit.sh -plt arrhenius_sigma
gpumdkit.sh -plt arrhenius_d
```

### Workflow 3: Structure Analysis Pipeline
```bash
# 1. Check composition
gpumdkit.sh -analyze_comp train.xyz

# 2. Check minimum distances
gpumdkit.sh -min_dist_pbc train.xyz

# 3. Filter structures
gpumdkit.sh  # Select: 5) Analyzer -> 506) Filter by distance

# 4. Check energy/force range
gpumdkit.sh -range train.xyz force
gpumdkit.sh -range train.xyz energy
```

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

## Troubleshooting

**Issue**: Command not found
**Solution**: Ensure GPUMDkit is installed and in PATH: `source ./install.sh`

**Issue**: Missing Python packages
**Solution**: Activate conda environment: `conda activate gpumdkit`

**Issue**: Script errors
**Solution**: Check input file format and required parameters with `-h` flag

## Detailed Documentation

- [Format Conversion Guide](../../docs/tutorials/format_conversion.md)
- [Calculators Guide](../../docs/tutorials/calculators.md)
- [Analyzer Guide](../../docs/tutorials/analyzer.md)
- [Visualization Guide](../../docs/tutorials/plot_scripts.md)
- [Workflow Guide](../../docs/tutorials/workflow.md)
- [Sampling Guide](../../docs/tutorials/sample_structures.md)
