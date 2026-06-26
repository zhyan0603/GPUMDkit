<div align="center">
  <h1>GPUMDkit Tutorials</h1>
  <p>
    <strong>English</strong> | <a href="../zh/README.md">简体中文</a>
  </p>
</div>

Welcome to **GPUMDkit** - a command-line toolkit for GPUMD and NEP.

## What is GPUMDkit?

GPUMDkit streamlines common tasks in computational materials science:

- **Format Conversion**: Convert between VASP, LAMMPS, CP2K, ABACUS, CIF, and extxyz formats
- **Visualization**: Plotting tools for NEP training, MD simulations, and analysis
- **Calculators**: Compute ionic conductivity, descriptors, MSD, NEB, and more
- **Analyzers**: Structure validation, filtering, composition analysis
- **Workflows**: Batch processing for DFT and MD simulations
- **Sampling**: Structure selection using uniform, random, or FPS methods

## Quick Start

### Installation

```bash
git clone https://github.com/zhyan0603/GPUMDkit.git
cd GPUMDkit
source ./install.sh
```

### Dependencies

```bash
conda create -n gpumdkit python=3.12
conda activate gpumdkit
pip install neptrain ase pymatgen dpdata numpy scipy matplotlib
```

### Usage

```bash
# Interactive mode
gpumdkit.sh

# Command-line mode
gpumdkit.sh -<option> [args...]

# Get help
gpumdkit.sh -h
```

## Tutorials

| Tutorial | Description |
|----------|-------------|
| [Quick Start](quickstart.md) | Installation and first steps |
| [Format Conversion](format_conversion.md) | Convert between file formats |
| [Calculators](calculators.md) | Compute material properties |
| [Analyzers](analyzers.md) | Structure analysis and validation |
| [Visualization](visualization.md) | Plotting and visualization |
| [Workflows](workflows.md) | Batch processing and automation |
| [Structure Sampling](sampling.md) | Structure selection methods |
| [NEP Training](nep_training.md) | Complete NEP training workflow |

## Example: NEP Training Pipeline

```bash
# Convert DFT data
gpumdkit.sh -out2xyz ./vasp_results/

# Add group labels
gpumdkit.sh -addgroup POSCAR Li Y Cl

# Sample structures
gpumdkit.sh  # Select: 2) Sample Structures -> 203

# Validate training
gpumdkit.sh -plt train
```

## Example: Ionic Conductivity

```bash
gpumdkit.sh -calc msd trajectory.xyz Li 10
gpumdkit.sh -calc ionic-cond Li 1
gpumdkit.sh -plt msd
```

## Module Overview

| Module | Menu | CLI | Description |
|--------|------|-----|-------------|
| Format Conversion | 1 | `-out2xyz`, `-pos2exyz` | File conversion |
| Sample Structures | 2 | `-pynep` | Structure sampling |
| Workflow | 3 | - | Batch processing |
| Calculators | 4 | `-calc <type>` | Property calculations |
| Analyzer | 5 | `-range`, `-min_dist` | Structure analysis |
| Visualization | 6 | `-plt <type>` | Plotting tools |
| Utilities | 7 | `-time`, `-clean` | Utility functions |

## Links

- GitHub: https://github.com/zhyan0603/GPUMDkit
- Documentation: https://zhyan0603.github.io/GPUMDkit/
