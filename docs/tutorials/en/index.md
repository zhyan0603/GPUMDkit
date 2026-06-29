<div align="center">
  <h1>📚 GPUMDkit Tutorials</h1>
  <p style="text-align: justify;">Welcome to <strong>GPUMDkit</strong>, a command-line toolkit for GPUMD and NEP.</p>
</div>

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
pip install neptrain dpdata calorine
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
| [Quick Start](quick_start.md) | Installation and first steps |
| [Command Reference](command_reference.md) | One-page CLI and menu reference |
| [Format Conversion](format_conversion.md) | Convert between file formats |
| [Calculator Scripts](calculator_scripts.md) | Compute material properties |
| [Analyzer Scripts](analyzer_scripts.md) | Structure analysis and validation |
| [Plot Scripts](plot_scripts.md) | Plotting and visualization |
| [Workflow Scripts](workflow_scripts.md) | Batch processing and automation |
| [Structure Sampling](structure_sampling.md) | Structure selection methods |
| [Custom Commands](custom_commands.md) | User-defined GPUMDkit shortcuts |
| [Active Learning Workflow](active_learning_workflow.md) | Batch active-learning workflow notes |
| [Polar Material Analysis](polar_material_analysis.md) | Ferroelectric and polarization tools |
| [Contributing to GPUMDkit](contributing_to_gpumdkit.md) | Development notes |

## Example: Convert VASP OUTCARs

```bash
gpumdkit.sh -out2xyz ./vasp_results/
```

## Example: Plot NEP Results

```bash
gpumdkit.sh -plt train
gpumdkit.sh -plt prediction
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
| Sample Structures | 2 | `-neptrain`, interactive menu | Structure sampling |
| Workflow | 3 | - | Batch processing |
| Calculators | 4 | `-calc <type>` | Property calculations |
| Analyzer | 5 | `-range`, `-min_dist` | Structure analysis |
| Visualization | 6 | `-plt <type>` | Plotting tools |
| Utilities | 7 | `-time`, `-clean` | Utility functions |

## Links

- GitHub: https://github.com/zhyan0603/GPUMDkit
- Documentation: https://zhyan0603.github.io/GPUMDkit/
