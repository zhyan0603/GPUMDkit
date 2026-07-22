<div align="center">
  <h1>📚 GPUMDkit Tutorials</h1>
  <p style="text-align: justify;">Welcome to <strong>GPUMDkit</strong>, a command-line toolkit for GPUMD and NEP.</p>
</div>

## What it does

GPUMDkit helps you perform common tasks in computational materials science without writing custom scripts. It supports format conversion, structure analysis, property calculation, visualization, and batch workflows — all from one command-line tool.

## Common Tasks

| I want to... | Tutorial |
|--------------|----------|
| Install GPUMDkit and run my first command | [Quick Start](quick_start.md) |
| Plan a GPUMD/NEP simulation and its post-processing | [Simulation and Post-processing](simulation_and_postprocessing.md) |
| Convert VASP, LAMMPS, CP2K, or CIF files to extxyz | [Format Conversion](format_conversion.md) |
| Check structure distances, filter datasets, or find outliers | [Analyzer Scripts](analyzer_scripts.md) |
| Calculate MSD, ionic conductivity, or descriptors | [Calculator Scripts](calculator_scripts.md) |
| Analyze polar materials, ferroelectrics, or ABO3 systems | [Polar Material Analysis](polar_material_analysis.md) |
| Plot NEP training results or thermodynamic data | [Plot Scripts](plot_scripts.md) |
| Run batch DFT or MD simulations | [Workflow Scripts](workflow_scripts.md) |
| Select representative structures from a dataset | [Structure Sampling](structure_sampling.md) |
| Add my own shortcuts to GPUMDkit | [Custom Commands](custom_commands.md) |

## Before you start

### Installation

#### Conda (Recommended)

```bash
conda create -n gpumdkit -c gpumdkit -c conda-forge gpumdkit
conda activate gpumdkit
```

Some features require optional packages:

```bash
pip install neptrain calorine
```

#### From Source

```bash
git clone https://github.com/zhyan0603/GPUMDkit.git
cd GPUMDkit
source ./install.sh
```

## Interactive mode

```bash
gpumdkit.sh
```

This opens a menu where you can select modules by number.

## CLI mode

```bash
gpumdkit.sh -<option> [args...]
```

For example:

```bash
gpumdkit.sh -pos2exyz POSCAR model.xyz
gpumdkit.sh -plt train
gpumdkit.sh -calc msd trajectory.xyz Li 10
```

Use `gpumdkit.sh -h` to see all available options.

## Module Overview

| Module | Menu | CLI | Description |
|--------|------|-----|-------------|
| Format Conversion | 1 | `-out2xyz`, `-pos2exyz` | File conversion |
| Sample Structures | 2 | interactive menu, `-frame_range`, `-pynep` (deprecated) | Structure sampling |
| Workflow | 3 | - | Batch processing |
| Calculators | 4 | `-calc <type>` | Property calculations |
| Analyzer | 5 | `-range`, `-min_dist` | Structure analysis |
| Visualization | 6 | `-plt <type>` | Plotting tools |
| Utilities | 7 | `-time`, `-clean` | Utility functions |

## All Tutorials

| Tutorial | Description |
|----------|-------------|
| [Quick Start](quick_start.md) | Installation and first steps |
| [Simulation and Post-processing](simulation_and_postprocessing.md) | End-to-end GPUMD/NEP workflow and Arrhenius example |
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

## Links

- GitHub: https://github.com/zhyan0603/GPUMDkit
- Documentation: https://gpumdkit.cn/
