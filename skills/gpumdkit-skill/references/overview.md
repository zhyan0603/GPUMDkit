# GPUMDkit Overview and Command Discovery

## Contents

- Quick start
- Module overview
- Common workflows
- Output files
- Utility commands
- Troubleshooting and documentation

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
gpumdkit.sh -<option> -h          # Python-backed command help when implemented
gpumdkit.sh -plt -h               # Plot help
gpumdkit.sh -calc -h              # Calculator help
```

For current Python-backed commands, detailed usage, type conversion, file checks, and errors live in the target Python script and are reached through `-h`. `gpumdkit.sh` mainly routes arguments; large dispatchers provide broad help. Legacy Shell/menu routes such as `-pynep` and `-time` do not implement the same per-command help contract, so inspect their module reference or interactive prompt instead of assuming `-h` support.

## Module Overview

| Module | Description | Reference |
|--------|-------------|-----------|
| Format Conversion | Convert between VASP, LAMMPS, CP2K, ABACUS, CIF, extxyz | `references/format-conversion.md` |
| Calculators | Compute ionic conductivity, descriptors, MSD, NEB, etc. | `references/calculators.md` |
| Analyzers | Structure validation, filtering, composition analysis | `references/analyzers.md` |
| Visualization | Plot training results, transport properties, structural data | `references/visualization.md` |
| Workflows | Batch processing for DFT and MD simulations | `references/workflows.md` |
| Sampling | Structure selection using uniform, random, or FPS methods | `references/sampling.md` |
| GPUMD/NEP | Plan simulation, training, validation, and post-processing | `references/gpumd.md`, `references/nep.md` |

## Utility Commands

| Command | Behavior and safety boundary |
|---|---|
| `gpumdkit.sh -skill` | Print the canonical unified Skill path and cross-client installation hints |
| `gpumdkit.sh -time <gpumd|nep>` | Legacy time-consumption analyzer; use only the supported `gpumd` or `nep` selector |
| `gpumdkit.sh -nep_modifier` | Launch the interactive NEP model modifier; inspect the source model and obtain overwrite authorization first |
| `gpumdkit.sh -clean` | Remove generated/extra files from the current directory; preview the cleanup implementation and get explicit deletion approval before running |
| `gpumdkit.sh -update` | Run GPUMDkit's networked self-update; inspect worktree changes and get explicit update authorization first |

Commands shown in the custom-command tutorials, such as `-greet`, `-batch_plot`, or `-prep_training`, are examples for extending GPUMDkit and are not built-in commands unless the user has installed those customizations.

## Common Workflows

### Workflow 1: NEP Model Training Pipeline
```bash
# 1. Convert DFT data to extxyz
gpumdkit.sh -out2xyz ./vasp_results/

# 2. Add group labels only if the downstream workflow uses GPUMD groups
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
# Path A: calculate MSD from a stored extxyz trajectory
gpumdkit.sh -calc msd trajectory.xyz Li 10

# Path B: use GPUMD compute_msd to generate msd.out directly
# Do not run both paths unless comparing implementations.

# Calculate ionic conductivity from the validated msd.out/model.xyz/run.in inputs
gpumdkit.sh -calc ionic-cond Li 1

# Plot results
gpumdkit.sh -plt msd
gpumdkit.sh -plt sdc

# Arrhenius analysis (multiple temperatures)
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

Read the matching page under `${GPUMDkit_path}/docs/tutorials/en/` or `/zh/` when user-facing tutorial detail is needed.
