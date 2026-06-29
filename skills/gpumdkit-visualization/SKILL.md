---
name: gpumdkit-visualization
description: >
  Use when plotting molecular dynamics simulation results, NEP training data, or transport properties.
  Supports 31+ plot types including training loss, parity plots, MSD, SDC, RDF, thermodynamic properties,
  Arrhenius plots, phonon DOS, thermal conductivity, and descriptor visualization.
  Use when user asks about: plotting, visualization, training results, MSD plot, RDF plot,
  thermal conductivity, Arrhenius plot, or NEP training visualization.
allowed-tools: Bash(gpumdkit *)
---

# GPUMDkit Visualization

## Quick Reference

```bash
gpumdkit.sh -plt <type>              # Display plot interactively
gpumdkit.sh -plt <type> save         # Save plot as PNG
gpumdkit.sh -plt <type> -h           # Get help for specific plot
```

## Plot Categories

### NEP Training & Evaluation (12 plot types)

| Command | Input Files | Description |
|---------|-------------|-------------|
| `train` | `loss.out`, `*_train.out` | Training loss curves and parity plots |
| `prediction` / `test` | `*_test.out` | Test set parity plots |
| `train_test` | `*_train.out`, `*_test.out` | Combined train/test parity plots |
| `parity_density` | `*_train.out` | Density-based parity plots for large datasets |
| `train_density` | `loss.out`, `*_train.out` | Training with density visualization |
| `force_errors` | `force_train.out` | Force error analysis |
| `restart` | `nep.restart` | NEP restart file visualization |
| `charge` | `charge_train.out` | Charge distribution (qNEP) |
| `born_charge` / `bec` | `bec_train.out`, `bec_test.out` | Born effective charges |
| `dimer` | NEP model | Dimer interaction curves |
| `des` | `descriptors.npy` | Descriptor PCA/UMAP visualization |
| `lr` | `loss.out` (gnep) | Learning rate decay |

```bash
# Training results
gpumdkit.sh -plt train
gpumdkit.sh -plt prediction
gpumdkit.sh -plt force_errors

# Descriptor visualization (requires prior calculation)
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li
gpumdkit.sh -plt des pca
gpumdkit.sh -plt des umap

# Dimer plot
gpumdkit.sh -plt dimer Li Li nep.txt
```

### Transport Properties (9 plot types)

| Command | Input Files | Description |
|---------|-------------|-------------|
| `msd` | `msd.out` | Mean square displacement |
| `msd_conv` | `msd_step*.out` | MSD convergence check |
| `msd_all` | `msd.out` (all_groups) | MSD per species |
| `sdc` | `msd.out` | Self-diffusion coefficient |
| `msd_sdc` | `msd.out` | MSD and SDC combined |
| `sigma` / `arrhenius_sigma` | `*K/` directories | Arrhenius ionic conductivity |
| `D` / `arrhenius_d` | `*K/` directories | Arrhenius diffusivity |
| `sigma_xyz` | `*K/` directories | Directional Arrhenius conductivity |
| `D_xyz` | `*K/` directories | Directional Arrhenius diffusivity |

```bash
# MSD and diffusion
gpumdkit.sh -plt msd
gpumdkit.sh -plt sdc
gpumdkit.sh -plt msd_sdc

# MSD for all species (requires all_groups in run.in)
gpumdkit.sh -plt msd_all Li P S

# Arrhenius plots (requires temperature-organized directories)
# Directory structure: 300K/, 350K/, 400K/, ... each containing msd.out
gpumdkit.sh -plt arrhenius_sigma
gpumdkit.sh -plt arrhenius_d
```

### Structural Analysis (10 plot types)

| Command | Input Files | Description |
|---------|-------------|-------------|
| `thermo` | `thermo.out` | Thermodynamic properties |
| `thermo2` | `thermo.out` | Alternative thermo style |
| `thermo3` | `thermo.out` | Third thermo style |
| `rdf` | `rdf.out` | Radial distribution function |
| `rdf_pmf` | `rdf.out` | RDF + potential of mean force |
| `vac` | `sdc.out` | Velocity autocorrelation |
| `cohesive` | `cohesive.out` | Cohesive energy curve |
| `net_force` | extxyz file | Net force distribution |
| `doas` | `doas.out` | Density of atomistic states |
| `plane-grid` | `model.xyz`, `displacements.dat` | Displacement grid visualization |

```bash
# Thermodynamic properties
gpumdkit.sh -plt thermo

# RDF analysis
gpumdkit.sh -plt rdf
gpumdkit.sh -plt rdf 2              # Specific column
gpumdkit.sh -plt rdf_pmf 300        # With PMF at 300K

# DOAS visualization (requires prior calculation)
gpumdkit.sh -plt doas doas.out Li

# Plane-grid displacement
gpumdkit.sh -plt plane-grid -i model.xyz -d displacements.dat -e Pb Sr
```

### Heat Transport (4 plot types)

| Command | Input Files | Description |
|---------|-------------|-------------|
| `emd` | EMD outputs | EMD thermal conductivity |
| `nemd` | NEMD outputs | NEMD thermal transport |
| `hnemd` | HNEMD outputs | HNEMD thermal transport |
| `viscosity` | `viscosity.out` | Viscosity components |

```bash
# EMD thermal conductivity
gpumdkit.sh -plt emd x

# NEMD thermal transport
# Parameters: real_length scale_eff_size cutoff_freq
gpumdkit.sh -plt nemd <real_length> <scale_eff_size> <cutoff_freq> save

# HNEMD thermal transport
gpumdkit.sh -plt hnemd <scale_eff_size> <cutoff_freq> save

# Viscosity
gpumdkit.sh -plt viscosity save
```

### Phonons (1 plot type)

| Command | Input Files | Description |
|---------|-------------|-------------|
| `pdos` | `model.xyz`, `run.in`, `dos.out`, `mvac.out` | Phonon DOS and heat capacity |

```bash
gpumdkit.sh -plt pdos save
```

## Common Workflows

### NEP Training Validation
```bash
# 1. Plot training loss
gpumdkit.sh -plt train
# 2. Check test predictions
gpumdkit.sh -plt prediction
# 3. Analyze force errors
gpumdkit.sh -plt force_errors
# 4. Visualize descriptors
gpumdkit.sh -plt des pca
```

### Diffusion Analysis
```bash
# 1. Plot MSD
gpumdkit.sh -plt msd
# 2. Plot self-diffusion coefficient
gpumdkit.sh -plt sdc
# 3. Combined MSD-SDC plot
gpumdkit.sh -plt msd_sdc
# 4. Arrhenius analysis (multi-temperature)
gpumdkit.sh -plt arrhenius_d
```

### Thermal Transport
```bash
# 1. Plot thermodynamic properties
gpumdkit.sh -plt thermo
# 2. Plot thermal conductivity
gpumdkit.sh -plt emd x
# 3. Or NEMD/HNEMD
gpumdkit.sh -plt nemd 10 1 60 save
```

## Output Files

| Plot Type | PNG Filename (with `save`) |
|-----------|---------------------------|
| `train` | `train.png` |
| `prediction` | `prediction.png` |
| `train_test` | `train_test.png` |
| `force_errors` | `force_errors.png` |
| `des` | `descriptors.png` |
| `msd` | `msd.png` |
| `sdc` | `sdc.png` |
| `msd_sdc` | `MSD_SDC.png` |
| `thermo` | `thermo.png` |
| `rdf` | `rdf.png` |
| `arrhenius_sigma` | `Arrhenius_sigma.png` |
| `arrhenius_d` | `Arrhenius_D.png` |
| `emd` | `emd.png` |
| `nemd` | `nemd.png` |
| `hnemd` | `hnemd.png` |
| `viscosity` | `viscosity.png` |
| `cohesive` | `Cohesive.png` |

## Dependencies

All plotting scripts require:
- `matplotlib`
- `numpy`

## Detailed Documentation

See [plot_scripts.md](../../docs/tutorials/en/plot_scripts.md) for comprehensive guide.
