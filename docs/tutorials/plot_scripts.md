# Plot Scripts

GPUMDkit provides comprehensive visualization tools for GPUMD and NEP data. These plotting scripts are primarily accessed via command-line mode, though some can also be accessed through the interactive menu.

## Command-Line Access

All plotting scripts are accessed via the `-plt` flag:

```bash
# General syntax
gpumdkit.sh -plt <plot_type> [options]

# View available plot types
gpumdkit.sh -plt

# Examples
gpumdkit.sh -plt thermo          # Plot thermodynamic properties
gpumdkit.sh -plt train           # Plot NEP training progress
gpumdkit.sh -plt msd             # Plot mean square displacement
```

## Interactive Mode Access

Some plotting functions are also available through interactive mode:

```bash
gpumdkit.sh
# Select: 6) Plot Scripts (if available in your version)
```

---

## Available Plots

### NEP Training & Prediction


### NEP Training & Prediction

**Training Results** (`plt_nep_train_results.py`)
```bash
gpumdkit.sh -plt train
```
Visualizes NEP training progress: loss curves, RMSE, parity plots.

![](../Gallery/train.png)

**Prediction Results** (`plt_nep_prediction_results.py`)
```bash
gpumdkit.sh -plt prediction
# or
gpumdkit.sh -plt test
```
Visualizes NEP predictions vs reference values.

![](../Gallery/prediction.png)

**Train + Test Comparison** (`plt_nep_train_test.py`)
```bash
gpumdkit.sh -plt train_test
```
Side-by-side parity plots for training and testing.

![](../Gallery/train_test.png)

**Parity Density Plots** (`plt_parity_density.py`)
```bash
gpumdkit.sh -plt parity_density
```
2D histogram density plots for large datasets.

![](../Gallery/parity_plot_density.png)

**Force Error Metrics** (`plt_force_errors.py`)
```bash
gpumdkit.sh -plt force_errors
```
Force error evaluation metrics from [Huang et al. (2023)](https://doi.org/10.1038/s41524-023-01123-3).

![](../Gallery/force_errors.png)

**Learning Rate** (`plt_learning_rate.py`)
```bash
gpumdkit.sh -plt lr
```
Visualizes learning rate schedule during training.

**Restart Parameters** (`plt_nep_restart.py`)
```bash
gpumdkit.sh -plt restart
```
Visualizes parameters in `nep.restart` file.

![](../Gallery/nep_restart.png)

### Thermodynamic Properties

**Thermo Evolution** (`plt_nep_thermo.py`)
```bash
gpumdkit.sh -plt thermo
```
Plots temperature, pressure, energy, lattice parameters, volume from `thermo.out`.

![](../Gallery/thermo.png)

Alternative layouts:
```bash
gpumdkit.sh -plt thermo2
gpumdkit.sh -plt thermo3
```

### Transport Properties

**Mean Square Displacement** (`plt_msd.py`)
```bash
gpumdkit.sh -plt msd
```
Plots MSD vs time from `msd.out`.

![](../Gallery/msd.png)

**MSD for All Species** (`plt_msd_all.py`)
```bash
gpumdkit.sh -plt msd_all Li Y Cl
```
Plots MSD for each species separately. Requires `all_groups` in `compute_msd`.

![](../Gallery/msd_all.png)

**MSD Convergence Check** (`plt_msd_convergence_check.py`)
```bash
gpumdkit.sh -plt msd_conv
```
Checks MSD convergence. Requires `save_every` in `compute_msd`.

![](../Gallery/msd_convergence.png)

**Self-Diffusion Coefficient** (`plt_sdc.py`)
```bash
gpumdkit.sh -plt sdc
```
Plots SDC calculated from MSD.

![](../Gallery/sdc.png)

**Arrhenius Plots**
```bash
# Diffusivity vs 1/T
gpumdkit.sh -plt arrhenius_d

# Conductivity vs 1/T
gpumdkit.sh -plt arrhenius_sigma
```

**Velocity Autocorrelation** (`plt_vac.py`)
```bash
gpumdkit.sh -plt vac
```
Plots VAC from `sdc.out` or `vac.out`.

**Thermal Conductivity**
```bash
# EMD method
gpumdkit.sh -plt emd

# NEMD method
gpumdkit.sh -plt nemd

# HNEMD method
gpumdkit.sh -plt hnemd
```

![](../Gallery/emd.png)

### Structural Analysis

**Radial Distribution Function** (`plt_rdf.py`)
```bash
# Plot all RDF pairs
gpumdkit.sh -plt rdf

# Plot specific column
gpumdkit.sh -plt rdf 2
```

![](../Gallery/rdf1.png)
![](../Gallery/rdf2.png)

**Charge Distribution** (`plt_charge.py`)
```bash
gpumdkit.sh -plt charge
```
Plots charge distribution from NEP charge predictions.

![](../Gallery/charge.png)

### Advanced Analysis

**Descriptors** (`plt_descriptors.py`)
```bash
# First calculate descriptors
gpumdkit.sh -calc des umap train.xyz descriptors.npy nep.txt Li

# Then plot
gpumdkit.sh -plt des umap descriptors.npy
```
Visualizes NEP descriptors using UMAP/t-SNE/PCA. See [arXiv:2504.15925](https://doi.org/10.48550/arXiv.2504.15925).

![](../Gallery/des_umap.png)

**Dimer Interaction** (`plt_dimer.py`)
```bash
gpumdkit.sh -plt dimer Li Li nep.txt
```
Plots dimer interaction curves.

![](../Gallery/dimer_nep.png)

**Density of Atomistic States** (`plt_doas.py`)
```bash
gpumdkit.sh -plt doas file1.dat file2.dat
```
Plots DOAS from [Wang et al.](https://doi.org/10.1002/anie.202215544)

![](../Gallery/doas.png)

**Net Force Distribution** (`plt_net_force.py`)
```bash
gpumdkit.sh -plt net_force [options]
```
Plots net force distribution for quality checking. See [arXiv:2510.19774](https://arxiv.org/abs/2510.19774).

![](../Gallery/net_force_distribution.png)

## Common Workflows

### Monitor NEP Training
```bash
# During training
gpumdkit.sh -plt train

# After completion
gpumdkit.sh -plt train
gpumdkit.sh -plt force_errors
gpumdkit.sh -plt lr
```

### Analyze MD Simulation
```bash
gpumdkit.sh -plt thermo
gpumdkit.sh -plt msd
gpumdkit.sh -plt rdf
```

### Validate NEP Model
```bash
gpumdkit.sh -plt prediction
gpumdkit.sh -plt train_test
gpumdkit.sh -plt dimer Li Li nep.txt
```

## Saving Plots

Add `save` argument to save plots as PNG:
```bash
gpumdkit.sh -plt thermo save
gpumdkit.sh -plt train save
```

Plots are saved with descriptive names (e.g., `thermo.png`, `train.png`).

## Quick Reference Table

| Command | Input File(s) | Description |
|---------|---------------|-------------|
| `thermo` | `thermo.out` | Thermodynamic properties |
| `train` | `loss.out`, `*_train.out` | Training progress |
| `prediction` | `*_test.out` | Prediction accuracy |
| `train_test` | `*_train.out`, `*_test.out` | Training vs testing |
| `msd` | `msd.out` | Mean square displacement |
| `sdc` | `msd.out` | Self-diffusion coefficient |
| `rdf` | `rdf.out` | Radial distribution function |
| `charge` | `charge_train.out` | Charge distribution |
| `force_errors` | `force_train.out` | Force error metrics |
| `lr` | `loss.out` | Learning rate |
| `des` | `descriptors.npy` | Descriptor visualization |
| `emd/nemd/hnemd` | Various | Thermal conductivity |

## Tips

- Use interactive matplotlib window to zoom and explore plots
- Save plots for reports and presentations
- Check convergence with `msd_conv` before extracting diffusion data
- Use density plots for large datasets (>10k structures)

---

For more details, see [Scripts/plt_scripts/README.md](../../Scripts/plt_scripts/README.md)
