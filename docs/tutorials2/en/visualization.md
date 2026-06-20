<div align="center">
  <h1>Visualization</h1>
  <p>
    <strong>English</strong> | <a href="../zh/visualization.md">简体中文</a>
  </p>
</div>

Plotting and visualization tools for simulation results.

## Quick Reference

```bash
gpumdkit.sh -plt <type>              # Display plot
gpumdkit.sh -plt <type> save         # Save as PNG
gpumdkit.sh -plt <type> -h           # Get help
```

## NEP Training Plots

| Command | Description |
|---------|-------------|
| `train` | Training loss curves and parity plots |
| `prediction` / `test` | Test set parity plots |
| `train_test` | Combined train/test parity plots |
| `force_errors` | Force error analysis |
| `des` | Descriptor PCA/UMAP visualization |
| `dimer` | Dimer interaction curves |
| `charge` | Charge distribution (qNEP) |
| `born_charge` / `bec` | Born effective charges |

```bash
gpumdkit.sh -plt train
gpumdkit.sh -plt prediction
gpumdkit.sh -plt force_errors
gpumdkit.sh -plt des pca
gpumdkit.sh -plt dimer Li Li nep.txt
```

## Transport Plots

| Command | Description |
|---------|-------------|
| `msd` | Mean square displacement |
| `sdc` | Self-diffusion coefficient |
| `msd_sdc` | MSD and SDC combined |
| `arrhenius_sigma` | Arrhenius ionic conductivity |
| `arrhenius_d` | Arrhenius diffusivity |

```bash
gpumdkit.sh -plt msd
gpumdkit.sh -plt sdc
gpumdkit.sh -plt arrhenius_sigma
gpumdkit.sh -plt arrhenius_d
```

## Structural Plots

| Command | Description |
|---------|-------------|
| `thermo` | Thermodynamic properties |
| `rdf` | Radial distribution function |
| `rdf_pmf` | RDF + potential of mean force |
| `vac` | Velocity autocorrelation |
| `doas` | Density of atomistic states |
| `net_force` | Net force distribution |

```bash
gpumdkit.sh -plt thermo
gpumdkit.sh -plt rdf
gpumdkit.sh -plt doas doas.out Li
gpumdkit.sh -plt net_force train.xyz
```

## Heat Transport Plots

| Command | Description |
|---------|-------------|
| `emd` | EMD thermal conductivity |
| `nemd` | NEMD thermal transport |
| `hnemd` | HNEMD thermal transport |
| `viscosity` | Viscosity components |

```bash
gpumdkit.sh -plt emd x
gpumdkit.sh -plt nemd <real_length> <scale_eff_size> <cutoff_freq> save
gpumdkit.sh -plt hnemd <scale_eff_size> <cutoff_freq> save
```

## Output Files

| Plot | PNG Filename |
|------|--------------|
| `train` | `train.png` |
| `prediction` | `prediction.png` |
| `msd` | `msd.png` |
| `sdc` | `sdc.png` |
| `thermo` | `thermo.png` |
| `rdf` | `rdf.png` |
| `arrhenius_sigma` | `Arrhenius_sigma.png` |
| `arrhenius_d` | `Arrhenius_D.png` |

## Dependencies

```bash
pip install matplotlib numpy
```
