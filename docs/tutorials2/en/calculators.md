# Calculators Guide

<div align="center">
  <p>
    <a href="../zh/calculators.md">中文</a> | <strong>English</strong>
  </p>
</div>

This guide covers all calculator tools in GPUMDkit for computing material properties.

## Available Calculators

| Calculator | Command | Description |
|------------|---------|-------------|
| Ionic Conductivity | `-calc ionic-cond` | Calculate ionic conductivity from MSD |
| NEP Properties | `-calc nep` | Compute energy/force/stress with NEP |
| Descriptors | `-calc des` | Calculate NEP descriptors for analysis |
| DOAS | `-calc doas` | Density of atomistic states |
| NEB | Direct Python | Nudged elastic band calculation |
| Neighbor List | `-calc nlist` | Build neighbor lists for analysis |
| Displacement | `-calc disp` | Calculate atomic displacements |
| Averaged Structure | `-calc avg-struct` | Time-averaged structure from trajectory |
| Octahedral Tilt | `-calc oct-tilt` | Perovskite octahedral tilt angles |
| Polarization | `-calc pol-abo3` | ABO3 local polarization |
| Minimization | `-calc minimize` | Structure relaxation with NEP |
| MSD | `-calc msd` | Mean square displacement from trajectory |

## Interactive Mode

```bash
gpumdkit.sh
# Select: 4) Calculators
```

You'll see:

```
+----------------------------------------------------------+
|                     CALCULATOR TOOLS                     |
+----------------------------------------------------------+
| 401) Calc ionic conductivity                             |
| 402) Calc properties by nep                              |
| 403) Calc descriptors of specific elements               |
| 404) Calc density of atomistic states (DOAS)             |
| 405) Calc nudged elastic band (NEB) by nep               |
| 406) Build neighbor list                                 |
| 407) Calc displacement from trajectory                   |
| 408) Calc averaged structure                             |
| 409) Calc octahedral tilt                                |
| 410) Calc polarization for ABO3                          |
| 411) Minimize structure by nep                           |
| 412) Calc mean square displacement (MSD) from trajectory |
+----------------------------------------------------------+
```

## Command Reference

### Ionic Conductivity

Calculate ionic conductivity from MSD data using the Nernst-Einstein equation.

```bash
gpumdkit.sh -calc ionic-cond <element> <charge>

# Examples
gpumdkit.sh -calc ionic-cond Li 1    # Lithium ion (Li+)
gpumdkit.sh -calc ionic-cond Na 1    # Sodium ion (Na+)
gpumdkit.sh -calc ionic-cond O -2    # Oxygen ion (O2-)
```

**Required files in current directory:**
- `msd.out` - From GPUMD `compute_msd`
- `thermo.out` - For temperature
- `model.xyz` - For volume
- `run.in` - For simulation parameters

**Output:** Ionic diffusivity and conductivity values

### NEP Property Prediction

Calculate energies, forces, and stresses using a NEP model as a DFT surrogate.

```bash
gpumdkit.sh -calc nep <input.xyz> <output.xyz> <nep.txt>

# Example
gpumdkit.sh -calc nep structures.xyz predictions.xyz nep.txt
```

**Note:** Clean input with `gpumdkit.sh -clean_xyz` first to remove existing properties.

### Descriptors

Calculate NEP descriptors for dimensionality reduction and structure analysis.

```bash
gpumdkit.sh -calc des <input.xyz> <output.npy> <nep.txt> <element>

# Example
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li
```

**Visualize with:**
```bash
gpumdkit.sh -plt des pca    # PCA visualization
gpumdkit.sh -plt des umap   # UMAP visualization
```

### Density of Atomistic States (DOAS)

Calculate per-atom energy distributions grouped by element type.

```bash
gpumdkit.sh -calc doas <input.xyz> <nep.txt> <output.txt>

# Example
gpumdkit.sh -calc doas structures.xyz nep.txt doas.out
```

**Visualize with:**
```bash
gpumdkit.sh -plt doas doas.out Li
```

### Mean Square Displacement (MSD)

Calculate directional MSD from an extxyz trajectory.

```bash
gpumdkit.sh -calc msd <trajectory.xyz> <element> <dt_fs> [max_corr_steps]

# Example: Li with 10 fs timestep
gpumdkit.sh -calc msd dump.xyz Li 10
```

**Output:** `msd.out` (Time/ps, MSD_x, MSD_y, MSD_z)

### Neighbor List

Build neighbor lists for perovskite analysis.

```bash
gpumdkit.sh -calc nlist -i <input> -c <cutoff> -n <num_neighbors> -C <center_elements> -E <neighbor_elements>

# Examples
gpumdkit.sh -calc nlist -i model.xyz -c 4 -n 12 -C Ti -E O
gpumdkit.sh -calc nlist -i model.xyz -c 4 -n 12 -C Pb Sr -E O
```

**Output:** `nl-<center>-<neighbor>.dat`

### Displacement

Calculate atomic displacements from trajectory using neighbor list.

```bash
gpumdkit.sh -calc disp -i <trajectory.xyz> -n <neighbor_list> -o <output>

# Example
gpumdkit.sh -calc disp -i movie.xyz -n nl-Pb-O.dat -o displacements.dat
```

**Optional frame slicing:** `-s <start> -t <stop> -p <step> -l <last_fraction>`

### Averaged Structure

Generate time-averaged structure from trajectory.

```bash
gpumdkit.sh -calc avg-struct -i <trajectory.xyz> -l <fraction> -o <output>

# Example: Average last 20% of frames
gpumdkit.sh -calc avg-struct -i movie.xyz -l 0.2 -o averaged.xyz
```

### Octahedral Tilt

Calculate octahedral tilt angles for perovskite systems.

```bash
gpumdkit.sh -calc oct-tilt -i <input.xyz> -n <B-O neighbor list> -o <output>

# Example
gpumdkit.sh -calc oct-tilt -i model.xyz -n nl-Ti-O.dat -o octahedral_tilt.dat
```

### Polarization (ABO3)

Calculate local polarization for ABO3 perovskites.

```bash
gpumdkit.sh -calc pol-abo3 -i <input.xyz> \
  --nl-ba <B-A neighbor list> \
  --nl-bo <B-O neighbor list> \
  --bec <Element=charge ...>

# Example
gpumdkit.sh -calc pol-abo3 -i model.xyz \
  --nl-ba nl-Ti-Pb.dat \
  --nl-bo nl-Ti-O.dat \
  --bec Pb=2.5 Sr=2.0 Ti=4.0 O=-2.0
```

### Structure Minimization

Relax structures using BFGS optimizer with NEP model.

```bash
gpumdkit.sh -calc minimize <structure> <nep.txt> [fmax] [max_steps]

# Example
gpumdkit.sh -calc minimize POSCAR nep.txt 0.01 1000
```

**Parameters:**
- `fmax`: Force convergence threshold (default: 0.01 eV/Å)
- `max_steps`: Maximum optimization steps (default: 1000)

**Output:** `minimize.xyz` (optimization trajectory)

### NEB Calculation

Perform nudged elastic band calculations for migration barriers.

```bash
# Direct Python execution
python Scripts/calculators/neb_calculation.py <initial.xyz> <final.xyz> <n_images> <nep.txt>

# Example
python Scripts/calculators/neb_calculation.py init.xyz fin.xyz 9 nep.txt
```

## Common Workflows

### Ionic Transport Analysis

```bash
# 1. Run MD with compute_msd in run.in
# 2. Calculate MSD
gpumdkit.sh -calc msd dump.xyz Li 10

# 3. Calculate conductivity
gpumdkit.sh -calc ionic-cond Li 1

# 4. Visualize
gpumdkit.sh -plt msd
gpumdkit.sh -plt sdc
```

### Descriptor Analysis

```bash
# 1. Calculate descriptors
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li

# 2. Visualize with PCA
gpumdkit.sh -plt des pca

# 3. Or visualize with UMAP
gpumdkit.sh -plt des umap
```

### Perovskite Analysis

```bash
# 1. Build neighbor lists
gpumdkit.sh -calc nlist -i model.xyz -c 4 -n 12 -C Ti -E O
gpumdkit.sh -calc nlist -i model.xyz -c 4 -n 12 -C Ti -E Pb

# 2. Calculate displacements
gpumdkit.sh -calc disp -i movie.xyz -n nl-Ti-O.dat -o disp.dat

# 3. Calculate tilt
gpumdkit.sh -calc oct-tilt -i movie.xyz -n nl-Ti-O.dat -o tilt.dat

# 4. Calculate polarization
gpumdkit.sh -calc pol-abo3 -i movie.xyz \
  --nl-ba nl-Ti-Pb.dat --nl-bo nl-Ti-O.dat \
  --bec Pb=2.5 Ti=4.0 O=-2.0
```

### Structure Relaxation

```bash
# 1. Minimize structure
gpumdkit.sh -calc minimize POSCAR nep.txt 0.01 1000

# 2. Use minimized structure for MD
cp minimize.xyz relaxed_model.xyz
```

## Dependencies

| Calculator | Required Packages |
|------------|------------------|
| All | `numpy`, `ase` |
| ionic-cond | `scipy` |
| nep, des, doas, minimize | `calorine` |
| nlist, disp, oct-tilt, pol-abo3 | `ferrodispcalc` |
| neb | `calorine`, `matplotlib` |

Install dependencies:
```bash
pip install numpy ase scipy calorine matplotlib tqdm
pip3 install git+https://github.com/MoseyQAQ/ferrodispcalc.git
```

## See Also

- [Visualization](visualization.md) - Plot calculation results
- [Analyzers](analyzers.md) - Structure analysis
- [NEP Training Guide](nep_training.md) - Complete training workflow
