---
name: gpumdkit-calculators
description: >
  Use when calculating material properties from molecular dynamics data.
  Provides ionic conductivity, NEP property prediction, descriptors, DOAS, NEB, MSD,
  neighbor lists, displacements, octahedral tilt, polarization, and structure minimization.
  Use when user asks about: ionic conductivity, diffusion coefficient, descriptors, NEB,
  mean square displacement, density of atomistic states, or structure relaxation.
allowed-tools: Bash(gpumdkit *) Bash(python3 *)
---

# GPUMDkit Calculators

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

## Command Reference

### Ionic Conductivity
```bash
# Calculate ionic conductivity from MSD data
gpumdkit.sh -calc ionic-cond <element> <charge>

# Example: Lithium ion (Li+)
gpumdkit.sh -calc ionic-cond Li 1

# Example: Oxygen ion (O2-)
gpumdkit.sh -calc ionic-cond O -2

# Required files in current directory:
# - msd.out (from GPUMD compute_msd)
# - thermo.out (for temperature)
# - model.xyz (for volume)
# - run.in (for simulation parameters)
```

### NEP Property Prediction
```bash
# Calculate properties using NEP model
gpumdkit.sh -calc nep <input.xyz> <output.xyz> <nep.txt>

# Example
gpumdkit.sh -calc nep structures.xyz predictions.xyz nep.txt

# Note: Clean input with gpumdkit.sh -clean_xyz first to remove existing properties
```

### Descriptors
```bash
# Calculate NEP descriptors for specific element
gpumdkit.sh -calc des <input.xyz> <output.npy> <nep.txt> <element>

# Example
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li

# Visualize with: gpumdkit.sh -plt des pca
# Or: gpumdkit.sh -plt des umap
```

### Density of Atomistic States (DOAS)
```bash
# Calculate DOAS
gpumdkit.sh -calc doas <input.xyz> <nep.txt> <output.txt>

# Example
gpumdkit.sh -calc doas structures.xyz nep.txt doas.out

# Visualize with: gpumdkit.sh -plt doas doas.out Li
```

### Mean Square Displacement
```bash
# Calculate MSD from trajectory
gpumdkit.sh -calc msd <trajectory.xyz> <element> <dt_fs> [max_corr_steps]

# Example: Li with 10 fs timestep
gpumdkit.sh -calc msd dump.xyz Li 10

# Output: msd.out (Time/ps, MSD_x, MSD_y, MSD_z)
```

### Neighbor List
```bash
# Build neighbor list for perovskite analysis
gpumdkit.sh -calc nlist -i <input> -c <cutoff> -n <num_neighbors> -C <center_elements> -E <neighbor_elements>

# Example: Ti-O neighbors in BaTiO3
gpumdkit.sh -calc nlist -i model.xyz -c 4 -n 12 -C Ti -E O

# Example: Pb/Sr-O neighbors in PZT
gpumdkit.sh -calc nlist -i model.xyz -c 4 -n 12 -C Pb Sr -E O

# Output: nl-<center>-<neighbor>.dat
```

### Displacement
```bash
# Calculate displacements from trajectory
gpumdkit.sh -calc disp -i <trajectory.xyz> -n <neighbor_list> -o <output>

# Example
gpumdkit.sh -calc disp -i movie.xyz -n nl-Pb-O.dat -o displacements.dat

# Optional frame slicing: -s <start> -t <stop> -p <step> -l <last_fraction>
```

### Averaged Structure
```bash
# Calculate time-averaged structure
gpumdkit.sh -calc avg-struct -i <trajectory.xyz> -l <fraction> -o <output>

# Example: Average last 20% of frames
gpumdkit.sh -calc avg-struct -i movie.xyz -l 0.2 -o averaged.xyz
```

### Octahedral Tilt
```bash
# Calculate octahedral tilt angles
gpumdkit.sh -calc oct-tilt -i <input.xyz> -n <B-O neighbor list> -o <output>

# Example
gpumdkit.sh -calc oct-tilt -i model.xyz -n nl-Ti-O.dat -o octahedral_tilt.dat
```

### Polarization (ABO3)
```bash
# Calculate local polarization for perovskites
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
```bash
# Minimize structure with NEP
gpumdkit.sh -calc minimize <structure> <nep.txt> [fmax] [max_steps]

# Example
gpumdkit.sh -calc minimize POSCAR nep.txt 0.01 1000

# Output: minimize.xyz (optimization trajectory)
# Default fmax: 0.01 eV/A, Default max_steps: 1000
```

### NEB Calculation
```bash
# Direct Python execution (no CLI shortcut)
python Scripts/calculators/neb_calculation.py <initial.xyz> <final.xyz> <n_images> <nep.txt>

# Example
python Scripts/calculators/neb_calculation.py init.xyz fin.xyz 9 nep.txt

# Alternative with NepTrainKit:
python Scripts/calculators/neb_calculation_neptrain.py init.xyz fin.xyz 9 nep.txt
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
gpumdkit.sh -calc pol-abo3 -i movie.xyz --nl-ba nl-Ti-Pb.dat --nl-bo nl-Ti-O.dat --bec Pb=2.5 Ti=4.0 O=-2.0
```

### Structure Relaxation
```bash
# 1. Minimize structure
gpumdkit.sh -calc minimize POSCAR nep.txt 0.01 1000
# 2. Check result
gpumdkit.sh -plt thermo  # If running MD after relaxation
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

## Detailed Documentation

See [calculator_scripts.md](../../docs/tutorials/en/calculator_scripts.md) for comprehensive guide.
