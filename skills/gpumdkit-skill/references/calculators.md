# Calculators

## Contents

- Available calculators
- Command reference
- Common workflows
- Dependencies

## Available Calculators

| Calculator | Command | Description |
|------------|---------|-------------|
| Ionic Conductivity | `-calc ionic-cond` | Calculate ionic conductivity from MSD |
| NEP Properties | `-calc nep` | Compute energy/force/stress with NEP |
| Descriptors | `-calc des` | Calculate NEP descriptors for analysis |
| DOAS | `-calc doas` | Density of atomistic states |
| NEB | `-calc neb` | Nudged elastic band calculation |
| Neighbor List | `-calc nlist` | Build neighbor lists for analysis |
| Displacement | `-calc disp` | Calculate atomic displacements |
| Averaged Structure | `-calc avg-struct` | Time-averaged structure from trajectory |
| Octahedral Tilt | `-calc oct-tilt` | Perovskite octahedral tilt angles |
| Polarization | `-calc pol-abo3` | ABO3 local polarization |
| Minimization | `-calc minimize` | Structure relaxation with NEP |
| MSD | `-calc msd` | Mean square displacement from trajectory |
| XRD | `4 -> 413` (interactive only) | X-ray diffraction from an extxyz trajectory |
| Phonon band structure | `4 -> 414` (interactive only) | Calculate a NEP phonon band structure along `QPOINTS` |

## Command Reference

### Ionic Conductivity
```bash
# Calculate ionic conductivity from MSD data
gpumdkit.sh -calc ionic-cond <element> <charge>

# Example: Lithium ion (Li+)
gpumdkit.sh -calc ionic-cond Li 1

# Example: Oxygen ion (O2-)
gpumdkit.sh -calc ionic-cond O -2

# Preferred non-interactive inputs in current directory:
# - msd.out (from GPUMD compute_msd, or the first four columns from calc_msd.py)
# - thermo.out (for temperature)
# - model.xyz (for volume)
# - run.in (optional; used to detect replicate factors)
```

### NEP Property Prediction
```bash
# Calculate properties using NEP model
gpumdkit.sh -calc nep <input.xyz> <output.xyz> <nep.txt>

# Example
gpumdkit.sh -calc nep structures.xyz predictions.xyz nep.txt

# Note: Clean input with gpumdkit.sh -clean_xyz first to remove existing properties
```

### NEP Prediction Output Files
```bash
# Write energy, force, stress, and virial prediction files for every frame.
# The output suffix comes from the input filename (test.xyz -> *_test.out).
gpumdkit.sh -prediction <input.xyz> <nep.txt> [workers]

# Single-core default
gpumdkit.sh -prediction train.xyz nep.txt

# Eight CPU workers
gpumdkit.sh -prediction train.xyz nep.txt 8
```

The output files are `energy_<input-stem>.out`, `force_<input-stem>.out`,
`stress_<input-stem>.out`, and `virial_<input-stem>.out` in the current
directory. Predicted columns precede target columns. Energy and virial are
eV/atom, force is eV/Angstrom, and stress is GPa. Tensor columns use NEP's
`xx yy zz xy yz xz` order and the same sign convention for stress and virial.
If only one of stress or virial is supplied, the other target is derived using
`stress = virial / volume` and the appropriate unit conversion. Only frames
missing both use the `-1e6` sentinel; energy and force targets are required. A
`tqdm` progress bar is shown during prediction.

### Descriptors
```bash
# Calculate NEP descriptors for specific element
gpumdkit.sh -calc des <input.xyz> <output.npy> <nep.txt> <element>

# Example
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li

# Visualize with: gpumdkit.sh -plt des pca descriptors.npy
# Or: gpumdkit.sh -plt des umap descriptors.npy
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

# Output: msd.out with four columns: Time/ps, MSD_x, MSD_y, MSD_z
# This output is compatible with: gpumdkit.sh -plt msd
```

GPUMD's native `compute_msd` writes a different `msd.out`. For one selected
group it contains seven columns: time, `MSD_x/y/z`, and `SDC_x/y/z`; this
single-group layout is required by `gpumdkit.sh -plt sdc` and
`gpumdkit.sh -plt msd_sdc`. With `all_groups` or multiple groups, GPUMD
appends additional group data, so do not assume that every `msd.out` has seven
columns. The separate `compute_sdc` command writes `sdc.out`, which is used by
`gpumdkit.sh -plt vac`.

### X-ray Diffraction (interactive only)

XRD is intentionally exposed only through the interactive calculator menu:

```text
gpumdkit.sh -> 4) Calculators -> 413) Calc XRD from extxyz trajectory
```

The Python page asks for the input extxyz trajectory, output file, wavelength
in Angstrom, a `2theta` minimum and maximum, the number of bins over that exact
interval, an element selection (`all` or comma/space-separated symbols), and a
CPU worker count (`0` means automatic). It uses the input `Lattice`/`pbc`, all
atoms by default, standard LAMMPS scattering factors, and the
Lorentz-polarization factor. Orthogonal cells are supported; triclinic cells
are not.

The output contains averaged `Count` and `Count/Total` columns at the centers
of the selected `2theta` bins. The implementation streams the trajectory once,
caches the reciprocal mesh for unchanged cells, groups atoms by element, and
uses ordered threaded workers without changing the XRD formula or bin range.

### Phonon Band Structure (interactive only)

The phonon calculator is exposed through the interactive calculator menu:

```text
gpumdkit.sh -> 4) Calculators -> 414) Calc phonon band structure
```

It prompts for a primitive-cell structure, a Calorine-compatible NEP model, a
line-mode `QPOINTS` file, diagonal supercell factors, a displacement, and an
output file. The q-point path, segment count, and labels are read from the
supplied `QPOINTS` file; do not replace them with a hard-coded path.

The output is line-oriented phonon data in THz, normally written as
`phonon_NEP.dat`. The calculation requires `ASE`, `Calorine`, `phonopy`, and
`NumPy`.

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
# CLI shortcut
gpumdkit.sh -calc neb <initial.xyz> <final.xyz> <n_images> <nep.txt>

# Example
gpumdkit.sh -calc neb init.xyz fin.xyz 9 nep.txt

# Direct Python execution (also works)
python Scripts/calculators/neb_calculation.py init.xyz fin.xyz 9 nep.txt

# Alternative with NepTrainKit:
python Scripts/calculators/neb_calculation_neptrain.py init.xyz fin.xyz 9 nep.txt
```

## Common Workflows

### Ionic Transport Analysis

Run each path from the working directory where its `msd.out` is written.

Path A derives a four-column `msd.out` from an extxyz trajectory and supports
the MSD plot:

```bash
gpumdkit.sh -calc msd dump.xyz Li 10
gpumdkit.sh -plt msd
```

Path B uses GPUMD `compute_msd`; for one selected group its `msd.out` has seven
columns (time, `MSD_x/y/z`, `SDC_x/y/z`) and supports the SDC plots:

```bash
# Use existing GPUMD compute_msd output
gpumdkit.sh -plt msd
gpumdkit.sh -plt sdc
gpumdkit.sh -plt msd_sdc
```

Choose one path; when comparing implementations, keep separate results to avoid
overwriting `msd.out`. From the same working
directory, either path can supply the validated `msd.out` required for ionic
conductivity:

```bash
# Conductivity also needs thermo.out, model.xyz, and optionally run.in
gpumdkit.sh -calc ionic-cond Li 1
```

`compute_sdc` writes a separate `sdc.out` for `gpumdkit.sh -plt vac`; it is not
the `msd.out` consumed by the SDC plotters.

### Descriptor Analysis
```bash
# 1. Calculate descriptors
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li
# 2. Visualize with PCA
gpumdkit.sh -plt des pca descriptors.npy
# 3. Or visualize with UMAP
gpumdkit.sh -plt des umap descriptors.npy
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
| xrd | `numpy` |
| phonon | `calorine`, `phonopy` |

Install dependencies:
```bash
pip install numpy ase scipy calorine matplotlib tqdm
pip3 install git+https://github.com/MoseyQAQ/ferrodispcalc.git
```

## Detailed Documentation

See `${GPUMDkit_path}/docs/tutorials/en/calculator_scripts.md` or the Chinese counterpart for the user-facing guide.
