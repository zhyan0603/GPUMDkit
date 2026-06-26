<div align="center">
  <h1>Calculators</h1>
  <p>
    <strong>English</strong> | <a href="../zh/calculators.md">简体中文</a>
  </p>
</div>

Compute material properties from molecular dynamics data.

## Available Calculators

| Calculator | Command | Description |
|------------|---------|-------------|
| Ionic Conductivity | `-calc ionic-cond` | Calculate ionic conductivity from MSD |
| NEP Properties | `-calc nep` | Compute energy/force/stress with NEP |
| Descriptors | `-calc des` | Calculate NEP descriptors for analysis |
| DOAS | `-calc doas` | Density of atomistic states |
| MSD | `-calc msd` | Mean square displacement from trajectory |
| Neighbor List | `-calc nlist` | Build neighbor lists for analysis |
| Displacement | `-calc disp` | Calculate atomic displacements |
| Averaged Structure | `-calc avg-struct` | Time-averaged structure from trajectory |
| Octahedral Tilt | `-calc oct-tilt` | Perovskite octahedral tilt angles |
| Polarization | `-calc pol-abo3` | ABO3 local polarization |
| Minimization | `-calc minimize` | Structure relaxation with NEP |
| NEB | Direct Python | Nudged elastic band calculation |

## Command Reference

### Ionic Conductivity

```bash
gpumdkit.sh -calc ionic-cond <element> <charge>

# Examples
gpumdkit.sh -calc ionic-cond Li 1    # Li+
gpumdkit.sh -calc ionic-cond O -2    # O2-
```

Required files: `msd.out`, `thermo.out`, `model.xyz`, `run.in`

### NEP Properties

```bash
gpumdkit.sh -calc nep <input.xyz> <output.xyz> <nep.txt>

# Example
gpumdkit.sh -calc nep structures.xyz predictions.xyz nep.txt
```

### Descriptors

```bash
gpumdkit.sh -calc des <input.xyz> <output.npy> <nep.txt> <element>

# Example
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li
```

Visualize with: `gpumdkit.sh -plt des pca` or `gpumdkit.sh -plt des umap`

### DOAS

```bash
gpumdkit.sh -calc doas <input.xyz> <nep.txt> <output.txt>

# Example
gpumdkit.sh -calc doas structures.xyz nep.txt doas.out
```

### MSD

```bash
gpumdkit.sh -calc msd <trajectory.xyz> <element> <dt_fs> [max_corr_steps]

# Example
gpumdkit.sh -calc msd dump.xyz Li 10
```

Output: `msd.out` (Time/ps, MSD_x, MSD_y, MSD_z)

### Neighbor List

```bash
gpumdkit.sh -calc nlist -i <input> -c <cutoff> -n <num> -C <center> -E <neighbor>

# Example
gpumdkit.sh -calc nlist -i model.xyz -c 4 -n 12 -C Ti -E O
```

### Displacement

```bash
gpumdkit.sh -calc disp -i <trajectory.xyz> -n <neighbor_list> -o <output>

# Example
gpumdkit.sh -calc disp -i movie.xyz -n nl-Pb-O.dat -o displacements.dat
```

### Averaged Structure

```bash
gpumdkit.sh -calc avg-struct -i <trajectory.xyz> -l <fraction> -o <output>

# Example: Average last 20% of frames
gpumdkit.sh -calc avg-struct -i movie.xyz -l 0.2 -o averaged.xyz
```

### Octahedral Tilt

```bash
gpumdkit.sh -calc oct-tilt -i <input.xyz> -n <B-O list> -o <output>

# Example
gpumdkit.sh -calc oct-tilt -i model.xyz -n nl-Ti-O.dat -o tilt.dat
```

### Polarization (ABO3)

```bash
gpumdkit.sh -calc pol-abo3 -i <input.xyz> --nl-ba <B-A list> --nl-bo <B-O list> --bec <Element=charge ...>

# Example
gpumdkit.sh -calc pol-abo3 -i model.xyz --nl-ba nl-Ti-Pb.dat --nl-bo nl-Ti-O.dat --bec Pb=2.5 Ti=4.0 O=-2.0
```

### Minimization

```bash
gpumdkit.sh -calc minimize <structure> <nep.txt> [fmax] [max_steps]

# Example
gpumdkit.sh -calc minimize POSCAR nep.txt 0.01 1000
```

Output: `minimize.xyz`

### NEB

```bash
python Scripts/calculators/neb_calculation.py <initial.xyz> <final.xyz> <n_images> <nep.txt>

# Example
python Scripts/calculators/neb_calculation.py init.xyz fin.xyz 9 nep.txt
```

## Dependencies

| Calculator | Packages |
|------------|----------|
| All | numpy, ase |
| ionic-cond | scipy |
| nep, des, doas, minimize | calorine |
| nlist, disp, oct-tilt, pol-abo3 | ferrodispcalc |
