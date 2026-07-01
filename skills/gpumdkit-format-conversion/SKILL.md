---
name: gpumdkit-format-conversion
description: >
  Use when converting structure files between computational materials science formats.
  Supports VASP (POSCAR/OUTCAR/XDATCAR), LAMMPS, CP2K, ABACUS, CIF, MTP, ASE trajectory, and extxyz.
  Use when user asks about: format conversion, file conversion, structure conversion,
  POSCAR to xyz, OUTCAR to extxyz, LAMMPS dump conversion, or adding group labels.
allowed-tools: Bash(gpumdkit *) Bash(python3 *)
---

# GPUMDkit Format Conversion

## Supported Formats

| Format | Extensions | Description |
|--------|-----------|-------------|
| VASP | POSCAR, OUTCAR, XDATCAR | Vienna Ab initio Simulation Package |
| LAMMPS | .data, dump.* | Large-scale Atomic/Molecular Massively Parallel Simulator |
| CP2K | .log, pos.xyz, frc.xyz, cell.cell | Quantum chemistry and solid state physics |
| ABACUS | running_scf.log, running_md.log | Atomic-orbital Based Ab-initio Computation at UStc |
| CIF | .cif | Crystallographic Information File |
| MTP | .cfg | Moment Tensor Potential format |
| ASE | .traj | Atomic Simulation Environment trajectory |
| extxyz | .xyz | Extended XYZ (primary working format) |

## Command Reference

### VASP Conversions

```bash
# OUTCAR to extxyz (directory, shell version)
gpumdkit.sh -out2xyz <directory>

# OUTCAR to extxyz (Python version)
gpumdkit.sh -out2exyz <directory>

# XDATCAR to extxyz
gpumdkit.sh -xdat2exyz XDATCAR output.xyz

# POSCAR to extxyz
gpumdkit.sh -pos2exyz POSCAR model.xyz

# extxyz to POSCAR (all frames)
gpumdkit.sh -exyz2pos structures.xyz
```

### LAMMPS Conversions

```bash
# LAMMPS dump to extxyz
# IMPORTANT: Element symbols must match atom type IDs in dump file
gpumdkit.sh -lmp2exyz dump.lammpstrj Li Y Cl

# POSCAR to LAMMPS data
gpumdkit.sh -pos2lmp POSCAR lammps.data
```

### CIF Conversions

```bash
# CIF to POSCAR
gpumdkit.sh -cif2pos input.cif POSCAR.vasp

# CIF to extxyz
gpumdkit.sh -cif2exyz input.cif model.xyz
```

### Other Conversions

```bash
# ASE trajectory to extxyz
gpumdkit.sh -traj2exyz input.traj output.xyz
```

### dp2xyz

Converts DeepMD npy datasets to extxyz format. Recursively scans a directory for datasets containing `type.raw`, `type_map.raw`, and `set.000/`.

**Usage:**
```
gpumdkit.sh -dp2xyz database train.xyz
```
Dependencies: `dpdata`, `ase`. Requires `conda activate gpumd` environment.
Author: Denan LI (lidenan@westlake.edu.cn)

### Structure Manipulation

```bash
# Add group labels (required for GPUMD/NEP)
gpumdkit.sh -addgroup POSCAR Li Y Cl

# Add weights to structures
gpumdkit.sh -addweight input.xyz output.xyz 2.0

# Replicate structure (by factors)
gpumdkit.sh -replicate POSCAR supercell.vasp 2 2 2

# Replicate structure (by target atom count)
gpumdkit.sh -replicate POSCAR supercell.vasp 256

# Extract specific frame (0-based index)
gpumdkit.sh -get_frame trajectory.xyz 1000

# Clean extxyz metadata
gpumdkit.sh -clean_xyz input.xyz clean.xyz
```

## Detailed CLI Flag Reference

| CLI Flag | Conversion | Syntax |
|----------|-----------|--------|
| `-out2xyz` | OUTCAR -> extxyz (shell) | `gpumdkit.sh -out2xyz <dir>` |
| `-out2exyz` | OUTCAR -> extxyz (python) | `gpumdkit.sh -out2exyz <dir>` |
| `-pos2exyz` | POSCAR -> extxyz | `gpumdkit.sh -pos2exyz <poscar> <xyz>` |
| `-exyz2pos` | extxyz -> POSCAR | `gpumdkit.sh -exyz2pos <xyz>` |
| `-pos2lmp` | POSCAR -> LAMMPS data | `gpumdkit.sh -pos2lmp <poscar> <lmp>` |
| `-lmp2exyz` | LAMMPS dump -> extxyz | `gpumdkit.sh -lmp2exyz <dump> <elem...>` |
| `-cif2pos` | CIF -> POSCAR | `gpumdkit.sh -cif2pos <cif> <output>` |
| `-cif2exyz` | CIF -> extxyz | `gpumdkit.sh -cif2exyz <cif> <output>` |
| `-xdat2exyz` | XDATCAR -> extxyz | `gpumdkit.sh -xdat2exyz XDATCAR dump.xyz` |
| `-traj2exyz` | ASE traj -> extxyz | `gpumdkit.sh -traj2exyz <traj> <xyz>` |
| `-dp2xyz` | DeepMD npy → extxyz (via dpdata) | `gpumdkit.sh -dp2xyz <input_dir/> [output.xyz]` |
| `-addgroup` | Add group labels | `gpumdkit.sh -addgroup <poscar> <elem...>` |
| `-addweight` | Add weight | `gpumdkit.sh -addweight <in> <out> <weight>` |
| `-replicate` | Replicate structure | `gpumdkit.sh -replicate <in> <out> a b c` |
| `-get_frame` | Extract frame | `gpumdkit.sh -get_frame <xyz> <index>` |
| `-clean_xyz` | Clean extxyz info | `gpumdkit.sh -clean_xyz <in> <out>` |

## Examples

### Example 1: Convert VASP MD Output
```bash
# Convert all OUTCAR files in current directory
gpumdkit.sh -out2xyz .

# Add group labels for NEP training
gpumdkit.sh -addgroup POSCAR Pb Ti O

# Result: model.xyz ready for NEP training
```

### Example 2: Prepare LAMMPS Simulation
```bash
# Convert POSCAR to LAMMPS data
gpumdkit.sh -pos2lmp POSCAR system.data

# After LAMMPS simulation, convert dump back
gpumdkit.sh -lmp2exyz dump.lammpstrj Li P S
```

### Example 3: Batch Conversion
```bash
# Convert multiple OUTCAR files
for dir in run_*; do
    gpumdkit.sh -out2xyz "$dir"
    mv "$dir"/model.xyz "$dir"/trajectory.xyz
done
```

### Example 4: Structure Replication
```bash
# Replicate to 2x2x2 supercell
gpumdkit.sh -replicate POSCAR supercell_222.vasp 2 2 2

# Replicate to target ~256 atoms
gpumdkit.sh -replicate POSCAR supercell_256.vasp 256
```

## Notes

1. **extxyz is the primary format**: Most GPUMDkit tools work with extxyz files
2. **Group labels are essential**: Required for GPUMD/NEP to identify atom types
3. **Frame indexing is 0-based**: First frame is index 0
4. **Element ordering matters for LAMMPS**: Must match atom type IDs in dump file
5. **exyz2pos exports all frames**: Creates separate POSCAR for each frame

## Dependencies

Most Python scripts require:
- `ase` (Atomic Simulation Environment)
- `numpy`

## Detailed Documentation

See [format_conversion.md](../../docs/tutorials/en/format_conversion.md) for comprehensive guide.
