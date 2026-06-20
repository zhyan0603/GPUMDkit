# Format Conversion Guide

<div align="center">
  <p>
    <a href="../zh/format_conversion.md">中文</a> | <strong>English</strong>
  </p>
</div>

This guide covers all format conversion tools in GPUMDkit.

## Supported Formats

| Format | Extensions | Description |
|--------|-----------|-------------|
| VASP | POSCAR, OUTCAR, XDATCAR, vasprun.xml | Vienna Ab initio Simulation Package |
| LAMMPS | .data, dump.* | Large-scale Atomic/Molecular Massively Parallel Simulator |
| CP2K | .log, pos.xyz, frc.xyz, cell.cell | Quantum chemistry and solid state physics |
| ABACUS | running_scf.log, running_md.log | Atomic-orbital Based Ab-initio Computation |
| CIF | .cif | Crystallographic Information File |
| MTP | .cfg | Moment Tensor Potential format |
| ASE | .traj | Atomic Simulation Environment trajectory |
| extxyz | .xyz | Extended XYZ (primary working format) |

## Interactive Mode

Launch the interactive menu and select option 1:

```bash
gpumdkit.sh
# Select: 1) Format Conversion
```

You'll see:

```
+-------------------------------------------------------------+
|                   FORMAT CONVERSION TOOLS                   |
+-------------------------------------------------------------+
| 101) VASP to extxyz            106) Add group labels        |
| 102) MTP to extxyz             107) Add weight to extxyz    |
| 103) CP2K to extxyz            108) Extract frame extxyz    |
| 104) ABACUS to extxyz          109) Clean XYZ info          |
| 105) extxyz to POSCAR          110) Replicate structure     |
+-------------------------------------------------------------+
| out2exyz) OUTCAR to extxyz     xdat2exyz) XDATCAR to extxyz |
| pos2exyz) POSCAR to extxyz     pos2lmp)   POSCAR to LAMMPS  |
| cif2pos)  CIF to POSCAR        lmp2exyz)  LAMMPS to extxyz  |
| cif2exyz) CIF to extxyz        traj2exyz) ASE traj to extxyz|
+-------------------------------------------------------------+
| 000) Return to main menu                                    |
+-------------------------------------------------------------+
```

## Command-Line Mode

### VASP Conversions

#### OUTCAR to extxyz

```bash
# Shell version (supports both OUTCAR and vasprun.xml)
gpumdkit.sh -out2xyz <directory>

# Python version (OUTCAR only)
gpumdkit.sh -out2exyz <directory>

# Example
gpumdkit.sh -out2xyz ./vasp_results/
```

#### XDATCAR to extxyz

```bash
gpumdkit.sh -xdat2exyz XDATCAR output.xyz
```

#### POSCAR to extxyz

```bash
gpumdkit.sh -pos2exyz POSCAR model.xyz
```

#### extxyz to POSCAR

```bash
# Converts ALL frames to separate POSCAR files
gpumdkit.sh -exyz2pos structures.xyz
```

### LAMMPS Conversions

#### LAMMPS dump to extxyz

```bash
# IMPORTANT: Element symbols must match atom type IDs in dump file
gpumdkit.sh -lmp2exyz dump.lammpstrj Li Y Cl

# Example
gpumdkit.sh -lmp2exyz trajectory.dump Na Cl
```

#### POSCAR to LAMMPS data

```bash
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

## Structure Manipulation

### Add Group Labels

Group labels are required for GPUMD/NEP to identify atom types:

```bash
gpumdkit.sh -addgroup POSCAR <element1> <element2> ...

# Example
gpumdkit.sh -addgroup POSCAR Li Y Cl
```

This creates `model.xyz` with group labels assigned to each atom.

### Add Weights

Add weight labels to structures for weighted training:

```bash
gpumdkit.sh -addweight input.xyz output.xyz <weight>

# Example
gpumdkit.sh -addweight train.xyz weighted.xyz 2.0
```

### Replicate Structure

Create supercells from unit cells:

```bash
# By replication factors
gpumdkit.sh -replicate POSCAR supercell.vasp 2 2 2

# By target atom count
gpumdkit.sh -replicate POSCAR supercell.vasp 256
```

### Extract Frame

Extract a specific frame from a trajectory:

```bash
# 0-based index
gpumdkit.sh -get_frame trajectory.xyz 1000
```

### Clean XYZ Metadata

Remove extra metadata from extxyz files:

```bash
gpumdkit.sh -clean_xyz input.xyz clean.xyz
```

## Complete CLI Reference

| Flag | Description | Syntax |
|------|-------------|--------|
| `-out2xyz` | OUTCAR to extxyz (shell) | `gpumdkit.sh -out2xyz <dir>` |
| `-out2exyz` | OUTCAR to extxyz (python) | `gpumdkit.sh -out2exyz <dir>` |
| `-pos2exyz` | POSCAR to extxyz | `gpumdkit.sh -pos2exyz <poscar> <xyz>` |
| `-exyz2pos` | extxyz to POSCAR | `gpumdkit.sh -exyz2pos <xyz>` |
| `-pos2lmp` | POSCAR to LAMMPS data | `gpumdkit.sh -pos2lmp <poscar> <lmp>` |
| `-lmp2exyz` | LAMMPS dump to extxyz | `gpumdkit.sh -lmp2exyz <dump> <elem...>` |
| `-cif2pos` | CIF to POSCAR | `gpumdkit.sh -cif2pos <cif> <output>` |
| `-cif2exyz` | CIF to extxyz | `gpumdkit.sh -cif2exyz <cif> <output>` |
| `-xdat2exyz` | XDATCAR to extxyz | `gpumdkit.sh -xdat2exyz XDATCAR dump.xyz` |
| `-traj2exyz` | ASE traj to extxyz | `gpumdkit.sh -traj2exyz <traj> <xyz>` |
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
# Create 2x2x2 supercell
gpumdkit.sh -replicate POSCAR supercell_222.vasp 2 2 2

# Create supercell with ~256 atoms
gpumdkit.sh -replicate POSCAR supercell_256.vasp 256
```

## Important Notes

1. **extxyz is the primary format**: Most GPUMDkit tools work with extxyz files
2. **Group labels are essential**: Required for GPUMD/NEP to identify atom types
3. **Frame indexing is 0-based**: First frame is index 0
4. **Element ordering matters for LAMMPS**: Must match atom type IDs in dump file
5. **exyz2pos exports all frames**: Creates separate POSCAR for each frame

## Dependencies

Most Python scripts require:
- `ase` (Atomic Simulation Environment)
- `numpy`

## Troubleshooting

### Issue: Missing element information
**Solution**: Ensure POSCAR has correct element names in the header

### Issue: Incorrect atom count
**Solution**: Check if POSCAR has selective dynamics enabled

### Issue: LAMMPS conversion wrong species
**Solution**: Verify element order matches atom type IDs in dump file

## See Also

- [Calculators](calculators.md) - Compute material properties
- [Analyzers](analyzers.md) - Structure analysis
- [NEP Training Guide](nep_training.md) - Complete training workflow
