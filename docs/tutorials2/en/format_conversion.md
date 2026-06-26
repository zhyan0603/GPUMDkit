<div align="center">
  <h1>Format Conversion</h1>
  <p>
    <strong>English</strong> | <a href="../zh/format_conversion.md">简体中文</a>
  </p>
</div>

Convert between computational materials science file formats.

## Supported Formats

| Format | Extensions | Software |
|--------|-----------|----------|
| VASP | POSCAR, OUTCAR, XDATCAR, vasprun.xml | VASP |
| LAMMPS | .data, dump.* | LAMMPS |
| CP2K | .log, pos.xyz, frc.xyz, cell.cell | CP2K |
| ABACUS | running_scf.log, running_md.log | ABACUS |
| CIF | .cif | Various |
| MTP | .cfg | MTP |
| ASE | .traj | ASE |
| extxyz | .xyz | Various (primary format) |

## Quick Reference

```bash
# VASP
gpumdkit.sh -out2xyz <dir>              # OUTCAR -> extxyz
gpumdkit.sh -pos2exyz POSCAR model.xyz  # POSCAR -> extxyz
gpumdkit.sh -exyz2pos structures.xyz    # extxyz -> POSCAR

# LAMMPS
gpumdkit.sh -lmp2exyz dump.lammpstrj Li Y Cl  # LAMMPS -> extxyz
gpumdkit.sh -pos2lmp POSCAR lammps.data        # POSCAR -> LAMMPS

# CIF
gpumdkit.sh -cif2pos input.cif POSCAR.vasp     # CIF -> POSCAR
gpumdkit.sh -cif2exyz input.cif model.xyz      # CIF -> extxyz

# Other
gpumdkit.sh -xdat2exyz XDATCAR output.xyz      # XDATCAR -> extxyz
gpumdkit.sh -traj2exyz input.traj output.xyz   # ASE traj -> extxyz
```

## Structure Manipulation

```bash
# Add group labels (required for NEP)
gpumdkit.sh -addgroup POSCAR Li Y Cl

# Add weights
gpumdkit.sh -addweight input.xyz output.xyz 2.0

# Replicate structure
gpumdkit.sh -replicate POSCAR supercell.vasp 2 2 2
gpumdkit.sh -replicate POSCAR supercell.vasp 256  # by atom count

# Extract frame (0-based)
gpumdkit.sh -get_frame trajectory.xyz 1000

# Clean metadata
gpumdkit.sh -clean_xyz input.xyz clean.xyz
```

## Examples

### Convert VASP Output

```bash
gpumdkit.sh -out2xyz .
gpumdkit.sh -addgroup POSCAR Pb Ti O
```

### Prepare LAMMPS Simulation

```bash
gpumdkit.sh -pos2lmp POSCAR system.data
# After simulation
gpumdkit.sh -lmp2exyz dump.lammpstrj Li P S
```

### Batch Conversion

```bash
for dir in run_*; do
    gpumdkit.sh -out2xyz "$dir"
done
```

## Notes

- extxyz is the primary format for GPUMDkit
- Group labels are required for GPUMD/NEP
- Frame indexing is 0-based
- LAMMPS element order must match atom type IDs

## Dependencies

```bash
pip install ase numpy
```
