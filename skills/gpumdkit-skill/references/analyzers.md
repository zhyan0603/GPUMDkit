# Analyzers

## Contents

- Available tools
- Command reference
- Common workflows
- CLI flags and dependencies

## Available Tools

| Tool | Command | Description |
|------|---------|-------------|
| Composition Analysis | `-analyze_comp` | Group structures by composition |
| Outlier Detection | Menu 502 | Find high-error structures |
| Chemical Species | `-chem_species` or Menu 503 | List unique elements |
| Charge Balance | `-cbc` | Check oxidation-state balance |
| Property Range | `-range` | Energy/force/virial statistics |
| Distance Filter | `-filter_dist` or Menu 506 | Filter by minimum distance (no PBC) |
| Distance Filter (PBC) | `-filter_dist_pbc` | Filter by minimum distance with PBC |
| Minimum Distance | `-min_dist` | Calculate min distances (no PBC) |
| Minimum Distance (PBC) | `-min_dist_pbc` | Calculate min distances (with PBC) |
| Probability Density | `-pda` or Menu 508 | 3D diffusion channel analysis |
| Energy Reference Shifting | `-shift_energy` or Menu 509 | Interactively shift extxyz energy references |

## Command Reference

### Composition Analysis
```bash
# Analyze composition of extxyz file
gpumdkit.sh -analyze_comp train.xyz

# Output: Table of unique compositions with counts
# Interactive: Select compositions to export as separate files
```

### Property Range Analysis
```bash
# Analyze energy range
gpumdkit.sh -range train.xyz energy

# Analyze force range
gpumdkit.sh -range train.xyz force

# Analyze virial range
gpumdkit.sh -range train.xyz virial

# With histogram plot
gpumdkit.sh -range train.xyz force hist
```

Energy and virial are read from `energy=`/`virial=` in either letter case,
including keys that ASE 3.23+ stores in the calculator. Frames without the
requested property are skipped with a note; the command fails only when no
frame has it. Energy is the total energy in the file (eV); force ranges use the
per-atom force magnitude; virial ranges use every stored component.

### Energy Reference Shifting
```bash
# Open the interactive energy reference shifter
gpumdkit.sh -shift_energy

# The same tool is available from menu 5 -> 509
```

`align_energy_reference.py` preserves atom lines and all other header fields
while updating `energy=` values. It supports reference-group alignment, zero
atomic baseline alignment, and DFT-to-NEP alignment. The latter requires the
`calorine` package and a compatible NEP model file.

### Minimum Distance
```bash
# Fast calculation (no PBC)
gpumdkit.sh -min_dist dump.xyz

# Accurate calculation (with PBC)
gpumdkit.sh -min_dist_pbc dump.xyz

# Output: Table of minimum distances for all element pairs
```

Both variants collect elements across all frames in first-appearance order and skip pairs absent from a frame.

### Charge Balance Check
```bash
# Check charge balance
gpumdkit.sh -cbc train.xyz

# Output:
# - balanced.xyz (charge-balanced structures)
# - unbalanced.xyz (unbalanced structures)
# - indices.txt (summary)
```

### Chemical Species Identification
```bash
# List unique elements
gpumdkit.sh -chem_species train.xyz

# Interactive mode also works:
# gpumdkit.sh  # Select: 5) Analyzer -> 503

# Output: Sorted list of all elements in the file
```

### Structure Filtering

#### By Minimum Distance (No PBC)
```bash
gpumdkit.sh -filter_dist dump.xyz <min_distance>
# Output: filtered_<input_name>, filtered_out_<input_name>
# e.g. filtered_dump.xyz and filtered_out_dump.xyz
```

#### By Minimum Distance (With PBC)
```bash
gpumdkit.sh -filter_dist_pbc dump.xyz <min_distance>
# More accurate for periodic systems; same output names as -filter_dist
```

#### By Element-Pair Distance Range
```bash
# Filter structures where Li-Li distance is between 1.9 and 2.0 Angstrom
gpumdkit.sh -filter_range dump.xyz Li Li 1.9 2.0

# Output: filtered_<elem1>_<elem2>_<min>_<max>.xyz
```

#### By Box Size
```bash
# Discard structures with any box edge longer than <max_edge> Angstrom
gpumdkit.sh -filter_box dump.xyz <max_edge>

# Output: filtered_by_box.xyz
```

#### By Property Value
```bash
# Keep structures whose values are all <= the threshold
gpumdkit.sh -filter_value train.xyz force <max_force_component>

# Output: filtered.xyz
```

A structure is kept when every value is at most the threshold: each force
component in magnitude (eV/Angstrom) for `force`, the total energy in the file
(eV, not per atom) for `energy`, and each virial component for `virial`.
Energy and virial are read in either letter case, including keys that ASE 3.23+
stores in the calculator. Frames without the chosen property cannot be checked
and are kept; the command reports how many. It fails without writing output
when no frame has the property, the property name is unknown, or the threshold
is not a number.

### Outlier Detection
```bash
# Find outlier structures based on RMSE
gpumdkit.sh  # Select: 5) Analyzer -> 502

# Non-interactive threshold mode also works
python3 ${GPUMDkit_path}/Scripts/analyzer/find_outliers.py <energy_rmse_meV_atom> <force_rmse_meV_A> <stress_rmse_GPa>

# Required files in current directory:
# - train.xyz
# - energy_train.out
# - force_train.out
# - stress_train.out

# Output: selected.xyz (high-error), remained.xyz (low-error)
```

The three optional thresholds are target RMSE values for energy (meV/atom),
force (meV/Angstrom), and stress (GPa). For each property the script removes
the highest-error structures one at a time until the RMSE of the remaining set
falls below the threshold; removed structures go to `selected.xyz`. Choose the
thresholds with the user.

Structures with a `-1e6` missing-stress marker or a stress component whose absolute value exceeds `1e6` are excluded from RMSE calculation and both output datasets. If none remain, the script reports an error without writing outputs.

### Probability Density Analysis
```bash
# Calculate 3D probability density of mobile ions
gpumdkit.sh -pda LLZO.vasp dump.xyz Li 0.25

# Parameters:
# - Reference structure (POSCAR)
# - Trajectory file (extxyz)
# - Mobile species (e.g., Li)
# - Grid interval (e.g., 0.25 Angstrom)

# Output: probability_density_<interval>.vasp (CHGCAR format)
# Can be visualized with VESTA or similar software
```

## Common Workflows

### Data Quality Check
```bash
# 1. Check composition
gpumdkit.sh -analyze_comp train.xyz

# 2. Check minimum distances
gpumdkit.sh -min_dist_pbc train.xyz

# 3. Check force range
gpumdkit.sh -range train.xyz force

# 4. Check for outliers
gpumdkit.sh  # Select: 5) Analyzer -> 502
```

### Structure Filtering Pipeline

Thresholds are system-specific; take them from the user.

```bash
# 1. Filter by distance
gpumdkit.sh -filter_dist_pbc dump.xyz <min_distance>

# 2. Filter by box size
gpumdkit.sh -filter_box filtered_dump.xyz <max_edge>

# 3. Filter by force value
gpumdkit.sh -filter_value filtered_by_box.xyz force <max_force_component>
```

### Diffusion Channel Analysis
```bash
# 1. Calculate probability density
gpumdkit.sh -pda LLZO.vasp dump.xyz Li 0.25

# 2. Visualize with VESTA or similar
# Open probability_density_0.25.vasp
```

## Additional Tools

### Time Monitoring
```bash
# Monitor GPUMD progress
gpumdkit.sh -time gpumd

# Monitor a WPE-enabled GPUMD run
gpumdkit.sh -time wpe

# Monitor NEP training progress
gpumdkit.sh -time nep
```

## CLI Flag Reference

| Flag | Description | Syntax |
|------|-------------|--------|
| `-analyze_comp` | Composition analysis | `gpumdkit.sh -analyze_comp <file>` |
| `-range` | Property range | `gpumdkit.sh -range <file> <prop>` |
| `-shift_energy` | Interactive energy reference shifting | `gpumdkit.sh -shift_energy` |
| `-min_dist` | Min distance (no PBC) | `gpumdkit.sh -min_dist <file>` |
| `-min_dist_pbc` | Min distance (PBC) | `gpumdkit.sh -min_dist_pbc <file>` |
| `-chem_species` | Chemical species list | `gpumdkit.sh -chem_species <file>` |
| `-cbc` | Charge balance | `gpumdkit.sh -cbc <file>` |
| `-filter_dist` | Filter by minimum distance | `gpumdkit.sh -filter_dist <file> <min_dist>` |
| `-filter_dist_pbc` | Filter by minimum distance with PBC | `gpumdkit.sh -filter_dist_pbc <file> <min_dist>` |
| `-filter_range` | Filter by distance range | `gpumdkit.sh -filter_range <file> <e1> <e2> <min> <max>` |
| `-filter_box` | Filter by box size | `gpumdkit.sh -filter_box <file> <limit>` |
| `-filter_value` | Filter by property | `gpumdkit.sh -filter_value <file> <prop> <thresh>` |
| `-time` | Time monitoring | `gpumdkit.sh -time <gpumd\|wpe\|nep\|gnep>` |

## Dependencies

| Tool | Required Packages |
|------|------------------|
| All Python scripts | `ase`, `numpy` |
| Distance calculations | `scipy` |
| Charge balance | `pymatgen`, `tqdm` |
| Property range | `matplotlib` |
| Probability density | `pymatgen` |
| DFT-to-NEP energy shifting | `calorine` |

## Detailed Documentation

See `${GPUMDkit_path}/docs/tutorials/en/analyzer_scripts.md` or the Chinese counterpart for the user-facing guide.
