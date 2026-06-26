---
name: gpumdkit-analyzers
description: >
  Use when analyzing, validating, or filtering molecular dynamics structures.
  Provides composition analysis, outlier detection, chemical species identification,
  charge balance check, distance calculations, property filtering, and probability density analysis.
  Use when user asks about: structure analysis, data quality, minimum distance, composition,
  charge balance, outlier detection, or structure filtering.
allowed-tools: Bash(gpumdkit *) Bash(python3 *)
---

# GPUMDkit Analyzers

## Available Tools

| Tool | Command | Description |
|------|---------|-------------|
| Composition Analysis | `-analyze_comp` | Group structures by composition |
| Outlier Detection | Menu 502 | Find high-error structures |
| Chemical Species | Menu 503 | List unique elements |
| Charge Balance | `-cbc` | Check oxidation-state balance |
| Property Range | `-range` | Energy/force/virial statistics |
| Distance Filter | Menu 506 | Filter by minimum distance (no PBC) |
| Distance Filter (PBC) | Menu 506b | Filter by minimum distance (with PBC) |
| Minimum Distance | `-min_dist` | Calculate min distances (no PBC) |
| Minimum Distance (PBC) | `-min_dist_pbc` | Calculate min distances (with PBC) |
| Probability Density | Menu 508 | 3D diffusion channel analysis |

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

### Minimum Distance
```bash
# Fast calculation (no PBC)
gpumdkit.sh -min_dist dump.xyz

# Accurate calculation (with PBC)
gpumdkit.sh -min_dist_pbc dump.xyz

# Output: Table of minimum distances for all element pairs
```

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
# List unique elements (interactive mode)
gpumdkit.sh  # Select: 5) Analyzer -> 503

# Output: Sorted list of all elements in the file
```

### Structure Filtering

#### By Minimum Distance (No PBC)
```bash
gpumdkit.sh  # Select: 5) Analyzer -> 506
# Input: extxyz file, distance threshold
# Output: filtered_<file>.xyz, filtered_out_<file>.xyz
```

#### By Minimum Distance (With PBC)
```bash
gpumdkit.sh  # Select: 5) Analyzer -> 506b
# More accurate for periodic systems
```

#### By Element-Pair Distance Range
```bash
# Filter structures where Li-Li distance is between 1.9 and 2.0 Angstrom
gpumdkit.sh -filter_range dump.xyz Li Li 1.9 2.0

# Output: filtered_<elem1>_<elem2>_<min>_<max>.xyz
```

#### By Box Size
```bash
# Filter structures with box edge > 20 Angstrom
gpumdkit.sh -filter_box dump.xyz 20

# Output: filtered_by_box.xyz
```

#### By Property Value
```bash
# Keep structures with force < 20 eV/Angstrom
gpumdkit.sh -filter_value train.xyz force 20

# Output: filtered.xyz
```

### Outlier Detection
```bash
# Find outlier structures based on RMSE
gpumdkit.sh  # Select: 5) Analyzer -> 502

# Required files in current directory:
# - train.xyz
# - energy_train.out
# - force_train.out
# - stress_train.out

# Output: selected.xyz (high-error), remained.xyz (low-error)
```

### Probability Density Analysis
```bash
# Calculate 3D probability density of mobile ions
gpumdkit.sh  # Select: 5) Analyzer -> 508

# Parameters (interactive):
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
```bash
# 1. Filter by distance
gpumdkit.sh  # Select: 5) Analyzer -> 506

# 2. Filter by box size
gpumdkit.sh -filter_box filtered.xyz 20

# 3. Filter by force value
gpumdkit.sh -filter_value filtered_by_box.xyz force 15
```

### Diffusion Channel Analysis
```bash
# 1. Calculate probability density
gpumdkit.sh  # Select: 5) Analyzer -> 508

# 2. Visualize with VESTA or similar
# Open probability_density_0.25.vasp
```

## Additional Tools

### Time Monitoring
```bash
# Monitor GPUMD progress
gpumdkit.sh -time gpumd

# Monitor NEP training progress
gpumdkit.sh -time nep
```

### Volume Analysis
```bash
# Get average volume per temperature
python Scripts/analyzer/get_volume.py
# Requires: */K/thermo.out directories
```

## CLI Flag Reference

| Flag | Description | Syntax |
|------|-------------|--------|
| `-analyze_comp` | Composition analysis | `gpumdkit.sh -analyze_comp <file>` |
| `-range` | Property range | `gpumdkit.sh -range <file> <prop>` |
| `-min_dist` | Min distance (no PBC) | `gpumdkit.sh -min_dist <file>` |
| `-min_dist_pbc` | Min distance (PBC) | `gpumdkit.sh -min_dist_pbc <file>` |
| `-cbc` | Charge balance | `gpumdkit.sh -cbc <file>` |
| `-filter_range` | Filter by distance range | `gpumdkit.sh -filter_range <file> <e1> <e2> <min> <max>` |
| `-filter_box` | Filter by box size | `gpumdkit.sh -filter_box <file> <limit>` |
| `-filter_value` | Filter by property | `gpumdkit.sh -filter_value <file> <prop> <thresh>` |
| `-time` | Time monitoring | `gpumdkit.sh -time <gpumd|nep>` |

## Dependencies

| Tool | Required Packages |
|------|------------------|
| All Python scripts | `ase`, `numpy` |
| Distance calculations | `scipy` |
| Charge balance | `pymatgen`, `tqdm` |
| Property range | `matplotlib` |
| Probability density | `pymatgen` |

## Detailed Documentation

See [analyzer.md](../../docs/tutorials/analyzer.md) for comprehensive guide.
