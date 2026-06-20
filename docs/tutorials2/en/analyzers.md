# Analyzers Guide

<div align="center">
  <p>
    <a href="../zh/analyzers.md">中文</a> | <strong>English</strong>
  </p>
</div>

This guide covers all analysis tools in GPUMDkit for structure validation, filtering, and quality control.

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

## Interactive Mode

```bash
gpumdkit.sh
# Select: 5) Analyzer
```

You'll see:

```
+------------------------------------------------------+
|                    ANALYZER TOOLS                    |
+------------------------------------------------------+
| 501) Analyze composition of extxyz                   |
| 502) Find outliers of extxyz                         |
| 503) Analyze chemical species of extxyz              |
| 504) Check charge balance of extxyz                  |
| 505) Analyze energy/force/virial range               |
| 506) Filter structures by minimum distance           |
| 507) Get minimum interatomic distance                |
| 508) Probability density analysis                    |
+------------------------------------------------------+
```

## Command Reference

### Composition Analysis

Analyze and group structures by chemical composition.

```bash
gpumdkit.sh -analyze_comp train.xyz
```

**Output:** Table of unique compositions with atom counts and structure counts. Interactive selection to export specific compositions.

### Property Range Analysis

Compute the range (min, max) of energy, force, or virial in an extxyz file.

```bash
gpumdkit.sh -range <filename> <property> [hist]

# Examples
gpumdkit.sh -range train.xyz energy
gpumdkit.sh -range train.xyz force
gpumdkit.sh -range train.xyz virial
gpumdkit.sh -range train.xyz force hist    # With histogram
```

**Output:** Min/max values. Optional histogram plot with `hist` argument.

### Minimum Distance

Calculate minimum interatomic distance for each element pair.

```bash
# Fast calculation (no PBC)
gpumdkit.sh -min_dist dump.xyz

# Accurate calculation (with PBC)
gpumdkit.sh -min_dist_pbc dump.xyz
```

**Output:** Table of minimum distances for all element pairs

### Charge Balance Check

Check oxidation-state balance for all structures.

```bash
gpumdkit.sh -cbc train.xyz
```

**Output:**
- `balanced.xyz` - Charge-balanced structures
- `unbalanced.xyz` - Unbalanced structures
- `indices.txt` - Summary

### Chemical Species Identification

List all unique chemical elements in the file.

```bash
# Interactive mode
gpumdkit.sh  # Select: 5) Analyzer -> 503
```

**Output:** Sorted list of all elements

### Structure Filtering

#### By Minimum Distance (No PBC)

```bash
gpumdkit.sh  # Select: 5) Analyzer -> 506
```

**Input:** extxyz file, distance threshold
**Output:** `filtered_<file>.xyz`, `filtered_out_<file>.xyz`

#### By Minimum Distance (With PBC)

```bash
gpumdkit.sh  # Select: 5) Analyzer -> 506b
```

More accurate for periodic systems.

#### By Element-Pair Distance Range

```bash
gpumdkit.sh -filter_range <file> <element1> <element2> <min_dist> <max_dist>

# Example: Filter Li-Li distances between 1.9 and 2.0 Å
gpumdkit.sh -filter_range dump.xyz Li Li 1.9 2.0
```

**Output:** `filtered_<elem1>_<elem2>_<min>_<max>.xyz`

#### By Box Size

```bash
gpumdkit.sh -filter_box <file> <edge_limit>

# Example: Filter structures with box edge > 20 Å
gpumdkit.sh -filter_box dump.xyz 20
```

**Output:** `filtered_by_box.xyz`

#### By Property Value

```bash
gpumdkit.sh -filter_value <file> <property> <threshold>

# Example: Keep structures with force < 20 eV/Å
gpumdkit.sh -filter_value train.xyz force 20
```

**Output:** `filtered.xyz`

### Outlier Detection

Find outlier structures based on RMSE thresholds.

```bash
gpumdkit.sh  # Select: 5) Analyzer -> 502
```

**Required files:**
- `train.xyz`
- `energy_train.out`
- `force_train.out`
- `stress_train.out`

**Output:** `selected.xyz` (high-error), `remained.xyz` (low-error)

### Probability Density Analysis

Calculate 3D probability density of mobile ions from AIMD trajectory.

```bash
gpumdkit.sh  # Select: 5) Analyzer -> 508
```

**Parameters:**
- Reference structure (POSCAR)
- Trajectory file (extxyz)
- Mobile species (e.g., Li)
- Grid interval (e.g., 0.25 Å)

**Output:** `probability_density_<interval>.vasp` (CHGCAR format)

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

## CLI Reference

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
| `-time` | Time monitoring | `gpumdkit.sh -time <gpumd\|nep>` |

## Dependencies

| Tool | Required Packages |
|------|------------------|
| All Python scripts | `ase`, `numpy` |
| Distance calculations | `scipy` |
| Charge balance | `pymatgen`, `tqdm` |
| Property range | `matplotlib` |
| Probability density | `pymatgen` |

## See Also

- [Format Conversion](format_conversion.md) - Convert file formats
- [Calculators](calculators.md) - Compute properties
- [Visualization](visualization.md) - Plot results
