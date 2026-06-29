<div align="center">
  <h1>🔍 Analyzer Scripts</h1>
  <p style="text-align: justify;">Analyzer scripts provide structure statistics, distance checks, filters, and composition analysis for extxyz datasets.</p>
</div>

**Script Location:** `Scripts/analyzer/`

## Overview

| Task | Command | Purpose |
|------|---------|---------|
| Composition | `gpumdkit.sh -analyze_comp train.xyz` | group structures by chemical composition |
| Chemical species | `gpumdkit.sh -chem_species train.xyz` | list species in a file |
| Property range | `gpumdkit.sh -range train.xyz force` | inspect energy/force/virial range |
| Minimum distance | `gpumdkit.sh -min_dist dump.xyz` | fast distance check without PBC |
| Minimum distance with PBC | `gpumdkit.sh -min_dist_pbc dump.xyz` | accurate distance check with PBC |
| Charge balance | `gpumdkit.sh -cbc train.xyz` | oxidation-state balance check |
| Distance filter | `gpumdkit.sh -filter_dist_pbc dump.xyz 1.0` | remove structures with short contacts |
| Box filter | `gpumdkit.sh -filter_box dump.xyz 13` | remove structures with too-large box edges |
| Property filter | `gpumdkit.sh -filter_value train.xyz force 20` | filter by energy/force/virial threshold |
| Pair-distance range | `gpumdkit.sh -filter_range dump.xyz Li Li 1.8 2.0` | extract structures by pair distance |

## Composition Analysis

```bash
gpumdkit.sh -analyze_comp train.xyz
```

Example output:

```text
Index    Compositions           N atoms      Count
---------------------------------------------------
1        Li56O96Zr16La24        192          51
---------------------------------------------------
Enter index to export (e.g., '1,2', '2-3', 'all'), or press Enter to skip:
```

This is useful when `train.xyz` contains structures from different systems or different cell sizes. You can export a subset by composition from the interactive prompt.

## Property Range Analysis

```bash
gpumdkit.sh -range train.xyz force
gpumdkit.sh -range train.xyz energy
gpumdkit.sh -range train.xyz virial
gpumdkit.sh -range train.xyz force hist
```

Use `hist` to show a histogram of the selected property.

<div align="center">
  <img src="../../Gallery/range_force.png" alt="Force range histogram" width="52%" />
</div>

## Minimum Distance Checks

Fast check without PBC:

```bash
gpumdkit.sh -min_dist dump.xyz
```

Accurate check with PBC:

```bash
gpumdkit.sh -min_dist_pbc dump.xyz
```

Typical output:

```text
Minimum interatomic distances (with PBC):
+---------------------------+
| Atom Pair |  Distance (A) |
+---------------------------+
|   Li-Li   |     1.696     |
|   Li-O    |     1.587     |
|   O-O     |     2.355     |
+---------------------------+
Overall min_distance: 1.587 A
```

The PBC version considers periodic boundary conditions.

## Filtering Structures

### Filter by Minimum Distance

```bash
gpumdkit.sh -filter_dist_pbc dump.xyz 1.0
```

This removes structures with any interatomic distance below `1.0 A`.

### Filter by Box Size

```bash
gpumdkit.sh -filter_box dump.xyz 13
```

This filters structures by box-edge length.

### Filter by Property Value

```bash
gpumdkit.sh -filter_value train.xyz force 20
gpumdkit.sh -filter_value train.xyz energy 5
```

### Filter by Pair-Distance Range

```bash
gpumdkit.sh -filter_range dump.xyz Li Li 1.8 2.0
```

This extracts structures where the Li-Li minimum distance is between `1.8 A` and `2.0 A`.

## Charge Balance Check

```bash
gpumdkit.sh -cbc train.xyz
```

Outputs may include:

- `balanced.xyz`
- `unbalanced.xyz`
- `indices.txt`

This is intended for systems where common oxidation states are meaningful.

## Example Commands

Property and distance checks:

```bash
gpumdkit.sh -analyze_comp train.xyz
gpumdkit.sh -range train.xyz force hist
gpumdkit.sh -range train.xyz energy hist
gpumdkit.sh -min_dist_pbc train.xyz
```

Filtering examples:

```bash
gpumdkit.sh -filter_dist_pbc dump.xyz 1.0
gpumdkit.sh -filter_box filtered_dist.xyz 13
```
