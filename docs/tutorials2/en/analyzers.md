<div align="center">
  <h1>Analyzers</h1>
  <p>
    <strong>English</strong> | <a href="../zh/analyzers.md">简体中文</a>
  </p>
</div>

Structure analysis, validation, and filtering tools.

## Available Tools

| Tool | Command | Description |
|------|---------|-------------|
| Composition Analysis | `-analyze_comp` | Group structures by composition |
| Property Range | `-range` | Energy/force/virial statistics |
| Minimum Distance | `-min_dist` | Calculate min distances (no PBC) |
| Minimum Distance (PBC) | `-min_dist_pbc` | Calculate min distances (with PBC) |
| Charge Balance | `-cbc` | Check oxidation-state balance |
| Distance Filter | `-filter_dist_pbc` | Filter by minimum distance |
| Box Filter | `-filter_box` | Filter by box size |
| Value Filter | `-filter_value` | Filter by property value |
| Distance Range | `-filter_range` | Filter by element-pair distance |
| Time Monitor | `-time` | Monitor GPUMD/NEP progress |

## Command Reference

### Composition Analysis

```bash
gpumdkit.sh -analyze_comp train.xyz
```

### Property Range

```bash
gpumdkit.sh -range <file> <property> [hist]

# Examples
gpumdkit.sh -range train.xyz force
gpumdkit.sh -range train.xyz energy
gpumdkit.sh -range train.xyz force hist
```

### Minimum Distance

```bash
# Fast (no PBC)
gpumdkit.sh -min_dist dump.xyz

# Accurate (with PBC)
gpumdkit.sh -min_dist_pbc dump.xyz
```

### Charge Balance

```bash
gpumdkit.sh -cbc train.xyz
```

Output: `balanced.xyz`, `unbalanced.xyz`

### Filtering

```bash
# By minimum distance
gpumdkit.sh -filter_dist_pbc train.xyz 1.0

# By box size
gpumdkit.sh -filter_box dump.xyz 20

# By property value
gpumdkit.sh -filter_value train.xyz force 20

# By element-pair distance range
gpumdkit.sh -filter_range dump.xyz Li Li 1.9 2.0
```

### Time Monitoring

```bash
gpumdkit.sh -time gpumd
gpumdkit.sh -time nep
```

## Dependencies

| Tool | Packages |
|------|----------|
| All | ase, numpy |
| Distance | scipy |
| Charge balance | pymatgen, tqdm |
| Property range | matplotlib |
