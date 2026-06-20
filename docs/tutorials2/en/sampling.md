# Structure Sampling Guide

<div align="center">
  <p>
    <a href="../zh/sampling.md">中文</a> | <strong>English</strong>
  </p>
</div>

This guide covers structure sampling and selection tools in GPUMDkit.

## Available Methods

| Method | Menu | Description |
|--------|------|-------------|
| Uniform/Random | 201 | Select evenly spaced or random frames |
| FPS by PyNEP | 202 | Farthest point sampling (deprecated) |
| FPS by NepTrain | 203 | Farthest point sampling (preferred) |
| Perturbation | 204 | Generate perturbed structures |
| Force Deviation | 205 | Select high-force-deviation structures |

## Interactive Mode

```bash
gpumdkit.sh
# Select: 2) Sample Structures
```

You'll see:

```
+------------------------------------------------------+
|                 SAMPLE STRUCTURE TOOLS               |
+------------------------------------------------------+
| 201) Sample structures from extxyz                   |
| 202) FPS sampling by PyNEP [deprecated]              |
| 203) FPS sampling by NepTrain [preferred]            |
| 204) Perturb structure                               |
| 205) Select max force deviation structs              |
+------------------------------------------------------+
```

## Command Reference

### Uniform/Random Sampling

Select structures using uniform or random sampling.

```bash
python Scripts/sample_structures/sample_structures.py <input.xyz> <method> <num_samples> [skip_initial]

# Examples
python Scripts/sample_structures/sample_structures.py dump.xyz uniform 50
python Scripts/sample_structures/sample_structures.py dump.xyz random 100 500
```

**Parameters:**
- `input.xyz`: Input trajectory file
- `method`: `uniform` or `random`
- `num_samples`: Number of frames to select
- `skip_initial`: (optional) Skip first N frames

**How it works:**
- **Uniform**: Uses `numpy.linspace` to select evenly spaced indices
- **Random**: Uses `numpy.random.choice` to select random indices without replacement

**Output:** `sampled_structures.xyz`

### Farthest Point Sampling (FPS) - NepTrain (Preferred)

Select diverse structures using FPS on NEP descriptors.

```bash
python Scripts/sample_structures/neptrain_select_structs.py <sample.xyz> <train.xyz> <nep.txt>

# Example
python Scripts/sample_structures/neptrain_select_structs.py dump.xyz train.xyz nep.txt
```

**Selection methods (interactive):**
1. **Minimum distance**: Select until max distance to existing training set drops below threshold
2. **Number of structures**: Select between `min_select` and `max_select` structures

**Output:**
- `selected.xyz`
- `select.png` (PCA visualization)
- `pca_sample.txt`, `pca_train.txt`, `pca_selected.txt`

**Citation**: Chen et al., Comput. Phys. Commun., 2025, 317, 109859

### Farthest Point Sampling (FPS) - PyNEP (Deprecated)

```bash
python Scripts/sample_structures/pynep_select_structs.py <sample.xyz> <train.xyz> <nep.txt>
```

**Note**: PyNEP package is no longer actively maintained. Use NepTrain (method 203) instead.

### Structure Perturbation

Generate perturbed structures for training data.

```bash
python Scripts/sample_structures/perturb_structure.py <input.vasp> <pert_num> <cell_pert> <atom_pert> <style>

# Example
python Scripts/sample_structures/perturb_structure.py POSCAR 20 0.03 0.2 uniform
```

**Parameters:**
- `input.vasp`: POSCAR/CONTCAR file
- `pert_num`: Number of perturbed structures (e.g., 20)
- `cell_pert`: Cell perturbation fraction (e.g., 0.03 = 3%)
- `atom_pert`: Atom perturbation distance in Angstrom (e.g., 0.2)
- `style`: `normal`, `uniform`, or `const`

**Output:** `POSCAR_01.vasp`, `POSCAR_02.vasp`, ..., `POSCAR_<N>.vasp`

### Force Deviation Selection

Select structures with high force deviations from active learning.

```bash
python Scripts/sample_structures/select_max_modev.py <top_n> <min_deviation>

# Example
python Scripts/sample_structures/select_max_modev.py 200 0.15
```

**Required files:**
- `active.out` - Max force deviation per structure (from GPUMD `active` command)
- `active.xyz` - Corresponding structures

**Output:** `selected.xyz`

### Frame Range Extraction

Extract a range of frames from a trajectory.

```bash
python Scripts/sample_structures/frame_range.py <input.xyz> <start_fraction> <end_fraction>

# Example: Extract first 80% of frames
python Scripts/sample_structures/frame_range.py dump.xyz 0 0.8
```

**Output:** `dump_0.00_0.80.xyz`

## Common Workflows

### Training Data Preparation

```bash
# 1. Run MD simulation
# ... generate trajectory ...

# 2. Filter by distance
gpumdkit.sh -min_dist_pbc dump.xyz
gpumdkit.sh -filter_dist_pbc dump.xyz 1.0

# 3. Sample diverse structures
gpumdkit.sh  # Select: 2) Sample Structures -> 203

# 4. Add to training set
cat selected.xyz >> train.xyz
```

### Active Learning Data Selection

```bash
# 1. Run MD with current NEP
# ... generate active.out and active.xyz ...

# 2. Select high-deviation structures
python Scripts/sample_structures/select_max_modev.py 100 0.1

# 3. Or use FPS for diversity
python Scripts/sample_structures/neptrain_select_structs.py active.xyz train.xyz nep.txt
```

### Structure Perturbation for Initial Data

```bash
# 1. Start from equilibrium structure
# 2. Generate perturbations
python Scripts/sample_structures/perturb_structure.py POSCAR 50 0.05 0.3 uniform

# 3. Run DFT on perturbed structures
# ... DFT calculations ...

# 4. Add to training set
```

## Dependencies

| Method | Required Packages |
|--------|------------------|
| Uniform/Random | `numpy`, `ase` |
| FPS (NepTrain) | `numpy`, `ase`, `matplotlib`, `scikit-learn`, `scipy`, `NepTrain` |
| FPS (PyNEP) | `numpy`, `ase`, `matplotlib`, `scikit-learn`, `pynep` |
| Perturbation | `dpdata` |
| Force Deviation | `numpy`, `ase` |

Install dependencies:
```bash
pip install numpy ase matplotlib scikit-learn scipy dpdata
pip install neptrain  # For FPS by NepTrain
```

## Important Notes

1. **FPS is preferred**: Farthest Point Sampling provides better diversity than uniform/random
2. **NepTrain over PyNEP**: NepTrain is actively maintained; PyNEP is deprecated
3. **Filter before sampling**: Remove unphysical structures (too close, too large) first
4. **Frame indexing is 0-based**: First frame is index 0
5. **Output is extxyz**: All sampling methods output extxyz format

## See Also

- [Format Conversion](format_conversion.md) - Convert file formats
- [Analyzers](analyzers.md) - Filter and analyze structures
- [NEP Training Guide](nep_training.md) - Complete training workflow
