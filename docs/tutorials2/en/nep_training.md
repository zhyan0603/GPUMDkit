# NEP Training Guide

<div align="center">
  <p>
    <a href="../zh/nep_training.md">中文</a> | <strong>English</strong>
  </p>
</div>

This guide provides a complete workflow for training NEP (Neuroevolution Potential) models using GPUMDkit.

## Overview

NEP training involves:

1. **Data Preparation** - Convert and organize training data
2. **Structure Sampling** - Select diverse structures
3. **Model Training** - Train the NEP model (external)
4. **Validation** - Check training quality
5. **Active Learning** - Iteratively improve the model

## Step 1: Data Preparation

### Convert DFT Data to extxyz

```bash
# From VASP OUTCAR
gpumdkit.sh -out2xyz ./vasp_results/

# From VASP XDATCAR
gpumdkit.sh -xdat2exyz XDATCAR trajectory.xyz

# From CP2K
gpumdkit.sh  # Select: 1) Format Conversion -> 103

# From LAMMPS
gpumdkit.sh -lmp2exyz dump.lammpstrj Li Y Cl
```

### Add Group Labels

Group labels are required for NEP to identify atom types:

```bash
gpumdkit.sh -addgroup POSCAR Li Y Cl
```

This creates `model.xyz` with group labels.

### Combine Training Data

```bash
# Combine multiple extxyz files
cat file1.xyz file2.xyz file3.xyz > train.xyz

# Or use composition analysis to select specific compositions
gpumdkit.sh -analyze_comp all_structures.xyz
```

## Step 2: Data Quality Check

### Check Composition

```bash
gpumdkit.sh -analyze_comp train.xyz
```

### Check Minimum Distances

```bash
# Fast check (no PBC)
gpumdkit.sh -min_dist train.xyz

# Accurate check (with PBC)
gpumdkit.sh -min_dist_pbc train.xyz
```

### Check Property Ranges

```bash
# Check force range
gpumdkit.sh -range train.xyz force

# Check energy range
gpumdkit.sh -range train.xyz energy

# With histogram
gpumdkit.sh -range train.xyz force hist
```

### Filter Structures

```bash
# Filter by minimum distance
gpumdkit.sh -filter_dist_pbc train.xyz 1.0

# Filter by box size
gpumdkit.sh -filter_box train.xyz 20

# Filter by force value
gpumdkit.sh -filter_value train.xyz force 20
```

## Step 3: Structure Sampling

### Sample Diverse Structures

```bash
# FPS sampling (preferred)
gpumdkit.sh  # Select: 2) Sample Structures -> 203

# Or uniform sampling
python Scripts/sample_structures/sample_structures.py train.xyz uniform 100
```

### Perturb Structures

Generate additional training data from existing structures:

```bash
python Scripts/sample_structures/perturb_structure.py POSCAR 20 0.03 0.2 uniform
```

## Step 4: NEP Model Training

### Prepare nep.in

Create `nep.in` with training parameters:

```
type 3 Li Y Cl
cutoff 6 4
n_max 4
l_max 4
neuron 100
batch 1000
generation 100000
```

### Run Training

```bash
# External to GPUMDkit
gpumd
```

## Step 5: Validation

### Plot Training Results

```bash
# Plot training loss and parity plots
gpumdkit.sh -plt train

# Plot test predictions
gpumdkit.sh -plt prediction

# Combined train/test plots
gpumdkit.sh -plt train_test

# Force error analysis
gpumdkit.sh -plt force_errors
```

### Check Convergence

```bash
# Plot loss curves
gpumdkit.sh -plt train
```

Look for:
- Smooth decrease in loss
- No oscillations
- Convergence to low values

### Analyze Descriptors

```bash
# Calculate descriptors
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li

# Visualize with PCA
gpumdkit.sh -plt des pca

# Or UMAP
gpumdkit.sh -plt des umap
```

## Step 6: Active Learning

### Run MD with Current NEP

```bash
# Setup MD simulation
gpumdkit.sh  # Select: 3) Workflow -> 302
```

### Select Diverse Structures

```bash
# Filter by distance
gpumdkit.sh -min_dist_pbc dump.xyz
gpumdkit.sh -filter_box dump.xyz 13

# Sample diverse structures
gpumdkit.sh  # Select: 2) Sample Structures -> 203
```

### Compute DFT Reference

```bash
# Prepare SCF calculations
gpumdkit.sh  # Select: 3) Workflow -> 301

# Run DFT (external)
```

### Retrain Model

```bash
# Add new data to training set
cat new_structures.xyz >> train.xyz

# Retrain NEP (external)
```

## Complete Workflow Example

```bash
# 1. Convert DFT data
gpumdkit.sh -out2xyz ./dft_results/

# 2. Add group labels
gpumdkit.sh -addgroup POSCAR Li Y Cl

# 3. Combine training data
cat */model.xyz > train.xyz

# 4. Check data quality
gpumdkit.sh -analyze_comp train.xyz
gpumdkit.sh -min_dist_pbc train.xyz
gpumdkit.sh -range train.xyz force

# 5. Filter structures
gpumdkit.sh -filter_dist_pbc train.xyz 1.0

# 6. Sample diverse structures
gpumdkit.sh  # Select: 2) Sample Structures -> 203

# 7. Train NEP (external)
# ... create nep.in and run gpumd ...

# 8. Validate training
gpumdkit.sh -plt train
gpumdkit.sh -plt prediction
gpumdkit.sh -plt force_errors

# 9. Analyze descriptors
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li
gpumdkit.sh -plt des pca
```

## Tips for Better Training

1. **Diverse Data**: Include various structures, temperatures, and compositions
2. **Balanced Dataset**: Ensure all element types are well-represented
3. **Quality Filter**: Remove structures with too-close atoms or high forces
4. **Sufficient Data**: At least 1000-5000 structures for good accuracy
5. **Validation Set**: Keep 10-20% of data for testing

## Troubleshooting

### Issue: High RMSE

**Solution**: 
- Add more diverse training data
- Increase `neuron` and `n_max` parameters
- Check for outliers in training data

### Issue: Poor Predictions

**Solution**:
- Verify training data quality
- Add structures from different phases/temperatures
- Increase training `generation`

### Issue: Overfitting

**Solution**:
- Use more training data
- Reduce model complexity
- Use regularization

## See Also

- [Format Conversion](format_conversion.md) - Convert file formats
- [Calculators](calculators.md) - Compute properties
- [Visualization](visualization.md) - Plot results
- [Structure Sampling](sampling.md) - Sample structures
