<div align="center">
  <h1>NEP Training Guide</h1>
  <p>
    <strong>English</strong> | <a href="../zh/nep_training.md">简体中文</a>
  </p>
</div>

Complete workflow for training NEP (Neuroevolution Potential) models.

## Overview

1. **Data Preparation** - Convert and organize training data
2. **Quality Check** - Validate and filter structures
3. **Sampling** - Select diverse structures
4. **Training** - Train the NEP model (external)
5. **Validation** - Check training quality
6. **Active Learning** - Iteratively improve the model

## Step 1: Data Preparation

```bash
# Convert DFT data
gpumdkit.sh -out2xyz ./vasp_results/

# Add group labels
gpumdkit.sh -addgroup POSCAR Li Y Cl

# Combine training data
cat */model.xyz > train.xyz
```

## Step 2: Quality Check

```bash
# Check composition
gpumdkit.sh -analyze_comp train.xyz

# Check distances
gpumdkit.sh -min_dist_pbc train.xyz

# Check force range
gpumdkit.sh -range train.xyz force

# Filter structures
gpumdkit.sh -filter_dist_pbc train.xyz 1.0
```

## Step 3: Sampling

```bash
# FPS sampling (recommended)
gpumdkit.sh  # Select: 2) Sample Structures -> 203

# Or uniform sampling
python Scripts/sample_structures/sample_structures.py train.xyz uniform 100
```

## Step 4: Training

Create `nep.in`:

```
type 3 Li Y Cl
cutoff 6 4
n_max 4
l_max 4
neuron 100
batch 1000
generation 100000
```

Run training:

```bash
gpumd
```

## Step 5: Validation

```bash
# Plot training results
gpumdkit.sh -plt train
gpumdkit.sh -plt prediction
gpumdkit.sh -plt force_errors

# Analyze descriptors
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li
gpumdkit.sh -plt des pca
```

## Step 6: Active Learning

```bash
# Run MD with current NEP
gpumdkit.sh  # Select: 3) Workflow -> 302

# Filter and sample
gpumdkit.sh -min_dist_pbc dump.xyz
gpumdkit.sh -filter_box dump.xyz 13

# Prepare DFT calculations
gpumdkit.sh  # Select: 3) Workflow -> 301

# Retrain with new data
cat new_structures.xyz >> train.xyz
```

## Tips

- Include diverse structures, temperatures, and compositions
- At least 1000-5000 structures for good accuracy
- Keep 10-20% of data for testing
- Check for outliers before training
