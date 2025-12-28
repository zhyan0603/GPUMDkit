# Function 2 - Sample Structures

This section covers the structure sampling tools in GPUMDkit (Interactive Mode - Function 2), which provide various methods for sampling and selecting structures from datasets. These tools are essential for creating diverse NEP training sets, implementing active learning workflows, and generating augmented datasets.

## Overview

Structure sampling is crucial for creating high-quality training datasets. GPUMDkit provides multiple sampling strategies:

- **Statistical sampling** - Uniform and random selection
- **Descriptor-based sampling** - FPS (Farthest Point Sampling) in descriptor space
- **Uncertainty-based sampling** - Select structures with high prediction uncertainty
- **Perturbation-based sampling** - Generate variations of existing structures

**Key concept:** Good sampling maximizes the diversity of atomic configurations while maintaining manageable dataset sizes.

## Interactive Mode Access

```bash
gpumdkit.sh
# Select: 2) Sample Structures
```

You'll see the following menu:

```
 ------------>>
 201) Sample structures from extxyz
 202) Sample structures by pynep
 203) Sample structures by neptrain
 204) Perturb structure
 205) Select max force deviation structs
 000) Return to the main menu
 ------------>>
```

## Command-Line Usage

For PyNEP sampling (the most commonly used method):

```bash
# PyNEP farthest point sampling
gpumdkit.sh -pynep <candidates.xyz> <train.xyz> <nep.txt>
```

Other sampling methods are accessed via interactive mode.

---

## Sampling Methods

### Uniform/Random Sampling (`sample_structures.py`)

**Option 201** in interactive mode.

Sample structures using statistical methods without requiring a NEP model.

**Algorithms:**

1. **Uniform Sampling**
   - Selects structures at evenly spaced intervals
   - Formula: `indices = range(0, total, step)` where `step = total / num_samples`
   - Ensures representation across entire dataset
   - Preserves temporal ordering (useful for trajectories)

2. **Random Sampling**
   - Randomly selects structures without replacement
   - Each structure has equal probability of selection
   - No bias toward any region of dataset
   - Better for unordered datasets

**Interactive Mode:**

Select option `201` from the sample structures menu:

```bash
201
```

You will see the following prompt:

```
>-------------------------------------------------<
| This function calls the script in Scripts       |
| Script: sample_structures.py                    |
| Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
>-------------------------------------------------<
Input <extxyz_file> <sampling_method> <num_samples>
Sampling_method: 'uniform' or 'random'
Examp: train.xyz uniform 50
------------>>
```

Enter the filename, method, and number of samples:

```bash
train.xyz uniform 50
```

**Direct execution:**
```bash
python ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py train.xyz uniform 100
python ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py train.xyz random 100
```

**Output:** `sampled_structures.xyz`

**When to use each method:**

| Method | Best For | Pros | Cons |
|--------|----------|------|------|
| **Uniform** | MD trajectories, ordered datasets | Preserves temporal information, evenly spaced | May miss outliers |
| **Random** | Unordered datasets, general reduction | Unbiased, simple | No guarantee of coverage |

**Example workflow:**
```bash
# Quick dataset reduction (10,000 → 1,000 structures)
python ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py \
    large_dataset.xyz uniform 1000

# Verify sampling quality
gpumdkit.sh -range sampled_structures.xyz force
gpumdkit.sh -analyze_comp sampled_structures.xyz
```

---

### PyNEP Farthest Point Sampling (`pynep_select_structs.py`)

**Option 202** in interactive mode.

Select maximally diverse structures using Farthest Point Sampling (FPS) in NEP descriptor space. This is the **recommended method** for creating high-quality training sets.

**Algorithm: Farthest Point Sampling (FPS)**

1. **Descriptor Calculation**
   ```
   For each structure:
       Calculate NEP descriptors for all atoms
       Aggregate to structure-level descriptor
   ```

2. **FPS Procedure**
   ```python
   1. Initialize: Select random first structure
   2. For i = 2 to n_select:
       a. Calculate distance from each remaining structure
          to all already-selected structures
       b. Find structure with maximum minimum distance
       c. Add to selected set
   3. Return selected structures
   ```

**Why FPS Works:**
- Maximizes diversity in descriptor space
- Ensures even coverage of configuration space
- No redundant similar structures
- Works in high-dimensional space (unlike physical space)

**Command-line:**
```bash
gpumdkit.sh -pynep candidates.xyz selected.xyz nep.txt
```

**Interactive:** Select option `202`

**Direct execution:**
```bash
python ${GPUMDkit_path}/Scripts/sample_structures/pynep_select_structs.py \
    candidates.xyz train.xyz nep.txt
```

**Parameters:**
- `candidates.xyz` - Pool of structures to select from
- `train.xyz` - Output file with selected structures  
- `nep.txt` - NEP model for descriptor calculation

**Interactive prompts:**
```
Input the number of structures to select:
>> 100

Select sampling method:
1) FPS (Farthest Point Sampling) - Recommended
2) Random sampling
>> 1
```

**Requirements:** PyNEP package installed
```bash
pip install pynep
```

**Use cases:**
- Active learning structure selection
- Training set creation from large pools
- Reducing redundant structures
- Ensuring diverse coverage

**Example:**
```bash
# After MD sampling, select 100 most diverse structures
gpumdkit.sh -pynep md_trajectory.xyz selected_100.xyz nep.txt

# Verify diversity
gpumdkit.sh -calc des umap selected_100.xyz desc.npy nep.txt Li
gpumdkit.sh -plt des umap desc.npy
# Should show good spread in descriptor space
```

**Performance:**
- Fast for < 10,000 candidates
- Moderate for 10,000-50,000 candidates  
- Slow for > 50,000 candidates (consider pre-filtering)

**Tips:**
- Pre-filter candidates with quality checks before FPS
- Use higher n_select for diverse systems
- Visualize descriptor space to verify sampling quality

---

### NEPtrain Selection (`neptrain_select_structs.py`)

Alternative FPS implementation using NEPtrain.

**Interactive:** Select option `203`

**Direct execution:**
```bash
python ${GPUMDkit_path}/Scripts/sample_structures/neptrain_select_structs.py candidates.xyz train.xyz nep.txt
```

**Requires:** `neptrain` package installed

### Structure Perturbation (`perturb_structure.py`)

**Option 204** in interactive mode.

Generate perturbed variants of structures for training data augmentation. This creates slightly modified copies to improve model robustness.

**Algorithm: Structure Perturbation**

1. **Cell Perturbation**
   ```python
   # For each lattice vector
   for i in range(3):
       perturbation = random_value × cell_fraction
       lattice_vector[i] += perturbation
   ```

2. **Atomic Perturbation**
   ```python
   # For each atom
   for atom in structure:
       displacement = random_vector × atom_distance
       atom.position += displacement
   ```

**Perturbation Styles:**

| Style | Distribution | When to Use |
|-------|--------------|-------------|
| `uniform` | Uniform random in [-δ, +δ] | General purpose |
| `normal` | Gaussian (mean=0, std=δ) | Physical thermal motion |
| `const` | Fixed displacement δ | Systematic exploration |

**Interactive Mode:**

Select option `204`:

```bash
204
```

You will see:
```
>-------------------------------------------------<
| This function calls the script in Scripts       |
| Script: perturb_structure.py                    |
| Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
>-------------------------------------------------<
Input <structure.vasp> <n_structures> <cell_pert_fraction> <atom_pert_distance> <style>
Examp: POSCAR 20 0.03 0.2 normal
------------>>
```

**Direct execution:**
```bash
python ${GPUMDkit_path}/Scripts/sample_structures/perturb_structure.py \
    structure.vasp 20 0.03 0.2 normal
```

**Parameters:**
- `structure.vasp` - Input structure (VASP POSCAR format)
- `n_structures` - Number of perturbed variants to generate
- `cell_pert_fraction` - Cell perturbation (fraction, e.g., 0.03 = 3%)
- `atom_pert_distance` - Atomic displacement (Angstroms)
- `style` - Perturbation type: `uniform`, `normal`, or `const`

**Output:** Multiple perturbed VASP files: `perturbed_001.vasp`, `perturbed_002.vasp`, ...

**Recommended Parameters by System Type:**

#### Rigid Crystals (e.g., diamond, Si)
```bash
cell_fraction: 0.01-0.02
atom_distance: 0.05-0.15 Å
style: normal
```

#### Soft Materials (e.g., molecular crystals, polymers)
```bash
cell_fraction: 0.03-0.05
atom_distance: 0.2-0.3 Å
style: normal
```

#### Ionic Conductors (e.g., LGPS, LLZO)
```bash
cell_fraction: 0.02-0.04
atom_distance: 0.15-0.25 Å
style: uniform
```

#### Phase Transitions
```bash
cell_fraction: 0.05-0.10
atom_distance: 0.3-0.5 Å
style: uniform
```

**Complete Workflow:**

```bash
# 1. Generate perturbed structures
python ${GPUMDkit_path}/Scripts/sample_structures/perturb_structure.py \
    POSCAR 20 0.03 0.2 normal

# 2. Set up DFT calculations
gpumdkit.sh
# Select: 3) Workflow → 301) SCF batch pretreatment

# 3. Run VASP calculations
for dir in struct_fp_*/; do
    cd $dir && sbatch submit.sh && cd ..
done

# 4. Collect results
gpumdkit.sh -out2xyz ./struct_fp_*/

# 5. Add to training set
cat train.xyz perturbed_results.xyz > train_augmented.xyz

# 6. Verify improvement
nep
gpumdkit.sh -plt train
```

**Use Cases:**
- **Data augmentation** - Increase training set size
- **Rare events** - Sample transition states
- **Thermal effects** - Simulate temperature effects
- **Robustness** - Improve model generalization

**Tips:**
- Start with small perturbations and increase if needed
- Always run DFT on perturbed structures
- Check that perturbed structures are physically reasonable
- Use `normal` style for realistic thermal motion
- Combine with FPS to avoid redundancy

---

### Max Force Deviation Selection (`select_max_modev.py`)

**Option 205** in interactive mode.

Select structures with highest force deviations from active learning simulations. This identifies configurations where the NEP model is most uncertain.

**Algorithm: Maximum Force Deviation**

1. **Calculate Force Deviation**
   ```
   For each structure in active.xyz:
       max_deviation = max|F_NEP - F_committee|
       where F_committee is average of multiple predictions
   ```

2. **Selection**
   ```python
   1. Read all structures with force deviations
   2. Sort by max_deviation (descending)
   3. Filter: keep structures with deviation > threshold
   4. Select top N structures
   5. Return selected structures
   ```

**What is Force Deviation?**

When GPUMD runs in active learning mode with `active` keyword, it:
1. Makes NEP predictions with slightly perturbed parameters
2. Calculates variance in force predictions
3. Identifies structures with high variance (uncertainty)
4. Saves these in `active.xyz` with deviation values

**Interactive Mode:**

Select option `205`:

```bash
205
```

You will see:
```
>-------------------------------------------------<
| This function calls the script in Scripts       |
| Script: select_max_modev.py                     |
| Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
>-------------------------------------------------<
Input <number_to_select> <force_deviation_threshold>
Examp: 100 0.15
------------>>
```

**Direct execution:**
```bash
python ${GPUMDkit_path}/Scripts/sample_structures/select_max_modev.py 100 0.15
```

**Parameters:**
- `number_to_select` - Maximum number of structures to select
- `force_deviation_threshold` - Minimum deviation (eV/Å)

**Input:** `active.xyz` from GPUMD active learning mode

**Output:** 
- `selected.xyz` - Structures with high uncertainty
- `report.txt` - Statistics on selected structures

**Setting up Active Learning in GPUMD:**

```bash
# In run.in
ensemble npt_ber 1200 1200 100 0 0 0 1000
dump_thermo 1000
dump_position 1000

# Enable active learning
# active <interval> <threshold>
active 100 0.15    # Check every 100 steps, threshold 0.15 eV/Å

run 100000
```

**Parameters:**
- `interval` - How often to check (steps)
- `threshold` - Force deviation threshold (eV/Å)

**Complete Active Learning Workflow:**

```bash
# Step 1: Run MD with active mode
cat > run.in << EOF
velocity        1200
ensemble        npt_ber 1200 1200 100 0 0 0 1000
dump_thermo     1000
dump_position   1000
active          100 0.15
run             100000
EOF

gpumd
# Generates active.xyz with uncertain structures

# Step 2: Select most uncertain structures
python ${GPUMDkit_path}/Scripts/sample_structures/select_max_modev.py 50 0.15
# Selects top 50 structures with deviation > 0.15 eV/Å

# Step 3: Run DFT on selected structures
gpumdkit.sh
# Select: 3) Workflow → 301) SCF batch pretreatment

# Step 4: Add to training set and retrain
cat train.xyz dft_results.xyz > train_new.xyz
nep
```

**Threshold Guidelines:**

| Threshold | Selectivity | When to Use |
|-----------|-------------|-------------|
| 0.10 eV/Å | Very selective | Late iterations, high-accuracy models |
| 0.15 eV/Å | Moderate | Standard active learning |
| 0.20 eV/Å | Permissive | Early iterations, exploration |
| 0.30 eV/Å | Very permissive | Initial sampling, diverse systems |

**Interpreting Results:**

```bash
# Check deviation distribution
grep "max_force_dev" selected.xyz | awk '{print $NF}' > deviations.dat

# Plot histogram (requires matplotlib)
python << EOF
import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt('deviations.dat')
plt.hist(data, bins=20)
plt.xlabel('Force Deviation (eV/Å)')
plt.ylabel('Count')
plt.title('Selected Structure Deviations')
plt.savefig('deviations.png')
EOF
```

**Advantages:**
- ✅ Targets model uncertainty directly
- ✅ Efficient - only DFT needed structures
- ✅ No manual selection required
- ✅ Works well for refinement

**Limitations:**
- ❌ Requires active mode in GPUMD
- ❌ May miss systematic errors
- ❌ Needs reasonable initial model

**Comparison with FPS:**

| Method | Best For | Requires |
|--------|----------|----------|
| **Max Force Dev** | Model refinement, targeted improvement | Active mode run, existing NEP |
| **FPS** | Diverse sampling, coverage | Large candidate pool, NEP model |

**Pro tip:** Combine both methods!
```bash
# Step 1: Active learning to find uncertain structures
# (generates active.xyz)

# Step 2: Select by force deviation
python select_max_modev.py 200 0.15
# Gets 200 uncertain structures

# Step 3: Further diversify with FPS
gpumdkit.sh -pynep selected.xyz final_selection.xyz nep.txt
# Selects 100 most diverse among uncertain ones

# Step 4: Run DFT only on final 100 structures
# Maximum efficiency!
```

---

## Interactive Mode

Access sampling tools through the interactive menu:

```bash
gpumdkit.sh
# Select: 2) Sample Structures
# Choose option 201-205
```

### Menu Options

```
 ------------>>
 201) Sample structures from extxyz
 202) Sample structures by pynep
 203) Sample structures by neptrain
 204) Perturb structure
 205) Select max force deviation structs
 000) Return to the main menu
 ------------>>
```

## Common Workflows

### Create Initial Training Set

```bash
# 1. Start with large dataset
# 2. Random sample to manageable size
python ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py large_dataset.xyz random 2000

# 3. Apply FPS for diversity
gpumdkit.sh -pynep sampled_structures.xyz train.xyz initial_nep.txt
```

### Active Learning Cycle

```bash
# 1. Run GPUMD with active learning mode
#    Add to run.in: active 100 0.15

# 2. After simulation, select uncertain structures
python ${GPUMDkit_path}/Scripts/sample_structures/select_max_modev.py 50 0.15

# 3. Run DFT on selected structures

# 4. Add to training set and retrain NEP
```

### Data Augmentation

```bash
# 1. Have equilibrium structures
# 2. Generate perturbed variants
python ${GPUMDkit_path}/Scripts/sample_structures/perturb_structure.py POSCAR 20 0.03 0.2 normal

# 3. Run DFT on perturbed structures
# 4. Add to training set
```

## Sampling Strategy Guide

Choosing the right sampling method depends on your goal and dataset characteristics.

### Decision Tree

```
┌─ Do you have a NEP model?
│
├─ NO → Use statistical sampling
│  │
│  ├─ Ordered dataset (trajectory)? → Uniform sampling
│  └─ Unordered dataset? → Random sampling
│
└─ YES → Do you want diversity or uncertainty?
   │
   ├─ DIVERSITY → Use FPS (PyNEP or NEPtrain)
   │  ├─ Large candidate pool (>10k) → Pre-filter, then FPS
   │  └─ Small candidate pool (<10k) → Direct FPS
   │
   └─ UNCERTAINTY → Use active learning
      ├─ Have active.xyz? → Max force deviation selection
      └─ No active.xyz? → Run GPUMD with active mode first
```

### Recommended Workflows by Task

#### Task 1: Initial Training Set Creation

**Scenario:** Starting NEP development with raw DFT data

```bash
# Step 1: Start with all available data
cat dft_calculations/*/*.xyz > all_data.xyz
# Total: 10,000 structures

# Step 2: Quick random sample to manageable size
python ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py \
    all_data.xyz random 2000
# Reduced to: 2,000 structures

# Step 3: Train initial NEP
mv sampled_structures.xyz train.xyz
nep
# Initial model ready

# Step 4: Use FPS for final diverse set
gpumdkit.sh -pynep train.xyz train_diverse_500.xyz nep.txt
# Final: 500 most diverse structures

# Step 5: Retrain with diverse set
mv train_diverse_500.xyz train.xyz
nep
```

---

#### Task 2: Active Learning Iteration

**Scenario:** Improving existing NEP model

```bash
# Step 1: Run MD with active mode
# (in run.in: active 100 0.15)
gpumd
# Generates: active.xyz with ~500 uncertain structures

# Step 2: Select top uncertain structures
python ${GPUMDkit_path}/Scripts/sample_structures/select_max_modev.py 100 0.15
# Selected: 100 structures with high force deviation

# Step 3: Run DFT calculations
gpumdkit.sh
# Workflow → SCF batch pretreatment

# Step 4: Add to training set
cat train.xyz dft_results.xyz > train_expanded.xyz

# Step 5: Retrain NEP
mv train_expanded.xyz train.xyz
nep
```

---

#### Task 3: Data Augmentation

**Scenario:** Limited DFT data, need more training structures

```bash
# Step 1: Identify key structures
gpumdkit.sh -analyze_comp train.xyz
# Find underrepresented compositions

# Step 2: Perturb these structures
python ${GPUMDkit_path}/Scripts/sample_structures/perturb_structure.py \
    important_structure.vasp 50 0.03 0.2 normal
# Generated: 50 perturbed variants

# Step 3: Run DFT on perturbed structures
# (Use workflow SCF batch pretreatment)

# Step 4: Add to training
cat train.xyz perturbed_dft_results.xyz > train_augmented.xyz

# Step 5: Verify no redundancy
gpumdkit.sh -pynep train_augmented.xyz train_final.xyz nep.txt
# Ensures diversity

# Step 6: Train
mv train_final.xyz train.xyz
nep
```

---

#### Task 4: Multi-Composition Dataset

**Scenario:** Training NEP for multiple compositions

```bash
# Step 1: Analyze compositions
gpumdkit.sh -analyze_comp all_structures.xyz
# Results:
#   Li8YCl9: 1000 structures
#   Li10YCl11: 500 structures
#   Li6YCl7: 200 structures  ← Underrepresented

# Step 2: Balance compositions
# Export each composition separately
# (Use interactive mode to export)

# Step 3: Sample proportionally
python sample_structures.py Li8YCl9.xyz random 300
python sample_structures.py Li10YCl11.xyz random 300  
python sample_structures.py Li6YCl7.xyz random 200  # Keep all

# Step 4: Combine balanced set
cat Li8YCl9_sampled.xyz Li10YCl11_sampled.xyz Li6YCl7.xyz > train_balanced.xyz

# Step 5: Apply FPS for final diversity
gpumdkit.sh -pynep train_balanced.xyz train_final.xyz nep.txt
```

---

## Advanced Techniques

### Hierarchical Sampling

Combine multiple methods for optimal results:

```bash
# Level 1: Coarse random sampling (100k → 10k)
python sample_structures.py huge_dataset.xyz random 10000

# Level 2: Filter by quality (10k → 5k)
gpumdkit.sh -filter_value sampled_structures.xyz force 30
gpumdkit.sh -filter_dist filtered_force.xyz 1.5

# Level 3: FPS for diversity (5k → 1k)
gpumdkit.sh -pynep filtered_dist.xyz train_hierarchical.xyz nep.txt

# Result: 1k highly diverse, high-quality structures
```

---

### Adaptive Sampling

Adjust sampling based on model performance:

```python
#!/bin/bash
# adaptive_sampling.sh

iteration=1
target_rmse=2.0  # meV/atom

while true; do
    echo "=== Iteration $iteration ==="
    
    # Train NEP
    nep
    
    # Check RMSE
    current_rmse=$(grep "Energy RMSE" loss.out | tail -1 | awk '{print $4}')
    echo "Current RMSE: $current_rmse meV/atom"
    
    # Check convergence
    if (( $(echo "$current_rmse < $target_rmse" | bc -l) )); then
        echo "Target accuracy reached!"
        break
    fi
    
    # Adaptive sampling based on RMSE
    if (( $(echo "$current_rmse > 5.0" | bc -l) )); then
        # High RMSE: aggressive sampling
        n_sample=100
        temp=1500
    elif (( $(echo "$current_rmse > 3.0" | bc -l) )); then
        # Medium RMSE: moderate sampling
        n_sample=50
        temp=1200
    else
        # Low RMSE: conservative sampling
        n_sample=30
        temp=1000
    fi
    
    echo "Sampling $n_sample structures at $temp K"
    
    # Run active learning with adaptive parameters
    # ... (run GPUMD, select, DFT, etc.)
    
    ((iteration++))
done
```

---

### Cross-Validation Sampling

Ensure training set quality:

```bash
# Split data for cross-validation
total=$(grep -c "^[0-9]" train.xyz)
train_size=$((total * 80 / 100))
test_size=$((total - train_size))

# Create 5-fold splits
for fold in {1..5}; do
    # Random shuffle and split
    python << EOF
import random
from ase.io import read, write

structures = read('train.xyz', ':')
random.shuffle(structures)

test_start = $((fold-1)) * $test_size
test_end = $fold * $test_size

test = structures[test_start:test_end]
train = structures[:test_start] + structures[test_end:]

write(f'train_fold{$fold}.xyz', train)
write(f'test_fold{$fold}.xyz', test)
EOF

    # Train on fold
    cp train_fold${fold}.xyz train.xyz
    cp test_fold${fold}.xyz test.xyz
    nep
    
    # Record RMSE
    grep "RMSE" loss.out >> cv_results.txt
done

# Analyze cross-validation results
python << EOF
import numpy as np
data = np.loadtxt('cv_results.txt', usecols=[3])
print(f"Mean RMSE: {np.mean(data):.2f} ± {np.std(data):.2f} meV/atom")
EOF
```

---

## Best Practices

### 1. Always Visualize Sampling Results

```bash
# After any sampling, check descriptor space coverage
gpumdkit.sh -calc des umap sampled.xyz desc.npy nep.txt Li
gpumdkit.sh -plt des umap desc.npy

# Good sampling shows:
# ✓ Even distribution across descriptor space
# ✓ No large gaps or clusters
# ✓ Coverage of extremes
```

---

### 2. Quality Before Quantity

```bash
# Don't just add more structures - add diverse ones
# Bad approach:
cat train.xyz more_structures.xyz > train_big.xyz  # May add redundancy

# Good approach:
cat train.xyz more_structures.xyz > combined.xyz
gpumdkit.sh -pynep combined.xyz train_diverse.xyz nep.txt  # Ensure diversity
```

---

### 3. Document Sampling Decisions

Keep a sampling log:

```bash
# sampling_log.txt
2025-01-15: Initial random sample 10000 → 2000 structures
2025-01-16: FPS selection 2000 → 500 structures  
2025-01-17: Active learning added 50 structures (threshold 0.15 eV/Å)
2025-01-18: Perturbation generated 20 variants (cell: 0.03, atom: 0.2)
```

---

### 4. Balance Composition

```bash
# For multi-element systems, ensure balanced representation
gpumdkit.sh -analyze_comp train.xyz

# If unbalanced:
# - Sample more from underrepresented compositions
# - Or use weighted training (gpumdkit.sh -addweight)
```

---

### 5. Validate Sampling Quality

```bash
# Check that sampled set maintains property distributions
gpumdkit.sh -range original.xyz force
gpumdkit.sh -range sampled.xyz force
# Ranges should be similar

gpumdkit.sh -min_dist_pbc original.xyz
gpumdkit.sh -min_dist_pbc sampled.xyz
# Min distances should be similar
```

---

### 6. Iterative Refinement

Don't try to get perfect sampling in one step:

```
Iteration 1: Coarse sampling (10k → 1k)
Iteration 2: Train initial NEP
Iteration 3: Active learning adds 100
Iteration 4: Train improved NEP
Iteration 5: Active learning adds 50
Iteration 6: Train final NEP ← Converged!
```

---

### 7. Start Conservative, Then Explore

```bash
# Early iterations: strict filtering
min_dist=1.5
max_force=25
n_select=30

# Middle iterations: balanced
min_dist=1.4
max_force=30
n_select=50

# Late iterations: targeted
min_dist=1.3
max_force=35
n_select=100  # More structures for refinement
```

---

## Troubleshooting

### Problem: FPS Takes Too Long

**Symptoms:**
```bash
gpumdkit.sh -pynep large_pool.xyz selected.xyz nep.txt
# Hangs for hours...
```

**Solutions:**

1. **Pre-filter before FPS**
   ```bash
   # Reduce pool size first
   python sample_structures.py large_pool.xyz random 5000
   gpumdkit.sh -pynep sampled_structures.xyz selected.xyz nep.txt
   ```

2. **Use hierarchical sampling**
   ```bash
   # Level 1: Random 50k → 5k
   # Level 2: FPS 5k → 500
   ```

3. **Upgrade PyNEP**
   ```bash
   pip install --upgrade pynep
   # Newer versions are faster
   ```

---

### Problem: Perturbed Structures Are Unphysical

**Symptoms:**
```bash
# After perturbation, minimum distances too small
gpumdkit.sh -min_dist_pbc perturbed_001.vasp
# Output: 0.8 Å ← Too small!
```

**Solutions:**

1. **Reduce perturbation**
   ```bash
   # Was: 0.05 cell, 0.4 Å atom
   # Try: 0.02 cell, 0.15 Å atom
   ```

2. **Filter after perturbation**
   ```bash
   # Convert to extxyz
   for f in perturbed_*.vasp; do
       gpumdkit.sh -pos2exyz $f ${f%.vasp}.xyz
   done
   cat perturbed_*.xyz > all_perturbed.xyz
   
   # Filter
   gpumdkit.sh -filter_dist all_perturbed.xyz 1.5
   ```

3. **Use normal distribution**
   ```bash
   # More realistic than uniform
   style=normal  # Instead of uniform
   ```

---

### Problem: Active Mode Selects Too Few/Many Structures

**Too few:**
```bash
# Only 5 structures in active.xyz
# Solution: Lower threshold or increase interval
active 50 0.10    # Was: active 100 0.15
```

**Too many:**
```bash
# 500+ structures in active.xyz (too expensive for DFT)
# Solution: Raise threshold or be more selective
active 100 0.20   # Was: active 100 0.15

# Or post-process with max force dev selection
python select_max_modev.py 50 0.20
```

---

### Problem: Sampled Set Has Poor Test Performance

**Symptoms:**
```bash
# Train RMSE: 2.0 meV/atom
# Test RMSE: 8.0 meV/atom ← Large gap!
```

**Causes & Solutions:**

1. **Biased sampling**
   ```bash
   # Check if test set different from training
   gpumdkit.sh -analyze_comp train.xyz
   gpumdkit.sh -analyze_comp test.xyz
   # Ensure similar compositions
   ```

2. **Overfitting to sampled structures**
   ```bash
   # Solution: More diverse sampling
   # Increase FPS selection number
   # Add perturbation for robustness
   ```

3. **Test set has outliers**
   ```bash
   # Check test set quality
   gpumdkit.sh -range test.xyz force
   gpumdkit.sh -plt force_errors
   ```

---

## Performance Optimization

### Memory-Efficient Sampling

For very large datasets (>100k structures):

```python
#!/usr/bin/env python3
# memory_efficient_sampling.py

from ase.io import read, write
import random

def chunked_sampling(input_file, output_file, n_select, chunk_size=1000):
    """Sample from large file without loading all into memory"""
    
    # First pass: count structures
    n_total = 0
    with open(input_file) as f:
        for line in f:
            if line.strip().isdigit():
                n_total += 1
    
    # Calculate sampling probability
    prob = n_select / n_total
    
    # Second pass: random sampling
    selected = []
    for atoms in read(input_file, ':'):
        if random.random() < prob:
            selected.append(atoms)
            if len(selected) >= n_select:
                break
    
    write(output_file, selected)
    print(f"Selected {len(selected)} from {n_total} structures")

if __name__ == "__main__":
    chunked_sampling('huge_dataset.xyz', 'sampled.xyz', 2000)
```

---

### Parallel FPS

For faster FPS with large datasets:

```bash
# Split dataset
python << EOF
from ase.io import read, write
import numpy as np

structures = read('large_pool.xyz', ':')
n_chunks = 4
chunk_size = len(structures) // n_chunks

for i in range(n_chunks):
    start = i * chunk_size
    end = (i + 1) * chunk_size if i < n_chunks - 1 else len(structures)
    write(f'chunk_{i}.xyz', structures[start:end])
EOF

# FPS on each chunk in parallel
for i in {0..3}; do
    gpumdkit.sh -pynep chunk_${i}.xyz selected_${i}.xyz nep.txt &
done
wait

# Combine and final FPS
cat selected_*.xyz > combined_selected.xyz
gpumdkit.sh -pynep combined_selected.xyz final_selected.xyz nep.txt
```

---

## Summary Table

| Sampling Method | Speed | Quality | Use Case | Best For |
|-----------------|-------|---------|----------|----------|
| **Uniform** | ⚡⚡⚡ | ⭐⭐ | Quick reduction | Trajectories |
| **Random** | ⚡⚡⚡ | ⭐⭐ | Unbiased sampling | General datasets |
| **PyNEP FPS** | ⚡⚡ | ⭐⭐⭐⭐ | Diverse selection | Training sets |
| **NEPtrain FPS** | ⚡⚡ | ⭐⭐⭐⭐ | Alternative FPS | Training sets |
| **Perturbation** | ⚡⚡ | ⭐⭐⭐ | Augmentation | Limited data |
| **Max Force Dev** | ⚡⚡⚡ | ⭐⭐⭐⭐⭐ | Uncertainty | Active learning |

**Legend:**
- Speed: ⚡ (slow) to ⚡⚡⚡ (fast)
- Quality: ⭐ (basic) to ⭐⭐⭐⭐⭐ (excellent)

---

---

For more details, see [Scripts/sample_structures/README.md](../../Scripts/sample_structures/README.md)
