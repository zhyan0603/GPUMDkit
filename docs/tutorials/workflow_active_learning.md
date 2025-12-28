# Active Learning Workflow

**Script Location:** `Scripts/workflow/workflow_active_learning_dev.sh`

This tutorial explains the automated active learning workflow for NEP model development.

## Overview

Active learning iteratively improves NEP models by:
1. Generating candidate structures via MD
2. Identifying uncertain configurations
3. Running DFT on selected structures
4. Adding to training set and retraining

This workflow automates the entire cycle.

## Prerequisites

- **GPUMD** installed and accessible
- **VASP** (or other DFT code) set up
- **GPUMDkit** installed
- **PyNEP** installed (for structure selection)
- Initial NEP model trained

## Workflow Steps

### 1. Preparation

Create a working directory with:
- `nep.txt` - Current NEP model
- `train.xyz` - Current training data
- `run.in` - GPUMD configuration
- VASP inputs (`POTCAR`, `INCAR`, `KPOINTS` in `fp/` directory)

### 2. Configuration

Edit `workflow_active_learning_dev.sh` to set:

```bash
# Basic settings
prefix_name=LiF_iter01        # Iteration identifier
min_dist=1.4                  # Minimum atomic distance filter
box_limit=13                  # Box size filter
max_force=30                  # Maximum force filter (eV/Å)

# Sampling
sample_number=100             # Structures to sample from MD
select_number=50              # Structures to select for DFT

# MD settings
md_steps=50000                # MD steps per structure
md_temp=1200                  # Temperature (K)

# DFT settings
dft_partition=intel-sc3       # Cluster partition
dft_queue=huge                # Queue name
dft_nodes=1                   # Nodes per calculation
```

### 3. Run the Workflow

**On SLURM cluster:**
```bash
sbatch workflow_active_learning_dev.sh
```

**Without SLURM:**
```bash
nohup bash workflow_active_learning_dev.sh &>workflow.log &
```

### 4. Monitor Progress

Check log file:
```bash
tail -f workflow.log
```

The workflow will:
1. Run MD sampling with current NEP
2. Filter structures by distance and force
3. Select diverse structures using PyNEP FPS
4. Set up DFT calculations
5. Submit DFT jobs
6. Collect results after DFT completes
7. Add to training set
8. Retrain NEP model

## Workflow Stages

### Stage 1: MD Sampling

Generates candidate structures by running MD at specified temperature.

### Stage 2: Structure Filtering

Applies quality filters:
- Minimum atomic distance > threshold
- Box size within limits
- Maximum force < threshold

### Stage 3: Structure Selection

Uses PyNEP FPS to select most diverse structures.

### Stage 4: DFT Calculations

Sets up and runs VASP calculations.

### Stage 5: Data Collection

Converts VASP results and adds to training.

### Stage 6: NEP Retraining

Retrains NEP with expanded dataset.

## Customization

### Adjust MD Parameters

For different exploration strategies:

**High temperature exploration:**
```bash
md_temp=1500    # More aggressive sampling
md_steps=100000 # Longer trajectories
```

**Multiple temperatures:**
```bash
# Run at several temperatures
for temp in 600 900 1200 1500; do
    # Modify run.in for each temperature
    # Run MD
done
```

### Adjust Selection Strategy

**More structures per iteration:**
```bash
sample_number=200
select_number=100
```

**Aggressive exploration:**
```bash
max_force=50    # Allow higher forces
min_dist=1.2    # Tighter distance threshold
```

## Multi-Element Systems

For systems with multiple species:

```bash
# Ensure proper group labels
gpumdkit.sh -addgroup POSCAR Li Y Cl

# Filter may need adjustment
min_dist=1.3    # Smaller for Li
box_limit=15    # Larger for multi-component
```

## Convergence Criteria

Stop active learning when:

1. **RMSE plateaus**: No significant improvement
2. **Few uncertain structures**: Most pass uncertainty threshold
3. **Validation metrics stable**: Test set performance stops improving

Example monitoring:
```bash
# Track RMSE over iterations
grep "RMSE" */loss.out
```

## Example Workflow

```bash
# Iteration 1
cd iteration_01
cp ../iter00/nep.txt ./
cp ../iter00/train.xyz ./
sbatch workflow_active_learning_dev.sh

# Monitor and wait for completion

# Iteration 2
cd ../iteration_02
cp ../iter01/nep.txt ./
cp ../iter01/train.xyz ./
# Adjust parameters if needed
sbatch workflow_active_learning_dev.sh

# Continue until converged
```

---

For script details, see `Scripts/workflow/workflow_active_learning_dev.sh`
