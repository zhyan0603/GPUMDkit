# Active Learning Workflow

This tutorial explains how to use the `workflow_active_learning_dev.sh` script for automated NEP model development through active learning.

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

Generates candidate structures by running MD at specified temperature:

```bash
# Multiple MD runs from sampled initial structures
for each structure in train.xyz:
    Run GPUMD MD simulation
    Collect trajectory
```

### Stage 2: Structure Filtering

Applies quality filters:
- Minimum atomic distance > threshold
- Box size within limits
- Maximum force < threshold

```bash
gpumdkit.sh -min_dist_pbc candidates.xyz
gpumdkit.sh -filter_box candidates.xyz $box_limit
gpumdkit.sh -filter_value candidates.xyz force $max_force
```

### Stage 3: Structure Selection

Uses PyNEP FPS to select most diverse structures:

```bash
gpumdkit.sh -pynep filtered_candidates.xyz selected.xyz nep.txt
```

### Stage 4: DFT Calculations

Sets up and runs VASP calculations:

```bash
# Use SCF batch workflow
# Submits jobs to cluster
# Waits for completion
```

### Stage 5: Data Collection

Converts VASP results and adds to training:

```bash
gpumdkit.sh -out2xyz ./struct_fp_*/
cat train.xyz dft_results.xyz > train_new.xyz
```

### Stage 6: NEP Retraining

Retrains NEP with expanded dataset:

```bash
nep  # Using updated train.xyz
```

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

### Cluster Settings

Adapt for your cluster:

```bash
# For different scheduler
#PBS -N workflow
#PBS -l nodes=1:ppn=48

# For different DFT code
module load quantum-espresso
pw.x < input.in > output.out
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

## Troubleshooting

### Issue: Too few structures selected

**Solution:** Relax filters
```bash
max_force=40    # Increase threshold
min_dist=1.3    # Decrease threshold
```

### Issue: DFT jobs fail

**Solution:** Check VASP inputs
```bash
# Verify POTCAR matches elements
# Check ENCUT is sufficient
# Ensure KPOINTS appropriate for cell size
```

### Issue: MD doesn't generate diverse structures

**Solution:** Adjust MD conditions
```bash
md_temp=1500    # Higher temperature
md_steps=100000 # Longer simulation
# Or use multiple starting configurations
```

## Example Complete Workflow

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

## Best Practices

1. **Start conservative**: Use strict filters initially
2. **Monitor carefully**: Check selected structures make sense
3. **Track provenance**: Keep detailed logs of each iteration
4. **Backup regularly**: Save NEP models and training data
5. **Validate frequently**: Test on independent validation set
6. **Document changes**: Note parameter adjustments between iterations

## Performance Tips

- Run multiple iterations in parallel if resources allow
- Use faster DFT settings for exploration, accurate for final
- Cache MD trajectories for analysis
- Profile bottlenecks (MD vs DFT vs selection)

---

For script details, see `Scripts/workflow/workflow_active_learning_dev.sh`
