# Active Learning Workflow

This tutorial explains how to implement automated active learning cycles for iterative NEP model development using GPUMDkit. Active learning is a powerful approach that systematically improves your NEP model by intelligently selecting the most informative structures for DFT calculations.

## Overview

Active learning iteratively improves NEP models through an automated cycle:

1. **Generate candidates** - Run MD simulations with current NEP to sample configuration space
2. **Identify uncertainty** - Select structures where NEP predictions are uncertain
3. **Validate with DFT** - Run accurate DFT calculations on selected structures
4. **Expand training** - Add new DFT data to training set
5. **Retrain model** - Train improved NEP with expanded dataset
6. **Repeat** - Continue until convergence

This workflow automates the entire cycle, significantly reducing manual effort and accelerating NEP development.

## Prerequisites

Before starting active learning, ensure you have:

**Software:**
- ✅ **GPUMD** installed and in PATH
- ✅ **VASP** (or CP2K/ABACUS) configured for DFT calculations
- ✅ **GPUMDkit** installed with proper environment variables
- ✅ **PyNEP** installed (`pip install pynep`) for structure selection

**Initial Data:**
- ✅ **Initial NEP model** (`nep.txt`) - Should have reasonable accuracy
- ✅ **Training dataset** (`train.xyz`) - Initial training structures
- ✅ **Test dataset** (`test.xyz`) - For validation
- ✅ **Model structure** (`model.xyz`) - Reference structure with group labels

**Cluster Setup:**
- ✅ **Job scheduler** (SLURM, PBS, or similar) configured
- ✅ **Submission scripts** tested and working
- ✅ **Sufficient compute allocation** for DFT calculations

## Active Learning Workflow Stages

The active learning cycle consists of six automated stages:

### Stage 1: MD Sampling

**Purpose:** Generate diverse candidate structures by running MD simulations at elevated temperature.

**Process:**
1. Sample initial structures from current training set
2. Run short MD simulations (~50k-100k steps)
3. Collect trajectory frames at regular intervals
4. Pool all MD snapshots as candidates

**Key Parameters:**
```bash
md_steps=50000          # MD steps per simulation
md_temp=1200           # Temperature (K) - higher = more exploration
sample_number=100      # Number of candidate structures
```

**Output:** `candidates.xyz` with ~100-500 structures

---

### Stage 2: Structure Filtering

**Purpose:** Remove unphysical or problematic structures before selection.

**Quality Filters Applied:**

1. **Minimum distance filter**
   ```bash
   gpumdkit.sh -min_dist_pbc candidates.xyz
   # Remove structures with overlapping atoms
   ```

2. **Box size filter**
   ```bash
   gpumdkit.sh -filter_box candidates.xyz $box_limit
   # Remove structures with cells too small
   ```

3. **Maximum force filter**
   ```bash
   gpumdkit.sh -filter_value candidates.xyz force $max_force
   # Remove high-force unphysical structures
   ```

**Key Parameters:**
```bash
min_dist=1.4           # Minimum atomic distance (Å)
box_limit=13           # Minimum box dimension (Å)
max_force=30           # Maximum force threshold (eV/Å)
```

**Output:** `filtered_candidates.xyz` with physically reasonable structures

---

### Stage 3: Structure Selection

**Purpose:** Select most diverse and informative structures using farthest point sampling (FPS) in descriptor space.

**Method:** PyNEP FPS

```bash
gpumdkit.sh -pynep filtered_candidates.xyz selected.xyz nep.txt
```

**How it works:**
1. Calculate NEP descriptors for each structure
2. Use FPS algorithm to maximize diversity
3. Select structures that cover unexplored regions

**Key Parameters:**
```bash
select_number=50       # Structures to select for DFT
```

**Output:** `selected.xyz` with most diverse ~50 structures

**Alternative:** Can also use uncertainty-based selection with GPUMD active mode

---

### Stage 4: DFT Calculation Setup

**Purpose:** Organize and submit DFT calculations for selected structures.

**Process:**

1. **Directory organization**
   ```bash
   # Creates structure directories:
   struct_fp_${prefix}_001/
   struct_fp_${prefix}_002/
   ...
   ```

2. **Link input files**
   ```bash
   fp/
   ├── POTCAR    # Pseudopotentials
   ├── INCAR     # DFT parameters
   └── KPOINTS   # K-point mesh
   ```

3. **Generate submission scripts**
   ```bash
   # Each directory gets submit.sh
   # Configured for cluster scheduler
   ```

4. **Submit jobs**
   ```bash
   for dir in struct_fp_*/; do
       cd $dir && sbatch submit.sh && cd ..
   done
   ```

**Key Parameters:**
```bash
prefix_name=LiF_iter01     # Iteration identifier
dft_partition=intel-sc3    # Cluster partition
dft_nodes=1                # Nodes per calculation
```

---

### Stage 5: Data Collection

**Purpose:** Collect DFT results and add to training set.

**Process:**

1. **Wait for DFT completion**
   ```bash
   # Monitor job status
   squeue -u $USER
   ```

2. **Convert VASP outputs**
   ```bash
   gpumdkit.sh -out2xyz ./struct_fp_${prefix}_*/
   # Generates train.xyz with new structures
   ```

3. **Merge with existing data**
   ```bash
   cat train_original.xyz train_new.xyz > train_expanded.xyz
   ```

4. **Quality check**
   ```bash
   gpumdkit.sh -range train_expanded.xyz force
   gpumdkit.sh -min_dist_pbc train_expanded.xyz
   ```

**Output:** Expanded training set ready for retraining

---

### Stage 6: NEP Retraining

**Purpose:** Train improved NEP model with expanded dataset.

**Process:**

1. **Update training files**
   ```bash
   cp train_expanded.xyz train.xyz
   cp test.xyz .
   cp model.xyz .
   cp nep.in .
   ```

2. **Train NEP**
   ```bash
   nep
   ```

3. **Evaluate improvement**
   ```bash
   # Check RMSE reduction
   grep "RMSE" loss.out
   
   # Plot results
   gpumdkit.sh -plt train
   gpumdkit.sh -plt prediction
   ```

4. **Prepare for next iteration**
   ```bash
   # If not converged, copy nep.txt to next iteration
   cp nep.txt ../iteration_02/
   ```

---

## Setup and Configuration

### Step 1: Prepare Working Directory

Create organized directory structure for active learning iterations:

```bash
# Create iteration directory
mkdir -p iteration_01
cd iteration_01

# Copy essential files
cp ../iteration_00/nep.txt .        # Current NEP model
cp ../iteration_00/train.xyz .      # Training data
cp ../iteration_00/test.xyz .       # Test data
cp ../iteration_00/model.xyz .      # Model with group labels
cp ../iteration_00/nep.in .         # NEP training parameters

# Create MD configuration
cat > run.in << EOF
velocity        ${md_temp}
time_step       1
ensemble        npt_ber ${md_temp} ${md_temp} 100 0 0 0 1000
dump_thermo     1000
dump_position   1000
run             ${md_steps}
EOF
```

---

### Step 2: Configure Workflow Script

The `workflow_active_learning_dev.sh` script contains all parameters. Key sections to edit:

#### Basic Settings
```bash
# Iteration identification
prefix_name="LiF_iter01"           # Unique identifier for this iteration

# Quality control thresholds
min_dist=1.4                       # Minimum atomic distance (Å)
                                   # Li-Li: 1.4-1.5, O-O: 2.0, Metal-Metal: 2.0
box_limit=13                       # Minimum box dimension (Å)
                                   # Should be > 2 × NEP cutoff
max_force=30                       # Maximum force threshold (eV/Å)
                                   # Typical: 20-50 eV/Å depending on system
```

#### Sampling Parameters
```bash
# Candidate generation
sample_number=100                  # MD candidate structures to generate
                                   # More = better coverage, longer runtime
                                   # Typical: 50-200

select_number=50                   # Structures to select for DFT
                                   # Balance: cost vs improvement
                                   # Typical: 20-100 per iteration
```

#### MD Parameters
```bash
# Molecular dynamics settings
md_steps=50000                     # Steps per MD run
                                   # Longer = better exploration
                                   # Typical: 50k-200k

md_temp=1200                       # Temperature (K)
                                   # Higher = more aggressive exploration
                                   # Typical: 1.2-2.0 × target temp
```

#### DFT Cluster Settings
```bash
# Cluster configuration (adapt to your system)
dft_partition="intel-sc3"          # Partition/queue name
dft_queue="huge"                   # Queue type
dft_nodes=1                        # Nodes per calculation
dft_cores=48                       # Cores per node
dft_time="24:00:00"               # Walltime limit
```

---

### Step 3: Prepare DFT Input Files

Create `fp/` directory with VASP input templates:

```bash
mkdir -p fp
cd fp
```

#### INCAR (DFT parameters)
```
# fp/INCAR
PREC = Accurate
ENCUT = 520
ALGO = Normal
ISMEAR = 0
SIGMA = 0.05
EDIFF = 1E-6
LREAL = Auto

# Single-point calculation
NSW = 0
IBRION = -1

# Output
LWAVE = .FALSE.
LCHARG = .FALSE.

# For forces and stress
ISIF = 2
```

#### KPOINTS (k-point mesh)
```
# fp/KPOINTS - Example for ~200 atom cell
Automatic mesh
0
Gamma
2 2 2
0 0 0
```

**Guidelines for k-points:**
- Small cells (< 100 atoms): 4×4×4 or denser
- Medium cells (100-200 atoms): 2×2×2
- Large cells (> 200 atoms): Gamma-point (1×1×1)

#### POTCAR
```bash
# Concatenate POTCAR files for each element
cat ~/VASP/potpaw_PBE/Li/POTCAR \
    ~/VASP/potpaw_PBE/Y/POTCAR \
    ~/VASP/potpaw_PBE/Cl/POTCAR > POTCAR
```

**Important:** Element order must match structure files!

---

### Step 4: Run the Workflow

#### Option A: Submit to Cluster (Recommended)

```bash
# Make script executable
chmod +x workflow_active_learning_dev.sh

# Submit with SLURM
sbatch workflow_active_learning_dev.sh

# Monitor progress
squeue -u $USER
tail -f workflow.log
```

#### Option B: Run Locally (Small systems)

```bash
# Run in background
nohup bash workflow_active_learning_dev.sh &> workflow.log &

# Monitor
tail -f workflow.log
```

#### Option C: Interactive Testing (Development)

```bash
# Run directly (for debugging)
bash workflow_active_learning_dev.sh
```

---

### Step 5: Monitor Progress

**Check workflow status:**
```bash
# View log file
tail -f workflow.log

# Check for errors
grep -i "error\|fail" workflow.log

# Check DFT job status
squeue -u $USER | grep struct_fp
```

**Monitor key metrics:**
```bash
# Number of structures at each stage
echo "Candidates: $(grep -c "^[0-9]" candidates.xyz)"
echo "Filtered: $(grep -c "^[0-9]" filtered_candidates.xyz)"
echo "Selected: $(grep -c "^[0-9]" selected.xyz)"

# DFT completion
completed=$(find struct_fp_* -name "OUTCAR" -exec grep -l "reached" {} \; | wc -l)
total=$(ls -d struct_fp_* | wc -l)
echo "DFT progress: $completed / $total"
```

---

## Customization Strategies

### Adjust Exploration Strategy

#### Conservative Approach (Early Iterations)
```bash
md_temp=1000               # Moderate temperature
md_steps=50000            # Standard length
sample_number=100         # Good coverage
select_number=30          # Smaller batches
max_force=25              # Strict filtering
```

**When to use:** First few iterations, when NEP is still learning basics

---

#### Aggressive Exploration (Middle Iterations)
```bash
md_temp=1500              # High temperature
md_steps=100000           # Longer sampling
sample_number=200         # Extensive candidates
select_number=80          # Larger batches
max_force=40              # Relaxed filtering
```

**When to use:** Mid-training, when exploring phase space

---

#### Refinement (Late Iterations)
```bash
md_temp=800               # Near-target temperature
md_steps=200000           # Long equilibration
sample_number=150         # Thorough sampling
select_number=50          # Moderate batches
max_force=20              # Strict quality
```

**When to use:** Final iterations, improving accuracy

---

### Multi-Temperature Sampling

Sample different temperature regimes in one iteration:

```bash
# Modify workflow to run multiple temperatures
temperatures=(600 900 1200 1500)

for temp in "${temperatures[@]}"; do
    # Update run.in with current temperature
    sed -i "s/velocity .*/velocity $temp/" run.in
    sed -i "s/ensemble npt_ber .*/ensemble npt_ber $temp $temp 100 0 0 0 1000/" run.in
    
    # Run MD with temperature-specific output
    gpumd
    mv dump.xyz dump_${temp}K.xyz
done

# Combine all temperatures
cat dump_*K.xyz > candidates_all_temp.xyz
```

**Benefits:** Better coverage of configuration space across temperatures

---

### System-Specific Adjustments

#### Li-ion Conductors
```bash
min_dist=1.3              # Li atoms can get close
md_temp=1200              # High mobility needed
max_force=35              # Allow flexible structures
```

#### Oxides
```bash
min_dist=1.8              # O-O repulsion
md_temp=1500              # Overcome barriers
max_force=30              # Standard filtering
```

#### Molecular Systems
```bash
min_dist=1.0              # Covalent bonds
md_temp=800               # Avoid dissociation
max_force=20              # Strict geometries
```

#### Alloys/Metals
```bash
min_dist=2.0              # Metal-metal distance
md_temp=1200              # Explore mixing
max_force=25              # Moderate filtering
```

---

### Cluster Optimization

#### For Large Cells (> 500 atoms)
```bash
# INCAR modifications
LREAL = Auto              # Real-space projection
NCORE = 12                # Core parallelization
KPAR = 4                  # K-point parallelization

# Increase resources
dft_nodes=2
dft_cores=96
dft_time="48:00:00"
```

#### For Many Small Calculations
```bash
# Run multiple jobs per node
dft_cores=12              # Reduce cores per job
dft_nodes=1
# Submit multiple simultaneously

# Or use job arrays
#SBATCH --array=1-50
```

---

## Complete Workflow Example

Here's a step-by-step example of a full active learning workflow for a Li-Y-Cl solid electrolyte:

### Iteration 0: Initial Training

```bash
# Start with small initial dataset (~500-1000 structures)
cd iteration_00

# Prepare training data
gpumdkit.sh -out2xyz ./vasp_calculations/
gpumdkit.sh -filter_value train.xyz force 30
gpumdkit.sh -min_dist_pbc train.xyz
gpumdkit.sh -addgroup POSCAR Li Y Cl

# Train initial NEP
nep

# Check initial accuracy
gpumdkit.sh -plt train
gpumdkit.sh -plt prediction
# Initial RMSE: Energy ~5 meV/atom, Force ~150 meV/Å
```

---

### Iteration 1: First Active Learning Cycle

```bash
cd ../iteration_01

# Copy from iteration 0
cp ../iteration_00/nep.txt .
cp ../iteration_00/train.xyz .
cp ../iteration_00/test.xyz .
cp ../iteration_00/model.xyz .

# Configure workflow (conservative first iteration)
cat > workflow_config.sh << 'EOF'
prefix_name="LiYCl_iter01"
min_dist=1.4
box_limit=15
max_force=30
sample_number=100
select_number=30
md_steps=50000
md_temp=1000
EOF

# Run active learning
sbatch workflow_active_learning_dev.sh

# Monitor progress
watch -n 60 'squeue -u $USER; echo "---"; \
             ls struct_fp_*/OUTCAR 2>/dev/null | wc -l'

# After completion (~ 1-2 days)
# New training set: 500 + 30 = 530 structures
# Retrain NEP
nep

# Check improvement
gpumdkit.sh -plt train
# Improved RMSE: Energy ~3 meV/atom, Force ~100 meV/Å
```

---

### Iteration 2: Aggressive Exploration

```bash
cd ../iteration_02

# Copy improved model
cp ../iteration_01/nep.txt .
cp ../iteration_01/train.xyz .
cp ../iteration_01/test.xyz .
cp ../iteration_01/model.xyz .

# More aggressive parameters
cat > workflow_config.sh << 'EOF'
prefix_name="LiYCl_iter02"
min_dist=1.3
box_limit=13
max_force=35
sample_number=150
select_number=50
md_steps=100000
md_temp=1200
EOF

# Run workflow
sbatch workflow_active_learning_dev.sh

# After completion
# Training set: 530 + 50 = 580 structures
nep

# Check improvement
gpumdkit.sh -plt train
# Further improved: Energy ~2 meV/atom, Force ~70 meV/Å
```

---

### Iteration 3: Multi-Temperature Refinement

```bash
cd ../iteration_03

cp ../iteration_02/nep.txt .
cp ../iteration_02/train.xyz .
cp ../iteration_02/test.xyz .
cp ../iteration_02/model.xyz .

# Sample multiple temperatures
for temp in 600 900 1200 1500; do
    mkdir md_${temp}K
    cd md_${temp}K
    
    cat > run.in << EOF
velocity        $temp
time_step       1
ensemble        npt_ber $temp $temp 100 0 0 0 1000
dump_thermo     1000
dump_position   1000
run             100000
EOF
    
    # Run sampling
    gpumd
    cd ..
done

# Combine and process
cat md_*/dump.xyz > candidates_all.xyz
gpumdkit.sh -filter_value candidates_all.xyz force 30
gpumdkit.sh -pynep filtered_force.xyz selected.xyz nep.txt

# Continue with DFT...
# Final: 580 + 60 = 640 structures
# RMSE: Energy ~1.5 meV/atom, Force ~50 meV/Å
```

---

### Iteration 4: Convergence Check

```bash
cd ../iteration_04

# Final refinement with strict quality control
# ... (similar to above)

# Training set: ~700 structures
# RMSE: Energy ~1.2 meV/atom, Force ~45 meV/Å

# Check convergence
echo "Iteration | Energy RMSE | Force RMSE | Structures"
echo "----------|-------------|------------|------------"
echo "    0     |     5.0     |    150     |    500"
echo "    1     |     3.0     |    100     |    530"
echo "    2     |     2.0     |     70     |    580"
echo "    3     |     1.5     |     50     |    640"
echo "    4     |     1.2     |     45     |    700"

# Convergence achieved! Stop iterations.
```

---

## Convergence Criteria

### Quantitative Metrics

**Stop active learning when:**

1. **RMSE Plateau** (Primary criterion)
   ```bash
   # Energy RMSE change < 10% between iterations
   ΔE_RMSE < 0.1 * E_RMSE_prev
   
   # Force RMSE change < 5%
   ΔF_RMSE < 0.05 * F_RMSE_prev
   ```

2. **Target Accuracy Reached**
   ```
   Energy RMSE < 2 meV/atom   (typical target)
   Force RMSE < 50 meV/Å      (typical target)
   Stress RMSE < 0.5 GPa      (if trained)
   ```

3. **Validation Performance**
   ```bash
   # Test set RMSE similar to training RMSE
   |E_RMSE_train - E_RMSE_test| < 0.5 meV/atom
   |F_RMSE_train - F_RMSE_test| < 10 meV/Å
   ```

---

### Qualitative Checks

**Also verify:**

1. **Structure Selection Rate**
   ```bash
   # Fewer structures needed as model improves
   # If < 5% of candidates selected → good coverage
   selected / total_candidates < 0.05
   ```

2. **Physical Properties**
   ```bash
   # Run test MD simulations
   # Check: lattice constant, density, diffusion
   # Compare with experiment or reference
   ```

3. **Descriptor Coverage**
   ```bash
   # Visualize descriptor space
   gpumdkit.sh -calc des umap train.xyz desc.npy nep.txt Li
   gpumdkit.sh -plt des umap desc.npy
   # Ensure good coverage, no gaps
   ```

---

### Monitoring Script

Create a convergence tracking script:

```bash
#!/bin/bash
# monitor_convergence.sh

echo "Iteration Analysis"
echo "=================="

for iter in iteration_*/; do
    cd $iter
    if [ -f loss.out ]; then
        e_rmse=$(grep "Energy" loss.out | tail -1 | awk '{print $4}')
        f_rmse=$(grep "Force" loss.out | tail -1 | awk '{print $4}')
        n_struct=$(grep -c "^[0-9]" train.xyz 2>/dev/null || echo "0")
        echo "$iter: E=$e_rmse meV/atom, F=$f_rmse meV/Å, N=$n_struct"
    fi
    cd ..
done
```

**Usage:**
```bash
bash monitor_convergence.sh
```

**Example output:**
```
Iteration Analysis
==================
iteration_00: E=5.2 meV/atom, F=148 meV/Å, N=500
iteration_01: E=3.1 meV/atom, F=102 meV/Å, N=530
iteration_02: E=2.0 meV/atom, F=68 meV/Å, N=580
iteration_03: E=1.4 meV/atom, F=49 meV/Å, N=640
iteration_04: E=1.2 meV/atom, F=46 meV/Å, N=700
Convergence achieved at iteration 4!
```

---

## Troubleshooting

### Issue: Too Few Structures Selected

**Symptoms:**
```bash
# Only 5-10 structures selected instead of 50
echo "Selected: $(grep -c "^[0-9]" selected.xyz)"
# Output: Selected: 8
```

**Causes & Solutions:**

1. **Filters too strict**
   ```bash
   # Relax thresholds
   max_force=40          # Increase from 30
   min_dist=1.3          # Decrease from 1.5
   box_limit=12          # Decrease from 15
   ```

2. **Temperature too low**
   ```bash
   # Increase exploration
   md_temp=1500          # Increase from 1000
   md_steps=150000       # Longer sampling
   ```

3. **NEP model too good**
   ```bash
   # Check if already converged
   gpumdkit.sh -plt prediction
   # If RMSE very low, may have converged!
   ```

---

### Issue: DFT Jobs Fail

**Symptoms:**
```bash
# Check for failures
grep -r "error\|ERROR" struct_fp_*/OUTCAR
```

**Common Causes:**

1. **POTCAR mismatch**
   ```bash
   # Solution: Verify element order
   head -1 struct_fp_001/POSCAR  # Check elements
   head -5 fp/POTCAR              # Check POTCAR order
   ```

2. **Insufficient memory**
   ```bash
   # Solution: Increase INCAR settings
   NCORE = 8              # Reduce from 12
   LREAL = Auto           # Add if missing
   ```

3. **k-points too dense**
   ```bash
   # Solution: For large cells, use Gamma
   # fp/KPOINTS
   Gamma
   0
   Gamma
   1 1 1
   ```

4. **Electronic convergence**
   ```bash
   # Solution: Relax INCAR
   ALGO = Fast            # Instead of Normal
   NELM = 200             # Increase from default
   ```

---

### Issue: MD Doesn't Generate Diverse Structures

**Symptoms:**
```bash
# All structures very similar
gpumdkit.sh -calc des umap candidates.xyz desc.npy nep.txt Li
gpumdkit.sh -plt des umap desc.npy
# Plot shows tight cluster
```

**Solutions:**

1. **Use multiple starting points**
   ```bash
   # Sample diverse initial structures
   python ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py \
       train.xyz uniform 20
   # Use these as starting points
   ```

2. **Increase temperature**
   ```bash
   md_temp=1800           # Very high temperature
   ```

3. **Try different conditions**
   ```bash
   # Vary pressure
   ensemble npt_ber 1200 1200 100 1 1 1 1000  # With pressure
   
   # Vary box shape
   ensemble npt_ber 1200 1200 100 1 1 0 1000  # Flexible XY
   ```

---

### Issue: Training Set Becomes Too Large

**Symptoms:**
```bash
# Training takes too long
wc -l train.xyz
# Output: 50000 lines (1000+ structures)
```

**Solutions:**

1. **Prune redundant structures**
   ```bash
   # Use FPS to reduce training set
   gpumdkit.sh -pynep train.xyz train_reduced.xyz nep.txt
   # Keep most diverse 70%
   ```

2. **Remove outliers**
   ```bash
   # Remove high-error structures
   gpumdkit.sh
   # Select: 5) Analyzer → 502
   # Use remained.xyz
   ```

3. **Split by composition**
   ```bash
   # Analyze and balance
   gpumdkit.sh -analyze_comp train.xyz
   # Export each composition separately
   # Keep equal numbers from each
   ```

---

### Issue: Workflow Script Hangs

**Symptoms:**
```bash
# No progress for hours
tail -f workflow.log
# Last line shows waiting...
```

**Diagnosis:**

1. **Check running jobs**
   ```bash
   squeue -u $USER
   # See if DFT jobs are stuck
   ```

2. **Check disk space**
   ```bash
   df -h .
   # Ensure sufficient space
   ```

3. **Check file permissions**
   ```bash
   ls -la struct_fp_*/
   # Verify write permissions
   ```

**Solutions:**
```bash
# Kill and restart
scancel <jobid>
# Fix issue
# Resubmit workflow
```

---

### Issue: Results Not Reproducible

**Symptoms:**
- Different results on re-run
- Structures change between iterations

**Causes & Solutions:**

1. **Random seed not set**
   ```bash
   # In run.in, add:
   random_seed 123456
   ```

2. **Different VASP versions**
   ```bash
   # Document version
   vasp_std --version
   # Use consistent version
   ```

3. **Floating point precision**
   ```bash
   # Normal for MD - accept small variations
   # Verify overall statistics match
   ```

---

### Issue: Poor Final Accuracy

**Symptoms:**
```bash
# After 5+ iterations, still high RMSE
# Force RMSE > 100 meV/Å
```

**Diagnosis Steps:**

1. **Check DFT quality**
   ```bash
   # Verify DFT convergence
   grep "reached" struct_fp_*/OUTCAR | wc -l
   # Should match number of calculations
   ```

2. **Check training parameters**
   ```bash
   # Review nep.in
   cat nep.in
   # Ensure appropriate cutoff, basis, etc.
   ```

3. **Check descriptor coverage**
   ```bash
   gpumdkit.sh -calc des umap train.xyz desc.npy nep.txt Li
   gpumdkit.sh -plt des umap desc.npy
   # Look for gaps or outliers
   ```

**Solutions:**

1. **Increase NEP complexity**
   ```
   # In nep.in
   n_max 10 8    # Increase from 8 6
   l_max 4 2 1   # Increase from 4 2 0
   ```

2. **Add targeted structures**
   ```bash
   # Identify problematic configurations
   # Add similar structures manually
   ```

3. **Clean training data**
   ```bash
   # Remove outliers
   gpumdkit.sh -plt force_errors
   # Identify and remove bad structures
   ```

---

## Best Practices

### 1. Data Organization

```bash
# Maintain clear directory structure
project/
├── iteration_00/          # Initial training
│   ├── train.xyz
│   ├── test.xyz
│   ├── model.xyz
│   ├── nep.in
│   └── nep.txt
├── iteration_01/          # First AL cycle
│   ├── workflow_active_learning_dev.sh
│   ├── workflow.log
│   ├── candidates.xyz
│   ├── selected.xyz
│   ├── struct_fp_iter01_001/
│   ├── ...
│   ├── train.xyz         # Expanded training
│   └── nep.txt           # Improved model
├── iteration_02/          # Second AL cycle
│   └── ...
└── analysis/              # Keep analysis separate
    ├── convergence.png
    └── descriptor_coverage.png
```

---

### 2. Documentation

**Keep a lab notebook** (can be markdown):

```markdown
# iteration_01/notes.md

## Date: 2025-01-15

### Configuration
- Temperature: 1000 K
- Selected: 30 structures
- DFT: VASP 6.3.0, PBE, 520 eV

### Results
- Energy RMSE: 3.1 meV/atom (was 5.2)
- Force RMSE: 102 meV/Å (was 148)
- Training time: 2.5 hours

### Observations
- Good improvement in force RMSE
- Some high-energy structures from MD
- Next iteration: increase temperature to 1200 K

### Files
- train.xyz: 530 structures (added 30)
- nep.txt: model saved
```

---

### 3. Validation

**Test on independent data:**

```bash
# After each iteration
# Run test MD simulations
mkdir validation
cd validation

# Test 1: Lattice constant
cat > run.in << EOF
velocity 300
ensemble npt_ber 300 300 100 0 0 0 1000
run 100000
EOF
gpumd
gpumdkit.sh -plt thermo
# Compare lattice with experiment

# Test 2: Diffusion
cat > run.in << EOF
velocity 600
ensemble npt_ber 600 600 100 0 0 0 1000
compute_msd 100 1
run 500000
EOF
gpumd
gpumdkit.sh -calc ionic-cond Li 1
# Compare with literature
```

---

### 4. Resource Management

**Optimize cluster usage:**

```bash
# Submit in batches, not all at once
# Prevents overwhelming scheduler
for i in {1..10}; do
    cd struct_fp_iter01_$(printf "%03d" $i)
    sbatch submit.sh
    cd ..
    sleep 2  # Small delay
done

# Monitor and resubmit failures
for dir in struct_fp_*/; do
    if ! grep -q "reached" $dir/OUTCAR 2>/dev/null; then
        echo "Resubmitting $dir"
        cd $dir && sbatch submit.sh && cd ..
    fi
done
```

---

### 5. Backup Strategy

```bash
# Regular backups
backup_dir="../backups/$(date +%Y%m%d)"
mkdir -p $backup_dir

# Backup essential files
cp train.xyz test.xyz nep.txt nep.in $backup_dir/
tar -czf $backup_dir/struct_fp.tar.gz struct_fp_*/POSCAR struct_fp_*/OUTCAR

# Keep logs
cp workflow.log $backup_dir/
```

---

### 6. Performance Tips

**Speed up iterations:**

1. **Parallel DFT**
   ```bash
   # Use job arrays
   #SBATCH --array=1-50
   # Submit all at once
   ```

2. **Fast DFT settings for exploration**
   ```
   # Early iterations - faster
   PREC = Normal
   ENCUT = 400
   EDIFF = 1E-5
   
   # Final iterations - accurate
   PREC = Accurate
   ENCUT = 520
   EDIFF = 1E-6
   ```

3. **Cache MD trajectories**
   ```bash
   # Save all MD dumps
   # Reuse if need different selection criteria
   ```

---

### 7. Quality Assurance

**Checklist before each iteration:**

- [ ] Previous NEP model tested and validated
- [ ] DFT inputs verified (POTCAR, INCAR, KPOINTS)
- [ ] Cluster resources available
- [ ] Disk space sufficient (> 100 GB recommended)
- [ ] Backup of previous iteration complete
- [ ] Workflow parameters documented
- [ ] Expected completion time estimated

---

---

For script details, see `Scripts/workflow/workflow_active_learning_dev.sh`
