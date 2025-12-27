# Workflow Automation

GPUMDkit provides workflow automation tools for batch processing of DFT calculations and MD simulations.

## Available Workflows

### SCF Batch Pretreatment (`scf_batch_pretreatment`)

Automates setup of VASP single-point energy calculations for multiple structures.

**Interactive:** Select option `301`

**What it does:**
1. Converts XYZ to POSCAR if needed
2. Organizes structures into individual directories
3. Creates symbolic links to VASP input files
4. Generates submission scripts

**Usage:**

1. Prepare structures (POSCAR files or single XYZ file)
2. Run workflow:
   ```bash
   gpumdkit.sh
   # Select: 3) Workflow
   # Select: 301
   ```
3. Enter directory prefix (e.g., `iter01_training`)
4. Place `POTCAR`, `INCAR`, `KPOINTS` in `fp/` directory
5. Submit calculations

**Output structure:**
```
work_directory/
├── struct_fp_iter01_001/
│   ├── POSCAR
│   ├── POTCAR -> ../fp/POTCAR
│   ├── INCAR -> ../fp/INCAR
│   ├── KPOINTS -> ../fp/KPOINTS
│   └── submit.sh
├── struct_fp_iter01_002/
│   └── ...
└── fp/
    ├── POTCAR
    ├── INCAR
    └── KPOINTS
```

**After calculations:**
```bash
# Convert results
gpumdkit.sh -out2xyz ./struct_fp_iter01_*/

# Output: train.xyz with all results
```

### MD Sample Batch Pretreatment - GPUMD (`md_sample_batch_pretreatment_gpumd`)

Automates setup of GPUMD MD simulations for multiple structures.

**Interactive:** Select option `302`

**What it does:**
1. Converts VASP to XYZ or splits XYZ file
2. Organizes structures into individual directories
3. Creates symbolic links to GPUMD input files
4. Generates submission scripts

**Usage:**

1. Prepare structures (POSCAR files or single XYZ file)
2. Run workflow:
   ```bash
   gpumdkit.sh
   # Select: 3) Workflow
   # Select: 302
   ```
3. Place `run.in` and `nep.txt` in `md/` directory
4. Submit simulations

**Output structure:**
```
work_directory/
├── struct_md_001/
│   ├── model.xyz
│   ├── run.in -> ../md/run.in
│   ├── nep.txt -> ../md/nep.txt
│   └── submit.sh
├── struct_md_002/
│   └── ...
└── md/
    ├── run.in
    └── nep.txt
```

### MD Sample Batch Pretreatment - LAMMPS (`md_sample_batch_pretreatment_lmp`)

Similar to GPUMD workflow but for LAMMPS simulations.

**Interactive:** Select option `303`

**Usage:**

1. Prepare structures
2. Run workflow and select option 303
3. Place LAMMPS input files in `lmp/` directory
4. Submit simulations

## Interactive Mode

```bash
gpumdkit.sh
# Select: 3) Workflow
```

### Menu Options

```
 ------------>>
 301) SCF batch pretreatment
 302) MD sample batch pretreatment (gpumd)
 303) MD sample batch pretreatment (lmp)
 000) Return to the main menu
 ------------>>
```

## Common Workflows

### Generate Training Data

```bash
# 1. Prepare candidate structures (multiple POSCAR files)

# 2. Set up SCF calculations
gpumdkit.sh
# Select: 3) Workflow
# Select: 301
# Enter prefix: training_set_01

# 3. Prepare VASP inputs
mkdir fp
# Copy POTCAR, INCAR, KPOINTS to fp/

# 4. Submit calculations
for dir in struct_fp_training_set_01_*/; do
    cd $dir
    sbatch submit.sh
    cd ..
done

# 5. After completion, collect results
gpumdkit.sh -out2xyz ./struct_fp_training_set_01_*/
```

### Validate NEP Model

```bash
# 1. Prepare test structures

# 2. Set up MD simulations
gpumdkit.sh
# Select: 3) Workflow
# Select: 302

# 3. Prepare GPUMD inputs
mkdir md
# Copy run.in and nep.txt to md/

# 4. Run simulations and analyze results
```

### Active Learning Iteration

```bash
# 1. Use select_max_modev to get uncertain structures
# 2. Run SCF batch workflow (301) on selected structures
# 3. Convert results and add to training set
# 4. Retrain NEP model
# 5. Repeat
```

## VASP Input Guidelines

### INCAR for SCF Calculations

```
PREC = Accurate
ENCUT = 520
ISMEAR = 0
SIGMA = 0.05
EDIFF = 1E-6
NSW = 0           # Single-point calculation
IBRION = -1       # No ionic relaxation
LWAVE = .FALSE.
LCHARG = .FALSE.
```

### KPOINTS

Adjust based on system size:
- Small cell (< 100 atoms): Dense mesh (e.g., 4×4×4)
- Large cell (> 100 atoms): Gamma-point or coarse mesh

## GPUMD Input Guidelines

### Example run.in

```
velocity        300
time_step       1
ensemble        npt_ber 300 300 100 0 0 0 1000
dump_thermo     1000
dump_position   1000
compute_msd     100 1 0
run             100000
```

## Tips

- **Unique prefixes**: Use descriptive prefixes with dates
- **Test first**: Run single calculation before batch
- **Resource management**: Don't overload cluster
- **Backup**: Keep copies of original structures
- **Track provenance**: Document workflow parameters

## Submission Script Customization

Edit the workflow scripts to match your cluster:

```bash
# Modify template in workflow script
cat > submit.sh << EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --time=24:00:00
#SBATCH --partition=your_partition

module load vasp
mpirun vasp_std
EOF
```

## Parallel Job Submission

```bash
# Submit all jobs
for dir in struct_fp_*/; do
    cd $dir
    sbatch submit.sh
    cd ..
done

# Or use GNU parallel
parallel -j 10 'cd {} && sbatch submit.sh && cd ..' ::: struct_fp_*/
```

## Monitoring Progress

```bash
# Check job status
squeue -u $USER

# Monitor specific calculations
for dir in struct_fp_*/; do
    echo "$dir: $(tail -1 $dir/OUTCAR 2>/dev/null | grep -o 'F=.*')"
done
```

---

For more details, see [Scripts/workflow/README.md](../../Scripts/workflow/README.md)
