# Workflow Scripts

This directory contains automation scripts for high-throughput computational workflows, particularly for NEP model development and active learning cycles.

## Overview

Workflow scripts automate repetitive tasks in computational materials research:
- **Batch preprocessing**: Prepare multiple structures for calculations
- **SCF calculations**: Set up VASP/CP2K single-point energy calculations
- **MD sampling**: Prepare molecular dynamics simulations with GPUMD/LAMMPS
- **Active learning**: Iterative NEP model improvement workflows
- **Submission management**: Generate job submission scripts

These tools bridge the gap between structure preparation and production calculations.

---

## Quick Command Reference

Access workflow tools through `gpumdkit.sh` interactive mode:

```bash
gpumdkit.sh
# Select: 3) Workflow
```

| Function | Option Number | Purpose |
|----------|---------------|---------|
| SCF batch pretreatment | 301 | Prepare VASP/CP2K calculations |
| MD batch (GPUMD) | 302 | Prepare GPUMD MD simulations |
| MD batch (LAMMPS) | 303 | Prepare LAMMPS MD simulations |
| Active learning (dev) | 3 | Automated active learning cycles |

---

## Scripts

### workflow_active_learning_dev.sh

Implements automated active learning cycles for iterative NEP model improvement.

**Purpose:** Automates the full active learning pipeline:
1. Generate candidate structures (MD sampling)
2. Evaluate with current NEP model
3. Select uncertain structures
4. Run DFT calculations
5. Add to training set
6. Retrain NEP model
7. Repeat until convergence

**Status:** Development version - under active development

**Documentation:** For detailed information on active learning workflows, see the [Active Learning tutorial](../../docs/tutorials/workflow_active_learning.md) or the [End-to-End Workflows](#end-to-end-workflow-examples) section in the main tutorials guide.

**Use Case:** Systematic NEP model development with minimal manual intervention

---

## General Usage Guidelines

### When to Use Workflows

Use workflow scripts when you have:
- **Multiple structures** requiring the same calculation setup
- **Repetitive tasks** that follow a standard procedure
- **Batch calculations** needing systematic organization
- **Active learning** requiring iteration

### Workflow Selection Guide

| Scenario | Recommended Workflow | Purpose |
|----------|---------------------|---------|
| Initial training data | SCF batch (301) | Single-point DFT energies/forces |
| NEP validation | MD batch GPUMD (302) | Test NEP in MD simulations |
| Explore configuration space | MD batch LAMMPS (303) | Generate candidate structures |
| Iterative training | Active learning (3) | Systematic model improvement |

### Common Workflows

#### Workflow 1: Generate Initial Training Set

```bash
# 1. Prepare structures (POSCAR or XYZ files)
# 2. Run batch preprocessing
gpumdkit.sh
# Select: 3) Workflow
# Select: 301) SCF batch pretreatment

# Enter prefix: initial_training_set
# Place POTCAR, INCAR, KPOINTS in fp/ directory

# 3. Submit calculations
cd struct_fp_*/
# Submit each directory's job

# 4. After completion, collect results
gpumdkit.sh -out2xyz .
```

#### Workflow 2: Validate NEP with MD

```bash
# 1. Have trained NEP model (nep.txt)
# 2. Prepare test structures
# 3. Set up MD batch
gpumdkit.sh
# Select: 3) Workflow  
# Select: 302) MD sample batch (GPUMD)

# 4. Place run.in and nep.txt in md/ directory
# 5. Run simulations
# 6. Analyze results with plotting tools
```

#### Workflow 3: Active Learning Cycle

```bash
# See detailed tutorial in docs/tutorials/
# Basic steps:
# 1. Initial NEP model
# 2. Generate candidates with MD
# 3. Select uncertain structures
# 4. Run DFT
# 5. Retrain model
# 6. Iterate
```

## File Organization

### SCF Workflow Structure

```
work_directory/
├── struct_fp_prefix_001/
│   ├── POSCAR
│   ├── POTCAR -> ../fp/POTCAR
│   ├── INCAR -> ../fp/INCAR
│   ├── KPOINTS -> ../fp/KPOINTS
│   └── submit.sh
├── struct_fp_prefix_002/
│   └── ...
└── fp/
    ├── POTCAR
    ├── INCAR
    └── KPOINTS
```

### MD Workflow Structure

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

## Input File Requirements

### For SCF Batch (VASP)

**Required in `fp/` directory:**
- `POTCAR`: Pseudopotential file (match elements in structures)
- `INCAR`: VASP input parameters
- `KPOINTS`: K-point mesh

**Recommended INCAR settings:**
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

### For MD Batch (GPUMD)

**Required in `md/` directory:**
- `run.in`: GPUMD run parameters
- `nep.txt`: NEP model file

**Example run.in:**
```
velocity        300
time_step       1
ensemble        npt_ber 300 300 100 0 0 0 1000
dump_thermo     1000
dump_position   1000
run             100000
```

### For MD Batch (LAMMPS)

**Required in `lmp/` directory:**
- `in.lammps`: LAMMPS input script
- Potential files (if applicable)

## Advanced Usage

### Custom Submission Scripts

The workflows generate `presub.sh` which creates individual `submit.sh` files. Customize the template in workflow scripts:

```bash
# Edit template in workflow script:
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

### Parallel Job Submission

After preprocessing:

```bash
# Submit all jobs in parallel
for dir in struct_fp_*/; do
    cd $dir
    sbatch submit.sh
    cd ..
done

# Or use GNU parallel
parallel -j 10 'cd {} && sbatch submit.sh && cd ..' ::: struct_fp_*/
```

### Monitoring Progress

```bash
# Check job status
squeue -u $USER

# Monitor specific calculations
for dir in struct_fp_*/; do
    echo "$dir: $(tail -1 $dir/OUTCAR 2>/dev/null | grep -o 'F=.*')"
done
```

## Best Practices

1. **Test first**: Run single calculation before batch to verify setup
2. **Unique prefixes**: Use descriptive prefixes with dates (e.g., `LiYCl_iter01_20250115`)
3. **Check symmetry**: Ensure POTCAR order matches POSCAR elements
4. **Resource management**: Don't overload cluster with too many simultaneous jobs
5. **Backup**: Keep copies of original structures before preprocessing
6. **Track provenance**: Document workflow parameters and iterations

## Integration with Other Tools

### With Format Conversion

```bash
# 1. Convert structures
gpumdkit.sh -cif2exyz structure.cif

# 2. Split into individual files
python split_single_xyz.py structures.xyz

# 3. Convert to POSCAR
for xyz in model_*.xyz; do
    gpumdkit.sh -exyz2pos $xyz
done

# 4. Run SCF workflow
gpumdkit.sh  # Select workflow option
```

### With Sampling

```bash
# 1. Sample from large MD trajectory
python sample_structures.py md.xyz uniform 100

# 2. Prepare sampled structures for DFT
gpumdkit.sh  # SCF batch pretreatment

# 3. Run calculations and add to training set
```

### With Analysis

```bash
# After calculations complete
# 1. Convert results
gpumdkit.sh -out2xyz ./struct_fp_*/

# 2. Analyze training set
gpumdkit.sh -analyze_comp train.xyz
gpumdkit.sh -range train.xyz force

# 3. Clean if needed
gpumdkit.sh -filter_value train.xyz force 30
```

## Troubleshooting

### Issue: "No .vasp files found"

**Problem**: Workflow expects VASP format but has XYZ

**Solution**:
- Workflow auto-converts XYZ if no VASP files present
- Or manually convert first: `gpumdkit.sh -exyz2pos structures.xyz`

### Issue: "POTCAR element mismatch"

**Problem**: POTCAR doesn't match elements in POSCAR

**Solution**:
1. Check POSCAR comment line for element symbols
2. Regenerate POTCAR with correct element order
3. Use `vaspkit` or similar tool to create proper POTCAR

### Issue: "Symbolic links broken"

**Problem**: Input files not found in calculation directories

**Solution**:
1. Check `fp/` or `md/` directory exists
2. Verify input files present before running workflow
3. Re-run preprocessing if links broken

### Issue: "presub.sh doesn't work"

**Problem**: Submission script template doesn't match cluster

**Solution**:
1. Edit workflow script to match your cluster's job scheduler
2. Modify resource requests (nodes, time, partition)
3. Update module loading commands

## Performance Tips

1. **Batch size**: Process 50-100 structures at a time for manageable submission
2. **Resource allocation**: Match resources to calculation requirements
3. **Parallel I/O**: Use separate directories to avoid I/O bottlenecks
4. **Checkpointing**: For long MD runs, enable restart capabilities
5. **Monitoring**: Set up automated monitoring for large batches

## Active Learning Considerations

### Convergence Criteria

Stop active learning when:
- Force RMSE plateaus
- No structures exceed uncertainty threshold
- Validation metrics stop improving
- Computational budget exhausted

### Structure Selection Strategy

Balance between:
- **Diversity**: Cover configuration space broadly
- **Uncertainty**: Focus on model weaknesses
- **Feasibility**: Structures must be DFT-calculable

### Iteration Planning

```
Iteration 1: ~500 diverse initial structures
Iteration 2: +100 structures from active selection
Iteration 3: +100 structures from refined selection
...
Continue until convergence
```

## Customization

### Adapting for Different Codes

The workflow scripts can be adapted for other codes:

**For Quantum ESPRESSO:**
```bash
# Modify SCF workflow to generate QE input files
# Update submission templates
# Adjust output parsing
```

**For CP2K:**
```bash
# CP2K workflow already implemented
# Use scf_batch_pretreatment_cp2k.py
```

**For ABACUS:**
```bash
# Similar structure to VASP workflow
# Adjust input file generation
```

## Dependencies

### Required
- **Bash**: Shell scripting
- **Python 3.x**: Preprocessing scripts
- **ASE**: Structure manipulation

### Job Scheduler
- SLURM (default templates)
- PBS/Torque (modify templates)
- LSF (modify templates)

### DFT Codes
- VASP (primary support)
- CP2K (supported)
- Others (customization needed)

## Contributing

To add new workflow capabilities:

1. **Follow structure**: Use existing workflows as templates
2. **Modular design**: Separate preprocessing, execution, and analysis
3. **Flexible parameters**: Allow user customization
4. **Error handling**: Validate inputs before processing
5. **Documentation**: Provide clear usage instructions and examples

See [CONTRIBUTING.md](../../CONTRIBUTING.md) for detailed guidelines.

---

Thank you for using GPUMDkit! For questions about workflows or to share your workflow adaptations, please open an issue on our [GitHub repository](https://github.com/zhyan0603/GPUMDkit/issues) or contact Zihan YAN (yanzihan@westlake.edu.cn).



#### scf_batch_pretreatment.sh

---

This script automates the preprocessing of `POSCAR` or `extxyz` files for *self-consistent field* (`SCF`) calculations. The script includes the following steps:

1. Converts a `.xyz` file to `POSCAR` format using `GPUMDkit` if no `.vasp` files are found in the current directory.
2. Renames and organizes `.vasp` files into a `struct_fp` directory.
3. Creates individual directories for each `POSCAR` file, setting up symbolic links to the necessary `VASP` input files.
4. Generates a `presub.sh` script to automate running `VASP` `SCF` calculations.

#### Usage

1. Prepare the environment:

   Ensure all `.vasp` files or a single `.xyz` file are in the current directory.

2. Enter:

   ```bash
   bash scf_batch_pretreatment.sh
   ```

3. You will see the following prompt: 

   ```sh
    Starting SCF batch pretreatment...
    Found 8 .vasp files.
    >-------------------------------------------------<
    | This function calls the script in Scripts       |
    | Script: scf_batch_pretreatment.sh               |
    | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
    >-------------------------------------------------<
   
    We recommend using the prefix to locate the structure.
    The folder name will be added to the second line of XYZ.
    config_type=<prefix>_<ID>
    ------------>>
    Please enter the prefix of directory (e.g. FAPBI3_iter01)
   ```

4. Enter the `prefix` of the folder name:

   ```sh
   FAPBI3_iter01
   ```

​		The script `scf_batch_pretreatment.sh` in the `Scripts` will be called to perform the pretreatment.

 5. You will see the following prompts:

    ```
     >-----------------------------------------------------<
     ATTENTION: Place POTCAR, KPOINTS and INCAR in 'fp' Dir.
     ATTENTION: Place POTCAR, KPOINTS and INCAR in 'fp' Dir.
     ATTENTION: Place POTCAR, KPOINTS and INCAR in 'fp' Dir.
     >-----------------------------------------------------<
    ```


You need to prepare the `POTCAR`, `KPOINTS`, and `INCAR` files and place them in a directory named `fp`.



#### md_sample_batch_pretreatment_gpumd.sh

---

This script automates the preprocessing of `POSCAR` or `extxyz` files for MD sampling using `GPUMD`. 

1. If `.vasp` files are found in the current directory, it will convert them to `extxyz` format to prepare the `model.xyz` file for `GPUMD`. If `.vasp` files are not found, the `.xyz` file will be read and all frames in it will be split into a individual sample.
2. Renames and organizes `.xyz` files into a `struct_md` directory.
3. Creates individual directories for each `model.xyz` file, setting up symbolic links to the necessary `GPUMD` input files.
4. Generates a `presub.sh` script to automate running MD simulations.

#### Usage

1. Prepare the environment:

   Ensure all `.vasp` files or a single `.xyz` file are in the current directory.

2. Enter:

   ```bash
   bash md_sample_batch_pretreatment_gpumd.sh
   ```

3. You will see the following prompt: 

   ```sh
    Starting MD sample batch pretreatment...
    No .vasp files found, but found one XYZ file.
    Converting it to model.xyz using GPUMDkit...
    All frames from "NEP-dataset.xyz" have been split into individual model files.
    20 model.xyz files were generated.
    >-------------------------------------------------<
    | This function calls the script in Scripts       |
    | Script: md_sample_batch_pretreatment.sh         |
    | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
    >-------------------------------------------------<
   ```

4. You will see the following prompts:

   ````
   ```
   >-----------------------------------------------<
   ATTENTION: Place run.in and nep.txt in 'md' Dir. 
   ATTENTION: Place run.in and nep.txt in 'md' Dir. 
   ATTENTION: Place run.in and nep.txt in 'md' Dir. 
   >-----------------------------------------------<
   ```
   ````

You need to prepare the `run.in` and`nep.txt` files and place them in a directory named `md`.



---

Thank you for using `GPUMDkit`! If you have any questions or need further assistance, feel free to open an issue on our GitHub repository or contact Zihan YAN (yanzihan@westlake.edu.cn).
