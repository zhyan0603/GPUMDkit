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

**Documentation:** See detailed tutorial in [docs/tutorials/](../../docs/tutorials/)

**Use Case:** Systematic NEP model development with minimal manual intervention

---

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
