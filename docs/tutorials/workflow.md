<div align="center">
  <h1>🔧 Workflow</h1>
    <p style="text-align: justify;">This directory (Scripts/workflow/) contains automation scripts for high-throughput computational workflows, particularly for NEP model development and active learning cycles.</p>
</div>

**Script Location:** `Scripts/workflow/`

This section covers the workflow automation tools in GPUMDkit (Interactive Mode - Function 3).

Workflow scripts automate repetitive tasks in computational materials research:

- **SCF Batch preprocessin**: Set up VASP/CP2K single-point energy calculations
- **MD sampling**: Prepare molecular dynamics simulations with GPUMD/LAMMPS

## Interactive Mode

```bash
gpumdkit.sh
# Select: 3) Workflow
```

You'll see the following menu:

```
 ------------>>
 301) SCF batch pretreatment
 302) MD sample batch pretreatment (gpumd)
 303) MD sample batch pretreatment (lmp)
 000) Return to the main menu
 ------------>>
```

---

### Option 301: SCF batch pretreatment

This script automates the preprocessing of `POSCAR` or `extxyz` files for *self-consistent field* (`SCF`) calculations. The script includes the following steps:

1. Converts a `.xyz` file to `POSCAR` format using `GPUMDkit` if no `.vasp` files are found in the current directory.
2. Renames and organizes `.vasp` files into a `struct_fp` directory.
3. Creates individual directories for each `POSCAR` file, setting up symbolic links to the necessary `VASP` input files.
4. Generates a `presub.sh` script to automate running `VASP` `SCF` calculations.

#### Usage

Ensure all `.vasp` files or a single `.xyz` file are in the current directory.

Select option `301` from the menu:

```
301
```

You will see the following prompt:

```
 ------------>>
 Input the function number:
 301
 ------------>>
 1) VASP scf batch pretreatment
 2) CP2K scf batch pretreatment
 ------------>>
 Input the function number:
```

choose `1`or `2`:

```
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

Enter the `prefix` of the folder name:

```
Li7La3Zr2O12_iter01
```

The script `scf_batch_pretreatment.sh` in the `Scripts` will be called to perform the pretreatment.

You will see the following prompts:

```
 >-----------------------------------------------------<
 ATTENTION: Place POTCAR, KPOINTS and INCAR in 'fp' Dir.
 ATTENTION: Place POTCAR, KPOINTS and INCAR in 'fp' Dir.
 ATTENTION: Place POTCAR, KPOINTS and INCAR in 'fp' Dir.
 >-----------------------------------------------------<
```

You need to prepare the `POTCAR`, `KPOINTS`, and `INCAR` files and place them in a directory named `fp`.

### Option 302: MD sample batch pretreatment (gpumd)

This script automates the preprocessing of `POSCAR` or `extxyz` files for MD sampling using `GPUMD`.

1. If `.vasp` files are found in the current directory, it will convert them to `extxyz` format to prepare the `model.xyz` file for `GPUMD`. If `.vasp` files are not found, the `.xyz` file will be read and all frames in it will be split into a individual sample.
2. Renames and organizes `.xyz` files into a `struct_md` directory.
3. Creates individual directories for each `model.xyz` file, setting up symbolic links to the necessary `GPUMD` input files.
4. Generates a `presub.sh` script to automate running MD simulations.

#### Usage

Ensure all `.vasp` files or a single `.xyz` file are in the current directory.

Select option `302` from the menu:

```
302
```

You will see the following prompt:

```
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

You will see the following prompts:

```
 >----------------------------------------------------<
 | ATTENTION: Place run_*.in and nep.txt in 'md' Dir. |
 >----------------------------------------------------<
 | You need to provide MD control parameter files in  |
 | the format run_*.in (e.g., run_1.in, run_2.in),    |
 | each corresponding to a sample (e.g., sample_1,    |
 | sample_2) for molecular dynamics simulations.      |
 >----------------------------------------------------<
 Thank you for using GPUMDkit. Have a great day!
```

You need to prepare the `run.in` and`nep.txt` files and place them in a directory named `md`.

### Option 303: MD sample batch pretreatment (lmp)

The same with option `302`.

## workflow_active_learning_dev.sh

Implements automated active learning cycles for iterative NEP model improvement.

**Purpose:** Automates the active learning pipeline:

1. Generate candidate structures (MD sampling)
2. Evaluate with current NEP model
3. Select uncertain structures
4. Run DFT calculations
5. Add to training set
6. Retrain NEP model

**Status:** Development version - **under active development**

**Documentation:** See detailed tutorial in [docs/tutorials/](tutorials.md)

---

## Contributing

To add new workflow capabilities, see [CONTRIBUTING.md](contributing.md) for detailed guidelines.

---

Thank you for using GPUMDkit! For questions about workflows or to share your workflow adaptations, please open an issue on our [GitHub repository](https://github.com/zhyan0603/GPUMDkit/issues) or contact Zihan YAN (yanzihan@westlake.edu.cn).
