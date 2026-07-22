<div align="center">
  <h1>🔧 Workflow Scripts</h1>
    <p style="text-align: justify;">This directory contains automation scripts for high-throughput computational workflows, particularly for NEP model development and active learning cycles.</p>
</div>

## Overview

Workflow scripts automate repetitive tasks in computational materials research:

- **SCF batch preprocessing**: Set up VASP/CP2K single-point energy calculations
- **MD sampling**: Prepare molecular dynamics simulations with GPUMD/LAMMPS
- **Active learning**: Iterative NEP model improvement workflows

---

## Via interactive mode

Access workflow tools through `gpumdkit.sh` interactive mode:

```text
 +---------------------------------------------------------+
 |                      WORKFLOW TOOLS                     |
 +---------------------------------------------------------+
 | 301) SCF batch pretreatment                             |
 | 302) MD sample batch pretreatment (gpumd)               |
 | 303) MD sample batch pretreatment (lmp)                 |
 +---------------------------------------------------------+
 | 000) Return to the main menu                            |
 +---------------------------------------------------------+
 Input the function number:
```

---

## Workflow Entries

### 301) SCF batch pretreatment

This entry prepares folders for single-point calculations. For VASP, keep either `.vasp` structures or one `.xyz` file in the current directory before running the script.

- If `.vasp` files are present, the script processes the `.vasp` files.
- If both `.vasp` and `.xyz` files are present, the script prints a notice and only processes `.vasp` files.
- If no `.vasp` file is present and multiple `.xyz` files are detected, the script asks which `.xyz` file to process.

After running, the script creates `struct_fp/`, `fp/`, and `fp_sample_*` directories. Put `INCAR`, `POTCAR`, and `KPOINTS` into the generated `fp/` directory; each calculation folder links to these files.

For CP2K, enter `3) Workflow`, then `301) SCF batch pretreatment`, and choose the CP2K branch. The CP2K script asks for:

```text
<input.xyz> <template.inp> <prefix>
```

The repository provides `Scripts/workflow/cp2k_template.inp` as a starting template.



### 302) MD sample batch pretreatment (GPUMD)

This entry prepares GPUMD MD sampling folders. Start from `.vasp` structures, or from one selected `.xyz` trajectory/structure file if no `.vasp` files are present.

- If `.vasp` files are present, they are converted to extxyz files and processed.
- If both `.vasp` and `.xyz` files are present, only `.vasp` files are processed.
- If multiple `.xyz` files are detected without `.vasp` files, the script asks which one to use.

After running, put `nep.txt` and `run_*.in` files into the generated `md/` directory. The sample folders link `run_1.in`, `run_2.in`, ... as their `run.in`.



### 303) MD sample batch pretreatment (LAMMPS)

This entry prepares LAMMPS MD sampling folders. The input-file selection rules are the same as function 302.

After running, put `lmprun.in` and `nep.txt` into the generated `md/` directory. The sample folders link `lammps.data`, `lmprun.in`, and `nep.txt`.


### workflow_active_learning_dev.sh

Implements automated active learning cycles for iterative NEP model improvement.

**Purpose:** Automates the active learning pipeline:

1. Generate candidate structures (MD sampling)
2. Evaluate with current NEP model
3. Select uncertain structures
4. Run DFT calculations
5. Add to training set
6. Retrain NEP model

**Status:** Development version - **under active development**

**Documentation:** See detailed tutorial in [docs/tutorials/](../../docs/tutorials/)

### Development active-learning templates: read before running

`workflow_active_learning_dev.sh` and
`workflow_active_learning_dev_multielement.sh` are scheduler-oriented development
templates, not one-command examples. They can create directories, move inputs,
submit jobs, and wait for those jobs to finish. Before using either file, make a
copy outside the repository and review the scheduler directives, executable
paths, input-file requirements, selection method, and all scientific thresholds
for your own system. Do not run the bundled values unchanged as a production
workflow.

---

## Contributing

To add new workflow capabilities, see [CONTRIBUTING.md](../../CONTRIBUTING.md) for detailed guidelines.

---

Thank you for using GPUMDkit! For questions about workflows or to share your workflow adaptations, please open an issue on our [GitHub repository](https://github.com/zhyan0603/GPUMDkit/issues) or contact Zihan YAN (yanzihan@westlake.edu.cn).
