<div align="center">
  <h1>🔧 Workflow Scripts</h1>
    <p style="text-align: justify;">This directory contains preparation tools for high-throughput SCF and MD workflows.</p>
</div>

## Overview

Workflow scripts automate repetitive tasks in computational materials research:

- **SCF batch preprocessing**: Set up VASP/CP2K single-point energy calculations
- **MD sampling**: Prepare molecular dynamics simulations with GPUMD/LAMMPS
- **Active learning support**: Reusable SCF/MD preparation steps for a manually reviewed iteration

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


### Manual active-learning support

For current active-learning work, use entries 301-303 as preparation steps in a
manually reviewed MD -> filtering -> NepTrain FPS -> DFT -> validation cycle.
Review every generated directory, and submit external calculations only after
explicit approval.

The staged workflow and validation gates are documented in the English and
Chinese active-learning tutorials and in the bilingual GPUMDkit skills.

---

## Contributing

To add new workflow capabilities, see [CONTRIBUTING.md](../../CONTRIBUTING.md) for detailed guidelines.

---

Thank you for using GPUMDkit! For questions about workflows or to share your workflow adaptations, please open an issue on our [GitHub repository](https://github.com/zhyan0603/GPUMDkit/issues) or contact Zihan YAN (yanzihan@westlake.edu.cn).
