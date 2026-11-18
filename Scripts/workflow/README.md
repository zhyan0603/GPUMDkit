<div align="center">
  <h1>🔧 Workflow Scripts</h1>
    <p style="text-align: justify;">This directory contains preparation tools for high-throughput SCF and MD workflows.</p>
</div>

## Overview

Workflow scripts automate repetitive tasks in computational materials research:

- **SCF batch preprocessing**: Set up VASP/CP2K/ABACUS single-point energy calculations
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

For extxyz input, structures with one element set keep the existing flat layout
and share `fp/POTCAR`. If multiple element sets are detected, POSCARs are grouped
under `struct_fp/<elements>/`; prepare one `fp/POTCAR_<elements>` file per listed
group. All groups share `fp/INCAR`, and each calculation folder links to its
matching POTCAR. POTCAR element order must match the POSCAR header. Each
calculation folder links to `fp/KPOINTS`; provide that file for explicit k-point
input. If it is absent, the link stays dangling and VASP uses `KSPACING` from
`INCAR`.

For CP2K, enter `3) Workflow`, then `301) SCF batch pretreatment`, and choose the CP2K branch. The CP2K script asks for:

```text
<input.xyz> <template.inp> <prefix>
```

The repository provides `Scripts/workflow/cp2k_template.inp` as a starting template.

For ABACUS, choose the ABACUS branch (`301 -> 3`). The script asks for:

```text
<extxyz_file> <prefix> <basis_dir>
```

For LCAO, `basis_dir` must contain exactly one `.upf` and one `.orb` per element
present in the extxyz, each named with the element symbol first (for
example `Li_ONCV_PBE-1.0.upf`). Each `<prefix>_<index>/` gets its own `STRU`
(only that frame's elements, Angstrom coordinates) and links to the shared
`fp/INPUT` and `fp/KPT`. Write `fp/INPUT` yourself with at least `calculation scf`,
`cal_force 1`, `cal_stress 1` and `basis_type lcao`. Set both `pseudo_dir` and
`orbital_dir` to the basis directory. Convert the finished runs with menu `104`.



### 302) MD sample batch pretreatment (GPUMD)

In a fresh working directory, provide one `.xyz` file (one or more frames) or one or more `.vasp` files. Mixed formats and multiple XYZ files are rejected before any files are moved.

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
