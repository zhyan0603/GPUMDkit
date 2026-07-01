<div align="center">
  <h1>🔧 Workflow Scripts</h1>
  <p style="text-align: justify;">Workflow scripts automate repetitive preparation steps for DFT, MD, and NEP development.</p>
</div>

**Script Location:** `Scripts/workflow/`

## Interactive Entry

```bash
gpumdkit.sh
```

Choose:

```text
3) Workflow
```

The workflow menu looks like:

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

Available entries:

| Menu | Workflow | Purpose |
|------|----------|---------|
| 301 | SCF batch pretreatment | Prepare single-point DFT calculations |
| 302 | MD sample batch pretreatment (GPUMD) | Prepare GPUMD sampling calculations |
| 303 | MD sample batch pretreatment (LAMMPS) | Prepare LAMMPS sampling calculations |

## SCF Batch Pretreatment

Use this when you have many structures and want to prepare single-point DFT calculations.

### VASP

Interactive mode:

```bash
gpumdkit.sh
# choose 3 -> 301 -> VASP
```

Before running, put structure files in the current directory. The script checks files in this order:

1. if one or more `.vasp` files exist, it processes `.vasp` files;
2. if both `.vasp` and `.xyz` files exist, it prints a notice and only processes `.vasp` files;
3. if no `.vasp` file exists and one `.xyz` file exists, it converts that extxyz file to POSCAR files;
4. if no `.vasp` file exists and multiple `.xyz` files exist, it asks you to choose one.

Typical input directory before running:

```text
current directory/
├── POSCAR_1.vasp
├── POSCAR_2.vasp
└── POSCAR_3.vasp
```

or:

```text
current directory/
└── sampled.xyz
```

After running, the script creates:

```text
current directory/
├── struct_fp/
│   ├── POSCAR_1.vasp
│   ├── POSCAR_2.vasp
│   └── ...
├── fp/
├── <prefix>_1/
├── <prefix>_2/
└── presub.sh
```

You must place `INCAR`, `POTCAR`, and `KPOINTS` in the generated `fp/` directory after the script creates it. The calculation folders link to files in `fp/`.

The `prefix` you input during the workflow becomes the calculation folder prefix. When `out2xyz` is later used to extract results, the folder name is stored as `config_type` in the extxyz file, which makes it easier to trace where each structure came from.

### CP2K

CP2K SCF batch pretreatment is available from the workflow menu and can also be called directly:

```bash
gpumdkit.sh
# choose 3 -> 301 -> CP2K
```

```bash
python Scripts/workflow/scf_batch_pretreatment_cp2k.py train.xyz template.inp calc
```

The interactive wrapper asks for:

```text
<extxyz_file> <template.inp> <prefix_name>
```

Example:

```text
dump.xyz template.inp H2O
```

The template should be configured to read coordinates from `pos.xyz`. A reference template is available at:

```text
Scripts/workflow/cp2k_template.inp
```

The script creates `<prefix>_<index>/` directories, writes `input.inp`, and writes `pos.xyz` for each structure.

## MD Sampling Pretreatment

### GPUMD

Interactive mode:

```bash
gpumdkit.sh
# choose 3 -> 302
```

Typical preparation:
Before running, put `.vasp` files or one extxyz file in the current directory. The detection order is the same as the VASP SCF workflow: `.vasp` files are preferred; if multiple `.xyz` files are found without `.vasp`, the script asks you to choose one.

Typical input directory:

```text
current directory/
├── POSCAR_1.vasp
└── POSCAR_2.vasp
```

or:

```text
current directory/
└── dump.xyz
```

After running, the script creates:

```text
current directory/
├── struct_md/
│   ├── model_1.xyz
│   ├── model_2.xyz
│   └── ...
├── md/
├── sample_1/
├── sample_2/
└── presub.sh
```

You must place `nep.txt` and `run_*.in` files in the generated `md/` directory. The script links `run_1.in` to `sample_1/run.in`, `run_2.in` to `sample_2/run.in`, and so on.

### LAMMPS

Interactive mode:

```bash
gpumdkit.sh
# choose 3 -> 303
```

Typical preparation:
Before running, put `.vasp` files or one extxyz file in the current directory. `.vasp` files are preferred when both formats exist.

After running, the script creates LAMMPS data files in `struct_md/`, sample folders, and an `md/` directory:

```text
current directory/
├── struct_md/
│   ├── lammps_1.data
│   ├── lammps_2.data
│   └── ...
├── md/
├── sample_1/
├── sample_2/
└── presub.sh
```

You must place `lmprun.in` and `nep.txt` in the generated `md/` directory.

## Active-Learning Style Workflow

A common NEP improvement loop is:

1. run MD with the current `nep.txt`;
2. collect `dump.xyz` or `active.xyz`;
3. remove unphysical structures;
4. select diverse or uncertain structures;
5. run DFT single-point calculations;
6. convert DFT results to extxyz;
7. append new structures to `train.xyz`;
8. retrain the NEP model.

For the development active-learning script, see [Active Learning](active_learning_workflow.md).

## Template-Based Temperature Series

For a simple temperature series, a shell loop is often enough:

```bash
for temp in 300 500 700 900; do
    mkdir -p "${temp}K"
    cp model.xyz "${temp}K/"
    cp nep.txt "${temp}K/"
    sed "s/TEMPERATURE/$temp/g" run_template.in > "${temp}K/run.in"
done
```

Then submit the generated folders according to your local cluster system.

## Practical Notes

- Workflow scripts are cluster- and template-dependent. Always test with a small number of structures first.
- Keep template files (`INCAR`, `KPOINTS`, `POTCAR`, `run.in`, CP2K inputs) under version control if possible.
- Check generated folders before submitting a large batch.
- Use analyzers and sampling scripts before expensive DFT calculations.
