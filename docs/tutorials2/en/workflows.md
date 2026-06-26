<div align="center">
  <h1>Workflows</h1>
  <p>
    <strong>English</strong> | <a href="../zh/workflows.md">简体中文</a>
  </p>
</div>

Batch processing and workflow automation tools.

## Available Workflows

| Workflow | Menu | Description |
|----------|------|-------------|
| SCF Batch (VASP) | 301.1 | VASP single-point batch setup |
| SCF Batch (CP2K) | 301.2 | CP2K single-point batch setup |
| MD Batch (GPUMD) | 302 | GPUMD MD sampling batch setup |
| MD Batch (LAMMPS) | 303 | LAMMPS MD sampling batch setup |

## SCF Batch Pretreatment

### VASP

```bash
gpumdkit.sh  # Select: 3) Workflow -> 301 -> VASP
```

Prerequisites:
1. POSCAR files in current directory
2. `fp/` directory with INCAR, POTCAR, KPOINTS

Output: `<prefix>_1/`, `<prefix>_2/`, ..., `presub.sh`

### CP2K

```bash
python Scripts/workflow/scf_batch_pretreatment_cp2k.py <extxyz> <template.inp> <prefix>
```

## MD Batch Pretreatment

### GPUMD

```bash
gpumdkit.sh  # Select: 3) Workflow -> 302
```

Prerequisites:
1. POSCAR files in current directory
2. `md/` directory with `nep.txt`, `run_*.in`

Output: `sample_1/`, `sample_2/`, ..., `presub.sh`

### LAMMPS

```bash
gpumdkit.sh  # Select: 3) Workflow -> 303
```

Prerequisites:
1. POSCAR files in current directory
2. `md/` directory with `nep.txt`, `lmprun.in`

## Active Learning

Iterative NEP model improvement:

```bash
# 1. MD sampling
gpumdkit.sh  # Select: 3) Workflow -> 302

# 2. Filter and sample
gpumdkit.sh -min_dist_pbc dump.xyz
gpumdkit.sh -filter_box dump.xyz 13

# 3. Prepare SCF
gpumdkit.sh  # Select: 3) Workflow -> 301

# 4. Run DFT (external)

# 5. Retrain NEP
```

## Examples

### Temperature-Dependent MD

```bash
for temp in 300 500 700 900; do
    mkdir -p "${temp}K"
    cp model.xyz "${temp}K/"
    cp nep.txt "${temp}K/"
    sed "s/TEMPERATURE/$temp/g" run_template.in > "${temp}K/run.in"
done
```

### Batch Structure Preparation

```bash
for struct in structures/*.vasp; do
    name=$(basename "$struct" .vasp)
    mkdir -p "calc_$name"
    cp "$struct" "calc_$name/POSCAR"
    ln -s ../POTCAR "calc_$name/"
    ln -s ../INCAR "calc_$name/"
done
```
