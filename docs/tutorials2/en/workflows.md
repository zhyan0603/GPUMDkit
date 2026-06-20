# Workflows Guide

<div align="center">
  <p>
    <a href="../zh/workflows.md">中文</a> | <strong>English</strong>
  </p>
</div>

This guide covers batch processing and workflow automation tools in GPUMDkit.

## Available Workflows

| Workflow | Menu | Description |
|----------|------|-------------|
| SCF Batch (VASP) | 301.1 | VASP single-point batch setup |
| SCF Batch (CP2K) | 301.2 | CP2K single-point batch setup |
| MD Batch (GPUMD) | 302 | GPUMD MD sampling batch setup |
| MD Batch (LAMMPS) | 303 | LAMMPS MD sampling batch setup |

## Interactive Mode

```bash
gpumdkit.sh
# Select: 3) Workflow
```

You'll see:

```
+---------------------------------------------------------+
|                      WORKFLOW TOOLS                     |
+---------------------------------------------------------+
| 301) SCF batch pretreatment                             |
| 302) MD sample batch pretreatment (gpumd)               |
| 303) MD sample batch pretreatment (lmp)                 |
+---------------------------------------------------------+
```

## SCF Batch Pretreatment

### VASP SCF Batch

Prepare batch single-point energy/force calculations for VASP.

```bash
gpumdkit.sh  # Select: 3) Workflow -> 301 -> VASP
```

**Prerequisites:**

1. Place POSCAR files (or single `.xyz`) in current directory
2. Prepare `fp/` directory with:
   - `INCAR`
   - `POTCAR`
   - `KPOINTS`

**Output Structure:**

```
struct_fp/
├── POSCAR_1.vasp
├── POSCAR_2.vasp
└── ...
fp/
├── INCAR
├── POTCAR
└── KPOINTS
<prefix>_1/
├── POSCAR -> ../struct_fp/POSCAR_1.vasp
├── POTCAR -> ../fp/POTCAR
├── INCAR -> ../fp/INCAR
└── KPOINTS -> ../fp/KPOINTS
<prefix>_2/
└── ...
presub.sh
```

### CP2K SCF Batch

Prepare batch single-point calculations for CP2K.

```bash
python Scripts/workflow/scf_batch_pretreatment_cp2k.py <extxyz_file> <template.inp> <prefix_name>

# Example
python Scripts/workflow/scf_batch_pretreatment_cp2k.py structures.xyz cp2k_template.inp H2O_batch
```

**Prerequisites:**

1. Prepare CP2K input template (`cp2k_template.inp`)
2. Template should read coordinates from `pos.xyz`

**Output Structure:**

```
<prefix>_1/
├── input.inp
└── pos.xyz
<prefix>_2/
└── ...
```

## MD Sample Batch Pretreatment

### GPUMD MD Batch

Set up batch MD simulations for GPUMD.

```bash
gpumdkit.sh  # Select: 3) Workflow -> 302
```

**Prerequisites:**

1. Place POSCAR files (or single `.xyz`) in current directory
2. Prepare `md/` directory with:
   - `nep.txt` (NEP model)
   - `run_1.in`, `run_2.in`, ... (GPUMD input files)

**Output Structure:**

```
struct_md/
├── model_1.xyz
├── model_2.xyz
└── ...
md/
├── nep.txt
├── run_1.in
├── run_2.in
└── ...
sample_1/
├── model.xyz -> ../struct_md/model_1.xyz
├── nep.txt -> ../md/nep.txt
└── run.in -> ../md/run_1.in
sample_2/
└── ...
presub.sh
```

### LAMMPS MD Batch

Set up batch MD simulations for LAMMPS.

```bash
gpumdkit.sh  # Select: 3) Workflow -> 303
```

**Prerequisites:**

1. Place POSCAR files (or single `.xyz`) in current directory
2. Prepare `md/` directory with:
   - `nep.txt` (NEP model)
   - `lmprun.in` (LAMMPS input file)

**Output Structure:**

```
struct_md/
├── lammps_1.data
├── lammps_2.data
└── ...
md/
├── nep.txt
└── lmprun.in
sample_1/
├── lammps.data -> ../struct_md/lammps_1.data
├── nep.txt -> ../md/nep.txt
└── lmprun.in -> ../md/lmprun.in
sample_2/
└── ...
presub.sh
```

## Active Learning Workflow

### Overview

Active learning iteratively improves NEP models by:

1. Running MD with current NEP model
2. Selecting diverse structures
3. Computing DFT reference data
4. Retraining NEP model

### Manual Active Learning Loop

```bash
# Step 1: MD sampling with current NEP
gpumdkit.sh  # Select: 3) Workflow -> 302

# Step 2: Filter and sample structures
gpumdkit.sh -min_dist_pbc dump.xyz
gpumdkit.sh -filter_box dump.xyz 13
gpumdkit.sh  # Select: 2) Sample Structures -> 203

# Step 3: Prepare SCF calculations
gpumdkit.sh  # Select: 3) Workflow -> 301

# Step 4: Run DFT calculations (external)

# Step 5: Convert and validate
gpumdkit.sh -out2xyz ./scf_results/
gpumdkit.sh -range new_data.xyz force

# Step 6: Retrain NEP model (external)
```

### Configurable Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `prefix_name` | Directory prefix | `LiF_iter01` |
| `min_dist` | Minimum distance filter | `1.4` |
| `box_limit` | Maximum box size filter | `13` |
| `max_fp_num` | Maximum DFT calculations | `50` |
| `sample_method` | Sampling method | `uniform`, `random`, or `pynep` |
| `pynep_sample_dist` | FPS distance threshold | `0.01` |

## Examples

### Example 1: Temperature-Dependent MD

```bash
# Setup directories for multiple temperatures
for temp in 300 500 700 900; do
    mkdir -p "${temp}K"
    cp model.xyz "${temp}K/"
    cp nep.txt "${temp}K/"
    sed "s/TEMPERATURE/$temp/g" run_template.in > "${temp}K/run.in"
done
```

### Example 2: High-Throughput Screening

```bash
# Prepare structures for multiple calculations
for struct in structures/*.vasp; do
    name=$(basename "$struct" .vasp)
    mkdir -p "calc_$name"
    cp "$struct" "calc_$name/POSCAR"
    ln -s ../POTCAR "calc_$name/"
    ln -s ../INCAR "calc_$name/"
    ln -s ../KPOINTS "calc_$name/"
done
```

### Example 3: NEP Training Data Generation

```bash
# 1. Start with initial structures
gpumdkit.sh -replicate POSCAR supercell.vasp 3 3 3

# 2. Perturb structures
python Scripts/sample_structures/perturb_structure.py supercell.vasp 20 0.03 0.2 uniform

# 3. Run MD sampling
gpumdkit.sh  # Select: 3) Workflow -> 302

# 4. Select diverse structures
gpumdkit.sh  # Select: 2) Sample Structures -> 203

# 5. Prepare DFT calculations
gpumdkit.sh  # Select: 3) Workflow -> 301
```

## Important Notes

1. **Directory structure matters**: Workflows expect specific directory layouts
2. **Symlinks are used**: Output directories use symlinks to save disk space
3. **presub.sh**: Generated submission script may need modification for your cluster
4. **Template files**: Ensure input templates are correct before batch processing
5. **Backup important files**: Always backup your data before running workflows

## See Also

- [Format Conversion](format_conversion.md) - Convert file formats
- [Structure Sampling](sampling.md) - Sample structures for training
- [NEP Training Guide](nep_training.md) - Complete training workflow
