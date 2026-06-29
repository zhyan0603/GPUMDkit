---
name: gpumdkit-workflows
description: >
  Use when setting up batch calculations or high-throughput simulations.
  Provides SCF batch pretreatment for VASP and CP2K, MD sample batch pretreatment for GPUMD and LAMMPS,
  and active learning workflow automation.
  Use when user asks about: batch calculation, workflow automation, SCF batch, MD sampling,
  high-throughput simulation, or active learning.
allowed-tools: Bash(gpumdkit *) Bash(python3 *)
---

# GPUMDkit Workflows

## Available Workflows

| Workflow | Menu | Description |
|----------|------|-------------|
| SCF Batch (VASP) | 301.1 | VASP single-point batch setup |
| SCF Batch (CP2K) | 301.2 | CP2K single-point batch setup |
| MD Batch (GPUMD) | 302 | GPUMD MD sampling batch setup |
| MD Batch (LAMMPS) | 303 | LAMMPS MD sampling batch setup |

## SCF Batch Pretreatment

### VASP SCF Batch

```bash
# Interactive mode
gpumdkit.sh  # Select: 3) Workflow -> 301 -> VASP
```

**Prerequisites:**
1. Place POSCAR files (or single `.xyz`) in current directory
2. Prepare `fp/` directory with:
   - `INCAR`
   - `POTCAR`
   - `KPOINTS`

**Output structure:**
```
struct_fp/POSCAR_1.vasp, POSCAR_2.vasp, ...
fp/POTCAR, INCAR, KPOINTS          (user-provided)
<prefix>_1/POSCAR -> struct_fp/POSCAR_1.vasp
<prefix>_1/POTCAR -> fp/POTCAR
<prefix>_2/...
presub.sh
```

### CP2K SCF Batch

```bash
# Direct Python execution
python Scripts/workflow/scf_batch_pretreatment_cp2k.py <extxyz_file> <template.inp> <prefix_name>

# Example
python Scripts/workflow/scf_batch_pretreatment_cp2k.py structures.xyz cp2k_template.inp H2O_batch
```

**Prerequisites:**
1. Prepare CP2K input template (`cp2k_template.inp`)
2. Template should read coordinates from `pos.xyz`

**Output structure:**
```
<prefix>_1/input.inp, pos.xyz
<prefix>_2/input.inp, pos.xyz
...
```

## MD Sample Batch Pretreatment

### GPUMD MD Batch

```bash
# Interactive mode
gpumdkit.sh  # Select: 3) Workflow -> 302
```

**Prerequisites:**
1. Place POSCAR files (or single `.xyz`) in current directory
2. Prepare `md/` directory with:
   - `nep.txt` (NEP model)
   - `run_1.in`, `run_2.in`, ... (GPUMD input files)

**Output structure:**
```
struct_md/model_1.xyz, model_2.xyz, ...
md/nep.txt, run_1.in, run_2.in, ...   (user-provided)
sample_1/model.xyz -> struct_md/model_1.xyz
sample_1/nep.txt -> md/nep.txt
sample_1/run.in -> md/run_1.in
presub.sh
```

### LAMMPS MD Batch

```bash
# Interactive mode
gpumdkit.sh  # Select: 3) Workflow -> 303
```

**Prerequisites:**
1. Place POSCAR files (or single `.xyz`) in current directory
2. Prepare `md/` directory with:
   - `nep.txt` (NEP model)
   - `lmprun.in` (LAMMPS input file)

**Output structure:**
```
struct_md/lammps_1.data, lammps_2.data, ...
md/lmprun.in, nep.txt                  (user-provided)
sample_1/lammps.data -> struct_md/lammps_1.data
sample_1/lmprun.in -> md/lmprun.in
sample_1/nep.txt -> md/nep.txt
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

### Automated Active Learning (Development)

```bash
# Scripts available but under development:
# Scripts/workflow/workflow_active_learning_dev.sh
# Scripts/workflow/workflow_active_learning_dev_multielement.sh
```

**Configurable parameters:**

| Parameter | Description | Example |
|-----------|-------------|---------|
| `prefix_name` | Directory prefix | `LiF_iter01` |
| `min_dist` | Minimum distance filter | `1.4` |
| `box_limit` | Maximum box size filter | `13` |
| `max_fp_num` | Maximum DFT calculations | `50` |
| `sample_method` | Sampling method | `uniform`, `random`, or `pynep` |
| `pynep_sample_dist` | FPS distance threshold | `0.01` |

## Examples

### Example 1: VASP High-Throughput Screening
```bash
# Prepare structures
for struct in structures/*.vasp; do
    name=$(basename "$struct" .vasp)
    mkdir -p "calc_$name"
    cp "$struct" "calc_$name/POSCAR"
    ln -s ../POTCAR "calc_$name/"
    ln -s ../INCAR "calc_$name/"
    ln -s ../KPOINTS "calc_$name/"
done
```

### Example 2: Temperature-Dependent MD
```bash
# Setup directories for multiple temperatures
for temp in 300 500 700 900; do
    mkdir -p "${temp}K"
    cp model.xyz "${temp}K/"
    cp nep.txt "${temp}K/"
    sed "s/TEMPERATURE/$temp/g" run_template.in > "${temp}K/run.in"
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

## Notes

1. **Directory structure matters**: Workflows expect specific directory layouts
2. **Symlinks are used**: Output directories use symlinks to save disk space
3. **presub.sh**: Generated submission script may need modification for your cluster
4. **Template files**: Ensure input templates are correct before batch processing

## Detailed Documentation

- [Workflow Guide](../../docs/tutorials/en/workflow_scripts.md)
- [Active Learning Guide](../../docs/tutorials/en/active_learning_workflow.md)
