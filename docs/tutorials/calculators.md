# Function 4 - Calculators

This section covers the calculator tools in GPUMDkit (Interactive Mode - Function 4), which provide computational tools for calculating various material properties using NEP models and analyzing simulation data.

## Interactive Mode Access

```bash
gpumdkit.sh
# Select: 4) Calculators
```

You'll see the following menu:

```
 ------------>>
 401) Calc ionic conductivity
 402) Calc properties by nep
 403) Calc descriptors
 404) Calc DOAS
 405) Calc NEB
 000) Return to the main menu
 ------------>>
```

## Command-Line Usage

For quick operations, you can also use command-line mode:

```bash
# Ionic conductivity
gpumdkit.sh -calc ionic-cond <element> <charge>

# NEP predictions
gpumdkit.sh -calc nep <input.xyz> <output.xyz> <nep.txt>

# Descriptors
gpumdkit.sh -calc des <method> <input.xyz> <output.npy> <nep.txt> <element>

# DOAS calculation
gpumdkit.sh -calc doas <input.xyz> <nep.txt> <output.txt>
```

---

## Available Calculators

### Ionic Conductivity (`calc_ion_conductivity.py`)

**Option 401** in interactive mode.

Calculates ionic diffusivity and conductivity from MSD data.

**Interactive Mode:**

Select option `401` from the calculator menu:

```bash
401
```

You will see the following prompt:

```
 >-------------------------------------------------<
 | This function calls the script in calculators   |
 | Script: calc_ion_conductivity.py                |
 | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
 >-------------------------------------------------<
 Input <element> <charge> (eg. Li 1)
 ------------>>
```

Enter the element and charge:

```bash
Li 1
```

**Command-line:**
```bash
gpumdkit.sh -calc ionic-cond Li 1
```

**Input files:**
- `msd.out` (required)
- `thermo.out` (optional - for temperature)
- `model.xyz` (optional - for volume/atom count)

**Output:**
```
Diffusivity (D):
  D_x: 4.153e-07 cm^2/s
  D_y: 4.174e-07 cm^2/s
  D_z: 2.610e-07 cm^2/s
  D_total: 3.646e-07 cm^2/s
------------------------------
Ionic Conductivity:
  Sigma_x: 2.576e-02 mS/cm
  Sigma_y: 2.589e-02 mS/cm
  Sigma_z: 1.619e-02 mS/cm
  Sigma_total: 2.261e-02 mS/cm
```

### NEP Property Calculation (`calc_properties_with_nep.py`)

Calculates energies, forces, and stresses using a NEP model.

**Command-line:**
```bash
gpumdkit.sh -calc nep input.xyz output.xyz nep.txt
```

**Interactive:** Select option `402`

**Use cases:**
- Generate predictions for test sets
- Calculate properties for new structures
- Validate NEP model performance

### Descriptor Calculation (`calc_descriptors.py`)

Calculates NEP descriptors for structure analysis and visualization.

**Command-line:**
```bash
gpumdkit.sh -calc des umap train.xyz descriptors.npy nep.txt Li
```

**Interactive:** Select option `403`

**Methods:**
- `umap` - UMAP dimensionality reduction
- `tsne` - t-SNE dimensionality reduction
- `pca` - PCA dimensionality reduction

**Workflow:**
```bash
# 1. Calculate descriptors
gpumdkit.sh -calc des umap train.xyz desc.npy nep.txt Li

# 2. Visualize
gpumdkit.sh -plt des umap desc.npy
```

### DOAS Calculation (`calc_doas.py`)

Calculates density of atomistic states by optimizing structures with NEP.

**Command-line:**
```bash
gpumdkit.sh -calc doas structures.xyz nep.txt doas_output.txt
```

**Interactive:** Select option `404`

**Workflow:**
```bash
# 1. Calculate DOAS
gpumdkit.sh -calc doas structures.xyz nep.txt doas.txt

# 2. Visualize
gpumdkit.sh -plt doas doas.txt
```

Reference: [Wang et al.](https://doi.org/10.1002/anie.202215544)

### NEB Calculation (`neb_calculation.py`)

Performs nudged elastic band calculations for transition state finding.

**Direct execution only:**
```bash
python ${GPUMDkit_path}/Scripts/calculators/neb_calculation.py initial.xyz final.xyz 9 nep.txt
```

**Interactive:** Select option `405`

**Output:**
- `neb.traj` - Trajectory file with all images
- `neb_pathway.png` - Energy vs reaction coordinate plot
- Migration barrier in terminal output

### RDF Calculator (`rdf_calculator_ovito.py`)

Calculates radial distribution function using OVITO.

**Direct execution:**
```bash
python ${GPUMDkit_path}/Scripts/calculators/rdf_calculator_ovito.py trajectory.xyz 6.0 400
```

**Parameters:**
- Input XYZ file
- Cutoff distance (Angstroms)
- Number of bins

**Output:** `rdf.txt` file

## Interactive Mode

Access calculators through the interactive menu:

```bash
gpumdkit.sh
# Select: 4) Calculators
# Choose option 401-405
```

### Menu Options

```
 ------------>>
 401) Calc ionic conductivity
 402) Calc properties by nep
 403) Calc descriptors
 404) Calc DOAS
 405) Calc NEB
 000) Return to the main menu
 ------------>>
```

## Common Workflows

### Analyze Ionic Transport
```bash
# 1. Run GPUMD with compute_msd
# 2. Calculate conductivity
gpumdkit.sh -calc ionic-cond Li 1

# 3. Plot diffusion
gpumdkit.sh -plt msd
gpumdkit.sh -plt sdc
```

### NEP Model Analysis
```bash
# 1. Test NEP on structures
gpumdkit.sh -calc nep test.xyz predictions.xyz nep.txt

# 2. Calculate descriptors
gpumdkit.sh -calc des umap train.xyz desc.npy nep.txt Li

# 3. Visualize
gpumdkit.sh -plt prediction
gpumdkit.sh -plt des umap desc.npy
```

### Find Migration Barrier
```bash
# 1. Prepare initial and final structures
# 2. Run NEB
python ${GPUMDkit_path}/Scripts/calculators/neb_calculation.py start.xyz end.xyz 9 nep.txt

# 3. Examine neb_pathway.png for barrier height
```

## Tips

- **Ionic conductivity**: Requires sufficient MSD sampling time
- **Descriptors**: Use for training set coverage analysis
- **DOAS**: Structures will be optimized before calculation
- **NEB**: Choose appropriate number of images (7-11 typical)
- **RDF**: Useful for offline analysis of trajectories

---

For more details, see [Scripts/calculators/README.md](../../Scripts/calculators/README.md)
