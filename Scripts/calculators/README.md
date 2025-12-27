# Calculator Scripts

This directory contains computational tools for calculating various material properties using NEP models and analyzing simulation data.

## Overview

The calculator scripts provide functionality for:
- Ionic conductivity and diffusivity calculations
- Property predictions using NEP models
- Descriptor calculations for machine learning analysis
- Density of atomistic states (DOAS) analysis
- Nudged elastic band (NEB) calculations
- Radial distribution function (RDF) calculations

All scripts can be accessed through `gpumdkit.sh` using the `-calc` flag or run directly.

---

## Scripts

### calc_ion_conductivity.py

Calculates ionic diffusivity and conductivity from mean square displacement (MSD) data.

#### Purpose
Determines transport properties of ionic species using the Nernst-Einstein equation, connecting MSD to diffusion coefficients and ionic conductivity.

#### Input Files
- `msd.out` - MSD data from GPUMD (required)
- `thermo.out` - Temperature data (optional, for automatic temperature detection)
- `model.xyz` - Structure file (optional, for automatic volume/atom count detection)

#### Usage

**Command-line mode:**
```bash
gpumdkit.sh -calc ionic-cond <element> <charge>
```

**Direct execution:**
```bash
python calc_ion_conductivity.py <element> <charge>
```

#### Parameters
- `<element>`: Chemical species symbol (e.g., Li, Na, O)
- `<charge>`: Formal charge of the ion (e.g., 1 for Li⁺, -2 for O²⁻)

#### Examples

**Automatic mode** (with thermo.out and model.xyz present):
```bash
gpumdkit.sh -calc ionic-cond Li 1
```

**Manual input mode** (if files not found):
```bash
gpumdkit.sh -calc ionic-cond Li 1
# Script will prompt:
Files 'thermo.out' and 'model.xyz' are not found.
Please provide the following values:
--------------------------->
Enter average temperature (in K): 800
Enter system volume (in Å^3): 16785
Enter number of ions: 448
```

#### Output

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

#### Theory
- **Diffusivity**: Calculated from MSD slope using Einstein relation: D = MSD/(6t)
- **Conductivity**: Nernst-Einstein equation: σ = (n·q²·D)/(kB·T)
  - n: number density of ions
  - q: ion charge
  - kB: Boltzmann constant
  - T: temperature

#### Requirements
- `msd.out` must contain at least 4 columns: time, MSD_x, MSD_y, MSD_z
- Sufficient sampling time for linear MSD regime
- Proper equilibration before MSD calculation

---

### calc_properties_with_nep.py

Calculates energies, forces, and stresses for structures using a NEP model.

#### Purpose
Performs predictions on structure datasets using trained NEP models, useful for:
- Generating predictions for test sets
- Creating training data from DFT-relaxed structures
- Validating NEP model performance
- Large-scale property screening

#### Input Files
- Structure file in extxyz format
- NEP model file (nep.txt)

#### Usage

**Command-line mode:**
```bash
gpumdkit.sh -calc nep <input.xyz> <output.xyz> <nep.txt>
```

**Direct execution:**
```bash
python calc_properties_with_nep.py <input.xyz> <output.xyz> <nep.txt>
```

#### Parameters
- `<input.xyz>`: Input structure file (extxyz format)
- `<output.xyz>`: Output file with NEP predictions
- `<nep.txt>`: NEP model file

#### Example

```bash
# Predict properties for a test set
gpumdkit.sh -calc nep test_structures.xyz predictions.xyz nep.txt

# Generate NEP properties for new structures
python calc_properties_with_nep.py new_configs.xyz nep_results.xyz nep.txt
```

#### Output Format
Output extxyz file contains:
- Original atomic positions
- NEP-predicted energies (per-atom and total)
- NEP-predicted forces (per atom)
- NEP-predicted stresses/virials (if applicable)
- Original lattice parameters

#### Use Cases
1. **Model validation**: Predict on test set and compare with DFT
2. **High-throughput screening**: Calculate properties for many structures
3. **Trajectory analysis**: Compute energies/forces for MD trajectories
4. **Active learning**: Generate predictions for structure selection

#### Requirements
- Python package: `calorine` (NEP interface)
- ASE (Atomic Simulation Environment)

```bash
pip install calorine ase
```

---

### calc_descriptors.py

Calculates NEP descriptors for structures, which can be used for dimensionality reduction and structure analysis.

#### Purpose
Extracts high-dimensional descriptors from NEP model for:
- Visualizing training set coverage (UMAP/t-SNE)
- Identifying similar structures
- Active learning structure selection
- Analyzing descriptor space distribution

#### Input Files
- Structure file (extxyz format)
- NEP model file (nep.txt)

#### Usage

**Command-line mode:**
```bash
gpumdkit.sh -calc des <method> <input.xyz> <output.npy> <element>
```

**Direct execution:**
```bash
python calc_descriptors.py <input.xyz> <output.npy> <nep.txt> <element>
```

#### Parameters
- `<method>`: Dimensionality reduction method (umap, tsne, pca)
- `<input.xyz>`: Input structure file
- `<output.npy>`: Output NumPy file with descriptors
- `<nep.txt>`: NEP model file
- `<element>`: Target element for descriptor extraction

#### Examples

```bash
# Calculate descriptors using UMAP
gpumdkit.sh -calc des umap train.xyz descriptors.npy nep.txt Li

# Calculate for multiple elements
gpumdkit.sh -calc des umap train.xyz desc_li.npy nep.txt Li
gpumdkit.sh -calc des umap train.xyz desc_o.npy nep.txt O

# Visualize after calculation
gpumdkit.sh -plt des umap descriptors.npy save
```

#### Output
- NumPy array file (.npy) containing descriptor vectors
- Can be visualized using `gpumdkit.sh -plt des`

#### Applications
- **Training set analysis**: Check if training data covers configuration space
- **Active learning**: Select diverse structures based on descriptor distance
- **Model comparison**: Compare descriptor distributions across different models
- **Outlier detection**: Identify unusual configurations

#### Reference
See [arXiv:2504.15925](https://doi.org/10.48550/arXiv.2504.15925) for methodology and applications.

#### Requirements
```bash
pip install calorine numpy scikit-learn
```

---

### calc_doas.py

Calculates the density of atomistic states (DOAS) by optimizing structures and computing per-atom energies.

#### Purpose
DOAS provides atomic-level electronic structure information, proposed by [Wang et al.](https://doi.org/10.1002/anie.202215544). Useful for:
- Understanding local atomic environments
- Analyzing electronic structure at atomic resolution
- Comparing different element contributions
- Identifying active sites in catalysts

#### Input Files
- Structure file (extxyz format)
- NEP model file (nep.txt)

#### Usage

**Command-line mode:**
```bash
gpumdkit.sh -calc doas <input.xyz> <nep.txt> <output.txt>
```

**Direct execution:**
```bash
python calc_doas.py <input.xyz> <nep.txt> <output.txt>
```

#### Parameters
- `<input.xyz>`: Input structure file
- `<nep.txt>`: NEP model file for structure optimization and energy calculation
- `<output.txt>`: Output file with grouped atomic energies by element

#### Example

```bash
# Calculate DOAS
gpumdkit.sh -calc doas structures.xyz nep.txt doas_results.txt

# Visualize results
gpumdkit.sh -plt doas doas_results.txt save
```

#### Process
1. Reads all structures from input file
2. Optimizes each structure using NEP (BFGS, fmax=0.05 eV/Å)
3. Extracts per-atom energies from optimized structures
4. Groups energies by element type
5. Outputs distribution for each element

#### Output Format
Text file with per-atom energies grouped by element:
```
Li: -5.23 -5.21 -5.25 ...
Y: -8.45 -8.43 -8.47 ...
Cl: -3.12 -3.15 -3.10 ...
```

#### Visualization
After calculation, visualize with:
```bash
gpumdkit.sh -plt doas element1.txt element2.txt
```

#### Requirements
```bash
pip install calorine ase tqdm
```

---

### neb_calculation.py

Performs nudged elastic band (NEB) calculations to find minimum energy pathways between initial and final states.

#### Purpose
NEB method identifies transition states and reaction pathways for:
- Ion migration barriers
- Phase transformations
- Chemical reactions
- Diffusion mechanisms

#### Input Files
- Initial structure (XYZ format)
- Final structure (XYZ format)
- NEP model file (nep.txt)

#### Usage

**Direct execution only:**
```bash
python neb_calculation.py <initial.xyz> <final.xyz> <n_images> <nep.txt>
```

#### Parameters
- `<initial.xyz>`: Initial structure file
- `<final.xyz>`: Final structure file
- `<n_images>`: Number of intermediate images (typically 5-15)
- `<nep.txt>`: NEP model for energy/force calculations

#### Example

```bash
# Calculate Li migration pathway with 9 intermediate images
python neb_calculation.py li_start.xyz li_end.xyz 9 nep.txt

# Interactive prompts will ask about atom constraints
```

#### Interactive Options
The script will prompt for:
1. **Constraint method**:
   - None: All atoms free to move
   - Index: Fix specific atom indices
   - Element: Fix all atoms of specific element
   - Layer: Fix bottom layer(s)

2. **Constraint details** based on chosen method

#### Output Files
- `neb.traj`: ASE trajectory file with all NEB images
- `neb_pathway.png`: Plot of energy vs reaction coordinate
- Terminal output with barrier height

#### Example Output
```
Created 11 images (1 initial + 9 intermediate + 1 final)
Select fixing method:
1) None (all atoms free)
2) Fix by indices
3) Fix by element
4) Fix bottom layer
Choice: 4
Bottom layer Z threshold: 5.0
NEB optimization converged
Forward barrier: 0.45 eV
Reverse barrier: 0.32 eV
```

#### Visualization
The script automatically generates:
- Energy vs reaction coordinate plot
- Shows initial, transition state, and final energies
- Indicates forward and reverse barriers

#### Best Practices
1. **Initial/final alignment**: Ensure structures are well-aligned
2. **Number of images**: 7-11 images usually sufficient for simple migrations
3. **Constraints**: Fix bottom layers for surface processes
4. **Convergence**: Check forces converge to <0.05 eV/Å
5. **Validation**: Verify path makes chemical sense

#### Requirements
```bash
pip install calorine ase matplotlib numpy
```

#### Authors
- Zhoulin LIU (1776627910@qq.com) - Original implementation
- Zihan YAN (yanzihan@westlake.edu.cn) - Modifications

---

### rdf_calculator_ovito.py

Calculates radial distribution function (RDF) using OVITO's analysis tools.

#### Purpose
Computes pair correlation functions for analyzing:
- Local atomic structure
- Nearest neighbor distances
- Coordination numbers
- Structural ordering

#### Input Files
- Structure file (extxyz format)

#### Usage

**Direct execution:**
```bash
python rdf_calculator_ovito.py <input.xyz> <cutoff> <bins>
```

#### Parameters
- `<input.xyz>`: Input structure file (can be single frame or trajectory)
- `<cutoff>`: Maximum distance for RDF calculation (Angstroms)
- `<bins>`: Number of histogram bins for RDF

#### Example

```bash
# Calculate RDF with 6 Å cutoff and 400 bins
python rdf_calculator_ovito.py trajectory.xyz 6.0 400

# Higher resolution
python rdf_calculator_ovito.py model.xyz 8.0 800
```

#### Output
- `rdf.txt`: Text file with RDF data
  - Column 1: Distance (Å)
  - Columns 2+: g(r) for each pair type

#### Advantages over GPUMD RDF
- **Offline analysis**: Calculate RDF from existing trajectories
- **Flexible binning**: Choose any number of bins
- **Post-processing**: Easy to compute RDF from saved structures
- **Multiple frames**: Average over trajectory frames

#### Use Cases
1. **Structure validation**: Verify structure quality before simulations
2. **Trajectory analysis**: Compute RDF from MD trajectories
3. **Phase identification**: Compare RDF with reference patterns
4. **Custom analysis**: Extract specific pair correlations

#### Visualization
After calculation:
```bash
gpumdkit.sh -plt rdf      # If saved as rdf.out
# Or plot manually with your preferred tool
```

#### Requirements
```bash
pip install ovito ase
```

**Note**: OVITO package can be large (~500MB). For simple RDF, consider using GPUMD's built-in `compute_rdf` command.

---

## General Usage Guidelines

### Command Syntax

Access calculators through `gpumdkit.sh`:
```bash
gpumdkit.sh -calc <calculator> [arguments]
```

Or run scripts directly:
```bash
python <script_name>.py [arguments]
```

### Quick Reference Table

| Calculator | Command | Primary Input | Output |
|------------|---------|---------------|--------|
| Ionic conductivity | `ionic-cond` | `msd.out` | Diffusivity, conductivity |
| NEP properties | `nep` | structures.xyz | predictions.xyz |
| Descriptors | `des` | structures.xyz | descriptors.npy |
| DOAS | `doas` | structures.xyz | doas.txt |
| NEB | Direct only | initial/final.xyz | pathway, barrier |
| RDF | Direct only | structures.xyz | rdf.txt |

### Common Workflows

#### Workflow 1: Analyze Ionic Conductivity

```bash
# 1. Run GPUMD with compute_msd
# 2. Calculate conductivity
gpumdkit.sh -calc ionic-cond Li 1

# 3. Plot diffusion properties
gpumdkit.sh -plt msd save
gpumdkit.sh -plt sdc save
```

#### Workflow 2: NEP Model Analysis

```bash
# 1. Calculate properties with NEP
gpumdkit.sh -calc nep test.xyz predictions.xyz nep.txt

# 2. Calculate descriptors
gpumdkit.sh -calc des umap train.xyz descriptors.npy nep.txt Li

# 3. Visualize
gpumdkit.sh -plt prediction save
gpumdkit.sh -plt des umap descriptors.npy save
```

#### Workflow 3: Migration Barrier Calculation

```bash
# 1. Prepare initial and final structures
# 2. Run NEB
python neb_calculation.py initial.xyz final.xyz 9 nep.txt

# 3. Analyze pathway
# (Automatically generates neb_pathway.png)
```

## Dependencies

### Core Requirements
- **Python 3.x**
- **NumPy**: Numerical computations
- **ASE**: Atomic structure handling

### Calculator-Specific
- **calorine**: NEP model interface (for NEP calculations)
- **OVITO**: RDF calculations (optional)
- **matplotlib**: NEB pathway visualization
- **tqdm**: Progress bars
- **scikit-learn**: Descriptor analysis

### Installation

**Basic installation:**
```bash
pip install numpy ase
```

**Full installation:**
```bash
pip install numpy ase calorine matplotlib tqdm scikit-learn
```

**Optional OVITO:**
```bash
pip install ovito
```

## Best Practices

1. **Verify inputs**: Check structure files are properly formatted before calculations
2. **Use appropriate models**: Ensure NEP model is trained for your system
3. **Check convergence**: Verify calculations converged properly
4. **Save outputs**: Always save important calculation results
5. **Validate results**: Compare with reference data when available

## Troubleshooting

### Issue: "calorine not found"

**Problem**: NEP calculations require calorine package

**Solution**:
```bash
pip install calorine
```

### Issue: "Descriptor calculation fails"

**Problem**: Element not in NEP model

**Solution**: Verify element exists in model and use correct symbol

### Issue: "MSD file format error"

**Problem**: `msd.out` not in expected format

**Solution**: Ensure MSD computed with correct GPUMD settings:
```
compute_msd 100 1 0
```

### Issue: "NEB doesn't converge"

**Problem**: Poor initial/final structure alignment or too few images

**Solution**:
1. Check structure alignment
2. Increase number of images
3. Relax initial/final structures first
4. Check constraints are reasonable

## Contributing

To add new calculator scripts:

1. **Follow naming**: `calc_<descriptive_name>.py`
2. **Add documentation**: Include docstring with usage
3. **Handle errors**: Validate inputs and provide helpful error messages
4. **Update README**: Add documentation to this file
5. **Test thoroughly**: Verify with multiple systems
6. **Consider integration**: Add command-line flag to `gpumdkit.sh` if appropriate

See [CONTRIBUTING.md](../../CONTRIBUTING.md) for detailed guidelines.

---

Thank you for using GPUMDkit! If you have questions or need assistance with calculator scripts, please open an issue on our [GitHub repository](https://github.com/zhyan0603/GPUMDkit/issues) or contact Zihan YAN (yanzihan@westlake.edu.cn).
