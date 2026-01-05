<div align="center">
  <h1>🧮 Calculators</h1>
    <p style="text-align: justify;">This directory (Scripts/calculators/) contains computational tools for calculating various material properties using NEP models and analyzing simulation data.</p>
</div>

**Script Location:** `Scripts/calculators/`

This section covers the calculator tools in GPUMDkit (Interactive Mode - Option 4).

The calculator scripts provide functionality for:

- Ionic conductivity and diffusivity calculations
- Property predictions using NEP models
- Descriptor calculations for PCA or UMAP analysis
- Density of atomistic states (DOAS)
- Nudged elastic band (NEB) calculations
- Radial distribution function (RDF) calculations

## Interactive Mode

```bash
gpumdkit.sh
# Select: 4) Calculators
```

You'll see the following menu:

```
 ------------>>
 401) Calc ionic conductivity
 402) Calc properties by nep
 403) Calc descriptors of specific elements
 404) Calc density of atomistic states (DOAS)
 405) Calc nudged elastic band (NEB) by nep
 000) Return to the main menu
 ------------>>
 Input the function number:
```

for `401`:

```
 >-------------------------------------------------<
 | This function calls the script in calculators   |
 | Script: calc_ion_conductivity.py                |
 | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
 >-------------------------------------------------<
 Input <element> <charge> (eg. Li 1)
 ------------>>
```

for `402`:

```
 >-------------------------------------------------<
 | This function calls the script in calculators   |
 | Script: calc_properties_with_nep.py             |
 | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
 >-------------------------------------------------<
 Input <input.xyz> <output.xyz> <nep_model>
 Examp: input.xyz output.xyz nep.txt
 ------------>>
```

for `403`:

```
 >-------------------------------------------------<
 | This function calls the script in calculators   |
 | Script: calc_descriptors.py                     |
 | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
 >-------------------------------------------------<
 Input <input.xyz> <output.npy> <nep_model> <element>
 Examp: train.xyz des_Li.npy nep.txt Li
 ------------>>
```

for `404`:

```
 >-------------------------------------------------<
 | This function calls the script in calculators   |
 | Script: calc_doas.py                            |
 | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
 >-------------------------------------------------<
 Input <input.xyz> <nep_model> <output_file>
 Examp: dump.xyz nep.txt doas.out
 ------------>>
```

for `405`:

```
 >-------------------------------------------------<
 | This function calls the script in calculators   |
 | Script: neb_calculation.py                      |
 | Developer: Zhoulin LIU (1776627910@qq.com)      |
 >-------------------------------------------------<
 Input <initial_struct> <final_struct> <n_image> <nep_model>
 Examp: IS.xyz FS.xyz 5 nep.txt
 ------------>>
```

Follow the prompts to complete the function.

## Command-Line Mode

```bash
gpumdkit.sh -calc ionic-cond <element> <charge>
gpumdkit.sh -calc nep <input.xyz> <output.xyz> <nep.txt>
gpumdkit.sh -calc des <input.xyz> <output.npy> <nep.txt> <element>
gpumdkit.sh -calc doas <input.xyz> <nep.txt> <output.txt>
```

---

## Scripts

### calc_ion_conductivity.py

Calculates ionic diffusivity and conductivity from mean square displacement (MSD) data.

#### Input Files
- `msd.out` - MSD data from GPUMD (required)
- `thermo.out` - Temperature data (optional, for automatic temperature detection)
- `model.xyz` - Structure file (optional, for automatic volume/atom count detection)
- `run.in` - MD input file (optional, for replicate detection)

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

---

### calc_properties_with_nep.py

Calculates energies, forces, and stresses for structures using a NEP model. You can provide a clean extxyz file and use a `nep.txt` to label the reference values for energies, forces, and stresses. 

***Note that you must know what you are doing and this is not a substitute for DFT calculations.***

Unless you have a well-trained NEP model and you want to use it as a DFT surrogate, you probably shouldn't use this feature. See [this article]([ https://doi.org/10.48550/arXiv.2504.15925](https://doi.org/10.48550/arXiv.2504.15925)) if you are interesting in it.

#### Input Files
- A clean structure file in `extxyz` format with energy, force, and virial/stress. 

  (use `gpumdkit.sh -clean_xyz train.xyz` to remove the original values)

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

#### Example

```bash
# Predict properties for a test set
gpumdkit.sh -calc nep structures.xyz predictions.xyz nep.txt
```

#### Requirements
- Python package: `calorine` 

```bash
pip install calorine
```

---

### calc_descriptors.py

Calculates descriptors for the specific species, which can be used for dimensionality reduction and structure analysis.

#### Usage

**Command-line mode:**
```bash
gpumdkit.sh -calc des <input.xyz> <output.npy> <nep.txt> <element>
```

**Direct execution:**
```bash
python calc_descriptors.py <input.xyz> <output.npy> <nep.txt> <element>
```

#### Parameters
- `<input.xyz>`: Input structure file
- `<output.npy>`: Output `npy` file with high-dimensional descriptors
- `<nep.txt>`: NEP model file
- `<element>`: Target element for descriptor calculations

#### Examples

```bash
# Calculate descriptors using UMAP
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li
```

#### Output
- NumPy array file (.npy) containing descriptor vectors
- Can be visualized using `gpumdkit.sh -plt des`

#### Applications
- **Training set analysis**: Check if training data covers configuration space
- **Model comparison**: Compare descriptor distributions across different models

#### Reference
See [arXiv:2504.15925](https://doi.org/10.48550/arXiv.2504.15925) if you are interesting in it.

---

### calc_doas.py

Calculates the [density of atomistic states (DOAS)](https://doi.org/10.1002/anie.202215544) by optimizing structures and computing per-atom energies.

#### Usage

**Command-line mode:**
```bash
gpumdkit.sh -calc doas <input.xyz> <nep.txt> <doas.out>
```

**Direct execution:**
```bash
python calc_doas.py <input.xyz> <nep.txt> <doas.out>
```

#### Parameters
- `<input.xyz>`: Input structure file
- `<nep.txt>`: NEP model file for structure optimization and atomistic-energies calculation
- `<output.txt>`: Output file with grouped atomic energies by element

#### Example

```bash
# Calculate DOAS
gpumdkit.sh -calc doas structures.xyz nep.txt doas.out

# Visualize results
gpumdkit.sh -plt doas doas.txt <species>
gpumdkit.sh -plt doas doas.txt Li
```

#### Process
1. Reads all structures from input file
2. Optimizes each structure using NEP (BFGS, fmax=0.05 eV/Å)
3. Extracts per-atom energies from optimized structures
4. Groups energies by element type
5. Outputs distribution for each element

You can also do the `minimize` and atomistic energy calculations in `gpumd`, which is prefered for large-scale structures.

---

### neb_calculation.py

Performs nudged elastic band (NEB) calculations to find minimum energy pathways between initial and final states.

#### Usage

**Direct execution only:**
```bash
python neb_calculation.py <initial.xyz> <final.xyz> <n_images> <nep.txt>
```

#### Parameters
- `<initial.xyz>`: Initial structure file
- `<final.xyz>`: Final structure file
- `<n_images>`: Number of intermediate images
- `<nep.txt>`: NEP model for energy/force calculations

#### Example

```bash
# Calculate Li migration pathway with 9 intermediate images
python neb_calculation.py init.xyz fin.xyz 9 nep.txt
```

#### Interactive Options
```
 Input the function number:
 4
 ------------>>
 401) Calc ionic conductivity
 402) Calc properties by nep
 403) Calc descriptors of specific elements
 404) Calc density of atomistic states (DOAS)
 405) Calc nudged elastic band (NEB) by nep
 000) Return to the main menu
 ------------>>
 Input the function number:
 405
 >-------------------------------------------------<
 | This function calls the script in calculators   |
 | Script: neb_calculation.py                      |
 | Developer: Zhoulin LIU (1776627910@qq.com)      |
 >-------------------------------------------------<
 Input <initial_struct> <final_struct> <n_image> <nep_model>
 Examp: IS.xyz FS.xyz 5 nep.txt
 ------------>>
```

---

### rdf_calculator_ovito.py

Calculates radial distribution function (RDF) using OVITO's analysis tools.

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
```

#### Visualization
After calculation:
```bash
gpumdkit.sh -plt rdf    
```

**Note**: It is recommended to use the `compute_rdf` command in `gpumd`.

---

## Contributing

See [CONTRIBUTING.md](contributing.md) for detailed guidelines.

---

Thank you for using GPUMDkit! If you have questions or need assistance with calculator scripts, please open an issue on our [GitHub repository](https://github.com/zhyan0603/GPUMDkit/issues) or contact Zihan YAN (yanzihan@westlake.edu.cn).
