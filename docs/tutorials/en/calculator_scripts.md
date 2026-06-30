<div align="center">
  <h1>🧮 Calculator Scripts</h1>
  <p style="text-align: justify;">The calculator module contains scripts for computing properties from MD trajectories, NEP models, and structure files.</p>
</div>

**Script Location:** `Scripts/calculators/`

## Overview

| Task | Command | Main Input |
|------|---------|------------|
| Ionic conductivity | `gpumdkit.sh -calc ionic-cond <element> <charge>` | `msd.out`, `thermo.out`, `model.xyz` |
| MSD from trajectory | `gpumdkit.sh -calc msd <trajectory.xyz> <element> <dt_fs>` | extxyz trajectory |
| NEP prediction | `gpumdkit.sh -calc nep <input.xyz> <output.xyz> <nep.txt>` | extxyz + NEP model |
| NEP descriptors | `gpumdkit.sh -calc des <input.xyz> <output.npy> <nep.txt> <element>` | extxyz + NEP model |
| DOAS | `gpumdkit.sh -calc doas <input.xyz> <nep.txt> <output.txt>` | extxyz + NEP model |
| NEB | `gpumdkit.sh -calc neb <initial.xyz> <final.xyz> <n_images> <nep.txt>` | initial/final structures |
| Minimization | `gpumdkit.sh -calc minimize <structure> <nep.txt> [fmax] [max_steps]` | structure + NEP model |
| Neighbor list | `gpumdkit.sh -calc nlist [args...]` | structure |
| Displacement | `gpumdkit.sh -calc disp [args...]` | trajectory + neighbor list |
| Average structure | `gpumdkit.sh -calc avg-struct [args...]` | trajectory |
| Octahedral tilt | `gpumdkit.sh -calc oct-tilt [args...]` | trajectory + B-O neighbor list |
| ABO3 polarization | `gpumdkit.sh -calc pol-abo3 [args...]` | trajectory + neighbor lists |

For a full command list:

```bash
gpumdkit.sh -calc -h
```

The command-line help table looks like:

```text
+-------------------------------------------------------------------------------------------------------+
|                                      CALCULATOR TOOLS                                                 |
+-------------------------------------------------------------------------------------------------------+
| Usage: gpumdkit.sh -calc <type> [args...]                                                             |
+-------------------------------------------------------------------------------------------------------+
| ionic-cond <element> <charge>                 Calculate ionic conductivity from MSD data              |
| nep <input.xyz> <output.xyz> <nep_model>      Calculate energy/force/virial with a NEP model          |
| des <input.xyz> <output.npy> <nep_model> <el> Calculate NEP descriptors for one element               |
| doas <input.xyz> <nep_model> <output.txt>     Calculate density of atomistic states                   |
| neb <initial.xyz> <final.xyz> <n_images> <nep> Run NEB calculation with a NEP model                   |
| minimize <structure> <nep_model> [fmax] [n]   Minimize a structure with a NEP model                   |
| msd <trajectory.xyz> <element> <dt_fs> [n]    Calculate MSD from an extxyz trajectory                 |
| nlist [script args...]                        Build neighbor lists                                    |
| disp [script args...]                         Calculate displacement from trajectory                  |
| avg-struct [script args...]                   Calculate averaged structure                            |
| oct-tilt [script args...]                     Calculate octahedral tilt                               |
| pol-abo3 [script args...]                     Calculate local polarization for ABO3                   |
+-------------------------------------------------------------------------------------------------------+
```

In interactive mode, choose `4) Calculators`. The menu is:

```text
+----------------------------------------------------------+
|                     CALCULATOR TOOLS                     |
+----------------------------------------------------------+
| 401) Calc ionic conductivity                             |
| 402) Calc properties by nep                              |
| 403) Calc descriptors of specific elements               |
| 404) Calc density of atomistic states (DOAS)             |
| 405) Calc nudged elastic band (NEB) by nep               |
| 406) Build neighbor list                                 |
| 407) Calc displacement from trajectory                   |
| 408) Calc averaged structure                             |
| 409) Calc octahedral tilt                                |
| 410) Calc polarization for ABO3                          |
| 411) Minimize structure by nep                           |
| 412) Calc mean square displacement (MSD) from trajectory |
+----------------------------------------------------------+
| 000) Return to the main menu                             |
+----------------------------------------------------------+
Input the function number:
```

## Ionic Conductivity

`calc_ion_conductivity.py` calculates ionic diffusivity and conductivity from `msd.out`.

### Required and Optional Files

| File | Role |
|------|------|
| `msd.out` | Required MSD data |
| `thermo.out` | Optional, used for automatic temperature detection |
| `model.xyz` | Optional, used for volume and ion-count detection |
| `run.in` | Optional, used to detect replication |

### Usage

```bash
gpumdkit.sh -calc ionic-cond Li 1
```

From interactive mode, choose `401`. You will see:

```text
>-------------------------------------------------<
| This function calls the script in calculators   |
| Script: calc_ion_conductivity.py                |
| Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
>-------------------------------------------------<
Input <element> <charge> (eg. Li 1)
------------>>
```

If automatic files are missing, the script will ask for temperature, volume, and ion count interactively.

Interactive prompts in manual mode look like:

```text
Files 'thermo.out' and 'model.xyz' are not found.
Please provide the following values:
--------------------------->
Enter average temperature (in K):
Enter system volume (in A^3):
Enter number of ions:
```

### Example Output

```text
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

## Mean Square Displacement

`calc_msd.py` computes MSD directly from an extxyz trajectory.

```bash
gpumdkit.sh -calc msd dump.xyz Li 10
gpumdkit.sh -calc msd dump.xyz Li 10 5000
```

From interactive mode, choose `412`. You will see:

```text
Input <extxyz_file> <element_symbol> <dt_fs> [max_corr_steps]
  Optional argument: max_corr_steps (default: frame number)
Example: dump.xyz Li 10
------------>>
```

Arguments:

| Argument | Meaning |
|----------|---------|
| `dump.xyz` | input trajectory |
| `Li` | target mobile species |
| `10` | time interval between frames, in fs |
| `5000` | optional maximum correlation steps |

Output:

- `msd.out`

The beginning of `msd.out` is a text table with time and MSD columns. After generating it, use the plot commands below.

You can then plot:

```bash
gpumdkit.sh -plt msd
gpumdkit.sh -plt sdc
```

<div align="center">
  <img src="../../Gallery/msd.png" alt="MSD plot" width="45%" />
  <img src="../../Gallery/sdc.png" alt="SDC plot" width="45%" />
</div>

## NEP Property Prediction

`calc_properties_with_nep.py` calculates energy, force, and stress for structures using a NEP model.

**Dependency:** `calorine`

```bash
pip install calorine
```

```bash
gpumdkit.sh -calc nep structures.xyz predictions.xyz nep.txt
```

This function is useful when you intentionally want to use a trained NEP model as a surrogate calculator. It should not be treated as a replacement for DFT unless the model quality has been carefully validated.

**Tip:** Before prediction, you may want to clean the extxyz metadata:

```bash
gpumdkit.sh -clean_xyz train.xyz clean_train.xyz
```

## NEP Descriptors

`calc_descriptors.py` extracts NEP descriptors for a selected element.

```bash
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li
```

Use cases:

- visualize chemical environments with PCA/UMAP;
- compare training and candidate structures;
- inspect whether new data expands descriptor space.

Plot descriptors with:

```bash
gpumdkit.sh -plt des pca
gpumdkit.sh -plt des umap
```

<div align="center">
  <img src="../../Gallery/des_umap.png" alt="Descriptor UMAP" width="52%" />
</div>

## Density of Atomistic States

`calc_doas.py` calculates density of atomistic states (DOAS), following the idea proposed by [Wang et al.](https://doi.org/10.1002/anie.202215544).

```bash
gpumdkit.sh -calc doas structures.xyz nep.txt doas.out
gpumdkit.sh -plt doas doas.out Li
```

The script:

1. reads all structures;
2. relaxes each structure with a NEP calculator;
3. extracts per-atom energies;
4. groups atomic energies by element;
5. writes the grouped values to the output file.

For very large systems, doing minimization and atomistic-energy extraction directly in GPUMD may be more efficient.

<div align="center">
  <img src="../../Gallery/doas.png" alt="Density of atomistic states" width="52%" />
</div>

## NEB with a NEP Model

```bash
gpumdkit.sh -calc neb initial.xyz final.xyz 9 nep.txt
```

This runs a NEB calculation with `9` intermediate images. During execution, the script asks how atoms should be fixed:

- `none`: no atoms fixed;
- `index`: fix atoms by index;
- `element`: fix all atoms of one element;
- `position`: fix atoms inside a coordinate range.

## Structure Minimization

`calc_minimize.py` minimizes a structure using a NEP model through `calorine`.

**Dependency:** `calorine`

```bash
pip install calorine
```

```bash
gpumdkit.sh -calc minimize POSCAR nep.txt 0.01 1000
```

Arguments:

| Argument | Meaning |
|----------|---------|
| `POSCAR` | input structure, POSCAR/CONTCAR or extxyz |
| `nep.txt` | NEP model |
| `0.01` | optional force convergence threshold in eV/Ang |
| `1000` | optional maximum optimization steps |

Output:

- `minimize.xyz`
- `minimize.log`

## RDF Calculation with OVITO

`rdf_calculator_ovito.py` calculates radial distribution function using OVITO's analysis tools.

**Dependency:** OVITO

```bash
pip install ovito
```

**Input file:** Structure file (single frame or trajectory)

```bash
python Scripts/calculators/rdf_calculator_ovito.py trajectory.xyz 6.0 400
```

**Parameters:**

| Argument | Meaning |
|----------|---------|
| `trajectory.xyz` | Input structure file |
| `6.0` | Maximum distance for RDF calculation (Å) |
| `400` | Number of histogram bins |

**Visualization:**

```bash
gpumdkit.sh -plt rdf
```

**Note:** It is recommended to use the `compute_rdf` command directly in GPUMD when possible.

---

## Ferroelectric and Polar Material Tools

The following tools are useful for perovskite and polar-material analysis. They require `ferrodispcalc`:

```bash
pip3 install git+https://github.com/MoseyQAQ/ferrodispcalc.git
```

`ferrodispcalc` currently has no associated paper to cite.

### Build Neighbor Lists

```bash
gpumdkit.sh -calc nlist -i model.xyz -c 4.0 -n 6 -C Ti -E O -o nl-Ti-O.dat
```

This creates a neighbor-list file where each row stores one center atom and its nearest neighbor atoms.

### Calculate Displacements

```bash
gpumdkit.sh -calc disp -i movie.xyz -n nl-Ti-O.dat -l 0.2 -o displacements.dat
```

Use `-l 0.2` to analyze the last 20% of frames. Use `-s`, `-t`, and `-p` if you want an explicit slice.

### Average a Structure

```bash
gpumdkit.sh -calc avg-struct -i movie.xyz -l 0.2 -o averaged_structure.xyz
```

This is useful before displacement or polarization analysis when a representative averaged structure is needed.

### Octahedral Tilt

```bash
gpumdkit.sh -calc oct-tilt -i movie.xyz -n nl-Ti-O.dat -l 0.2 -o octahedral_tilt.dat
```

The neighbor list should usually contain six B-O neighbors for each B-site center in an `ABO3`-like structure.

### ABO3 Polarization

```bash
gpumdkit.sh -calc pol-abo3 -i movie.xyz --nl-ba nl-Ti-Pb.dat --nl-bo nl-Ti-O.dat \
  --bec Pb=2 Ti=4.0 O=-2.0 -o polarization.dat
```

The `--bec` values should be chosen according to your model or reference data. The tool is designed for `ABO3`-type local polarization analysis.
