### calc_ion_conductivity.py

---

This script will read the `msd.out` file and then calculate the diffusivity and ionic conductivity.

#### Usage

```
python calc_ion_conductivity.py <element> <charge>
```

- `<element>`: Chemical species (e.g., Li).
- `<charge>`: Ion charge (e.g., 1 for Li⁺).

#### Example

```sh
python calc_ion_conductivity.py Li 1
```

If required files (`thermo.out` and `model.xyz`) are not found, the script will prompt for manual input of structure `volume`, `number of ions`, and `temperature`.

```
Files 'thermo.out' and 'model.xyz' are not found.
Please provide the following values:
--------------------------->
Enter average temperature (in K): 800
Enter system volume (in Å^3): 16785
Enter number of ions: 448
```

#### Command-Line Mode Example

```
gpumdkit.sh -calc ion-cond Li 1
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

### calc_properties_with_nep.py

---

This script will calculate energies, forces, and stresses by using the `calorine` package.

#### Usage

```
python calc_properties_with_nep.py <input.xyz> <output.xyz> <nep_model>
```

#### Example

```sh
python calc_properties_with_nep.py input.xyz output.xyz nep.txt
```

### rdf_calculator_ovito.py

---

This script will read the `dump.xyz` file and then calculate the RDF by using the `ovito` package.

#### Usage

```
python rdf_calculator_ovito.py <extxyz_file> <cutoff> <bins>
```

- `<extxyz_file>`: The path to the `extxyz` file.
- `<cutoff>`: The cutoff used in the RDF calculation.
- `<bins>`: The bins used in the RDF calculation.

#### Example

```sh
python rdf_calculator_ovito.py dump.xyz 6 400
```

This command will read the `dump.xyz` file and calculate the RDF by using the `ovito` package and output the `rdf.txt` file.



---

Thank you for using `GPUMDkit`! If you have any questions or need further assistance, feel free to open an issue on our GitHub repository or contact Zihan YAN (yanzihan@westlake.edu.cn).
