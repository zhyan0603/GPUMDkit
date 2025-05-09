# Function4 - Calculators

This script provides a menu-driven interface to perform various tasks related to property calculations.

### Menu Options

```sh
------------>>
401) Calc ionic conductivity
402) Calc properties by nep
403) Developing ... 
000) Return to the main menu
------------>>
Input the function number:
```

### Option 401: Calc ionic conductivity

This script will read the `msd.out` file and then calculate the diffusivity and ionic conductivity.

Select option `401` from the menu:

```sh
401
```

You will see the following prompt:

```sh
 >-------------------------------------------------<
 | This function calls the script in calculators   |
 | Script: calc_ion_conductivity.py                |
 | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
 >-------------------------------------------------<
 Input <element> <charge> (eg. Li 1)
 ------------>>
```

Enter the `<element>` and `<charge>`:

```sh
Li 1
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

You will see the following output:
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

The script `calc_ion_conductivity.py` in the `Scripts/calculators` will be called to perform the calculations.


### Option 402: Calc properties by nep

This script will calculate energies, forces, and stresses by using the `calorine` package.

Select option `402` from the menu:

```sh
402
```

You will see the following prompt:

```sh
 >-------------------------------------------------<
 | This function calls the script in calculators   |
 | Script: calc_properties_with_nep.py             |
 | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
 >-------------------------------------------------<
 Input <input.xyz> <output.xyz> <nep_model> 
 Examp: input.xyz outpt.xyz nep.txt
 ------------>>
```

Enter the `<input.xyz>` `<output.xyz>` and `<nep_model>`:

```sh
input.xyz output.xyz nep.txt
```

The script `calc_properties_with_nep.py` in the `Scripts/calculators` will be called to perform the calculations.

### Option 403: Calc descriptors of specific elements

This script will calculate the descriptors of specific elements by using the `calorine` package.

Select option `403` from the menu:

```sh
403
```

You will see the following prompt:

```sh
 >-------------------------------------------------<
 | This function calls the script in calculators   |
 | Script: calc_descriptors.py                     |
 | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
 >-------------------------------------------------<
 Input <input.xyz> <output.npy> <nep_model> <element>
 Examp: train.xyz des_Li.npy nep.txt Li
 ------------>>
```

Enter the `<input.xyz>` `<output.npy>`, `<nep_model>` and `<element>`:

```sh
train.xyz des_Li.npy nep.txt Li
```

The script `calc_descriptors.py` in the `Scripts/calculators` will be called to perform the calculations, and output `des_Li.npy` for the analysis.
