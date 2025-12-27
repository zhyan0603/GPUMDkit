# Format Conversion

GPUMDkit provides comprehensive tools for converting between various file formats used in computational materials science. This guide covers both interactive and command-line modes for format conversion. 

### Menu Options

```sh
 ------------>>
 101) Convert OUTCAR to extxyz
 102) Convert mtp to extxyz
 103) Convert cp2k to extxyz
 104) Convert castep to extxyz
 105) Convert extxyz to POSCAR
 106) Developing ...
 000) Return to the main menu
 ------------>>
 Input the function number:
```

### Option 101: Convert OUTCAR to extxyz

This option allows you to convert VASP `OUTCAR` files to `extxyz` format.

Select option `101` from the menu:

```
101
```

You will see the following prompt:

```
 >-------------------------------------------------<
 | This function calls the script in GPUMD's tools |
 | Script: multipleFrames-outcars2nep-exyz.sh      |
 | Developer: Yanzhou WANG (yanzhowang@gmail.com ) |
 >-------------------------------------------------<
 Input the directory containing OUTCARs
 ------------>>
```

Enter the directory containing your `OUTCAR` files:

```
/path/to/your/outcars
```

The script `multipleFrames-outcars2nep-exyz.sh` in GPUMD's tools will be called to perform the conversion.

### Option 102: Convert mtp to extxyz

This option allows you to convert `cfg` files to `extxyz` format.

Select option `102` from the menu:

```
102
```

You will see the following prompt:

```
 >-------------------------------------------------<
 | This function calls the script in GPUMD's tools |
 | Script: mtp2xyz.py                              |
 | Developer: Ke XU (kickhsu@gmail.com)            |
 >-------------------------------------------------<
 Input <filename.cfg> <Symbol1 Symbol2 Symbol3 ...>
 Examp: train.cfg Pd Ag
 ------------>>
```

Enter the `<filename.cfg>` `<Symbol1 Symbol2 Symbol3 ...>` :

```
train.cfg Pd Ag
```

The script `mtp2xyz.py` in GPUMD's tools will be called to perform the conversion.

### Option 103: Convert cp2k to extxyz

This option allows you to convert cp2k's output to `extxyz` format.

Select option `103` from the menu:

```
103
```

You will see the following prompt:

```
 >-------------------------------------------------<
 | This function calls the script in GPUMD's tools |
 | Script: cp2k2xyz.py                             |
 | Developer: Ke XU (kickhsu@gmail.com)            |
 >-------------------------------------------------<
 Input <dir_cp2k>
 Examp: ./cp2k
 ------------>>
```

Enter the `<dir_cp2k>` :

```
./cp2k
```

The script `cp2k2xyz.py` in GPUMD's tools will be called to perform the conversion.

### Option 104: Convert castep to extxyz

This option allows you to convert castep's output to `extxyz` format.

Select option `104` from the menu:

```
104
```

You will see the following prompt:

```
 >-------------------------------------------------<
 | This function calls the script in GPUMD's tools |
 | Script: castep2nep-exyz.sh                      |
 | Developer: Yanzhou WANG (yanzhowang@gmail.com ) |
 >-------------------------------------------------<
 Input <dir_castep>
 Examp: ./castep
 ------------>>
```

Enter the `<dir_castep>` :

```
./castep
```

The script `castep2nep-exyz.sh` in GPUMD's tools will be called to perform the conversion.

### Option 105: Convert extxyz to POSCAR

This option allows you to convert `extxyz` file to POSCAR.

Select option `105` from the menu:

```
105
```

You will see the following prompt:

```
 >-------------------------------------------------<
 | This function calls the script in Scripts       |
 | Script: exyz2pos.py                             |
 | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
 >-------------------------------------------------<
 Input the name of extxyz
 Examp: ./train.xyz 
 ------------>>
```

Enter the `<extxyz_filename>` :

```
./train.xyz
```

The script `exyz2pos.py` in Scripts will be called to perform the conversion.



---

## Command-Line Usage

For quick operations, use command-line mode:

```bash
# Convert VASP OUTCAR to extxyz
gpumdkit.sh -out2xyz ./path/to/outcars/

# Convert POSCAR to extxyz
gpumdkit.sh -pos2exyz POSCAR model.xyz

# Convert CIF to extxyz
gpumdkit.sh -cif2exyz structure.cif

# Convert CIF to POSCAR
gpumdkit.sh -cif2pos structure.cif

# Convert POSCAR to LAMMPS data
gpumdkit.sh -pos2lmp POSCAR lammps.data Li Y Cl

# Convert LAMMPS dump to extxyz
gpumdkit.sh -lmp2exyz dump.lammpstrj Li Y Cl

# Add group labels for NEP training
gpumdkit.sh -addgroup POSCAR Li Y Cl

# Add weight to structures
gpumdkit.sh -addweight train.xyz train_weighted.xyz 2.0

# Extract specific frame
gpumdkit.sh -get_frame dump.xyz 100
```

## Supported Formats

| Source Format | Target Format | Interactive | Command-Line | Script |
|--------------|---------------|-------------|--------------|---------|
| VASP OUTCAR | extxyz | 101 | `-out2xyz` | `out2xyz.sh` |
| VASP vasprun.xml | extxyz | 101 | `-out2xyz` | `xml2xyz.py` |
| MTP .cfg | extxyz | 102 | - | `mtp2xyz.py` |
| CP2K output | extxyz | 103 | `-cp2k2exyz` | `cp2k2xyz.py` |
| ABACUS output | extxyz | 104 | `-abacus2xyz` | `abacus2xyz_*.sh` |
| extxyz | POSCAR | 105 | `-exyz2pos` | `exyz2pos.py` |
| POSCAR | extxyz | - | `-pos2exyz` | `pos2exyz.py` |
| POSCAR | LAMMPS | - | `-pos2lmp` | `pos2lmp.py` |
| LAMMPS dump | extxyz | - | `-lmp2exyz` | - |
| CIF | extxyz | - | `-cif2exyz` | `cif2exyz.py` |
| CIF | POSCAR | - | `-cif2pos` | - |

## Common Workflows

### Prepare NEP Training Data from VASP

```bash
# 1. Convert all OUTCAR files in subdirectories
gpumdkit.sh -out2xyz ./dft_calculations/

# 2. Output: train.xyz with all structures combined

# 3. Add group labels (for multi-element systems)
gpumdkit.sh -addgroup POSCAR Li Y Cl
# Output: model.xyz with group information

# 4. (Optional) Add weights to specific structures
gpumdkit.sh -addweight train.xyz train_weighted.xyz 2.0
```

### Convert CIF Database to GPUMD Format

```bash
# Convert multiple CIF files
for cif in *.cif; do
    gpumdkit.sh -cif2exyz "$cif"
done

# Combine into single file if needed
cat *.xyz > combined_structures.xyz
```

### Prepare LAMMPS Simulation from VASP

```bash
# 1. Convert POSCAR to LAMMPS data file
gpumdkit.sh -pos2lmp POSCAR lammps.data Li La Zr O

# 2. Use in LAMMPS simulation
# Note: Specify element order matching your LAMMPS potential
```

## File Format Details

### Extended XYZ (extxyz)

The extended XYZ format used by GPUMD includes:

```
<number_of_atoms>
Lattice="a1 a2 a3 b1 b2 b3 c1 c2 c3" Properties=species:S:1:pos:R:3:forces:R:3 energy=<value> virial="v1 v2 ... v6"
<element> <x> <y> <z> <fx> <fy> <fz>
...
```

**Key fields:**
- `Lattice`: 3x3 lattice vectors (row-major, in Angstroms)
- `energy`: Total energy (eV)
- `forces`: Forces on atoms (eV/Angstrom)
- `virial`: Stress tensor components (eV)

### Group Labels

For NEP training with multiple species, add group labels:

```
<number_of_atoms>
Lattice="..." Properties=species:S:1:pos:R:3 group=<group_string>
Li <x> <y> <z>
Y <x> <y> <z>
Cl <x> <y> <z>
```

The `group` string assigns atom types, e.g., `group=0 0 0 1 1 2 2 2` means:
- First 3 atoms (Li) are type 0
- Next 2 atoms (Y) are type 1  
- Last 3 atoms (Cl) are type 2

## Tips

- **Batch conversion**: Use shell loops for multiple files
- **Check results**: Always verify a few structures after conversion
- **Element order**: Keep consistent element ordering across conversions
- **Group labels**: Required for NEP training with multiple species
- **Weights**: Use to emphasize important structures in training

---

Thank you for using `GPUMDkit`! If you have any questions or need further assistance, feel free to open an issue on our GitHub repository or contact Zihan YAN (yanzihan@westlake.edu.cn).