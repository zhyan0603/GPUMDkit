# Format Conversion Scripts

This directory contains utilities for converting between different atomic structure file formats commonly used in computational materials science.

## Overview

The format conversion scripts provide seamless interconversion between:
- **VASP** (POSCAR, OUTCAR) ↔ extxyz
- **CP2K** output → extxyz
- **ABACUS** output → extxyz  
- **LAMMPS** dump → extxyz
- **CIF** → POSCAR/extxyz
- Adding metadata (group labels, weights)
- Frame extraction and manipulation

All scripts can be accessed through `gpumdkit.sh` using various flags or run directly.

---

## Quick Command Reference

| Source Format | Target Format | Command |
|---------------|---------------|---------|
| OUTCAR | extxyz | `gpumdkit.sh -out2xyz <dir>` |
| POSCAR | extxyz | `gpumdkit.sh -pos2exyz <poscar> <xyz>` |
| extxyz | POSCAR | `gpumdkit.sh -exyz2pos <xyz>` |
| POSCAR | LAMMPS | `gpumdkit.sh -pos2lmp <poscar> <lmp> <elem...>` |
| LAMMPS dump | extxyz | `gpumdkit.sh -lmp2exyz <dump> <elem...>` |
| CIF | extxyz | `gpumdkit.sh -cif2exyz <cif>` |
| CIF | POSCAR | `gpumdkit.sh -cif2pos <cif>` |
| Add groups | - | `gpumdkit.sh -addgroup <poscar> <elem...>` |
| Add weight | - | `gpumdkit.sh -addweight <in> <out> <weight>` |

---

## Scripts

### add_groups.py

---

This script adds group labels to structures based on specified elements.

#### Usage

```
python add_groups.py <filename> <Symbols>
```

- `<filename>`: The path to the input file (e.g., POSCAR).
- `<Symbols>`: Space-separated list of element symbols to group (e.g., Li Y Cl).

#### Example

```sh
python add_groups.py POSCAR Li Y Cl
```

#### Command-Line Mode Example

```
gpumdkit.sh -addgroup POSCAR Li Y Cl
```

This command will read the `POSCAR` file and add group labels for the elements `Li`, `Y`, and `Cl`. The output will be saved to a file named `model.xyz`.



### add_weight.py

---

This script adds weight labels to structures.

#### Usage

```
python add_weight.py <input_file> <output_file> <new_weight>
```

- `<inputfile>`: The path to the input file (e.g., train.xyz).
- `<outputfile>`: The path to the input file (e.g., train_weighted.xyz).
- `<new_weight>`: The `weight` you need to change.

#### Example

```sh
python add_weight.py train.xyz train_weighted.xyz 5
```

#### Command-Line Mode Example

```
gpumdkit.sh -addweight train.xyz train_weighted.xyz 5
```

This command will read the `train.xyz` file and add `Weight=5` labels for all structures. The output will be saved to a file named `train_weighted.xyz`.



### exyz2pos.py

---

This script converts all frames in an `extxyz` file to `POSCAR` format.

#### Usage

```sh
python exyz2pos.py [extxyz_filename]
```

- `[extxyz_filename]`: (Optional) The path to the input `extxyz` file. If not specified, the default is `train.xyz`.

#### Example

```sh
python exyz2pos.py my_structures.xyz
```

#### Command-Line Mode Example

```
gpumdkit.sh -exyz2pos my_structures.xyz
```

This command will convert all frames in `my_structures.xyz` to `POSCAR` files.



### pos2exyz.py

---

This script converts a `POSCAR` file to `extxyz` format.

#### Usage

```
python pos2exyz.py <POSCAR_filename> <extxyz_filename>
```

- `<POSCAR_filename>`: The path to the input `POSCAR` file.
- `<extxyz_filename>`: The desired name for the output `extxyz` file.

#### Example

```
python pos2exyz.py POSCAR model.xyz
```

#### Command-Line Mode Example

```
gpumdkit.sh -pos2exyz POSCAR model.xyz
```

This command will read the `POSCAR` file and convert it to `model.xyz` in `extxyz` format.



### pos2lmp.py

---

This script converts a `POSCAR` file to `lammps-data` format.

#### Usage

```
python pos2lmp.py <poscar_file> <lammps_data_file> <elements_order>
```

- `<poscar_file>`: The path to the input `POSCAR` file.
- `<lammps_data_file>`: The desired name for the output `lammps-data` file.
- `<elements_order>`: The desired order of elements in `lammps-data` file.

#### Example

```
python pos2lmp.py POSCAR lammps.data Li La Zr O
```

#### Command-Line Mode Example

```
gpumdkit.sh -pos2lmp POSCAR lammps.data Li La Zr O
```

This command will read the `POSCAR` file and convert it to `lammps.data` in `lammps-data` format.



### split_single_xyz.py

---

This script splits an `extxyz` file into individual frames, each written to a separate file.

#### Usage

```
python split_single_xyz.py <extxyz_filename>
```

- `<extxyz_filename>`: The path to the input `extxyz` file.

#### Example

```sh
python split_single_xyz.py train.extxyz
```

This command will split all frames in `train.extxyz` into separate files named `model_${i}.xyz`, where `${i}` is the frame index.



### lmp2exyz.py

---

This script will convert the `lammps-dump` to `extxyz` format.

#### Usage

```
python lmp2exyz.py <dump_file> <element1> <element2> ...
```

- `<dump_file>`: The path to the input `lammps-dump` file.
- `<element>`: The order of the specified elements.

#### Example

```sh
python lmp2exyz.py dump.lammps Li Y Cl
```

#### Command-Line Mode Example

```
gpumdkit.sh -lmp2exyz dump.lammps Li Y Cl
```



### get_frame.py

---

This script will read the `extxyz` file and return the specified frame by index..

#### Usage

```
python get_frame.py <extxyz_file> <frame_index>
```

- `<extxyz_file>`: The path to the input `extxyz` file.
- `<frame_index>`: The index of the specified frame.

#### Example

```sh
python get_frame.py dump.xyz 1000
```

#### Command-Line Mode Example

```
gpumdkit.sh -get_frame dump.xyz 1000
```

You will get the `frame_1000.xyz` file after perform the script.



---

## General Usage Guidelines

### Command Syntax

Access format conversion tools through `gpumdkit.sh`:
```bash
gpumdkit.sh -<conversion_flag> [arguments]
```

Or run scripts directly:
```bash
python <script_name>.py [arguments]
```

### Common Workflows

#### Workflow 1: VASP → NEP Training

```bash
# 1. Convert VASP outputs to extxyz
gpumdkit.sh -out2xyz ./VASP_calculations/

# 2. Add group labels for NEP
gpumdkit.sh -addgroup POSCAR Li Y Cl

# 3. (Optional) Add weights for specific structures
gpumdkit.sh -addweight train.xyz train_weighted.xyz 2.0
```

#### Workflow 2: Structure Preparation

```bash
# 1. Convert CIF to POSCAR
gpumdkit.sh -cif2pos structure.cif

# 2. Convert to extxyz with groups
gpumdkit.sh -pos2exyz POSCAR model.xyz
gpumdkit.sh -addgroup POSCAR Li La Zr O

# 3. Use in GPUMD
# (Now have model.xyz ready for simulation)
```

#### Workflow 3: Multi-Format Pipeline

```bash
# 1. Start with CIF from database
gpumdkit.sh -cif2exyz structure.cif

# 2. Extract single frame
gpumdkit.sh -get_frame structure.xyz 0

# 3. Convert to LAMMPS format
gpumdkit.sh -exyz2pos frame_0.xyz
gpumdkit.sh -pos2lmp POSCAR lammps.data Li La Zr O
```

#### Workflow 4: Batch OUTCAR Processing

```bash
# Directory structure:
# calculations/
# ├── struct_001/OUTCAR
# ├── struct_002/OUTCAR
# └── struct_003/OUTCAR

# Convert all OUTCARs
gpumdkit.sh -out2xyz ./calculations/

# Output: combined train.xyz with all structures
```

## Format Specifications

### extxyz Format

Extended XYZ format used by GPUMD and NEP:

```
<number_of_atoms>
Lattice="a1 a2 a3 b1 b2 b3 c1 c2 c3" Properties=species:S:1:pos:R:3:forces:R:3 energy=<total_energy> virial="v1 v2 v3 v4 v5 v6" pbc="T T T"
<element> <x> <y> <z> <fx> <fy> <fz>
...
```

**Key fields:**
- `Lattice`: 3x3 lattice vectors (row-major)
- `Properties`: Atomic properties definition
- `energy`: Total energy (eV)
- `virial`: Stress tensor (eV)
- `pbc`: Periodic boundary conditions

### Group Labels

Groups specify atom types for NEP training:

```
<number_of_atoms>
Lattice="..." Properties=species:S:1:pos:R:3 group=<group_string>
Li <x> <y> <z>
Y <x> <y> <z>
Cl <x> <y> <z>
```

**Group string**: Space-separated integers, one per atom
- Example: `group=0 0 0 1 1 2 2 2` (3 Li as type 0, 2 Y as type 1, 3 Cl as type 2)

### Weight Labels

Weights adjust structure importance in NEP training:

```
<number_of_atoms>
Lattice="..." Properties=... Weight=<weight_value>
```

**Typical weights:**
- `Weight=1.0` - Normal importance (default)
- `Weight=2.0` - Double importance
- `Weight=0.5` - Half importance

## File Format Details

### VASP POSCAR

```
System name
1.0
a1x a1y a1z
a2x a2y a2z
a3x a3y a3z
Li Y Cl
4 2 8
Direct
fractional_coords...
```

### VASP OUTCAR

Large file with detailed calculation output. Converters extract:
- Final or all ionic positions
- Energies
- Forces
- Stress tensor

### LAMMPS Dump

```
ITEM: TIMESTEP
<step>
ITEM: NUMBER OF ATOMS
<N>
ITEM: BOX BOUNDS
<xlo> <xhi>
<ylo> <yhi>
<zlo> <zhi>
ITEM: ATOMS id type x y z fx fy fz
<atom_data>...
```

### CIF (Crystallographic Information File)

Standard crystallographic format. Converters handle:
- Symmetry operations
- Fractional coordinates
- Space groups
- Unit cell parameters

## Advanced Usage

### Batch Conversion with Shell Loops

```bash
# Convert multiple POSCAR files
for poscar in POSCAR_*; do
    output="${poscar//POSCAR/model}.xyz"
    gpumdkit.sh -pos2exyz $poscar $output
done

# Merge multiple xyz files
cat model_*.xyz > combined.xyz
```

### Conditional Weights

```python
# Custom script for conditional weighting
import ase.io as io

atoms_list = io.read('train.xyz', ':')
for i, atoms in enumerate(atoms_list):
    if atoms.info.get('energy_per_atom', 0) > -5.0:
        atoms.info['Weight'] = 2.0  # Higher weight for high-energy configs
    else:
        atoms.info['Weight'] = 1.0
io.write('train_weighted.xyz', atoms_list)
```

### Element Order Consistency

For NEP training, maintain consistent element ordering:

```bash
# Define element order
ELEMENTS="Li Y Cl"

# Use same order in all conversions
gpumdkit.sh -addgroup POSCAR $ELEMENTS
gpumdkit.sh -pos2lmp POSCAR lmp.data $ELEMENTS
```

## Dependencies

### Required
- **Python 3.x**
- **ASE** (Atomic Simulation Environment)
- **NumPy**

### Optional
- **pymatgen**: Enhanced CIF handling
- **MDAnalysis**: Additional format support

### Installation

```bash
pip install ase numpy
# Optional:
pip install pymatgen
```

## Best Practices

1. **Verify conversions**: Always check first few structures after conversion
2. **Consistent units**: Ensure energy (eV), distance (Å), force (eV/Å)
3. **Check lattice**: Verify lattice vectors after conversion
4. **Element mapping**: Keep element order consistent across files
5. **Backup originals**: Keep original files before conversion

## Troubleshooting

### Issue: "File format not recognized"

**Problem**: Input file has unexpected format

**Solution**:
1. Check file has correct extension
2. Verify file isn't corrupted
3. Try opening in text editor to inspect format

### Issue: "Energy/forces not found in OUTCAR"

**Problem**: VASP calculation didn't complete or write output

**Solution**:
1. Check VASP completed successfully
2. Verify OUTCAR isn't empty
3. Check INCAR has appropriate output settings

### Issue: "Element order mismatch"

**Problem**: Converted file has different element order than expected

**Solution**:
1. Explicitly specify element order in command
2. Check POSCAR comment line for element symbols
3. Manually reorder if needed

### Issue: "Group labels incorrect"

**Problem**: Wrong number of group values or incorrect mapping

**Solution**:
1. Count atoms per element: should match group counts
2. Verify element list matches POSCAR
3. Check for extra spaces in element specification

### Issue: "Lattice vectors seem wrong"

**Problem**: Unusual cell shape or negative volumes after conversion

**Solution**:
1. Check original file lattice definition
2. Verify conversion preserves cell parameters
3. Visualize structure with VESTA or ASE to confirm

## Performance Tips

1. **Batch operations**: Process multiple files in parallel when possible
2. **Large trajectories**: Use frame slicing for huge OUTCAR files
3. **Memory management**: Process in chunks if RAM limited
4. **Disk I/O**: Use SSD for faster file operations

## Integration with Other Tools

### With Visualization

```bash
# Convert and visualize
gpumdkit.sh -pos2exyz POSCAR model.xyz
ase gui model.xyz  # Visualize with ASE
```

### With Analysis

```bash
# Convert and analyze
gpumdkit.sh -out2xyz ./
gpumdkit.sh -analyze_comp train.xyz
gpumdkit.sh -range train.xyz force
```

### With Training

```bash
# Prepare training set
gpumdkit.sh -out2xyz ./DFT_calcs/
gpumdkit.sh -addgroup POSCAR Li Y Cl
# Now ready for NEP training with train.xyz
```

## File Organization

Recommended directory structure:

```
project/
├── raw_data/
│   ├── struct_001/
│   │   ├── POSCAR
│   │   └── OUTCAR
│   └── struct_002/
│       ├── POSCAR
│       └── OUTCAR
├── converted/
│   ├── train.xyz
│   └── test.xyz
└── models/
    └── model.xyz
```

## Contributing

To add new format converters:

1. **Follow naming**: `<source>2<target>.py`
2. **Handle errors**: Validate input format before processing
3. **Preserve metadata**: Keep all relevant information
4. **Test thoroughly**: Verify with various structure types
5. **Document**: Add usage to this README
6. **Update gpumdkit.sh**: Add command-line flag if appropriate

See [CONTRIBUTING.md](../../CONTRIBUTING.md) for detailed guidelines.

---

Thank you for using GPUMDkit! If you have questions about format conversion, please open an issue on our [GitHub repository](https://github.com/zhyan0603/GPUMDkit/issues) or contact Zihan YAN (yanzihan@westlake.edu.cn).
