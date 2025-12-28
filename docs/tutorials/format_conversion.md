<div align="center">
  <h1>🔄 Format Conversion</h1>
    <p style="text-align: justify;">This directory (Scripts/format_conversion/) contains utilities for converting between different file formats commonly used in computational materials science.</p>
</div>

**Script Location:** `Scripts/format_conversion/`

This section covers the format conversion tools in GPUMDkit (Interactive Mode - Option 1).

## Interactive Mode

```bash
gpumdkit.sh
# Select: 1) Format Conversion
```

You'll see the following menu:

```
 ------------>>
 101) Convert VASP to extxyz
 102) Convert mtp to extxyz
 103) Convert CP2K to extxyz
 104) Convert ABACUS to extxyz
 105) Convert extxyz to POSCAR
 000) Return to the main menu
 ------------>>
 Input the function number:
```

for `101`:

```
 >-------------------------------------------------<
 | Calling the script in Scripts/format_conversion |
 | Script: out2xyz.sh                              |
 | Developer: Yanzhou WANG (yanzhowang@gmail.com)  |
 >-------------------------------------------------<
 Choose the type of conversion:
 1) OUTCAR to extxyz
 2) vasprun.xml to extxyz
 ------------>>
```

for `102`:

```
 >-------------------------------------------------<
 | Calling the script in Scripts/format_conversion |
 | Script: mtp2xyz.py                              |
 | Developer: Ke XU (kickhsu@gmail.com)            |
 >-------------------------------------------------<
 Input <filename.cfg> <Symbol1 Symbol2 Symbol3 ...>
 Examp: train.cfg Pd Ag
 ------------>>
```

for `103`:

```
 >-------------------------------------------------<
 | Calling the script in Scripts/format_conversion |
 | Script: cp2k2xyz.py                             |
 | Developer: Ke XU (kickhsu@gmail.com)            |
 >-------------------------------------------------<
 Input [pos.xyz] [frc.xyz] [cell.cell] [-shifted yes/no]
 ------------>>
```

for `104`:

```
 >-------------------------------------------------<
 | Calling the script in Scripts/format_conversion |
 | Script: abacus2xyz-scf.sh/abacus2xyz-md.sh      |
 | Developer: Benrui TANG (tang070205@proton.me)   |
 >-------------------------------------------------<
 Choose the type of ABACUS calculation:
 1) SCF calculation
 2) MD calculation
 ------------>>
```

for `105`:

```
 >-------------------------------------------------<
 | Calling the script in Scripts/format_conversion |
 | Script: exyz2pos.py                             |
 | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
 >-------------------------------------------------<
 Input the name of extxyz
 Examp: ./train.xyz
 ------------>>
```

Follow the prompts to complete the function.

---

## Command-line Mode

| Source Format | Target Format | Command |
|---------------|---------------|---------|
| OUTCAR | extxyz | `gpumdkit.sh -out2xyz <dir>` |
| vasprun.xml | extxyz | `gpumdkit.sh -xml2xyz <dir>` |
| POSCAR | extxyz | `gpumdkit.sh -pos2exyz <poscar> <xyz>` |
| extxyz | POSCAR | `gpumdkit.sh -exyz2pos <xyz>` |
| POSCAR | LAMMPS | `gpumdkit.sh -pos2lmp <poscar> <lmp> <elem...>` |
| LAMMPS dump | extxyz | `gpumdkit.sh -lmp2exyz <dump> <elem...>` |
| CIF | extxyz | `gpumdkit.sh -cif2exyz <cif>` |
| CIF | POSCAR | `gpumdkit.sh -cif2pos <cif>` |
| Add groups | - | `gpumdkit.sh -addgroup <poscar> <elem...>` |
| Add weight | - | `gpumdkit.sh -addweight <in> <out> <weight>` |
| Replicate1 | - | `gpumdkit.sh -replicate input.vasp output.vasp 2 2 2` |
| Replicate2 | - | `gpumdkit.sh -replicate input.vasp output.vasp <target_num>` |
| Get frame | - | `gpumdkit.sh -ge_frame <extxyz> <index>` |

---

## Scripts

### add_groups.py

---

This script adds group labels to structures based on specified elements.

#### Usage

```
python add_groups.py <filename> <Symbols1> <Symbols2> ...
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
python exyz2pos.py <extxyz_file>
```

#### Example

```sh
python exyz2pos.py structs.xyz
```

#### Command-Line Mode Example

```
gpumdkit.sh -exyz2pos structs.xyz
```

This command will convert all frames in `structs.xyz` to `POSCAR_*.vasp` files.



### pos2exyz.py

---

This script converts a `POSCAR` file to `extxyz` format.

#### Usage

```
python pos2exyz.py <POSCAR_file> <extxyz_file>
```

- `<POSCAR_file>`: The path to the input `POSCAR` file.
- `<extxyz_file>`: The desired name for the output `extxyz` file.

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
python pos2lmp.py <poscar_file> <lammps_data_file>
```

- `<poscar_file>`: The path to the input `POSCAR` file.
- `<lammps_data_file>`: The desired name for the output `lammps-data` file.

#### Example

```
python pos2lmp.py POSCAR lammps.data
```

#### Command-Line Mode Example

```
gpumdkit.sh -pos2lmp POSCAR lammps.data
```

This command will read the `POSCAR` file and convert it to `lammps.data` in `lammps-data` format.



### split_single_xyz.py

---

This script splits an `extxyz` file into individual frames, each written to a separate file.

#### Usage

```
python split_single_xyz.py <extxyz_file>
```

This command will split all frames in `extxyz_file` into separate files named `model_*.xyz`.



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
python lmp2exyz.py dump.data Li Y Cl
```

#### Command-Line Mode Example

```
gpumdkit.sh -lmp2exyz dump.data Li Y Cl
```

It will convert the `dump.data` to `dump.xyz` file



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



## Contributing

To add new format converters:

1. **Follow naming**: `<source>2<target>.py`
2. **Handle errors**: Validate input format before processing
3. **Document**: Add usage to this README
6. **Update gpumdkit.sh**: Add command-line flag if appropriate

See [CONTRIBUTING.md](Contributing.md) for detailed guidelines.

---

Thank you for using GPUMDkit! If you have questions about format conversion, please open an issue on our [GitHub repository](https://github.com/zhyan0603/GPUMDkit/issues) or contact Zihan YAN (yanzihan@westlake.edu.cn).
