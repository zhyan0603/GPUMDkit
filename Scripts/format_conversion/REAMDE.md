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

This command will read the `POSCAR` file and add group labels for the elements `Li`, `Y`, and `Cl`. The output will be saved to a file named `model.xyz`.



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

This command will read the `POSCAR` file and convert it to `model.xyz` in `extxyz` format.



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

