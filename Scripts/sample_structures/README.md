#### sample_structures.py

---

This script samples structures from an `extxyz` file using either '`uniform`' or '`random`' sampling methods. The sampled structures are then written to the '`sampled_structures.xyz`' file.

`<extxyz_file>`: The path to the input extxyz file containing structures.

`<sampling_method>`: The sampling method to use. Can be '`uniform`' or '`random`'.

`<num_samples>`: The number of samples to extract from the input file.

#### Usage

```sh
python sample_structures.py <extxyz_file> <sampling_method> <num_samples>
```

#### Example

```sh
python sample_structures.py train.xyz uniform 10
```

This command will sample 10 structures uniformly from `train.xyz` and save them to `sampled_structures.xyz`.



#### get_min_dist.py

---

This script calculates the minimum atomic distance within a given system. The input is an `extxtz` file, and the script outputs the minimum distance between any two atoms in the structure.

#### Usage

```bash
python get_min_dist.py <extxyz_file>
```

