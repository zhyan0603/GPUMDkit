<div align="center">
  <h1>📦 Structure Sampling</h1>
    <p style="text-align: justify;">This directory (Scripts/sample_structures/) contains tools for sampling, selecting, and generating atomic structures for NEP training and molecular dynamics simulations.</p>
</div>

**Script Location:** `Scripts/sample_structures/`

This section covers the structure sampling tools in GPUMDkit (Interactive Mode - Option 2).

Structure sampling is crucial for creating diverse configurations. These scripts provide:

- **Uniform and random sampling**: Basic sampling methods
- **FPS (Farthest Point Sampling)**: Maximizing structural diversity using PyNEP or NEPtrain
- **Structure perturbation**: Generating variations for training set augmentation
- **maximum force deviations selection**: Identifying structures with maximum force deviations
- **Frame extraction**: Selecting specific range from trajectories

## Interactive Mode

```bash
gpumdkit.sh
# Select: 2) Sample Structures
```

You'll see the following menu:

```
 ------------>>
 201) Sample structures from extxyz
 202) Sample structures by pynep
 203) Sample structures by neptrain
 204) Perturb structure
 205) Select max force deviation structs
 000) Return to the main menu
 ------------>>
 Input the function number:
```

### Option 201: Sample Structures from extxyz

This option allows you to sample structures from an `extxyz` file using a specified method.

Select option `201` from the menu:

```
201
```

You will see the following prompt:

```
 >-------------------------------------------------<
 | Calling the script in Scripts/sample_structures |
 | Script: sample_structures.py                    |
 | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
 >-------------------------------------------------<
 Input <extxyz_file> <sampling_method> <num_samples> [skip_num]
 [skip_num]: number of initial frames to skip, default value is 0
 Sampling_method: 'uniform' or 'random'
 Examp: train.xyz uniform 50
 ------------>>
```

Enter the `extxyz` file name, sampling method, and number of samples:

```
train.xyz uniform 50
```

The script `sample_structures.py` in the `Scripts/sample_strcutures` will be called to perform the sampling.

### Option 202: Sample structures by pynep

This function calls the `pynep_select_structs.py` in the `Scripts/sample_structures` to sampling the structures by `pynep`.

Select option `202` from the menu:

```
202
```

You will see the following prompt:

```
 +-------------------------------------------------+
 | PyNEP package is no longer actively maintained  |
 |       Recommend using function 203 instead      |
 +-------------------------------------------------+
 +-------------------------------------------------+
 |     To use parallel version, please use:        |
 |            gpumdkit.sh -pynep                   |
 +-------------------------------------------------+
 >-------------------------------------------------<
 | Calling the script in Scripts/sample_structures |
 | Script: pynep_select_structs.py                 |
 | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
 >-------------------------------------------------<
 Input <sample.xyz> <train.xyz> <nep_model>
 Examp: dump.xyz train.xyz nep.txt
 ------------>>
```

`<samle.xyz>`: extxyz file

`<train.xyz>`: train.xyz

`<nep_model>`: nep.txt

Enter the following parameters:

```
dump.xyz train.xyz nep.txt
```

**Note:** pynep package is no longer actively maintained, this function will be removed soon.

### Option 203: Sample structures by neptrain

This function calls the `neptrain_select_structs.py` script in the `Scripts/sample_structures` to sampling the structures by `neptrain`

Select option `203` from the menu:

```
203
```

You will see the following prompt:

```
 >-------------------------------------------------<
 | Calling the script in Scripts/sample_structures |
 | Script: neptrain_select_structs.py              |
 | Developer: Benrui TANG (tang070205@proton.me)   |
 >-------------------------------------------------<
 Input <sample.xyz> <train.xyz> <nep_model>
 Examp: dump.xyz train.xyz nep.txt
 ------------>>
```

It is the same with option `202`.

### Option 204: Perturb structure

This function calls the `perturb_structure.py` script  in the `Scripts/sample_structures`  to generate the perturbed structures.

Select option `204` from the menu:

```
204
```

You will see the following prompt:

```
>-------------------------------------------------<
| This function calls the script in Scripts       |
| Script: perturb_structure.py                    |
| Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
>-------------------------------------------------<
Input <input.vasp> <pert_num> <cell_pert_fraction> <atom_pert_distance> <atom_pert_style>
The default paramters for perturb are 20 0.03 0.2 uniform
Examp: POSCAR 20 0.03 0.2 normal
------------>>
```

`<input.vasp>`: filename.vasp

`<pert_num>`: number of perturbed structures

`<cell_pert_fraction>`: A fraction determines how much (relatively) will cell deform

`<atom_pert_distance>`: A distance determines how far atoms will move (in angstrom).

`<atom_pert_style>`: `<uniform>`, `<normal>`, `<const>`

Enter your parameters like:

```
POSCAR 20 0.03 0.2 uniform
```

The script `perturb_structure.py` in the `Scripts/sample_strcutures` will be called to perform the perturbation.

### Option 205: Select max force deviation structs from active.xyz

This option allows you to select max force deviation structures from `active.xyz` generated by the command `active` in `gpumd`.

Select option `205` from the menu:

```
205
```

You will see the following prompt:

```
 >-------------------------------------------------<
 | Calling the script in Scripts/sample_structures |
 | Script: select_max_modev.py                     |
 | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
 >-------------------------------------------------<
 +----------------------------------------------------+
 | Select max force deviation structs from active.xyz |
 |     generated by the active command in gpumd.      |
 +----------------------------------------------------+
 Input <structs_num> <threshold> (eg. 200 0.15)
 ------------>>
```

Enter the number of structures and the `threshold` used in your `run.in`:

```
200 0.15
```

The script `select_max_modev.py` in the `Scripts/sample_strcutures` will be called to perform the sampling.

## Contributing

See [CONTRIBUTING.md](contributing.md) for detailed guidelines.

---

Thank you for using GPUMDkit! If you have questions about structure sampling, please open an issue on our [GitHub repository](https://github.com/zhyan0603/GPUMDkit/issues) or contact Zihan YAN (yanzihan@westlake.edu.cn).
