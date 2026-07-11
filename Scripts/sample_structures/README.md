<div align="center">
  <h1>📦 Structure Sampling Scripts</h1>
    <p style="text-align: justify;">This directory contains tools for sampling, selecting, and generating atomic structures for NEP training and molecular dynamics simulations.</p>
</div>

## Overview

Structure sampling is crucial for creating diverse configurations. These scripts provide:
- **Uniform and random sampling**: Basic sampling methods
- **FPS (Farthest Point Sampling)**: Maximizing structural diversity using PyNEP or NEPtrain
- **Structure perturbation**: Generating variations for training set augmentation
- **maximum force deviations selection**: Identifying structures with maximum force deviations
- **Training/test splitting**: Creating train and test extxyz files using uniform, random, or FPS selection
- **Frame extraction**: Selecting specific range from trajectories

---

## Via interactive mode

Run `gpumdkit.sh` and choose `2) Sample Structures`. The current submenu is:

```text
 +------------------------------------------------------+
 |                 SAMPLE STRUCTURE TOOLS               |
 +------------------------------------------------------+
 | 201) Sample structures from extxyz                   |
 | 202) PyNEP sampling [deprecated]                     |
 | 203) FPS sampling by NepTrain [preferred]            |
 | 204) Perturb structure                               |
 | 205) Select max force deviation structs              |
 | 206) Split training and test sets                    |
 +------------------------------------------------------+
 | 000) Return to the main menu                         |
 +------------------------------------------------------+
 Input the function number:
```

---

## Scripts

### sample_structures.py

This script samples frames from an extxyz trajectory using either uniform or random sampling.

#### Usage

```sh
python sample_structures.py <extxyz_file> <sampling_method> <num_samples> [skip_initial]
```

- `<sampling_method>`: `uniform` or `random`
- `[skip_initial]`: optional number of initial frames to skip; default is `0`

#### Example

```sh
python sample_structures.py train.xyz uniform 50
```

#### Interactive Mode Example

```
201) Sample structures from extxyz
Input <extxyz_file> <sampling_method> <num_samples> [skip_num]
[skip_num]: number of initial frames to skip, default value is 0
Sampling_method: 'uniform' or 'random'
Example: train.xyz uniform 50
```

The output file is `sampled_structures.xyz`.



### pynep_select_structs.py / parallel_pynep_select_structs.py

PyNEP-based FPS sampling is kept for compatibility, but it is only executed through `gpumdkit.sh -pynep`. The interactive `202` entry only prints a notice and does not run PyNEP. The interactive sampling menu uses `203) FPS sampling by NepTrain` for descriptor-based FPS.

#### Direct Usage

```sh
python pynep_select_structs.py <sampledata_file> <traindata_file> <nep_model_file>
python parallel_pynep_select_structs.py <sampledata_file> <traindata_file> <nep_model_file> <threads>
```

#### GPUMDkit Entry

```sh
gpumdkit.sh -pynep
```

Then input:

```sh
dump.xyz train.xyz nep.txt 8
```



### neptrain_select_structs.py

This script selects diverse structures from candidate structures using NepTrain descriptors and farthest-point sampling.

#### Usage

```sh
python neptrain_select_structs.py <sampledata_file> <traindata_file> <nep_model_file>
```

#### Example

```sh
python neptrain_select_structs.py dump.xyz train.xyz nep.txt
```

#### Interactive Mode Example

```
203) FPS sampling by NepTrain [preferred]
Input <sample.xyz> <train.xyz> <nep_model>
Example: dump.xyz train.xyz nep.txt
```

This function requires the `NepTrain` package. If you use this function, we recommend citing the NepTrain paper.



### perturb_structure.py

This script generates perturbed structures from an input VASP structure.

#### Usage

```sh
python perturb_structure.py <input.vasp> <pert_num> <cell_pert_fraction> <atom_pert_distance> <atom_pert_style>
```

#### Example

```sh
python perturb_structure.py POSCAR 20 0.03 0.2 uniform
```

#### Interactive Mode Example

```
204) Perturb structure
Input <input.vasp> <pert_num> <cell_pert_fraction> <atom_pert_distance> <atom_pert_style>
The default parameters for perturb are 20 0.03 0.2 uniform
Example: POSCAR 20 0.03 0.2 uniform
```

This function requires the `dpdata` package. If you use this function, we recommend citing dpdata.



### select_max_modev.py

This script selects structures with large maximum force deviation from GPUMD active-learning output files.

#### Required Files

- `active.out`
- `active.xyz`

#### Usage

```sh
python select_max_modev.py <top_n> <min_deviation>
```

#### Example

```sh
python select_max_modev.py 200 0.15
```

#### Interactive Mode Example

```
205) Select max force deviation structs
Input <structs_num> <threshold> (eg. 200 0.15)
```

The output file is `selected.xyz`.



### split_train_test.py

This script inspects an extxyz dataset and creates complementary training and
test files. Before asking for any split settings, it reports:

- the input filename;
- all chemical elements found across the dataset;
- the total number of frames;
- the minimum and maximum number of atoms per frame.

#### Usage

```sh
python split_train_test.py <input.xyz>
```

#### Interactive Mode Example

```text
 >-------------------------------------------------<
 | Calling the script in Scripts/sample_structures |
 | Script: split_train_test.py                     |
 | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
 >-------------------------------------------------<
 Input <extxyz_file>
 Example: data.xyz
 The dataset summary and split options will be shown next.
 ------------>>
```

After the dataset summary, enter the test-set size in one of two forms:

| Input | Interpretation |
|-------|----------------|
| decimal strictly between `0` and `1` | fraction of all frames; `0.1` means 10% |
| positive integer | exact number of test frames; `100` means 100 frames |

Fractional frame counts are rounded to the nearest integer with halves rounded
up, and at least one frame is selected. The requested test set must be smaller
than the complete dataset so that the training set is not empty.

The selected method determines which frames form the test set:

| Method | Selection behavior | Extra input |
|--------|--------------------|-------------|
| Uniform | evenly spaced frame indices across the complete dataset | none |
| Random | random selection without replacement | optional integer seed |
| FPS | NepTrain mean atomic descriptors; start from frame 0 and repeatedly add the farthest remaining frame | compatible NEP model, e.g. `nep.txt` |

A fixed random seed makes the split reproducible; pressing Enter uses no fixed
seed. FPS uses the same NepTrain descriptor approach as function `203` and also
requires `scipy`.

For `data.xyz`, the outputs are `data_train.xyz` and `data_test.xyz`. Frames in
both outputs retain their original input order. The frames selected by the
chosen method go to the test file, and every other frame goes to the training
file.

Dependencies:

- Uniform/Random: `numpy`, `ase`;
- FPS: `numpy`, `ase`, `scipy`, `NepTrain`, and a compatible NEP model.

## Contributing

See [CONTRIBUTING.md](../../CONTRIBUTING.md) for detailed guidelines.

---

Thank you for using GPUMDkit! If you have questions about structure sampling, please open an issue on our [GitHub repository](https://github.com/zhyan0603/GPUMDkit/issues) or contact Zihan YAN (yanzihan@westlake.edu.cn).
