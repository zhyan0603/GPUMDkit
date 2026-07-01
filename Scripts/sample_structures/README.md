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
- **Frame extraction**: Selecting specific range from trajectories

---

## Via interactive mode

---

```
          ____ ____  _   _ __  __ ____  _    _ _
        / ___|  _ \| | | |  \/  |  _ \| | _(_) |_
       | |  _| |_) | | | | |\/| | | | | |/ / | __|
       | |_| |  __/| |_| | |  | | |_| |   <| | |_
        \____|_|    \___/|_|  |_|____/|_|\_\_|\__|

        GPUMDkit Version 1.5.6 (dev) (2026-06-17)
  Core Developer: Zihan YAN (yanzihan@westlake.edu.cn)
 Main Contributors: Denan LI, Xin WU, Zhoulin LIU & Chen HUA

 ----------------------- GPUMD -----------------------
  1) Format Conversion          2) Sample Structures
  3) Workflow                   4) Calculators
  5) Analyzer                   6) Visualization
  7) Utilities                  8) Help                
  0) Exit
 ------------>>
 Input the function number:
 2
 ------------>>
 201) Sample structures from extxyz
 202) PyNEP sampling [deprecated]
 203) FPS sampling by NepTrain [preferred]
 204) Perturb structure
 205) Select max force deviation structs
 000) Return to the main menu
 ------------>>
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

## Contributing

See [CONTRIBUTING.md](../../CONTRIBUTING.md) for detailed guidelines.

---

Thank you for using GPUMDkit! If you have questions about structure sampling, please open an issue on our [GitHub repository](https://github.com/zhyan0603/GPUMDkit/issues) or contact Zihan YAN (yanzihan@westlake.edu.cn).
