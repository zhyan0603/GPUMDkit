<div align="center">
  <h1>📦 Structure Sampling</h1>
  <p style="text-align: justify;">Structure sampling is used to select, generate, or filter candidate configurations for NEP training and molecular dynamics workflows.</p>
</div>

**Script Location:** `Scripts/sample_structures/`

## Interactive Mode

Sampling tasks are usually easier to run from the interactive menu:

```bash
gpumdkit.sh
```

Choose:

```text
2) Sample Structures
```

The menu is:

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

Available entries:

| Menu | Method | Use Case |
|------|--------|----------|
| 201 | Uniform/random sampling | quick frame selection from a trajectory |
| 202 | PyNEP notice | prints a notice; use `gpumdkit.sh -pynep` to run PyNEP |
| 203 | NepTrain FPS | descriptor-based FPS with NepTrain |
| 204 | Perturb structure | generate perturbed structures from POSCAR/CONTCAR |
| 205 | Force-deviation selection | select structures with high model deviation |
| 206 | Train/test split | inspect and split an extxyz dataset using uniform, random, or FPS selection |

## Uniform and Random Sampling

From interactive mode, choose `201`. You will see:

```text
Input <extxyz_file> <sampling_method> <num_samples> [skip_num]
[skip_num]: number of initial frames to skip, default value is 0
Sampling_method: 'uniform' or 'random'
Example: train.xyz uniform 50
------------>>
```

Example input:

```text
dump.xyz uniform 50
```

or:

```text
dump.xyz random 100 500
```

Arguments:

| Argument | Meaning |
|----------|---------|
| `dump.xyz` | input trajectory |
| `uniform` / `random` | sampling method |
| `50` / `100` | number of selected structures |
| `500` | optional number of initial frames to skip |

Output:

- `sampled_structures.xyz`

## Farthest Point Sampling with NepTrain

This entry uses NepTrain descriptors for FPS. Choose `203` in interactive mode, or run the script directly:

```bash
python Scripts/sample_structures/parallel_neptrain_select_structs.py dump.xyz train.xyz nep.txt [threads]
```

Descriptor calculation uses one CPU worker by default. Set the optional fourth
argument to a positive integer to enable parallel calculation, for example
`4`. Each worker loads its own copy of the NEP model and uses one native OpenMP
thread, so the requested worker count is the total descriptor parallelism.

From interactive mode, choose `203`. You will see:

```text
Input <sample.xyz> <train.xyz> <nep_model> [threads]
[threads]: descriptor threads, default value is 1
Example: dump.xyz train.xyz nep.txt 4
------------>>
```

Inputs:

| File | Meaning |
|------|---------|
| `dump.xyz` | candidate structures |
| `train.xyz` | existing training set |
| `nep.txt` | current NEP model used to compute descriptors |
| `threads` | optional number of descriptor workers; default: `1` |

Outputs:

- `selected.xyz`
- `select.png`
- `pca_sample.txt`
- `pca_train.txt`
- `pca_selected.txt`

During execution, choose one of two selection modes:

1. select until the descriptor distance is below a threshold;
2. select a specified number of structures.

The script then prints a selection prompt:

```text
Choose selection method:
1) Select structures based on minimum distance
2) Select structures based on number of structures
------------>>
```

This function requires the `NepTrain` package. If you use this function, we recommend citing the NepTrain paper printed by the script.

## Deprecated PyNEP FPS

PyNEP FPS is kept for compatibility, but it is only exposed through the direct `-pynep` entry.

If you choose `202` from the interactive menu, it prints:

```text
+-------------------------------------------------+
| Function 202 is no longer supported here.       |
| PyNEP package is no longer actively maintained. |
| Please use 203) NepTrain sampling instead.      |
| If you still need PyNEP compatibility, run:     |
|                 gpumdkit.sh -pynep              |
+-------------------------------------------------+
```

### Serial PyNEP

```bash
python Scripts/sample_structures/pynep_select_structs.py dump.xyz train.xyz nep.txt
```

### Parallel PyNEP

```bash
python Scripts/sample_structures/parallel_pynep_select_structs.py dump.xyz train.xyz nep.txt 8
```

The last argument is the number of CPU threads.

### GPUMDkit Entry

```bash
gpumdkit.sh -pynep
```

Then input:

```text
dump.xyz train.xyz nep.txt 8
```

This entry requires `pynep` and is kept for compatibility.

## Structure Perturbation

Use perturbation when you want to generate initial structures around a known configuration. Choose `204` in interactive mode, or run the script directly:

```bash
python Scripts/sample_structures/perturb_structure.py POSCAR 20 0.03 0.2 uniform
```

From interactive mode, choose `204`. You will see:

```text
Input <input.vasp> <pert_num> <cell_pert_fraction> <atom_pert_distance> <atom_pert_style>
The default parameters for perturb are 20 0.03 0.2 uniform
Example: POSCAR 20 0.03 0.2 uniform
------------>>
```

Arguments:

| Argument | Meaning |
|----------|---------|
| `POSCAR` | input VASP structure |
| `20` | number of perturbed structures |
| `0.03` | cell perturbation fraction |
| `0.2` | atom perturbation distance in Angstrom |
| `uniform` | perturbation style: `normal`, `uniform`, or `const` |

Output:

- `POSCAR_01.vasp`, `POSCAR_02.vasp`, ...

This function requires `dpdata`. If you use this function, we recommend citing the dpdata package.

## Force-Deviation Selection

Use this together with the GPUMD `active` command. The `active` command uses a committee model approach: multiple potentials predict forces for the same structure, and GPUMD records the maximum force deviation. `select_max_modev.py` selects structures with large force deviations from `active.out` and `active.xyz`, then writes them to `selected.xyz`.

Choose `205` in interactive mode, or run the script directly:

```bash
python Scripts/sample_structures/select_max_modev.py 200 0.15
```

From interactive mode, choose `205`. You will see:

```text
+----------------------------------------------------+
| Select max force deviation structs from active.xyz |
|     generated by the active command in gpumd.      |
+----------------------------------------------------+
Input <structs_num> <threshold> (eg. 200 0.15)
------------>>
```

Arguments:

| Argument | Meaning |
|----------|---------|
| `200` | maximum number of structures to keep |
| `0.15` | minimum force deviation threshold |

Required files:

- `active.out`
- `active.xyz`

Output:

- `selected.xyz`

## Training and Test Set Split

Choose `206` to inspect an extxyz dataset and split it into training and test
sets. The script first reports the element list, total number of frames, and the
range of atom counts per frame. It then asks for the test-set size and selection
method.

From interactive mode:

```text
gpumdkit.sh -> 2) Sample Structures -> 206
```

Enter the input filename when prompted:

```text
Input <extxyz_file>
Example: data.xyz
------------>>
```

The test-set size accepts two forms:

| Input | Meaning |
|-------|---------|
| `0.1` | select approximately 10% of all frames for testing |
| `100` | select exactly 100 frames for testing |

For a fractional input, the number of test frames is rounded to the nearest
integer with halves rounded up, with a minimum of one test frame. The test set
must remain smaller than the complete dataset so that the training set is not
empty.

Three selection methods are available:

| Method | Behavior | Additional input |
|--------|----------|------------------|
| Uniform | selects evenly spaced frame indices across the dataset | none |
| Random | selects without replacement | optional integer random seed |
| FPS | selects diverse structures using mean atomic NEP descriptors | a compatible `nep.txt` model |

A fixed random seed makes a random split reproducible. Leaving the seed empty
uses a non-fixed seed. FPS follows the NepTrain descriptor approach used by
function `203`: it starts from the first frame and repeatedly selects the frame
farthest from the structures already selected.

For an input named `data.xyz`, the outputs are:

- `data_train.xyz`: all frames not selected for testing;
- `data_test.xyz`: the selected test frames.

Both output files preserve the input order of their retained frames. The script
can also be launched directly, but the remaining choices are still interactive:

```bash
python Scripts/sample_structures/split_train_test.py data.xyz
```

Uniform and random splitting require `numpy` and `ase`. FPS additionally
requires `NepTrain`, `scipy`, and a NEP model compatible with the elements in the
dataset.

## Frame Range Extraction

`frame_range.py` is an independent CLI tool and is not part of the interactive sampling menu `201–206`. Use it when you want to keep only a fraction of a trajectory before sampling, for example after equilibration.

```bash
gpumdkit.sh -frame_range dump.xyz 0 0.8
python Scripts/sample_structures/frame_range.py dump.xyz 0 0.8
```

This writes frames from 0% to 80% of the trajectory. The range arguments are trajectory fractions from `0` to `1`.

## Example Commands

An example preprocessing-and-sampling sequence is shown below. The first two commands are analyzer tools; the final step is structure sampling.

```bash
gpumdkit.sh -min_dist_pbc dump.xyz
gpumdkit.sh -filter_box dump.xyz 13
python Scripts/sample_structures/parallel_neptrain_select_structs.py dump.xyz train.xyz nep.txt 4
```

The final step can also be run from interactive mode by choosing `2) Sample Structures -> 203`.

## Common Mistakes

| Problem | Recommendation |
|---------|----------------|
| FPS selects too many similar structures | increase the distance threshold or select fewer structures |
| Perturbed structures are unphysical | reduce `cell_pert` and `atom_pert`, then check minimum distances |
| `active.out` and `active.xyz` do not match | regenerate both from the same GPUMD active run |
| PCA plot looks strange | check whether `nep.txt` matches the chemical species in both datasets |
| FPS split cannot load the model | check that NepTrain is installed and that `nep.txt` supports every dataset element |
| Random split changes between runs | enter a fixed integer random seed |
