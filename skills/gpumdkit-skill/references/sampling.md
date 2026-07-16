# Structure Sampling

## Contents

- Available methods
- Command reference
- Common workflows
- Dependencies and notes

## Available Methods

| Method | Menu | Description |
|--------|------|-------------|
| Uniform/Random | 201 | Select evenly spaced or random frames |
| PyNEP notice | 202 | Prints a notice only; run PyNEP through `gpumdkit.sh -pynep` |
| FPS by PyNEP | `-pynep` | Farthest point sampling (deprecated compatibility entry) |
| FPS by NepTrain | 203 | Farthest point sampling (preferred) |
| Perturbation | 204 | Generate perturbed structures |
| Force Deviation | 205 | Select high-force-deviation structures |
| Train/Test Split | 206 | Split extxyz data using uniform, random, or NepTrain FPS selection |

## Command Reference

### Uniform/Random Sampling

```bash
# Direct Python execution
python3 ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py <input.xyz> <method> <num_samples> [skip_initial]

# Examples
python3 ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py dump.xyz uniform 50
python3 ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py dump.xyz random 100 500

# Parameters:
# - input.xyz: Input trajectory file
# - method: 'uniform' or 'random'
# - num_samples: Number of frames to select
# - skip_initial: (optional) Skip first N frames

# Output: sampled_structures.xyz
```

**How it works:**
- **Uniform**: Uses `numpy.linspace` to select evenly spaced indices across the trajectory
- **Random**: Uses `numpy.random.choice` to select random indices without replacement

### Farthest Point Sampling (FPS)

#### Using NepTrain (Preferred)

```bash
# Interactive mode
gpumdkit.sh  # Select: 2) Sample Structures -> 203

# Semi-interactive direct Python execution
python3 ${GPUMDkit_path}/Scripts/sample_structures/parallel_neptrain_select_structs.py <sample.xyz> <train.xyz> <nep.txt> [threads]

# Example
python3 ${GPUMDkit_path}/Scripts/sample_structures/parallel_neptrain_select_structs.py dump.xyz train.xyz nep.txt 4

# Selection methods (interactive):
# 1. Minimum distance: Select until max distance < threshold
# 2. Number of structures: Select N structures

# Output:
# - selected.xyz
# - select.png (PCA visualization)
# - pca_sample.txt, pca_train.txt, pca_selected.txt
```

`threads` is optional and defaults to `1`. Set a positive integer to enable
parallel descriptor calculation. Each worker loads an independent NEP model,
uses one native OpenMP thread, and returns results in input frame order.

**Algorithm**: Farthest Point Sampling using NEP descriptors (mean descriptor per structure)

**Citation**: Chen et al., Comput. Phys. Commun., 2025, 317, 109859

#### Using PyNEP (Deprecated)

```bash
# GPUMDkit entry
gpumdkit.sh -pynep

# Note: PyNEP package is no longer actively maintained
# Use NepTrain (method 203) instead
```

### Structure Perturbation

```bash
# Interactive mode
gpumdkit.sh  # Select: 2) Sample Structures -> 204

# Direct Python execution
python3 ${GPUMDkit_path}/Scripts/sample_structures/perturb_structure.py <input.vasp> <pert_num> <cell_pert> <atom_pert> <style>

# Example
python3 ${GPUMDkit_path}/Scripts/sample_structures/perturb_structure.py POSCAR 20 0.03 0.2 uniform

# Parameters:
# - input.vasp: POSCAR/CONTCAR file
# - pert_num: Number of perturbed structures (e.g., 20)
# - cell_pert: Cell perturbation fraction (e.g., 0.03 = 3%)
# - atom_pert: Atom perturbation distance in Angstrom (e.g., 0.2)
# - style: 'normal', 'uniform', or 'const'

# Output: POSCAR_01.vasp, POSCAR_02.vasp, ..., POSCAR_<N>.vasp
```

**How it works**: Uses `dpdata.System.perturb()` to generate perturbed variants of the input structure, varying both cell vectors and atomic positions.

### Force Deviation Selection

```bash
# Interactive mode
gpumdkit.sh  # Select: 2) Sample Structures -> 205

# Direct Python execution
python3 ${GPUMDkit_path}/Scripts/sample_structures/select_max_modev.py <top_n> <min_deviation>

# Example
python3 ${GPUMDkit_path}/Scripts/sample_structures/select_max_modev.py 200 0.15

# Required files in current directory:
# - active.out (from GPUMD active command)
# - active.xyz (corresponding structures)

# Output: selected.xyz
```

**How it works**:
1. Reads `active.out` and filters structures with force deviation above `min_deviation`
2. Sorts by deviation (descending) and extracts the top N structures
3. Reads corresponding frames from `active.xyz` and writes them to output

### Training and Test Set Split

```bash
# Interactive mode
gpumdkit.sh  # Select: 2) Sample Structures -> 206

# Direct script entry; the split settings are prompted interactively
python3 ${GPUMDkit_path}/Scripts/sample_structures/split_train_test.py <input.xyz>

# Example
python3 ${GPUMDkit_path}/Scripts/sample_structures/split_train_test.py data.xyz
```

The script reports the dataset elements, frame count, and atom-count range before
asking for the test-set size:

- A value strictly between `0` and `1` is a fraction, e.g. `0.1` means 10%.
- A positive integer is an exact frame count, e.g. `100` means 100 frames.
- Fractional counts are rounded to the nearest integer with halves rounded up,
  with a minimum of one test frame.
- The test set must contain fewer frames than the complete dataset.

Selection methods:

- **Uniform**: evenly spaced indices over the full dataset.
- **Random**: sampling without replacement; an optional integer seed makes the
  result reproducible.
- **FPS**: mean atomic NEP descriptors calculated by NepTrain, starting from the
  first frame and repeatedly adding the farthest remaining frame. This requires
  a compatible NEP model such as `nep.txt`.

For `data.xyz`, the outputs are `data_train.xyz` and `data_test.xyz`. Structures
in each output retain their original input order.

### Frame Range Extraction

```bash
# GPUMDkit CLI
gpumdkit.sh -frame_range <input.xyz> <start_fraction> <end_fraction>

# Direct Python execution also works
python3 ${GPUMDkit_path}/Scripts/sample_structures/frame_range.py <input.xyz> <start_fraction> <end_fraction>

# Example: Extract first 80% of frames
gpumdkit.sh -frame_range dump.xyz 0 0.8

# Output: dump_0.00_0.80.xyz
```

## Common Workflows

### Training Data Preparation

```bash
# 1. Run MD simulation
# ... generate trajectory ...

# 2. Filter by distance
gpumdkit.sh -min_dist_pbc dump.xyz
gpumdkit.sh -filter_dist_pbc dump.xyz 1.0

# 3. Sample diverse structures
gpumdkit.sh  # Select: 2) Sample Structures -> 203

# 4. Add to training set
cat selected.xyz >> train.xyz
```

### Active Learning Data Selection

```bash
# 1. Run MD with current NEP
# ... generate active.out and active.xyz ...

# 2. Select high-deviation structures
python3 ${GPUMDkit_path}/Scripts/sample_structures/select_max_modev.py 100 0.1

# 3. Or use FPS for diversity
python3 ${GPUMDkit_path}/Scripts/sample_structures/parallel_neptrain_select_structs.py active.xyz train.xyz nep.txt 4
```

### Structure Perturbation for Initial Data

```bash
# 1. Start from equilibrium structure
# 2. Generate perturbations
python3 ${GPUMDkit_path}/Scripts/sample_structures/perturb_structure.py POSCAR 50 0.05 0.3 uniform

# 3. Run DFT on perturbed structures
# ... DFT calculations ...

# 4. Add to training set
```

### PyNEP Compatibility Path

```bash
# gpumdkit.sh -pynep calls the parallel PyNEP compatibility script.
gpumdkit.sh -pynep

# Example input after launching:
dump.xyz train.xyz nep.txt 8
```

## Dependencies

| Method | Required Packages |
|--------|------------------|
| Uniform/Random | `numpy`, `ase` |
| FPS (NepTrain) | `numpy`, `ase`, `matplotlib`, `scikit-learn`, `scipy`, `NepTrain` |
| FPS (PyNEP) | `numpy`, `ase`, `matplotlib`, `scikit-learn`, `pynep` |
| Perturbation | `dpdata` |
| Force Deviation | `numpy`, `ase` |
| Train/Test Split (Uniform/Random) | `numpy`, `ase` |
| Train/Test Split (FPS) | `numpy`, `ase`, `scipy`, `NepTrain` |

Install dependencies:
```bash
pip install numpy ase matplotlib scikit-learn scipy dpdata
pip install neptrain  # For FPS by NepTrain
```

## Notes

1. **FPS is preferred**: Farthest Point Sampling provides better diversity than uniform/random
2. **NepTrain over PyNEP**: NepTrain is actively maintained; PyNEP is deprecated
3. **Filter before sampling**: Remove unphysical structures (too close, too large) first
4. **Frame ranges use fractions**: `-frame_range` uses start/end fractions between 0.0 and 1.0
5. **Output is extxyz**: All sampling methods output extxyz format
6. **Split selection targets the test set**: Function 206 selects test frames;
   every unselected frame is written to the training set

## Detailed Documentation

See `${GPUMDkit_path}/docs/tutorials/en/structure_sampling.md` or the Chinese counterpart for the user-facing guide.
