# Sample Structures

GPUMDkit provides various methods for sampling and selecting structures from datasets, essential for creating diverse NEP training sets and active learning workflows.

## Command-Line Usage

```bash
# PyNEP farthest point sampling
gpumdkit.sh -pynep <candidates.xyz> <train.xyz> <nep.txt>

# Other sampling methods via interactive mode only
```

## Sampling Methods

### Uniform/Random Sampling (`sample_structures.py`)

Sample structures using statistical methods.

**Interactive:** Select option `201`

**Direct execution:**
```bash
python ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py train.xyz uniform 100
python ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py train.xyz random 100
```

**Methods:**
- `uniform` - Evenly spaced samples across the dataset
- `random` - Randomly selected samples

**Output:** `sampled_structures.xyz`

### PyNEP Farthest Point Sampling (`pynep_select_structs.py`)

Select maximally diverse structures using FPS in descriptor space.

**Command-line:**
```bash
gpumdkit.sh -pynep candidates.xyz selected.xyz nep.txt
```

**Interactive:** Select option `202`

**Direct execution:**
```bash
python ${GPUMDkit_path}/Scripts/sample_structures/pynep_select_structs.py candidates.xyz train.xyz nep.txt
```

**Use case:** Create diverse training set from large candidate pool

**Requires:** `pynep` package installed

### NEPtrain Selection (`neptrain_select_structs.py`)

Alternative FPS implementation using NEPtrain.

**Interactive:** Select option `203`

**Direct execution:**
```bash
python ${GPUMDkit_path}/Scripts/sample_structures/neptrain_select_structs.py candidates.xyz train.xyz nep.txt
```

**Requires:** `neptrain` package installed

### Structure Perturbation (`perturb_structure.py`)

Generate perturbed variants of structures for training data augmentation.

**Interactive:** Select option `204`

**Direct execution:**
```bash
python ${GPUMDkit_path}/Scripts/sample_structures/perturb_structure.py structure.vasp 20 0.03 0.2 normal
```

**Parameters:**
- Input structure file (VASP format)
- Number of perturbed structures
- Cell perturbation fraction
- Atom perturbation distance (Angstroms)
- Perturbation style: `uniform`, `normal`, or `const`

**Output:** Multiple perturbed VASP files

**Use case:** Create training data variations from equilibrium structures

### Max Force Deviation Selection (`select_max_modev.py`)

Select structures with highest force deviations from active learning.

**Interactive:** Select option `205`

**Direct execution:**
```bash
python ${GPUMDkit_path}/Scripts/sample_structures/select_max_modev.py 100 0.15
```

**Parameters:**
- Number of structures to select
- Force deviation threshold (eV/Å)

**Input:** `active.xyz` from GPUMD active learning mode

**Output:** Selected structures with high uncertainty

## Interactive Mode

Access sampling tools through the interactive menu:

```bash
gpumdkit.sh
# Select: 2) Sample Structures
# Choose option 201-205
```

### Menu Options

```
 ------------>>
 201) Sample structures from extxyz
 202) Sample structures by pynep
 203) Sample structures by neptrain
 204) Perturb structure
 205) Select max force deviation structs
 000) Return to the main menu
 ------------>>
```

## Common Workflows

### Create Initial Training Set

```bash
# 1. Start with large dataset
# 2. Random sample to manageable size
python ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py large_dataset.xyz random 2000

# 3. Apply FPS for diversity
gpumdkit.sh -pynep sampled_structures.xyz train.xyz initial_nep.txt
```

### Active Learning Cycle

```bash
# 1. Run GPUMD with active learning mode
#    Add to run.in: active 100 0.15

# 2. After simulation, select uncertain structures
python ${GPUMDkit_path}/Scripts/sample_structures/select_max_modev.py 50 0.15

# 3. Run DFT on selected structures

# 4. Add to training set and retrain NEP
```

### Data Augmentation

```bash
# 1. Have equilibrium structures
# 2. Generate perturbed variants
python ${GPUMDkit_path}/Scripts/sample_structures/perturb_structure.py POSCAR 20 0.03 0.2 normal

# 3. Run DFT on perturbed structures
# 4. Add to training set
```

## Sampling Strategy Guide

| Method | When to Use | Pros | Cons |
|--------|-------------|------|------|
| **Uniform** | Initial exploration | Simple, evenly spaced | May miss important regions |
| **Random** | Quick reduction | Unbiased | May cluster |
| **PyNEP FPS** | Training set creation | Maximizes diversity | Requires NEP model |
| **NEPtrain FPS** | Alternative FPS | Also maximizes diversity | Different implementation |
| **Perturbation** | Data augmentation | Creates variations | Needs good initial |
| **Max force dev** | Active learning | Targets uncertainty | Requires active mode |

## Tips

- **Start diverse**: Begin with varied initial structures
- **Iterative refinement**: Use active learning for continuous improvement
- **Validate samples**: Check that sampled structures make sense
- **Track provenance**: Keep records of sampling parameters
- **Balance composition**: Ensure all compositions represented

## Perturbation Guidelines

For different materials:

**Rigid crystals:**
- Cell perturbation: 0.01-0.02
- Atom perturbation: 0.1-0.2 Å

**Soft materials:**
- Cell perturbation: 0.03-0.05
- Atom perturbation: 0.2-0.3 Å

**Phase transitions:**
- Cell perturbation: 0.05-0.10
- Atom perturbation: 0.3-0.5 Å

## Dependencies

- **PyNEP**: For PyNEP sampling (`pip install pynep`)
- **NEPtrain**: For NEPtrain sampling
- **ASE**: For structure manipulation
- **NumPy**: For numerical operations

---

For more details, see [Scripts/sample_structures/README.md](../../Scripts/sample_structures/README.md)
