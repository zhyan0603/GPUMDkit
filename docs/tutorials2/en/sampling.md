<div align="center">
  <h1>Structure Sampling</h1>
  <p>
    <strong>English</strong> | <a href="../zh/sampling.md">简体中文</a>
  </p>
</div>

Structure selection and sampling tools.

## Available Methods

| Method | Menu | Description |
|--------|------|-------------|
| Uniform/Random | 201 | Select evenly spaced or random frames |
| FPS by NepTrain | 203 | Farthest point sampling (preferred) |
| Perturbation | 204 | Generate perturbed structures |
| Force Deviation | 205 | Select high-force-deviation structures |

## Uniform/Random Sampling

```bash
python Scripts/sample_structures/sample_structures.py <input.xyz> <method> <num> [skip]

# Examples
python Scripts/sample_structures/sample_structures.py dump.xyz uniform 50
python Scripts/sample_structures/sample_structures.py dump.xyz random 100 500
```

Output: `sampled_structures.xyz`

## Farthest Point Sampling (FPS)

```bash
python Scripts/sample_structures/neptrain_select_structs.py <sample.xyz> <train.xyz> <nep.txt>

# Example
python Scripts/sample_structures/neptrain_select_structs.py dump.xyz train.xyz nep.txt
```

Output: `selected.xyz`, `select.png`

## Structure Perturbation

```bash
python Scripts/sample_structures/perturb_structure.py <input.vasp> <num> <cell_pert> <atom_pert> <style>

# Example
python Scripts/sample_structures/perturb_structure.py POSCAR 20 0.03 0.2 uniform
```

Parameters:
- `num`: Number of perturbed structures
- `cell_pert`: Cell perturbation fraction (e.g., 0.03 = 3%)
- `atom_pert`: Atom perturbation distance (Angstrom)
- `style`: `normal`, `uniform`, or `const`

Output: `POSCAR_01.vasp`, `POSCAR_02.vasp`, ...

## Force Deviation Selection

```bash
python Scripts/sample_structures/select_max_modev.py <top_n> <min_deviation>

# Example
python Scripts/sample_structures/select_max_modev.py 200 0.15
```

Required files: `active.out`, `active.xyz`

Output: `selected.xyz`

## Frame Range Extraction

```bash
python Scripts/sample_structures/frame_range.py <input.xyz> <start> <end>

# Example: First 80% of frames
python Scripts/sample_structures/frame_range.py dump.xyz 0 0.8
```

## Dependencies

| Method | Packages |
|--------|----------|
| Uniform/Random | numpy, ase |
| FPS | numpy, ase, matplotlib, scikit-learn, scipy, NepTrain |
| Perturbation | dpdata |
| Force Deviation | numpy, ase |
