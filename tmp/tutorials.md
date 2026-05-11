# GPUMDkit Tutorials

Welcome to **GPUMDkit**.

This page is designed to help new users get productive quickly, then move to deeper topics when needed.

---

## What is GPUMDkit?

GPUMDkit is a command-line toolkit for GPUMD/NEP workflows, covering:

- **Format conversion** (VASP/CP2K/LAMMPS/CIF/extxyz...)
- **Data analysis** (quality checks, composition, distance filters)
- **Calculators** (ionic conductivity, descriptors, DOAS, etc.)
- **Plotting** (training curves, thermo/MSD/RDF, transport properties)
- **Workflow automation** (batch preparation and active learning)

---

## 5-Minute Quick Start

### 1) Installation

```bash
git clone https://github.com/zhyan0603/GPUMDkit.git
export GPUMDkit_path=/path/to/GPUMDkit
export PATH=${GPUMDkit_path}:${PATH}
source ${GPUMDkit_path}/Scripts/utils/completion.sh
source ~/.bashrc
cd ${GPUMDkit_path}
chmod +x gpumdkit.sh
```

### 2) Check whether installation is ready

```bash
gpumdkit.sh -h
```

If this prints the help table, GPUMDkit is available.

### 3) First practical task (recommended)

Convert a structure file to `extxyz`:

```bash
gpumdkit.sh -pos2exyz POSCAR model.xyz
```

---

## Two Usage Modes

### A) Interactive mode (recommended for beginners)

```bash
gpumdkit.sh
```

Use this mode if you are still exploring available functions.

### B) Command-line mode (recommended for scripts)

```bash
gpumdkit.sh -h
```

Use this mode for reproducible workflows and batch automation.

---

## Suggested Learning Path

### Path 1: New user (start here)

1. Read **Format Conversion** and run one conversion.
2. Run one **Analyzer** command on your dataset.
3. Generate one **Plot** from training or MD output.

### Path 2: Model development

1. Prepare/clean data with **Format Conversion + Analyzer**.
2. Use **Sample Structures** tools to reduce redundancy.
3. Use **Workflow/Active Learning** for iterative improvement.

### Path 3: Production workflow

1. Use command-line options only.
2. Wrap commands into shell scripts.
3. Keep each pipeline step atomic and logged.

---

## Common Tasks (copy-and-run)

### Format conversion

```bash
gpumdkit.sh -out2xyz .
gpumdkit.sh -out2exyz .
gpumdkit.sh -cif2exyz input.cif model.xyz
gpumdkit.sh -lmp2exyz dump.lammpstrj Li Y Cl
```

### Data analysis

```bash
gpumdkit.sh -range train.xyz force
gpumdkit.sh -analyze_comp train.xyz
gpumdkit.sh -min_dist_pbc train.xyz
```

### Plotting

```bash
gpumdkit.sh -plt train
gpumdkit.sh -plt thermo
gpumdkit.sh -plt msd
```

---

## Documentation Map

- [Format Conversion](format_conversion.md)
- [Sample Structures](sample_structures.md)
- [Workflow](workflow.md)
- [Active Learning](workflow_active_learning.md)
- [Calculators](calculators.md)
- [Analyzer](analyzer.md)
- [Plot Scripts](plot_scripts.md)
- [Custom Commands](custom_commands.md)

---

## Troubleshooting (quick checks)

1. **`gpumdkit.sh: command not found`**
   - Check whether `PATH` includes `${GPUMDkit_path}`.
2. **Script exists but fails with permission error**
   - Run `chmod +x gpumdkit.sh`.
3. **Unsure about parameters**
   - Use `gpumdkit.sh -h` and `gpumdkit.sh -<option> -h`.

---

## Next Step

If you're starting now, continue with:

➡️ **[Format Conversion Guide](format_conversion.md)**

