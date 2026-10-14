<div align="center">
  <h1>📚 GPUMDkit Documentation</h1>
  <p style="text-align: justify;">Bilingual tutorials for GPUMDkit: installation, interactive menus, command-line tools, plotting, and analysis.</p>
</div>

## Language

| Language | Start Here | Command Reference |
|----------|------------|-------------------|
| English | [GPUMDkit Tutorials](en/index.md) | [Command Reference](en/command_reference.md) |
| 简体中文 | [GPUMDkit 教程](zh/index.md) | [命令参考](zh/命令参考.md) |

## What You Can Learn Here

The pages are organized by task. New users should start with Quick Start, then open the topic that matches their work.

| Topic | What it Covers |
|-------|----------------|
| Quick Start | Installation, environment setup, interactive mode, direct command mode |
| Remote Web Console | Browse and manage remote files through an SSH tunnel |
| Format Conversion | VASP, LAMMPS, CP2K, ABACUS, CIF, extxyz, group labels, weights |
| Sampling | Uniform/random sampling, NepTrain FPS, perturbation, force-deviation selection |
| Calculators | MSD, ionic conductivity, XRD, NEP prediction, descriptors, DOAS, NEB, polarization tools |
| Analyzers | Composition checks, distance checks, filters, property range analysis |
| Visualization | NEP training, MD, diffusion, and thermal transport plots, with example figures |
| Workflows | Batch DFT/MD processing and active-learning style workflows |

Use the interactive menu for guided, multi-step tasks:

```bash
gpumdkit.sh
```

Use command-line mode when the arguments are known:

```bash
gpumdkit.sh -h
gpumdkit.sh -plt -h
gpumdkit.sh -calc -h
```
