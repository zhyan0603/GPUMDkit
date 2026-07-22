<div align="center">
  <h1>📚 GPUMDkit Documentation</h1>
  <p style="text-align: justify;">A practical bilingual tutorial for GPUMDkit, covering installation, command-line tools, interactive workflows, plotting, and analysis utilities.</p>
</div>

## Language

| Language | Start Here | Command Reference |
|----------|------------|-------------------|
| English | [GPUMDkit Tutorials](en/index.md) | [Command Reference](en/command_reference.md) |
| 简体中文 | [GPUMDkit 教程](zh/index.md) | [命令参考](zh/命令参考.md) |

## What You Can Learn Here

This documentation is organized as a practical tutorial rather than only a command list. New users are encouraged to start from the overview and quick start pages, then move to the topic that matches their task.

| Topic | What it Covers |
|-------|----------------|
| Quick Start | Installation, environment setup, interactive mode, direct command mode |
| Format Conversion | VASP, LAMMPS, CP2K, ABACUS, CIF, extxyz, group labels, weights |
| Sampling | Uniform/random sampling, NepTrain FPS, perturbation, force-deviation selection |
| Calculators | MSD, ionic conductivity, NEP prediction, descriptors, DOAS, NEB, polarization tools |
| Analyzers | Composition checks, distance checks, filters, property range analysis |
| Visualization | NEP training plots, MD plots, transport plots, thermal transport plots, example figures |
| Workflows | Batch DFT/MD processing and active-learning style workflows |

For complex operations, use the interactive menu:

```bash
gpumdkit.sh
```

For common operations with fixed arguments, use direct command-line mode:

```bash
gpumdkit.sh -h
gpumdkit.sh -plt -h
gpumdkit.sh -calc -h
```
