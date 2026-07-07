# GPUMDkit

GPUMDkit is a Bash + Python command-line toolkit for [GPUMD](https://gpumd.org/) and [NEP](https://gitlab.com/brucefan1983/nep) — molecular dynamics and neuroevolution potential programs. It provides format conversion, structure sampling, workflows, property calculation, plotting, and analysis utilities.

**Entry point**: `gpumdkit.sh` — dual-mode (interactive menu + CLI flags):

```bash
gpumdkit.sh                                # Interactive mode
gpumdkit.sh -h                             # Show all CLI flags
gpumdkit.sh -plt train                     # Plot NEP training
gpumdkit.sh -calc msd dump.xyz Li 10       # Calculate MSD
gpumdkit.sh -dp2xyz database train.xyz     # DeepMD npy → extxyz
gpumdkit.sh -skill                         # Show agent skill paths and install hints
```

**Documentation**: https://gpumdkit.cn/

## Agent Skills

When the user asks about GPUMDkit features, use the GPUMDkit skills. The canonical skill source is `${GPUMDkit_path}/skills`; agent tools may use global or project-local symlinks to those skills. Run `gpumdkit.sh -skill` to print the local skill path, available skills, supported install targets, and install hints.

| Skill | Use When |
|-------|----------|
| `gpumdkit-main` | General usage, navigation, entry point |
| `gpumdkit-format-conversion` | Converting structure files between formats |
| `gpumdkit-calculators` | Computing material properties |
| `gpumdkit-analyzers` | Analyzing structures and datasets |
| `gpumdkit-visualization` | Plotting simulation data |
| `gpumdkit-workflows` | Batch processing and automation |
| `gpumdkit-sampling` | Sampling and selecting structures |

## If You Need to Modify Code

**Read [skills/gpumdkit-contributing/SKILL.md](skills/gpumdkit-contributing/SKILL.md) first.** It contains all coding conventions, function templates, shell/Python patterns, and explicit rejections (things NOT to propose). The project follows a strict minimum-invasion principle — do not touch working code unless fixing a confirmed bug.

Key facts:
- Tutorials are bilingual (`docs/tutorials/en/` + `zh/`). Feature changes must update both languages.
- After Python edits: `python3 -m py_compile <file>` then clean `__pycache__`.
- After doc edits: `mkdocs build -f docs/mkdocs.yml`.
