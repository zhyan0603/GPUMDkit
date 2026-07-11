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

## Agent Skill

Use `gpumdkit-skill` whenever the user asks about GPUMDkit, GPUMD simulation, NEP training, atomistic-data conversion/sampling/analysis, plotting, workflows, or post-processing. Its `SKILL.md` routes each task to focused reference files and combines references for cross-module workflows.

The canonical source is `${GPUMDkit_path}/skills/gpumdkit-skill`. Run `gpumdkit.sh -skill` for the local path and cross-client installation hints.

A Chinese version is available at `${GPUMDkit_path}/skills/gpumdkit-skill-zh`.

The skill requires agents to ask the user about unresolved scientific or execution choices rather than inventing values. Running simulations, training, DFT, or scheduler jobs requires an explicit user request.

## If You Need to Modify Code

**Read [skills/gpumdkit-skill/SKILL.md](skills/gpumdkit-skill/SKILL.md), then [references/contributing.md](skills/gpumdkit-skill/references/contributing.md), before changing code.** They contain routing rules, coding conventions, shell/Python patterns, and explicit rejections. The project follows a strict minimum-invasion principle — do not touch working code unless fixing a confirmed bug.

Key facts:
- Tutorials are bilingual (`docs/tutorials/en/` + `zh/`). Feature changes must update both languages.
- After Python edits: `python3 -m py_compile <file>` then clean `__pycache__`.
- After doc edits: `mkdocs build -f docs/mkdocs.yml`.
