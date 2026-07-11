# GPUMDkit

GPUMDkit is a Bash + Python CLI toolkit for GPUMD molecular dynamics and NEP potentials: data conversion/sampling, workflows, calculations, analysis, plotting, simulation, training, and post-processing. Docs: <https://gpumdkit.cn/>.

## Architecture

```text
gpumdkit.sh -> src/fN_*.sh menu wrappers -> Scripts/<module>/ implementations
```

- `gpumdkit.sh`: interactive and CLI entry point (`-h` lists commands; `-skill` locates skills).
- `Scripts/`: conversion, sampling, workflow, calculator, analyzer, plot, and utility implementations.
- `docs/tutorials/{en,zh}/`: bilingual sources; `docs/mkdocs.yml`: site config.
- `skills/gpumdkit-skill/`: canonical agent skill; `gpumdkit-skill-zh/`: Chinese version.

## Agent Routing

For any GPUMDkit, GPUMD, NEP, atomistic-data, workflow, plotting, or post-processing task:

1. Read `skills/gpumdkit-skill/SKILL.md`.
2. Classify the request and load only its routed references; combine them for cross-module work.
3. Inspect local help/files, then prefer documented `gpumdkit.sh` commands over direct scripts.

Key references under `skills/gpumdkit-skill/references/`:

- Navigation: `overview.md`.
- Data/tools: `format-conversion.md`, `sampling.md`, `workflows.md`, `calculators.md`, `analyzers.md`.
- Plots: `visualization.md`, `plotting-style.md`.
- GPUMD: start with `gpumd.md`; it routes to `gpumd-{inputs,setup,ensembles,computes,outputs}.md`.
- NEP: start with `nep.md`; it routes to `nep-{data,parameters,outputs}.md`.
- Diffusion/conductivity fits: `arrhenius.md` plus relevant simulation, calculator, and plot references.
- Code, docs, CLI, scripts, or skills: `contributing.md` plus the affected module reference.

## Rules

- Never invent consequential scientific or execution choices; ask the user.
- Simulations, training, DFT, scheduler jobs, destructive operations, and expensive runs require explicit authorization.
- Preserve unrelated changes. Before repository edits, read `contributing.md`; keep changes minimal, update both documentation languages for user-visible features, and run its validation checklist.
