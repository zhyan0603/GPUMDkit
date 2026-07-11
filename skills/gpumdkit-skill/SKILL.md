---
name: gpumdkit-skill
description: Use GPUMDkit, GPUMD, and NEP end to end. Use when installing or navigating GPUMDkit; converting, sampling, filtering, analyzing, calculating, or plotting atomistic data; preparing batch or active-learning workflows; creating model.xyz, run.in, or nep.in; setting up, running, validating, or post-processing GPUMD simulations; training, fine-tuning, or predicting with NEP; performing diffusion, ionic-conductivity, Arrhenius, thermal-transport, or related workflows; or modifying GPUMDkit code and documentation.
---

# GPUMDkit

Use one entry point for GPUMDkit utilities, GPUMD simulation, NEP training, validation, and post-processing. Load only the references required by the current task.

## Mandatory rules

- Read the relevant reference file before giving commands, editing inputs, running tools, or interpreting results. Do not rely on memory when a local reference or help command is available.
- Ask the user whenever an unresolved choice can affect scientific meaning, numerical stability, data selection, units, cost, files, or interpretation. Never invent a plausible parameter.
- Do not choose potentials, species/type mappings, charges, temperatures, time steps, ensembles, pressures, constraints, run lengths, sampling intervals, fit windows, cutoffs, seeds, or convergence criteria without user input or an explicitly identified existing project setting.
- Distinguish documented defaults, project conventions, recommendations, and user decisions. Do not silently turn a default into a production choice.
- Do not overwrite inputs, discard structures, launch GPUMD/NEP/DFT, submit scheduler jobs, or start long/expensive calculations unless the user explicitly authorizes that action.
- Stop on parser errors, NaN/Inf values, missing outputs, unstable behavior, or unexplained warnings. Report the evidence and ask before changing scientific parameters or retrying.
- Preserve original data and unrelated worktree changes. Record transformations, filters, exclusions, commands, versions, and assumptions.

## Reference router

Read one or more references according to the task. Do not load every file by default.

| Task | Required reference |
|---|---|
| Installation, command discovery, module selection, common CLI usage | [overview.md](references/overview.md) |
| Structure/data conversion, labels, weights, replication, frame extraction | [format-conversion.md](references/format-conversion.md) |
| Uniform/random/FPS sampling, train/test splitting, perturbation, force-deviation selection | [sampling.md](references/sampling.md) |
| Batch SCF/MD preparation and active-learning workflows | [workflows.md](references/workflows.md) |
| MSD, ionic conductivity, descriptors, NEB, minimization, polarization | [calculators.md](references/calculators.md) |
| Composition, ranges, distances, filtering, outliers, probability density | [analyzers.md](references/analyzers.md) |
| NEP/MD/transport/structure plots and fit outputs | [visualization.md](references/visualization.md) |
| Creating or restyling GPUMDkit plots; project plotting style and review checklist | [plotting-style.md](references/plotting-style.md) |
| Any GPUMD simulation task; workflow and bundled-reference routing | [gpumd.md](references/gpumd.md) |
| GPUMD `model.xyz`, `run.in`, groups, units, auxiliary inputs | [gpumd-inputs.md](references/gpumd-inputs.md) |
| GPUMD potential/setup, box controls, forces, constraints, minimization | [gpumd-setup.md](references/gpumd-setup.md) |
| GPUMD NVE/NVT/NPT, MTTK, QTB, heat, PIMD, TI, shock ensembles | [gpumd-ensembles.md](references/gpumd-ensembles.md) |
| GPUMD compute keywords, signatures, outputs, compatibility | [gpumd-computes.md](references/gpumd-computes.md) |
| GPUMD dumps, active/observer modes, output write behavior | [gpumd-outputs.md](references/gpumd-outputs.md) |
| Any NEP training/prediction task; workflow and bundled-reference routing | [nep.md](references/nep.md) |
| NEP `train.xyz`/`test.xyz` schemas, targets, units, audits | [nep-data.md](references/nep-data.md) |
| All current `nep.in` parameters, defaults, ranges, and modes | [nep-parameters.md](references/nep-parameters.md) |
| NEP loss/model/restart/parity outputs and validation | [nep-outputs.md](references/nep-outputs.md) |
| Diffusion or conductivity temperature series and Arrhenius activation energy | [arrhenius.md](references/arrhenius.md) |
| GPUMDkit code, CLI, scripts, docs, tests, or skill maintenance | [contributing.md](references/contributing.md) |

For cross-module work, load every relevant reference. Examples:

- Arrhenius study: read `gpumd.md`, `gpumd-inputs.md`, `gpumd-ensembles.md`, `gpumd-computes.md`, `gpumd-outputs.md`, `arrhenius.md`, `calculators.md`, and `visualization.md`; also read `format-conversion.md` if preparing `model.xyz`.
- NEP training pipeline: read `nep.md`, `nep-data.md`, `nep-parameters.md`, `nep-outputs.md`, `format-conversion.md`, `analyzers.md`, `sampling.md`, and `visualization.md`.
- GPUMD batch sampling for active learning: read `gpumd.md`, the GPUMD parameter references used by the protocol, `workflows.md`, `sampling.md`, and `analyzers.md`.
- New GPUMDkit command: read `contributing.md` plus the reference for the affected module.
- New or restyled plot: read `visualization.md`, `plotting-style.md`, and `contributing.md`; also read the scientific reference for the quantity being plotted.

## Source priority

Use evidence in this order:

1. The self-contained GPUMD and NEP references bundled with this skill.
2. Local executable help, parser messages, and version output.
3. GPUMDkit `gpumdkit.sh -h`, module help, and target Python script `-h` output.
4. Existing project inputs confirmed to work with the same executable version.
5. Version-specific documentation supplied by the user when local behavior differs from this bundled snapshot.

If sources disagree, show the exact conflict, executable version, and parser evidence, then ask which software version governs the task. Do not guess version-specific syntax. Treat old standalone GPUMD/NEP skills as non-authoritative.

## Operating workflow

1. Classify the request and load the smallest sufficient set of references.
2. Inspect available files, executable versions, local help, and existing project conventions.
3. List missing decisions as concise questions. Continue safe read-only inspection, but do not fill scientific gaps yourself.
4. Present the resolved file plan, commands, expected outputs, and validation criteria before expensive execution.
5. Prefer `gpumdkit.sh` CLI commands when documented. Use direct scripts under `${GPUMDkit_path}/Scripts/` only when the reference identifies a menu-only/debug path or no CLI route exists.
6. When execution is authorized, capture the exact command, working directory, version, exit status, and warnings. Use a user-approved smoke test first when appropriate.
7. Validate file structure and physical behavior before downstream analysis. Separate equilibration from production and report exclusions.
8. Deliver results with units, provenance, convergence/fit evidence, uncertainty, and limitations.

## Portability

- Resolve all links relative to this skill directory. The skill does not require a product-specific activation tool.
- Use standard Markdown, YAML frontmatter, shell commands, and repository-relative resources so any Agent Skills-compatible client can load it.
- Treat `agents/openai.yaml` as optional interface metadata; never depend on it for instructions or behavior.
- Use `${GPUMDkit_path}` for repository resources when available. If it is unset, locate the repository with the user or from the current workspace instead of guessing a path.
