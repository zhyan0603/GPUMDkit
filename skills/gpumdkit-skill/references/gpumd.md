# GPUMD Simulation Router and Workflow

Read this file for every GPUMD simulation task. Then load the parameter reference that covers the requested command. Do not construct `model.xyz` or `run.in` from memory.

## Contents

- Bundled reference scope
- Route the task
- Build a simulation safely
- Verify before and after execution
- Known version-sensitive points

## Bundled reference scope

The GPUMD files in this skill form a self-contained, categorized summary of the GPUMD input and output specification bundled for GPUMDkit. They define the supported keyword families, current signatures, defaults, constraints, lifecycle, and expected outputs without requiring a GPUMD source tree.

Use the local executable's parser and version output as runtime evidence. If it rejects syntax in these references, stop and report the command, input line, exact error, and executable version. Ask the user for the version-specific behavior or documentation; do not silently switch to remembered or older syntax.

## Route the task

Load these files directly from the main `SKILL.md` router:

| Need | Reference |
|---|---|
| `model.xyz`, `run.in`, groups, units, auxiliary inputs | `references/gpumd-inputs.md` |
| Potential loading, initialization, box changes, forces, constraints, minimization, run | `references/gpumd-setup.md` |
| NVE/NVT/NPT, MTTK, QTB, heat, PIMD, TI, shock and TTM integrators | `references/gpumd-ensembles.md` |
| Static calculations, MSD, RDF, transport, phonons, viscosity, LSQT | `references/gpumd-computes.md` |
| Dump/active keywords, output filenames, append/overwrite behavior | `references/gpumd-outputs.md` |
| NEP potential creation or prediction | `references/nep.md` and its parameter references |
| Diffusion/conductivity temperature series | `references/arrhenius.md` |

## Build a simulation safely

1. Confirm the observable, structure, PBC, potential file, species order, GPUMD build, and deliverables.
2. Load `gpumd-inputs.md`; validate `model.xyz` atom count, lattice, properties, units, and every group referenced by `run.in`.
3. Load the parameter references for every keyword that will appear in `run.in`.
4. Ask the user for all unresolved scientific choices: time step, temperature/pressure path, ensemble, coupling periods, constraints, stage lengths, seeds, sample intervals, correlation lengths, and fit/convergence criteria.
5. Separate minimization, equilibration, and production into explicit stages. Reissue non-propagating ensemble, compute, dump, and control commands for each `run` block where needed.
6. Calculate physical duration and output cadence explicitly:

```text
duration_ps = run_steps * time_step_fs / 1000
frame_dt_fs = dump_interval_steps * time_step_fs
```

7. List expected outputs, write mode, estimated size, and downstream GPUMDkit command before running.
8. Do not launch GPUMD until the user explicitly authorizes execution.

## Verify before and after execution

Before execution:

- ensure all referenced files exist and the potential supports every species;
- preserve old append-mode outputs or use a clean directory;
- check cell/PBC compatibility with the potential and pressure control;
- check that correlation windows fit inside the production length;
- confirm GPU/runtime and exact executable path.

After execution:

- inspect exit status and logs for parser errors, warnings, NaN/Inf, and instability;
- verify expected files, row counts, column counts, and output cadence;
- separate equilibration from production averages;
- inspect temperature, energy, pressure, cell, and trajectory behavior;
- validate observable-specific convergence before fitting or reporting.

## Known version-sensitive points

- Current `dump_xyz` starts with grouping method and group ID; do not use the obsolete interval/frame-count signature.
- `compute_phonon` currently takes only `<displacement>` and requires `kpoints.in`; older guides may show a cutoff argument.
- HNEMD requires temperature control; use Nose-Hoover chain and do not use Langevin for this purpose.
- `compute_hac` belongs in the EMD production run; it is not a one-step post-processing command for prior `compute ... jp jk` data.
- Many keywords are run-scoped and do not propagate. Reissue ensembles, computes, dumps, constraints, deformations, and driving controls for each `run` block unless their bundled entry explicitly says otherwise.
