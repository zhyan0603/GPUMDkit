# Simulation and Post-processing

The unified `gpumdkit-skill` connects GPUMDkit utilities, GPUMD and NEP input preparation, controlled execution, validation, and post-processing. Its main file selects focused references for each task, including requests such as "run an Arrhenius diffusion study" or "train a NEP model and validate it in GPUMD."

## Decision boundary

An agent must ask before filling any unresolved scientific or execution choice. This includes the potential, species/type mapping, charge, temperatures, time step, ensemble, pressure, run length, output interval, fit policy, convergence criteria, scheduler resources, and overwrite behavior.

Preparing inputs does not authorize running GPUMD or NEP. Long calculations and scheduler submissions require an explicit request.

The skill bundles a portable, categorized summary of `model.xyz`, `run.in`, `nep.in`, GPUMD/NEP parameters, and output formats; users do not need a local GPUMD documentation tree. Local GPUMDkit and Python-script help govern post-processing commands. If parser behavior differs from the bundled snapshot, the agent must report the executable version and exact error, then request version-specific information instead of guessing.

## One skill, focused references

Agents activate only `gpumdkit-skill`, then load the smallest relevant set from `references/`:

| Task | References |
|---|---|
| GPUMD simulation | `gpumd.md`, then the required `gpumd-inputs.md`, `gpumd-setup.md`, `gpumd-ensembles.md`, `gpumd-computes.md`, and `gpumd-outputs.md` catalogs |
| NEP training or prediction | `nep.md`, then `nep-data.md`, `nep-parameters.md`, and `nep-outputs.md` as needed |
| Arrhenius workflow | `gpumd.md`, `arrhenius.md`, `calculators.md`, and `visualization.md` |
| GPUMDkit utility | The matching conversion, sampling, workflow, calculator, analyzer, or visualization reference |
| GPUMDkit development | `contributing.md` plus the affected module reference |

This progressive-loading structure keeps discovery simple while avoiding the context cost of loading every module for every request.

## End-to-end workflow

1. Define the observable, deliverables, and acceptance criteria.
2. Confirm the structure, potential/model, software version, and scientific parameters.
3. Prepare inputs and check species, cell, PBC, groups, units, and command syntax.
4. Present the resolved protocol and unresolved questions.
5. Run only when explicitly requested, beginning with a short approved smoke test when appropriate.
6. Check completion, thermodynamic stability, trajectory integrity, and observable-specific convergence.
7. Run GPUMDkit post-processing and report units, fit choices, uncertainty, exclusions, and limitations.

## Arrhenius example

Before setting up an Arrhenius study, confirm whether the target is diffusivity, ionic conductivity, or both; the mobile species and charge; temperature points; ensemble sequence; time step; production length; replicas; MSD sampling; and fit/uncertainty policy.

Use one integer-temperature directory per run:

```text
project/
  500K/{model.xyz,nep.txt,run.in}
  600K/{model.xyz,nep.txt,run.in}
  700K/{model.xyz,nep.txt,run.in}
```

MSD can come directly from GPUMD `compute_msd` or from an extxyz trajectory:

```bash
cd 500K
gpumdkit.sh -calc msd dump.xyz Li <frame_dt_fs> [max_corr_steps]
gpumdkit.sh -plt msd save
gpumdkit.sh -plt sdc save
```

`frame_dt_fs` is the time between stored trajectory frames:

```text
frame_dt_fs = integration_time_step_fs * dump_interval_steps
```

After validating every temperature point, run from the parent directory:

```bash
gpumdkit.sh -plt arrhenius_d save
gpumdkit.sh -plt D_xyz save
```

The current Arrhenius scripts fit 40%-80% of each MSD series. This is an established GPUMDkit choice and should not be changed without maintainer approval. The fit must still be checked against the observed diffusive regime.

For one temperature, ionic conductivity can be calculated with a user-confirmed species and charge:

```bash
gpumdkit.sh -calc ionic-cond Li 1
```

The automatic conductivity Arrhenius scripts have narrower assumptions:

```bash
gpumdkit.sh -plt arrhenius_sigma save
gpumdkit.sh -plt sigma_xyz save
```

They count only Li/Na, assume charge magnitude one, use the first temperature's composition/replication, and extrapolate the fit to 300 K. Do not use them as a general multivalent or mixed-ion workflow. Confirm the Nernst-Einstein approximation and any extrapolation with the user.

## Skill discovery

```bash
gpumdkit.sh -skill
```

The command prints the canonical skill path and installation hints. The single portable skill is located at `${GPUMDkit_path}/skills/gpumdkit-skill/`; simulation details live in its `references/gpumd.md`, `nep.md`, and `arrhenius.md` files.
