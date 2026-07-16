# Arrhenius Diffusion and Ionic-Conductivity Workflow

Use this reference when the user asks for an Arrhenius task, activation energy, temperature-dependent diffusivity, or ionic conductivity. This workflow joins GPUMD simulation setup with GPUMDkit analysis; it does not authorize choosing scientific parameters or launching simulations.

## Contents

- Clarify the target
- Plan the temperature series
- Prepare GPUMD inputs
- Produce and validate MSD
- Run GPUMDkit analysis
- Interpret and report
- Current tool limitations

## Clarify the target

Ask the user to resolve:

- diffusivity `D`, ionic conductivity `sigma`, or both;
- mobile species, group definition, and formal/effective charge;
- isotropic total or directional result;
- temperature list and whether a single Arrhenius regime is expected;
- structure, potential, pressure/volume condition, and phase stability;
- ensemble sequence, time step, equilibration length, production length, and replicas;
- MSD sampling interval/correlation length or trajectory frame interval;
- fit acceptance, uncertainty method, and whether extrapolation (for example to 300 K) is allowed.

Never infer these from a folder name, chemical formula, or common practice. If the user only says "do an Arrhenius calculation," ask these questions before creating production inputs.

## Plan the temperature series

Use one independent directory per integer temperature because the current plotting scripts discover names matching `<integer>K`, for example:

```text
project/
  500K/
    model.xyz
    nep.txt
    run.in
  600K/
  700K/
```

Use identical composition and a controlled, user-approved protocol across temperatures unless the study design says otherwise. Record seeds and replicas. Two valid temperatures are the mathematical minimum for a line, but they are weak scientific evidence; ask the user for an adequate design rather than selecting a count yourself.

Do not assume all temperatures belong to one linear regime. Check phase changes, mechanism changes, nonlinearity, and unstable trajectories before fitting.

## Prepare GPUMD inputs

For direct GPUMD MSD, add a group mapping to `model.xyz` when only one species/subset is mobile, then use a verified command such as:

```text
compute_msd <sample_interval> <Nc> group <group_method> <group_id> [save_every <interval>]
```

The exact ensemble and run lengths must come from the user-approved design. Reissue non-propagating compute/dump commands in the correct production run.

For trajectory-based GPUMDkit MSD, output extxyz frames with enough metadata for periodic unwrapping. Pass the time between stored frames to `-calc msd`; calculate it as:

```text
frame_dt_fs = integration_time_step_fs * dump_interval_steps
```

Confirm this value from `run.in` and the actual trajectory cadence. Do not pass the integration time step when frames were dumped less frequently.

## Produce and validate MSD

Two supported paths are:

### GPUMD native MSD

Run GPUMD with `compute_msd`. The native `msd.out` contains time, three MSD columns, and three SDC columns for one group. GPUMDkit Arrhenius scripts consume the first four columns.

### GPUMDkit trajectory MSD

From each temperature directory:

```bash
gpumdkit.sh -calc msd <trajectory.xyz> <element> <frame_dt_fs> [max_corr_steps]
```

This writes a four-column `msd.out` in the current directory.

Before fitting any point:

- inspect the trajectory and thermodynamic stability;
- verify the selected atoms and periodic unwrapping;
- confirm a physically meaningful diffusive linear regime rather than ballistic, caged, or saturated behavior;
- inspect `gpumdkit.sh -plt msd`, `-plt sdc`, and `-plt msd_conv` where applicable;
- compare independent replicas or block estimates according to the user's uncertainty plan;
- reject or retain questionable points only with an explicit documented decision.

## Run GPUMDkit analysis

### Diffusivity Arrhenius fit

From the parent containing `<integer>K/msd.out`:

```bash
gpumdkit.sh -plt arrhenius_d save
gpumdkit.sh -plt D_xyz save
```

The current scripts fit the middle 40%-80% of each MSD series, convert the slope from Angstrom^2/ps to cm^2/s, fit `log10(D)` versus `1000/T`, and report activation energy in eV. `D_xyz` reports total and directional values.

The 40%-80% range is an established GPUMDkit project choice. Do not change it arbitrarily. If it does not match the observed diffusive regime, show the evidence and ask the user how to proceed.

### Per-temperature ionic conductivity

In a temperature directory containing `msd.out`, `thermo.out`, `model.xyz`, and optionally `run.in`:

```bash
gpumdkit.sh -calc ionic-cond <element> <charge>
```

This uses the Nernst-Einstein relation, reads temperature/volume, counts the requested species, accounts for `replicate` when found, and prints directional/total diffusivity and conductivity. Confirm that independent-ion Nernst-Einstein conductivity is the intended observable; it does not include correlation effects unless the method is extended.

### Conductivity Arrhenius fit

For monovalent Li/Na systems that satisfy the current script assumptions:

```bash
gpumdkit.sh -plt arrhenius_sigma save
gpumdkit.sh -plt sigma_xyz save
```

These commands fit `ln(sigma*T)` versus `1000/T` and report activation energies. Inspect the current limitations below before using them.

## Interpret and report

Report at minimum:

- temperatures included/excluded and why;
- ensemble, production duration, sampling, species/group, and number of replicas;
- MSD fit interval (40%-80% for current GPUMDkit scripts) and evidence that it is diffusive;
- `D` or `sigma` values with units and uncertainty;
- Arrhenius equation, transformed axes, activation energy, fit quality, and uncertainty;
- whether any value is interpolated or extrapolated;
- Nernst-Einstein assumptions for conductivity;
- phase/mechanism changes or non-Arrhenius behavior.

Do not report a precise activation energy from an unvalidated fit merely because the plotting command completed.

## Current tool limitations

- `plt_arrhenius_d.py` and directional variant discover only integer `<temperature>K` directories and require usable `msd.out` files.
- Arrhenius scripts use a fixed 40%-80% MSD fit interval.
- `plt_arrhenius_sigma.py` and `plt_arrhenius_sigma_xyz.py` count only Li and Na from the first temperature's `model.xyz` and hard-code charge magnitude `z = 1`.
- Conductivity plotting reads replication from the first temperature's `run.in` and assumes the same ion count across the series.
- The conductivity scripts extrapolate to 300 K. Ask the user whether that extrapolation is scientifically acceptable before presenting it as a result.
- For other species, multivalent ions, mixed mobile species, changing compositions, or correlated conductivity, do not use the automatic conductivity Arrhenius scripts as though they were general. Use per-temperature verified calculations and a user-approved fitting method, or ask for a code change.

