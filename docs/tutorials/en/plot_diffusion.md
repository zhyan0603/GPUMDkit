<div align="center">
  <h1>🚶 Diffusion and Ionic Transport</h1>
  <p style="text-align: justify;">Plots for MSD, self-diffusion coefficients, and Arrhenius fits of diffusivity and ionic conductivity.</p>
</div>

← Back to the [Plot Scripts overview](plot_scripts.md)

## Diffusion and Ionic Transport

### msd (`plt_msd.py`)

Plots mean square displacement (MSD) for all directions.

**Input File:** `msd.out` with time and `MSD_x/y/z` in its first four columns.
This accepts the four-column output from `gpumdkit.sh -calc msd` and the first
four columns of GPUMD `compute_msd` output.

The slope annotations use the middle 40%-80% of the MSD series.

```bash
gpumdkit.sh -plt msd
```

<div align="center">
  <img src="../../Gallery/msd.png" alt="Mean square displacement" width="48%" />
</div>

---

### msd_all (`plt_msd_all.py`)

Plots MSD for all atomic species separately when using `all_groups` in GPUMD.

**Input File:** `msd.out` (computed with `all_groups` option)

**Requirements:** Must use `all_groups` in the `compute_msd` command in `run.in`.
For multiple groups, GPUMD appends group data; inspect that layout rather than
assuming the file has the seven columns of a single-group result.

```bash
gpumdkit.sh -plt msd_all msd.out Li P S
```

<div align="center">
  <img src="../../Gallery/msd_all.png" alt="MSD for all species" width="48%" />
</div>

---

### msd_conv (`plt_msd_convergence_check.py`)

Checks convergence of MSD calculations across different time windows.

**Input File:** `msd_step*.out` (computed with `save_every` option)

**Requirements:** Use `save_every` in the `compute_msd` command.

**Purpose:** Verify MSD has converged sufficiently for accurate diffusion coefficient calculation.

```bash
gpumdkit.sh -plt msd_conv
```

<div align="center">
  <img src="../../Gallery/msd_convergence.png" alt="MSD convergence" width="48%" />
</div>

---

### sdc (`plt_sdc.py`)

Plots self-diffusion coefficient (SDC) vs time.

**Input File:** a single-group seven-column `msd.out` from GPUMD `compute_msd`:
time, `MSD_x/y/z`, and `SDC_x/y/z`. The four-column file from `-calc msd` is
not sufficient for this plot.

```bash
gpumdkit.sh -plt sdc
```

<div align="center">
  <img src="../../Gallery/sdc.png" alt="Self-diffusion coefficient" width="48%" />
</div>

---

### msd_sdc (`plt_msd_sdc.py`)

Plots MSD and self-diffusion coefficient (SDC) side by side. The slope
annotations use the middle 40%-80% of the MSD series; the inset shows the
last 80% of the SDC data with a moving-average overlay.

**Input File:** a single-group seven-column `msd.out` from GPUMD `compute_msd`:
time, `MSD_x/y/z`, and `SDC_x/y/z`. The four-column file from `-calc msd` is
not sufficient for this plot.

```bash
gpumdkit.sh -plt msd_sdc
```

---

### D (`plt_arrhenius_d.py`)

Creates Arrhenius plot for diffusivity (log10 D vs 1000/T).

**Input Files:** `*K/msd.out` files (each temperature subdirectory should contain an `msd.out`)

```bash
gpumdkit.sh -plt arrhenius_d
gpumdkit.sh -plt D             # Alternative command
```

**Activation-energy output format:**

```text
Activation Energy: <Ea> eV, R2 = <R2>
```

<div align="center">
  <img src="../../Gallery/Arrhenius_D.png" alt="Arrhenius diffusivity" width="58%" />
</div>

---

### D_PT (`plt_arrhenius_d_PT.py`)

Creates a piecewise Arrhenius diffusivity plot around a user-specified phase-transition temperature.
The transition-temperature point is included in both the `LowT` (`T <= Tc`) and `HighT` (`T >= Tc`) fits.

**Input Files:** `*K/msd.out` files

```bash
gpumdkit.sh -plt D_PT 380
gpumdkit.sh -plt D_PT 380 save
```

The legend reports the activation energy of each branch as `HighT (xx eV)` and `LowT (xx eV)`.

<div align="center">
  <img src="../../Gallery/Arrhenius_D_PT.png" alt="Phase-transition Arrhenius diffusivity" width="58%" />
</div>

---

### D_xyz (`plt_arrhenius_d_xyz.py`)

Calculates diffusion coefficients from `msd.out` in temperature folders for the
x/y/z directions, generates a directional Arrhenius plot, and extracts the
activation energy of each component.

**Input Files:** `*K/msd.out` files

```bash
gpumdkit.sh -plt D_xyz
```

---

!!! note "Assumptions of the conductivity plots (sigma, sigma_PT, sigma_xyz)"
    These conductivity plots count the atoms of the given mobile element in the
    first temperature's `model.xyz`, use the given charge number, and apply the
    same ion count and replication across temperatures. Use them for one mobile
    species with matching MSD columns; mixed carriers or changing composition
    require a separate calculation. They fit the 40%–80% MSD interval and
    extrapolate conductivity to 300 K using the Nernst–Einstein relation.
    The temperature is taken from the folder name and the volume is the mean over
    the whole `thermo.out`, whereas `gpumdkit.sh -calc ionic-cond` averages the
    temperature and volume over the last half of `thermo.out`, so the two
    conductivities can differ slightly. Validate the diffusive interval and the
    extrapolation range before interpreting the results.

### sigma (`plt_arrhenius_sigma.py`)

Creates Arrhenius plot for ionic conductivity (log10(σ·T) vs 1000/T).

**Input Files:** `thermo.out` and `msd.out` in each `*K/` directory;
`model.xyz` and optional `run.in` in the first temperature directory

The mobile ion species and its charge number `z` are required by the Nernst-Einstein
conversion. Give them after the command (for example `Li 1`, `Na 1`, `H 1`, or `Mg 2`);
if they are omitted, the script asks for both on one line (for example `Li 1`). Only atoms
of that element in the first `model.xyz` are counted, multiplied by any `replicate` factor.
When `run.in` contains `compute_msd`, the script warns if its MSD group does not contain
exactly that species.

```bash
gpumdkit.sh -plt arrhenius_sigma Li 1
gpumdkit.sh -plt sigma Li 1 save   # Alternative command; omit "Li 1" to be prompted
```

<div align="center">
  <img src="../../Gallery/Arrhenius_sigma.png" alt="Arrhenius ionic conductivity" width="58%" />
</div>

---

### sigma_PT (`plt_arrhenius_sigma_PT.py`)

Creates a piecewise Arrhenius ionic-conductivity plot around a user-specified phase-transition temperature.
The transition-temperature point is included in both the `LowT` (`T <= Tc`) and `HighT` (`T >= Tc`) fits.
The 300 K conductivity is extrapolated from the branch on the same side of the transition as 300 K.

**Input Files:** `thermo.out` and `msd.out` in each `*K/` directory;
`model.xyz` and optional `run.in` in the first temperature directory

The element and charge come before the transition temperature, as for `sigma`;
if they are omitted (`gpumdkit.sh -plt sigma_PT 380`), the script asks for them.

```bash
gpumdkit.sh -plt sigma_PT Li 1 380
gpumdkit.sh -plt sigma_PT Li 1 380 save
```

The legend reports the activation energy of each branch as `HighT (xx eV)` and `LowT (xx eV)`.

<div align="center">
  <img src="../../Gallery/Arrhenius_sigma_PT.png" alt="Phase-transition Arrhenius ionic conductivity" width="58%" />
</div>

---

### sigma_xyz (`plt_arrhenius_sigma_xyz.py`)

Calculates the ionic conductivity from MSD and thermo data in temperature
folders for the x/y/z directions, plots the directional Arrhenius relationship,
and extracts the activation energy of each direction using the Nernst-Einstein
relation.

**Input Files:** `thermo.out` and `msd.out` in each `*K/` directory;
`model.xyz` and optional `run.in` in the first temperature directory

The element and charge arguments work as for `sigma`.

```bash
gpumdkit.sh -plt sigma_xyz Li 1
gpumdkit.sh -plt sigma_xyz Li 1 save
```
