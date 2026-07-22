# GPUMD Computation Keywords

Use this self-contained reference to find a computation's signature, output, compatibility rules, and parameter roles. Ask before choosing physical or convergence parameters.

## Contents

- General and structural calculations
- Dynamics and transport
- Phonons and modal transport
- Electronic and order calculations
- Compatibility and convergence

## General and structural calculations

| Keyword | Current signature | Main output |
|---|---|---|
| `compute` | `compute <grouping_method> <sample_interval> <output_interval> {<quantity>}` | `compute.out` |
| `compute_chunk` | `compute_chunk <sample_interval> <output_interval> bin/1d|bin/2d|bin/3d <bin parameters> {<quantity>}` | `compute_chunk.out` |
| `compute_adf` | `compute_adf <interval> <num_bins> <rc_min> <rc_max>` or typed multi-triplet form | `adf.out` |
| `compute_rdf` | `compute_rdf <cutoff> <num_bins> <interval>` | `rdf.out` |
| `compute_angular_rdf` | `compute_angular_rdf <cutoff> <r_bins> <angle_bins> <interval> [atom <i> <j> ...]` | `angular_rdf.out` |
| `compute_cohesive` | `compute_cohesive <scale_start> <scale_end> <direction>` | `cohesive.out` |
| `compute_elastic` | `compute_elastic <strain_value>` | `elastic.out` |
| `compute_phonon` | `compute_phonon <displacement>` | `D.out`, `omega2.out` |

For `compute_cohesive`, direction values 0-6 mean x, y, z, xy, yz, zx, and xyz scaling; the number of points is `(scale_end-scale_start)*1000+1`. `compute_elastic` uses a dimensionless finite strain. Both must follow `potential`.

`compute_phonon` requires `kpoints.in`; place `replicate` at the head of `run.in` when replication is needed, then load the potential before this immediate calculation. Displacement is in Angstrom. The current signature has no cutoff argument.

Common sampling rules:

- `sample_interval` is the number of MD steps between samples. `output_interval` is the number of samples averaged per write, so one block is written every `sample_interval * output_interval` steps.
- `compute` accepts one or more distinct quantities: `temperature`, `potential`, `force`, `virial`, `jp` (potential heat current), `jk` (kinetic heat current), and `momentum`. Output order follows command order and is resolved per group in the selected grouping method.
- `compute_chunk` bins atoms from their current positions. For each axis use `<dim> lower <delta>`, where dimension is `x`, `y`, or `z`, origin is currently only `lower`, and positive `delta` is the bin width in Angstrom. Axes must be distinct in 2D/3D.
- `compute_chunk` quantities are `temperature` (K), `density/number` (Angstrom^-3), `density/mass` (amu/Angstrom^3), `vx vy vz` (Angstrom/fs), and `fx fy fz` (eV/Angstrom). Each output block has one line per chunk: zero-based ID, bin center coordinate(s), average atom count, then requested values.
- Structural cutoffs, bin counts, strain amplitudes, phonon displacements, and element/type selections alter resolution or physics. Ask the user instead of copying example values.

## Dynamics and transport

| Keyword | Current signature | Main output |
|---|---|---|
| `compute_msd` | `compute_msd <sample_interval> <Nc> [group <method> <id> | all_groups <method>] [save_every <interval>]` | `msd.out`, snapshots |
| `compute_sdc` | `compute_sdc <sample_interval> <Nc> [group <method> <id>]` | `sdc.out` |
| `compute_ic` | `compute_ic <sample_interval> <Nc> <type_index> <charge>` | `ic.out` |
| `compute_hac` | `compute_hac <sampling_interval> <correlation_steps> <output_interval>` | `hac.out` |
| `compute_hnemd` | `compute_hnemd <output_interval> <Fe_x> <Fe_y> <Fe_z>` | `kappa.out` |
| `compute_hnemdec` | `compute_hnemdec <drive_type> <output_interval> <Fe_x> <Fe_y> <Fe_z>` | `onsager.out` |
| `compute_shc` | `compute_shc <sample_interval> <Nc> <0|1|2> <num_omega> <max_omega> [group <method> <id>]` | `shc.out` |
| `compute_viscosity` | `compute_viscosity <sampling_interval> <correlation_steps>` | `viscosity.out` |

Interpretation rules:

- `sample_interval` is in MD steps; correlation time is `sample_interval * Nc * time_step`.
- `Nc` is the maximum number of correlation samples. The first complete correlation output requires at least `sample_interval * Nc` MD steps.
- `compute_msd` output contains time, directional MSD, and directional SDC for each selected group.
- `compute_ic` uses a zero-based potential type index and an ionic charge supplied by the user; verify species/type order and state the ionic-conductivity definition used in the analysis.
- `compute_hac` runs during an EMD production trajectory, normally NVE after equilibration.
- `compute_hnemd` requires temperature control; use one nonzero driving-force component unless the method explicitly requires otherwise.
- HNEMD/HNEMA driving-force components are in Angstrom^-1. For `compute_hnemdec`, `drive_type=0` selects thermal driving with Angstrom^-1 force; a positive integer `i` selects diffusive driving for the i-th species in the `nep.txt` header with force in eV/Angstrom. Langevin thermostats are incompatible with both HNEMD and HNEMDEC dynamics.
- `compute_shc` requires `1 <= sample_interval <= 10`, `100 <= Nc <= 1000`, and direction 0/1/2 for x/y/z. `max_omega` is THz. `group <method> -1` calculates every nonzero group and can be expensive.
- Driving force, correlation length, and production duration require convergence tests chosen with the user.

## Phonons and modal transport

| Keyword | Current signature | Main output |
|---|---|---|
| `compute_dos` | `compute_dos <sample_interval> <Nc> <omega_max> [group <method> <id>] [num_dos_points <n>]` | `mvac.out`, `dos.out` |
| `compute_gkma` | `compute_gkma <sample_interval> <first_mode> <last_mode> <bin_size|f_bin_size> <size>` | `heatmode.out` |
| `compute_hnema` | `compute_hnema <sample_interval> <output_interval> <Fe_x> <Fe_y> <Fe_z> <first_mode> <last_mode> <bin_option> <size>` | `kappamode.out` |

GKMA/HNEMA require a matching `eigenvector.in`; mode indices and binning must match that file. `compute_gkma` and `compute_hnema` cannot be active in the same run; the last one wins. Modal calculations can require substantial memory and disk.
For HNEMA, `sample_interval` must divide `output_interval`. `bin_size` groups a fixed number of modes; `f_bin_size` groups modes by frequency width.

## Electronic, dipole, and order calculations

| Keyword | Current signature | Main output |
|---|---|---|
| `compute_dpdt` | `compute_dpdt <sampling_interval>` | `dpdt.out` |
| `compute_lsqt` | `compute_lsqt <direction> <num_moments> <num_energies> <E_1> <E_2> <E_max>` | `lsqt_dos.out`, `lsqt_velocity.out`, `lsqt_sigma.out` |
| `compute_orientorder` | `compute_orientorder <interval> <cutoff|nnn> <mode_value> <ndegrees> <l...> [<average> <wl> <wlhat>]` | `orientorder.out` |

For `compute_orientorder`, `cutoff` uses a distance and `nnn` uses a neighbor count. Optional `average`, `wl`, and `wlhat` flags are 0/1 and control neighbor averaging and third-order invariants.

LSQT accepts transport direction `x`, `y`, or `z`; energies are in eV. It is preliminary and currently supports only the hard-coded carbon tight-binding model. `E_max` must be slightly larger than the model's maximum absolute energy and requires a user-approved convergence procedure.

## Compatibility and convergence

- `compute_dos` and `compute_sdc` cannot be used in the same run.
- For `compute_dos`, angular-frequency maximum is in THz, `num_dos_points` defaults to `Nc`, and `group <method> -1` requests all groups.
- Compute keywords are generally run-scoped/non-propagating; reissue them for each intended production block.
- Ensure `Nc * sample_interval` is comfortably shorter than the production run and verify convergence with saved windows/replicas where applicable.
- Confirm output append behavior before reruns.
- If output columns or parser behavior differ from this snapshot, stop and resolve the executable version before interpreting data.
