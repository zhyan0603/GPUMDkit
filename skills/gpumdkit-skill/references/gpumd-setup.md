# GPUMD Setup, Controls, and Actions

Use this self-contained catalog for non-ensemble setup commands, external controls, minimization, and `run`. Advanced physical parameters still require explicit user decisions.

## Contents

- Core setup
- Box and external controls
- Constraints and hybrid MD/MC
- Immediate actions
- Version handling

## Core setup

| Keyword | Current signature | Important meaning/constraint |
|---|---|---|
| `replicate` | `replicate <n_a> <n_b> <n_c>` | Positive integer replication counts; apply before loading the potential |
| `potential` | `potential <potential_filename>` | Relative/absolute file; format comes from its content; validate species/type order |
| `velocity` | `velocity <initial_temperature> [seed <seed_number>]` | Temperature in K; seed is optional. If velocities are absent and this command is omitted, GPUMD initializes at 300 K; prefer an explicit user-approved value |
| `correct_velocity` | `correct_velocity <interval> [<group_method>]` | Removes translational/rotational drift periodically; interval must be at least 10 steps |
| `time_step` | `time_step <dt_in_fs> [<max_distance_per_step>]` | Step in fs; documented default is 1 fs; optional positive displacement limiter is in Angstrom |
| `compute_extrapolation` | `compute_extrapolation asi_file <file> gamma_low <v> gamma_high <v> check_interval <n> dump_interval <n>` | ASI monitoring. Defaults: `gamma_low=0`, effectively unbounded `gamma_high`, and both intervals 1; user approval is required before relying on them in production |
| `dftd3` | `dftd3 <functional> <potential_cutoff> <coordination_number_cutoff>` | Cutoffs are in Angstrom; functional and cutoffs must match the intended D3 parameterization |
| `kspace` | `kspace <ewald|pppm>` | Reciprocal-space electrostatics method; default is `pppm`; use only for a compatible charge model |

Recognized potential families include Tersoff-1988/1989/mini, EAM, FCP, LJ, NEP, NEP+ILP, SW+ILP, Tersoff+ILP, and Deep Potential. Potential file schemas are family-specific and are not interchangeable. Do not create or rewrite a potential file from this summary; use a user-supplied validated model or obtain its exact format specification.

## Box and external controls

| Keyword | Current signature | Important meaning/constraint |
|---|---|---|
| `change_box` | `change_box <delta>`; `change_box <delta_xx> <delta_yy> <delta_zz>`; or the 6-parameter form followed by `<epsilon_yz> <epsilon_xz> <epsilon_xy>` | Immediate box change. Length increments are in Angstrom; shear terms are dimensionless; the 6-parameter form requires a triclinic box |
| `deform` | `deform <A_per_step> <deform_x> <deform_y> <deform_z>` or component-wise rates followed by three flags | Run-scoped deformation; if the parser rejects the chosen form, resolve the executable version |
| `add_force` | `add_force <group_method> <group_id> <Fx> <Fy> <Fz>` or `add_force <group_method> <group_id> <add_force_file>` | Adds force in eV/Angstrom to every atom in the selected group; file form supplies a periodic force series |
| `add_efield` | `add_efield <method> <group> <Ex> <Ey> <Ez> [charge|bec]` or `add_efield <method> <group> <field_file> [charge|bec]` | Field is in V/Angstrom. Default force uses qNEP BEC or `model.xyz` charge; explicit `charge` uses predicted/supplied charge, while `bec` requires qNEP |
| `add_spring` | `add_spring ghost_com|ghost_atom|com_com ... couple|decouple ...` | Use one of the six complete signatures below; spring constants and references are user inputs |
| `electron_stop` | `electron_stop <file>` | Loads an electronic stopping-power table; validate table units/schema before use |
| `plumed` | `plumed <plumed_file> <interval> <if_restart>` | Requires a PLUMED-enabled GPUMD build; interval is in MD steps and restart flag is 0/1 |

`add_spring` signatures:

```text
add_spring ghost_com  <method> <group> <gvx> <gvy> <gvz> couple   <k> <R0> <ox> <oy> <oz>
add_spring ghost_com  <method> <group> <gvx> <gvy> <gvz> decouple <kx> <ky> <kz> <ox> <oy> <oz>
add_spring ghost_atom <method> <group> <gvx> <gvy> <gvz> couple   <k> <R0> <ox> <oy> <oz>
add_spring ghost_atom <method> <group> <gvx> <gvy> <gvz> decouple <kx> <ky> <kz> <ox> <oy> <oz>
add_spring com_com <method> <group1> <group2> couple   <k> <R0>
add_spring com_com <method> <group1> <group2> decouple <kx> <ky> <kz>
```

An `electron_stop` table starts with `N E_min E_max`, followed by `N` evenly spaced energy rows. Each row has one stopping-power column per potential species, in the same species order as the potential file.

## Constraints and hybrid MD/MC

| Keyword | Current signature | Important meaning/constraint |
|---|---|---|
| `fix` | `fix <group_label>` or `fix <grouping_method> <group_label>` | Freezes selected atoms; one-argument form uses method 0 |
| `move` | `move <group_id> <vx> <vy> <vz>` or `move <method> <group_id> <vx> <vy> <vz>` | Constant velocity in Angstrom/fs; supported with `nvt_ber`, `nvt_nhc`, `nvt_bdp`, and `heat_lan`, not NPT; when combined with `fix`, both must use the same grouping method |
| `mc canonical` | `mc canonical <md_steps> <mc_trials> <T_i> <T_f> [group <method> <id>]` | Canonical MC/MD; temperatures are K |
| `mc sgc` | `mc sgc <md_steps> <mc_trials> <T_i> <T_f> <num_species> {<species> <mu> ...} [group ...]` | Semi-grand canonical; chemical potentials are species-specific physical inputs |
| `mc vcsgc` | `mc vcsgc <md_steps> <mc_trials> <T_i> <T_f> <num_species> {<species> <phi> ...} kappa [group ...]` | Variance-constrained SGC; `phi` and `kappa` require a defined thermodynamic protocol |

## Immediate actions

| Keyword | Current signature | Important meaning/constraint |
|---|---|---|
| `minimize` | `minimize <sd|fire> <force_tolerance> <max_steps> [<box_change> [<hydrostatic_strain>]]` | Tolerance in eV/Angstrom; box optimization currently requires FIRE |
| `run` | `run <number_of_steps>` | Positive step count; executes current run-scoped settings |

To output a minimized structure, use a zero-time-step, one-step NVE block with `dump_xyz -1 0 1 relaxed.xyz`. Load `gpumd-outputs.md` before using it.

## Version handling

These signatures are bundled with the skill. If a local parser rejects one, capture the executable version and exact error, then ask for the matching version behavior. Ask the user before selecting any physical value regardless of parser acceptance.
