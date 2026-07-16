# GPUMD Ensembles and Integrators

Use this self-contained reference to choose the correct command family and interpret its parameters. Ensemble choice and coupling values are scientific decisions; ask the user.

## Contents

- Standard ensembles
- Pressure-control forms
- MTTK and QTB
- Heat transport and TTM
- PIMD, thermodynamic integration, and shock
- Selection and validation rules

## Standard ensembles

```text
ensemble nve
ensemble nvt_ber <T_1> <T_2> <T_coup>
ensemble nvt_nhc <T_1> <T_2> <T_coup>
ensemble nvt_bdp <T_1> <T_2> <T_coup>
ensemble nvt_lan <T_1> <T_2> <T_coup>
ensemble nvt_bao <T_1> <T_2> <T_coup>
ensemble npt_ber <T_1> <T_2> <T_coup> {<pressure_control_parameters>}
ensemble npt_scr <T_1> <T_2> <T_coup> {<pressure_control_parameters>}
```

`T_1` and `T_2` are initial/final target temperatures in K and ramp linearly during the run. `T_coup` controls thermostat coupling in timestep units. Do not select it without the user's protocol or a source-backed method.

## Pressure-control forms

The standard `npt_ber` and `npt_scr` families accept one of:

```text
# Isotropic; orthogonal cell; all directions periodic
<p_hydro> <C_hydro> <p_coup>

# Orthorhombic; periodic directions are controlled independently
<p_xx> <p_yy> <p_zz> <C_xx> <C_yy> <C_zz> <p_coup>

# Triclinic; all directions periodic; six cell components
<p_xx> <p_yy> <p_zz> <p_yz> <p_xz> <p_xy> <C_xx> <C_yy> <C_zz> <C_yz> <C_xz> <C_xy> <p_coup>
```

Pressures and elastic moduli are in GPa. The elastic constants are rough conversion parameters, not values the agent may invent. A component above 2000 GPa disables coupling for that cell component in the documented standard barostat.

## MTTK and QTB

```text
ensemble nvt_mttk temp <T_1> <T_2> [tperiod <tau_temp>]
ensemble npt_mttk temp <T_1> <T_2> <direction> <p_1> <p_2> [tperiod <tau_temp>] [pperiod <tau_press>]
ensemble nph_mttk <direction> <p_1> <p_2> [pperiod <tau_press>]

ensemble nvt_qtb <T_1> <T_2> <T_coup> [f_max <value>] [N_f <value>]
ensemble npt_qtb <direction> <p_1> <p_2> temp <T_1> <T_2> tperiod <tau_T> pperiod <tau_p> [f_max <value>] [N_f <value>]
```

MTTK directions include `iso`, `aniso`, `tri`, and individual `x`, `y`, `z`, `xy`, `yz`, `xz` components. Documented default MTTK periods are 100 timesteps for temperature and 1000 for pressure; pressure period should be at least 200 timesteps. Treat these as documented defaults, not automatic production choices. QTB `f_max` is the maximum colored-noise frequency in ps^-1 (default 200 and larger than the highest system phonon frequency); `N_f` is the frequency-point count (default 100, with `2*N_f` filter points). QTB pressure period must be at least 200 timesteps.

## Heat transport and TTM

```text
ensemble heat_nhc <T> <T_coup> <delta_T> <label_source> <label_sink>
ensemble heat_bdp <T> <T_coup> <delta_T> <label_source> <label_sink>
ensemble heat_lan <T> <T_coup> <delta_T> <label_source> <label_sink>

ensemble ttm <grouping_method> <group_id> <Ce> <rho_e> <kappa_e> <gamma_p> <gamma_s> <v_0> <nx> <ny> <nz> <T_e_init> [{optional_args}]
ensemble heat_ttm <T> <T_coup> <delta_T> <label_source> <label_sink> <grouping_method> <group_id> <Ce> <rho_e> <kappa_e> <gamma_p> <gamma_s> <v_0> <nx> <ny> <nz> <T_e_init> [{optional_args}]
```

For heat source/sink ensembles, labels refer to grouping method 0 and target temperatures are `T + delta_T` and `T - delta_T`. TTM units are: `Ce` eV/K per electron, `rho_e` Angstrom^-3, `kappa_e` eV/(ps K Angstrom), `gamma_p` and `gamma_s` amu/ps, `v_0` Angstrom/ps, and `T_e_init` K. The product `Ce*rho_e` is volumetric heat capacity in eV/(K Angstrom^3); `nx ny nz` are positive electron-grid counts.

TTM optional arguments:

| Key | Meaning/schema |
|---|---|
| `ttm_out_interval <n>` | Electron-temperature snapshot interval; default 1 step |
| `ttm_infile <file>` | Initial grid temperature; one `ix iy iz T_e` row per cell with 1-based indices |
| `ttm_properties_file <file>` | One `ix iy iz C_vol kappa_e gamma_p eta` row per cell; overrides uniform heat capacity, conductivity, and coupling |
| `ttm_source <value>` | Volumetric source in eV/(ps Angstrom^3); multiplied by per-cell `eta` when supplied |
| `ttm_active_x|y|z <range>` | `all`, one 1-based index, or inclusive `3:10`/`3-10`; inactive cells remain at zero electron temperature |

HNEMD is not an ensemble keyword. Use a temperature-controlling ensemble with `compute_hnemd`; Nose-Hoover chain is recommended and Langevin is excluded for this purpose.

## PIMD, thermodynamic integration, and shock

| Family | Current signature | Parameter meaning |
|---|---|---|
| PIMD | `ensemble pimd <num_beads> <T_1> <T_2> <T_coup> [{pressure parameters}]` | Bead count must be a positive even integer no larger than 128; set it in the first PIMD-related run and do not change it later |
| RPMD | `ensemble rpmd <num_beads>` | Real-time ring-polymer dynamics using the requested bead count |
| TRPMD | `ensemble trpmd <num_beads>` | Thermostatted ring-polymer dynamics using the requested bead count |
| Liquid TI | `ensemble ti_liquid temp <T> [tperiod <tau_T>] [tequil <n>] [tswitch <n>] [press <p>] sigmasqrd <v> p <v>` | Thermostat period defaults to 100; omitted equilibration/switch lengths are assigned in a 1:4 ratio; implemented final `p` values are 1, 25, 50, 75, 100 |
| Spring TI | `ensemble ti_spring temp <T> [tperiod <tau_T>] [tequil <n>] [tswitch <n>] [press <p>] [spring <element> <k> ...]` | Thermostat period defaults to 100; spring constants are eV/Angstrom^2 and must come last; if absent they are estimated from atomic MSD |
| Adiabatic switching | `ensemble ti_as temp <T> [tperiod <tau_T>] <pressure_control> <pmin> <pmax> [pperiod <tau_p>] [tswitch <n>] [tequil <n>]` | Isothermal pressure path; omitted switching/equilibration lengths use a 4:1 ratio |
| Reversible scaling | `ensemble ti_rs temp <tmin> <tmax> [tperiod <tau_T>] <pressure_control> <p> [pperiod <tau_p>] [tswitch <n>] [tequil <n>]` | Isobaric temperature-scaling path; omitted switching/equilibration lengths use a 4:1 ratio |
| Fixed-lambda TI | `ensemble ti lambda <lambda> temp <T> [tperiod <tau_T>] spring <element> <k> ...` | Testing integrator with fixed coupling parameter; thermostat period defaults to 100 |
| Wall piston | `ensemble wall_piston vp <vp> [thickness <thickness>]` | Piston velocity and optional wall thickness in Angstrom; thickness default 20 |
| Mirror wall | `ensemble wall_mirror vp <vp>` | Moving reflecting-wall velocity |
| Harmonic wall | `ensemble wall_harmonic vp <vp> [k <k>]` | Wall velocity and optional harmonic force constant in eV/Angstrom^2; `k` default 10 |
| MSST | `ensemble msst <x|y|z> <shock_velocity> qmass <q> mu <mu> [tscale <v>] [p0 <p0>] [v0 <v0>] [e0 <e0>]` | Shock velocity is km/s; `qmass` and `mu` are required; `tscale` defaults to 0; omitted initial-state values are calculated on the first step |
| NPHug | `ensemble nphug <direction> <p_1> <p_2> [tperiod <tau_T>] [pperiod <tau_p>] [p0 <p0>] [v0 <v0>] [e0 <e0>]` | Target pressures are in GPa and should be equal; pressure period defaults to 1000; omitted initial-state values are calculated on the first step |

The table records syntax and roles, not recommended settings. These methods require user-approved physical parameters and version-specific file schemas where auxiliary files are involved.

TI output is method-specific: liquid and spring TI write CSV data plus a YAML free-energy summary; adiabatic switching writes pressure/volume CSV data; reversible scaling writes lambda/dlambda/enthalpy CSV data. Preserve screen output and both forward/backward branches when assessing convergence.

## Selection and validation rules

- Ask whether cell shape/volume must be fixed, isotropic, anisotropic, or fully triclinic.
- Check PBC and cell-shape restrictions before pressure control.
- Ask whether temperature/pressure ramps are intentional.
- Keep ensemble commands in the same `run` block they govern; they are not propagating.
- Validate achieved temperature/pressure and cell behavior, not only parser success.
- For transport, preserve the dynamics required by the chosen method; do not substitute thermostats.
