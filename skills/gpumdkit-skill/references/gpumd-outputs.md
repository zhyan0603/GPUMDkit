# GPUMD Output and Dump Keywords

Use this self-contained reference to select outputs, estimate storage, and prevent accidental append contamination.

## Contents

- Trajectory and restart dumps
- Thermodynamic and model-specific dumps
- Active/observer modes
- Output files and write modes
- Output audit

## Current trajectory and restart dumps (GPUMD 5.7+)

| Keyword | Current signature | Output | Behavior |
|---|---|---|---|
| `dump_xyz` | `dump_xyz <interval> <filename> [group <grouping_method> <group_id>] [precision <single|double>] {<property> ...}` | User filename(s) | Appends extended XYZ frames; a filename ending in `*` writes one file per frame; precision defaults to single |
| `dump_netcdf` | `dump_netcdf <interval> <filename> [{<optional_args> ...}]` | User NetCDF filename | Requires the NetCDF package; positions are always written; optional arguments include properties, `group`, `precision <single|double>`, and `compression none|deflate [<level>]` |
| `dump_restart` | `dump_restart <interval>` | `restart.xyz` | Overwrites latest restart state; interval should be close to total run length (e.g. write once at the end or a few times) |
| `dump_beads` | no parameters | `beads_dump_<k>.xyz` | PIMD bead data only |

`dump_xyz` always includes wrapped positions. Optional properties currently include `mass`, `velocity`, `force`, `potential`, `virial`, `charge`, `bec`, `group_labels`, and `unwrapped_position`. Without `group`, it writes the whole system; with `group`, it restricts output to the selected group. A filename ending in `*` produces one file per frame.

## Legacy dump keywords (GPUMD before 5.7 only)

GPUMD 5.7 removed `dump_position`, `dump_velocity`, `dump_force`, and `dump_exyz`. Do not generate these keywords for GPUMD 5.7 or later. They may appear in older input files from GPUMD versions before 5.7; when the executable version is unknown, ask the user before reusing or translating them.

Typical pre-5.7 to 5.7+ replacements are:

```text
dump_position <i>                   -> dump_xyz <i> movie.xyz
dump_position <i> group <m> <g>     -> dump_xyz <i> movie.xyz group <m> <g>
dump_velocity <i>                   -> dump_xyz <i> dump.xyz velocity
dump_force <i>                      -> dump_xyz <i> dump.xyz force
dump_exyz <i> <v> <f> <p> <s>       -> dump_xyz <i> dump.xyz[*] velocity force potential
```

These are migration examples, not current commands. The `*` is needed only when separate per-frame files are wanted. Extended XYZ remains a supported data format, and `dump_observer` still has an `interval_exyz` parameter; do not remove those format/parameter references.

For trajectory-based GPUMDkit MSD, request unwrapped positions or verify that downstream unwrapping is valid for the cell and frame cadence.

## Thermodynamic and model-specific dumps

| Keyword | Current signature | Output | Constraint |
|---|---|---|---|
| `dump_thermo` | `dump_thermo <interval>` | `thermo.out` | Positive interval in MD steps |
| `dump_dipole` | `dump_dipole <interval>` | `dipole.out` | Requires a dipole-capable NEP setup |
| `dump_polarizability` | `dump_polarizability <interval>` | `polarizability.out` | Requires a polarizability model |
| `dump_shock_nemd` | `dump_shock_nemd interval <n> [bin_size <Angstrom>]` | shock spatial thermo output | Bin-size default is 10 Angstrom |

In the bundled GPUMD snapshot, `thermo.out` has 18 columns: `T K U Pxx Pyy Pzz Pyz Pxz Pxy ax ay az bx by bz cx cy cz`. Temperature is K, energies are eV, pressure components are GPa, and cell-vector components are Angstrom. Older executables can differ; resolve their version before interpreting a different layout.

## Active and observer modes

```text
active <interval> <has_velocity> <has_force> <has_uncertainty> <threshold>
dump_observer <observe|average> <interval_thermo> <interval_exyz> <has_velocity> <has_force>
```

`active` requires a committee of NEP potentials and propagates MD with the first potential. It writes `active.out` every check and appends structures above the force-uncertainty threshold to `active.xyz`. Inspect selected structures for explosions/unphysical geometry.

`dump_observer observe` propagates with the first potential and writes `observer0.out/.xyz`, `observer1.out/.xyz`, and so on for every model. `average` evaluates all supplied NEP potentials each step, propagates on their average, and writes the unnumbered `observer.out` and `observer.xyz`. All models must use the same species order.

## Output files and write modes

| Output | Generator | Typical write mode |
|---|---|---|
| `thermo.out` | `dump_thermo` | Append |
| User-named `.xyz` file | `dump_xyz` | Append |
| `restart.xyz` | `dump_restart` | Overwrite |
| `observer*.out`, `observer*.xyz` | `dump_observer` | Append |
| `active.out`, `active.xyz` | `active` | Append |
| `compute.out` | `compute` | Append |
| `ttm_electron_temperature.out` | `ensemble ttm` or `heat_ttm` | Overwrite |
| `hac.out`, `kappa.out`, `shc.out` | transport computes | Append |
| `msd.out`, `sdc.out`, `ic.out` | diffusion computes | Append |
| `rdf.out`, `adf.out`, `angular_rdf.out` | structure computes | Append |
| `mcmd.out` | `mc` | Append |
| `elastic.out`, `cohesive.out`, `D.out`, `omega2.out` | immediate/static computes | Check matching page; commonly overwrite |

`ttm_electron_temperature.out` begins with grid/active-range/source metadata. Each snapshot starts with a step marker followed by one `ix iy iz T_e` row per cell; indices are 1-based, temperature is K, and ordering is x-fastest, then y, then z. `mcmd.out` columns are MD step, MC acceptance ratio, then species concentrations in command order.

## Output audit

- Use a clean run directory or archive append-mode outputs before rerunning.
- Compute expected number of frames/rows and estimated storage before production.
- Confirm dump properties, precision, groups, and cadence match downstream analysis.
- Verify files are nonempty and have the documented columns/units.
- Detect incomplete final blocks, duplicated appended data, and mixed runs before fitting.
- Record which `run` block generated each output.
