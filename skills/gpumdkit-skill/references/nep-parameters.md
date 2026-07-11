# NEP Input Parameters

Use this self-contained catalog for current `nep.in` syntax, constraints, and documented defaults. Ask before changing a parameter that affects the model or fit objective.

## Contents

- File rules and minimal input
- Mode and species parameters
- Descriptor parameters
- Loss parameters
- Optimization and output parameters
- Fine-tuning and conflict checks

## File rules and minimal input

Lines use `keyword parameter...`; blank lines and text after `#` are ignored. Keywords can appear in any order except `type_weight` must follow `type`. `type` is mandatory; all others have documented defaults.

A minimal parser-level training input may contain only `type`, but a production model requires explicit user review of architecture, loss, optimizer, checkpoint, and special-physics settings.

## Mode and species parameters

| Keyword | Syntax/values | Documented default/constraint |
|---|---|---|
| `version` | `version <integer>` | Only `4` in this bundled snapshot |
| `prediction` | `prediction 0|1` | `0` optimization; `1` prediction of `train.xyz`; default 0 |
| `model_type` | `model_type 0|1|2` | `0` potential, `1` dipole, `2` polarizability; default 0 |
| `type` | `type <number_of_species> <species...>` | Mandatory; case-sensitive periodic-table symbols |
| `type_weight` | `type_weight <weight_for_each_species...>` | Exactly one non-negative force-loss weight per species in `type` order; default 1.0 each; must follow `type` |
| `charge_mode` | `charge_mode 0|1|2` | 0 original; 1 real+reciprocal qNEP; 2 reciprocal-only qNEP; default 0 |
| `atomic_v` | `atomic_v 0|1` | 0 global default; 1 atomic for dipole/polarizability; potential virial is global only |

## Descriptor parameters

| Keyword | Syntax | Documented default/constraint |
|---|---|---|
| `zbl` | `zbl <outer_cutoff>` | Absent means off; cutoff 1-3 Angstrom; inner cutoff is half |
| `use_typewise_cutoff_zbl` | `use_typewise_cutoff_zbl [<factor>]` | Absent means off; dimensionless factor default 0.7; enabled mode sets inner ZBL cutoff to zero |
| `cutoff` | `cutoff <radial> <angular>` or one radial/angular pair per species | Default 8,4 Angstrom; `3 <= angular <= radial <= 10`; cross-species cutoff is the arithmetic mean |
| `n_max` | `n_max <radial> <angular>` | Default 6,6; radial 0-12, angular 0-8 |
| `basis_size` | `basis_size <radial> <angular>` | Default 6,6; each 0-8 |
| `l_max` | `l_max <l3> [q222 q1111 q112 q123 q233 q134]` | Default `4 1` = `4 1 0 0 0 0 0`; `l3` 2-8; switches 0/1 |
| `neuron` | `neuron <number>` | Default 30; range 1-120 |

ZBL notes:

- `zbl.in` enables flexible pair-specific parameters; it needs `n(n+1)/2` pair-ordered lines and 10 values per pair.
- The `zbl` cutoff argument remains syntactically required even when `zbl.in` supplies its own cutoffs.
- `q1111` is backward-compatibility-only and strongly discouraged; newer optional four-body switches are not comprehensively tested.

## Loss parameters

| Keyword | Syntax | Documented default/constraint |
|---|---|---|
| `lambda_1` | `lambda_1 <weight>` | Non-negative; default `sqrt(N_parameters)/1000` |
| `lambda_2` | `lambda_2 <weight>` | Non-negative; default `sqrt(N_parameters)/1000` |
| `lambda_e` | `lambda_e <weight>` | Non-negative; default 1.0 |
| `lambda_f` | `lambda_f <weight>` | Non-negative; default 1.0 |
| `lambda_v` | `lambda_v <weight>` | Non-negative; default 0.1 |
| `lambda_q` | `lambda_q <weight>` | qNEP only; non-negative; default 0.1 |
| `lambda_z` | `lambda_z <weight>` | qNEP only; non-negative; default 0.5 |
| `lambda_shear` | `lambda_shear <weight>` | Extra shear-virial multiplier; effective shear weight is `lambda_v*lambda_shear`; non-negative; default 1 |
| `force_delta` | `force_delta <delta>` | Non-negative eV/Angstrom; default 0; positive values emphasize errors on smaller target forces |

Loss weights change the fitted objective and are not generic tuning knobs. Confirm available targets, intended observables, and normalization before changing them.

## Optimization and output parameters

| Keyword | Syntax | Documented default/constraint |
|---|---|---|
| `batch` | `batch <batch_size>` | Integer >=1; default 1000 |
| `population` | `population <size>` | Integer 10-100; default 50 |
| `generation` | `generation <count>` | Integer 0 to 10^7; default 100000 |
| `save_potential` | `save_potential <interval> <format> <save_restart>` | Interval default 100000; format 0 uses generation name, 1 uses timestamp/extended name (default); restart flag is 0/1 |
| `output_descriptor` | `output_descriptor 0|1|2` | Prediction only; 0 off, 1 per-structure, 2 per-atom; default 0 |

For `save_potential`, `save_restart=1` requests matching restart checkpoints. A saved `nep.restart` is required for restart/foundation workflows. If checkpoint names differ locally, resolve the executable version before automating collection.

## Fine-tuning and conflict checks

```text
fine_tune <nep_model_file> <nep_restart_file>
type <number_of_types> <species...>
```

Fine-tuning fixes `version`, `zbl`, `cutoff`, `n_max`, `basis_size`, `l_max`, and `neuron` to the source model architecture. Inspect model/restart provenance; do not substitute values from another model.

Before execution:

- compare all species and architecture values with the model/restart header;
- ensure prediction/fine-tune/restart modes are not accidentally mixed;
- ensure qNEP settings and k-space deployment requirements are understood;
- record all explicit values and implicit defaults;
- ask before using a documented default in production.
