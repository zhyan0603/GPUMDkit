# NEP Outputs, Restart, and Model Validation

Use this reference to interpret training/prediction files and decide whether a model is ready for deployment.

## Contents

- Write behavior and cadence
- Loss and model files
- Prediction/parity files
- Special-model outputs
- Validation checklist

## Write behavior and cadence

Bundled write behavior:

- `loss.out` is append mode;
- other outputs are continuously overwritten;
- energy/force/virial/stress/dipole/polarizability train/test outputs update every 1000 steps;
- other output files update every 100 steps.

Use a clean directory or archive old `loss.out` before a new training run.

## Loss and model files

| File | Meaning |
|---|---|
| `loss.out` | Generation, total/regularization loss, and train/test RMSE terms |
| `nep.txt` | Current trained model parameters; deployable in GPUMD |
| `nep.restart` | Continuously updated optimization state |
| checkpointed `nep_*.txt`/restart | Periodic model/state files controlled by `save_potential` |

Potential-model `loss.out` columns:

```text
gen L_t L_1 L_2 L_e_train L_f_train L_v_train L_e_test L_f_test L_v_test
```

Dipole model:

```text
gen L_t L_1 L_2 L_mu_train L_mu_test
```

Polarizability model:

```text
gen L_t L_1 L_2 L_alpha_train L_alpha_test
```

Energy/virial RMSE values are documented per atom; force RMSE is eV/Angstrom. Interpret dipole/polarizability RMSE in the dataset's recorded units.

If `nep.restart` is present, optimization starts from it. Descriptor-related hyperparameters must match the state that produced it. Do not rename or reuse a restart without provenance checks.

## Prediction and parity files

| File | Columns/meaning |
|---|---|
| `energy_train.out`, `energy_test.out` | 2 columns: predicted then target energy, eV/atom |
| `force_train.out`, `force_test.out` | 6 columns: predicted xyz then target xyz, eV/Angstrom; one row per atom |
| `virial_train.out`, `virial_test.out` | 12 columns: predicted then target xx yy zz xy yz zx, eV/atom |
| `stress_train.out`, `stress_test.out` | 12 columns: predicted then target xx yy zz xy yz zx, GPa |
| `descriptor.out` | Normalized descriptors in prediction mode when requested |

For structures without target virial/stress, the output target sentinel is `-1e6`. Exclude sentinel values explicitly and report the exclusion; do not treat them as physical labels.

Prediction mode evaluates only `train.xyz`. File row order follows the source configuration/atom order when evaluated as a full batch; preserve a mapping to source frames and verify row counts.

## Special-model outputs

| File | Meaning |
|---|---|
| `dipole_train.out`, `dipole_test.out` | 6 columns: predicted xyz then target xyz; values are normalized per atom |
| `polarizability_train.out`, `polarizability_test.out` | 12 columns: predicted then target xx yy zz xy yz zx; normalized per atom |
| `charge_train.out`, `charge_test.out` | Predicted qNEP charges |
| `bec_train.out`, `bec_test.out` | 18 columns: predicted then target 11 12 13 21 22 23 31 32 33; one row per atom |

Component ordering must be kept consistent with the training data convention. If a local executable emits a different layout, stop interpretation and resolve its version first.

## Validation checklist

- Confirm training completed or was intentionally stopped; preserve logs and restart.
- Plot loss and parity outputs with `references/visualization.md`.
- Compare train/test errors and inspect outliers by species, composition, phase, and configuration class.
- Check energy offsets, force tails, virial/stress sentinel handling, and label units.
- Evaluate errors for the intended downstream observable, not only global RMSE.
- Test extrapolation/stability in the target MD domain and inspect short-range behavior/ZBL if relevant.
- Run independent seeds/models when uncertainty or committee methods require them.
- Record model file hash, training data revision, `nep.in`, executable/source revision, and acceptance decision.
- Do not claim a universal RMSE threshold or production readiness without user-approved criteria.
