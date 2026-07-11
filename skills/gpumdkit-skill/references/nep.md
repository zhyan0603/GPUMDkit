# NEP Training Router and Workflow

Read this file for every NEP training, fine-tuning, prediction, or deployment task. Then load the data, parameter, and output references required by the task.

## Contents

- Bundled reference scope
- Route the task
- Select the operating mode
- Training workflow
- Prediction and GPUMD deployment
- Version-sensitive points

## Bundled reference scope

The NEP files in this skill are a self-contained, categorized summary of the input data schema, all current `nep.in` keywords, defaults, constraints, output files, and deployment checks. They do not require a GPUMD source tree.

Use parser messages and executable version output as runtime evidence. If local behavior differs, stop and report the exact line, error, and version, then ask the user for the governing version-specific behavior. Do not substitute syntax from memory.

## Route the task

| Need | Reference |
|---|---|
| `train.xyz`/`test.xyz` fields, units, model-type targets, validation | `references/nep-data.md` |
| Every current `nep.in` keyword, syntax, constraints, and defaults | `references/nep-parameters.md` |
| `loss.out`, model/restart files, parity outputs, prediction outputs | `references/nep-outputs.md` |
| Conversion, filtering, sampling, plotting | Matching GPUMDkit module reference |
| Using `nep.txt` in MD | `references/gpumd.md` and GPUMD parameter references |

## Select the operating mode

Ask the user to choose and define acceptance criteria:

| Mode | Key settings/files |
|---|---|
| Potential training | `model_type 0`, `prediction 0`, `type`, `train.xyz`, optional `test.xyz` |
| Dipole training | `model_type 1`, dipole targets, global/atomic choice |
| Polarizability training | `model_type 2`, polarizability targets, global/atomic choice |
| qNEP potential | `model_type 0`, `charge_mode 1` or `2`, charge-specific loss/settings |
| Prediction | `prediction 1`, `nep.txt`, structures in `train.xyz` |
| Foundation fine-tuning | `fine_tune <model> <restart>`, compatible fixed architecture, new `type` |
| Restart interrupted training | compatible `nep.in` plus `nep.restart` in working directory |

Do not infer model type from the dataset filename. Confirm the intended downstream observable and labels.

## Training workflow

1. Preserve the raw reference data and record DFT code/settings, units, pseudopotentials, functional, and provenance.
2. Load `nep-data.md`; convert to validated extxyz and audit composition, cell, labels, ranges, minimum distances, duplicates, and train/test leakage.
3. Confirm model type, species count/order, dataset split, weighting, and missing targets.
4. Load `nep-parameters.md`; resolve every nontrivial architecture, loss, optimizer, checkpoint, charge, and ZBL decision with the user.
5. Write `nep.in`, `train.xyz`, and `test.xyz` as applicable. Record whether each value is user-selected, inherited, or a documented default.
6. Preflight working-directory files, executable path/build, GPU/runtime, expected duration, and output collision behavior.
7. Run only with explicit authorization. Capture command, working directory, source revision/build, environment, exit status, and warnings.
8. Load `nep-outputs.md`; monitor training/test losses and inspect parity outputs by property, species, composition, and configuration class.
9. Evaluate stability/extrapolation in the user's target regime before approving the model for production MD.

Do not declare success from training RMSE alone. Accuracy thresholds, dataset size, number of generations, and stopping policy are system- and application-specific.

## Prediction and GPUMD deployment

Prediction mode evaluates only structures in `train.xyz` and requires `nep.txt`. Use a clean prediction directory to avoid confusing training and prediction outputs.

Before using a model in GPUMD:

- verify species and type order against `model.xyz`;
- verify potential format/version and charge/k-space requirements;
- confirm the training domain covers composition, phase, temperature, strain, defects, surfaces, and short distances to be sampled;
- define user-approved extrapolation/stability checks;
- perform an authorized short smoke test before production.

The GPUMD command is normally:

```text
potential <path-to-nep.txt>
```

Confirm the model header, species order, and any charge/k-space requirements before use.

## Version-sensitive points

- This bundled snapshot accepts only `version 4`. If a local parser behaves differently, report the executable version and stop for version resolution.
- Current `l_max` is a three-body limit followed by optional 0/1 descriptor switches, not the old three-number NEP3 form.
- Current fine-tuning syntax is `fine_tune <nep_model_file> <nep_restart_file>`, not `fine_tune 1`.
- `nep.restart` (dot) is the current documented restart filename; older text may use `nep_restart`.
- Do not copy hyperparameter defaults from old skills. Use `nep-parameters.md` and distinguish defaults from user-approved production choices.
