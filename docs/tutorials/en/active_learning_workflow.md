<div align="center">
  <h1>🔁 Active Learning Workflow</h1>
  <p style="text-align: justify;">A manually reviewed NEP data iteration using GPUMDkit sampling, filtering, batch preparation, conversion, and validation tools.</p>
</div>

## Scope

This page breaks one NEP data iteration into checkable stages: generate candidate
structures, filter them, prepare reference calculations, convert the results,
update the dataset, and validate the model. GPUMDkit prepares and analyzes files;
GPUMD, DFT, and NEP calculations use the project's own execution procedure.

## Before starting an iteration

Create a dedicated working directory, preserve the previous dataset and model, and
record these decisions:

| Item | Required decision |
|---|---|
| Model | Current `nep.txt` revision and intended deployment domain |
| Candidate generation | Starting structures and GPUMD conditions |
| Selection | MD trajectory filtering followed by NepTrain FPS |
| Criteria | Distance, box, FPS stopping, and structure-count limits |
| DFT | Code, functional, pseudopotentials, convergence settings, and job budget |
| Dataset | Merge, weighting, exclusion, split, and leakage policy |
| Validation | NEP accuracy and target-domain stability acceptance criteria |

## 1. Generate candidate structures

When there are multiple starting structures or `run.in` files, first prepare GPUMD
batch directories with function 302:

~~~bash
gpumdkit.sh  # Select: 3) Workflow -> 302
~~~

This function creates directories and symlinks only. After preparation, put the
reviewed `nep.txt` and the matching `run_<index>.in` file for each sample in the
generated `md/` directory. Inspect every mapping before running GPUMD through the
project's procedure.

Use a clean output directory before a rerun, or archive GPUMD outputs that use
append mode. Before selecting structures, check exit status, logs, temperature,
energy, pressure, cell, and trajectory, and confirm that `dump.xyz` is complete
and contains no obvious unphysical structures.

## 2. Select candidate structures

For a stable MD trajectory, inspect the geometry first and then apply only the
approved filters:

~~~bash
gpumdkit.sh -min_dist_pbc dump.xyz
gpumdkit.sh -filter_dist_pbc dump.xyz <minimum_distance>
# Run this only when a box-edge criterion has been approved:
gpumdkit.sh -filter_box filtered_dump.xyz <maximum_box_edge>
gpumdkit.sh  # Select: 2) Sample Structures -> 203) FPS by NepTrain
~~~

`-min_dist_pbc` reports distances but does not remove structures;
`-filter_dist_pbc` performs the distance filter. The box filter is optional and
should be skipped when no applicable, approved box-edge criterion exists. When
used, it writes `filtered_by_box.xyz`; provide this file as the candidate input at
the 203 prompt. Record all thresholds, rejected structures, the FPS stopping rule,
and the mapping from selected structures to source frames.

### Optional committee-uncertainty alternative

Use committee uncertainty only when it is part of the project. GPUMD `active`
requires compatible committee models and produces `active.out` and `active.xyz`
from the same run. Function 205 can then rank or cap high-force-deviation
structures. Confirm the model order, check interval, output options, threshold,
`top_n`, and `min_deviation` before using it.

## 3. Prepare reference calculations

Use function 301 with the reviewed `selected.xyz`:

~~~bash
gpumdkit.sh  # Select: 3) Workflow -> 301 -> VASP or CP2K
~~~

The VASP branch creates `struct_fp/`, `fp/`, calculation directories, and
`presub.sh`. After preparation, put the reviewed `INCAR`, `POTCAR`, and `KPOINTS`
files in `fp/`. The CP2K branch, through menu 301 -> 2, accepts an extxyz file,
a reviewed template, and a directory prefix.

Before submission, check the structure count, atom and type order, cells, DFT
inputs, resource request, executable, and expected outputs. GPUMDkit prepares files;
it does not submit DFT jobs.

## 4. Convert and audit reference data

After confirming that every accepted DFT calculation completed successfully, use
the converter for the DFT code:

For VASP results:

~~~bash
gpumdkit.sh -out2xyz <scf_results_directory>
~~~

For CP2K logs and structure files, first enter the CP2K results root. The following command recursively scans `.log` files and their accompanying `.xyz`/`.inp` files under the current directory. Successful conversion writes `cp2k_exyz.xyz` there, with processing details in `Logfile.txt`.

~~~bash
gpumdkit.sh -cp2k2xyz
~~~

The VASP `-out2xyz`/menu 101 route writes `NEPdataset/train.xyz` in the current
directory by default; the CP2K route above writes `cp2k_exyz.xyz`. Replace `<new_reference.xyz>` in the following commands with the extxyz
file actually produced by the selected converter:

~~~bash
gpumdkit.sh -range <new_reference.xyz> energy
gpumdkit.sh -range <new_reference.xyz> force
gpumdkit.sh -min_dist_pbc <new_reference.xyz>
gpumdkit.sh -analyze_comp <new_reference.xyz>
~~~

Check energy, force, and virial labels, units, cells, species and type order,
composition coverage, duplicates, extreme structures, and failed or incomplete
calculations. Preserve the raw DFT outputs and record every exclusion. See [Format Conversion](format_conversion.md) for output locations and overwrite risks.

## 5. Update the dataset and validate the model

Keep the previous dataset revision before merging new reference structures. Record
which configurations were added, rejected, weighted, or assigned to training and
test sets. Where independence matters, prevent neighboring trajectory frames from
leaking across the split.

After checking the dataset, prepare or update `nep.in` using the confirmed
settings. Run NEP training through the project's execution procedure, then inspect:

- training and test losses;
- energy, force, and virial/stress parity outputs;
- errors by species, composition, phase, and configuration class;
- outliers, short-range behavior, and target-domain MD stability;
- the project's model-acceptance criteria.

Aggregate RMSE alone is not evidence that a model is ready for production. Preserve
the model and dataset revisions and record the changes before using the new model
for another MD iteration.

## Iteration record

For reproducibility, retain model and data revisions, resolved decisions, exact
commands, working directories, executable versions, exit status, warnings,
generated files, rejected structures, validation results, and known limitations.
Stop the iteration when required outputs are missing or a validation issue remains
unresolved.

## AI Assistance and Authorization Boundaries

Sampling conditions, selection criteria, DFT settings, scheduler resources,
dataset policy, and model-acceptance criteria depend on the scientific system.
An agent must not copy thresholds or simulation parameters from another material
system, or infer the potential, species/type mapping, charge, temperature, time
step, ensemble, pressure, run length, fit window, resources, or convergence
criteria.

Preparing directories, templates, or input files does not authorize execution.
GPUMD, DFT, and NEP calculations and scheduler submissions require explicit user
approval of the inputs, cost, and scientific settings. `presub.sh` is a template,
not permission to submit. Model architecture, loss, optimizer, checkpoint, ZBL,
and other nontrivial `nep.in` settings also require approval.

At each stage, an agent should check parser and log output, geometry and labels,
units, cells, species/type order, duplicates, leakage, failed outputs, and
target-domain stability. Stop and report missing outputs, parser errors, NaN/Inf
values, unphysical structures, and unexplained warnings.

Relevant skill references:

- `skills/gpumdkit-skill/references/workflows.md`
- `skills/gpumdkit-skill/references/sampling.md`
- `skills/gpumdkit-skill/references/format-conversion.md`
- `skills/gpumdkit-skill/references/nep-data.md`
- `skills/gpumdkit-skill/references/nep-parameters.md`
- `skills/gpumdkit-skill/references/nep-outputs.md`

This page does not execute a complete active-learning loop or choose the candidate,
selection, DFT, dataset, or model-acceptance procedure for the user.
