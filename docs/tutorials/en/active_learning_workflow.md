<div align="center">
  <h1>🔁 Active Learning Workflow</h1>
  <p style="text-align: justify;">A manually reviewed NEP data-iteration protocol built from GPUMDkit sampling, filtering, batch-preparation, conversion, and validation tools.</p>
</div>

## Scope

This page intentionally documents a staged, manually reviewed iteration.
Sampling conditions, selection criteria, DFT settings, scheduler resources,
dataset policy, and model acceptance criteria depend on the scientific system.

GPUMDkit supports the individual preparation and analysis stages. GPUMD, DFT,
and NEP calculations are run separately only after their inputs and cost have
been approved.

## Before starting an iteration

Define and record:

| Item | Required decision |
|---|---|
| Model | Current `nep.txt` revision and intended deployment domain |
| Candidate generation | Starting structures and GPUMD conditions |
| Selection | MD trajectory filtering followed by NepTrain FPS |
| Criteria | Distance, box, FPS stopping, and structure-count limits |
| DFT | Code, functional, pseudopotentials, convergence settings, and job budget |
| Dataset | Merge, weighting, exclusion, split, and leakage policy |
| Validation | NEP accuracy and target-domain stability acceptance criteria |

Do not copy thresholds or simulation parameters from another material system.

## 1. Generate candidate structures

Prepare GPUMD batch directories with function 302 when multiple starting
structures or `run.in` files are needed:

```bash
gpumdkit.sh  # Select: 3) Workflow -> 302
```

The command prepares directories and symlinks. Put the approved `nep.txt` and
`run_<index>.in` files in the generated `md/` directory, inspect every mapping,
and run GPUMD through your own approved execution procedure.

Use a clean output directory or archive append-mode GPUMD outputs before a
rerun. Check the exit status, log, temperature, energy, pressure, cell, and
trajectory before selecting structures.

## 2. Select candidates

The recommended route is direct MD sampling followed by geometry filtering and
NepTrain FPS.

### MD trajectory and NepTrain FPS

For an ordinary stable MD trajectory, first inspect periodic distances, then
apply only approved filters:

```bash
gpumdkit.sh -min_dist_pbc dump.xyz
gpumdkit.sh -filter_dist_pbc dump.xyz <minimum_distance>
gpumdkit.sh -filter_box filtered_dump.xyz <maximum_box_edge>
gpumdkit.sh  # Select: 2) Sample Structures -> 203) FPS by NepTrain
```

`-min_dist_pbc` reports distances but does not remove structures. Use
`-filter_dist_pbc` for actual filtering. The box filter is optional and should
be skipped when no approved box-edge criterion applies. Record all thresholds,
rejected structures, and the NepTrain FPS stopping rule.

### Optional committee-uncertainty alternative

Use committee uncertainty only when it is explicitly part of the project.
GPUMD `active` requires compatible committee models and produces `active.out`
and `active.xyz`; function 205 can then rank or cap high-force-deviation
structures. Confirm all uncertainty settings and keep the two files from the
same run.

## 3. Prepare reference calculations

Use function 301 with the reviewed `selected.xyz`:

```bash
gpumdkit.sh  # Select: 3) Workflow -> 301 -> VASP or CP2K
```

For VASP, the tool creates `struct_fp/`, `fp/`, calculation directories, and a
`presub.sh` template. Put the approved `INCAR`, `POTCAR`, and `KPOINTS` files in
`fp/` after preparation. For CP2K, provide an extxyz file, a reviewed template,
and a prefix through menu 301 -> 2.

Before submission, verify the structure count, atom/type order, cells, DFT
inputs, resource request, executable, and expected output. GPUMDkit does not
submit the DFT jobs.

## 4. Convert and audit reference data

After confirming that every accepted DFT calculation completed successfully,
use the converter matching the DFT code:

```bash
# VASP results
gpumdkit.sh -out2xyz <scf_results_directory>

# CP2K results (interactive converter)
gpumdkit.sh -cp2k2xyz

# Continue with the extxyz file produced by the selected converter
gpumdkit.sh -range <new_reference.xyz> energy
gpumdkit.sh -range <new_reference.xyz> force
gpumdkit.sh -min_dist_pbc <new_reference.xyz>
gpumdkit.sh -analyze_comp <new_reference.xyz>
```

Validate energy/force/virial labels, units, cells, species/type order,
composition coverage, duplicates, extreme structures, and failed or incomplete
calculations. Preserve the raw DFT outputs and record all exclusions.

## 5. Update the dataset and validate the model

Keep the previous dataset revision before merging new reference structures.
Record which configurations were added, rejected, weighted, or assigned to the
training/test sets. Prevent closely related trajectory frames from leaking
across an independence-sensitive split.

Prepare or update `nep.in` only after its architecture, loss, optimizer,
checkpoint, ZBL, and other nontrivial choices are approved. Run NEP separately,
then inspect:

- training and test losses;
- energy, force, and virial/stress parity outputs;
- errors by species, composition, phase, and configuration class;
- outliers, short-range behavior, and target-domain MD stability;
- the project-specific acceptance criteria.

Aggregate RMSE alone is not evidence that a model is ready for production.

## Iteration record

For reproducibility, retain the model/data revisions, approved decisions, exact
commands, working directories, executable versions, exit status, warnings,
generated files, rejected structures, validation results, and known
limitations. Stop the iteration whenever required outputs are missing or a
validation issue remains unresolved.

---

Thank you for using `GPUMDkit`! For questions or feedback, open an issue on the
GPUMDkit repository or contact Zihan YAN (yanzihan@westlake.edu.cn).
