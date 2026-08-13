# Batch and Manual Active-Learning Workflows

Use this reference for SCF/MD batch preparation and for manually reviewed NEP
active-learning iterations. GPUMDkit prepares directories and data; it does not
choose scientific settings or submit an active-learning cycle automatically.

## Contents

- Available workflows
- SCF batch preparation
- MD sample batch preparation
- Manual active-learning protocol
- Validation and handoff record

## Available Workflows

| Workflow | Menu | Description |
|---|---|---|
| SCF Batch (VASP) | 301 -> 1 | Prepare VASP single-point directories |
| SCF Batch (CP2K) | 301 -> 2 | Prepare CP2K single-point directories |
| MD Batch (GPUMD) | 302 | Prepare GPUMD MD sampling directories |
| MD Batch (LAMMPS) | 303 | Prepare LAMMPS MD sampling directories |

Run these tools in a dedicated working copy. The VASP/XYZ branches can move
input structures into generated `struct_fp/` or `struct_md/` directories.
Preserve the source structures elsewhere before proceeding.

## SCF Batch Preparation

### VASP SCF Batch

1. Put the structures to process in a dedicated working directory. The tool
   accepts `.vasp` files or one selected extxyz file; `.vasp` files take
   precedence when both formats are present.
2. Run the menu and provide a directory prefix when prompted:

```bash
gpumdkit.sh  # Select: 3) Workflow -> 301 -> 1) VASP
```

3. After preparation, place the user-approved `INCAR`, `POTCAR`, and `KPOINTS`
   files in the generated `fp/` directory.

```text
struct_fp/POSCAR_1.vasp, POSCAR_2.vasp, ...
fp/INCAR, POTCAR, KPOINTS                 (user-provided after preparation)
<prefix>_1/POSCAR -> ../struct_fp/POSCAR_1.vasp
<prefix>_1/INCAR  -> ../fp/INCAR
<prefix>_2/...
presub.sh
```

`presub.sh` is a template, not an authorization to run VASP. Review the
executable, resources, scheduler integration, and every DFT input before use.

### CP2K SCF Batch

Prefer the supported menu entry:

```bash
gpumdkit.sh  # Select: 3) Workflow -> 301 -> 2) CP2K
```

The prompt expects:

```text
<extxyz_file> <template.inp> <prefix_name>
```

The CP2K template must read coordinates from `pos.xyz`. The tool writes
`<prefix>_<index>/input.inp` and `pos.xyz`. It prepares files only and does not
submit CP2K.

## MD Sample Batch Preparation

### GPUMD MD Batch

1. Put `.vasp` structures or one selected extxyz file in a dedicated working
   directory.
2. Run:

```bash
gpumdkit.sh  # Select: 3) Workflow -> 302
```

3. After preparation, place the approved `nep.txt` and one matching
   `run_<index>.in` per sample in the generated `md/` directory.

```text
struct_md/model_1.xyz, model_2.xyz, ...
md/nep.txt, run_1.in, run_2.in, ...       (user-provided after preparation)
sample_1/model.xyz -> ../struct_md/model_1.xyz
sample_1/nep.txt   -> ../md/nep.txt
sample_1/run.in    -> ../md/run_1.in
presub.sh
```

### LAMMPS MD Batch

Run:

```bash
gpumdkit.sh  # Select: 3) Workflow -> 303
```

After preparation, place the approved `lmprun.in` and `nep.txt` in the
generated `md/` directory. The tool writes `struct_md/lammps_<index>.data`,
`sample_<index>/`, and `presub.sh`.

For both MD workflows, inspect every symlink and input mapping before execution.
Do not run the generated template until the user explicitly authorizes the MD
calculation.

## Manual Active-Learning Protocol

Use a staged, manually reviewed iteration so system-specific scientific and
scheduler choices remain explicit. The recommended route is ordinary MD
sampling followed by geometry filtering and NepTrain FPS selection.

### Resolve the iteration before execution

Ask the user to provide or approve:

- the model revision, target compositions/phases, and intended deployment domain;
- the candidate-generation structures and GPUMD conditions;
- the MD trajectory, distance/box filters, NepTrain FPS rule, and structure-count criteria;
- the DFT code, settings, pseudopotentials, convergence criteria, and job budget;
- dataset merge/split rules and leakage controls;
- NEP training configuration and model-acceptance criteria.

Do not replace any unresolved item with the examples or defaults from another
system.

### Recommended path: MD trajectory and NepTrain FPS

Run authorized ordinary MD across the user-approved starting structures and
conditions. Inspect stability first, then apply only user-approved filters:

```bash
gpumdkit.sh -min_dist_pbc dump.xyz
gpumdkit.sh -filter_dist_pbc dump.xyz <minimum_distance>
gpumdkit.sh -filter_box filtered_dump.xyz <maximum_box_edge>
gpumdkit.sh  # Select: 2) Sample Structures -> 203) FPS by NepTrain
```

`-min_dist_pbc` reports distances but does not filter. Use
`-filter_dist_pbc` when structures must actually be removed. The box filter is
optional; skip it unless a box-edge criterion is meaningful and approved for
the dataset. Record every rejected structure and filter threshold. NepTrain FPS
then selects structures that add descriptor-space diversity relative to the
current training set. Record the FPS stopping rule and selected frame mapping.

### Optional alternative: committee uncertainty

Use committee uncertainty only when the project explicitly chooses it. Read
`gpumd-outputs.md`: GPUMD `active` requires compatible NEP committee models and
writes `active.out` plus `active.xyz`. Confirm the model order, check interval,
output flags, threshold, `top_n`, and `min_deviation`; after an authorized run,
validate the logs and geometry before using menu 205. Keep `active.out` and
`active.xyz` from the same run.

### Reference calculation and dataset update

1. Use menu 301 to prepare DFT directories from the approved `selected.xyz`.
2. Review the generated structures, DFT inputs, resource request, and job count.
3. Run DFT only after separate explicit authorization and through the user's
   scheduler procedure.
4. Confirm every calculation completed and exclude failed/incomplete outputs
   with recorded reasons.
5. Convert the reference outputs with the route matching the DFT code, then
   inspect the converter's generated dataset:

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

6. Load `nep-data.md`; validate labels, units, cells, species/type order,
   duplicates, train/test leakage, and configuration coverage before merging.
7. Preserve the previous dataset and record exactly which structures were added,
   rejected, weighted, or reassigned.
8. Load `nep.md`, `nep-parameters.md` when configuration changes, and
   `nep-outputs.md`. Train only after authorization, then evaluate train/test
   behavior, outliers, target-domain stability, and the user-approved acceptance
   criteria. Do not approve a model from aggregate RMSE alone.

## Validation and Handoff Record

For each iteration, report at least:

| Field | Record |
|---|---|
| Inputs | Model hash/revision, structures, GPUMD and DFT inputs |
| Decisions | Candidate path and every approved scientific threshold |
| Authorization | Which MD, DFT, and NEP runs were explicitly approved |
| Commands | Exact commands, working directories, versions, and exit status |
| Outputs | Candidate, selected, reference, dataset, and model files |
| Validation | Parser/log checks, geometry/data audits, errors, and stability tests |
| Limitations | Failed jobs, exclusions, uncovered regimes, and unresolved risks |

Do not advance to the next stage when required outputs are missing, logs contain
unresolved errors, structures are unphysical, or the validation gate has not
been approved.

## Detailed Documentation

- `${GPUMDkit_path}/docs/tutorials/en/workflow_scripts.md`
- `${GPUMDkit_path}/docs/tutorials/en/active_learning_workflow.md`
