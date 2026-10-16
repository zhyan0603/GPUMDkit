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
   requires one batch of `.vasp` files or one multi-frame extxyz file named
   `*.xyz` for this protocol. Do not mix formats or leave multiple structure
   XYZ files. The current script prioritizes `.vasp` if mixed; that fallback
   is not the recommended preparation procedure.
2. Run the menu and provide a directory prefix when prompted:

```bash
gpumdkit.sh  # Select: 3) Workflow -> 301 -> 1) VASP
```

3. After preparation, place the user-approved VASP inputs in `fp/`. Extxyz
   structures with one element set keep the flat layout and share `fp/POTCAR`.
   If multiple element sets are detected, structures go under
   `struct_fp/<elements>/`, and each calculation directory links to a matching
   `fp/POTCAR_<elements>` file. Different atom ratios with the same elements
   share one group. Detection uses source symbols; the converter preserves the
   trajectory-wide first-seen element order and verifies each written POSCAR's
   order and counts by reading it back. All groups share `fp/INCAR`; POTCAR
   order must match the POSCAR header. Every calculation directory links to
   `../fp/KPOINTS`; provide the target for explicit k-point input. If it is
   absent, the link stays dangling and VASP uses `KSPACING` from `INCAR`.

```text
struct_fp/Al_O/POSCAR_1.vasp, ...          (only when element sets differ)
fp/INCAR, POTCAR_Al_O, ...                 (user-provided after preparation)
<prefix>_1/POSCAR -> ../struct_fp/Al_O/POSCAR_1.vasp
<prefix>_1/POTCAR -> ../fp/POTCAR_Al_O
<prefix>_1/KPOINTS -> ../fp/KPOINTS
<prefix>_1/INCAR  -> ../fp/INCAR
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

### Directory convention and stage handoffs

Use `iter00`, `iter01`, etc.; do not use `round_00`. Keep each stage's inputs,
outputs, logs, frame mapping, and checks inside that stage. The following is a
workflow convention, not a directory tree automatically created by GPUMDkit:

```text
active_learning/
+-- iter00/
|   +-- perturbation/
|   +-- nep89md/
|   +-- 01.fps/
|   +-- 02.scf/
|   +-- 03.neptrain/
|   +-- 04.diffusion/
+-- iter01/
|   +-- 00.md/
|   +-- 01.fps/
|   +-- 02.scf/
|   +-- 03.neptrain/
|   +-- 04.diffusion/
+-- iter02/
    +-- ...
```

`iter00` prepares initial data with perturbations, NEP89 MD, or both, as chosen
by the user. Later iterations normally start from `00.md` using the accepted
model from the previous iteration's `03.neptrain/nep.txt`. `04.xx` is a
placeholder: replace `xx` with the actual validation task, such as `diffusion`;
number additional property stages consecutively (`05.elastic`, etc.). Create
only stages needed for the project. Tool-generated `md/`, `fp/`, `struct_md/`,
`struct_fp/`, and calculation directories remain inside their owning stage;
do not make parallel candidate/conversion/training directories with overlapping
purposes. Do not move an existing project into this layout without authorization.

The names below describe handoffs. Preserve actual tool output names, or record
any copy/rename; they are not promises of automatic output naming.

| Stage | Inputs and model source | Outputs and next consumer |
|---|---|---|
| `iter00/perturbation` | Approved initial structures and perturbation settings; no model required | `POSCAR_*.vasp`, `perturb.txt`; candidate geometries to `01.fps` |
| `iter00/nep89md` | Approved initial structures, NEP89 model/revision and MD inputs | Raw trajectories and MD logs; candidate geometries to `01.fps` |
| `iterNN/00.md` (`NN >= 01`) | Starting structures, accepted `iter(NN-1)/03.neptrain/nep.txt`, approved MD inputs | Raw trajectories and MD logs; candidate geometries to `01.fps` |
| `01.fps` | Upstream geometries, approved descriptor model and baseline structures | Raw combined candidates, filtered candidates, `selected.xyz`, FPS outputs and source-frame mapping; selected geometries to `02.scf` |
| `02.scf` | Working copy of selected geometries and approved DFT settings | Calculation inputs/results, completion/convergence audit and converted DFT reference data (e.g. `new_reference.xyz`); accepted labels to `03.neptrain` |
| `03.neptrain` | Accepted new DFT labels, previous dataset (if any), merge decisions and `nep.in` | Final cumulative `train.xyz`, optional user-approved `test.xyz`, `nep.txt`, training logs and dataset membership record; model to validation and the next iteration |
| `04.<property>` | Current model and approved property protocol/reference | Property results, checks and acceptance decision; feedback for the next iteration |

### Decisions and authorization

Directory naming, copying without overwriting, and routine organization within
an authorized task may proceed directly. Reuse confirmed project settings;
for missing scientific choices, present concrete suggestions and their reasons
for approval before writing executable inputs. Resolve:

- deployment compositions/phases, structures, NEP89 applicability or previous model;
- perturbation count, atomic displacement (Angstrom), cell perturbation fraction
  and distribution (`normal`, `uniform`, or `const`);
- MD temperature, time step, duration, sampling interval and ensemble; for NPT,
  pressure, `iso/aniso/tri`, thermostat/barostat coupling periods and units;
- geometry filters, descriptor model, FPS baseline and count/stopping rule;
- DFT functional, pseudopotentials, cutoffs, k-point treatment, electronic
  convergence, labels and job budget;
- dataset inclusion, weighting, split choice, NEP architecture (version,
  cutoffs, basis, orders and neurons), loss weights, ZBL and training generations;
- validation protocols, acceptance and iteration stopping criteria.

After a user changes a parameter, update affected inputs and the stage record,
mark results from old settings as superseded, and reassess dependent stages.
Never carry the old plan forward silently. Material-specific values and cited
papers' settings are examples, not universal defaults.

Preparing files does not authorize MD, DFT, NEP training or scheduler submission.
If the user asks only for preparation, stop after input checks. If execution of
a stage is already authorized, complete that stage within its approved scope
without asking for the same authorization again. Ask only for missing decisions
or changes outside that scope. Resume from the last completed and checked stage;
do not overwrite results or rerun calculations just to recreate the layout.

### Stage 0: initial candidates or subsequent MD

For perturbations, use menu `2 -> 204` and load `sampling.md`; explicitly record
`<input.vasp> <pert_num> <cell_pert> <atom_pert> <style>`. Keep generated
`POSCAR_*.vasp` and `perturb.txt` in `perturbation/`. The report describes changes;
it does not decide which geometries are acceptable. Convert geometries for the
FPS input through a supported route in `format-conversion.md`, inside `01.fps`.

For NEP89 or later MD, load the relevant GPUMD input/ensemble/output references,
verify species coverage, model provenance and the installed version, then
prepare `model.xyz`, `run.in`, and the model in `nep89md/` or `00.md/`. Menu 302
is available for batches; its generated directories stay inside that stage.
Run only if authorized. Keep raw trajectories and logs; inspect stability and
separate equilibration from sampling before handing off geometries.

NEP89 and later NEP trajectory energies/forces are predictions, not DFT labels.
Keep them out of reference data and final training inputs; replacing labels
requires actual, checked DFT results for the corresponding configurations.

### Stage 1: geometry audit and NepTrain FPS

Work in `01.fps`. Combine the chosen upstream candidates with source IDs (source
file, run, frame, iteration); preserve raw candidates and record filters and
rejections separately. Apply only approved, physically meaningful thresholds:

```bash
gpumdkit.sh -min_dist_pbc <candidates.xyz>
gpumdkit.sh -filter_dist_pbc <candidates.xyz> <minimum_distance>
# Optional: only with an approved box criterion, using the actual filtered file
gpumdkit.sh -filter_box <filtered_candidates.xyz> <maximum_box_edge>
gpumdkit.sh  # Select: 2) Sample Structures -> 203) FPS by NepTrain
```

`-min_dist_pbc` reports distances; `-filter_dist_pbc` removes structures. The box
filter writes `filtered_by_box.xyz`. Feed the actual retained file to FPS.
Load `sampling.md` for the prompts: candidates, baseline structures, `nep.txt`,
and a distance threshold or structure-count bounds. Prefer NepTrain. Identify
the descriptor model's iteration/revision and calculate candidate and baseline
descriptors with the same model. Later baselines normally use the previous
cumulative training geometries; `iter00` needs an explicitly chosen seed/baseline
file because the current tool requires one. Do not invent DFT labels for it.

The current implementation uses mean atomic NEP descriptors per structure.
Record threshold or count bounds and actual selection count; descriptor distances
are model-dependent, not a calibrated energy/force error. Keep `selected.xyz`,
`select.png`, `pca_sample.txt`, `pca_train.txt`, `pca_selected.txt` in `01.fps`.
Maintain an explicit selected-frame mapping: PCA files alone are not provenance.
Changing descriptor models requires recalculating both descriptor sets and
reviewing the distance scale. PCA plots are diagnostics, not proof of coverage.

### Stage 2: DFT preparation, labeling and conversion

1. Copy only selected geometries to `02.scf`. Before running **301**, its working
   directory must contain either one batch of `.vasp` structure files or one
   multi-frame extxyz file named `*.xyz`, never both and never multiple structure
   XYZ files. This protocol is stricter than the current VASP script, which
   prioritizes `.vasp` when formats are mixed; do not rely on that fallback.
2. Use `gpumdkit.sh`, menu `3 -> 301`, then the appropriate DFT branch. For VASP,
   generate calculation directories, `struct_fp/` and `fp/` first; only then
   supply approved `INCAR`, `POTCAR` (or grouped POTCARs) and other inputs in
   `fp/`. With approved `KSPACING` in `INCAR`, `KPOINTS` is not required; check
   that no unintended existing KPOINTS overrides it. For CP2K, use the required
   approved template at the 301 prompt; see the CP2K section above.
3. Check input/species order, links, selected-frame-to-calculation mapping,
   resources and job count. Run DFT only within the user's execution authorization.
4. Audit each selected configuration: expected task, completion, electronic
   convergence, finite labels and matching geometry. Report missing, failed or
   unconverged tasks explicitly with their IDs. Do not silently treat a partial
   batch as final; resolve retries or user-approved exclusions before finalizing.
5. Convert only accepted outputs using the corresponding route, inside `02.scf`:

```bash
gpumdkit.sh -out2xyz <accepted_scf_results_directory>
gpumdkit.sh -cp2k2xyz  # CP2K alternative: run in the results root; scans logs recursively, writes cp2k_exyz.xyz
# Inspect the actual extxyz output from the chosen converter
gpumdkit.sh -range <new_reference.xyz> energy
gpumdkit.sh -range <new_reference.xyz> force
gpumdkit.sh -min_dist_pbc <new_reference.xyz>
gpumdkit.sh -analyze_comp <new_reference.xyz>
```

Load `nep-data.md` and `format-conversion.md`. Confirm energy in eV (total per
configuration), force in eV/Angstrom, and virial convention/units when present;
check cells, PBC, species/counts, finite labels and one-to-one configuration
mapping. File existence or successful conversion does not establish convergence.
Raw candidates, selected geometries and accepted DFT reference labels are distinct
products; only the last supplies new training labels.

### Stage 3: cumulative dataset and NEP training

Work in `03.neptrain`. Merge accepted new references with the previous checked
dataset without modifying that previous dataset. Audit duplicates and reference
setting consistency; record source frame, DFT task and inclusion/weight/split
for each final configuration. `train.xyz` is the final cumulative training
input, not an unreviewed concatenation of trajectories or converter outputs.

Do not create a test set by default. Ask whether a split is wanted and how it
should be made. If approved, keep related trajectory frames/parent structures
from leaking across splits and record membership. If no test set is requested,
omit `test.xyz`; never copy training data as testing data. Verify optional-test
behavior for the installed NEP version and report that independent held-out
error has not been measured.

Load `nep.md`, `nep-parameters.md` and `nep-outputs.md`; write only confirmed
`nep.in` settings. After authorized training, inspect losses, energy/force/virial
errors where applicable, outliers and per-species/environment errors. Preserve
model provenance and training outputs in this stage; pass the current model to
property validation, not directly to acceptance based on aggregate RMSE.

### Stage 4 and later: model validation and the next iteration

Name each directory for its actual property (`04.diffusion`, etc.). Use approved
protocols and references to check target-domain dynamics, geometry, stability,
and relevant properties; examples include equation of state, elastic response,
relative phase/defect energies, migration barriers or diffusion. Choose only
checks relevant to the intended use; literature does not supply universal
thresholds or authorize expensive validation runs.

Report training fit, independent validation (if performed), property agreement,
and uncovered regimes separately. A training fit or stable MD alone does not
establish transferability. Record the accepted model revision and pass it with
the checked cumulative dataset to the next iteration. If criteria fail, retain
the evidence and propose targeted candidates. Stop only against the user's
agreed coverage, error/property and budget criteria; small FPS distances alone
are insufficient.

### Optional selection route: committee uncertainty

Use only when chosen by the user. Load `gpumd-outputs.md`: GPUMD `active` writes
`active.out` and `active.xyz` using compatible NEP committee models. Confirm
model order, check interval, output flags, threshold, `top_n` and `min_deviation`;
after an authorized run, check logs and geometry before menu 205. Keep paired
outputs in the sampling stage. Committee disagreement is a selection signal,
not a guaranteed error bound; correlated models can share blind spots.

### Method guidance and literature

These sources motivate methods; the project must approve their application.

- **Coverage:** descriptor diversity complements MD sampling. Check composition,
  phase and local-environment coverage as well as trajectory size. See
  [Training data selection for accuracy and transferability of interatomic potentials](https://doi.org/10.1038/s41524-022-00872-x).
- **Rare environments:** ordinary MD can undersample defects and transition
  states. Propose targeted starting structures or enhanced sampling when needed;
  these require approved protocols. See
  [Hyperactive learning for data-driven interatomic potentials](https://doi.org/10.1038/s41524-023-01104-6).
- **Validation:** low mean energy/force errors can coexist with errors in dynamics,
  defects and rare events. Inspect relevant local errors and properties. See
  [Discrepancies and error evaluation metrics for machine learning interatomic potentials](https://doi.org/10.1038/s41524-023-01123-3).

For mean-descriptor FPS, separately inspect minority species, defects and unusual
local environments: structure averaging can hide them. This is a practical
limitation of the current representation, not evidence that every such candidate
needs labeling. Selection, reference labeling, retraining and validation remain
separate stages.

## Validation and Handoff Record

Keep a concise record in each owning stage; do not add a parallel QA directory.
Record input paths/model revision, confirmed parameters and changes,
authorization scope, command/cwd/version/status, output paths and frame mapping,
checks and acceptance, failures/exclusions and remaining limits. On resumption,
verify these records against the actual files and identify the next incomplete
stage. Missing required outputs or unresolved errors block dependent work;
report them rather than declaring the iteration complete.

## Detailed Documentation

- `${GPUMDkit_path}/docs/tutorials/en/workflow_scripts.md`
- `${GPUMDkit_path}/docs/tutorials/en/active_learning_workflow.md`
