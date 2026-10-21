# NEP Model Modifier

> **English** | [简体中文](README_zh-CN.md)

The NEP model modifier is an interactive GPUMDkit front end to the model
modification functions provided by
[calorine](https://calorine.materialsmodeling.org/dev/get_started/modifying_models.html).
It turns the Python methods `augment()`, `prune()`, `add_species()`,
`remove_species()`, and `keep_species()` into a guided terminal workflow with
input validation, architecture inspection, change review, and package export.

This feature is under testing. Validate the generated model, restart, and
`nep.in` together before starting further training.

## When to use this tool

A trained model does not always need to be rebuilt from scratch when its scope
changes. This tool is useful when you want to:

- increase the capacity of a model that appears too small for the target data;
- introduce higher-body descriptors or a charge output head;
- make a smaller model as a starting point for speed-oriented retraining;
- extract one element or a small element subset from a multi-element model;
- extend an existing model to additional chemical species while retaining the
  parameters learned for its original species.

Modification is only a starting point. Adding or removing parameters changes
the optimization problem, and adding a species creates untrained parameters for
the new element and its descriptor pairs. Retrain and validate the resulting
model against data representative of the intended application.

## Requirements

- `calorine >= 3.4`
- A NEP4 model in `nep.txt` format
- A matching `nep.restart` for model expansion, reduction, or adding species
- An optional source `nep.in` whose training-specific settings should be kept

Install calorine with:

```bash
pip install 'calorine>=3.4'
```

## Usage

Start with default files in the current directory:

```bash
gpumdkit.sh -nep_modifier
```

Pass the model package explicitly:

```bash
gpumdkit.sh -nep_modifier path/nep.txt path/nep.restart path/nep.in
```

Use `-` to skip the restart or training input:

```bash
gpumdkit.sh -nep_modifier path/nep.txt - path/nep.in
gpumdkit.sh -nep_modifier path/nep.txt path/nep.restart -
```

Show help without importing calorine:

```bash
gpumdkit.sh -nep_modifier -h
```

When paths are omitted, the program prompts separately for `nep.txt`,
`nep.restart`, and `nep.in`. Defaults are resolved beside the selected model,
not from an unrelated working directory.

## Main menu

```text
+-------------------------------------------------------------------------------+
|                              NEP MODEL MODIFIER                               |
+-------------------------------------------------------------------------------+
| 1) Expand model capacity             4) Remove chemical species               |
| 2) Reduce model capacity             5) Keep selected species                 |
| 3) Add chemical species              6) Inspect current model                 |
| 7) Review pending changes            8) Export model package                  |
| 0) Exit                                                                       |
+-------------------------------------------------------------------------------+
```

The program warns before exiting if accepted changes have not been exported.

## Operations

| Operation | Calorine method | Restart needed | Further training |
|---|---|---:|---:|
| Expand capacity | `augment()` | Yes | Required |
| Reduce capacity | `prune()` | Yes | Required |
| Add species | `add_species()` | Yes | Required |
| Remove species | `remove_species()` | No | Optional |
| Keep selected species | `keep_species()` | No | Optional |

### Expand model capacity

Multiple changes can be selected and passed to one atomic `augment()` call:

- increase `n_neuron`;
- increase `l_max_4b` or `l_max_5b`;
- enable `q_112`, `q_123`, `q_233`, or `q_134` descriptor terms;
- add a charge output head.

When the installed calorine exposes `charge_mode`, the user must select mode 1
or 2 explicitly. The program shows calorine's documented SNES defaults and asks
the user whether to accept them or enter `sigma_new`, `sigma_factor`, and
`sigma_floor` manually.

Existing trained means (`mu`) are retained. New parameters start inactive or
newly initialized according to calorine, while the associated SNES exploration
widths are initialized from the selected sigma settings. This lets subsequent
training adapt the enlarged architecture without discarding the learned part of
the model.

### Reduce model capacity

One `prune()` operation may combine:

- reducing `n_neuron`;
- reducing or removing `l_max_4b` and `l_max_5b` terms;
- disabling enabled higher-body descriptor terms;
- removing a charge output head from a charge-aware model.

Calorine ranks neurons by their model weights when reducing the neuron count.
The reduced model must be trained and validated before production use.
Turning a descriptor family fully off removes its stored descriptor columns;
merely lowering a non-zero `l_max_4b` can change the header without reducing the
stored descriptor dimension.

### Add species

`add_species()` is available only for NEP4 models with restart data. The program:

- rejects duplicate or existing species;
- requires a user-provided random seed;
- requests radial and angular cutoffs when the source model uses typewise
  cutoffs;
- records the seed, cutoff values, and accepted SNES settings.

New species parameters are not trained, so the exported model is not ready for
production until it has been trained and validated on appropriate reference
data.

### Remove or keep species

These operations can be performed without restart data. If restart data are
loaded, the program shows the documented adaptive-sigma settings before the
operation. At least one species must remain.

`remove_species()` is convenient when the elements to discard are few.
`keep_species()` is the inverse interface and is usually easier for extracting
a small subset from a large foundation model. Both remove the species-specific
ANN subnetworks and descriptor-weight pairs associated with discarded species.

## Guided examples

The examples below show the important responses rather than every status line.
The exact current values depend on the loaded model.

### Increase neurons and enable a descriptor term

```text
Input the function number:
------------>>
1
Input one or more choices, separated by spaces:
------------>>
1 4
Input target neuron count (current: 50):
------------>>
60
Use these SNES initialization defaults? (Y/n)
------------>>
y
Apply these changes? (y/N)
------------>>
y
```

This performs one combined `augment()` operation: the neuron count becomes 60
and `q_112` is enabled. Export the package and use the generated `.in`, `.txt`,
and `.restart` together for continued training.

### Extract a Li-O submodel

```text
Input the function number:
------------>>
5
Input species to keep, separated by spaces:
------------>>
Li O
Keep only these species? (y/N)
------------>>
y
```

This is easier than listing every other species for removal. The species order
follows the retained order returned by calorine; inspect it before preparing
new training data or a GPUMD `model.xyz`.

### Add carbon reproducibly

```text
Input the function number:
------------>>
3
Input species to add, separated by spaces:
------------>>
C
Input random seed for reproducible initialization:
------------>>
20260712
Use these SNES initialization defaults? (Y/n)
------------>>
y
Add these species? (y/N)
------------>>
y
```

For a typewise-cutoff model, radial and angular cutoffs for carbon are requested
before the seed. These are scientific model choices and should come from the
intended training setup, not from this example. The selected seed and cutoffs
are written to the change summary.

### Review and export

Choose `7` to inspect every accepted operation, then choose `8` and provide an
output directory. A successful export prints all generated paths. If pending
changes remain when `0` is selected, the tool asks whether to leave without
exporting them.

## Model inspection and change history

The inspection page reports:

- NEP and model type;
- species order;
- radial and angular cutoffs;
- `n_max`, `basis_size`, and `l_max`;
- higher-body descriptor switches;
- descriptor, ANN, neuron, and total parameter counts;
- restart availability, ZBL, and charge mode when available.

`Review pending changes` shows accepted operations, arguments, changed model
fields, and whether each operation has already been exported.

## Safe export

The tool never overwrites the source package. It chooses a suffix that is free
for every output file, writes all files in a temporary directory, verifies that
they exist, and then moves them into the requested output directory.

For an input named `nep.txt`, the first package is:

| File | Purpose |
|---|---|
| `nep_modified.txt` | Modified model |
| `nep_modified.restart` | Matching restart, when loaded |
| `nep_modified.in` | Updated architecture plus retained training settings |
| `nep_modified.changes.txt` | Source paths, calorine version, parameters, seed, and operation history |

If any member already exists, the complete package uses the next common suffix,
such as `nep_modified_1.*`.

When a source `nep.in` is supplied, its training settings are preserved and the
architecture fields are replaced with `model.training_parameters`. Without a
source `nep.in`, the generated file contains architecture fields only.

The generated `nep.in` must be checked with the NEP executable version that will
perform the next training run. Do not combine the exported model or restart with
an older, incompatible `nep.in`.

## Error handling

- Missing restart data disables operations that require it without terminating
  the session.
- Invalid numbers, choices, species, and confirmations are reported without a
  Python traceback.
- EOF and `Ctrl-C` exit cleanly.
- Loading or export failures are shown explicitly; they are not reported as
  successful partial exports.

## Citation

This feature directly uses calorine's NEP model modification implementation.
If you use it in research, cite calorine as requested by the calorine project;
also cite GPUMDkit when appropriate:

- E. Lindgren, J. M. Rahm, E. Fransson, F. Eriksson, N. Österbacka, Z. Fan,
  and P. Erhart, “calorine: A Python package for constructing and sampling
  neuroevolution potential models,” *Journal of Open Source Software* **9**(95),
  6264 (2024), <https://doi.org/10.21105/joss.06264>.
- Z. Yan et al., *MGE Advances* **4**, e70074 (2026),
  <https://doi.org/10.1002/mgea.70074>.
