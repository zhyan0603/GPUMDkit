<div align="center">
  <h1>🔁 Active Learning Workflow</h1>
  <p style="text-align: justify;">One round of NEP data iteration with GPUMDkit: sample with MD, select structures, label them with DFT, and retrain.</p>
</div>

## Overview

One iteration has five steps: generate candidate structures with MD, select the informative ones, prepare DFT calculations, convert the results, and update the dataset and model. GPUMDkit prepares inputs and analyzes outputs; you run GPUMD, DFT, and NEP yourself (or through your scheduler).

## Before you start

Work in a new directory and keep the previous dataset and `nep.txt`. Decide these points up front, since they shape every later step:

| Item | What to decide |
|---|---|
| Model | Which `nep.txt` to start from and which systems or conditions it should cover |
| MD sampling | Starting structures and GPUMD conditions |
| Selection | Distance and box-size limits, the FPS stopping rule, and how many structures to keep |
| DFT | Code, functional, pseudopotentials, convergence settings, and compute budget |
| Dataset | How new data is merged, weighted, and split into training and test sets |
| Acceptance | Which accuracy and MD-stability results make the new model good enough |

## 1. Generate candidate structures

With several starting structures or `run.in` files, prepare the GPUMD batch directories with function 302:

~~~bash
gpumdkit.sh  # Select: 3) Workflow -> 302
~~~

Function 302 only creates directories and symlinks. Afterwards, put `nep.txt` and one `run_<index>.in` per structure into the generated `md/` directory, check that each directory points to the right files, and run GPUMD.

GPUMD appends to some output files, so rerun in a clean directory or archive old outputs first. Before selecting structures, check that each run finished normally, that temperature, energy, pressure, and cell behave as expected, and that `dump.xyz` is complete and free of obviously unphysical structures.

## 2. Select candidate structures

Check distances first, then filter and run farthest-point sampling (FPS):

~~~bash
gpumdkit.sh -min_dist_pbc dump.xyz                            # report only
gpumdkit.sh -filter_dist_pbc dump.xyz <minimum_distance>      # writes filtered_dump.xyz
gpumdkit.sh -filter_box filtered_dump.xyz <maximum_box_edge>  # optional, writes filtered_by_box.xyz
gpumdkit.sh  # Select: 2) Sample Structures -> 203) FPS by NepTrain
~~~

`-min_dist_pbc` only reports distances; `-filter_dist_pbc` removes structures. Use the box filter only if a box-size limit applies to your system; in that case give `filtered_by_box.xyz` to function 203. Note the thresholds and FPS stopping rule you used, and which source frames were selected.

### Alternative: committee uncertainty

If you train a committee of NEP models, GPUMD's `active` keyword writes `active.out` and `active.xyz` during MD. Function 205 then ranks the structures by force deviation and can cap their number. Check the model order, check interval, threshold, `top_n`, and `min_deviation` before using it.

## 3. Prepare the DFT calculations

Pass the selected structures (`selected.xyz`) to function 301:

~~~bash
gpumdkit.sh  # Select: 3) Workflow -> 301 -> VASP or CP2K
~~~

The VASP branch creates `struct_fp/`, `fp/`, one directory per structure, and `presub.sh`; then put your `INCAR`, `POTCAR`, and `KPOINTS` into `fp/`. The CP2K branch (301 -> 2) asks for an extxyz file, an input template, and a directory prefix.

Before submitting, check the number of structures, atom and species order, cells, DFT inputs, and the job script. GPUMDkit does not submit jobs.

## 4. Convert and check the reference data

Once the DFT calculations have finished, convert them.

VASP (writes `NEPdataset/train.xyz` in the current directory):

~~~bash
gpumdkit.sh -out2xyz <scf_results_directory>
~~~

If the OUTCAR files contain complete `BORN EFFECTIVE CHARGES` blocks and you need per-atom BEC labels, use `gpumdkit.sh -out2xyz_bec <scf_results_directory>` instead.

CP2K: run the command from the root of the CP2K results. It scans `.log` files and their `.xyz`/`.inp` files recursively and writes `cp2k_exyz.xyz`, with a processing summary in `Logfile.txt`:

~~~bash
gpumdkit.sh -cp2k2xyz
~~~

Then check the new file (replace `<new_reference.xyz>` with the file written above):

~~~bash
gpumdkit.sh -range <new_reference.xyz> energy
gpumdkit.sh -range <new_reference.xyz> force
gpumdkit.sh -min_dist_pbc <new_reference.xyz>
gpumdkit.sh -analyze_comp <new_reference.xyz>
~~~

Look for wrong units or labels, unexpected species order, missing compositions, duplicates, extreme structures, and failed calculations. Keep the raw DFT outputs. See [Format Conversion](format_conversion.md) for output locations and when files are overwritten.

## 5. Update the dataset and retrain

Keep a copy of the previous dataset, then merge the new structures. If neighboring trajectory frames end up in both the training and test sets, the test error will look better than it is, so split by trajectory or time when independence matters.

Update `nep.in`, train NEP, and check:

- training and test losses;
- energy, force, and virial/stress parity plots (`gpumdkit.sh -plt train`, `-plt test`);
- errors by species, composition, phase, and structure type;
- outliers, short-range behavior, and MD stability in the target conditions.

A low overall RMSE alone does not show that the model is ready. Once it meets your acceptance criteria, use it for the next round of MD.

## Keep a record

For each iteration, keep the model and dataset versions, the commands and thresholds used, rejected structures, and validation results. Stop and investigate if outputs are missing or a validation problem is unresolved.
