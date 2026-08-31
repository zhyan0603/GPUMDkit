<div align="center">
  <h1>🔧 Workflow Scripts</h1>
  <p style="text-align: justify;">Workflow scripts prepare batch directories, input files, and links for DFT, MD, and NEP data work.</p>
</div>

The workflow scripts prepare directories and files. They do not submit calculations.

**Script Location:** `Scripts/workflow/`

## Interactive Entry

Use the interactive menu as the normal user entry point:

~~~bash
gpumdkit.sh
~~~

Choose:

~~~text
3) Workflow
~~~

The workflow menu looks like this (the text matches the program output):

~~~text
+---------------------------------------------------------+
|                      WORKFLOW TOOLS                     |
+---------------------------------------------------------+
| 301) SCF batch pretreatment                             |
| 302) MD sample batch pretreatment (gpumd)               |
| 303) MD sample batch pretreatment (lmp)                 |
+---------------------------------------------------------+
| 000) Return to the main menu                            |
+---------------------------------------------------------+
Input the function number:
~~~

| Menu | Workflow | Use it for |
|---|---|---|
| 301 | SCF batch pretreatment | Prepare DFT single-point directories for structures |
| 302 | MD sample batch pretreatment (gpumd) | Prepare multiple GPUMD sampling directories |
| 303 | MD sample batch pretreatment (lmp) | Prepare multiple LAMMPS sampling directories |

## 1. Prepare DFT single-point calculations {#scf-batch-pretreatment}

This section covers the input files and single-point directory preparation.

Work in a dedicated directory and keep a separate copy of the source structures.
The scripts create `struct_fp/`, `struct_md/`, calculation directories, and
`presub.sh` in the current directory.

### VASP SCF preparation (301 -> 1)

The input can be several `.vasp` files:

~~~text
current directory/
├── POSCAR_1.vasp
├── POSCAR_2.vasp
└── POSCAR_3.vasp
~~~

or one extxyz file:

~~~text
current directory/
└── sampled.xyz
~~~

The input selection order is:

1. If one or more `.vasp` files exist, process `.vasp` files first.
2. If both `.vasp` and `.xyz` files exist, print a notice and process only `.vasp`.
3. If no `.vasp` file exists and one `.xyz` file exists, convert that extxyz to POSCAR files.
4. If no `.vasp` file exists and multiple `.xyz` files exist, ask you to choose one.

Run the menu:

~~~bash
gpumdkit.sh
# choose 3 -> 301 -> VASP
~~~

The usual output is:

~~~text
current directory/
├── struct_fp/
│   ├── POSCAR_1.vasp
│   └── ...
├── fp/
├── <prefix>_1/
├── <prefix>_2/
└── presub.sh
~~~

The script creates `fp/`. After preparation, place the reviewed `INCAR`, `POTCAR`,
and `KPOINTS` files in `fp/`; the calculation directories link to them. The
`prefix` becomes the calculation-directory prefix. When `out2xyz` is used later,
the directory name is written as `config_type` in the extxyz output.

### CP2K SCF preparation (301 -> 2)

Prepare an extxyz file and a template. The template must read coordinates from
`pos.xyz`; the bundled example is `Scripts/workflow/cp2k_template.inp`. Run
through the menu:

~~~bash
gpumdkit.sh
# choose 3 -> 301 -> CP2K
~~~

The prompt is:

~~~text
<extxyz_file> <template.inp> <prefix_name>
~~~

For example:

~~~text
dump.xyz template.inp H2O
~~~

The script creates `<prefix>_<index>/` for each structure and writes `input.inp`
and `pos.xyz`.

Source users can also invoke the script directly. This relative script path
assumes the GPUMDkit checkout root as the working directory; replace the input
filenames with actual data paths. Alternatively, use the script's absolute path
from your data directory. Conda users should prefer the menu above.

~~~bash
python Scripts/workflow/scf_batch_pretreatment_cp2k.py train.xyz template.inp calc
~~~

## 2. Prepare MD sampling directories {#md-sampling-pretreatment}

### GPUMD (302)

Put `.vasp` files or one extxyz file in a dedicated working directory. The input
selection rules are the same as for VASP preparation: `.vasp` files take
precedence, and multiple `.xyz` files without `.vasp` require a choice.

~~~bash
gpumdkit.sh
# choose 3 -> 302
~~~

Typical inputs are:

~~~text
current directory/
├── POSCAR_1.vasp
└── POSCAR_2.vasp
~~~

or:

~~~text
current directory/
└── dump.xyz
~~~

The usual output is:

~~~text
current directory/
├── struct_md/
│   ├── model_1.xyz
│   └── ...
├── md/
├── sample_1/
├── sample_2/
└── presub.sh
~~~

After preparation, place the reviewed `nep.txt` and the matching `run_*.in` file
for each sample in `md/`. For example, `run_1.in` links to `sample_1/run.in`
and `run_2.in` links to `sample_2/run.in`.

### LAMMPS (303)

Put `.vasp` files or one extxyz file in a dedicated working directory. When both
formats are present, `.vasp` files take precedence.

~~~bash
gpumdkit.sh
# choose 3 -> 303
~~~

The script creates LAMMPS data files in `struct_md/`, sample directories, `md/`,
and `presub.sh`:

~~~text
current directory/
├── struct_md/
│   ├── lammps_1.data
│   └── ...
├── md/
├── sample_1/
├── sample_2/
└── presub.sh
~~~

After preparation, place the reviewed `lmprun.in` and `nep.txt` files in `md/`.

## 3. Confirm before running and inspect results

Before submitting through your own execution procedure, check:

- the structure count, atom and type order, cells, and generated-directory count;
- `INCAR`, `POTCAR`, and `KPOINTS` in `fp/`, or `nep.txt`, `run_*.in`, and
  `lmprun.in` in `md/`;
- every sample-directory link and file name;
- the executable, resources, and scheduler settings in `presub.sh`.

`presub.sh` is a template. It does not validate the inputs and is not permission
to submit a job. Inspect the directories and links first, then run through the
project's procedure. Afterward, check exit status, logs, temperature, energy,
pressure, cells, and trajectories, and identify failed or incomplete samples.

## 4. Post-process and continue

After DFT results have completed and passed review, use the converter for the DFT
code:

For VASP results:

~~~bash
gpumdkit.sh -out2xyz <scf_results_directory>
~~~

For CP2K logs and structure files, first enter the CP2K results root. The following command recursively scans `.log` files and their accompanying `.xyz`/`.inp` files under the current directory. Successful conversion writes `cp2k_exyz.xyz` there, with processing details in `Logfile.txt`.

~~~bash
gpumdkit.sh -cp2k2xyz
~~~

The VASP `-out2xyz`/menu 101 route creates `NEPdataset/train.xyz` in the current
directory; the Python `-out2exyz` route writes `train.xyz` there. Continue with
the file actually produced by the selected route and preserve the raw DFT outputs.

Use analysis and sampling tools to inspect the resulting extxyz and record excluded
directories.

## Active-Learning Style Workflow

For a NEP data iteration, proceed in this order: MD trajectory, geometry filtering,
NepTrain FPS, DFT, conversion, and dataset checks. See [Active Learning Workflow](active_learning_workflow.md)
for the manually reviewed protocol.

## Practical Notes

- Workflow scripts depend on the cluster, templates, and software versions; check a small number of structures first.
- Preserve templates and source structures, and inspect generated directories, links, and file names before continuing.

## AI Assistance and Authorization Boundaries

Workflow scripts depend on the cluster, templates, and software versions. An
agent must not infer the potential, species/type mapping, temperatures, time
steps, ensembles, pressures, run lengths, selection thresholds, resources, or
convergence criteria from examples or another material system. `INCAR`,
`POTCAR`, `KPOINTS`, `run.in`, CP2K templates, and `nep.txt` are not automatically
approved files. Preparing directories does not authorize DFT, GPUMD, LAMMPS, or
NEP execution; execution and scheduler submission require an explicit user request.

For direct source calls, first confirm the checkout root, interpreter, and
dependencies, then use the target script's own help. Ordinary users should use
the existing menu. Do not invent a CLI route for a menu-only function. Stop and
report parser errors, missing files, NaN/Inf values, failed outputs, and
unexplained warnings.

Relevant skill references:

- `skills/gpumdkit-skill/references/workflows.md`
- `skills/gpumdkit-skill/references/sampling.md`
- `skills/gpumdkit-skill/references/format-conversion.md`
- `skills/gpumdkit-skill/references/nep-data.md`
- `skills/gpumdkit-skill/references/nep-parameters.md`

For a reproducible iteration, record model and data revisions, approved decisions,
exact commands, working directories, versions, exit status, warnings, generated
files, exclusions, validation results, and known limitations.
