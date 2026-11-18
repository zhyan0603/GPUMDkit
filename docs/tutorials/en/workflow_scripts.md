<div align="center">
  <h1>🔧 Workflow Scripts</h1>
  <p style="text-align: justify;">Workflow scripts prepare batch directories, input files, and links for DFT, MD, and NEP data work.</p>
</div>

The workflow scripts prepare directories and files. They do not submit calculations.

**Script Location:** `Scripts/workflow/`

## Interactive Entry

Run `gpumdkit.sh` and choose `3) Workflow`. The entries are described in the table below.

<details markdown="1">
<summary>Show the full menu</summary>

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

</details>

| Menu | Workflow | Use it for |
|---|---|---|
| 301 | SCF batch pretreatment | Prepare DFT single-point directories for structures |
| 302 | MD sample batch pretreatment (gpumd) | Prepare multiple GPUMD sampling directories |
| 303 | MD sample batch pretreatment (lmp) | Prepare multiple LAMMPS sampling directories |

## 1. Prepare DFT single-point calculations {#scf-batch-pretreatment}

Work in a dedicated directory and keep a separate copy of the source structures.
Function 301 creates `struct_fp/`, `fp/`, the calculation directories, and
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

With one element set, the output stays flat:

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

For extxyz input, POSCAR species follow the trajectory-wide first-seen element
order. If all structures contain the same element set, 301 keeps the flat
layout and uses the shared `fp/POTCAR`. Different atom ratios with the same
element set do not create separate groups. If multiple element sets are
present, the POSCARs are stored under `struct_fp/<elements>/`, and each
calculation directory links to the corresponding `fp/POTCAR_<elements>` file:

~~~text
struct_fp/
├── Al_O/POSCAR_1.vasp
└── Al_O_Te_Bi/POSCAR_2.vasp
fp/
├── INCAR
├── POTCAR_Al_O
└── POTCAR_Al_O_Te_Bi
<prefix>_1/POTCAR -> ../fp/POTCAR_Al_O
<prefix>_1/KPOINTS -> ../fp/KPOINTS
<prefix>_2/POTCAR -> ../fp/POTCAR_Al_O_Te_Bi
<prefix>_2/KPOINTS -> ../fp/KPOINTS
~~~

Prepare each listed POTCAR in `fp/`; its element order must match the POSCAR
header. All groups share `fp/INCAR`, and every calculation directory links to
`../fp/KPOINTS`. Provide `fp/KPOINTS` for explicit k-point input. If the target
file is absent, the link stays dangling and VASP uses `KSPACING` in `INCAR`.
The `prefix` becomes the calculation-directory prefix. When `out2xyz` is used
later, the directory name is written as `config_type` in the extxyz output.

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

With a source installation you can also call the script directly from your data
directory:

~~~bash
python "${GPUMDkit_path}/Scripts/workflow/scf_batch_pretreatment_cp2k.py" train.xyz template.inp calc
~~~

### ABACUS SCF preparation (301 -> 3)

#### Before running 301

1. **Structure file**. Put the structures into one extxyz, for example `frames.xyz`. Every frame needs a valid three-dimensional cell. If you start from a POSCAR, convert it with `gpumdkit.sh -pos2exyz POSCAR frame.xyz`.
2. **Basis directory**. Put exactly one `.upf` pseudopotential and one compatible `.orb` numerical orbital for every element in the batch directly in one directory. Choose the pairs from a compatible library. Filenames must start with the case-sensitive element symbol. A missing file or multiple matches for an element stops preparation before any calculation directories are created.

For a Li–La–Zr–O batch, the files can be arranged as follows:

~~~text
work/
├── frames.xyz
└── basis/
    ├── Li.upf
    ├── Li_gga_9au_100Ry_4s1p.orb
    ├── La_ONCV_PBE-1.2.upf
    ├── La_gga_9au_100Ry_4s2p2d1f.orb
    ├── Zr_ONCV_PBE-1.0.upf
    ├── Zr_gga_10au_100Ry_4s2p2d1f.orb
    ├── O.upf
    └── O_gga_6au_100Ry_2s2p1d.orb
~~~

Start `gpumdkit.sh` and choose **3 → 301 → 3**. Enter the structure file, the output directory prefix and the basis directory.

~~~text
frames.xyz LLZO ./basis
~~~

The script creates:

~~~text
current directory/
├── fp/                 # put the shared INPUT (and KPT) here
├── LLZO_1/
│   ├── STRU            # this frame; species, masses, and file names filled in
│   ├── INPUT -> ../fp/INPUT
│   └── KPT -> ../fp/KPT
├── LLZO_2/
└── presub.sh
~~~

Each `STRU` lists only the elements present in that frame, so structures with
different compositions need no grouping. Lattice vectors and Cartesian
positions are written in Angstrom with `LATTICE_CONSTANT 1.889726125457828`.

## 2. Prepare MD sampling directories {#md-sampling-pretreatment}

### GPUMD (302)

In a fresh working directory, provide one `.xyz` file or one or more `.vasp`
files. The script rejects mixed formats and multiple XYZ files before moving
any inputs.

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

After preparation, put `nep.txt` and one `run_*.in` per sample into `md/`. Each
sample directory links to them: `sample_1/run.in` points to `md/run_1.in`,
`sample_2/run.in` to `md/run_2.in`, and so on.

### LAMMPS (303)

Use the same input rules as for GPUMD (302).

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

After preparation, put `lmprun.in` and `nep.txt` into `md/`; every sample directory links to them.

## A Worked Example: 302 on Two Structures

Starting from a folder that only contains `a.vasp` (Li) and `b.vasp` (Si), run `gpumdkit.sh`, choose `3` and then `302`. The script renames the structures, converts each to extxyz and creates one sample folder per structure:

```text
 Progress: [##################################################] 100% (2/2)
        2 model.xyz files were generated.
 |    ATTENTION: Place run_*.in and nep.txt in 'md' Dir.   |
```

```text
.
├── md/                      # empty: put nep.txt and run_1.in, run_2.in here
├── presub.sh                # template loop over sample_*; adapt it to your scheduler
├── struct_md/
│   ├── POSCAR_1.vasp        # the input structures, renamed (moved out of the top folder)
│   ├── POSCAR_2.vasp
│   ├── model_1.xyz
│   └── model_2.xyz
├── sample_1/
│   ├── model.xyz -> ../struct_md/model_1.xyz
│   ├── nep.txt   -> ../md/nep.txt
│   └── run.in    -> ../md/run_1.in
└── sample_2/
    └── ...                  # the same links, pointing to *_2 files
```

Only links are created, so `sample_1/nep.txt` and `sample_1/run.in` stay dangling until you place `nep.txt` and the `run_*.in` files in `md/`. That is the moment to check that each `run_N.in` really belongs to structure `N`. The input `.vasp` files are **moved** into `struct_md/`, so keep your own copy of anything you cannot regenerate.

## 3. Check before running

Before submitting, check:

- the number of structures, atom and species order, cells, and number of generated directories;
- `INCAR`, `POTCAR`, and `KPOINTS` in `fp/`, or `nep.txt`, `run_*.in`, and `lmprun.in` in `md/`;
- the links and file names in each directory;
- the executable, resources, and scheduler settings in `presub.sh`, which is only a template.

Workflow scripts depend on your cluster, templates, and software versions, so try a few structures first. After the runs, look for failed or incomplete samples before continuing.

## 4. Convert the results

After the DFT runs finish, convert them with `gpumdkit.sh -out2xyz <dir>` (VASP, writes `NEPdataset/train.xyz`), `-out2xyz_bec <dir>` (VASP with Born effective charges), `-cp2k2xyz` (CP2K, run from the results root, writes `cp2k_exyz.xyz`), or menu `104` -> SCF (ABACUS, reads each `OUT.*/running_scf.log`, writes `NEPdataset/train.xyz`). See [Format Conversion](format_conversion.md) for details, and [Active Learning Workflow](active_learning_workflow.md) for how these steps fit into a full NEP data iteration.
