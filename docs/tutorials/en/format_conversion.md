<div align="center">
  <h1>🔄 Format Conversion</h1>
  <p style="text-align: justify;">GPUMDkit contains utilities for converting common atomistic simulation formats, with extra handling for metadata such as group labels, weights, and frame extraction.</p>
</div>

## What it does

This module converts structure and trajectory files between common formats used in computational materials science. It supports VASP, LAMMPS, CP2K, ABACUS, CIF, ASE trajectories, and extxyz. It also provides tools for adding group labels, weights, extracting frames, and replicating structures.

## Before you start

**Script location:** `Scripts/format_conversion/`

Make sure GPUMDkit is installed. See [Quick Start](quick_start.md) for installation instructions.

## Supported formats

The format conversion module covers:

- **VASP**: `POSCAR`, `CONTCAR`, `OUTCAR`, `XDATCAR`
- **LAMMPS**: data files and dump trajectories
- **CP2K**: output logs, position, force, and cell files
- **ABACUS**: SCF/MD logs
- **CIF**: crystallographic structure files
- **ASE trajectory**: `.traj`
- **DeepMD**: `npy` datasets
- **extxyz**: a common structure format for GPUMD and NEP

## Interactive Mode

Open GPUMDkit:

```bash
gpumdkit.sh
```

Choose:

```text
1) Format Conversion
```

The format conversion menu is:

```text
+-------------------------------------------------------------+
|                   FORMAT CONVERSION TOOLS                   |
+-------------------------------------------------------------+
| 101) VASP to extxyz            106) Add group labels        |
| 102) MTP to extxyz             107) Add weight to extxyz    |
| 103) CP2K to extxyz            108) Extract frame extxyz    |
| 104) ABACUS to extxyz          109) Clean XYZ info          |
| 105) extxyz to POSCAR          110) Replicate structure     |
+-------------------------------------------------------------+
| out2exyz) OUTCAR to extxyz     xdat2exyz) XDATCAR to extxyz |
| pos2exyz) POSCAR to extxyz     pos2lmp)   POSCAR to LAMMPS  |
| cif2pos)  CIF to POSCAR        lmp2exyz)  LAMMPS to extxyz  |
| cif2exyz) CIF to extxyz        traj2exyz) ASE traj to extxyz|
| dp2xyz)   DeepMD to extxyz                                  |
+-------------------------------------------------------------+
| 000) Return to main menu                                    |
+-------------------------------------------------------------+
Input the function number or converter keyword:
```

### What Each Entry Does

| Entry | Called script | Function | When to Use |
|-------|---------------|----------|-------------|
| `101` | `out2xyz.sh` | VASP OUTCAR to extxyz, shell version | convert VASP calculation directories |
| `102` | `mtp2xyz.py` | MTP cfg to extxyz | convert MTP training data |
| `103` | `cp2k_log2xyz.py` / `cp2k2xyz.py` | CP2K to extxyz | choose CP2K log/inp route or pos/frc/cell route |
| `104` | `abacus2xyz_scf.sh` / `abacus2xyz_md.sh` | ABACUS to extxyz | convert ABACUS SCF or MD output |
| `105` | `exyz2pos.py` | extxyz to POSCAR | write each extxyz frame as a POSCAR-style file |
| `106` | `add_groups.py` | add group labels | add atom group labels for GPUMD-related workflows |
| `107` | `add_weight.py` | add weights | assign training weights in extxyz |
| `108` | `get_frame.py` | extract frame | export one frame from an extxyz trajectory |
| `109` | `clean_xyz.py` | clean XYZ info | remove extra extxyz properties |
| `110` | `replicate.py` | replicate structure | build supercells by factors or target atom count |
| `out2exyz` | `out2exyz.py` | OUTCAR to extxyz, Python version | alternative VASP OUTCAR converter |
| `pos2exyz` | `pos2exyz.py` | POSCAR to extxyz | convert a single structure |
| `cif2pos` | `cif2pos.py` | CIF to POSCAR | prepare VASP input from CIF |
| `cif2exyz` | `cif2exyz.py` | CIF to extxyz | prepare GPUMDkit input from CIF |
| `xdat2exyz` | `xdatcar2exyz.py` | XDATCAR to extxyz | convert VASP MD trajectory |
| `pos2lmp` | `pos2lmp.py` | POSCAR to LAMMPS data | prepare LAMMPS input |
| `lmp2exyz` | `lmp2exyz.py` | LAMMPS dump to extxyz | convert LAMMPS trajectory |
| `traj2exyz` | `traj2exyz.py` | ASE traj to extxyz | convert ASE trajectory |
| `dp2xyz` | `dp2xyz.py` | DeepMD npy to extxyz | convert DeepMD datasets for GPUMD/NEP workflows |

## Quick Command Reference

| Source | Target | Command |
|--------|--------|---------|
| OUTCAR directory | extxyz | `gpumdkit.sh -out2xyz <dir>` |
| OUTCAR directory | extxyz | `gpumdkit.sh -out2exyz <dir>` |
| POSCAR | extxyz | `gpumdkit.sh -pos2exyz <POSCAR> <output.xyz>` |
| extxyz | POSCAR files | `gpumdkit.sh -exyz2pos <input.xyz>` |
| XDATCAR | extxyz | `gpumdkit.sh -xdat2exyz XDATCAR dump.xyz` |
| POSCAR | LAMMPS data | `gpumdkit.sh -pos2lmp POSCAR lammps.data` |
| LAMMPS dump | extxyz | `gpumdkit.sh -lmp2exyz dump.lammpstrj Li Y Cl` |
| CIF | POSCAR | `gpumdkit.sh -cif2pos input.cif POSCAR.vasp` |
| CIF | extxyz | `gpumdkit.sh -cif2exyz input.cif model.xyz` |
| ASE traj | extxyz | `gpumdkit.sh -traj2exyz input.traj output.xyz` |
| DeepMD npy directory | extxyz | `gpumdkit.sh -dp2xyz database train.xyz` |
| extxyz | clean extxyz | `gpumdkit.sh -clean_xyz input.xyz clean.xyz` |

## Understand input and output names

For a first conversion, use the POSCAR route because its two arguments make the
data flow explicit:

```text
$ gpumdkit.sh -pos2exyz -h
 Usage: gpumdkit.sh -pos2exyz <POSCAR> <output.xyz>

 Arguments:
   POSCAR       One or more VASP POSCAR files (supports wildcards)
   output.xyz   Output extxyz file
```

`POSCAR` is the input structure name; `output.xyz` is the extxyz file created by
the command. Shell wildcards are expanded by your shell, so check the matched
filenames before converting a group of files. After conversion, a lightweight
format check is:

```bash
head -n 2 model.xyz
```

The first line is the atom count. The second line contains extxyz metadata such
as the cell and available properties. It does not by itself validate that the
underlying calculation is physically appropriate; use the analyzer tutorials to
check a dataset before training or simulation.

## Common Examples

### Convert a DeepMD dataset to extxyz

Use `-dp2xyz` to recursively scan a DeepMD dataset directory and write an extxyz file:

```bash
gpumdkit.sh -dp2xyz database train.xyz
```

The input directory must contain DeepMD dataset files such as `type.raw`, `type_map.raw`, and `set.000/`. This command requires `dpdata` and `ase`.

### Convert VASP calculations to extxyz

**What it does:** Searches a directory for VASP OUTCAR files and converts them into a single extxyz file for NEP training or analysis.

**CLI mode:**

```bash
gpumdkit.sh -out2xyz ./vasp_results/
```

The shell version searches the target directory and converts VASP results into an extxyz file. If you prefer the Python implementation:

```bash
gpumdkit.sh -out2exyz ./vasp_results/
```

**Interactive mode:** Choose `101` from the format conversion menu. You will see:

```text
>-------------------------------------------------<
| Calling the script in Scripts/format_conversion |
| Script: out2xyz.sh                              |
| Developer: Yanzhou WANG (yanzhowang@gmail.com)  |
>-------------------------------------------------<
Input the directory containing OUTCARs
Example: ./
------------>>
```

**Output:** An extxyz file containing the converted structures, suitable for NEP training or further analysis.

### Add group labels

**What it does:** Adds atom group labels to a structure file. Group labels are required by some GPUMD-related workflows, such as species-specific MSD or diffusion calculations.

**CLI mode:**

```bash
gpumdkit.sh -addgroup POSCAR Li Y Cl
```

This command reads the input structure and writes an extxyz file with group information.

**Interactive mode:** Choose `106` from the format conversion menu. You will see:

```text
>-------------------------------------------------<
| Calling the script in Scripts/format_conversion |
| Script: add_groups.py                           |
| Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
>-------------------------------------------------<
Input <POSCAR> <element1> <element2> ...
Example: POSCAR Li Y Cl
------------>>
```

**Output:** An extxyz file with group labels added.

## Script Details

### POSCAR to extxyz

```bash
gpumdkit.sh -pos2exyz POSCAR model.xyz
```

Use this when you have a single VASP structure and want an extxyz output.

Interactive keyword: `pos2exyz`

```text
Input <POSCAR> <output.xyz>
Example: POSCAR model.xyz
------------>>
```

### extxyz to POSCAR

```bash
gpumdkit.sh -exyz2pos structures.xyz
```

This converts all frames in an extxyz file into `POSCAR_*.vasp` files. Frame indices are 0-based in most GPUMDkit scripts, but output filenames are meant for direct inspection and batch calculations.

Interactive entry: `105`

```text
Input the name of extxyz
Example: ./train.xyz
------------>>
```

### LAMMPS dump to extxyz

```bash
gpumdkit.sh -lmp2exyz dump.lammpstrj Li Y Cl
```

The element symbols must match the LAMMPS atom type IDs. For example, if type `1` is Li, type `2` is Y, and type `3` is Cl, the order should be `Li Y Cl`.

Interactive keyword: `lmp2exyz`

```text
Input <dump_file> <element1> <element2> ...
Example: dump.lammpstrj Li O
------------>>
```

### CIF conversion

```bash
gpumdkit.sh -cif2pos input.cif POSCAR.vasp
gpumdkit.sh -cif2exyz input.cif model.xyz
```

Use `-cif2pos` if you want a VASP-style structure, and `-cif2exyz` if the next step is GPUMDkit analysis.

### Add weights

```bash
gpumdkit.sh -addweight train.xyz train_weighted.xyz 5
```

This is useful when you want some structures to have a different training weight in a NEP dataset.

### Replicate structures

```bash
gpumdkit.sh -replicate POSCAR supercell.vasp 2 2 2
gpumdkit.sh -replicate POSCAR supercell.vasp 256
```

The first form uses explicit replication factors. The second form tries to build a supercell close to a target atom count.

### Extract one frame

```bash
gpumdkit.sh -get_frame dump.xyz 1000
```

This extracts frame index `1000` from an extxyz trajectory.

### Split multi-frame extxyz

`split_single_xyz.py` splits an extxyz file into individual frames, each written to a separate file.

```bash
python Scripts/format_conversion/split_single_xyz.py dump.xyz
```

This creates `model_0.xyz`, `model_1.xyz`, ... for each frame in the trajectory.

### MTP conversion

Convert MTP `.cfg` format to extxyz:

```bash
gpumdkit.sh    # Select: 1) Format Conversion → 102
```

Interactive prompt:

```text
Input <filename.cfg> <Symbol1 Symbol2 Symbol3 ...>
Example: train.cfg Pd Ag
------------>>
```

### ABACUS conversion

Convert ABACUS output to extxyz:

```bash
gpumdkit.sh    # Select: 1) Format Conversion → 104
```

The menu offers two options:

1. SCF output (`running_scf.log`)
2. MD output (`running_md.log`)

## Common Mistakes

| Problem | What to Check |
|---------|---------------|
| LAMMPS elements are wrong | Check the element order passed after the dump file |
| A trajectory has strange metadata | Try `-clean_xyz input.xyz clean.xyz` |
| A converted structure looks shifted | Inspect PBC/cell information in the source file |
| Frame extraction gives the wrong structure | Remember that frame indices are 0-based |

## Notes

If a Python package required by a specific converter is missing, Python will report it when that converter is used.
