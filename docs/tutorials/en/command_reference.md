<div align="center">
  <h1>📋 Command Reference</h1>
  <p style="text-align: justify;">This page lists the stable GPUMDkit command-line shortcuts. Interactive mode is available for workflows with many choices.</p>
</div>

The source table is maintained in `docs/command_reference.tsv`.

## `gpumdkit.sh -h` Output

```text
+-------------------------------------------------------------------------------------------------------+
|                          GPUMDkit 1.5.6 (dev) (2026-07-10) Command Help                               |
+-------------------------------------------------------------------------------------------------------+
|                                          MAIN FUNCTIONS                                               |
+-------------------------------------------------------------------------------------------------------+
| -h            Show this help table            | -plt <type>        Plot and visualization tools       |
| -calc <type>  Calculator tools                | -time <gpumd|nep>  Time-consuming analyzer            |
| -update       Update GPUMDkit                 | -clean             Clean extra files in current dir   |
| -skill        Show GPUMDkit agent skill info  | -doctor           Check Python environment          |
+-------------------------------------------------------------------------------------------------------+
|                                         FORMAT CONVERSION                                             |
+-------------------------------------------------------------------------------------------------------+
| -out2xyz      OUTCAR -> extxyz (shell)        | -out2exyz          OUTCAR -> extxyz (python)          |
| -cp2k2xyz     CP2K log -> xyz                 | -xdat2exyz         XDATCAR -> extxyz                  |
| -cif2pos      cif -> POSCAR                   | -cif2exyz          cif -> extxyz                      |
| -pos2exyz     POSCAR -> extxyz                | -exyz2pos          extxyz -> POSCAR                   |
| -pos2lmp      POSCAR -> LAMMPS data           | -lmp2exyz          LAMMPS dump -> extxyz              |
| -traj2exyz    ASE traj -> extxyz              | -replicate         Replicate structure                |
| -addgroup     Add group labels                | -addweight         Add structure weight in extxyz     |
| -clean_xyz    Clean extra info in extxyz      | -get_frame         Extract specific frame             |
| -frame_range  Extract frames by range         | -dp2xyz            DeepMD npy -> extxyz               |
+-------------------------------------------------------------------------------------------------------+
|                                            ANALYSIS                                                   |
+-------------------------------------------------------------------------------------------------------+
| -range        Energy/force/virial statistics  | -analyze_comp      Analyze composition                |
| -chem_species Analyze chemical species        | -cbc               Charge balance check               |
| -min_dist     Min distance (no PBC)           | -min_dist_pbc      Min distance with PBC              |
| -filter_dist  Filter by min_dist (no PBC)     | -filter_dist_pbc   Filter by min_dist (PBC)           |
| -pda          Probability density analysis    | -filter_box        Filter by box-edge length          |
| -pynep        Deprecated PyNEP sampling       | -nep_modifier      Modify NEP model interactively     |
+-------------------------------------------------------------------------------------------------------+
| Python option help: gpumdkit.sh -<option> -h    Plot list: gpumdkit.sh -plt -h                     |
+-------------------------------------------------------------------------------------------------------+
```

## Main

| Command | Syntax | Description |
|---|---|---|
| `-h` | `gpumdkit.sh -h` | Show general help |
| `-doctor` | `gpumdkit.sh -doctor` | Check Python and GPUMDkit package availability |
| `-update` | `gpumdkit.sh -update` | Update GPUMDkit |
| `-clean` | `gpumdkit.sh -clean` | Clean extra files in the current directory |

## Format Conversion

| Command | Syntax | Description |
|---|---|---|
| `-out2xyz` | `gpumdkit.sh -out2xyz <dir>` | OUTCAR to extxyz, shell version |
| `-out2exyz` | `gpumdkit.sh -out2exyz <dir>` | OUTCAR to extxyz, Python version |
| `-cp2k2xyz` | `gpumdkit.sh -cp2k2xyz` | CP2K output to xyz/extxyz |
| `-xdat2exyz` | `gpumdkit.sh -xdat2exyz <XDATCAR> <output.xyz>` | XDATCAR to extxyz |
| `-cif2pos` | `gpumdkit.sh -cif2pos <input.cif> <output.vasp>` | CIF to POSCAR/VASP |
| `-cif2exyz` | `gpumdkit.sh -cif2exyz <input.cif> <output.xyz>` | CIF to extxyz |
| `-pos2exyz` | `gpumdkit.sh -pos2exyz <POSCAR> <output.xyz>` | POSCAR to extxyz |
| `-exyz2pos` | `gpumdkit.sh -exyz2pos <input.xyz>` | extxyz frames to POSCAR files |
| `-pos2lmp` | `gpumdkit.sh -pos2lmp <POSCAR> <output.data>` | POSCAR to LAMMPS data |
| `-lmp2exyz` | `gpumdkit.sh -lmp2exyz <dump> <element...>` | LAMMPS dump to extxyz |
| `-traj2exyz` | `gpumdkit.sh -traj2exyz <input.traj> <output.xyz>` | ASE trajectory to extxyz |
| `-replicate` | `gpumdkit.sh -replicate <input> <output> a b c` | Replicate by cell factors |
| `-replicate` | `gpumdkit.sh -replicate <input> <output> <target_num>` | Replicate toward a target atom count |
| `-addgroup` | `gpumdkit.sh -addgroup <POSCAR> <element...>` | Add GPUMD group labels |
| `-addweight` | `gpumdkit.sh -addweight <input.xyz> <output.xyz> <weight>` | Add structure weights |
| `-get_frame` | `gpumdkit.sh -get_frame <input.xyz> <frame_index>` | Extract one frame |
| `-clean_xyz` | `gpumdkit.sh -clean_xyz <input.xyz> <output.xyz>` | Remove extra extxyz properties |
| `-frame_range` | `gpumdkit.sh -frame_range <input.xyz> <start_frac> <end_frac>` | Extract frames by fractional range |
| `-dp2xyz` | `gpumdkit.sh -dp2xyz <input_dir/> [output.xyz]` | DeepMD npy datasets to extxyz |

## Calculators

| Command | Syntax | Description |
|---|---|---|
| `-calc ionic-cond` | `gpumdkit.sh -calc ionic-cond <element> <charge>` | Ionic conductivity |
| `-calc nep` | `gpumdkit.sh -calc nep <input.xyz> <output.xyz> <nep.txt>` | NEP property prediction |
| `-calc des` | `gpumdkit.sh -calc des <input.xyz> <output.npy> <nep.txt> <element>` | NEP descriptors |
| `-calc doas` | `gpumdkit.sh -calc doas <input.xyz> <nep.txt> <output.txt>` | Density of atomistic states |
| `-calc neb` | `gpumdkit.sh -calc neb <initial.xyz> <final.xyz> <n_images> <nep.txt>` | NEB with a NEP model |
| `-calc minimize` | `gpumdkit.sh -calc minimize <structure> <nep.txt> [fmax] [max_steps]` | Structure minimization |
| `-calc msd` | `gpumdkit.sh -calc msd <trajectory.xyz> <element> <dt_fs> [max_corr_steps]` | MSD from trajectory |
| `-calc nlist` | `gpumdkit.sh -calc nlist [args...]` | Neighbor lists |
| `-calc disp` | `gpumdkit.sh -calc disp [args...]` | Displacements |
| `-calc avg-struct` | `gpumdkit.sh -calc avg-struct [args...]` | Averaged structure |
| `-calc oct-tilt` | `gpumdkit.sh -calc oct-tilt [args...]` | Octahedral tilt |
| `-calc pol-abo3` | `gpumdkit.sh -calc pol-abo3 [args...]` | ABO3 local polarization |

## Analyzers

| Command | Syntax | Description |
|---|---|---|
| `-range` | `gpumdkit.sh -range <input.xyz> <energy\|force\|virial> [hist]` | Property range analysis |
| `-analyze_comp` | `gpumdkit.sh -analyze_comp <input.xyz>` | Composition analysis |
| `-chem_species` | `gpumdkit.sh -chem_species <input.xyz>` | Chemical species list |
| `-cbc` | `gpumdkit.sh -cbc <input.xyz>` | Charge-balance check |
| `-min_dist` | `gpumdkit.sh -min_dist <input.xyz>` | Minimum distance without PBC |
| `-min_dist_pbc` | `gpumdkit.sh -min_dist_pbc <input.xyz>` | Minimum distance with PBC |
| `-filter_dist` | `gpumdkit.sh -filter_dist <input.xyz> <min_dist>` | Distance filtering |
| `-filter_dist_pbc` | `gpumdkit.sh -filter_dist_pbc <input.xyz> <min_dist>` | PBC-aware distance filtering |
| `-filter_box` | `gpumdkit.sh -filter_box <input.xyz> <edge_limit>` | Box-edge filtering |
| `-filter_value` | `gpumdkit.sh -filter_value <input.xyz> <property> <threshold>` | Property threshold filtering |
| `-filter_range` | `gpumdkit.sh -filter_range <input.xyz> <element1> <element2> <min_dist> <max_dist>` | Element-pair distance range filtering |
| `-pda` | `gpumdkit.sh -pda <ref_struct> <trajectory.xyz> <species> <interval>` | Probability density analysis |

## Visualization

Use:

```bash
gpumdkit.sh -plt <type> [options]
gpumdkit.sh -plt -h
```

Common types include `train`, `prediction` (alias: `test`), `thermo`, `msd`, `sdc`, `rdf`, `emd`, `nemd`, `hnemd`, `pdos`, and `plane-grid`.

## Utilities

| Command | Syntax | Description |
|---|---|---|
| `-time` | `gpumdkit.sh -time <gpumd\|nep>` | Monitor GPUMD or NEP progress |
| `-nep_modifier` | `gpumdkit.sh -nep_modifier [nep.txt] [nep.restart\|-] [nep.in\|-]` | Safely modify and export a NEP4 model package |
| `-pynep` | `gpumdkit.sh -pynep` | Deprecated PyNEP FPS sampling |

### NEP model modifier

`-nep_modifier` starts a guided two-column editor built on calorine's NEP model
modification API. It is intended for developing an existing NEP4 model further:
for example, increasing capacity, extracting a chemical subset from a foundation
model, or adding a new species without discarding the learned parameters for the
original species.

#### Requirements and startup

The command requires `calorine >= 3.4`. Expansion, reduction, and adding species
also require the `nep.restart` that matches the model because it contains the
SNES parameter means and exploration widths. A source `nep.in` is recommended so
that non-architecture training settings can be retained.

```bash
# Prompt for files; defaults are resolved beside the selected nep.txt
gpumdkit.sh -nep_modifier

# Load a complete model package directly
gpumdkit.sh -nep_modifier models/nep.txt models/nep.restart models/nep.in

# Load without restart; species removal/retention remains available
gpumdkit.sh -nep_modifier models/nep.txt - models/nep.in

# Display command help without importing calorine
gpumdkit.sh -nep_modifier -h
```

#### What each operation does

| Menu operation | Purpose and model effect |
|---|---|
| **Expand model capacity** | Increases neurons, enables 4-/5-body or `q_*` descriptor terms, or adds a charge head. Existing trained parameter means are retained and new parameters are initialized for continued optimization. |
| **Reduce model capacity** | Keeps the highest-ranked neurons while discarding lower-ranked ones, disables descriptor terms, or removes the charge head. This can provide a smaller starting model, but accuracy and speed must be measured after retraining. |
| **Add chemical species** | Adds a species-specific ANN subnetwork and every descriptor-weight pair involving the new species. The new parameters are untrained; an explicit seed makes their initialization reproducible. |
| **Remove chemical species** | Removes selected species together with their ANN subnetworks and descriptor-weight pairs. This is convenient when only a few species should be discarded. |
| **Keep selected species** | Retains the listed species and removes all others. This is the more convenient inverse operation when extracting a small subset from a large model. |
| **Inspect current model** | Shows species order, cutoffs, descriptor switches and dimensions, neuron and parameter counts, restart state, ZBL, and charge mode. |
| **Review pending changes** | Lists accepted operations, their arguments, changed architecture fields, and export state before files are written. |
| **Export model package** | Writes a common, collision-free package suffix for the model, optional restart, updated input, and provenance summary. |

#### Example: expand one model in a reproducible workflow

After loading the package, enter `1`, select the desired expansion fields, and
enter their target values. Multiple fields may be selected together; calorine
applies them in one `augment()` call. For example:

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

This example changes the neuron count from 50 to 60 and enables `q_112`. Choose
`7` to review the recorded arguments and architecture changes, then choose `8`
to export. Continued training must use the generated `.in`, `.txt`, and
`.restart` from the same package; using a stale `nep.in` can make the restart
layout inconsistent with the modified parameter count.

#### Example: extract or extend a chemical model

To extract a Li-O submodel, choose `5`, enter `Li O`, review the reported species
order, and export. To add carbon, choose `3`, enter `C`, then provide a deliberate
random seed. A typewise-cutoff model additionally asks for carbon's radial and
angular cutoffs. These cutoffs are scientific choices that must be taken from
the intended training design; the tool does not choose them.

The export contains `*_modified.txt`, `*_modified.restart` when restart data are
loaded, `*_modified.in`, and `*_modified.changes.txt`. Check the generated input
with the exact NEP executable version intended for continued training, then
retrain and validate the modified model on representative reference data before
using it in production simulations.

#### Required calorine citation

This GPUMDkit feature directly uses calorine's model modification
implementation. Research that uses this feature should cite calorine as
requested by its developers:

E. Lindgren, J. M. Rahm, E. Fransson, F. Eriksson, N. Österbacka, Z. Fan, and
P. Erhart, “calorine: A Python package for constructing and sampling
neuroevolution potential models,” *Journal of Open Source Software* **9**(95),
6264 (2024), <https://doi.org/10.21105/joss.06264>.

The detailed operation guide is available in the
[NEP modifier README](https://github.com/zhyan0603/GPUMDkit/tree/main/Scripts/utils/nep_modifier)
and the [official calorine model-modification tutorial](https://calorine.materialsmodeling.org/dev/get_started/modifying_models.html).
