<div align="center">
  <h1>📋 Command Reference</h1>
  <p style="text-align: justify;">This page lists the stable GPUMDkit command-line shortcuts. Interactive mode is available for workflows with many choices.</p>
</div>

The source table is maintained in `docs/command_reference.tsv`.

If you are new to GPUMDkit, start with [Quick Start](quick_start.md). For
agent-assisted work, use the [GPUMDkit Agent Skill](simulation_and_postprocessing.md).
For remote browser access, see [Remote Web Console](remote_web_console.md).
The `-h` form is available only for commands that expose their own option help;
otherwise use `gpumdkit.sh -h`, the relevant module help, or interactive mode.

## `gpumdkit.sh -h` Output

```text
 +-------------------------------------------------------------------------------------------------------+
 |                          GPUMDkit 1.5.8 (2026-09-29)       Command Help                               |
 +-------------------------------------------------------------------------------------------------------+
 |                                          MAIN FUNCTIONS                                               |
 +-------------------------------------------------------------------------------------------------------+
 | -h            Show this help table            | -plt <type>        Plot and visualization tools       |
 | -calc <type>  Calculator tools                | -time <gpumd|wpe|nep|gnep>  Time-consuming analyzer   |
 | -update       Update GPUMDkit                 | -clean             Clean extra files in current dir   |
 | -skill        Show GPUMDkit agent skill info  | -doctor            Check Python environment           |
 | -prediction   Write NEP prediction .out files | -prediction_dpa    Write DPA prediction .out files    |
 | -server       Start the GPUMDkit web server   |                                                       |
 +-------------------------------------------------------------------------------------------------------+
 |                                         FORMAT CONVERSION                                             |
 +-------------------------------------------------------------------------------------------------------+
 | -out2xyz      OUTCAR -> extxyz (shell)        | -out2xyz_bec       OUTCAR -> extxyz with BEC          |
 | -out2exyz     OUTCAR -> extxyz (python)       | -xyz2dp            extxyz -> DeepMD npy               |
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
 | -shift_energy Interactive energy shift        |                                                       |
 +-------------------------------------------------------------------------------------------------------+
 | Python option help: gpumdkit.sh -<option> -h    Plot list: gpumdkit.sh -plt -h                        |
 +-------------------------------------------------------------------------------------------------------+
```

## Main

| Command | Syntax | Description |
|---|---|---|
| `-h` | `gpumdkit.sh -h` | Show general help |
| `-doctor` | `gpumdkit.sh -doctor` | Check Python and GPUMDkit package availability |
| `-update` | `gpumdkit.sh -update` | Update GPUMDkit |
| `-clean` | `gpumdkit.sh -clean` | Clean extra files in the current directory |
| `-server` | `gpumdkit.sh -server [port] [-b <address>] [-pw <password>]` | Start the Web Console for the current directory; use an SSH tunnel for remote access. |
| `-prediction` | `gpumdkit.sh -prediction <input.xyz> <nep.txt> [workers]` | Write NEP-compatible energy, force, stress, and virial prediction files |
| `-prediction_dpa` | `gpumdkit.sh -prediction_dpa <input.xyz> <dpa_model>` | Write DPA training-set prediction files in the current directory |

## Format Conversion

| Command | Syntax | Description |
|---|---|---|
| `-out2xyz` | `gpumdkit.sh -out2xyz <dir>` | OUTCAR to extxyz, shell version |
| `-out2xyz_bec` | `gpumdkit.sh -out2xyz_bec <dir>` | OUTCAR to extxyz with BEC labels |
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
| `-get_frame` | `gpumdkit.sh -get_frame <input.xyz> <frame_index>` | Extract one frame (1-based index) |
| `-clean_xyz` | `gpumdkit.sh -clean_xyz <input.xyz> <output.xyz>` | Remove extra extxyz properties |
| `-frame_range` | `gpumdkit.sh -frame_range <input.xyz> <start_frac> <end_frac>` | Extract frames by fractional range |
| `-dp2xyz` | `gpumdkit.sh -dp2xyz <input_dir/> [output.xyz]` | DeepMD npy datasets to extxyz |
| `-xyz2dp` | `gpumdkit.sh -xyz2dp` | Convert labeled extxyz to DeepMD npy data through the interactive prompts |

## Calculators

| Command | Syntax | Description |
|---|---|---|
| `-calc ionic-cond` | `gpumdkit.sh -calc ionic-cond <element> <charge>` | Ionic conductivity |
| `-calc nep` | `gpumdkit.sh -calc nep <input.xyz> <output.xyz> <nep.txt>` | NEP property prediction |
| `-prediction` | `gpumdkit.sh -prediction <input.xyz> <nep.txt> [workers]` | NEP-compatible prediction output files |
| `-prediction_dpa` | `gpumdkit.sh -prediction_dpa <input.xyz> <dpa_model>` | DPA training-set prediction output files |
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

`-calc nep` writes an output extxyz containing properties predicted by a NEP
model. `-prediction` writes the separate NEP prediction files used by training
and validation workflows.

## Analyzers

| Command | Syntax | Description |
|---|---|---|
| `-range` | `gpumdkit.sh -range <input.xyz> <energy\|force\|virial> [hist]` | Property range analysis |
| `-shift_energy` | `gpumdkit.sh -shift_energy` | Interactive energy reference shifting |
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

Common types include `train`, `prediction` (alias: `test`), `thermo`, `msd`, `sdc`, `rdf`, `emd`, `emd2`, `nemd`, `hnemd`, `pdos`, and `plane-grid`.

## Utilities

| Command | Syntax | Description |
|---|---|---|
| `-time` | `gpumdkit.sh -time <gpumd\|wpe\|nep\|gnep>` | Monitor GPUMD, WPE, NEP, or GNEP progress |
| `-nep_modifier` | `gpumdkit.sh -nep_modifier [nep.txt] [nep.restart\|-] [nep.in\|-]` | Safely modify and export a NEP4 model package |
| `-pynep` | `gpumdkit.sh -pynep` | Deprecated PyNEP FPS sampling |

### NEP model modifier

`-nep_modifier` adjusts NEP4 model capacity, adds or removes species, and exports
model packages. It requires `calorine >= 3.4`. Expansion, reduction, and adding
species also require a matching `nep.restart`; providing the source `nep.in`
is recommended.

```bash
gpumdkit.sh -nep_modifier
gpumdkit.sh -nep_modifier models/nep.txt models/nep.restart models/nep.in
gpumdkit.sh -nep_modifier -h
```

Use `-` in place of the restart path to omit it; species removal and retention
remain available. Continue training with the `.in`, `.txt`, and `.restart`
from the same exported package, then retrain and validate before production use.

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
