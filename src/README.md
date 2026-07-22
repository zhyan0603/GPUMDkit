# Functions for GPUMDkit Interactive Mode

This directory contains shell script functions that implement the interactive menu system for `gpumdkit.sh`.

## Files

### f1_format_conversions.sh
Contains functions for format conversion operations (Interactive mode `option 1`).

**Functions:**

- `f101_out2xyz` - Convert VASP OUTCAR to extxyz
- `f102_mtp2xyz` - Convert MTP .cfg format to extxyz
- `cp2k2xyz_chenhua` - Convert CP2K to extxyz (method by Chen HUA)
- `cp2k2xyz_kexu` - Convert CP2K to extxyz (method by Ke XU)
- `f103_cp2k2xyz` - Convert CP2K to extxyz (main function)
- `f104_abacus2xyz` - Convert ABACUS output to extxyz
- `f105_extxyz2poscar` - Convert extxyz to POSCAR format
- `f106_add_group_labels` - Add group labels to structures
- `f107_add_weight` - Add weight labels to structures
- `f108_get_frame` - Get a single frame from extxyz
- `f109_clean_xyz` - Clean extra info in XYZ file
- `f110_replicate_structure` - Replicate unit cell to supercell
- `out2exyz` - Convert OUTCAR to extxyz (Python version)
- `pos2extxyz` - Convert POSCAR to extxyz
- `cif2poscar` - Convert CIF to POSCAR
- `cif2extxyz` - Convert CIF to extxyz
- `xdatcar2extxyz` - Convert XDATCAR to extxyz
- `poscar2lammps` - Convert POSCAR to LAMMPS data
- `lammps2extxyz` - Convert LAMMPS dump to extxyz
- `traj2extxyz` - Convert ASE trajectory to extxyz

### f2_sample_structures.sh
Contains functions for structure sampling (Interactive mode `option 2`).

**Functions:**
- `f201_sample_structures` - Sample structures using uniform or random methods
- `parallel_pynep_sample_structures` - Compatibility PyNEP sampling path used by `gpumdkit.sh -pynep`
- `f203_neptrain_sample_structures` - Sample structures using farthest point sampling (powered by NepTrain)
- `f204_perturb_structure` - Generate perturbed structures
- `f205_select_max_force_deviation_structs` - Select structures with high force deviations

### f3_workflows.sh
Contains workflow automation functions (Interactive mode `option 3`).

**Functions:**

- `cp2k_batch_pretreatment` - CP2K SCF batch pretreatment
- `f301_scf_batch_pretreatment` - SCF batch pretreatment (main function, supports VASP and CP2K)
- `f302_md_sample_batch_pretreatment_gpumd` - Batch setup for GPUMD MD simulations (source `Scripts/workflow/md_sample_batch_pretreatment_gpumd.sh`)
- `f303_md_sample_batch_pretreatment_lmp` - Batch setup for LAMMPS MD simulations (source `Scripts/workflow/md_sample_batch_pretreatment_lmp.sh`)

### f4_calculators.sh
Contains calculator functions (Interactive mode `option 4`).

**Functions:**
- `f401_calc_ionic_conductivity` - Calculate ionic conductivity from MSD data
- `f402_calc_properties_with_nep` - Calculate properties using NEP model
- `f403_calc_descriptors` - Calculate descriptors
- `f404_calc_doas` - Calculate density of atomistic states proposed by [Wang et al](https://doi.org/10.1002/anie.202215544)
- `f405_calc_neb` - Calculate nudged elastic band (NEB) by NEP
- `f406_calc_neighbor_list` - Build neighbor list for displacement analysis
- `f407_calc_displacement` - Calculate displacement from trajectory and neighbor list
- `f408_calc_averaged_structure` - Calculate averaged structure from trajectory
- `f409_calc_oct_tilt` - Calculate octahedral tilt from trajectory and B-O neighbor list
- `f410_calc_polarization_abo3` - Calculate local polarization for ABO3 from neighbor lists
- `f411_minimize_structure_by_nep` - Minimize structure using NEP
- `f412_calc_msd_from_trajectory` - Calculate MSD from trajectory

### f5_analyzers.sh
Contains analysis functions (Interactive mode `option 5`).

**Functions:**

- `f501_analyze_composition` - Analyze composition of extxyz file
- `f502_find_outliers` - Find outlier structures based on RMSE thresholds
- `f503_analyze_chem_species` - Analyze chemical species in extxyz file
- `f504_charge_balance_check` - Check charge balance of structures
- `f505_energy_force_virial_analyzer` - Analyze energy/force/virial range
- `f506_filter_structures_by_distance` - Filter structures by minimum distance (without PBC)
- `f507_get_min_dist` - Get minimum interatomic distance (with optional PBC)
- `f508_probability_density_analysis` - Probability density analysis for diffusion channels

### f6_plots.sh
Contains plotting functions (Interactive mode `option 6`).

**Functions:**

- `f6_plots_one_column` - Display plot menu in single-column layout
- `f6_plots_two_column` - Display plot menu in two-column layout

### f7_utilities.sh
Contains utility functions (Interactive mode `option 7`).

**Functions:**

- `f701_time_consuming_analyzer` - Real-time monitoring of GPUMD/NEP/gnep progress

## Usage

These functions are sourced by `gpumdkit.sh` during interactive mode operation. They are not intended to be called directly.

**Interactive mode:**

```
gpumdkit.sh
```

then you will see:

```bash
         ____ ____  _   _ __  __ ____  _    _ _
        / ___|  _ \| | | |  \/  |  _ \| | _(_) |_
       | |  _| |_) | | | | |\/| | | | | |/ / | __|
       | |_| |  __/| |_| | |  | | |_| |   <| | |_
        \____|_|    \___/|_|  |_|____/|_|\_\_|\__|

        GPUMDkit Version 1.5.6 (dev) (2026-07-10)
  Core Developer: Zihan YAN (yanzihan@westlake.edu.cn)
 Main Contributors: Denan LI, Xin WU, Zhoulin LIU & Chen HUA

 ----------------------- GPUMD -----------------------
  1) Format Conversion          2) Sample Structures
  3) Workflow                   4) Calculators
  5) Analyzer                   6) Visualization
  7) Utilities                  8) Help                
  0) Exit
 ------------>>
 Input the function number:
```

Each function:
1. Displays information about the script being called
2. Prompts user for required inputs
3. Calls the corresponding Python/Bash script from `Scripts/` directory
4. Shows the path to the called script

## Development

When adding new interactive functions:

1. Add the function to the appropriate `f*.sh` file
2. Update the menu in `gpumdkit.sh`
3. Ensure prompts clearly explain expected inputs
4. Display the script path after execution
5. Follow existing function naming pattern: `f<category><number>_<description>`

## Related

- **Scripts/**: Contains the actual implementation scripts
- **gpumdkit.sh**: Main entry point that sources these files
- **Command-line mode**: Direct script invocation bypasses these functions

---

For detailed usage of each script, see the README files in the corresponding `Scripts/` subdirectories.
