# Source Functions for GPUMDkit Interactive Mode

This directory contains shell script functions that implement the interactive menu system for `gpumdkit.sh`.

## Files

### f1_format_conversions.sh
Contains functions for format conversion operations (Interactive mode option 1).

**Functions:**
- `f101_out2xyz()` - Convert VASP OUTCAR or vasprun.xml to extxyz
- `f102_mtp2xyz()` - Convert MTP .cfg format to extxyz
- `f103_cp2k2xyz()` - Convert CP2K output to extxyz
- `f104_abacus2xyz()` - Convert ABACUS output to extxyz
- `f105_extxyz2poscar()` - Convert extxyz to POSCAR format

### f2_sample_structures.sh
Contains functions for structure sampling (Interactive mode option 2).

**Functions:**
- `f201_sample_structures()` - Sample structures using uniform or random methods
- `f202_pynep_sample_structures()` - Sample using PyNEP farthest point sampling
- `f203_neptrain_sample_structures()` - Sample using NEPtrain selection
- `f204_perturb_structure()` - Generate perturbed structures for training
- `f205_select_max_force_deviation_structs()` - Select structures with high force deviations

### f3_workflows.sh
Contains workflow automation functions (Interactive mode option 3).

**Functions:**
- `f301_scf_batch_pretreatment()` - Batch preprocessing for DFT SCF calculations
- `f302_md_sample_batch_pretreatment_gpumd()` - Batch setup for GPUMD MD simulations
- `f303_md_sample_batch_pretreatment_lmp()` - Batch setup for LAMMPS MD simulations

### f4_calculators.sh
Contains calculator functions (Interactive mode option 4).

**Functions:**
- `f401_calc_ionic_conductivity()` - Calculate ionic conductivity from MSD data
- `f402_calc_properties_with_nep()` - Calculate properties using NEP model
- `f403_calc_descriptors()` - Calculate NEP descriptors
- `f404_calc_doas()` - Calculate density of atomistic states
- `f405_calc_neb()` - Nudged elastic band calculations

### f5_analyzer.sh
Contains analysis functions (Interactive mode option 5).

**Functions:**
- `f501_analyze_composition()` - Analyze composition of structure files

## Usage

These functions are sourced by `gpumdkit.sh` during interactive mode operation. They are not intended to be called directly.

**Interactive mode workflow:**
```bash
gpumdkit.sh
# Menu appears
# Select category (1-6)
# Select specific function
# Follow prompts
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
