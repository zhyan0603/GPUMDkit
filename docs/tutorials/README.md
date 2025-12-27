# GPUMDkit Tutorials

Welcome to **GPUMDkit** tutorials! This guide helps you learn how to use GPUMDkit for GPUMD and NEP workflows.

## What is GPUMDkit?

GPUMDkit is a command-line toolkit that simplifies working with GPUMD (Graphics Processing Units Molecular Dynamics) and NEP (neuroevolution potential). It provides:

- **Easy-to-use interface**: Both interactive menu and command-line modes
- **Format conversion tools**: Convert between VASP, LAMMPS, CP2K, and other formats
- **Analysis tools**: Analyze structures, forces, energies, and more
- **Visualization scripts**: Plot NEP training, MD results, and material properties
- **Workflow automation**: Batch processing for DFT and MD calculations

## Getting Started

### Installation

```bash
# Clone the repository
git clone --depth 1 https://github.com/zhyan0603/GPUMDkit.git

# Add to your ~/.bashrc
export GPUMDkit_path=/path/to/GPUMDkit
export PATH=${GPUMDkit_path}:${PATH}
source ${GPUMDkit_path}/Scripts/utils/completion.sh

# Make executable
cd ${GPUMDkit_path}
chmod +x gpumdkit.sh
```

### Basic Usage

**Interactive Mode** (recommended for beginners):
```bash
gpumdkit.sh
# Follow the menu to select functions
```

**Command-Line Mode** (for quick operations):
```bash
# Show help
gpumdkit.sh -h

# Plot thermodynamic data
gpumdkit.sh -plt thermo

# Convert VASP to extxyz
gpumdkit.sh -out2xyz ./

# Analyze composition
gpumdkit.sh -analyze_comp train.xyz
```

## Command-Line Reference

### Format Conversion
```bash
gpumdkit.sh -out2xyz <dir>              # OUTCAR to extxyz
gpumdkit.sh -pos2exyz <P> <xyz>         # POSCAR to extxyz
gpumdkit.sh -cif2exyz <cif>             # CIF to extxyz
gpumdkit.sh -lmp2exyz <dump> <elements> # LAMMPS to extxyz
gpumdkit.sh -addgroup <P> <elements>    # Add group labels
gpumdkit.sh -addweight <in> <out> <w>   # Add weights
```

### Analysis Tools
```bash
gpumdkit.sh -range <xyz> <property>     # Show range of property
gpumdkit.sh -min_dist <xyz>             # Calculate min distance
gpumdkit.sh -min_dist_pbc <xyz>         # Min distance with PBC
gpumdkit.sh -filter_value <xyz> <p> <t> # Filter by threshold
gpumdkit.sh -filter_dist <xyz> <dist>   # Filter by distance
gpumdkit.sh -analyze_comp <xyz>         # Analyze composition
```

### Plotting
```bash
gpumdkit.sh -plt thermo       # Plot thermodynamic properties
gpumdkit.sh -plt train        # Plot NEP training results
gpumdkit.sh -plt prediction   # Plot NEP predictions
gpumdkit.sh -plt msd          # Plot mean square displacement
gpumdkit.sh -plt rdf          # Plot radial distribution function
gpumdkit.sh -plt sdc          # Plot self-diffusion coefficient
gpumdkit.sh -plt charge       # Plot charge distribution
gpumdkit.sh -plt des <method> # Plot descriptors (umap/tsne/pca)
```

### Calculators
```bash
gpumdkit.sh -calc ionic-cond <elem> <charge>          # Ionic conductivity
gpumdkit.sh -calc nep <in> <out> <model>              # NEP predictions
gpumdkit.sh -calc des <method> <in> <out> <model> <elem>  # Descriptors
```

### Utilities
```bash
gpumdkit.sh -clean      # Clean extra files
gpumdkit.sh -update     # Update GPUMDkit
gpumdkit.sh -time gpumd # Estimate GPUMD time
gpumdkit.sh -time nep   # Estimate NEP time
```

## Quick Examples

### Example 1: Convert VASP Data to NEP Training Set

```bash
# 1. Convert VASP OUTCAR files to extxyz
gpumdkit.sh -out2xyz ./vasp_calculations/

# 2. Check data quality
gpumdkit.sh -range train.xyz force
gpumdkit.sh -min_dist_pbc train.xyz

# 3. Filter outliers
gpumdkit.sh -filter_value train.xyz force 30

# 4. Add group labels (for multi-element systems)
gpumdkit.sh -addgroup POSCAR Li Y Cl
```

### Example 2: Analyze NEP Training

```bash
# After NEP training completes:

# Plot training progress
gpumdkit.sh -plt train

# Plot prediction accuracy
gpumdkit.sh -plt prediction

# Plot force errors
gpumdkit.sh -plt force_errors

# Visualize learning rate
gpumdkit.sh -plt lr
```

### Example 3: Analyze MD Simulation

```bash
# After GPUMD simulation:

# Plot thermodynamic properties
gpumdkit.sh -plt thermo

# Plot diffusion properties
gpumdkit.sh -plt msd
gpumdkit.sh -plt sdc

# Calculate ionic conductivity
gpumdkit.sh -calc ionic-cond Li 1

# Plot structure
gpumdkit.sh -plt rdf
```

## Detailed Tutorials

For detailed tutorials on specific topics, see:

- [Format Conversion](format_conversion.md) - Converting between file formats
- [Sample Structures](sample_structures.md) - Sampling and selecting structures
- [Calculators](calculators.md) - Using calculator tools
- [Plot Scripts](plot_scripts.md) - Visualization and plotting
- [Workflow Development](workflow_dev.md) - Batch processing workflows
- [Active Learning](workflow_active_learning.md) - NEP active learning cycles
- [Custom Commands](custom_commands.md) - Creating custom commands

## Script Documentation

For comprehensive documentation of all scripts:

- [Scripts/plt_scripts/](../../Scripts/plt_scripts/README.md) - Plotting scripts
- [Scripts/analyzer/](../../Scripts/analyzer/README.md) - Analysis tools
- [Scripts/calculators/](../../Scripts/calculators/README.md) - Calculator scripts
- [Scripts/format_conversion/](../../Scripts/format_conversion/README.md) - Format converters
- [Scripts/sample_structures/](../../Scripts/sample_structures/README.md) - Sampling tools
- [Scripts/workflow/](../../Scripts/workflow/README.md) - Workflow automation
- [Scripts/utils/](../../Scripts/utils/README.md) - Utility scripts

## Tips

- Use **tab completion** for faster command entry
- Add `save` at the end of plot commands to save figures (e.g., `-plt thermo save`)
- Check help for any command: `gpumdkit.sh -<command> -h`
- Use interactive mode when learning, command-line mode for automation
- Keep training data organized in separate directories

## Getting Help

- **Documentation**: [https://github.com/zhyan0603/GPUMDkit](https://github.com/zhyan0603/GPUMDkit)
- **Issues**: [GitHub Issues](https://github.com/zhyan0603/GPUMDkit/issues)
- **Contact**: Zihan YAN (yanzihan@westlake.edu.cn)

---

Thank you for using GPUMDkit!
