# Quick Start Guide

<div align="center">
  <p>
    <a href="../zh/quickstart.md">中文</a> | <strong>English</strong>
  </p>
</div>

This guide will help you get started with GPUMDkit quickly.

## Prerequisites

Before installing GPUMDkit, ensure you have:

- **Linux or macOS** (Windows via WSL)
- **Bash** shell
- **Python 3.8+** (3.12 recommended)
- **Git** (for cloning the repository)

## Installation

### Step 1: Clone the Repository

```bash
git clone https://github.com/zhyan0603/GPUMDkit.git
cd GPUMDkit
```

### Step 2: Run Installation Script

```bash
source ./install.sh
```

This script will:
1. Add GPUMDkit to your PATH
2. Set up bash completion
3. Configure environment variables

### Step 3: Reload Shell

```bash
source ~/.bashrc
# or
source ~/.zshrc
```

### Step 4: Verify Installation

```bash
gpumdkit.sh -h
```

You should see the help message with available commands.

## Python Dependencies

Some features require Python packages. Create a conda environment:

```bash
# Create environment
conda create -n gpumdkit python=3.12
conda activate gpumdkit

# Install core packages
pip install numpy ase scipy matplotlib

# Install optional packages (for full functionality)
pip install neptrain pymatgen dpdata

# For perovskite analysis (optional)
pip3 install git+https://github.com/MoseyQAQ/ferrodispcalc.git
```

## Two Usage Modes

### Interactive Mode

Launch the interactive menu:

```bash
gpumdkit.sh
```

You'll see:

```
       ____ ____  _   _ __  __ ____  _    _ _
      / ___|  _ \| | | |  \/  |  _ \| | _(_) |_
     | |  _| |_) | | | | |\/| | | | | |/ / | __|
     | |_| |  __/| |_| | |  | | |_| |   <| | |_
      \____|_|    \___/|_|  |_|____/|_|\_\_|\__|

      GPUMDkit Version 1.5.5 (dev) (2026-05-10)
Core Developer: Zihan YAN (yanzihan@westlake.edu.cn)

 ---------------------- GPUMD ------------------------
 1) Format Conversion          2) Sample Structures
 3) Workflow                   4) Calculators
 5) Analyzer                   6) Visualization
 7) Utilities                  8) Developing...
 0) Exit
 ------------>>
 Input the function number:
```

Navigate by entering the number of your desired function.

### Command-Line Mode

Execute commands directly:

```bash
# Show help
gpumdkit.sh -h

# Convert POSCAR to extxyz
gpumdkit.sh -pos2exyz POSCAR model.xyz

# Plot training results
gpumdkit.sh -plt train

# Calculate ionic conductivity
gpumdkit.sh -calc ionic-cond Li 1
```

## First Steps Example

### Example 1: Convert a Structure

```bash
# Convert POSCAR to extxyz format
gpumdkit.sh -pos2exyz POSCAR model.xyz

# Add group labels (required for NEP training)
gpumdkit.sh -addgroup POSCAR Li Y Cl
```

### Example 2: Analyze Training Data

```bash
# Check force range in training data
gpumdkit.sh -range train.xyz force

# Check minimum interatomic distances
gpumdkit.sh -min_dist_pbc train.xyz

# Analyze composition
gpumdkit.sh -analyze_comp train.xyz
```

### Example 3: Plot Results

```bash
# Plot NEP training results
gpumdkit.sh -plt train

# Plot mean square displacement
gpumdkit.sh -plt msd

# Plot radial distribution function
gpumdkit.sh -plt rdf
```

## Getting Help

### General Help

```bash
gpumdkit.sh -h
```

### Module-Specific Help

```bash
# Plot help
gpumdkit.sh -plt -h

# Calculator help
gpumdkit.sh -calc -h

# Specific command help
gpumdkit.sh -pos2exyz -h
gpumdkit.sh -calc ionic-cond -h
```

### Interactive Mode Help

In interactive mode, each submenu displays available options with descriptions.

## Common Issues

### Issue: Command not found

**Solution**: Ensure GPUMDkit is in your PATH:

```bash
# Check if gpumdkit.sh is in PATH
which gpumdkit.sh

# If not found, add to PATH manually
export GPUMDkit_path=/path/to/GPUMDkit
export PATH=${GPUMDkit_path}:${PATH}
```

### Issue: Missing Python packages

**Solution**: Activate the conda environment:

```bash
conda activate gpumdkit
```

### Issue: Permission denied

**Solution**: Make the script executable:

```bash
chmod +x /path/to/GPUMDkit/gpumdkit.sh
```

## Next Steps

Now that you have GPUMDkit installed, explore these tutorials:

- [Format Conversion](format_conversion.md) - Learn to convert between file formats
- [Calculators](calculators.md) - Compute material properties
- [Visualization](visualization.md) - Create publication-ready plots
- [NEP Training Guide](nep_training.md) - Complete training workflow

## Updating GPUMDkit

To update to the latest version:

```bash
gpumdkit.sh -update
```

Or manually:

```bash
cd /path/to/GPUMDkit
git pull origin main
```
