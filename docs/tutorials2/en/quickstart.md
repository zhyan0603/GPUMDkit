<div align="center">
  <h1>Quick Start</h1>
  <p>
    <strong>English</strong> | <a href="../zh/quickstart.md">简体中文</a>
  </p>
</div>

Get up and running with GPUMDkit in 5 minutes.

## Prerequisites

- Linux or macOS (Windows via WSL)
- Python 3.8+ (3.12 recommended)
- Git

## Installation

### Step 1: Clone

```bash
git clone https://github.com/zhyan0603/GPUMDkit.git
cd GPUMDkit
```

### Step 2: Install

```bash
source ./install.sh
```

### Step 3: Reload Shell

```bash
source ~/.bashrc  # or source ~/.zshrc
```

### Step 4: Verify

```bash
gpumdkit.sh -h
```

## Python Dependencies

```bash
conda create -n gpumdkit python=3.12
conda activate gpumdkit
pip install neptrain ase pymatgen dpdata numpy scipy matplotlib
```

> Make sure the `gpumdkit` environment is activated before using GPUMDkit.

## Usage Modes

### Interactive Mode

```bash
gpumdkit.sh
```

```
           ____ ____  _   _ __  __ ____  _    _ _
          / ___|  _ \| | | |  \/  |  _ \| | _(_) |_
         | |  _| |_) | | | | |\/| | | | | |/ / | __|
         | |_| |  __/| |_| | |  | | |_| |   <| | |_
          \____|_|    \___/|_|  |_|____/|_|\_\_|\__|

          GPUMDkit Version 1.5.5 (dev) (2026-05-10)

 ---------------------- GPUMD ------------------------
  1) Format Conversion          2) Sample Structures
  3) Workflow                   4) Calculators
  5) Analyzer                   6) Visualization
  7) Utilities                  0) Exit
 ------------>>
```

### Command-Line Mode

```bash
gpumdkit.sh -h                    # Show help
gpumdkit.sh -pos2exyz POSCAR model.xyz  # Direct command
```

## First Examples

### Convert Structure

```bash
gpumdkit.sh -pos2exyz POSCAR model.xyz
gpumdkit.sh -addgroup POSCAR Li Y Cl
```

### Analyze Data

```bash
gpumdkit.sh -range train.xyz force
gpumdkit.sh -min_dist_pbc train.xyz
gpumdkit.sh -analyze_comp train.xyz
```

### Plot Results

```bash
gpumdkit.sh -plt train
gpumdkit.sh -plt msd
gpumdkit.sh -plt rdf
```

## Getting Help

```bash
gpumdkit.sh -h              # General help
gpumdkit.sh -plt -h         # Plot help
gpumdkit.sh -calc -h        # Calculator help
gpumdkit.sh -pos2exyz -h    # Specific command help
```

## Troubleshooting

**Command not found**
```bash
export GPUMDkit_path=/path/to/GPUMDkit
export PATH=${GPUMDkit_path}:${PATH}
```

**Missing packages**
```bash
conda activate gpumdkit
```

## Next Steps

- [Format Conversion](format_conversion.md)
- [Calculators](calculators.md)
- [Visualization](visualization.md)
- [NEP Training](nep_training.md)
