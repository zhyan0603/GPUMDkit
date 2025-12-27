# Utility Scripts

This directory contains utility scripts for GPUMDkit maintenance, convenience features, and helper functions.

## Overview

The utilities in this folder help with GPUMDkit's operation and user experience, including:
- Cleaning temporary files from working directories
- Updating GPUMDkit to the latest version
- Bash command completion for improved CLI experience
- Atom renumbering for LAMMPS dump files

---

## Scripts

### clean_extra_files.sh

Cleans the current working directory by removing temporary and output files while preserving essential input files.

#### Purpose
After running GPUMD or NEP training, many intermediate and output files are generated. This script helps maintain a clean workspace by removing these files while keeping important inputs.

#### Files Preserved
- `run.in` - GPUMD run configuration
- `nep.in` - NEP training configuration
- `model.xyz` - Structure file
- `nep.txt` - NEP model file
- `train.xyz` - Training dataset
- `test.xyz` - Testing dataset
- All files matching patterns: `*sub*`, `*.sh`, `*slurm` (submission scripts)

#### Usage

**Direct execution:**
```bash
bash clean_extra_files.sh
```

**Command-line mode:**
```bash
gpumdkit.sh -clean
```

#### Example Output
```
Cleaning working directory...
Keeping: run.in, nep.in, model.xyz, train.xyz, test.xyz
Removing: thermo.out, loss.out, energy_train.out, force_train.out, ...
Cleanup complete!
```

#### Use Cases
- After completing a training run
- Before archiving a project directory
- When disk space is limited
- To prepare for a fresh calculation

---

### update_gpumdkit.sh

Updates GPUMDkit to the latest version from the GitHub repository.

#### Purpose
Automatically pulls the latest changes from the GPUMDkit repository, ensuring you have access to new features, bug fixes, and improvements.

#### Requirements
- Git must be installed
- GPUMDkit must be in a Git repository
- Internet connection to GitHub
- `GPUMDkit_path` environment variable must be set

#### Usage

**Command-line mode:**
```bash
gpumdkit.sh -update
```

**Alternative:**
```bash
gpumdkit.sh -U
```

#### Example Output
```
Checking for updates...
Current version: 1.4.1
Latest version: 1.4.2
Pulling latest changes from GitHub...
Successfully updated to version 1.4.2
```

#### Behavior
1. Verifies you are in a Git repository
2. Checks current branch and remote status
3. Pulls latest changes from origin
4. Displays update information from `docs/updates.info`
5. Reports any conflicts or errors

#### Manual Update Alternative
If automatic update fails or you don't have GitHub access:
```bash
cd $GPUMDkit_path
git pull origin main
```

Or download manually:
```bash
wget https://github.com/zhyan0603/GPUMDkit/archive/refs/heads/main.zip
unzip main.zip
```

---

### completion.sh

Provides Bash command completion for GPUMDkit commands.

#### Purpose
Enables tab-completion for `gpumdkit.sh` commands, making the CLI more user-friendly and reducing typing errors.

#### Features
- Completes main command flags (e.g., `-plt`, `-calc`, `-h`)
- Completes subcommands (e.g., `thermo`, `train` after `-plt`)
- Completes file paths and options
- Context-aware suggestions

#### Setup

This should already be configured if you followed the installation instructions. It's sourced in your `~/.bashrc`:

```bash
source ${GPUMDkit_path}/Scripts/utils/completion.sh
```

#### Usage

Type `gpumdkit.sh -` and press **Tab** to see available options:
```bash
$ gpumdkit.sh -<Tab>
-calc       -clean      -filter_box  -h          -min_dist    -plt
-range      -time       -update      ...
```

Type `gpumdkit.sh -plt ` and press **Tab** to see plotting options:
```bash
$ gpumdkit.sh -plt <Tab>
thermo      train       prediction   msd         rdf         sdc
dimer       charge      doas         ...
```

#### Benefits
- **Faster command entry**: Reduce typing with tab completion
- **Discover commands**: See available options without consulting docs
- **Avoid typos**: Completion ensures correct command syntax
- **Learn as you go**: Explore GPUMDkit features interactively

---

### renumber_atoms.py

Renumbers atom IDs in LAMMPS dump files sequentially.

#### Purpose
LAMMPS dump files may have non-sequential or discontinuous atom IDs. This script renumbers atoms sequentially starting from 1, which is useful for certain analysis tools or visualization software.

#### Usage

```bash
python renumber_atoms.py <input_file> <output_file>
```

#### Parameters
- `<input_file>`: Path to the input LAMMPS dump file
- `<output_file>`: Path for the output renumbered dump file

#### Example

```bash
python renumber_atoms.py dump.lammps dump_renumbered.lammps
```

#### Output Format
The script:
- Preserves all dump file formatting
- Maintains atom coordinates and properties
- Only modifies atom ID column
- Shows progress bar for large files (using `tqdm`)

#### Progress Display
```
Processing atoms: |████████████████████| 100% [50000/50000]
Renumbering complete!
```

#### Use Cases
- Preparing LAMMPS output for visualization tools requiring sequential IDs
- Post-processing dumps after atom deletion/addition
- Standardizing atom numbering across multiple trajectories
- Debugging issues related to atom indexing

#### Requirements
- Python 3.x
- `tqdm` package (for progress bar)

Install dependencies:
```bash
pip install tqdm
```

#### Author
Dian HUANG (huangdian@stu.xjtu.edu.cn)

---

## Integration with gpumdkit.sh

Most utilities can be accessed through the main `gpumdkit.sh` interface:

| Script | Command | Description |
|--------|---------|-------------|
| `clean_extra_files.sh` | `gpumdkit.sh -clean` | Clean working directory |
| `update_gpumdkit.sh` | `gpumdkit.sh -update` or `-U` | Update GPUMDkit |
| `completion.sh` | Auto-loaded | Enable tab completion |
| `renumber_atoms.py` | Direct call only | Renumber LAMMPS atoms |

## Best Practices

1. **Regular updates**: Run `gpumdkit.sh -update` periodically to get latest features
2. **Clean workspace**: Use `-clean` after completing calculations to save disk space
3. **Use completion**: Take advantage of tab completion to explore available commands
4. **Backup important data**: Run `-clean` carefully - review what will be deleted first

## Contributing Utilities

To add new utility scripts:

1. **Keep utilities general-purpose**: Should benefit multiple workflows
2. **Add proper documentation**: Include usage examples and parameter descriptions
3. **Follow naming conventions**: Use descriptive snake_case names
4. **Update this README**: Document your new utility
5. **Consider integration**: Decide if it should be accessible via `gpumdkit.sh`

## Troubleshooting

### Update fails
**Problem**: `git pull` fails or reports conflicts
**Solution**: 
```bash
cd $GPUMDkit_path
git stash  # Save local changes
git pull origin main
git stash pop  # Restore local changes
```

### Completion not working
**Problem**: Tab completion doesn't show GPUMDkit commands
**Solution**:
1. Verify `completion.sh` is sourced in `~/.bashrc`
2. Reload bash configuration: `source ~/.bashrc`
3. Check `GPUMDkit_path` is set correctly: `echo $GPUMDkit_path`

### Clean removes important files
**Problem**: Accidentally deleted necessary files
**Solution**:
- Review `keep_files` array in `clean_extra_files.sh` before running
- Consider modifying the script to preserve additional files
- Keep backups of critical data before cleaning

---

Thank you for using GPUMDkit! If you have suggestions for new utilities or encounter issues, please open an issue on our [GitHub repository](https://github.com/zhyan0603/GPUMDkit/issues) or contact Zihan YAN (yanzihan@westlake.edu.cn).
