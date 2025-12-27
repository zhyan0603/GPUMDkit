# GPUMDkit Tutorials

Welcome to the **GPUMDkit** tutorials! This comprehensive guide will help you master GPUMDkit, a powerful toolkit for GPUMD (Graphics Processing Units Molecular Dynamics) and NEP (neuroevolution potential) workflows.

## Table of Contents

- [Introduction](#introduction)
- [Getting Started](#getting-started)
- [Quick Links to Script Documentation](#quick-links-to-script-documentation)
- [End-to-End Workflow Examples](#end-to-end-workflow-examples)
- [Interactive vs Command-Line Mode](#interactive-vs-command-line-mode)
- [Best Practices](#best-practices)
- [Contributing Tutorials](#contributing-tutorials)
- [Additional Resources](#additional-resources)

---

## Introduction

GPUMDkit is designed to streamline your computational materials research workflow by providing:

- **Unified Interface**: Access all tools through a single `gpumdkit.sh` command
- **Two Modes**: Interactive menu-driven mode OR fast command-line mode
- **Comprehensive Tools**: Format conversion, analysis, visualization, and automation
- **NEP Integration**: Complete support for NEP machine learning potential development
- **GPUMD Support**: Tools for setting up and analyzing molecular dynamics simulations

### Why Use GPUMDkit?

- ✅ **Save Time**: Automate repetitive tasks with workflow scripts
- ✅ **Reduce Errors**: Standardized procedures ensure consistency
- ✅ **Easy Learning**: Interactive mode guides you through options
- ✅ **Flexible**: Use command-line flags for scripting and automation
- ✅ **Well-Documented**: Comprehensive README files for every component

---

## Getting Started

### Installation

If you haven't installed GPUMDkit yet, see the [main README](../../README.md) for installation instructions.

**Quick setup:**
```bash
# Clone repository
git clone https://github.com/zhyan0603/GPUMDkit.git

# Set environment variable in ~/.bashrc
export GPUMDkit_path=/path/to/GPUMDkit
export PATH=${GPUMDkit_path}:${PATH}
source ${GPUMDkit_path}/Scripts/utils/completion.sh

# Make executable
cd ${GPUMDkit_path}
chmod +x gpumdkit.sh
```

### First Steps

1. **Check installation**:
   ```bash
   gpumdkit.sh -h
   ```

2. **Try interactive mode**:
   ```bash
   gpumdkit.sh
   ```

3. **Test a simple command**:
   ```bash
   gpumdkit.sh -plt -h
   ```

---

## Quick Links to Script Documentation

Each Scripts/ subdirectory has comprehensive documentation. Click to explore:

### 📊 [Plotting Scripts](../../Scripts/plt_scripts/README.md)
Visualization tools for NEP training, MD simulations, and materials properties.

**Key capabilities:**
- NEP training and prediction visualization
- Thermodynamic property evolution
- Transport properties (MSD, diffusion, conductivity)
- Structural analysis (RDF, charge distribution)
- Advanced analysis (descriptors, DOAS, force errors)

**Quick example:**
```bash
gpumdkit.sh -plt thermo        # Plot temperature, pressure, energy evolution
gpumdkit.sh -plt train         # Visualize NEP training progress
gpumdkit.sh -plt msd save      # Plot and save mean square displacement
```

### 🔍 [Analyzer Scripts](../../Scripts/analyzer/README.md)
Data analysis and quality control tools for structure files and datasets.

**Key capabilities:**
- Energy/force/virial range analysis
- Minimum distance calculations (with/without PBC)
- Structure filtering by various criteria
- Composition analysis
- Time estimation for running calculations

**Quick example:**
```bash
gpumdkit.sh -range train.xyz force     # Analyze force distribution
gpumdkit.sh -min_dist_pbc model.xyz    # Calculate min distance with PBC
gpumdkit.sh -filter_value train.xyz force 30  # Remove high-force outliers
```

### 🧮 [Calculator Scripts](../../Scripts/calculators/README.md)
Computational tools for property calculations using NEP models.

**Key capabilities:**
- Ionic conductivity calculations
- NEP model predictions
- Descriptor calculations for ML analysis
- Density of atomistic states (DOAS)
- NEB (migration barrier) calculations

**Quick example:**
```bash
gpumdkit.sh -calc ionic-cond Li 1              # Calculate Li⁺ conductivity
gpumdkit.sh -calc nep test.xyz pred.xyz nep.txt  # NEP predictions
gpumdkit.sh -calc des umap train.xyz desc.npy nep.txt Li  # Descriptors
```

### 🔄 [Format Conversion Scripts](../../Scripts/format_conversion/README.md)
Convert between various structure file formats.

**Key capabilities:**
- VASP (POSCAR/OUTCAR) ↔ extxyz
- CP2K, ABACUS output → extxyz
- LAMMPS dump ↔ extxyz
- CIF → POSCAR/extxyz
- Adding group labels and weights

**Quick example:**
```bash
gpumdkit.sh -out2xyz ./VASP_calcs/     # Convert OUTCAR files
gpumdkit.sh -pos2exyz POSCAR model.xyz # POSCAR to extxyz
gpumdkit.sh -addgroup POSCAR Li Y Cl   # Add group labels
```

### 📦 [Structure Sampling Scripts](../../Scripts/sample_structures/README.md)
Sample and select structures for NEP training.

**Key capabilities:**
- Uniform and random sampling
- FPS (Farthest Point Sampling) with PyNEP
- Structure perturbation for data augmentation
- Active learning structure selection

**Quick example:**
```bash
python sample_structures.py large.xyz uniform 1000  # Uniform sampling
gpumdkit.sh -pynep candidates.xyz train.xyz nep.txt  # FPS sampling
```

### 🔧 [Workflow Scripts](../../Scripts/workflow/README.md)
Automation for high-throughput calculations and active learning.

**Key capabilities:**
- Batch SCF preprocessing (VASP, CP2K)
- MD sampling batch preprocessing (GPUMD, LAMMPS)
- Active learning automation
- Job submission management

**Quick example:**
```bash
gpumdkit.sh
# Select: 3) Workflow
# Select: 301) SCF batch pretreatment
```

### 🛠️ [Utility Scripts](../../Scripts/utils/README.md)
Helper tools for maintenance and convenience.

**Key capabilities:**
- Cleaning working directories
- Updating GPUMDkit
- Bash completion support
- Atom renumbering for LAMMPS

**Quick example:**
```bash
gpumdkit.sh -clean    # Clean extra files from work directory
gpumdkit.sh -update   # Update GPUMDkit to latest version
```

### 📚 [Scripts Overview](../../Scripts/README.md)
Main documentation page for all Scripts/ with quick navigation.

---

## End-to-End Workflow Examples

### Workflow 1: Building a NEP Model from VASP Calculations

This complete example shows how to create a NEP model from DFT calculations.

#### Step 1: Prepare Training Data from VASP

```bash
# You have multiple directories with VASP calculations
# project/
# ├── struct_001/OUTCAR
# ├── struct_002/OUTCAR
# └── struct_003/OUTCAR

# Convert all OUTCAR files to extxyz
cd project/
gpumdkit.sh -out2xyz .

# Output: train.xyz with all structures combined

# Check the data quality
gpumdkit.sh -analyze_comp train.xyz
gpumdkit.sh -range train.xyz energy
gpumdkit.sh -range train.xyz force hist

# Filter problematic structures
gpumdkit.sh -filter_value train.xyz force 30
gpumdkit.sh -min_dist_pbc filtered_force.xyz
```

#### Step 2: Add Group Labels for NEP

```bash
# NEP requires group labels to identify atom types
# Assuming LiYCl system with POSCAR from one calculation

gpumdkit.sh -addgroup struct_001/POSCAR Li Y Cl

# This creates model.xyz with proper group labels
# Use this as a template for group definition
```

#### Step 3: Prepare Training and Test Sets

```bash
# Split filtered data into training and test sets
# Use 80-20 split

total_structs=$(grep -c "Lattice" filtered_force.xyz)
train_size=$((total_structs * 80 / 100))

# Sample for training set
python ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py \
    filtered_force.xyz uniform $train_size

mv sampled_structures.xyz train.xyz

# Remaining for test set (requires custom script or manual extraction)
# For simplicity, use all data for training in first iteration
```

#### Step 4: Train Initial NEP Model

```bash
# Create nep.in file
cat > nep.in << EOF
version        4
type           3 Li Y Cl
cutoff         8 4
n_max          8 4
basis_size     12 8
l_max          4 2 0
neuron         100
batch          1000
generation     1000000
population     50
lambda_1       0.1
lambda_2       0.01
lambda_e       1.0
lambda_f       1.0
lambda_v       0.1
EOF

# Run NEP training (requires nep executable from GPUMD)
nep
```

#### Step 5: Analyze Training Results

```bash
# Visualize training progress
gpumdkit.sh -plt train save
gpumdkit.sh -plt train_test save
gpumdkit.sh -plt force_errors save

# Check learning rate evolution
gpumdkit.sh -plt lr save

# Analyze model parameters
gpumdkit.sh -plt restart save
```

#### Step 6: Validate Model with MD Simulations

```bash
# Create test structure (if not already have one)
gpumdkit.sh -pos2exyz struct_001/POSCAR test_model.xyz
gpumdkit.sh -addgroup struct_001/POSCAR Li Y Cl

# Create run.in for GPUMD
cat > run.in << EOF
velocity        800
time_step       1
ensemble        npt_ber 800 800 100 0 0 0 1000
dump_thermo     1000
dump_position   1000
compute_msd     100 1 0
run             100000
EOF

# Run GPUMD simulation
gpumd

# Analyze results
gpumdkit.sh -plt thermo save
gpumdkit.sh -plt msd save
gpumdkit.sh -calc ionic-cond Li 1
```

---

### Workflow 2: Active Learning Cycle for NEP Improvement

Iteratively improve your NEP model using active learning.

#### Overview

```
Initial Model → MD Sampling → Uncertainty → DFT → Retrain → Improved Model
                     ↑____________________________________________|
```

#### Step 1: Generate Candidate Structures with MD

```bash
# Use current NEP model for MD sampling at various conditions

# High temperature exploration
cat > run_explore.in << EOF
velocity        1200
time_step       1
ensemble        nvt_ber 1200 100
dump_position   100
run             50000
EOF

gpumd
mv dump.xyz candidates_highT.xyz

# Multiple temperatures, pressures, etc.
# Combine all candidate structures
cat candidates_*.xyz > all_candidates.xyz
```

#### Step 2: Select Uncertain Structures

```bash
# Use PyNEP FPS for diversity
gpumdkit.sh -pynep all_candidates.xyz diverse_candidates.xyz nep.txt

# Or use GPUMD's active learning mode
# In run.in, add:
# active 100 0.15  # Check every 100 steps, threshold 0.15 eV/Å

gpumd  # With active mode enabled

# Select structures with high force uncertainty
python ${GPUMDkit_path}/Scripts/sample_structures/select_max_modev.py 50 0.15

# Output: selected_structures.xyz with 50 most uncertain structures
```

#### Step 3: Prepare DFT Calculations

```bash
# Use batch workflow for DFT
# First convert to POSCAR files if needed
gpumdkit.sh -exyz2pos selected_structures.xyz

# Then use SCF batch workflow
gpumdkit.sh
# Select: 3) Workflow
# Select: 301) SCF batch pretreatment
# Enter prefix: iter02_active

# Prepare VASP inputs
mkdir fp
# Copy POTCAR, INCAR, KPOINTS to fp/

# Submit all calculations
for dir in struct_fp_iter02_active_*/; do
    cd $dir
    sbatch submit.sh
    cd ..
done
```

#### Step 4: Collect Results and Retrain

```bash
# After DFT calculations complete
gpumdkit.sh -out2xyz ./struct_fp_iter02_active_*/

# Merge with previous training data
cat train_iter01.xyz train.xyz > train_iter02.xyz

# Check data quality
gpumdkit.sh -analyze_comp train_iter02.xyz
gpumdkit.sh -range train_iter02.xyz force hist

# Retrain NEP
cp train_iter02.xyz train.xyz
nep  # Using same or updated nep.in

# Compare models
gpumdkit.sh -plt train save
# Examine if RMSE decreased
```

#### Step 5: Validate and Iterate

```bash
# Test improved model
gpumdkit.sh -calc nep test.xyz predictions_iter02.xyz nep.txt
gpumdkit.sh -plt prediction save

# Compare metrics:
# - Training RMSE (lower is better)
# - Test RMSE (should also decrease)
# - Force error metrics

# If not converged, repeat cycle:
# - Generate new candidates with improved model
# - Select uncertain structures
# - Run DFT
# - Retrain

# Stop when:
# - RMSE plateaus
# - No structures exceed uncertainty threshold
# - Satisfactory accuracy achieved
```

---

### Workflow 3: High-Throughput Structure Screening

Screen many structures for specific properties.

#### Step 1: Prepare Structure Library

```bash
# Download CIF files from materials databases
# (e.g., Materials Project, ICSD)

# Convert all CIF to extxyz
for cif in *.cif; do
    name=$(basename "$cif" .cif)
    gpumdkit.sh -cif2exyz "$cif"
    mv *.xyz "${name}.xyz"
done

# Combine if needed
cat *.xyz > structure_library.xyz
```

#### Step 2: Calculate Properties with NEP

```bash
# Use trained NEP model for fast screening
gpumdkit.sh -calc nep structure_library.xyz nep_results.xyz nep.txt

# Extract energies
grep "energy=" nep_results.xyz | sed 's/.*energy=\([^ ]*\).*/\1/' > energies.dat

# Screen for stability (example: formation energy)
# (Requires reference energies and composition analysis)
```

#### Step 3: Analyze Structural Features

```bash
# Calculate descriptors for all structures
gpumdkit.sh -calc des umap structure_library.xyz descriptors.npy nep.txt Li

# Visualize descriptor space
gpumdkit.sh -plt des umap descriptors.npy save

# Find similar structures
# (Use descriptor distances to identify structural motifs)
```

#### Step 4: MD Simulations for Selected Structures

```bash
# Identify promising candidates (low energy, interesting features)
# Run MD to test stability and dynamics

# Batch MD setup
gpumdkit.sh
# Select: 3) Workflow
# Select: 302) MD batch (GPUMD)

# After simulations, analyze results
for dir in struct_md_*/; do
    cd $dir
    echo "Analyzing $dir"
    gpumdkit.sh -plt thermo
    gpumdkit.sh -plt msd
    cd ..
done
```

---

## Interactive vs Command-Line Mode

### Interactive Mode

**Best for:**
- Learning GPUMDkit features
- Infrequent, complex tasks
- Guided workflows

**How to use:**
```bash
gpumdkit.sh
```

**Example session:**
```
 ----------------------- GPUMD -----------------------
 1) Format Conversion          2) Sample Structures
 3) Workflow                   4) Calculators
 5) Analyzer                   6) Developing ...
 0) Quit!
 ------------>>
 Input the function number: 1

 ============ Format Conversion ============
 101) OUTCAR to XYZ          102) MTP to XYZ
 103) CP2K to XYZ            104) ABACUS to XYZ
 105) extXYZ to POSCAR       ...
 Input the function number: 101

 # Follow prompts...
```

### Command-Line Mode

**Best for:**
- Quick one-off commands
- Scripting and automation
- Batch processing

**How to use:**
```bash
gpumdkit.sh -<flag> [arguments]
```

**Common patterns:**
```bash
# Single command
gpumdkit.sh -plt thermo save

# In a script
#!/bin/bash
for temp in 300 600 900 1200; do
    cd temp_${temp}
    gpumdkit.sh -plt thermo save
    gpumdkit.sh -plt msd save
    cd ..
done

# Piped with other commands
gpumdkit.sh -range train.xyz force | grep "range:"
```

### Choosing the Right Mode

| Task | Recommended Mode | Why |
|------|------------------|-----|
| First time using a feature | Interactive | Guided prompts |
| Quick visualization | Command-line | Fast, no menus |
| Batch processing | Command-line | Scriptable |
| Complex workflow | Interactive | Step-by-step guidance |
| Automation | Command-line | Integration with scripts |
| Exploring options | Interactive | See all choices |

---

## Best Practices

### Data Management

1. **Organize by project**:
   ```
   project_name/
   ├── 00_structures/      # Initial structures
   ├── 01_dft_calcs/       # DFT calculations
   ├── 02_training_data/   # Processed training sets
   ├── 03_nep_models/      # NEP models by iteration
   ├── 04_validation/      # Validation simulations
   └── 05_analysis/        # Analysis results and figures
   ```

2. **Version control for models**:
   ```bash
   nep_v1.0_initial/
   ├── nep.txt
   ├── nep.in
   ├── train.xyz
   └── results/
   
   nep_v1.1_iter01/
   ├── nep.txt
   ├── nep.in
   ├── train.xyz
   └── results/
   ```

3. **Keep logs**:
   ```bash
   # Document your workflow
   cat > workflow.log << EOF
   Date: 2025-01-15
   Task: Initial NEP training
   Training set: 500 structures from VASP
   Command: nep
   Results: RMSE_E=5.2 meV/atom, RMSE_F=45 meV/Å
   EOF
   ```

### Quality Control

1. **Always check data before training**:
   ```bash
   gpumdkit.sh -range train.xyz force hist
   gpumdkit.sh -analyze_comp train.xyz
   gpumdkit.sh -min_dist_pbc train.xyz
   ```

2. **Validate conversions**:
   ```bash
   # After format conversion, spot check
   ase gui converted_structure.xyz
   # Verify lattice, atoms, energy/forces
   ```

3. **Monitor training**:
   ```bash
   # During NEP training, periodically check
   gpumdkit.sh -plt train
   gpumdkit.sh -time nep
   ```

### Performance Tips

1. **Use command-line for batch operations**:
   ```bash
   # Process multiple files efficiently
   for xyz in *.xyz; do
       gpumdkit.sh -range "$xyz" force
   done
   ```

2. **Parallel processing when possible**:
   ```bash
   # GNU parallel for independent tasks
   parallel gpumdkit.sh -plt thermo save ::: temp_*/
   ```

3. **Filter early**:
   ```bash
   # Remove outliers before expensive operations
   gpumdkit.sh -filter_value train.xyz force 30
   # Then use filtered_force.xyz for training
   ```

### Documentation

1. **Comment your workflow scripts**:
   ```bash
   #!/bin/bash
   # NEP Model Development Pipeline
   # Author: Your Name
   # Date: 2025-01-15
   # Purpose: Automated NEP training with active learning
   
   # Step 1: Prepare data
   gpumdkit.sh -out2xyz ./dft_results/
   ...
   ```

2. **Track parameters**:
   ```bash
   # Save nep.in with each model version
   cp nep.in nep_models/v1.2/nep.in
   
   # Document sampling parameters
   echo "Sampled 1000 structures uniformly" > sampling.log
   ```

---

## Contributing Tutorials

We welcome contributions! Help make GPUMDkit better for everyone.

### How to Contribute

1. **Share your workflows**:
   - Document a workflow you've developed
   - Explain the problem it solves
   - Provide example commands and results

2. **Improve existing tutorials**:
   - Fix errors or unclear explanations
   - Add more detailed examples
   - Update for new features

3. **Create specialized guides**:
   - System-specific guides (e.g., "NEP for superionic conductors")
   - Code-specific integration (e.g., "Using GPUMDkit with Quantum ESPRESSO")
   - Advanced techniques (e.g., "Optimizing NEP hyperparameters")

### Contribution Process

1. **Fork the repository**:
   ```bash
   # On GitHub, click "Fork"
   git clone https://github.com/YOUR_USERNAME/GPUMDkit.git
   ```

2. **Create your tutorial**:
   ```bash
   cd GPUMDkit/docs/tutorials/
   # Add your tutorial as a new markdown file
   vi my_awesome_workflow.md
   ```

3. **Follow the template**:
   ```markdown
   # Tutorial Title
   
   ## Overview
   Brief description
   
   ## Prerequisites
   What you need to know/have
   
   ## Step-by-Step Instructions
   Detailed steps with commands
   
   ## Expected Results
   What should happen
   
   ## Troubleshooting
   Common issues
   
   ## References
   Links to papers, docs
   ```

4. **Test your tutorial**:
   - Follow your own instructions
   - Ensure commands work
   - Verify results are reproducible

5. **Submit pull request**:
   ```bash
   git add docs/tutorials/my_awesome_workflow.md
   git commit -m "Add tutorial: My Awesome Workflow"
   git push origin main
   # Create pull request on GitHub
   ```

### Tutorial Guidelines

- **Clear objectives**: State what the tutorial achieves
- **Complete examples**: Include all necessary commands
- **Explain concepts**: Don't just list commands, explain why
- **Realistic data**: Use representative examples
- **Troubleshooting**: Anticipate and address common issues
- **References**: Cite papers, documentation, external resources

### Topics We'd Love to See

- **Materials-specific workflows**: Different material classes
- **Advanced analysis**: Novel analysis techniques
- **Integration guides**: Using GPUMDkit with other codes
- **Optimization strategies**: Best practices for specific scenarios
- **Case studies**: Real research examples
- **Video tutorials**: Screen recordings demonstrating workflows

---

## Additional Resources

### Documentation

- **Main README**: [../../README.md](../../README.md) - Installation and overview
- **Contributing Guide**: [../../CONTRIBUTING.md](../../CONTRIBUTING.md) - Development guidelines
- **Scripts Documentation**: [../../Scripts/README.md](../../Scripts/README.md) - All scripts overview

### External Resources

- **GPUMD Documentation**: [https://gpumd.org](https://gpumd.org)
- **NEP Formalism**: [https://gpumd.org/potentials/nep.html](https://gpumd.org/potentials/nep.html)
- **ASE Documentation**: [https://wiki.fysik.dtu.dk/ase/](https://wiki.fysik.dtu.dk/ase/)
- **Materials Project**: [https://materialsproject.org](https://materialsproject.org)

### Related Tutorials

Existing tutorial files in this directory:

- [Format Conversion](format_conversion.md) - Detailed format conversion guide
- [Sample Structures](sample_structures.md) - Structure sampling guide
- [Workflow Dev](workflow_dev.md) - Workflow development guide
- [Calculators](calculators.md) - Calculator tools guide
- [Plot Scripts](plot_scripts.md) - Visualization guide
- [Active Learning](workflow_active_learning.md) - Active learning tutorial
- [Custom Commands](custom_commands.md) - Creating custom commands

### Getting Help

- **GitHub Issues**: [https://github.com/zhyan0603/GPUMDkit/issues](https://github.com/zhyan0603/GPUMDkit/issues)
- **Email**: Zihan YAN (yanzihan@westlake.edu.cn)
- **Documentation**: [https://zhyan0603.github.io/GPUMDkit](https://zhyan0603.github.io/GPUMDkit)

### Community

- Share your GPUMDkit workflows on GitHub
- Contribute to documentation improvements
- Report bugs and request features
- Help other users on GitHub Discussions (if enabled)

---

## Quick Reference Card

### Essential Commands

```bash
# Help
gpumdkit.sh -h                      # General help
gpumdkit.sh -plt -h                 # Plotting help

# Format Conversion
gpumdkit.sh -out2xyz <dir>          # VASP to extxyz
gpumdkit.sh -pos2exyz <P> <xyz>     # POSCAR to extxyz

# Analysis
gpumdkit.sh -range <xyz> force      # Force range
gpumdkit.sh -analyze_comp <xyz>     # Composition

# Plotting
gpumdkit.sh -plt thermo save        # Plot thermodynamics
gpumdkit.sh -plt train save         # NEP training progress

# Calculators
gpumdkit.sh -calc ionic-cond Li 1   # Ionic conductivity

# Utilities
gpumdkit.sh -clean                  # Clean directory
gpumdkit.sh -update                 # Update GPUMDkit
```

### File Formats

| Format | Extension | Usage |
|--------|-----------|-------|
| Extended XYZ | `.xyz` | GPUMD, NEP training |
| VASP POSCAR | `POSCAR` | VASP input structure |
| VASP OUTCAR | `OUTCAR` | VASP output data |
| LAMMPS dump | `.lammps`, `.dump` | LAMMPS trajectory |
| CIF | `.cif` | Crystallographic data |
| NEP model | `nep.txt` | Trained NEP model |

---

## Conclusion

This tutorial guide has covered:
- ✅ Overview of GPUMDkit capabilities
- ✅ Quick links to all script documentation
- ✅ Three comprehensive end-to-end workflows
- ✅ Interactive vs command-line mode guidance
- ✅ Best practices for data management and quality control
- ✅ How to contribute tutorials

**Next Steps:**

1. **Beginners**: Start with interactive mode to explore features
2. **Intermediate**: Try the end-to-end workflows with your data
3. **Advanced**: Develop custom workflows and contribute back

**Remember:**
- Documentation is your friend - each Scripts/ subdirectory has detailed README
- Start simple, build complexity gradually
- Test workflows on small datasets first
- Keep good records of your procedures
- Don't hesitate to ask for help!

---

Thank you for using GPUMDkit! We're excited to see what you'll discover and create. If you have questions, suggestions, or success stories to share, please reach out via [GitHub Issues](https://github.com/zhyan0603/GPUMDkit/issues) or email Zihan YAN (yanzihan@westlake.edu.cn).

**Happy researching!** 🚀

## Introduction

`GPUMDkit` offers two main modes of operation:

1. **Interactive Mode**: Run `gpumdkit.sh` and follow the menu prompts for a guided experience.
2. **Command-Line Mode**: Directly pass arguments to `gpumdkit.sh` for quick and streamlined command execution.

## Interactive Mode

### Getting Started

1. Open your terminal.

2. Execute the `gpumdkit.sh` script:
    ```sh
    ./gpumdkit.sh
    ```
    
3. Follow the on-screen prompts to interactively select and run the desired script.

    ```
             ____ ____  _   _ __  __ ____  _    _ _
            / ___|  _ \| | | |  \/  |  _ \| | _(_) |_
           | |  _| |_) | | | | |\/| | | | | |/ / | __|
           | |_| |  __/| |_| | |  | | |_| |   <| | |_
            \____|_|    \___/|_|  |_|____/|_|\_\_|\__|
    
            GPUMDkit Version 1.3.2 (dev) (2025-07-21)
         Developer: Zihan YAN (yanzihan@westlake.edu.cn)
    
     ----------------------- GPUMD -----------------------
     1) Format Conversion          2) Sample Structures
     3) Workflow (dev)             4) Calculators
     5) Developing ...             6) Developing ...
     0) Quit!
     ------------>>
     Input the function number:
    ```

    

This mode is useful for new users or for tasks that require step-by-step guidance.

## Command-Line Mode

### Quick Commands

For users familiar with the `GPUMDkit` , the command-line mode allows for faster execution by directly passing arguments to `gpumdkit.sh`. Here are some examples:

#### Example 1: View help information

```
gpumdkit.sh -h
```

the help information:

```
+==================================================================================================+
|                              GPUMDkit 1.3.0 (dev) (2025-07-11) Usage                             |
|                                                                 --- by Zihan YAN                 |
+======================================== Conversions =============================================+
| -outcar2exyz   Convert OUTCAR to extxyz       | -pos2exyz     Convert POSCAR to extxyz           |
| -castep2exyz   Convert castep to extxyz       | -pos2lmp      Convert POSCAR to LAMMPS           |
| -cp2k2exyz     Convert cp2k output to extxyz  | -lmp2exyz     Convert LAMMPS-dump to extxyz      |
| -addgroup      Add group label                | -addweight    Add weight to the struct in extxyz |
| Developing...                                 | Developing...                                    |
+========================================= Analysis ===============================================+
| -range         Print range of energy etc.     | -max_rmse     Get max RMSE from XYZ              |
| -min_dist      Get min_dist between atoms     | -min_dist_pbc Get min_dist considering PBC       |
| -filter_box    Filter struct by box limits    | -filter_value Filter struct by value (efs)       |
| -filter_dist   Filter struct by min_dist      | Developing...                                    |
+=========================================    Misc  ==============+================================+
| -plt           Plot scripts                   | -get_frame     Extract the specified frame       |
| -calc          Calculators                    | -clear_xyz     Clear extra info in XYZ file      |
| -clean         Clear files for work_dir       | -time          Time consuming Analyzer           |
| Developing...                                 | Developing...                                    |
+==================================================================================================+
| For detailed usage and examples, use: gpumdkit.sh -<option> -h                                   |
+==================================================================================================+
```

#### Example 2: Convert VASP OUTCARs to extxyz
To convert a `VASP` `OUTCARs` to an extended XYZ format (`extxyz`) file, use the following command:
```sh
gpumdkit.sh -outcar2exyz <dir_of_OUTCARs>
gpumdkit.sh -outcar2exyz .
```

#### Example 3: Plot thermo evolution

To visualize `thermo` evolution from `thermo.out` :

```sh
gpumdkit.sh -plt thermo
```



## Detailed Tutorials

For more detailed tutorials on specific functionalities, refer to the following documents:

1. [Format Conversion](format_conversion.md): Detailed guide on `1) Format Conversion`.
2. [Sample Structures](sample_structures.md): Detailed guide on `2) Sample Structures`.
3. [Workflow (dev)](workflow_dev.md): Detailed guide on `3) Workflow (dev)`.
4. [Calculators](calculators.md): Detailed guide on `4) Calculators`.
5. [Active Learning](workflow_active_learning.md): Detailed guide on `workflow_active_learning_dev.sh`.


---

Thank you for using `GPUMDkit`! If you have any questions or need further assistance, feel free to open an issue on our GitHub repository or contact Zihan YAN (yanzihan@westlake.edu.cn).
