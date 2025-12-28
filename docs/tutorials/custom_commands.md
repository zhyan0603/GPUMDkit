# Custom Commands in GPUMDkit

GPUMDkit allows users to extend the command-line interface with **custom commands**. This powerful feature lets you define shortcuts for frequently used operations, create custom workflows, and automate repetitive tasks by editing your `~/.gpumdkit.in` file.

## Overview

Custom commands enable you to:
- 🔧 Create shortcuts for complex command sequences
- ⚡ Automate repetitive workflows
- 🎯 Build custom analysis pipelines
- 📦 Package frequently used operations
- 🔄 Share workflows with collaborators

All custom functions must start with the `custom_` prefix to be recognized by GPUMDkit.

---

## Setup

Create or edit `~/.gpumdkit.in`:
```bash
vi ~/.gpumdkit.in
```

Define custom functions with the `custom_` prefix.

## Examples

### Simple Greeting

Write the following codes in your `~/.gpumdkit.in` 

```
custom_hello() {
    echo "Hello, GPUMDkit user! This is a custom command."
    echo "Current GPUMDkit path: ${GPUMDkit_path}"
}
```

and run `gpumdkit.sh -hello`, you will see:

```
Hello, GPUMDkit user! This is a custom command.
Current GPUMDkit path: /home/yanzihan/software/GPUMDkit
```

### Command with Required Arguments

For the following function:

```
custom_greet() {
    name=$1
    echo "Hello, ${name}!"
}
```

and you can add a argument, for example: `gpumdkit.sh -greet Zihan`, then you will see:

```
Hello, Zihan!
```

**Passing arguments**: Use `"$@"` to forward all arguments safely. Use `${@:1}` when you need to refer to arguments starting from the first one. 

### Calling an External Script

If you want to call an external script like `gpumdkit.sh`.

```
custom_nepanalyse() {
    python ~/my_gpumd_tools/analyse_nep.py
}
```

and run `gpumdkit.sh -nepanalyse`



## Advanced Examples

### Comprehensive Data Quality Pipeline

Create a complete quality control workflow:

```bash
custom_qc_pipeline() {
    local xyz_file=$1
    
    if [ -z "$xyz_file" ]; then
        echo "Usage: gpumdkit.sh -qc_pipeline <train.xyz>"
        return 1
    fi
    
    echo "╔═══════════════════════════════════════╗"
    echo "║  GPUMDkit Quality Control Pipeline    ║"
    echo "╚═══════════════════════════════════════╝"
    
    # 1. Property ranges
    echo ""
    echo "=== Step 1: Analyzing Property Ranges ==="
    gpumdkit.sh -range "$xyz_file" energy
    gpumdkit.sh -range "$xyz_file" force
    gpumdkit.sh -range "$xyz_file" virial
    
    # 2. Distance check
    echo ""
    echo "=== Step 2: Checking Minimum Distances ==="
    gpumdkit.sh -min_dist_pbc "$xyz_file"
    
    # 3. Composition analysis
    echo ""
    echo "=== Step 3: Analyzing Composition ==="
    gpumdkit.sh -analyze_comp "$xyz_file"
    
    # 4. Summary report
    echo ""
    echo "╔═══════════════════════════════════════╗"
    echo "║  Quality Control Complete!            ║"
    echo "╠═══════════════════════════════════════╣"
    echo "║  Review the results above.            ║"
    echo "║  Use filtering commands if needed:    ║"
    echo "║  - gpumdkit.sh -filter_value          ║"
    echo "║  - gpumdkit.sh -filter_dist           ║"
    echo "║  - gpumdkit.sh -filter_box            ║"
    echo "╚═══════════════════════════════════════╝"
}
```

**Usage:**
```bash
gpumdkit.sh -qc_pipeline train.xyz
```

---

### Automated Training Set Preparation

Complete workflow from VASP to NEP-ready data:

```bash
custom_prep_nep_training() {
    local dft_dir=${1:-.}
    local elements=${@:2}
    
    if [ -z "$elements" ]; then
        echo "Usage: gpumdkit.sh -prep_nep_training <dft_dir> <element1> <element2> ..."
        echo "Example: gpumdkit.sh -prep_nep_training ./vasp_calcs Li Y Cl"
        return 1
    fi
    
    echo "╔═══════════════════════════════════════╗"
    echo "║  NEP Training Data Preparation        ║"
    echo "╚═══════════════════════════════════════╝"
    
    # Step 1: Convert VASP outputs
    echo ""
    echo "=== Step 1/6: Converting VASP outputs ==="
    gpumdkit.sh -out2xyz "$dft_dir"
    
    # Step 2: Initial quality checks
    echo ""
    echo "=== Step 2/6: Initial quality checks ==="
    gpumdkit.sh -range train.xyz force
    gpumdkit.sh -min_dist_pbc train.xyz
    
    # Step 3: Filter high forces
    echo ""
    echo "=== Step 3/6: Filtering high forces (>30 eV/Å) ==="
    gpumdkit.sh -filter_value train.xyz force 30
    mv filtered_force.xyz train_filtered.xyz
    
    # Step 4: Filter close distances
    echo ""
    echo "=== Step 4/6: Filtering close distances (<1.5 Å) ==="
    gpumdkit.sh -filter_dist train_filtered.xyz 1.5
    mv filtered_dist.xyz train_clean.xyz
    
    # Step 5: Add group labels
    echo ""
    echo "=== Step 5/6: Adding group labels ==="
    if [ -f "POSCAR" ]; then
        gpumdkit.sh -addgroup POSCAR $elements
        echo "✓ model.xyz created with group labels"
    else
        echo "⚠ Warning: No POSCAR found, skipping group labels"
        echo "  You can add them later with: gpumdkit.sh -addgroup POSCAR $elements"
    fi
    
    # Step 6: Final verification
    echo ""
    echo "=== Step 6/6: Final verification ==="
    n_original=$(grep -c "^[0-9]" train.xyz 2>/dev/null || echo "0")
    n_filtered=$(grep -c "^[0-9]" train_clean.xyz 2>/dev/null || echo "0")
    n_removed=$((n_original - n_filtered))
    
    mv train_clean.xyz train.xyz
    
    echo ""
    echo "╔═══════════════════════════════════════╗"
    echo "║  Preparation Complete!                ║"
    echo "╠═══════════════════════════════════════╣"
    echo "║  Original structures: $n_original"
    echo "║  Filtered structures: $n_filtered"
    echo "║  Removed: $n_removed"
    echo "╠═══════════════════════════════════════╣"
    echo "║  Files ready for NEP training:        ║"
    echo "║  - train.xyz (training data)          ║"
    echo "║  - model.xyz (with group labels)      ║"
    echo "╠═══════════════════════════════════════╣"
    echo "║  Next steps:                          ║"
    echo "║  1. Split train/test if needed        ║"
    echo "║  2. Create nep.in                     ║"
    echo "║  3. Run: nep                          ║"
    echo "╚═══════════════════════════════════════╝"
}
```

**Usage:**
```bash
gpumdkit.sh -prep_nep_training ./vasp_calculations Li Y Cl
```

---

### Batch Visualization Suite

Generate all relevant plots at once:

```bash
custom_plot_all() {
    local save_flag=""
    if [ "$1" == "save" ]; then
        save_flag="save"
    fi
    
    echo "╔═══════════════════════════════════════╗"
    echo "║  Generating All Plots...              ║"
    echo "╚═══════════════════════════════════════╝"
    
    # NEP training plots
    if [ -f "loss.out" ]; then
        echo ""
        echo "=== NEP Training ==="
        gpumdkit.sh -plt train $save_flag
        gpumdkit.sh -plt force_errors $save_flag
        [ -f "energy_test.out" ] && gpumdkit.sh -plt prediction $save_flag
        [ -f "loss.out" ] && grep -q "learning_rate" loss.out && gpumdkit.sh -plt lr $save_flag
    fi
    
    # MD simulation plots
    if [ -f "thermo.out" ]; then
        echo ""
        echo "=== MD Simulation ==="
        gpumdkit.sh -plt thermo $save_flag
    fi
    
    if [ -f "msd.out" ]; then
        echo ""
        echo "=== Transport Properties ==="
        gpumdkit.sh -plt msd $save_flag
        gpumdkit.sh -plt sdc $save_flag
    fi
    
    if [ -f "rdf.out" ]; then
        echo ""
        echo "=== Structure ==="
        gpumdkit.sh -plt rdf $save_flag
    fi
    
    echo ""
    echo "╔═══════════════════════════════════════╗"
    echo "║  All Plots Generated!                 ║"
    echo "╚═══════════════════════════════════════╝"
}
```

**Usage:**
```bash
gpumdkit.sh -plot_all          # Display plots
gpumdkit.sh -plot_all save     # Save as PNG files
```

---

### Multi-Temperature Analysis

Analyze simulations across multiple temperatures:

```bash
custom_temp_scan() {
    echo "╔═══════════════════════════════════════╗"
    echo "║  Multi-Temperature Analysis           ║"
    echo "╚═══════════════════════════════════════╝"
    
    # Find all temperature directories
    temp_dirs=$(ls -d *K 2>/dev/null | sort -V)
    
    if [ -z "$temp_dirs" ]; then
        echo "Error: No temperature directories found (e.g., 600K, 900K)"
        return 1
    fi
    
    echo ""
    echo "Found temperature directories:"
    echo "$temp_dirs"
    echo ""
    
    # Create summary file
    echo "# Temperature Analysis Summary" > temp_summary.txt
    echo "# Generated: $(date)" >> temp_summary.txt
    echo "Temperature(K) | E_avg(eV) | T_avg(K) | P_avg(GPa) | D_total(cm²/s)" >> temp_summary.txt
    echo "---------------|-----------|----------|------------|----------------" >> temp_summary.txt
    
    for dir in $temp_dirs; do
        temp=$(echo $dir | sed 's/K//')
        echo ""
        echo "=== Processing $dir ==="
        
        cd $dir
        
        # Extract average values
        if [ -f "thermo.out" ]; then
            e_avg=$(awk 'NR>1 {sum+=$2; n++} END {if(n>0) printf "%.3f", sum/n}' thermo.out)
            t_avg=$(awk 'NR>1 {sum+=$3; n++} END {if(n>0) printf "%.1f", sum/n}' thermo.out)
            p_avg=$(awk 'NR>1 {sum+=$4; n++} END {if(n>0) printf "%.2f", sum/n}' thermo.out)
        fi
        
        # Calculate diffusivity if MSD available
        if [ -f "msd.out" ]; then
            # Simple linear fit to last 50% of MSD
            d_total=$(awk 'NR>1 {
                if(NR==2) start=NR;
                if(NR>1+int((NR-1)/2)) {
                    n++; sum_t+=$1; sum_msd+=$5;
                    sum_t2+=$1*$1; sum_t_msd+=$1*$5;
                }
            } END {
                if(n>0) {
                    slope=(n*sum_t_msd-sum_t*sum_msd)/(n*sum_t2-sum_t*sum_t);
                    D=slope/6*1e4;  # Convert to cm²/s
                    printf "%.2e", D;
                }
            }' msd.out)
        fi
        
        # Add to summary
        echo "$temp | ${e_avg:--} | ${t_avg:--} | ${p_avg:--} | ${d_total:--}" >> ../temp_summary.txt
        
        cd ..
    done
    
    echo ""
    echo "╔═══════════════════════════════════════╗"
    echo "║  Analysis Complete!                   ║"
    echo "╠═══════════════════════════════════════╣"
    echo "║  Summary saved to: temp_summary.txt   ║"
    echo "╚═══════════════════════════════════════╝"
    echo ""
    cat temp_summary.txt
}
```

**Usage:**
```bash
# In directory with subdirectories: 600K/ 900K/ 1200K/ etc.
gpumdkit.sh -temp_scan
```

---

### Automated Convergence Checker

Monitor NEP training convergence:

```bash
custom_check_convergence() {
    local max_iterations=${1:-10}
    local check_interval=${2:-60}
    
    echo "╔═══════════════════════════════════════╗"
    echo "║  NEP Training Convergence Monitor     ║"
    echo "╚═══════════════════════════════════════╝"
    echo ""
    echo "Checking every $check_interval seconds (max $max_iterations checks)"
    echo ""
    
    prev_energy_rmse=""
    prev_force_rmse=""
    iter=0
    
    while [ $iter -lt $max_iterations ]; do
        ((iter++))
        
        if [ ! -f "loss.out" ]; then
            echo "[$iter] Waiting for loss.out..."
            sleep $check_interval
            continue
        fi
        
        # Extract current RMSE
        energy_rmse=$(grep "Energy" loss.out | tail -1 | awk '{print $4}')
        force_rmse=$(grep "Force" loss.out | tail -1 | awk '{print $4}')
        generation=$(grep -c "Generation" loss.out)
        
        # Display current status
        echo "[$iter] Gen: $generation | E: $energy_rmse meV/atom | F: $force_rmse meV/Å"
        
        # Check convergence
        if [ -n "$prev_energy_rmse" ]; then
            e_change=$(echo "scale=4; ($prev_energy_rmse - $energy_rmse) / $prev_energy_rmse * 100" | bc)
            f_change=$(echo "scale=4; ($prev_force_rmse - $force_rmse) / $prev_force_rmse * 100" | bc)
            
            echo "      Change: E: ${e_change}% | F: ${f_change}%"
            
            # Check if converged (< 1% change)
            e_small=$(echo "$e_change < 1.0 && $e_change > -1.0" | bc)
            f_small=$(echo "$f_change < 1.0 && $f_change > -1.0" | bc)
            
            if [ "$e_small" -eq 1 ] && [ "$f_small" -eq 1 ]; then
                echo ""
                echo "╔═══════════════════════════════════════╗"
                echo "║  ✓ Training Converged!                ║"
                echo "╠═══════════════════════════════════════╣"
                echo "║  Energy RMSE: $energy_rmse meV/atom"
                echo "║  Force RMSE: $force_rmse meV/Å"
                echo "╚═══════════════════════════════════════╝"
                return 0
            fi
        fi
        
        prev_energy_rmse=$energy_rmse
        prev_force_rmse=$force_rmse
        
        sleep $check_interval
    done
    
    echo ""
    echo "╔═══════════════════════════════════════╗"
    echo "║  Monitoring timeout reached           ║"
    echo "║  Training still in progress           ║"
    echo "╚═══════════════════════════════════════╝"
}
```

**Usage:**
```bash
# Check every 60 seconds, up to 10 times
gpumdkit.sh -check_convergence

# Custom: check every 120 seconds, up to 20 times
gpumdkit.sh -check_convergence 20 120
```

---

## Best Practices and Patterns

### 1. Always Add Help Messages

```bash
custom_my_function() {
    if [ $# -eq 0 ]; then
        echo "Usage: gpumdkit.sh -my_function <arg1> <arg2>"
        echo ""
        echo "Description:"
        echo "  Does something useful with args"
        echo ""
        echo "Examples:"
        echo "  gpumdkit.sh -my_function file.xyz 100"
        return 1
    fi
    
    # Function logic here
}
```

---

### 2. Validate Input Files

```bash
custom_process_data() {
    local input_file=$1
    
    # Check if file exists
    if [ ! -f "$input_file" ]; then
        echo "Error: File '$input_file' not found!"
        return 1
    fi
    
    # Check if file is not empty
    if [ ! -s "$input_file" ]; then
        echo "Error: File '$input_file' is empty!"
        return 1
    fi
    
    # Check file format (extxyz should have numbers on first lines)
    if ! head -1 "$input_file" | grep -q "^[0-9]"; then
        echo "Error: '$input_file' doesn't look like extxyz format!"
        return 1
    fi
    
    # Process file...
}
```

---

### 3. Use Meaningful Variable Names

```bash
# Bad
custom_func() {
    a=$1
    b=$2
    c=$(cat $a | grep $b)
}

# Good
custom_analyze_structure() {
    local xyz_file=$1
    local element=$2
    local element_count=$(grep "$element" "$xyz_file" | wc -l)
    echo "Found $element_count $element atoms in $xyz_file"
}
```

---

### 4. Provide Progress Feedback

```bash
custom_long_process() {
    local files=("$@")
    local total=${#files[@]}
    local current=0
    
    echo "Processing $total files..."
    
    for file in "${files[@]}"; do
        ((current++))
        echo "[$current/$total] Processing $file..."
        
        # Do processing
        gpumdkit.sh -some_operation "$file"
        
        echo "  ✓ Complete"
    done
    
    echo "All files processed!"
}
```

---

### 5. Handle Errors Gracefully

```bash
custom_safe_conversion() {
    local input=$1
    local output=$2
    
    echo "Converting $input..."
    
    if gpumdkit.sh -pos2exyz "$input" "$output" 2>/dev/null; then
        echo "✓ Success: $output created"
        return 0
    else
        echo "✗ Error: Conversion failed"
        echo "  Check that $input is a valid POSCAR file"
        return 1
    fi
}
```

---

### 6. Create Backup Before Destructive Operations

```bash
custom_filter_and_backup() {
    local xyz_file=$1
    local threshold=$2
    
    # Create backup
    local backup="${xyz_file}.backup_$(date +%Y%m%d_%H%M%S)"
    echo "Creating backup: $backup"
    cp "$xyz_file" "$backup"
    
    # Perform filtering
    echo "Filtering..."
    gpumdkit.sh -filter_value "$xyz_file" force "$threshold"
    
    echo "Original saved to: $backup"
    echo "Filtered data in: filtered_force.xyz"
}
```

---

### 7. Log Operations

```bash
custom_logged_workflow() {
    local logfile="workflow_$(date +%Y%m%d_%H%M%S).log"
    
    echo "Starting workflow, logging to $logfile"
    
    {
        echo "=== Workflow Started: $(date) ==="
        echo "Arguments: $@"
        echo ""
        
        # Run operations
        echo "Step 1: Converting files..."
        gpumdkit.sh -out2xyz .
        
        echo "Step 2: Filtering..."
        gpumdkit.sh -filter_value train.xyz force 30
        
        echo ""
        echo "=== Workflow Completed: $(date) ==="
    } | tee "$logfile"
    
    echo "Log saved to: $logfile"
}
```

---

## Common Patterns Library

### Pattern 1: Iterate Over Directories

```bash
custom_process_all_dirs() {
    for dir in */; do
        echo "=== Processing $dir ==="
        cd "$dir"
        
        # Do something
        gpumdkit.sh -plt thermo save
        
        cd ..
    done
}
```

---

### Pattern 2: Conditional Execution

```bash
custom_smart_plot() {
    # Plot training if files exist
    [ -f "loss.out" ] && gpumdkit.sh -plt train save
    
    # Plot MD if files exist
    [ -f "thermo.out" ] && gpumdkit.sh -plt thermo save
    [ -f "msd.out" ] && gpumdkit.sh -plt msd save
    
    # Plot RDF if exists
    [ -f "rdf.out" ] && gpumdkit.sh -plt rdf save
    
    echo "Generated all available plots"
}
```

---

### Pattern 3: Parallel Processing

```bash
custom_parallel_convert() {
    local max_jobs=4
    local running=0
    
    for dir in struct_*/; do
        # Convert in background
        (
            cd "$dir"
            gpumdkit.sh -pos2exyz POSCAR structure.xyz
        ) &
        
        ((running++))
        
        # Wait if too many jobs running
        if [ $running -ge $max_jobs ]; then
            wait -n  # Wait for any job to finish
            ((running--))
        fi
    done
    
    wait  # Wait for all remaining jobs
    echo "All conversions complete!"
}
```

---

### Pattern 4: Data Aggregation

```bash
custom_collect_rmse() {
    echo "Collecting RMSE from all iterations..."
    echo "Iteration | Energy RMSE | Force RMSE" > rmse_summary.txt
    
    for dir in iteration_*/; do
        iter=$(echo $dir | grep -o "[0-9]*")
        cd "$dir"
        
        if [ -f "loss.out" ]; then
            e_rmse=$(grep "Energy" loss.out | tail -1 | awk '{print $4}')
            f_rmse=$(grep "Force" loss.out | tail -1 | awk '{print $4}')
            echo "$iter | $e_rmse | $f_rmse" >> ../rmse_summary.txt
        fi
        
        cd ..
    done
    
    cat rmse_summary.txt
}
```

---

### Pattern 5: Interactive Selection

```bash
custom_interactive_filter() {
    local xyz_file=$1
    
    # Show current statistics
    echo "Current statistics for $xyz_file:"
    gpumdkit.sh -range "$xyz_file" force
    
    # Ask user for threshold
    read -p "Enter force threshold (eV/Å): " threshold
    
    if [ -n "$threshold" ]; then
        echo "Filtering with threshold: $threshold eV/Å"
        gpumdkit.sh -filter_value "$xyz_file" force "$threshold"
        echo "Filtered data saved to: filtered_force.xyz"
    else
        echo "No threshold provided, skipping filter"
    fi
}
```

---

## Tips for Writing Custom Commands

### 💡 Tip 1: Test Incrementally

Build complex commands step by step:

```bash
# Start simple
custom_test() {
    echo "Hello, GPUMDkit!"
}

# Add argument handling
custom_test() {
    echo "Hello, $1!"
}

# Add validation
custom_test() {
    if [ -z "$1" ]; then
        echo "Usage: gpumdkit.sh -test <name>"
        return 1
    fi
    echo "Hello, $1!"
}

# Add logic
custom_test() {
    if [ -z "$1" ]; then
        echo "Usage: gpumdkit.sh -test <name>"
        return 1
    fi
    echo "Hello, $1!"
    echo "Time: $(date)"
    echo "Directory: $(pwd)"
}
```

---

### 💡 Tip 2: Use Functions for Reusability

```bash
# Helper function (not a custom command)
print_header() {
    local title=$1
    echo "╔═══════════════════════════════════════╗"
    echo "║  $title"
    echo "╚═══════════════════════════════════════╝"
}

# Use in multiple custom commands
custom_analyze() {
    print_header "Structure Analysis"
    # analysis code...
}

custom_process() {
    print_header "Data Processing"
    # processing code...
}
```

---

### 💡 Tip 3: Document in Code

```bash
custom_complex_workflow() {
    # Purpose: Automate NEP training data preparation
    # Author: Your Name
    # Date: 2025-01-15
    # Dependencies: gpumdkit.sh, VASP tools
    
    # Step 1: Convert VASP outputs to extxyz
    # This processes all OUTCAR files in subdirectories
    gpumdkit.sh -out2xyz ./vasp_calculations
    
    # Step 2: Apply quality filters
    # Remove structures with forces > 30 eV/Å
    gpumdkit.sh -filter_value train.xyz force 30
    
    # Continue workflow...
}
```

---

### 💡 Tip 4: Make Commands Composable

Design commands that work well together:

```bash
# Each command does one thing well
custom_convert_all() {
    for d in struct_*/; do
        gpumdkit.sh -out2xyz "$d"
    done
}

custom_merge_all() {
    cat struct_*/train.xyz > combined.xyz
}

custom_filter_merged() {
    gpumdkit.sh -filter_value combined.xyz force "$1"
}

# Use together:
# gpumdkit.sh -convert_all
# gpumdkit.sh -merge_all
# gpumdkit.sh -filter_merged 30
```

---

### 💡 Tip 5: Version Your Custom Commands

```bash
# In ~/.gpumdkit.in, add a header:
# GPUMDkit Custom Commands
# Version: 2.0
# Last updated: 2025-01-15
# Author: Your Name

custom_version() {
    echo "Custom commands version 2.0"
    echo "Available commands:"
    grep "^custom_" ~/.gpumdkit.in | sed 's/custom_/  - /' | sed 's/().*$//'
}
```

---

## Troubleshooting Custom Commands

### Issue: Command Not Recognized

**Symptoms:**
```bash
gpumdkit.sh -my_function
# Error: Unknown option: -my_function
```

**Solutions:**

1. **Check function name starts with `custom_`**
   ```bash
   # Wrong
   my_function() { ... }
   
   # Correct
   custom_my_function() { ... }
   ```

2. **Reload shell configuration**
   ```bash
   source ~/.bashrc
   # Or restart terminal
   ```

3. **Verify file location**
   ```bash
   ls -la ~/.gpumdkit.in
   # Should exist and be readable
   ```

---

### Issue: Function Doesn't Execute

**Symptoms:**
```bash
gpumdkit.sh -my_function
# Nothing happens or immediate return
```

**Solutions:**

1. **Check for syntax errors**
   ```bash
   bash -n ~/.gpumdkit.in
   # Reports syntax errors without executing
   ```

2. **Add debug output**
   ```bash
   custom_my_function() {
       echo "Function started"
       # Your code
       echo "Function completed"
   }
   ```

3. **Check return codes**
   ```bash
   custom_my_function() {
       gpumdkit.sh -some_operation || {
           echo "Operation failed!"
           return 1
       }
   }
   ```

---

### Issue: Arguments Not Passed Correctly

**Symptoms:**
```bash
gpumdkit.sh -my_function arg1 arg2
# Only receives arg1
```

**Solutions:**

```bash
# Use "$@" for all arguments
custom_my_function() {
    echo "All arguments: $@"
    echo "Number of arguments: $#"
    echo "First argument: $1"
    echo "Second argument: $2"
}

# For operations on multiple files
custom_process_files() {
    for file in "$@"; do
        echo "Processing: $file"
    done
}
```

---

## Sharing Custom Commands

### Option 1: Share File

```bash
# Export your commands
cp ~/.gpumdkit.in my_gpumdkit_commands.sh

# Others can use:
cat my_gpumdkit_commands.sh >> ~/.gpumdkit.in
source ~/.bashrc
```

---

### Option 2: Create Library

```bash
# Create a shared library
cat > gpumdkit_library.sh << 'EOF'
# GPUMDkit Command Library
# Source this file in ~/.gpumdkit.in

# Quality control suite
custom_qc() { ... }

# Batch processing
custom_batch() { ... }

# Add more commands...
EOF

# In ~/.gpumdkit.in:
source /path/to/gpumdkit_library.sh
```

---

### Option 3: Documentation

Create a README for your custom commands:

```markdown
# My GPUMDkit Custom Commands

## Available Commands

### Quality Control
- `gpumdkit.sh -qc_pipeline <file>` - Run complete QC
- `gpumdkit.sh -qc_quick <file>` - Quick checks only

### Analysis
- `gpumdkit.sh -analyze_all` - Analyze all MD outputs
- `gpumdkit.sh -temp_scan` - Multi-temperature analysis

### Utilities
- `gpumdkit.sh -clean_backups` - Remove old backup files
- `gpumdkit.sh -archive <dir>` - Archive completed work

## Installation

```bash
cat my_commands.sh >> ~/.gpumdkit.in
source ~/.bashrc
```

## Examples

...
```

---

---

Thank you for using `GPUMDkit`! If you have any questions or need further assistance, feel free to open an issue on our GitHub repository or contact Zihan YAN (yanzihan@westlake.edu.cn).