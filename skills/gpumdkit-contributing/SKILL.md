---
name: gpumdkit-contributing
description: >
  Use when adding new features, commands, or scripts to GPUMDkit.
  Provides coding conventions, file structure patterns, and step-by-step guides
  for implementing new interactive menu items, CLI flags, calculators, analyzers,
  format converters, plotting scripts, and workflow tools.
  Use when user asks about: adding a new command, creating a new script,
  contributing to GPUMDkit, code conventions, or how to extend the toolkit.
allowed-tools: Bash(gpumdkit *) Bash(python3 *) Bash(python *) Bash(ls *) Bash(mkdir *)
---

# GPUMDkit Contributing Guide

## Project Structure

```
GPUMDkit/
├── gpumdkit.sh              # Main entry: CLI routing + interactive menu
├── install.sh               # Installation script
├── src/                     # Shell menu modules (interactive mode wrappers)
│   ├── f1_format_conversions.sh
│   ├── f2_sample_structures.sh
│   ├── f3_workflows.sh
│   ├── f4_calculators.sh
│   ├── f5_analyzers.sh
│   ├── f6_plots.sh
│   └── f7_utilities.sh
├── Scripts/                 # Implementation scripts (Python + Bash)
│   ├── format_conversion/
│   ├── sample_structures/
│   ├── workflow/
│   ├── calculators/
│   ├── analyzer/
│   ├── plt_scripts/
│   └── utils/
├── skills/                  # AI agent skill definitions
├── docs/                    # Documentation
│   ├── tutorials/en/        # English tutorials
│   ├── tutorials/zh/        # Chinese tutorials
│   ├── mkdocs.yml           # MkDocs config
│   └── htmls/               # Generated HTML
```

## How Routing Works

`gpumdkit.sh` has two modes:

1. **Interactive mode**: `gpumdkit.sh` (no arguments) → shows numbered menu
   - User selects a category (1-7) → sources the corresponding `src/fN_*.sh`
   - User selects a function (e.g., 403) → calls the function (e.g., `f403_calc_descriptors`)
   - The function prompts for input, calls a Python/Bash script from `Scripts/`, prints the code path

2. **CLI mode**: `gpumdkit.sh -flag [args]` → routes via `case $1 in`
   - Each flag maps to a script in `Scripts/`
   - Help messages go in `help_info_table()` and `calculator_help_table()`

## Adding a New Interactive Mode Feature

### Step 1: Create the implementation script

Place it in the appropriate `Scripts/` subdirectory:

```bash
# Python script example: Scripts/calculators/calc_new_feature.py
```

### Step 2: Add a wrapper function in `src/fN_*.sh`

Follow this exact pattern:

```bash
function f413_new_feature(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/calculators       |"
echo " | Script: calc_new_feature.py                     |"
echo " | Developer: Your Name (your@email.com)           |"
echo " >-------------------------------------------------<"
echo " Input <param1> <param2> [optional_param3]"
echo " Example: input.xyz nep.txt"
echo " ------------>>"
read -r -a feature_args
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/calc_new_feature.py "${feature_args[@]}"
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_new_feature.py"
echo " ---------------------------------------------------"
}
```

Key conventions:
- Banner uses `>-----<` and `| ... |` format
- `echo " Input <...>"` describes parameters
- `echo " Example: ..."` shows a concrete example (not "Examp:")
- `read -r -a varname` reads into an array
- `"${varname[@]}"` passes all arguments properly quoted
- Function name format: `f<Category><Number>_<descriptive_name>`

### Step 3: Register in `gpumdkit.sh`

1. Add the choice number to `array_choice` (line ~57)
2. Add the case in the `main()` function's nested case statement:
   ```bash
   "413") f413_new_feature ;;
   ```

### Step 4: Update the menu display in `src/f6_plots.sh` or the relevant `fN_*.sh`

If it's a calculator, update the calculator menu display in `f4_calculators.sh`.

## Adding a New CLI Flag

### Step 1: Create the implementation script

Same as above — place in `Scripts/`.

### Step 2: Add the CLI handler in `gpumdkit.sh`

Add a new case in the main `case $1 in` block:

```bash
-new_flag)
    if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ]; then
        python ${analyzer_path}/new_script.py $2 $3
    else
        echo " Usage: -new_flag <param1> <param2>"
        echo " Example: gpumdkit.sh -new_flag input.xyz 3.0"
        echo " Code path: ${analyzer_path}/new_script.py"
    fi ;;
```

Key conventions:
- Check `$2 != "-h"` to handle help requests
- No quotes around `$2`, `$3` in the python call (paths with spaces not supported)
- No `exit 1` in the else branch (just print usage and fall through)
- Place alphabetically among existing flags for readability

### Step 3: Add to help table

Update `help_info_table()` or `calculator_help_table()` in `gpumdkit.sh`:

```bash
echo "| -new_flag  Description of new flag           | -existing_flag    Description                   |"
```

### Step 4: Update tab completion

Edit `Scripts/utils/completion.sh`:

1. Add `-new_flag` to the `opts` string (line 16)
2. If it needs file completion, add it to the `|`-separated case (line 33):
   ```bash
   -out2xyz|-...|-new_flag)
       COMPREPLY=($(compgen -f -- "$cur")) ;;
   ```
3. If it needs secondary options, add a new case block

## Python Script Conventions

### File header

```python
"""
Script:     new_script.py
Category:   Calculator Scripts
Purpose:    Brief description of what this script does.
Usage:      python new_script.py <param1> <param2>
Arguments:
  param1    Description of param1
  param2    Description of param2
Author:     Your Name (your@email.com)
"""
```

### Argument parsing

Use `sys.argv` for simple scripts (most common pattern):

```python
import sys

if len(sys.argv) != 3:
    print("Usage: python new_script.py <param1> <param2>")
    print("Example: python new_script.py input.xyz nep.txt")
    sys.exit(1)

param1 = sys.argv[1]
param2 = sys.argv[2]
```

Use `argparse` only if the script has many optional flags with complex defaults.

### Main guard (optional)

GPUMDkit scripts run directly — a `main()` function with `if __name__ == "__main__":` guard is optional. Most existing scripts execute at module level. Use whichever style you prefer, but do not wrap the entire script in a `main()` function just for the sake of it.

### Dependency notices

If your script requires a heavy/special package (`NepTrain`, `calorine`, `dpdata`), add a notice:

```python
def print_dependency_notice():
    print("This function requires the 'calorine' package.")
    print("If you use this function, please cite:")
    print("  Bochkov et al., npj Comput Mater 6, 170 (2020)")

# Call it early in the script
print_dependency_notice()
```

Do NOT add notices for common packages like `ase`, `numpy`, `pymatgen`.

### Error messages

```python
# Good
print("Error: File 'input.xyz' not found.")
sys.exit(1)

# Avoid
raise FileNotFoundError("input.xyz")
```

## Shell Script Conventions

### File header (in `src/`)

```bash
# ============================================================
# GPUMDkit <module name> module
# Repository: https://github.com/zhyan0603/GPUMDkit
# Author: Your Name (your@email.com)
# ============================================================
```

### Variable checks

```bash
# Non-empty check (project convention)
if [ ! -z "$var" ]; then
    ...
fi

# Not the -h flag check
if [ ! -z "$2" ] && [ "$2" != "-h" ]; then
    ...
fi
```

### Interactive menu pattern

```bash
function f413_new_feature(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/calculators       |"
echo " | Script: calc_new_feature.py                     |"
echo " | Developer: Your Name (your@email.com)           |"
echo " >-------------------------------------------------<"
echo " Input <param1> <param2>"
echo " Example: input.xyz nep.txt"
echo " ------------>>"
read -r -a feature_args
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/calc_new_feature.py "${feature_args[@]}"
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_new_feature.py"
echo " ---------------------------------------------------"
}
```

## Validation Checklist

Before committing, run these checks:

```bash
# Shell syntax
bash -n gpumdkit.sh
find src Scripts -name '*.sh' -exec bash -n {} +

# Python syntax
python3 -m py_compile Scripts/path/to/your_script.py

# Clean up __pycache__ after py_compile
find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null

# MkDocs build (if docs changed)
mkdocs build -f docs/mkdocs.yml

# Whitespace check
git diff --check
```
