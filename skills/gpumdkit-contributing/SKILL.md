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
read_menu_choice feature_args || return 1
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/calc_new_feature.py "${feature_args[@]}"
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_new_feature.py"
echo " ---------------------------------------------------"
}
```

Key conventions:
- Banner uses `>-----<` and `| ... |` format (49 dashes inner, 53 chars total)
- `echo " Input <...>"` describes parameters
- `echo " Example: ..."` shows a concrete example
- Use `read_menu_choice var || return 1` for single-value input, `read_menu_array arr || return 1` for array input
- **Never use bare `read -p` or `read -r -a`** — these hang on closed stdin. The helpers are defined in `gpumdkit.sh`
- `"${varname[@]}"` passes all arguments properly quoted
- Function name format: `f<Category><Number>_<descriptive_name>`
- **Script banner must name the EXACT filename** (not an abbreviation or wrong suffix)
- All banners use uniform width: 49 dashes inner (lines 82-86 style). Inconsistent widths (51/52/53/56) should be normalized when touching a file.

### Step 3: Register in `gpumdkit.sh`

1. Add the choice number to `array_choice` (line ~64)
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
2. If it needs file completion, add it to the `|`-separated case (line 38):
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
# Good (leading space, clear message, sys.exit)
print(" Error: File 'input.xyz' not found.")
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
read_menu_choice feature_args || return 1   # EOF-safe; see project conventions
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/calc_new_feature.py "${feature_args[@]}"
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_new_feature.py"
echo " ---------------------------------------------------"
}
```

## Project Conventions & Personal Preferences

These are the maintainer's conventions. Follow them strictly.

### Echo / Print Leading Space

**Every user-facing `echo`/`print` must start with a space.**

```bash
# ✓ Correct
echo " Error: File not found."
echo " Deleting the files..."
echo " Operation canceled."

# ✗ Wrong
echo "Error: File not found."
```

```python
# ✓ Correct
print(" Error: 'input.xyz' not found.")
print(" Deleting the files...")

# ✗ Wrong
print("Error: 'input.xyz' not found.")
```

This applies to ALL script types: interactive wrappers, workflow scripts, utility messages, tutorial transcripts, and Python print() output. The only exception is the `" ------------>>"` prompt marker which starts with a space by design.

### ASCII Only Box Drawing

All border characters must be pure ASCII: `|`, `+`, `-`, `/`, `\`, `>`, `<`. **Never use** Unicode box-drawing characters like `│` (U+2502) — they render as mojibake on non-UTF-8 terminals.

### EOF-Safe Input

All interactive input must go through the helpers defined in `gpumdkit.sh`:

- `read_menu_choice var || return 1` — reads a single line into a variable
- `read_menu_array arr || return 1` — reads whitespace-separated tokens into an array

**Never use bare `read -p` or `read -r -a`.** When stdin is closed (pipe, redirect, /dev/null), bare read silently succeeds with empty/unset values or hangs forever. The helpers print `" Input closed. Exiting."` and return 1.

### Menu Validation

- Validation array: always named `valid_menu_choices`
- Validation loop: `while ! echo "${valid_menu_choices[@]}" | grep -wq "$num_choice"; do ...; done`
- Return-to-main-menu: `case ... "000") menu; main ;;`

### Script Banners

The info banner displayed before running a script follows this exact format:

```bash
echo " >-------------------------------------------------<"   # exactly 49 dashes
echo " | Calling the script in Scripts/<category>       |"
echo " | Script: <exact_filename>                         |"   # MUST match actual filename
echo " | Developer: Name (email)                          |"
echo " >-------------------------------------------------<"
```

- The script name in the banner **must be the exact filename** (e.g., `scf_batch_pretreatment_vasp.sh`, not `scf_batch_pretreatment.sh`)
- Width: 49 dashes inner (53 chars total between `>` and `<`). Keep uniform across projects.

### The Prompt Arrow

User input follows the marker pattern: `echo " ------------>>"` before reading input. This is consistent across all interactive functions.

### Error Message Format

Shell: `echo " Error: <message>."` — leading space, capital E, colon, period at end.
Python: `print(" Error: <message>."); sys.exit(1)` — same format, always call sys.exit after.

### What NOT to Propose or Do

These changes have been **explicitly rejected** by the project maintainer. Do not suggest them.

| Rejected | Reason |
|---|---|
| Add `if __name__ == "__main__":` to Python scripts | Scripts run standalone; main guard is unnecessary |
| Replace `from pylab import *` with explicit imports | Keep as-is |
| Unify `-plt` "save" argument position | Different scripts have different arg counts; keep flexible |
| Unify DPI (150 → 300) across plotting scripts | Keep per-script DPI settings |
| Change `exit` → `return` in workflow sourced scripts | Keep as-is |
| Add `-filter_value` / `-filter_range` / `-get_volume` / `-re_atoms` to the help table | Intentionally undocumented |
| Make interactive mode loop back to menu after one function | Keeps shots single-shot; designed this way |
| Add color/ANSI escape codes | User explicitly rejected; keep monochrome |
| Add debug examples to every tutorial page | Unnecessary; keep tutorials lean |
| Change MSD fitting range arbitrarily | Scientific choice; ask maintainer first |
| Add `des_compare` to CLI | Script exists but user chose not to expose it |

### macOS / Cross-Platform Notes

- macOS ships **bash 3.2** (no `local -n` nameref). Use `read -a "$varname"` with eval or herestring instead of namerefs.
- macOS **sed is BSD**, not GNU. Cross-platform sed: use `perl -i -pe` instead.
- **zsh does not support `read -a`** (bash-only). All scripts run under bash.
- `timeout` command is not available on macOS. Use background jobs or rely on EOF-safe read.

### Specific Design Decisions

- **calc_ion_conductivity.py**: MSD fitting range is 40%-80% of data (lines 107-108). The Arrhenius plotting scripts use the same range. Do not change without maintainer approval.
- **clean_extra_files.sh**: Keep-substring deletion uses **whole-word matching**, NOT substring removal. Do not regress to `${var/pattern/}` substring removal (it corrupts filenames sharing substrings).
- **troubleshooting.md**: Only 5 FAQ items. Do not expand without asking.
- **select_max_modev.py line 111**: The index mapping `atoms_list[filtered_indices[i]]` was recently fixed from a bug (`atoms_list[i]`). Do not undo.
- **charge_balance_check.py lines 118-131**: Now uses manual loop instead of `dict()` to handle 3-tuple error returns. Do not revert to `dict()`.

---

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
