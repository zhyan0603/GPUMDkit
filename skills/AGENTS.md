# GPUMDkit — AI Agent Knowledge File

> This document records architecture, conventions, preferences, and known pitfalls of the GPUMDkit project for future AI agents. Last updated: 2026-07-01 (project VERSION="1.5.6 (dev)").

## Project Overview

GPUMDkit is a Bash + Python command-line toolkit for GPUMD and NEP. It provides format conversion, structure sampling, workflows, property calculation, plotting, and analysis utilities. The entry point is `gpumdkit.sh` (dual-mode: interactive menu + CLI).

**Repository**: https://github.com/zhyan0603/GPUMDkit
**Primary Author**: Zihan YAN (yanzihan@westlake.edu.cn)
**Citation**: MGE Advances, 2026, 4, e70074

## Architecture

```
gpumdkit.sh          ← Main entry point: menu rendering + CLI dispatch
  │
  ├── src/f1_format_conversions.sh  ← Interactive menu wrappers (f101-f110 + keyword converters)
  ├── src/f2_sample_structures.sh   ← f201-f205
  ├── src/f3_workflows.sh           ← f301-f303
  ├── src/f4_calculators.sh         ← f401-f412
  ├── src/f5_analyzers.sh           ← f501-f508
  ├── src/f6_plots.sh               ← f6_plots_one_column / f6_plots_two_column (dynamic layout)
  │     └── Layout is selected based on terminal width
  ├── src/f7_utilities.sh           ← f701
  │
  └── Scripts/                      ← Implementation scripts
        ├── format_conversion/      ← Python/Bash converters
        ├── sample_structures/      ← Sampling utilities
        ├── calculators/            ← Calculators (MSD, ionic conductivity, NEB, descriptors, ...)
        ├── analyzer/               ← Analyzers (composition, distance, filter, ...)
        ├── plt_scripts/            ← Plotting (35 scripts)
        ├── workflow/               ← Workflows (SCF pretreatment, active learning)
        └── utils/                  ← Update, clean, completion, NEP modifier
```

## Key Files

| File | Purpose |
|------|---------|
| `gpumdkit.sh` | Main entry: defines `GPUMDkit_path`, `menu()`, `main()`, `help_info_table()`, `calculator_help_table()`, `citation()`, CLI `case`, `read_menu_choice()`, `read_menu_array()`, `run_python_script()` |
| `src/fN_*.sh` | Interactive menu modules. Each defines `fN_category()` display functions and leaf `fNXX_*` wrapper functions with `valid_menu_choices` validation loops. |
| `Scripts/utils/completion.sh` | Bash TAB completion |
| `Scripts/utils/update_gpumdkit.sh` | Update command |
| `Scripts/utils/clean_extra_files.sh` | Clean command |
| `install.sh` | Installation (writes `GPUMDkit_path` to bashrc/zshrc) |
| `docs/tutorials/en/*.md` + `zh/*.md` | Bilingual tutorials (13 pairs + index + troubleshooting) |
| `skills/*/SKILL.md` | AI agent skill definitions |
| `~/.gpumdkit.in` | User custom commands |

## Variable Reference

```bash
GPUMDkit_path          # Set by user; points to GPUMDkit root directory
VERSION                # "1.5.6 (dev) (2026-06-17)" (only in gpumdkit.sh)
plt_path               # ${GPUMDkit_path}/Scripts/plt_scripts
analyzer_path          # ${GPUMDkit_path}/Scripts/analyzer
calc_path              # ${GPUMDkit_path}/Scripts/calculators
workflow_path          # ${GPUMDkit_path}/Scripts/workflow
format_conv_path       # ${GPUMDkit_path}/Scripts/format_conversion
utils_path             # ${GPUMDkit_path}/Scripts/utils
sample_path            # ${GPUMDkit_path}/Scripts/sample_structures
```

## Shell Script Conventions

### Entry Point / Dispatch

**Interactive**: `gpumdkit.sh` (no args) → `menu()` → `main()` reads choice → first character of choice selects module → nested `case` calls `fNXX_*` functions → after function returns, `citation()` runs → script exits (no loop back to menu).

**CLI**: `gpumdkit.sh -flag [args]` → dispatched via `case $1 in`.

### Read Conventions (EOF Safety)

**Always use the EOF-safe helpers defined in `gpumdkit.sh`:**

- `read_menu_choice var || return 1` — reads a single line into a variable
- `read_menu_array arr || return 1` — reads whitespace-separated tokens into an array

**Never use bare `read -p " "` or `read -r -a`.** On closed stdin (pipe, redirect, /dev/null), bare read silently hangs or returns empty values without detection. The helpers print `" Input closed. Exiting."` and return 1.

### Menu Conventions

- Validation array: always named `valid_menu_choices`
- Validation loop: `while ! echo "${valid_menu_choices[@]}" | grep -wq "$num_choice"; do ...; done`
- Return to main menu: `case ... "000") menu; main ;;` (recursive)

### Box Drawing Styles

- **Menu boxes**: `+---+` style, constructed with sequential `echo` lines.
- **Script banners**: `>----<` style, inner width 49 dashes (53 chars total including `>` and `<`). Content uses ` | ... |` wrapper.
- **Input prompt marker**: `echo " ------------>>"` before any `read_menu_choice` or `read_menu_array`.

### Error Message Style

**Every user-facing `echo` must start with a space character:**

```bash
echo " Error: GPUMDkit_path is not set."        # ✓ Correct
echo " Error: Unable to determine current branch."  # ✓ Correct
echo "Error: ..."   # ✗ Wrong — missing leading space

# Normal informational messages also use leading space:
echo " Deleting the files..."   # ✓
echo " Operation canceled."     # ✓
```

### ASCII Only Box Drawing

All border characters must be pure ASCII: `|`, `+`, `-`, `/`, `\`, `>`, `<`. **Never use** Unicode box-drawing characters like `│` (U+2502) — they render as mojibake on non-UTF-8 terminals.

### Function Naming

- Interactive entry functions: `f<Module><Number>_<descriptive_name>` (e.g., `f101_out2xyz`, `f302_md_sample_batch_pretreatment_gpumd`)
- Menu display functions: `f<Module>_<module_name>` (e.g., `f1_format_conversion`, `f4_calculators`)
- Utility functions: lowercase snake_case (e.g., `clean_extra_files`, `update_gpumdkit`)
- `run_python_script` (helper in gpumdkit.sh): prints author + path banner, then calls `python "$@"`

### Other Shell Conventions

- Non-empty check: `[ ! -z "$var" ]` (project convention)
- CLI help flag check: `[ "$2" != "-h" ]`
- `"${@:2}"` / `"${@:3}"` patterns for argument forwarding
- `run_python_script "Author" /path/to/script.py "${@:2}"` used in CLI branches
- Functions are `source`d into the current shell; `type -t` can verify definition

## Python Script Conventions

- **Scripts are designed to run standalone, not to be imported.** Module-level `sys.exit(1)`, `input()`, `plt.show()` etc. are all acceptable.
- **Arguments**: prefer `sys.argv` positional (88% of scripts), `argparse` only for complex flag-heavy scripts (12%).
- **`np.loadtxt()` default `comments='#'`** handles GPUMD's newer output files that have `#` header lines automatically. No extra handling needed.
- **`np.savetxt(header=...)`** auto-prefixes the header string with `# `. Writers and readers are pair-compatible.
- **File header docstring**: all scripts over 50 lines should have a module-level docstring (Script / Author / Purpose / Usage).
- **Error exits**: `print(" Error: ..."); sys.exit(1)` — do not raise raw exceptions.
- `ase.io.read` provides reasonably clear errors when files are missing; this can serve as a minimal defense.

## Documentation Conventions

- **Tutorials**: `docs/tutorials/en/` + `zh/` bilingual. Any feature change should update the corresponding tutorial page in both languages.
- **SKILL.md**: `skills/<category>/SKILL.md` with YAML frontmatter (the `name` field must match the directory name exactly).
- **Help tables**: `help_info_table()` and `calculator_help_table()` in `gpumdkit.sh` — CLI flag tables must be kept in sync manually.

## Contributor Rules

### Must Do

1. All CLI help text and tutorial `echo`/`print` messages start with a **space** character.
2. All box-drawing characters must be **pure ASCII** (no `│` or other Unicode).
3. All interactive input reads must go through `read_menu_choice` or `read_menu_array` (**never** bare `read`).
4. Menu validation arrays must be named `valid_menu_choices`.
5. Script banners must name the **exact** filename (no abbreviations or wrong suffixes).
6. Changes must not break existing functionality — **minimum-invasion principle**.
7. After any change, validate with `bash -n` and `python3 -m py_compile`.

### Do NOT Do (explicitly rejected by maintainer)

| Rejected | Rationale |
|---|---|
| Force `if __name__ == "__main__":` on Python scripts | Scripts are standalone; main guards are unnecessary churn |
| Replace `from pylab import *` with explicit imports | Keep as-is |
| Unify `-plt` "save" argument position | Different scripts have different arg counts; inflexible |
| Unify DPI (e.g., 150 → 300) across plotting scripts | Keep per-script DPI |
| Change `exit` → `return` in workflow sourced scripts | Keep as-is |
| Make interactive mode loop back to menu after one function | Intentionally single-shot |
| Add `-filter_value` / `-filter_range` / `-get_volume` / `-re_atoms` to help table | Intentionally undocumented |
| Add color / ANSI escape codes | User rejected; keep monochrome |
| Add debug examples to every tutorial page | Keep tutorials lean |
| Change MSD fitting range without maintainer approval | 40%-80% is the settled range; changing it is a scientific decision |
| Add `des_compare` to CLI | Script exists but not exposed |
| Add extra troubleshooting FAQ items | Keep the 5-item set |

### macOS / Cross-Platform Notes

- macOS ships **bash 3.2** (no `local -n` nameref). Use `read -a "$varname"` for array reads.
- macOS **sed is BSD**, not GNU. For cross-platform sed, use `perl -i -pe` instead.
- **zsh does not support `read -a`** (bash-only). All scripts run under bash.
- `timeout` command is unavailable on macOS. Rely on EOF-safe `read_menu_choice` / `read_menu_array` to avoid hangs.

## Key Data Flows

### Interactive Mode

```
User runs: gpumdkit.sh
         → menu() prints menu
         → main() reads choice
         → ${choice:0:1} selects module (1-7)
         → source src/fN_*.sh
         → nested case calls fNXX_* function
         → Function: echo banner → read_menu_choice/read_menu_array for args → python/bash runs script
         → Return to main → citation() → script exits
```

### CLI Mode

```
User runs: gpumdkit.sh -flag args
         → case $1 in
         → source / python / bash dispatches to script
         → ${@:2} or ${@:3} forwards arguments
         → exit
```

### Custom Commands (~/.gpumdkit.in)

```
gpumdkit.sh -mycommand args
         → alias_key="${1#-}"        → "mycommand"
         → func_name="custom_${alias_key}" → "custom_mycommand"
         → source ~/.gpumdkit.in
         → custom_mycommand "$@"
```

## GPUMD/NEP Output File Behavior

- **`#` comment lines**: GPUMD now prints `#`-prefixed comment lines at the start of `.out` files (e.g., msd.out). `np.loadtxt()` default `comments='#'` handles these automatically — all 35+ reading scripts are SAFE.
- **`np.savetxt(header=...)`**: auto-prefixes the header with `# `. Writers (like `calc_msd.py`) and readers (like `plt_msd.py`) are pair-compatible.
- **`.out` column counts vary by compute command**: msd.out may have 4 or 7 columns depending on whether `compute_sdc` was enabled. Scripts with hardcoded column assumptions (`plt_sdc.py`, `plt_msd_sdc.py` assume 7) will crash on 4-column files. See "Known Pitfalls" below.

## Known Pitfalls

1. **`plt_sdc.py` / `plt_msd_sdc.py`**: hardcode 7-column msd.out. Crashes with `IndexError` when only `compute_msd` (4 columns) was used. These are the most common script crashes for GPUMD users.
2. **`plt_arrhenius_sigma.py` / `_xyz.py`**: hardcode `z = 1` (species charge). Produces **wrong conductivity values** (off by z²) for multivalent ions (Mg²⁺, Ca²⁺, Al³⁺, etc.).
3. **`plt_emd.py` / `plt_nemd.py` / `plt_hnemd.py` / `plt_pdos.py`**: call `plt.show()` without Agg backend detection. Crash with `TclError` on headless servers/SSH. The other 34 plotting scripts all have proper Agg detection.
4. **`plt_rdf_pmf.py`**: calls `input()` when run without arguments — hangs in non-interactive mode. `gpumdkit.sh -plt rdf_pmf` (no extra args) triggers this.
5. **~30 plotting scripts** use `np.loadtxt()` / `open()` without file-existence checks. Users running from the wrong directory get cryptic tracebacks instead of friendly "file not found" messages.
6. **macOS bash 3.2**: no `local -n` nameref. Use alternative patterns (read -a "$varname", eval, or herestring redirection).

## Design Decisions (recent fixes — do not revert)

- **`calc_ion_conductivity.py` lines 107-108**: MSD fitting range is now 40%-80% (was 10%-40%). The Arrhenius plotting scripts also use 40%-80%. This is a settled scientific choice.
- **`select_max_modev.py` line 111**: The index mapping `atoms_list[filtered_indices[i]]` was recently fixed from a bug (`atoms_list[i]`). The old version silently selected the wrong structures.
- **`charge_balance_check.py` lines 118-131**: Now uses a manual loop instead of `dict()` to handle 3-tuple error returns from `check_oxidation_state()`. The `dict()` constructor crashes on 3-tuples.
- **`clean_extra_files.sh`**: Uses whole-word matching for the keep-exclusion logic, not substring removal. Substring removal (`${var/pattern/}`) corrupts filenames that share substrings with the keep input.
- **`clean_extra_files.sh`**: `read user_input` changed to EOF-safe read with `|| return 1`. On EOF (closed stdin), the old code would proceed to `rm -f` the entire deletion list.
- **`src/f6_plots.sh`**: `doas` is now in the "Diffusion & Transport" section of `f6_plots_two_column` (previously was in "MD & Structural Analysis"). Both one-column and two-column layouts are now consistent.
- **All workflow ATTENTION blocks**: standardized to uniform width (57 dashes inner), with consistent `|` bar edges.

## Verification Checklist

```bash
# Shell syntax
bash -n gpumdkit.sh
find src Scripts -name '*.sh' -exec bash -n {} +

# Python syntax
python3 -m py_compile <modified_file>

# Interactive menu paths
echo "0" | bash gpumdkit.sh                         # Exit
printf '8\n0\n' | bash gpumdkit.sh                   # About/Help → Exit

# CLI paths
bash gpumdkit.sh -h >/dev/null                       # Help table
bash gpumdkit.sh -plt badtype 2>&1                    # Should print "Unknown plot type"
bash gpumdkit.sh -plt -h >/dev/null                   # Plot list

# Documentation build
mkdocs build -f docs/mkdocs.yml
```

## Quick Reference: Skill Files

| Skill | Directory | Purpose |
|-------|-----------|---------|
| Main | `skills/gpumdkit-main/SKILL.md` | Entry point and project overview |
| Format Conversion | `skills/gpumdkit-format-conversion/SKILL.md` | Structure file conversion |
| Calculators | `skills/gpumdkit-calculators/SKILL.md` | Material property calculations |
| Analyzers | `skills/gpumdkit-analyzers/SKILL.md` | Structure analysis and validation |
| Visualization | `skills/gpumdkit-visualization/SKILL.md` | Plotting and visualization |
| Workflows | `skills/gpumdkit-workflows/SKILL.md` | Batch processing and automation |
| Sampling | `skills/gpumdkit-sampling/SKILL.md` | Structure sampling and selection |
| Contributing | `skills/gpumdkit-contributing/SKILL.md` | Conventions and adding new features |
| Agents | `skills/AGENTS.md` | This file — AI agent project knowledge |
