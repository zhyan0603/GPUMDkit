<div align="center">
  <h1>🤝 Contributing to GPUMDkit</h1>
  <p style="text-align: justify;">Thank you for your interest in contributing to <strong>GPUMDkit</strong>! We appreciate your time and effort in helping improve this toolkit.</p>
</div>

<p align="center">
  <strong>English</strong>
  &nbsp;·&nbsp;
  <a href="../zh/贡献指南.md">简体中文</a>
</p>

`GPUMDkit` is an open-source package, and we welcome contributions from the community, whether you're fixing bugs, adding new features, improving documentation, or suggesting enhancements.

> **Note**: The authoritative version of this document is [`CONTRIBUTING.md`](https://github.com/zhyan0603/GPUMDkit/blob/main/CONTRIBUTING.md) in the repository root. This tutorial page mirrors its content.

---

## Table of Contents

- [Development Restrictions and Guidelines](#development-restrictions-and-guidelines)
- [Reporting Bugs](#reporting-bugs)
- [Suggesting Features](#suggesting-features)
- [Contributing Code](#contributing-code)
  - [Fork and Clone the Repository](#fork-and-clone-the-repository)
  - [Create a Feature Branch](#create-a-feature-branch)
  - [Making Changes](#making-changes)
    - [Interactive Mode Contributions](#interactive-mode-contributions)
    - [Command-Line Mode Contributions](#command-line-mode-contributions)
    - [Updating Tab Completion](#updating-tab-completion)
  - [Code Style and Best Practices](#code-style-and-best-practices)
  - [Testing Your Changes](#testing-your-changes)
  - [Commit Messages](#commit-messages)
  - [Push and Open a Pull Request](#push-and-open-a-pull-request)

---

## Development Restrictions and Guidelines

To maintain code quality and consistency across the project, please adhere to these guidelines:

### Language Requirements

- **All code should be written in English**: This includes variable names, function names, comments, docstrings, commit messages, and documentation.
- While we understand that contributors come from diverse backgrounds, using English ensures that the codebase is accessible to the widest possible audience.

### Modularity and Reusability

- **Write modular and reusable code**: Functions and scripts should be designed with flexibility in mind.
- **Prefer passing parameters over hardcoding values**: Use function arguments, command-line parameters, or configuration files instead of hardcoded values. However, some filenames are required to be fixed by GPUMD and NEP programs (e.g., `train.xyz`, `thermo.out`), and these can be hardcoded as needed.
- **Example**: 
  ```bash
  # Good: accepts parameters
  function convert_file() {
      local input_file=$1
      local output_file=$2
      # processing logic
  }
  
  # Avoid: hardcoded values (unless required by GPUMD/NEP)
  function convert_file() {
      local input_file="hardcoded_input.xyz"
      # processing logic
  }
  ```

### Code Style

- Write code that is clear and easy for others to maintain and modify.
- Use meaningful variable and function names.
- Add comments where necessary to explain complex logic.

### Documentation

- Ensure that help messages (e.g., `-h` flags) are clear and accurate.
- Tutorial documentation lives in `docs/tutorials/en/` (English) and `docs/tutorials/zh/` (Chinese).
  If you add a new feature, consider updating the relevant tutorial page.
- After editing tutorial markdown files, rebuild the HTML with `mkdocs build -f docs/mkdocs.yml`.

---

## Reporting Bugs

If you encounter a bug, please help us fix it by reporting it through GitHub Issues:

1. **Search existing issues** to see if the bug has already been reported.
2. **Open a new issue** if it hasn't been reported yet: [Create an Issue](https://github.com/zhyan0603/GPUMDkit/issues/new)
3. **Describe the problem in detail** and provide test files if possible.

You can also reach out to us through:
- **QQ Group**: 825696376
- **Email**: yanzihan@westlake.edu.cn
- **Script Developer**: Contact the developer listed in the script header

---

## Suggesting Features

We welcome feature suggestions! To propose a new feature:

1. **Search existing issues** to see if someone has already suggested it.
2. **Open a new feature request**: [Create a Feature Request](https://github.com/zhyan0603/GPUMDkit/issues/new) and mention @zhyan0603
3. **Describe what you need** clearly and explain how it would be useful.

---

## Contributing Code

### Fork and Clone the Repository

1. **Fork the repository** on GitHub by clicking the "Fork" button at the top right of the [repository page](https://github.com/zhyan0603/GPUMDkit).

2. **Clone your fork** to your local machine:
   
   ```bash
   git clone https://github.com/YOUR_USERNAME/GPUMDkit.git
   cd GPUMDkit
   ```
   
3. **Set up the upstream remote** to keep your fork in sync:
   ```bash
   git remote add upstream https://github.com/zhyan0603/GPUMDkit.git
   ```

### Create a Feature Branch

Create a new branch for your changes:

```bash
# Create and checkout a new branch from main
git checkout -b your-branch-name
```

Name your branch as you prefer.

### Making Changes

The structure of `GPUMDkit` consists of:
- **`gpumdkit.sh`**: Main entry point (Bash script) handling both interactive menu mode and command-line mode
- **`install.sh`**: Installation script that sets up environment variables and shell configuration
- **`Scripts/`**: Python and Bash implementation scripts organized by functionality
  - `plt_scripts/`: Plotting scripts
  - `calculators/`: Calculation utilities
  - `format_conversion/`: Format conversion tools
  - `workflow/`: Workflow automation
  - `sample_structures/`: Structure sampling utilities
  - `analyzer/`: Analysis tools
  - `utils/`: Utility functions (including `completion.sh` for tab completion)
- **`src/`**: Shell scripts containing menu functions for interactive mode
  - `f1_format_conversions.sh`
  - `f2_sample_structures.sh`
  - `f3_workflows.sh`
  - `f4_calculators.sh`
  - `f5_analyzers.sh`
  - `f6_plots.sh`
  - `f7_utilities.sh`
- **`skills/`**: Portable Agent Skills definition (`gpumdkit-skill/SKILL.md`) and on-demand references
- **`docs/`**: Documentation files
  - `tutorials/en/` and `tutorials/zh/`: Bilingual tutorial pages
  - `mkdocs.yml`: MkDocs configuration for building tutorial HTML
  - `command_reference.tsv`: Machine-readable command reference
  - `htmls/`: Generated HTML output from MkDocs

#### Interactive Mode Contributions

To add a new feature accessible through the interactive menu:

1. **Create your implementation script** in the appropriate `Scripts/` subdirectory:
   
   ```bash
   # Example: Add a new format conversion script
   touch Scripts/format_conversion/new_converter.py
   ```
   
2. **Implement your functionality** following the modularity guidelines (accept parameters, avoid hardcoding).

3. **Add a wrapper function** in the corresponding `src/` file:
   ```bash
   # Example: Edit src/f1_format_conversions.sh
   vim src/f1_format_conversions.sh
   ```
   
   Add a new function:
   ```bash
   function f111_new_converter(){
   echo " >-------------------------------------------------<"
   echo " | Calling the script in Scripts/format_conversion |"
   echo " | Script: new_converter.py                        |"
   echo " | Developer: Your Name (your@email.com)           |"
   echo " >-------------------------------------------------<"
   echo " Input <required_param1> <required_param2>"
   echo " Example: input.xyz output.lmp"
   echo " ------------>>"
   read_menu_array converter_args || return 1
   echo " ---------------------------------------------------"
   python ${GPUMDkit_path}/Scripts/format_conversion/new_converter.py "${converter_args[@]}"
   echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/new_converter.py"
   echo " ---------------------------------------------------"
   }
   ```

4. **Update the menu display** in `gpumdkit.sh`:
   ```bash
   # Find the menu() function and update it if needed
   # Find the array_choice array and add your new choice number
   array_choice=(
       "0" "1" "101" "102" "103" "104" "105" "106" "107" "108" "109" "110" "111"  # Added "111"
       # ... rest of choices
   )
   ```

5. **Add the case statement** in `gpumdkit.sh`:
   ```bash
   # In the main() function, find the appropriate case statement
   case "${choice:0:1}" in
       # ...
       "1")
           source ${GPUMDkit_path}/src/f1_format_conversions.sh
           case $choice in
               "1") f1_format_conversion ;;
                "101") f101_out2xyz ;;
                # ... existing cases ...
                "111") f111_new_converter ;;  # Add your case
           esac ;;
       # ...
   esac
   ```

#### Command-Line Mode Contributions

To add a new command-line flag or subcommand:

1. **Create your implementation script** in the appropriate `Scripts/` subdirectory:
   ```bash
   # Example: Add a new analyzer script
   touch Scripts/analyzer/analyze_bonds.py
   ```

2. **Implement your functionality** with clear parameter requirements. Python scripts
   should own detailed argument parsing, help text, type conversion, file checks, and
   user-facing error messages:
   ```python
   # Example: Scripts/analyzer/analyze_bonds.py
   import sys
   
   args = sys.argv[1:]
   if len(args) < 2 or args[0] in ("-h", "--help"):
       print(" Usage: gpumdkit.sh -analyze_bonds <input.xyz> <cutoff_distance>")
       print("    or: python3 analyze_bonds.py <input.xyz> <cutoff_distance>")
       print("")
       print(" Arguments:")
       print("   input.xyz        Input extxyz file")
       print("   cutoff_distance  Bond cutoff distance")
       print("")
       print(" Example: gpumdkit.sh -analyze_bonds structure.xyz 3.0")
       print("")
       sys.exit(0 if args and args[0] in ("-h", "--help") else 1)
   
   input_file = args[0]
   cutoff = float(args[1])
   
   # Your implementation here
   ```
   
   > **Note**: GPUMDkit Python scripts are designed to run directly (not imported as modules).
   > You may use a `main()` function with `if __name__ == "__main__":` guard if you prefer,
   > but it is not required. Most existing scripts execute at module level.

3. **Add the command-line flag handler** in `gpumdkit.sh`. Keep this layer as a
   lightweight router. Detailed parameter validation should stay in the Python script,
   not in the main shell entry point.
   ```bash
   # Find the command-line parsing section (the large "case $1 in" block)
   # Add your new flag in the appropriate location
   
   case $1 in
       # ... existing cases ...
       -analyze_bonds)
           run_python_script "Your Name (your@email.com)" "${analyzer_path}/analyze_bonds.py" "${@:2}" ;;
       # ... rest of cases ...
   esac
   ```

4. **Update the help information** in `gpumdkit.sh`:
   ```bash
   # Find the help_info_table() function and add your command
   function help_info_table(){
       echo "+====================================== Analysis ===============================================+"
       echo "| -analyze_bonds Analyze bond lengths in struct | -analyze_comp Analyze composition of extxyz      |"
       # ... rest of help table ...
   }
   ```

#### Updating Tab Completion

When adding new command-line flags, update the tab completion script:

```bash
# Edit Scripts/utils/completion.sh
vim Scripts/utils/completion.sh
```

Add your new flag to the `opts` variable:
```bash
# Find the line with local opts=...
local opts="-h -help -update -U -clean -time -plt -calc ... -your_new_flag ..."
#                                                        Add your flag here
```

If your flag requires file arguments, add it to the existing file-completion case:
```bash
# Find the case for file-requiring flags and add yours with |
-out2xyz|-out2exyz|-...|-your_new_flag)
    COMPREPLY=($(compgen -f -- "$cur")) ;;
```

If your flag accepts secondary options (like `-plt` or `-calc`), add a new case:
```bash
case "$prev" in
    # ... existing cases ...
    -your_new_flag)
        COMPREPLY=($(compgen -W "option1 option2 option3" -- "$cur")) ;;
    # ... rest of cases ...
esac
```

### Code Style and Best Practices

- **Shell Scripts**:
  - Use `${variable}` for variable expansion
  - Use `[ ! -z "$var" ]` for shell control flow and broad dispatch checks (project convention)
  - Interactive functions should follow the banner format:
    ```bash
    echo " >-------------------------------------------------<"
    echo " | Calling the script in Scripts/<category>        |"
    echo " | Script: <script_name>.py                        |"
    echo " | Developer: <Name> (<email>)                     |"
    echo " >-------------------------------------------------<"
    ```
  - Use `read_menu_choice varname || return 1` or `read_menu_array varname || return 1` for interactive input, then pass arrays with `"${varname[@]}"`
  - Shell scripts in `src/` should start with a file header block:
    ```bash
    # ============================================================
    # GPUMDkit <module name> module
    # Repository: https://github.com/zhyan0603/GPUMDkit
    # Author: <Name> (<email>)
    # ============================================================
    ```

- **Python Scripts**:
  - Write clear, maintainable code
  - Use meaningful variable names
  - Use `args = sys.argv[1:]` for simple positional argument parsing, or `argparse` for complex option sets
  - Provide clear usage messages with leading-space prints such as `print(" Usage: ...")` and `print(" Example: ...")`
  - Put help and missing-argument checks before heavy optional imports so `-h` works without optional packages installed
  - If your script uses heavy/special packages (`NepTrain`, `calorine`, `dpdata`), add a
    `print_dependency_notice()` function to inform users about citation recommendations

### Testing Your Changes

Before submitting your contribution, run the relevant validation commands:

```bash
# Shell syntax checks (always run these)
bash -n gpumdkit.sh
find src Scripts -name '*.sh' -exec bash -n {} +

# Python syntax checks (for modified Python files)
python3 -m py_compile path/to/modified_script.py

# MkDocs build (if you modified documentation)
mkdocs build -f docs/mkdocs.yml

# Check for trailing whitespace issues
git diff --check
```

Then test functionality:

1. **Test interactive mode** (if applicable):
   ```bash
   gpumdkit.sh
   # Navigate to your new feature and test it thoroughly
   ```

2. **Test command-line mode** (if applicable):
   ```bash
   gpumdkit.sh -your_new_flag [arguments]
   gpumdkit.sh -your_new_flag -h  # Test help message
   ```

3. **Test tab completion** (if you modified `completion.sh`):
   ```bash
   source Scripts/utils/completion.sh
   gpumdkit.sh -<TAB>  # Should show your new flag
   ```

4. **Test with various inputs**:
   - Test with typical inputs
   - Test with edge cases (empty files, large files, etc.)
   - Test error handling (missing files, invalid parameters)

5. **Verify no regressions**: Ensure existing functionality still works correctly.

### Commit Messages

Write clear, descriptive commit messages that explain what you changed. There are no strict format requirements - just make sure your message is understandable.

### Push and Open a Pull Request

1. **Commit your changes**:
   ```bash
   git add .
   git commit -m "your descriptive message"
   ```

2. **Keep your branch updated** with the latest `dev` branch:
   ```bash
   git fetch origin
   git rebase origin/dev
   ```

3. **Push your branch** to your fork:
   ```bash
   git push origin your-branch-name
   ```

4. **Open a Pull Request**:
   - Go to the [GPUMDkit repository](https://github.com/zhyan0603/GPUMDkit)
   - Click "Pull Requests" → "New Pull Request"
   - Set the base branch to `dev` for code review
   - Set the compare branch to your feature branch
   - Fill in the PR description with:
     - **Title**: Brief description of changes
     - **Description**: Detailed explanation of what changed and why
     - **Testing**: Describe how you tested the changes

5. **Respond to review feedback**:
   - Be open to suggestions and constructive criticism
   - Make requested changes promptly
   - Push additional commits to your branch (they'll automatically appear in the PR)

6. **After approval**: A maintainer will merge your PR.

---

## Questions or Need Help?

If you have questions about contributing or need help with your contribution:

- **Open a discussion**: Use [GitHub Discussions](https://github.com/zhyan0603/GPUMDkit/discussions) for general questions
- **Ask in an issue**: Comment on related issues for specific questions
- **Contact maintainers**: Reach out to the core developers listed in the README

---

Thank you again for contributing to `GPUMDkit`! Your efforts help make this toolkit better for the entire GPUMD and NEP community.
