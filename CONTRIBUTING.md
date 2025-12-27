# Contributing to GPUMDkit

Thank you for your interest in contributing to **GPUMDkit**! We appreciate your time and effort in helping improve this toolkit for GPUMD and NEP. This document provides guidelines to make the contribution process smooth and effective for everyone involved.

GPUMDkit is an open-source scientific software project, and we welcome contributions from the community, whether you're fixing bugs, adding new features, improving documentation, or suggesting enhancements.

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

- **All code must be written in English**: This includes variable names, function names, comments, docstrings, commit messages, and documentation.
- While we understand that contributors come from diverse backgrounds, using English ensures that the codebase is accessible to the widest possible audience.

### Modularity and Reusability

- **Write modular and reusable code**: Functions and scripts should be designed with flexibility in mind.
- **Prefer passing parameters over hardcoding values**: Use function arguments, command-line parameters, or configuration files instead of hardcoded values.
- **Example**: 
  ```bash
  # Good: accepts parameters
  function convert_file() {
      local input_file=$1
      local output_file=$2
      # processing logic
  }
  
  # Avoid: hardcoded values
  function convert_file() {
      local input_file="hardcoded_input.xyz"
      # processing logic
  }
  ```

### Code Style

- **Shell scripts**: Follow consistent indentation (4 spaces) and use descriptive function names.
- **Python scripts**: Follow PEP 8 style guidelines. Use meaningful variable names and add docstrings to functions and classes.
- **Comments**: Write clear, concise comments explaining *why* something is done, not just *what* is done (the code should be self-explanatory for the "what").

### Documentation

- Update documentation when adding new features or changing existing behavior.
- Ensure that help messages (e.g., `-h` flags) are clear and accurate.
- If your change affects user workflows, consider updating the relevant files in the `docs/` directory.

---

## Reporting Bugs

If you encounter a bug, please help us fix it by reporting it through GitHub Issues:

1. **Search existing issues** to see if the bug has already been reported.
2. **Open a new issue** if it hasn't been reported yet: [Create an Issue](https://github.com/zhyan0603/GPUMDkit/issues/new)
3. **Provide detailed information**:
   - **Description**: Clearly describe the bug and its impact.
   - **Steps to reproduce**: List the exact steps to reproduce the issue.
   - **Expected behavior**: What you expected to happen.
   - **Actual behavior**: What actually happened.
   - **Environment**: Your operating system, shell version, Python version, and GPUMDkit version.
   - **Error messages**: Include any error messages or stack traces.
   - **Screenshots**: If applicable, add screenshots to illustrate the problem.

**Example Issue Title**: `Bug: -out2xyz fails with OUTCAR files containing magnetic moments`

---

## Suggesting Features

We welcome feature suggestions! To propose a new feature:

1. **Search existing issues** to see if someone has already suggested it.
2. **Open a new feature request**: [Create a Feature Request](https://github.com/zhyan0603/GPUMDkit/issues/new)
3. **Provide detailed information**:
   - **Description**: Clearly describe the feature and its purpose.
   - **Use case**: Explain why this feature would be useful and how it would be used.
   - **Proposed implementation**: If you have ideas about how to implement it, share them.
   - **Examples**: Provide examples or mockups if applicable.

**Example Feature Request Title**: `Feature: Add support for converting LAMMPS trajectory to extxyz format`

---

## Contributing Code

### Fork and Clone the Repository

1. **Fork the repository** on GitHub by clicking the "Fork" button at the top right of the [repository page](https://github.com/zhyan0603/GPUMDkit).

2. **Clone your fork** to your local machine:
   ```bash
   git clone https://github.com/YOUR_USERNAME/GPUMDkit.git
   cd GPUMDkit
   ```

3. **Checkout the `refactor` branch** (this is our main development branch):
   ```bash
   git checkout refactor
   ```

4. **Set up the upstream remote** to keep your fork in sync:
   ```bash
   git remote add upstream https://github.com/zhyan0603/GPUMDkit.git
   ```

### Create a Feature Branch

Always create a new branch for your changes. Branch names should be descriptive:

```bash
# Create and checkout a new branch from refactor
git checkout -b feature/add-xyz-to-lammps-converter refactor
```

**Branch naming conventions**:
- Feature additions: `feature/description-of-feature`
- Bug fixes: `fix/description-of-bug`
- Documentation: `docs/description-of-change`
- Refactoring: `refactor/description-of-change`

### Making Changes

The structure of GPUMDkit consists of:
- **`gpumdkit.sh`**: Main entry point (Bash script) handling both interactive menu mode and command-line mode
- **`Scripts/`**: Python utility scripts organized by functionality
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
- **`docs/`**: Documentation files

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
   function f106_new_converter(){
   echo " >-------------------------------------------------<"
   echo " | Calling the script in Scripts/format_conversion |"
   echo " | Script: new_converter.py                        |"
   echo " | Developer: Your Name (your@email.com)           |"
   echo " >-------------------------------------------------<"
   echo " Input the required parameters:"
   echo " Examp: input.xyz output.lmp"
   echo " ------------>>"
   read -p " " converter_params
   echo " ---------------------------------------------------"
   python ${GPUMDkit_path}/Scripts/format_conversion/new_converter.py ${converter_params}
   echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/new_converter.py"
   echo " ---------------------------------------------------"
   }
   ```

4. **Update the menu display** in `gpumdkit.sh`:
   ```bash
   # Find the menu() function and update it if needed
   # Find the array_choice array and add your new choice number
   array_choice=(
       "0" "1" "101" "102" "103" "104" "105" "106"  # Added "106"
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
               "106") f106_new_converter ;;  # Add your case
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

2. **Implement your functionality** with proper argument parsing:
   ```python
   # Example: Scripts/analyzer/analyze_bonds.py
   import sys
   import argparse
   
   def main():
       parser = argparse.ArgumentParser(description='Analyze bond lengths in structure')
       parser.add_argument('input_file', help='Input XYZ file')
       parser.add_argument('--cutoff', type=float, default=3.0, help='Bond cutoff distance')
       args = parser.parse_args()
       
       # Your implementation here
       
   if __name__ == "__main__":
       main()
   ```

3. **Add the command-line flag handler** in `gpumdkit.sh`:
   ```bash
   # Find the command-line parsing section (after line ~167)
   # Add your new flag in the appropriate location
   
   case $1 in
       # ... existing cases ...
       -analyze_bonds)
           if [ ! -z "$2" ]; then
               python ${analyzer_path}/analyze_bonds.py "$2" "${@:3}"
           else
               echo " Usage: gpumdkit.sh -analyze_bonds <input.xyz> [--cutoff 3.0]"
               echo " Code path: ${analyzer_path}/analyze_bonds.py"
               exit 1
           fi ;;
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
local opts="-h -update -U -help -clean -time -plt -calc -analyze_bonds -range -out2xyz ..."
#                                                      ^^^^^^^^^^^^^^ Add your flag here
```

If your flag accepts secondary options (like `-plt` or `-calc`), add a case for it:
```bash
case "$prev" in
    # ... existing cases ...
    -analyze_bonds)
        COMPREPLY=($(compgen -f -- "$cur"))  # Complete with filenames
        ;;
    # ... rest of cases ...
esac
```

### Code Style and Best Practices

- **Shell Scripts**:
  - Use `#!/bin/bash` shebang
  - Use 4-space indentation
  - Quote variables: `"$variable"` not `$variable`
  - Use `[[` instead of `[` for conditionals
  - Add error checking: `set -e` or check return codes

- **Python Scripts**:
  - Follow PEP 8 style guide
  - Use type hints where appropriate
  - Add docstrings to functions and classes:
    ```python
    def calculate_distance(coord1, coord2):
        """
        Calculate Euclidean distance between two coordinates.
        
        Args:
            coord1 (list): First coordinate [x, y, z]
            coord2 (list): Second coordinate [x, y, z]
            
        Returns:
            float: Distance between the two coordinates
        """
        # implementation
    ```
  - Use meaningful variable names (avoid single letters except for common iterators)

### Testing Your Changes

Before submitting your contribution:

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

Write clear, descriptive commit messages following conventional commit format:

**Format**:
```
<type>(<scope>): <short summary>

<detailed description (optional)>

<footer (optional)>
```

**Types**:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `style`: Code style changes (formatting, no functional changes)
- `refactor`: Code refactoring
- `test`: Adding or updating tests
- `chore`: Maintenance tasks

**Examples**:
```bash
feat(format_conversion): add LAMMPS trajectory to XYZ converter

Implement new converter that reads LAMMPS dump files and converts
them to extended XYZ format. Supports custom atom types and periodic
boundary conditions.

Closes #123
```

```bash
fix(calculators): correct ionic conductivity calculation for multi-species systems

The previous implementation incorrectly averaged conductivities.
Now properly weights by species concentration.
```

```bash
docs(README): update installation instructions for Python dependencies
```

**Guidelines**:
- Use present tense: "add feature" not "added feature"
- Keep the first line under 72 characters
- Reference related issues: `Fixes #123`, `Closes #456`, `Related to #789`
- Write in English
- Be descriptive but concise

### Push and Open a Pull Request

1. **Commit your changes**:
   ```bash
   git add .
   git commit -m "feat(scope): your descriptive message"
   ```

2. **Keep your branch updated** with the upstream `refactor` branch:
   ```bash
   git fetch upstream
   git rebase upstream/refactor
   ```

3. **Push your branch** to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```

4. **Open a Pull Request**:
   - Go to the [GPUMDkit repository](https://github.com/zhyan0603/GPUMDkit)
   - Click "Pull Requests" → "New Pull Request"
   - Set the base branch to `dev` (not `refactor` directly)
   - Set the compare branch to your feature branch
   - Fill in the PR template with:
     - **Title**: Brief description of changes (following commit message format)
     - **Description**: Detailed explanation of what changed and why
     - **Related Issues**: Link any related issues
     - **Testing**: Describe how you tested the changes
     - **Checklist**: Complete the checklist items

5. **Respond to review feedback**:
   - Be open to suggestions and constructive criticism
   - Make requested changes promptly
   - Push additional commits to your branch (they'll automatically appear in the PR)

6. **After approval**: A maintainer will merge your PR into the `dev` branch, and it will eventually be merged into `refactor`.

---

## Questions or Need Help?

If you have questions about contributing or need help with your contribution:

- **Open a discussion**: Use [GitHub Discussions](https://github.com/zhyan0603/GPUMDkit/discussions) for general questions
- **Ask in an issue**: Comment on related issues for specific questions
- **Contact maintainers**: Reach out to the core developers listed in the README

---

## Code of Conduct

We are committed to providing a welcoming and inclusive environment. Please:

- Be respectful and considerate
- Welcome newcomers and help them get started
- Give and accept constructive feedback gracefully
- Focus on what is best for the community and the project
- Show empathy towards other community members

---

Thank you again for contributing to GPUMDkit! Your efforts help make this toolkit better for the entire GPUMD and NEP community. 🚀
