# Function 3 - Workflow

**Script Location:** `Scripts/workflow/`

This section covers the workflow automation tools in GPUMDkit (Interactive Mode - Function 3).

## Interactive Mode Access

```bash
gpumdkit.sh
# Select: 3) Workflow
```

You'll see the following menu:

```
 ------------>>
 301) SCF batch pretreatment
 302) MD sample batch pretreatment (gpumd)
 303) MD sample batch pretreatment (lmp)
 000) Return to the main menu
 ------------>>
```

---

<div align="center">
  <h1>🔧 Workflow Scripts</h1>
    <p style="text-align: justify;">This directory contains automation scripts for high-throughput computational workflows, particularly for NEP model development and active learning cycles.</p>
</div>

## Overview

Workflow scripts automate repetitive tasks in computational materials research:

- **SCF Batch preprocessin**: Set up VASP/CP2K single-point energy calculations
- **MD sampling**: Prepare molecular dynamics simulations with GPUMD/LAMMPS
- **Active learning**: Iterative NEP model improvement workflows

---

## Via interactive mode

Access workflow tools through `gpumdkit.sh` interactive mode:

```bash
         ____ ____  _   _ __  __ ____  _    _ _
        / ___|  _ \| | | |  \/  |  _ \| | _(_) |_
       | |  _| |_) | | | | |\/| | | | | |/ / | __|
       | |_| |  __/| |_| | |  | | |_| |   <| | |_
        \____|_|    \___/|_|  |_|____/|_|\_\_|\__|

        GPUMDkit Version 1.4.2 (dev) (2025-12-17)
  Core Developer: Zihan YAN (yanzihan@westlake.edu.cn)

 ----------------------- GPUMD -----------------------
 1) Format Conversion          2) Sample Structures
 3) Workflow                   4) Calculators
 5) Analyzer                   6) Developing ...
 0) Quit!
 ------------>>
 Input the function number:
 3
 ------------>>
 301) SCF batch pretreatment
 302) MD sample batch pretreatment (gpumd)
 303) MD sample batch pretreatment (lmp)
 000) Return to the main menu
 ------------>>
 Input the function number:
```

---

### workflow_active_learning_dev.sh

Implements automated active learning cycles for iterative NEP model improvement.

**Purpose:** Automates the active learning pipeline:

1. Generate candidate structures (MD sampling)
2. Evaluate with current NEP model
3. Select uncertain structures
4. Run DFT calculations
5. Add to training set
6. Retrain NEP model

**Status:** Development version - **under active development**

**Documentation:** See detailed tutorial in [docs/tutorials/](../../docs/tutorials/)

---

## Contributing

To add new workflow capabilities, see [CONTRIBUTING.md](../../CONTRIBUTING.md) for detailed guidelines.

---

Thank you for using GPUMDkit! For questions about workflows or to share your workflow adaptations, please open an issue on our [GitHub repository](https://github.com/zhyan0603/GPUMDkit/issues) or contact Zihan YAN (yanzihan@westlake.edu.cn).
