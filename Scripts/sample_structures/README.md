<div align="center">
  <h1>📦 Structure Sampling Scripts</h1>
    <p style="text-align: justify;">This directory contains tools for sampling, selecting, and generating atomic structures for NEP training and molecular dynamics simulations.</p>
</div>

## Overview

Structure sampling is crucial for creating diverse configurations. These scripts provide:
- **Uniform and random sampling**: Basic sampling methods
- **FPS (Farthest Point Sampling)**: Maximizing structural diversity using PyNEP or NEPtrain
- **Structure perturbation**: Generating variations for training set augmentation
- **maximum force deviations selection**: Identifying structures with maximum force deviations
- **Frame extraction**: Selecting specific range from trajectories

---

## Via interactive mode

---

```
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
 2
 ------------>>
 201) Sample structures from extxyz
 202) Sample structures by pynep
 203) Sample structures by neptrain
 204) Perturb structure
 205) Select max force deviation structs
 000) Return to the main menu
 ------------>>
 Input the function number:
```

## Contributing

See [CONTRIBUTING.md](../../CONTRIBUTING.md) for detailed guidelines.

---

Thank you for using GPUMDkit! If you have questions about structure sampling, please open an issue on our [GitHub repository](https://github.com/zhyan0603/GPUMDkit/issues) or contact Zihan YAN (yanzihan@westlake.edu.cn).
