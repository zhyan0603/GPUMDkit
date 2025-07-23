# GPUMDkit Tutorials

Welcome to the `GPUMDkit` tutorials! This guide will help you get started with using `GPUMDkit`, a shell interface designed to streamline your `GPUMD` workflows.

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
