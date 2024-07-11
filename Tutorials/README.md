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
 GPUMDkit 0.0.0 (dev) (2024-07-11)
 Usage: GPUMDkit -[options]
 Options:
    -outcar2exyz    Convert OUTCAR to nep-exyz file
                    Usage: -outcar2exyz dir_name
                      Examp: gpumdkit.sh -outcar2exyz .

    -castep2exyz    Convert castep to nep-exyz file
                    Usage: -castep2exyz dir_name
                      Examp: gpumdkit.sh -castep2exyz .

    -cp2k2exyz    Convert cp2k output to nep-exyz file
                    Usage: -cp2k2exyz dir_name
                      Examp: gpumdkit.sh -cp2k2exyz .

    -max_rmse         get_max_rmse_xyz
                    Usage: -getmax|-get_max_rmse_xyz train.xyz force_train.out 13

    -h,-help    Show this help message
```

#### Example 2: Converting VASP OUTCARs to extxyz
To convert a VASP OUTCARs to an extended XYZ format (extxyz) file, use the following command:
```sh
gpumdkit.sh -outcar2exyz <dir_of_OUTCARs>
gpumdkit.sh -outcar2exyz .
```



