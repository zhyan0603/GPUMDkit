# GPUMDkit

**`GPUMDkit`** is a shell interface for the [GPUMD](https://github.com/brucefan1983/GPUMD) (*Graphics Processing Units Molecular Dynamics*) program. It provides a set of tools to streamline the use of common scripts in GPUMD and  [NEP](https://gpumd.org/potentials/nep.html#nep-formalism) (neuroevolution potential), simplifying workflows and enhancing efficiency.

## Features
- **Simplified Script Invocation**: Easily run scripts for GPUMD and NEP.
- **Workflow Automation**: Automate common tasks to save time and reduce manual intervention.
- **User-Friendly Interface**: Intuitive shell commands designed to enhance user experience.

## Installation
To install `GPUMDkit`, follow these steps:

1. Clone the repository or download the whole project.

2. Set the `GPUMD_path` and `GPUMDkit_path` variables in your `~/.bashrc` file, for example:
   
    ```
    vi ~/.bashrc
    ```
    
    add these two variables
    
    ```sh
    export GPUMD_path=/your_dir_of_GPUMD
    export GPUMDkit_path=/your_dir_of_GPUMDkit
    ```
    
    add `GPUMDkit_path` to the PATH

    ```sh
    export PATH=/your_dir_of_GPUMDkit:${PATH}
    ```

    then
    
    ```sh
    source ~/.bashrc
    ```
    
3. Add executable permissions to the `gpumdkit.sh` file:
    ```sh
    chmod +x gpumdkit.sh
    ```
    
4. Move the `gpumdkit.sh` file to a directory in your PATH, for example:
    ```sh
    mv gpumdkit.sh ~/bin/
    ```

## Usage
Refer to the [Usage Instructions](./Tutorials/README.md) for detailed examples and command options.



---

Thank you for using `GPUMDkit`! If you have any questions or need further assistance, feel free to open an issue on our GitHub repository or contact Zihan YAN (yanzihan@westlake.edu.cn).
