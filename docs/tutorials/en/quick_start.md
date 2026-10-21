<div align="center">
  <h1>🚀 Quick Start</h1>
  <p style="text-align: justify;">Install GPUMDkit, check that it works, and run a first conversion.</p>
</div>

## Install {#before-you-start}

### With Conda (recommended) {#install-with-conda-recommended}

```bash
conda create -n gpumdkit -c gpumdkit -c conda-forge gpumdkit
conda activate gpumdkit
```

This installs GPUMDkit and its standard dependencies. No repository clone or environment variables are needed.

### From source {#install-from-source-optional}

Use this route if you want to read or modify the code:

```bash
git clone https://github.com/zhyan0603/GPUMDkit.git
cd GPUMDkit
source ./install.sh
```

The installer writes `GPUMDkit_path` and `PATH` to `~/.zshrc` if `$SHELL` is zsh, otherwise to `~/.bashrc`, and loads them into the current shell. It backs up the rc file first, and asks before replacing an existing GPUMDkit path.

## Check the installation {#verify-the-installation}

```bash
gpumdkit.sh -h
gpumdkit.sh -doctor
```

`-h` lists all options. `-doctor` checks the Python environment and reports common and optional packages. `MISS` for an optional package is not an error; install it only when you need the feature that uses it.

## What a healthy `-doctor` looks like

`gpumdkit.sh -doctor` prints one line per package. Your versions will differ; what matters is that nothing in **Environment** or **Common packages** is missing:

```text
 Environment
   [OK]   GPUMDkit_path    /path/to/GPUMDkit
   [OK]   Python           3.9.19
   [OK]   Bash             3.2.57(1)-release

 Common packages
   [OK]   numpy            1.26.4
   [OK]   scipy            1.13.1
   [OK]   matplotlib       3.9.2
   [OK]   ase              3.23.0
   [OK]   pandas           2.2.3
   [OK]   pymatgen         2024.8.9

 Optional packages
   [MISS] dpdata           DeepMD data conversion and perturbation
   [MISS] NepTrain         NepTrain structure sampling
   [OK]   seaborn          0.13.2
```

`[MISS]` under **Optional packages** only means that the functions listed next to the name are unavailable until you install the package. A `[MISS]` under **Common packages** should be fixed first.

## First practice: convert a POSCAR {#hello-world-example}

Run this in a new directory so no existing files are overwritten:

```bash
mkdir gpumdkit-first-check && cd gpumdkit-first-check
cat > POSCAR << 'EOF'
Si
1.0
0 2.715 2.715
2.715 0 2.715
2.715 2.715 0
Si
1
direct
0 0 0
EOF
gpumdkit.sh -pos2exyz POSCAR model.xyz
head model.xyz
```

`model.xyz` should list one Si atom together with the lattice.

`extxyz` (extended XYZ) is the training-data format used by NEP: line 1 is the atom count, line 2 holds the lattice and structure-level properties, and each following line holds one atom's element, coordinates, and per-atom properties.

## Two ways to run GPUMDkit {#two-modes}

### Interactive mode {#interactive-mode}

Run `gpumdkit.sh` and choose a module by number; each module has its own sub-menu. At file-name prompts, `Tab` completes names and the arrow keys edit the line.

### Command-line mode {#cli-mode}

Pass an option and its positional arguments directly, which also works in scripts:

```bash
gpumdkit.sh -<option> [args...]
```

`gpumdkit.sh -h` lists every option and `gpumdkit.sh -plt -h` lists the plot types. For Python-backed options, `gpumdkit.sh -<option> -h` shows the arguments of that command.

## Reading the examples {#read-a-command-before-running-it}

| Notation | Meaning | Example |
|---|---|---|
| `<required>` | A value you must provide; do not type the angle brackets. | `<output.xyz>` → `model.xyz` |
| `[optional]` | Can be omitted; the script then uses its documented default. | `[max_corr_steps]` |
| `...` | One or more values. | `<element...>` |

Time parameters such as `dt_fs` must match how your own trajectory was written; GPUMDkit does not read them from the MD input file.

## Optional packages and agent skills {#optional-packages}

Some features need extra packages:

```bash
pip install neptrain calorine
```

GPUMDkit also ships English and Chinese agent skills. To find where they are installed and how to add them to your agent client, run:

```bash
gpumdkit.sh -skill
```

## Next steps {#common-examples}

If you already have structures or trajectories, pick a task from the [tutorial index](index.md): format conversion, sampling, analysis, calculators, plotting, or workflows. Each page lists its inputs, outputs, and limits.
