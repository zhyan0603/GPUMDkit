<div align="center">
  <h1>🔧 Custom Commands</h1>
    <p style="text-align: justify;">Add your own shortcuts to <code>gpumdkit.sh</code> by defining Bash functions in <strong>~/.gpumdkit.in</strong>.</p>
</div>

## How it works

When `gpumdkit.sh -<name>` is not a built-in option, GPUMDkit loads `~/.gpumdkit.in` and runs the function `custom_<name>`, passing it every argument after `-<name>`. Built-in options always take precedence, so choose names that do not clash with them (see `gpumdkit.sh -h`).

Create or edit the file with any editor:

```bash
vi ~/.gpumdkit.in
```

## Basic examples

### A first command

```bash
custom_hello() {
    echo "Hello, GPUMDkit user! This is a custom command."
    echo "Current GPUMDkit path: ${GPUMDkit_path}"
}
```

Running `gpumdkit.sh -hello` prints:

```text
Hello, GPUMDkit user! This is a custom command.
Current GPUMDkit path: /home/yanzihan/software/GPUMDkit
```

### Using arguments

Arguments are available as `$1`, `$2`, ... inside the function:

```bash
custom_greet() {
    echo "Hello, $1!"
}
```

`gpumdkit.sh -greet Zihan` prints `Hello, Zihan!`. To pass all arguments on to another program, use `"$@"`, which keeps spaces and special characters intact.

### Calling your own script

```bash
custom_nepanalyse() {
    python ~/my_gpumd_tools/analyse_nep.py "$@"
}
```

Run it with `gpumdkit.sh -nepanalyse`, followed by any arguments your script accepts.

## Combining GPUMDkit commands

### Plot several directories

```bash
custom_batch_plot() {
    for dir in "$@"; do
        (
            cd "$dir" || exit 1
            gpumdkit.sh -plt thermo save
            gpumdkit.sh -plt msd save
        ) || return 1
    done
}
```

Usage: `gpumdkit.sh -batch_plot dir1 dir2 dir3`

### Check a training set

```bash
custom_analyze_training() {
    xyz_file=$1
    echo "Analyzing $xyz_file..."
    gpumdkit.sh -range "$xyz_file" force
    gpumdkit.sh -min_dist_pbc "$xyz_file"
    gpumdkit.sh -analyze_comp "$xyz_file"
    echo "Analysis complete!"
}
```

Usage: `gpumdkit.sh -analyze_training train.xyz`

### Prepare training data

```bash
custom_prep_training() {
    gpumdkit.sh -out2xyz ./                                  # writes NEPdataset/train.xyz
    gpumdkit.sh -filter_value NEPdataset/train.xyz force 30   # writes filtered.xyz
    gpumdkit.sh -range filtered.xyz energy
}
```

Usage: `gpumdkit.sh -prep_training`. The force threshold `30` is only an example; choose a value suited to your system.

## Tips

- Use descriptive names: `custom_analyze_nep_training` is clearer than `custom_ant`.
- Print a usage line when required arguments are missing.
- Check that input files exist before processing them.

---

Thank you for using `GPUMDkit`! If you have any questions or need further assistance, feel free to open an issue on our GitHub repository or contact Zihan YAN (yanzihan@westlake.edu.cn).
