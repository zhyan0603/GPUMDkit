<div align="center">
  <h1>🔧 Custom Commands</h1>
    <p style="text-align: justify;">GPUMDkit allows users to extend the command-line interface with <strong>custom commands</strong>. Define shortcuts for frequently used scripts in your <strong>~/.gpumdkit.in</strong> file.</p>
</div>

## Setup

Create or edit `~/.gpumdkit.in`:
```bash
vi ~/.gpumdkit.in
```

Define custom functions with the `custom_` prefix.

## Examples

### Simple Greeting

Write the following codes in your `~/.gpumdkit.in` 

```bash
custom_hello() {
    echo "Hello, GPUMDkit user! This is a custom command."
    echo "Current GPUMDkit path: ${GPUMDkit_path}"
}
```

and run `gpumdkit.sh -hello`, you will see:

```
Hello, GPUMDkit user! This is a custom command.
Current GPUMDkit path: /home/yanzihan/software/GPUMDkit
```

### Command with Required Arguments

For the following function:

```bash
custom_greet() {
    name=$1
    echo "Hello, ${name}!"
}
```

and you can add a argument, for example: `gpumdkit.sh -greet Zihan`, then you will see:

```
Hello, Zihan!
```

**Passing arguments**: Use `"$@"` to forward all arguments safely. Use `${@:1}` when you need to refer to arguments starting from the first one. 

### Calling an External Script

If you want to call an external script like `gpumdkit.sh`.

```bash
custom_nepanalyse() {
    python ~/my_gpumd_tools/analyse_nep.py
}
```

and run `gpumdkit.sh -nepanalyse`

## Advanced Examples

### Batch Processing

Process multiple files:
```bash
custom_batch_plot() {
    for dir in "$@"; do
        cd "$dir"
        gpumdkit.sh -plt thermo save
        gpumdkit.sh -plt msd save
        cd ..
    done
}
```

Usage: `gpumdkit.sh -batch_plot dir1 dir2 dir3`

### Analysis Pipeline

Combine multiple operations:
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

### Custom Workflow

Automate common workflow:
```bash
custom_prep_training() {
    echo "Preparing training data..."
    
    # Convert VASP outputs
    gpumdkit.sh -out2xyz ./
    
    # Filter outliers
    gpumdkit.sh -filter_value train.xyz force 30
    
    # Check quality
    gpumdkit.sh -range filtered_force.xyz energy
    
    echo "Training data ready!"
}
```

Usage: `gpumdkit.sh -prep_training`

## Tips

- **Use descriptive names**: `custom_analyze_nep_training` is clearer than `custom_ant`
- **Add help messages**: Echo usage info when no arguments provided
- **Error handling**: Check if required files exist before processing
- **Forward arguments safely**: Use `"$@"` to preserve spaces and special characters

---

Thank you for using `GPUMDkit`! If you have any questions or need further assistance, feel free to open an issue on our GitHub repository or contact Zihan YAN (yanzihan@westlake.edu.cn).
