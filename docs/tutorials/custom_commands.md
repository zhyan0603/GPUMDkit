# Custom Commands in GPUMDkit

This feature allows users to extend GPUMDkit's command-line interface with **personal custom commands**. You can define shortcuts for frequently used scripts.

## Detailed Examples

### Simple Greeting

Write the following codes in your `~/.bashrc` 

```
custom_hello() {
    echo "Hello, GPUMDkit user! This is a custom command."
    echo "Current GPUMDkit path: $GPUMDkit_path"
}
```

and run `gpumdkit.sh -hello`, you will see:

```
Hello, GPUMDkit user! This is a custom command.
Current GPUMDkit path: /home/yanzihan/software/GPUMDkit
```

### Command with Required Arguments

For the following function:

```
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

If you want to call an external script like gpumdkit.sh.

```
custom_nepanalyse() {
    python ~/my_gpumd_tools/analyse_nep.py
}
```

and run `gpumdkit.sh -nepanalyse`



---

Thank you for using `GPUMDkit`! If you have any questions or need further assistance, feel free to open an issue on our GitHub repository or contact Zihan YAN (yanzihan@westlake.edu.cn).