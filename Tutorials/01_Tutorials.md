# function1 - Sample Structures

This tutorial will guide you through using the `sample_structures` function in `GPUMDkit`. For now, this function offers two main features:
1. Converting VASP `OUTCAR` files to `extxyz` format.
2. Sampling structures from an `extxyz` file using a specified method.

## Introduction

When you run the `sample_structures` function, you will be presented with a menu to select one of the available operations. Follow the steps below to use each feature.

### Menu Options

```sh
echo " ------------>>"
echo " 101) Convert OUTCAR to extxyz"
echo " 102) Sample structures from extxyz"
echo " 000) Return to the main menu"
echo " ------------>>"
echo " Input the function number:"
```



### Option 101: Convert OUTCAR to EXYZ

This option allows you to convert VASP OUTCAR files to EXYZ format.

1. Select option `101` from the menu:

   ```sh
   101
   ```

2. You will see the following prompt:

   ```sh
   echo " Input the directory containing OUTCARs"
   ```

3. Enter the directory containing your `OUTCAR` files:

   ```sh
   /path/to/your/outcars
   ```

4. The script   `multipleFrames-outcars2nep-exyz.sh`  in GPUMD's tools will be called to perform the conversion.

**Script Information:**

```sh
echo " >-------------------------------------------------<"
echo " | This function calls the script in GPUMD's tools |"
echo " | Script: multipleFrames-outcars2nep-exyz.sh      |"
echo " | Developer: Yanzhou WANG (yanzhowang@gmail.com ) |"
echo " >-------------------------------------------------<"
```



### Option 102: Sample Structures from extxyz

This option allows you to sample structures from an `extxyz` file using a specified method.

1. Select option `102` from the menu:

   ```sh
   102
   ```

2. You will see the following prompt:

   ```sh
   echo " Input <extxyz_file> <sampling_method> <num_samples>"
   echo " Sampling_method: 'uniform' or 'random'"
   echo " Examp: train.xyz uniform 50"
   ```

3. Enter the `extxyz` file name, sampling method, and number of samples:

   ```sh
   train.xyz uniform 50
   ```

4. The script `sample_structures.py` in the `Scripts` will be called to perform the sampling.

**Script Information:**

```sh
echo " >-------------------------------------------------<"
echo " | This function calls the script in Scripts       |"
echo " | Script: sample_structures.py                    |"
echo " | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |"
echo " >-------------------------------------------------<"
```



### Option 000: Return to the Main Menu

If you wish to return to the main menu, select option `000`:

```sh
000
```



Thank you for using `GPUMDkit`! If you have any questions or need further assistance, feel free to open an issue on our GitHub repository or contact Zihan YAN (yanzihan@westlake.edu.cn).