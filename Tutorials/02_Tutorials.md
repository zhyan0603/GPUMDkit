# Function2 - Sample Structures

This script provides a menu-driven interface to perform various tasks related to structure sampling.

### Menu Options

```sh
------------>>
201) Sample structures from extxyz
202) Find the outliers in training set
000) Return to the main menu
------------>>
Input the function number:
```

#### Option 201: Sample Structures from extxyz

This option allows you to sample structures from an `extxyz` file using a specified method.

1. Select option `201` from the menu:

   ```sh
   201
   ```

2. You will see the following prompt:

   ```sh
   >-------------------------------------------------<
   | This function calls the script in Scripts       |
   | Script: sample_structures.py                    |
   | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
   >-------------------------------------------------<
   Input <extxyz_file> <sampling_method> <num_samples>
   Sampling_method: 'uniform' or 'random'
   Examp: train.xyz uniform 50
   ------------>>
   ```

3. Enter the `extxyz` file name, sampling method, and number of samples:

   ```sh
   train.xyz uniform 50
   ```

4. The script `sample_structures.py` in the `Scripts` will be called to perform the sampling.

#### Option 202: Find the outliers in training set

This function calls the `get_max_rmse_xyz.py` script to find outliers in a training set.

1. Select option `202` from the menu:

   ```sh
   202
   ```

2. You will see the following prompt:

   ```sh
   >-------------------------------------------------<
   | This function calls the script in GPUMD's tools |
   | Script: get_max_rmse_xyz.py                     |
   | Developer: Ke XU (kickhsu@gmail.com)            |
   >-------------------------------------------------<
   Input <extxyz_file> <*_train.out> <num_outliers>
   Examp: train.xyz energy_train.out 13 
   ------------>>
   ```

   `<extxyz_file>`: extxyz file

   `<*_train.out>`: `energy_train.out`/`force_train.out`/`virial_train.out``

   `<num_outliers>`: number of outliers

3. Enter the `extxyz` file name, <*_train.out>, and number of outliers:

```sh
train.xyz energy_train.out 13 
```

4. The script `sample_structures.py` in the `Scripts` will be called to perform the sampling.



Thank you for using `GPUMDkit`! If you have any questions or need further assistance, feel free to open an issue on our GitHub repository or contact Zihan YAN (yanzihan@westlake.edu.cn).

