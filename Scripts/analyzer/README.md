# Analyzer Scripts

This directory contains analysis tools for structure files, simulation data, and dataset quality control. These scripts help validate, filter, and understand molecular dynamics data and training sets.

## Overview

The analyzer scripts provide functionality for:
- Energy, force, and virial range analysis
- Minimum distance calculations (with/without periodic boundary conditions)
- Structure filtering by various criteria (distance, box size, property values)
- Composition analysis of multi-component systems
- Time estimation for GPUMD and NEP calculations
- Dataset quality checks and outlier detection

Access analyzers through `gpumdkit.sh` using various flags or run scripts directly.

---

## Scripts

### energy_force_virial_analyzer.py

---

This script calculates and visualizes the range of properties (such as energy, forces, and virial) from the `extxyz` file.

#### Usage

The script requires at least two arguments:

- The filename of the `extxyz` file.
- The name of the property to analyze (`energy`, `force`, or `virial`).

An optional third argument (`hist`) can be provided to generate a histogram plot of the values.

```
python energy_force_virial_analyzer.py <filename> <property> [hist]
```

#### Example

```sh
python energy_force_virial_analyzer.py dump.xyz force
```

#### Command-Line Mode Example

```
gpumdkit.sh -range dump.xyz force
```

#### Output

```
Force range: 0.03210566767721861 to 9.230115912468435
```

If you add the `[hist]` option, it will calculate the range of forces and display a histogram:

```sh
python energy_force_virial_analyzer.py dump.xyz force hist
```

<div align="center">
<img src="../../docs/Gallery/range_force.png" width = "50%" />
</div>



### filter_dist_range.py

---

This script is used to extract the structures with specified `min_dist`.

#### Usage

```
python filter_dist_range.py <input.xyz> <element1> <element2> <min_dist> <max_dist>
```

#### Example

```sh
python filter_dist_range.py dump.xyz Li Li 1.8 2.0
```

#### Command-Line Mode Example

```
gpumdkit.sh -filter_range dump.xyz Li Li 1.9 2.0
```

This means you need to extract the structures with the `min_dist` of Li-Li in the range of 1.8-2.0 in `dump.xyz` file. Finally, the `filtered_Li_Li_1.8_2.0.xyz` file will be generated.



### get_min_dist.py

---

This script is used to calculate the min_dist of the structures.

#### Usage

```
python get_min_dist.py <extxyz_file>
```

#### Example

```sh
python get_min_dist.py dump.xyz
```

#### Command-Line Mode Example

```
gpumdkit.sh -min_dist dump.xyz
```

#### Output

```
Minimum interatomic distance: 1.478098603206159 Å
```

NOTE: This script is fast because it does not take into account periodic boundary conditions (PBC), but in some cases it can be problematic.



### get_min_dist_pbc.py

---

This script is used to calculate the min_dist of the structures considering the PBC.

#### Usage

```
python get_min_dist_pbc.py <extxyz_file>
```

#### Example

```sh
python get_min_dist_pbc.py dump.xyz
```

#### Command-Line Mode Example

```
gpumdkit.sh -min_dist_pbc dump.xyz
```

#### Output

```
Minimum interatomic distance: 1.478098603206159 Å
```



### filter_exyz_by_value.py

---

This script filter the structures by min_dist.

#### Usage

```sh
python filter_exyz_by_value.py <extxyz_file> <min_dist>
```

#### Example

```sh
python filter_exyz_by_value.py dump.xyz 1.4
```

#### Command-Line Mode Example

```
gpumdkit.sh -filter_value dump.xyz 1.4
```

#### 

### filter_exyz_by_box.py

---

This script filter the structures by box limit.

#### Usage

```
python filter_exyz_by_box.py <extxyz_file> <min_dist>
```

#### Example

```
python filter_exyz_by_box.py dump.xyz 20
```

#### Command-Line Mode Example

```
gpumdkit.sh -filter_box dump.xyz 20
```



### filter_exyz_by_value.py

---

This script filter the structures by specified value.

#### Usage

```
python filter_exyz_by_value.py <extxyz_file> <property> <threshold>
```

- `<extxyz_file>`: The path to the input `extxyz` file.
- `<property>`: Filtering property: `energy`, `force`, or `virial`
- `<threshold>`: Threshold value for filtering

#### Example

```
python filter_exyz_by_value.py train.xyz force 20
```

#### Command-Line Mode Example

```
gpumdkit.sh -filter_value train.xyz force 20
```

This command will filter out the structure in `train.xyz` with a force greater than 20 eV/angstrom.



### time_consuming_gpumd.sh

---

This script calculates the remaining time for GPUMD.

#### Usage

```
bash time_consuming_gpumd.sh <logfile>
```

- `<logfile>`: The path to your `log` file.

#### Example

```
bash time_consuming_gpumd.sh log
```

#### Command-Line Mode Example

```
gpumdkit.sh -time gpumd <logfile>
```

#### Output

```
------------------ Time Consuming Results ------------------
num of atoms: 7168
atom*step/s : 4.85401e+06
timesteps/s : 677.178
total frames: 1050000
total time  : 0h 25min 50s
time left   : 0h 0min 0s
Progress Bar: [########################################] 100%
```



### time_consuming_nep.sh

---

This script calculates the remaining time for nep.

#### Usage

```
bash time_consuming_gpumd.sh
```

#### Example

```
bash time_consuming_nep.sh
```

#### Command-Line Mode Example

```
gpumdkit.sh -time nep
```

#### Output

```
+-----------------+-----------+-----------------+---------------------+
|       Step      | Time Diff |    Time Left    |    Finish Time      |
+-----------------+-----------+-----------------+---------------------+
| 6700            | 1 s       | 0 h 15 m 33 s   | 2025-10-23 15:34:11 |
| 6800            | 2 s       | 0 h 31 m 4 s    | 2025-10-23 15:49:44 |
| 6900            | 2 s       | 0 h 31 m 2 s    | 2025-10-23 15:49:44 |
| 7000            | 2 s       | 0 h 31 m 0 s    | 2025-10-23 15:49:44 |
| 7100            | 3 s       | 0 h 46 m 27 s   | 2025-10-23 16:05:14 |
```



---

## File Naming Conventions

Scripts automatically generate output files with descriptive names:

- `filtered_<property>_<threshold>.xyz` - Filtered by property
- `filtered_<elem1>_<elem2>_<min>_<max>.xyz` - Filtered by distance range
- `<composition>.xyz` - Exported by composition

Keep track of filtering history in a log file for reproducibility.

## Integration with Other Tools

### With Sample Structures
```bash
# 1. Filter dataset
gpumdkit.sh -filter_value train.xyz force 30

# 2. Sample from filtered set
python sample_structures.py filtered_force.xyz uniform 1000
```

### With NEP Training
```bash
# 1. Clean training set
gpumdkit.sh -filter_dist train.xyz 1.5
gpumdkit.sh -filter_value filtered_dist.xyz force 25

# 2. Analyze final set
gpumdkit.sh -analyze_comp filtered_force.xyz
gpumdkit.sh -range filtered_force.xyz force

# 3. Use in NEP training
# (Use filtered files as training input)
```

### With Plotting
```bash
# 1. Analyze with histogram
gpumdkit.sh -range train.xyz force hist

# 2. Further analysis in Python
# (Use generated images or data)
```

## Contributing

To add new analyzer scripts:

1. **Follow naming**: `<descriptive_name>.py`
2. **Handle errors**: Validate inputs gracefully
3. **Progress indicators**: Use tqdm for long operations
4. **Document thoroughly**: Include docstrings and usage examples
5. **Update README**: Add to this documentation
6. **Test comprehensively**: Try various input formats and edge cases

See [CONTRIBUTING.md](../../CONTRIBUTING.md) for detailed guidelines.

---

Thank you for using GPUMDkit! If you have questions or need assistance with analyzer scripts, please open an issue on our [GitHub repository](https://github.com/zhyan0603/GPUMDkit/issues) or contact Zihan YAN (yanzihan@westlake.edu.cn).
