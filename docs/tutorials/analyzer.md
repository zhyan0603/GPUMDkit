# Function 5 - Analyzer Tools

This section covers the analyzer tools in GPUMDkit (Interactive Mode - Function 5), which provide comprehensive analysis capabilities for structure files, simulation data, and dataset quality control.

## Interactive Mode Access

```bash
gpumdkit.sh
# Select: 5) Analyzer
```

You'll see the following menu:

```
 ------------>>
 501) Analyze composition of extxyz
 502) Find outliers of extxyz
 000) Return to the main menu
 ------------>>
```

## Command-Line Usage

Most analyzer functions are accessible via command-line:

```bash
# Range analysis
gpumdkit.sh -range <file.xyz> <property>

# Minimum distance
gpumdkit.sh -min_dist <file.xyz>
gpumdkit.sh -min_dist_pbc <file.xyz>

# Filtering
gpumdkit.sh -filter_value <file.xyz> <property> <threshold>
gpumdkit.sh -filter_dist <file.xyz> <min_dist>
gpumdkit.sh -filter_box <file.xyz> <box_limit>

# Composition analysis
gpumdkit.sh -analyze_comp <file.xyz>

# Time estimation
gpumdkit.sh -time gpumd
gpumdkit.sh -time nep
```

---

## Available Analyzers

### Property Range Analysis (`energy_force_virial_analyzer.py`)

**Purpose:** Analyze the distribution and range of energy, force, or virial values in training datasets.

**Command-line:**
```bash
# Show range only
gpumdkit.sh -range train.xyz force
gpumdkit.sh -range train.xyz energy
gpumdkit.sh -range train.xyz virial

# Show range with histogram
python ${GPUMDkit_path}/Scripts/analyzer/energy_force_virial_analyzer.py train.xyz force hist
```

**Example Output:**
```
Force range: 0.032 to 9.230 eV/Å
```

**With histogram option:**

<div align="center">
<img src="../Gallery/range_force.png" width="50%" />
</div>

**Use Cases:**
- Verify training data quality
- Identify outliers before training
- Check for unphysical values
- Compare datasets

---

### Minimum Distance Calculation

#### Without PBC (`get_min_dist.py`)

Fast calculation ignoring periodic boundary conditions.

**Command-line:**
```bash
gpumdkit.sh -min_dist structure.xyz
```

**Output:**
```
+---------------------------+
|   PBC ignored for speed   |
| use -min_dist_pbc for PBC |
+---------------------------+
Minimum interatomic distances:
+---------------------------+
| Atom Pair |  Distance (Å) |
+---------------------------+
|   Li-Li   |     1.696     |
|   Li-La   |     2.498     |
|   Li-Zr   |     2.506     |
|   Li-O    |     1.587     |
|   La-La   |     3.463     |
|   La-Zr   |     3.243     |
|   La-O    |     2.043     |
|   Zr-Zr   |     5.060     |
|   Zr-O    |     1.867     |
|   O-O     |     2.480     |
+---------------------------+
Overall min_distance: 1.587 Å
```

**Advantages:** Very fast for large systems
**Limitations:** May give incorrect results for small periodic cells

#### With PBC (`get_min_dist_pbc.py`)

Accurate calculation considering periodic boundary conditions.

**Command-line:**
```bash
gpumdkit.sh -min_dist_pbc structure.xyz
```

**Output:**
```
Minimum interatomic distances (with PBC):
+---------------------------+
| Atom Pair |  Distance (Å) |
+---------------------------+
|   Li-Li   |     1.696     |
|   Li-La   |     2.498     |
|   Li-Zr   |     2.477     |
|   Li-O    |     1.587     |
|   La-La   |     3.463     |
|   La-Zr   |     3.210     |
|   La-O    |     2.043     |
|   Zr-Zr   |     5.060     |
|   Zr-O    |     1.867     |
|   O-O     |     2.355     |
+---------------------------+
Overall min_distance: 1.587 Å
```

**Recommendation:** Always use `-min_dist_pbc` for periodic systems to ensure accuracy.

---

### Structure Filtering

#### Filter by Property Value (`filter_exyz_by_value.py`)

Remove structures with property values exceeding thresholds.

**Command-line:**
```bash
# Filter structures with forces > 30 eV/Å
gpumdkit.sh -filter_value train.xyz force 30

# Filter structures with energy per atom > 0.5 eV
gpumdkit.sh -filter_value train.xyz energy 0.5

# Filter structures with virial > 100 GPa
gpumdkit.sh -filter_value train.xyz virial 100
```

**Output:** `filtered_<property>.xyz` containing only structures within threshold

**Use Cases:**
- Remove unphysical high-force structures
- Filter convergence issues from DFT
- Clean training datasets before NEP training

---

#### Filter by Minimum Distance (`filter_exyz_by_dist.py`)

Remove structures with atomic distances below threshold.

**Command-line:**
```bash
# Remove structures with any interatomic distance < 1.5 Å
gpumdkit.sh -filter_dist train.xyz 1.5
```

**Output:** `filtered_dist.xyz` containing only structures passing the distance check

**Use Cases:**
- Remove overlapping atoms
- Filter DFT errors or unrelaxed structures
- Quality control for training data

---

#### Filter by Box Size (`filter_exyz_by_box.py`)

Remove structures with simulation cells smaller than threshold.

**Command-line:**
```bash
# Remove structures with box dimensions < 20 Å
gpumdkit.sh -filter_box train.xyz 20
```

**Output:** `filtered_box.xyz` containing only large enough structures

**Use Cases:**
- Ensure sufficient cutoff for MD
- Remove small test structures
- Prepare data for specific simulations

---

#### Filter by Distance Range (`filter_dist_range.py`)

Extract structures with specific interatomic distances.

**Command-line:**
```bash
# Extract structures with Li-Li distance between 1.8 and 2.0 Å
gpumdkit.sh -filter_range dump.xyz Li Li 1.8 2.0
```

**Output:** `filtered_Li_Li_1.8_2.0.xyz`

**Use Cases:**
- Study specific bonding configurations
- Analyze diffusion pathways
- Sample transition states

---

### Composition Analysis (`analyze_composition.py`)

**Option 501** in interactive mode.

Analyze and export structures by chemical composition.

**Interactive Mode:**

Select option `501` from the analyzer menu:

```bash
501
```

You will see:

```
>-------------------------------------------------<
| This function calls the script in analyzer      |
| Script: analyze_composition.py                  |
| Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
>-------------------------------------------------<
Input the name of extxyz file
------------>>
```

Enter filename:
```bash
train.xyz
```

**Command-line:**
```bash
gpumdkit.sh -analyze_comp train.xyz
```

**Example Output:**
```
 Index    Compositions           N atoms      Count
 ---------------------------------------------------
 1        Li56O96Zr16La24        192          51
 2        Li48O80Zr20La24        172          23
 3        Li64O112Zr12La24       212          15
 ---------------------------------------------------
 Enter index to export (e.g., '1,2', '2-3', 'all'), or press Enter to skip:
```

**Options:**
- Enter single index: `1`
- Multiple indices: `1,2,3`
- Range: `2-5`
- All: `all`
- Skip: Press Enter

**Output:** Separate XYZ files for each selected composition

**Use Cases:**
- Organize multi-composition datasets
- Balance training data composition
- Analyze composition-property relationships

---

### Outlier Detection (`find_outliers.py`)

**Option 502** in interactive mode.

Identify and remove outlier structures based on NEP prediction errors.

**Interactive Mode:**

Select option `502`:

```bash
502
```

You will be prompted:

```
Input the threshold of RMSE to identify outliers
---------------------------------------------------
Enter energy RMSE threshold (meV/atom): 1
Enter force RMSE threshold (meV/Å): 60
Enter stress RMSE threshold (GPa): 0.03
```

**Requirements:**
- NEP training results (`*_train.out` files)
- `train.xyz` in current directory

**Output:**
- `selected.xyz` - Outlier structures
- `remained.xyz` - Good structures
- `selected_remained.png` - Visualization of outliers vs retained structures

**Workflow:**
```bash
# 1. Train NEP model
nep

# 2. Find outliers
gpumdkit.sh
# Select: 5) Analyzer → 502

# 3. Review outliers
# 4. Decide: retrain with only remained.xyz or investigate outliers
```

**Use Cases:**
- Remove problematic training structures
- Identify DFT calculation errors
- Improve NEP model quality
- Iterative data cleaning

---

### Charge Balance Check (`charge_balance_check.py`)

Verify charge neutrality in ionic systems.

**Command-line:**
```bash
gpumdkit.sh -cbc train.xyz
```

**Process:**
```
Computing compositions: 100%|██████████| 51/51 [00:03<00:00, 13.03it/s]
Checking oxidation states: 100%|██████| 1/1 [00:03<00:00,  3.67s/it]
```

**Output:**
- `balanced.xyz` - Charge-balanced structures
- `indices.txt` - Indices of balanced structures

**Use Cases:**
- Validate ionic system charges
- Filter unphysical charge states
- Prepare data for polarizable models

---

### Time Estimation

#### GPUMD Simulation (`time_consuming_gpumd.sh`)

Monitor running GPUMD simulations and estimate completion time.

**Command-line:**
```bash
gpumdkit.sh -time gpumd
```

**Example Output:**
```
 ----------------- System Information ----------------
 total frames: 1050000
 -----------------------------------------------------
 Current Frame  Speed (steps/s)   Total Time       Time Left       Estimated End
 -------------   -------------   -------------   -------------   -----------------
     13000          499.86         0h 35m 0s      0h 34m 34s    2025-12-27 18:12:04
     14000          199.93        1h 27m 31s      1h 26m 21s    2025-12-27 19:03:56
     15000          199.94        1h 27m 31s      1h 26m 16s    2025-12-27 19:03:56
```

**Use Cases:**
- Monitor long MD simulations
- Estimate when to check results
- Optimize simulation parameters

---

#### NEP Training (`time_consuming_nep.sh`)

Monitor NEP training progress and estimate completion.

**Command-line:**
```bash
gpumdkit.sh -time nep
```

**Example Output:**
```
+-----------------+-----------+-----------------+---------------------+
|       Step      | Time Diff |    Time Left    |    Finish Time      |
+-----------------+-----------+-----------------+---------------------+
| 6700            | 1 s       | 0 h 15 m 33 s   | 2025-10-23 15:34:11 |
| 6800            | 2 s       | 0 h 31 m 4 s    | 2025-10-23 15:49:44 |
| 6900            | 2 s       | 0 h 31 m 2 s    | 2025-10-23 15:49:44 |
| 7000            | 2 s       | 0 h 31 m 0 s    | 2025-10-23 15:49:44 |
| 7100            | 3 s       | 0 h 46 m 27 s   | 2025-10-23 16:05:14 |
+-----------------+-----------+-----------------+---------------------+
```

**Use Cases:**
- Monitor NEP training progress
- Plan when to check convergence
- Estimate computational cost

---

## Common Workflows

### Pre-Training Data Quality Check

```bash
# 1. Check property ranges
gpumdkit.sh -range train.xyz energy
gpumdkit.sh -range train.xyz force
gpumdkit.sh -range train.xyz virial

# 2. Check atomic distances
gpumdkit.sh -min_dist_pbc train.xyz

# 3. Check composition
gpumdkit.sh -analyze_comp train.xyz

# 4. Filter outliers
gpumdkit.sh -filter_value train.xyz force 30
gpumdkit.sh -filter_dist filtered_force.xyz 1.5
```

### Post-Training Outlier Removal

```bash
# 1. Train NEP
nep

# 2. Identify outliers
gpumdkit.sh
# Select: 5) Analyzer → 502
# Set thresholds

# 3. Retrain with clean data
mv remained.xyz train.xyz
nep
```

### Dataset Organization

```bash
# 1. Analyze composition
gpumdkit.sh -analyze_comp combined_data.xyz

# 2. Export each composition
# (Follow interactive prompts)

# 3. Balance training set
# Use equal numbers from each composition
```

### Monitor Long Calculations

```bash
# During GPUMD simulation
gpumdkit.sh -time gpumd

# During NEP training
gpumdkit.sh -time nep
```

---

## Quality Control Guidelines

### Energy Thresholds

**Typical ranges (per atom):**
- Cohesive energy: -10 to 0 eV
- Formation energy: -2 to 2 eV
- High-temperature MD: wider range acceptable

**Action:** Filter structures with unusual energies

### Force Thresholds

**Typical ranges:**
- Relaxed structures: < 0.05 eV/Å
- MD sampling: < 10 eV/Å
- High-temperature: < 30 eV/Å

**Action:** Remove very high-force structures (> 30 eV/Å)

### Distance Thresholds

**Element-dependent:**
- Li-Li: > 1.5 Å
- O-O: > 2.0 Å
- Metal-metal: > 2.0 Å

**Action:** Filter structures with unphysical overlaps

### Box Size Recommendations

**Minimum dimensions:**
- Short-range interactions: > 10 Å
- Medium-range: > 15 Å
- Long-range: > 20 Å

**Action:** Ensure box large enough for cutoff radius

---

## Tips and Best Practices

💡 **Always use PBC**: Use `-min_dist_pbc` for periodic systems to avoid incorrect distances

💡 **Check before training**: Run all quality checks before starting expensive NEP training

💡 **Document thresholds**: Keep records of filtering parameters for reproducibility

💡 **Visualize distributions**: Use histogram option to understand data distribution

💡 **Iterative cleaning**: Multiple rounds of filtering may be necessary

💡 **Save originals**: Keep unfiltered data as backup before applying filters

💡 **Composition balance**: Ensure balanced representation of all compositions in training

💡 **Monitor progress**: Use time estimation tools for long calculations

---

## Troubleshooting

### Issue: Minimum distances seem wrong

**Solution:** Use `-min_dist_pbc` instead of `-min_dist` for periodic systems

### Issue: Too many structures filtered out

**Solution:** 
- Check if thresholds are too strict
- Verify DFT calculations converged properly
- Consider if structures are physically meaningful

### Issue: Outlier detection removes too many structures

**Solution:**
- Relax RMSE thresholds
- Check if NEP model is undertrained
- Verify training data quality

### Issue: Time estimation jumps around

**Solution:**
- Normal for adaptive time steps
- Wait for more data points
- Average over multiple readings

---

For more details, see [Scripts/analyzer/README.md](../../Scripts/analyzer/README.md)
