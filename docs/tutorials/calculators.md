# Function 4 - Calculators

This section covers the calculator tools in GPUMDkit (Interactive Mode - Function 4), which provide computational tools for calculating various material properties using NEP models and analyzing simulation data.

## Interactive Mode Access

```bash
gpumdkit.sh
# Select: 4) Calculators
```

You'll see the following menu:

```
 ------------>>
 401) Calc ionic conductivity
 402) Calc properties by nep
 403) Calc descriptors
 404) Calc DOAS
 405) Calc NEB
 000) Return to the main menu
 ------------>>
```

## Command-Line Usage

For quick operations, you can also use command-line mode:

```bash
# Ionic conductivity
gpumdkit.sh -calc ionic-cond <element> <charge>

# NEP predictions
gpumdkit.sh -calc nep <input.xyz> <output.xyz> <nep.txt>

# Descriptors
gpumdkit.sh -calc des <method> <input.xyz> <output.npy> <nep.txt> <element>

# DOAS calculation
gpumdkit.sh -calc doas <input.xyz> <nep.txt> <output.txt>
```

---

## Available Calculators

### Ionic Conductivity (`calc_ion_conductivity.py`)

**Option 401** in interactive mode.

Calculates ionic diffusivity and conductivity from MSD data.

**Interactive Mode:**

Select option `401` from the calculator menu:

```bash
401
```

You will see the following prompt:

```
 >-------------------------------------------------<
 | This function calls the script in calculators   |
 | Script: calc_ion_conductivity.py                |
 | Developer: Zihan YAN (yanzihan@westlake.edu.cn) |
 >-------------------------------------------------<
 Input <element> <charge> (eg. Li 1)
 ------------>>
```

Enter the element and charge:

```bash
Li 1
```

**Command-line:**
```bash
gpumdkit.sh -calc ionic-cond Li 1
```

**Input files:**
- `msd.out` (required)
- `thermo.out` (optional - for temperature)
- `model.xyz` (optional - for volume/atom count)

**Output:**
```
Diffusivity (D):
  D_x: 4.153e-07 cm^2/s
  D_y: 4.174e-07 cm^2/s
  D_z: 2.610e-07 cm^2/s
  D_total: 3.646e-07 cm^2/s
------------------------------
Ionic Conductivity:
  Sigma_x: 2.576e-02 mS/cm
  Sigma_y: 2.589e-02 mS/cm
  Sigma_z: 1.619e-02 mS/cm
  Sigma_total: 2.261e-02 mS/cm
```

### NEP Property Calculation (`calc_properties_with_nep.py`)

Calculates energies, forces, and stresses using a NEP model.

**Command-line:**
```bash
gpumdkit.sh -calc nep input.xyz output.xyz nep.txt
```

**Interactive:** Select option `402`

**Use cases:**
- Generate predictions for test sets
- Calculate properties for new structures
- Validate NEP model performance

### Descriptor Calculation (`calc_descriptors.py`)

Calculates NEP descriptors for structure analysis and visualization.

**Command-line:**
```bash
gpumdkit.sh -calc des umap train.xyz descriptors.npy nep.txt Li
```

**Interactive:** Select option `403`

**Methods:**
- `umap` - UMAP dimensionality reduction
- `tsne` - t-SNE dimensionality reduction
- `pca` - PCA dimensionality reduction

**Workflow:**
```bash
# 1. Calculate descriptors
gpumdkit.sh -calc des umap train.xyz desc.npy nep.txt Li

# 2. Visualize
gpumdkit.sh -plt des umap desc.npy
```

### DOAS Calculation (`calc_doas.py`)

Calculates density of atomistic states by optimizing structures with NEP.

**Command-line:**
```bash
gpumdkit.sh -calc doas structures.xyz nep.txt doas_output.txt
```

**Interactive:** Select option `404`

**Workflow:**
```bash
# 1. Calculate DOAS
gpumdkit.sh -calc doas structures.xyz nep.txt doas.txt

# 2. Visualize
gpumdkit.sh -plt doas doas.txt
```

Reference: [Wang et al.](https://doi.org/10.1002/anie.202215544)

### NEB Calculation (`neb_calculation.py`)

Performs nudged elastic band calculations for transition state finding.

**Direct execution only:**
```bash
python ${GPUMDkit_path}/Scripts/calculators/neb_calculation.py initial.xyz final.xyz 9 nep.txt
```

**Interactive:** Select option `405`

**Output:**
- `neb.traj` - Trajectory file with all images
- `neb_pathway.png` - Energy vs reaction coordinate plot
- Migration barrier in terminal output

### RDF Calculator (`rdf_calculator_ovito.py`)

Calculates radial distribution function using OVITO.

**Direct execution:**
```bash
python ${GPUMDkit_path}/Scripts/calculators/rdf_calculator_ovito.py trajectory.xyz 6.0 400
```

**Parameters:**
- Input XYZ file
- Cutoff distance (Angstroms)
- Number of bins

**Output:** `rdf.txt` file

## Interactive Mode

Access calculators through the interactive menu:

```bash
gpumdkit.sh
# Select: 4) Calculators
# Choose option 401-405
```

### Menu Options

```
 ------------>>
 401) Calc ionic conductivity
 402) Calc properties by nep
 403) Calc descriptors
 404) Calc DOAS
 405) Calc NEB
 000) Return to the main menu
 ------------>>
```

## Common Workflows

### Analyze Ionic Transport
```bash
# 1. Run GPUMD with compute_msd
# 2. Calculate conductivity
gpumdkit.sh -calc ionic-cond Li 1

# 3. Plot diffusion
gpumdkit.sh -plt msd
gpumdkit.sh -plt sdc
```

### NEP Model Analysis
```bash
# 1. Test NEP on structures
gpumdkit.sh -calc nep test.xyz predictions.xyz nep.txt

# 2. Calculate descriptors
gpumdkit.sh -calc des umap train.xyz desc.npy nep.txt Li

# 3. Visualize
gpumdkit.sh -plt prediction
gpumdkit.sh -plt des umap desc.npy
```

### Find Migration Barrier
```bash
# 1. Prepare initial and final structures
# 2. Run NEB
python ${GPUMDkit_path}/Scripts/calculators/neb_calculation.py start.xyz end.xyz 9 nep.txt

# 3. Examine neb_pathway.png for barrier height
```

## Tips and Best Practices

### Ionic Conductivity Calculations

**Requirements:**
- ✅ Long enough MD simulation (>500 ps recommended)
- ✅ Sufficient statistics (multiple independent runs preferred)
- ✅ Equilibrated system before MSD calculation

**Best practices:**
```bash
# 1. Run long simulation
cat > run.in << EOF
velocity 600
ensemble npt_ber 600 600 100 0 0 0 1000
dump_thermo 1000
compute_msd 100 1  # Start after equilibration
run 500000         # 500 ps at 1 fs timestep
EOF

# 2. Check MSD linearity
gpumdkit.sh -plt msd
# Should be linear in diffusive regime

# 3. Check convergence
gpumdkit.sh -plt msd_conv
# If available, verify convergence

# 4. Calculate conductivity
gpumdkit.sh -calc ionic-cond Li 1

# 5. Compare with experiment
# For LGPS at 300K: ~10-12 mS/cm expected
```

**Troubleshooting:**
- If conductivity seems too high → Check temperature, may not be equilibrated
- If MSD not linear → Simulation too short or not diffusive
- If D_x, D_y, D_z very different → Check if anisotropic (normal) or simulation issue

---

### Descriptor Calculations

**Purpose:** Understand training set coverage in descriptor space

**Workflow:**
```bash
# 1. Calculate descriptors for all data
gpumdkit.sh -calc des umap train.xyz train_desc.npy nep.txt Li
gpumdkit.sh -calc des umap test.xyz test_desc.npy nep.txt Li
gpumdkit.sh -calc des umap active.xyz active_desc.npy nep.txt Li

# 2. Visualize overlaid
# (Requires custom plotting script)
python << 'EOF'
import numpy as np
import matplotlib.pyplot as plt

train = np.load('train_desc.npy')
test = np.load('test_desc.npy')
active = np.load('active_desc.npy')

plt.figure(figsize=(10, 8))
plt.scatter(train[:, 0], train[:, 1], alpha=0.3, s=20, label='Training', c='blue')
plt.scatter(test[:, 0], test[:, 1], alpha=0.5, s=30, label='Test', c='green')
plt.scatter(active[:, 0], active[:, 1], alpha=0.7, s=50, label='Active', c='red')
plt.xlabel('Descriptor 1')
plt.ylabel('Descriptor 2')
plt.legend()
plt.title('Descriptor Space Coverage')
plt.savefig('descriptor_coverage.png', dpi=300)
EOF
```

**Interpretation:**
- ✅ Good: Test points within training cloud
- ✅ Good: Active points at training boundaries (exploring new regions)
- ⚠️ Warning: Test points far from training (extrapolation)
- ⚠️ Warning: Large gaps in training cloud (need more sampling)

**Use for training set curation:**
```bash
# If descriptor analysis shows gaps:
# 1. Identify undersampled regions
# 2. Run targeted MD to sample those regions
# 3. Add to training set
# 4. Recalculate descriptors to verify coverage
```

---

### DOAS Calculations

**Density of Atomistic States** - Energy distribution analysis

**Setup:**
```bash
# 1. Prepare diverse structures (different compositions/phases)
cat phase_A/*.xyz phase_B/*.xyz > structures.xyz

# 2. Calculate DOAS (will optimize each structure first)
gpumdkit.sh -calc doas structures.xyz nep.txt doas.out

# 3. Visualize
gpumdkit.sh -plt doas doas.out Li
```

**Advanced: Compare with DFT**
```bash
# Calculate DOAS with DFT for reference
# Then compare NEP vs DFT DOAS distributions

python << 'EOF'
import numpy as np
import matplotlib.pyplot as plt

doas_nep = np.loadtxt('doas_nep.out')
doas_dft = np.loadtxt('doas_dft.out')

plt.figure(figsize=(12, 5))

plt.subplot(121)
plt.hist(doas_nep, bins=50, alpha=0.7, label='NEP')
plt.hist(doas_dft, bins=50, alpha=0.7, label='DFT')
plt.xlabel('Energy (eV)')
plt.ylabel('Frequency')
plt.legend()
plt.title('DOAS Comparison')

plt.subplot(122)
plt.hist2d(doas_dft, doas_nep, bins=50, cmap='viridis')
plt.xlabel('DFT Energy (eV)')
plt.ylabel('NEP Energy (eV)')
plt.title('Correlation')
plt.colorbar(label='Count')

plt.tight_layout()
plt.savefig('doas_comparison.png', dpi=300)
EOF
```

---

### NEB Calculations

**Nudged Elastic Band** - Transition state finding

**Preparing Structures:**
```bash
# 1. Start with relaxed initial and final states
gpumd  # Or DFT relaxation
cp dump_last_frame.xyz initial.xyz

# 2. Prepare final state (manual or from trajectory)
# Move atoms to final position
cp target_structure.xyz final.xyz

# 3. Verify endpoints are relaxed
# Check forces are small
```

**Parameter Selection:**

| System Size | N_images | When to Use |
|-------------|----------|-------------|
| < 50 atoms | 11-15 | Small systems, simple barriers |
| 50-200 atoms | 7-11 | Medium systems, typical use |
| > 200 atoms | 5-7 | Large systems, computational cost |

**Running NEB:**
```bash
# Standard NEB
python ${GPUMDkit_path}/Scripts/calculators/neb_calculation.py \
    initial.xyz final.xyz 9 nep.txt

# Outputs:
# - neb.traj: Full NEB trajectory
# - neb_pathway.png: Energy vs reaction coordinate
# - Terminal: Migration barrier value
```

**Interpreting Results:**
```
# Example output:
Migration barrier: 0.45 eV

# Check neb_pathway.png:
# ✓ Smooth curve: Good convergence
# ✓ Single peak: Simple barrier
# ✗ Multiple peaks: Complex pathway, may need more images
# ✗ Jagged: Poor convergence, adjust parameters
```

**Advanced: Temperature-dependent barriers**
```bash
# Calculate at multiple temperatures
for temp in 300 600 900; do
    # Thermalize structures
    cat > run.in << EOF
velocity $temp
ensemble nvt_ber $temp $temp 100
run 10000
EOF
    
    gpumd
    cp dump_last.xyz initial_${temp}K.xyz
    
    # Run NEB
    python neb_calculation.py initial_${temp}K.xyz final_${temp}K.xyz 9 nep.txt
    mv neb_pathway.png neb_${temp}K.png
done

# Compare barriers at different temperatures
```

---

## Dependencies and Installation

### Required Packages

**Core (always needed):**
```bash
pip install numpy ase
```

**For specific calculators:**

| Calculator | Required Packages | Install Command |
|-----------|-------------------|-----------------|
| `calc_ion_conductivity.py` | numpy | `pip install numpy` |
| `calc_properties_with_nep.py` | calorine | `pip install calorine` |
| `calc_descriptors.py` | pynep, numpy | `pip install pynep numpy` |
| `calc_doas.py` | calorine, ase | `pip install calorine ase` |
| `neb_calculation.py` | ase, calorine | `pip install ase calorine` |
| `rdf_calculator_ovito.py` | ovito | `pip install ovito` |

**Complete installation:**
```bash
# For all calculator functionality
pip install numpy ase calorine pynep

# Optional (for RDF calculator)
pip install ovito
```

### Verifying Installation

```bash
# Test each calculator
python << 'EOF'
try:
    import numpy
    print("✓ NumPy")
except:
    print("✗ NumPy - install with: pip install numpy")

try:
    import ase
    print("✓ ASE")
except:
    print("✗ ASE - install with: pip install ase")

try:
    import calorine
    print("✓ Calorine")
except:
    print("✗ Calorine - install with: pip install calorine")

try:
    from pynep.calculate import NEP
    print("✓ PyNEP")
except:
    print("✗ PyNEP - install with: pip install pynep")

try:
    from ovito.io import import_file
    print("✓ OVITO")
except:
    print("✗ OVITO - install with: pip install ovito")
EOF
```

---

## Troubleshooting

### Issue: Ionic Conductivity Gives Zero or Very Small Values

**Symptoms:**
```
Diffusivity (D):
  D_total: 1.234e-15 cm^2/s  ← Too small!
```

**Causes & Solutions:**

1. **Simulation too short**
   ```bash
   # Check MSD plot
   gpumdkit.sh -plt msd
   # If not linear → run longer
   
   # Solution: Increase run steps
   run 1000000  # Instead of 100000
   ```

2. **Temperature too low (no diffusion)**
   ```bash
   # Solution: Increase temperature
   velocity 800  # Instead of 300
   ```

3. **Wrong ion species or charge**
   ```bash
   # Solution: Verify species
   gpumdkit.sh -calc ionic-cond Li 1  # Correct charge
   ```

---

### Issue: NEP Property Calculation Fails

**Symptoms:**
```
Error: Could not load NEP model
```

**Solutions:**

1. **Calorine not installed**
   ```bash
   pip install calorine
   ```

2. **NEP model file corrupted**
   ```bash
   # Verify nep.txt is complete
   tail -5 nep.txt
   # Should show "potential_name"
   ```

3. **Structure file format issue**
   ```bash
   # Verify extxyz format
   head -3 input.xyz
   # Check lattice and properties are correct
   ```

---

### Issue: Descriptor Calculation Memory Error

**Symptoms:**
```
MemoryError: Unable to allocate array
```

**Solutions:**

1. **Process in chunks**
   ```python
   # Instead of processing all at once
   for i in range(0, n_structures, chunk_size):
       chunk = structures[i:i+chunk_size]
       descriptors_chunk = calculate_descriptors(chunk)
       np.save(f'desc_chunk_{i}.npy', descriptors_chunk)
   ```

2. **Use fewer structures**
   ```bash
   # Sample before descriptor calculation
   python sample_structures.py train.xyz random 1000
   gpumdkit.sh -calc des umap sampled_structures.xyz desc.npy nep.txt Li
   ```

---

### Issue: DOAS Takes Very Long

**Symptoms:**
```bash
# Optimization of structures is slow
```

**Solutions:**

1. **Use GPUMD minimize instead**
   ```bash
   # In GPUMD:
   minimize sd 1e-4     # Steepest descent
   minimize cg 1e-5     # Conjugate gradient
   dump_exyz            # Output optimized structures
   ```

2. **Pre-filter structures**
   ```bash
   # Remove high-energy outliers before DOAS
   gpumdkit.sh -filter_value structures.xyz energy 0.5
   ```

3. **Parallelize**
   ```bash
   # Split and run in parallel
   split -n 4 structures.xyz chunk_
   for f in chunk_*; do
       gpumdkit.sh -calc doas $f nep.txt doas_$f.out &
   done
   wait
   cat doas_*.out > doas_all.out
   ```

---

### Issue: NEB Doesn't Converge

**Symptoms:**
```
# neb_pathway.png shows jagged curve
# Or very high barrier (> 5 eV)
```

**Solutions:**

1. **Check endpoints are relaxed**
   ```bash
   # Relax initial and final states first
   # Forces should be < 0.05 eV/Å
   ```

2. **Increase number of images**
   ```bash
   # Was: 7 images
   # Try: 11 or 15 images
   ```

3. **Adjust NEB parameters**
   ```python
   # In neb_calculation.py, modify:
   neb = NEB(..., climb=True, k=0.1)  # Adjust spring constant
   opt = FIRE(neb, maxstep=0.2)       # Adjust max step
   ```

---

For more details, see [Scripts/calculators/README.md](../../Scripts/calculators/README.md)

