<div align="center">
  <h1>🛠️ Troubleshooting</h1>
  <p style="text-align: justify;">This page collects the most common problems users hit when installing or running GPUMDkit, with a quick diagnostic command and fix for each. If your issue is not listed here, see <a href="#getting-more-help">Getting More Help</a> at the bottom.</p>
</div>

## Common Issues

### 1. After install, `gpumdkit.sh: command not found` in a new terminal

**Symptom:** Opening a fresh terminal and typing `gpumdkit.sh` gives `command not found`, even though installation appeared to succeed.

**Diagnostic:**

```bash
echo $GPUMDkit_path
```

This should print the install directory. If it is empty or wrong, your shell did not load the environment variables.

**Fix:**

```bash
source ~/.bashrc        # on Linux / bash
source ~/.zshrc         # on macOS / zsh
```

If the variable is still empty after sourcing, run the installer again from the GPUMDkit directory:

```bash
cd /path/to/GPUMDkit
bash install.sh
```

---

### 2. Plot crashes with `FileNotFoundError: loss.out` (or `thermo.out`, `msd.out`, etc.)

**Symptom:** Running `gpumdkit.sh -plt <type>` fails with a `FileNotFoundError` for an `.out` file.

**Diagnostic:**

```bash
ls *.out
```

If nothing is listed, you are not in the directory that holds the GPUMD/NEP output files.

**Fix:**

```bash
cd /path/to/your/simulation/run
gpumdkit.sh -plt <type>
```

Always run `-plt` from the directory containing the output files, or pass the correct path.

---

### 3. `No module named 'calorine'` (or `'neptrain'`, `'dpdata'`) when running certain calculators

**Symptom:** A calculator or converter aborts with `ModuleNotFoundError: No module named 'calorine'` (or `neptrain`, `dpdata`).

**Diagnostic:**

```bash
pip list | grep -iE 'calorine|neptrain|dpdata'
```

If the package is missing, install it.

**Fix:**

```bash
pip install calorine neptrain dpdata
```

Or, if you use conda, run the equivalent `conda install` / `pip install` inside your activated `gpumdkit` environment. Make sure the environment is activated before running GPUMDkit.

---

### 4. MSD plot oscillates wildly instead of showing linear diffusion

**Symptom:** The MSD curve is noisy or oscillates instead of growing roughly linearly with time, giving an unreliable diffusion coefficient.

**Diagnostic:**

Check the actual frame interval used in the simulation:

```bash
head thermo.out
```

Compare the sampling interval there with the `dt_fs` value you passed to `-calc msd`.

**Fix:**

- Ensure `dt_fs` passed to `gpumdkit.sh -calc msd` equals the real MD timestep multiplied by the sample (output) interval — not just the MD timestep alone.
- Confirm the trajectory is long enough for the diffusion timescale of your system.
- Note that solids at low temperature may genuinely not show a linear MSD regime (the atoms are localized); this is physical, not a bug.

---

### 5. `-out2xyz` or `-pos2exyz` says "No files found"

**Symptom:** A format-conversion command reports that no input files were found, even though you expected them to be there.

**Diagnostic:**

```bash
ls *.vasp *.xyz *.OUTCAR
```

If nothing matches, you are either in the wrong directory or the file extensions do not match what the command expects.

**Fix:**

- `cd` into the directory that actually contains your structure files before running the conversion.
- Verify the file extensions match the expected ones (e.g., `POSCAR`, `*.vasp`, `OUTCAR`, `*.xyz`).

---

## Getting More Help

If none of the above solves your problem:

- **Issue tracker:** Open an issue at <https://github.com/zhyan0603/GPUMDkit/issues> (please include the full command you ran and the error message).
- **Email:** Zihan YAN — <yanzihan@westlake.edu.cn>
- **QQ group:** 825696376 — for community discussion and quick help.
