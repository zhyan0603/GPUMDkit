# 🔄 Format Conversion

This guide helps you quickly choose and run GPUMDkit format-conversion tools.

**Script directory:** `Scripts/format_conversion/`  
**Interactive menu:** `gpumdkit.sh` → `1) Format Conversion`

---

## Quick Start

### Interactive mode

```bash
gpumdkit.sh
# choose: 1
```

### Command-line mode (recommended for reproducibility)

```bash
gpumdkit.sh -h
```

---

## Most Common Recipes

### VASP/OUTCAR → extxyz

```bash
gpumdkit.sh -out2xyz .
# or
gpumdkit.sh -out2exyz .
```

### POSCAR → extxyz

```bash
gpumdkit.sh -pos2exyz POSCAR model.xyz
```

### extxyz → POSCAR series

```bash
gpumdkit.sh -exyz2pos train.xyz
```

### CIF → extxyz / POSCAR

```bash
gpumdkit.sh -cif2exyz input.cif model.xyz
gpumdkit.sh -cif2pos input.cif POSCAR.vasp
```

### LAMMPS dump → extxyz

```bash
gpumdkit.sh -lmp2exyz dump.lammpstrj Li Y Cl
```

### Extract frames from extxyz

```bash
gpumdkit.sh -get_frame dump.xyz 1000
gpumdkit.sh -frame_range dump.xyz 0.2 0.5
```

---

## Command Reference

| Task | Command |
|---|---|
| OUTCAR → extxyz (shell) | `gpumdkit.sh -out2xyz <dir>` |
| OUTCAR → extxyz (python) | `gpumdkit.sh -out2exyz <dir>` |
| CP2K log → xyz | `gpumdkit.sh -cp2k2xyz <log_file>` |
| XDATCAR → extxyz | `gpumdkit.sh -xdat2exyz <xdatcar>` |
| POSCAR → extxyz | `gpumdkit.sh -pos2exyz <poscar> <xyz>` |
| extxyz → POSCAR | `gpumdkit.sh -exyz2pos <xyz>` |
| POSCAR → LAMMPS data | `gpumdkit.sh -pos2lmp <poscar> <lmp_data> <elem...>` |
| LAMMPS dump → extxyz | `gpumdkit.sh -lmp2exyz <dump> <elem...>` |
| CIF → POSCAR | `gpumdkit.sh -cif2pos <cif> <poscar>` |
| CIF → extxyz | `gpumdkit.sh -cif2exyz <cif> <xyz>` |
| Add group labels | `gpumdkit.sh -addgroup <poscar> <elem...>` |
| Add structure weights | `gpumdkit.sh -addweight <in.xyz> <out.xyz> <weight>` |
| Replicate structure | `gpumdkit.sh -replicate <in> <out> 2 2 2` |
| Replicate to target atom count | `gpumdkit.sh -replicate <in> <out> <target_num>` |
| Clean extra extxyz tags | `gpumdkit.sh -clean_xyz <in.xyz> <out.xyz>` |
| Extract one frame | `gpumdkit.sh -get_frame <xyz> <index>` |
| Extract frame range | `gpumdkit.sh -frame_range <xyz> <start_frac> <end_frac>` |

---

## Interactive Menu Mapping

Menu items map to frequently used scripts:

- `101` → `out2xyz.sh`
- `102` → `mtp2xyz.py`
- `103` → `cp2k2xyz.py`
- `104` → `abacus2xyz-*.sh`
- `105` → `exyz2pos.py`
- `106` → `add_groups.py`
- `107` → `add_weight.py`
- `108` → `get_frame.py`
- `109` → clean extxyz info
- `110` → replicate structure

If unsure which tool to use, start from interactive mode first.

---

## Output Validation Checklist

After conversion, always check:

1. **Atom count** is unchanged (unless intentionally modified).
2. **Element order** is correct (especially LAMMPS dump conversions).
3. **Cell vectors / PBC** are reasonable.
4. **Energy/force fields** exist when expected.

You can use analyzer commands for quick checks:

```bash
gpumdkit.sh -analyze_comp model.xyz
gpumdkit.sh -range model.xyz force
```

---

## Common Pitfalls

1. **Wrong element order in `-lmp2exyz`**
   - The element list must match dump type ordering.
2. **Missing output filename**
   - Some converters require explicit output paths.
3. **Mixed-quality datasets**
   - Run analyzer tools after conversion before model training.

---

## Related Guides

- [Tutorials Overview](tutorials.md)
- [Analyzer](analyzer.md)
- [Sample Structures](sample_structures.md)
- [Workflow](workflow.md)

