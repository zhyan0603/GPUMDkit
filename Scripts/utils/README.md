<div align="center">
  <h1>🛠️ Utility Scripts</h1>
  <p style="text-align: justify;">Utility scripts provide maintenance helpers for GPUMDkit itself and small data-cleanup tools used by command-line workflows.</p>
</div>

## Overview

| Tool | Entry | Purpose |
|------|-------|---------|
| `clean_extra_files.sh` | `gpumdkit.sh -clean` | Interactively remove temporary files in the current directory |
| `update_gpumdkit.sh` | `gpumdkit.sh -update` or `gpumdkit.sh -U` | Pull the latest GPUMDkit code from the current Git branch |
| `completion.sh` | sourced by the installer | Bash completion for `gpumdkit.sh` commands |
| `renumber_atoms.py` | `gpumdkit.sh -re_atoms <input> <output>` | Renumber atom IDs in LAMMPS dump files |
| `nep_modifier/` | `gpumdkit.sh -nep_modifier` | Interactive NEP model editing utilities |

---

## Clean Extra Files

```bash
gpumdkit.sh -clean
```

This helper scans the current directory and proposes files to remove. It keeps common GPUMDkit input files such as `run.in`, `nep.in`, `model.xyz`, `nep.txt`, `train.xyz`, and `test.xyz`, plus shell/slurm submission files.

Before deleting anything, it prints the candidate file list and asks for confirmation:

```text
Do you want to delete all these files?
y/yes to delete, n/no to cancel, or input files to keep (separated by spaces):
```

Use this command only after checking that the current directory is the calculation folder you intend to clean.

---

## Update GPUMDkit

```bash
gpumdkit.sh -update
gpumdkit.sh -U
```

The updater checks the current Git branch against the remote GPUMDkit repository and runs `git pull` when updates are available. If the command is launched outside a Git repository, it tries to use the configured `GPUMDkit_path`.

---

## Shell Completion

`completion.sh` provides Bash completion for common GPUMDkit options, plot types, and calculator types. It is normally configured by the installer, so users do not need to run it manually.

---

## Renumber LAMMPS Atom IDs

```bash
gpumdkit.sh -re_atoms dump.lammpstrj renumbered.lammpstrj
```

This command rewrites atom IDs in each `ITEM: ATOMS` section of a LAMMPS dump file so that IDs start from 1 in every frame.

---

## NEP Modifier

```bash
gpumdkit.sh -nep_modifier
```

The NEP modifier is an interactive utility for editing NEP model files. See:

- `Scripts/utils/nep_modifier/README.md`
- `Scripts/utils/nep_modifier/README_zh-CN.md`
