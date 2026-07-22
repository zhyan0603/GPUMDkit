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
| `doctor.py` | `gpumdkit.sh -doctor` | Report the configured path, runtime versions, and available Python packages |
| `skill_info.sh` | `gpumdkit.sh -skill` | Show the English and Chinese Agent Skill directories and installation hints |
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

## Check the environment

```bash
gpumdkit.sh -doctor
```

This read-only check reports the configured `GPUMDkit_path`, Python and Bash
versions, and common or feature-specific Python packages. A package marked
`[MISS]` is only needed for the feature named beside it; it does not prevent
unrelated commands from running. Use this before installing an optional
dependency or reporting an environment issue.

---

## Discover the Agent Skill

```bash
gpumdkit.sh -skill
```

This command only prints the available English and Chinese skill locations and
cross-client installation examples. It does not create links or change an agent
configuration. See [`skills/README.md`](../../skills/README.md) before choosing
whether a skill should be installed globally or for one project.

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
gpumdkit.sh -nep_modifier path/nep.txt path/nep.restart path/nep.in
```

The NEP modifier provides a guided interface to calorine's `augment()`,
`prune()`, `add_species()`, `remove_species()`, and `keep_species()` methods. It
can change neuron or descriptor capacity, add a charge head, extend the chemical
space, or extract a species subset. It also inspects the current architecture,
records accepted operations and seeds, updates a matching `nep.in`, and exports
a collision-free model package. Expansion, reduction, and adding species require
a matching restart file; modified models must be retrained and validated before
production use.

This implementation depends directly on calorine. Research use should cite
E. Lindgren et al., *Journal of Open Source Software* **9**(95), 6264 (2024),
<https://doi.org/10.21105/joss.06264>. Detailed usage, interaction examples, and
the complete citation are provided in:

- `Scripts/utils/nep_modifier/README.md`
- `Scripts/utils/nep_modifier/README_zh-CN.md`
