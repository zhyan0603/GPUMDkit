# GPUMDkit Skills

This directory contains GPUMDkit skill definitions for AI agent tools such as opencode and Claude Code.

After installing GPUMDkit, make sure `GPUMDkit_path` is set, then run:

```bash
gpumdkit.sh -skill
```

The command prints the skill source path, available GPUMDkit skills, global and project install targets, and copyable install examples for common agent tools.

## Skill Source

The canonical skill source is:

```bash
${GPUMDkit_path}/skills
```

Use this absolute environment-variable path when linking skills into an agent tool. This is more robust than relative paths such as `../../skills/...` because it does not depend on the current working directory.

## Available Skills

| Skill | Use When |
|-------|----------|
| [gpumdkit-main](gpumdkit-main/SKILL.md) | Installation, entry point, overview, and module navigation |
| [gpumdkit-format-conversion](gpumdkit-format-conversion/SKILL.md) | Converting structure files between VASP, LAMMPS, CP2K, ABACUS, CIF, DeepMD, MTP, ASE, and extxyz formats |
| [gpumdkit-calculators](gpumdkit-calculators/SKILL.md) | Computing material properties such as MSD, ionic conductivity, NEP predictions, descriptors, DOAS, NEB, and polarization |
| [gpumdkit-analyzers](gpumdkit-analyzers/SKILL.md) | Analyzing, validating, filtering, and checking structures or datasets |
| [gpumdkit-visualization](gpumdkit-visualization/SKILL.md) | Plotting NEP training, MD outputs, transport properties, RDF, thermal conductivity, and descriptors |
| [gpumdkit-workflows](gpumdkit-workflows/SKILL.md) | Preparing batch SCF/MD calculations and active-learning workflows |
| [gpumdkit-sampling](gpumdkit-sampling/SKILL.md) | Sampling, perturbing, and selecting structures for training data |
| [gpumdkit-contributing](gpumdkit-contributing/SKILL.md) | Adding features, commands, scripts, or documentation to GPUMDkit |

## Install Examples

Agents should ask the user whether to install GPUMDkit skills globally or only for the current project.

Global install targets:

| Tool mode | Target directory |
|-----------|------------------|
| OpenCode | `~/.config/opencode/skills` |
| Claude-compatible | `~/.claude/skills` |
| Agents-compatible | `~/.agents/skills` |

Project install targets:

| Tool mode | Target directory |
|-----------|------------------|
| OpenCode | `.opencode/skills` |
| Claude-compatible | `.claude/skills` |
| Agents-compatible | `.agents/skills` |

Install globally for OpenCode:

```bash
target_dir="${HOME}/.config/opencode/skills"
mkdir -p "${target_dir}"
for skill_dir in ${GPUMDkit_path}/skills/gpumdkit-*; do
    ln -s "${skill_dir}" "${target_dir}/$(basename "${skill_dir}")"
done
```

Install for the current project only:

```bash
target_dir=".opencode/skills"
mkdir -p "${target_dir}"
for skill_dir in ${GPUMDkit_path}/skills/gpumdkit-*; do
    ln -s "${skill_dir}" "${target_dir}/$(basename "${skill_dir}")"
done
```

Use one of the other target directories above for Claude-compatible or agents-compatible installs.

## Agent Usage Notes

- Prefer `gpumdkit.sh -skill` as the discovery command.
- Prefer specialized GPUMDkit skills for specific tasks instead of `gpumdkit-main`.
- When modifying code, read [gpumdkit-contributing](gpumdkit-contributing/SKILL.md) first.
- For direct script calls, use `${GPUMDkit_path}/Scripts/...` so the command does not depend on the current directory.

## References

- [Agent Skills Standard](https://agentskills.io)
- [opencode Documentation](https://opencode.ai)
- [Claude Code Skills](https://docs.anthropic.com/en/docs/claude-code/skills)
- [GPUMDkit Documentation](https://gpumdkit.cn/)
