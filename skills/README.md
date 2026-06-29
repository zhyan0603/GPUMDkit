# GPUMDkit Skills

This directory contains skill definitions for AI agent tools (opencode, Claude Code, etc.)
to effectively use GPUMDkit functionality.

## Available Skills

| Skill | Description |
|-------|-------------|
| [gpumdkit-main](gpumdkit-main/SKILL.md) | Main entry point and overview |
| [gpumdkit-format-conversion](gpumdkit-format-conversion/SKILL.md) | Structure file format conversion |
| [gpumdkit-calculators](gpumdkit-calculators/SKILL.md) | Material property calculations |
| [gpumdkit-analyzers](gpumdkit-analyzers/SKILL.md) | Structure analysis and validation |
| [gpumdkit-visualization](gpumdkit-visualization/SKILL.md) | Plotting and visualization |
| [gpumdkit-workflows](gpumdkit-workflows/SKILL.md) | Batch processing and automation |
| [gpumdkit-sampling](gpumdkit-sampling/SKILL.md) | Structure sampling and selection |
| [gpumdkit-contributing](gpumdkit-contributing/SKILL.md) | Code conventions and adding new features |

## Usage with Different Agent Tools

### opencode (Recommended)

Configure in `opencode.json` to load skills directly from the skills directory:

```json
{
  "skills": {
    "paths": ["./skills"]
  }
}
```

Or create symlinks manually:

```bash
# Run from GPUMDkit root directory
mkdir -p .opencode/skills
cd .opencode/skills
ln -s ../../skills/gpumdkit-main .
ln -s ../../skills/gpumdkit-format-conversion .
ln -s ../../skills/gpumdkit-calculators .
ln -s ../../skills/gpumdkit-analyzers .
ln -s ../../skills/gpumdkit-visualization .
ln -s ../../skills/gpumdkit-workflows .
ln -s ../../skills/gpumdkit-sampling .
ln -s ../../skills/gpumdkit-contributing .
cd ../..
```

### Claude Code

Create symlinks from GPUMDkit root directory:

```bash
# Run from GPUMDkit root directory
mkdir -p .claude/skills
cd .claude/skills
ln -s ../../skills/gpumdkit-main .
ln -s ../../skills/gpumdkit-format-conversion .
ln -s ../../skills/gpumdkit-calculators .
ln -s ../../skills/gpumdkit-analyzers .
ln -s ../../skills/gpumdkit-visualization .
ln -s ../../skills/gpumdkit-workflows .
ln -s ../../skills/gpumdkit-sampling .
ln -s ../../skills/gpumdkit-contributing .
cd ../..
```

Or use a one-liner from GPUMDkit root:

```bash
mkdir -p .claude/skills && cd .claude/skills && for skill in ../../skills/gpumdkit-*/; do ln -s "$skill" .; done && cd ../..
```

### Other Agent Tools

Any agent tool supporting the [Agent Skills](https://agentskills.io) standard can use these skills.
Each skill directory contains:
- `SKILL.md` - Main skill file with YAML frontmatter

## Skill Design Principles

1. **Modular**: Each skill focuses on a specific functionality
2. **Self-contained**: Skills include all necessary information
3. **Cross-platform**: Works with multiple agent tools
4. **Well-documented**: Clear descriptions and examples
5. **Maintainable**: Easy to update as GPUMDkit evolves

## Skill Structure

Each skill follows the Agent Skills standard:

```markdown
---
name: skill-name
description: >
  Description of what the skill does and when to use it.
allowed-tools: Bash(gpumdkit *)
---

# Skill Title

## Content...
```

## Contributing

When adding new features to GPUMDkit:

1. Update the relevant skill file
2. Add examples and usage instructions
3. Update this README if adding new skills
4. Test with at least one agent tool

## References

- [Agent Skills Standard](https://agentskills.io)
- [opencode Documentation](https://opencode.ai)
- [Claude Code Skills](https://docs.anthropic.com/en/docs/claude-code/skills)
- [GPUMDkit Documentation](https://zhyan0603.github.io/GPUMDkit/)
