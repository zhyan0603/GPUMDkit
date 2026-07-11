# GPUMDkit Agent Skill

GPUMDkit provides two portable Agent Skills-compatible entry points (English and Chinese):

```text
skills/gpumdkit-skill/     (English)
  SKILL.md
  references/
  agents/openai.yaml

skills/gpumdkit-skill-zh/  (Chinese)
  SKILL.md
  references/
  agents/openai.yaml
```

The single `gpumdkit-skill` (and its Chinese counterpart `gpumdkit-skill-zh`) handles GPUMDkit commands, GPUMD simulation, NEP training, validation, post-processing, and GPUMDkit development. Its `SKILL.md` routes agents to focused files under `references/`, so detailed module content is loaded only when required.

GPUMD and NEP references are self-contained summaries organized by inputs, setup/actions, ensembles, computations, outputs, training data, and model parameters. They do not point to or require a local GPUMD documentation source tree. When an executable differs from the bundled snapshot, agents record the version and parser evidence and ask for version-specific information.

## Discover the skill

After installing GPUMDkit and setting `GPUMDkit_path`, run:

```bash
gpumdkit.sh -skill
```

Canonical skill directory:

```bash
${GPUMDkit_path}/skills/gpumdkit-skill
${GPUMDkit_path}/skills/gpumdkit-skill-zh
```

## Install for Agent Skills-compatible clients

The cross-client convention is `.agents/skills/`. Ask the user whether installation should be global or project-local before creating a link.

Global:

```bash
target_dir="${HOME}/.agents/skills"
mkdir -p "${target_dir}"
ln -s "${GPUMDkit_path}/skills/gpumdkit-skill" "${target_dir}/gpumdkit-skill"
ln -s "${GPUMDkit_path}/skills/gpumdkit-skill-zh" "${target_dir}/gpumdkit-skill-zh"
```

Current project only:

```bash
target_dir=".agents/skills"
mkdir -p "${target_dir}"
ln -s "${GPUMDkit_path}/skills/gpumdkit-skill" "${target_dir}/gpumdkit-skill"
ln -s "${GPUMDkit_path}/skills/gpumdkit-skill-zh" "${target_dir}/gpumdkit-skill-zh"
```

Common client-specific alternatives:

| Client mode | Global directory | Project directory |
|---|---|---|
| Agent Skills standard | `~/.agents/skills` | `.agents/skills` |
| OpenCode | `~/.config/opencode/skills` | `.opencode/skills` |
| Claude-compatible | `~/.claude/skills` | `.claude/skills` |

Use the same `ln -s` command with the appropriate target directory. Clients that support a custom skill path can point directly to `${GPUMDkit_path}/skills/gpumdkit-skill`.

## Portability notes

- `SKILL.md` uses standard YAML frontmatter with a directory-matching `name` and a trigger-focused `description`.
- Core instructions remain concise; module details use one-level `references/` for progressive disclosure.
- The skill does not rely on experimental `allowed-tools` metadata or a product-specific activation function.
- `agents/openai.yaml` is optional UI metadata. Clients that do not recognize it can ignore it.
- Commands and resources use repository-relative paths or `${GPUMDkit_path}`.
- Any unresolved scientific or execution choice must be asked of the user rather than guessed.

## References

- [Agent Skills specification](https://agentskills.io/specification)
- [Agent Skills client implementation guide](https://agentskills.io/client-implementation/adding-skills-support)
- [GPUMDkit documentation](https://gpumdkit.cn/)
