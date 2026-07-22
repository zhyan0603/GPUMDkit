# GPUMDkit Agent Skill

GPUMDkit 提供英文和中文两个可移植、兼容 Agent Skills 规范的入口：

```text
skills/gpumdkit-skill/     (English)
  SKILL.md
  references/
  agents/openai.yaml

skills/gpumdkit-skill-zh/  (中文)
  SKILL.md
  references/
  agents/openai.yaml
```

`gpumdkit-skill` 是英文技术基准，`gpumdkit-skill-zh` 是逐文件对应并经过中文润色的版本。两者都覆盖 GPUMDkit 命令、GPUMD 模拟、NEP 训练与验证、结果后处理和 GPUMDkit 开发。`SKILL.md` 根据任务把 Agent 路由到 `references/` 中的专题文件，避免一次加载全部细节。

GPUMD 和 NEP 参考内容自成体系，分别按输入、初始化与操作、系综、计算、输出、训练数据和模型参数组织，不依赖本地 GPUMD 文档源码。当本地可执行文件与技能内置的版本快照不一致时，Agent 应记录版本和解析器证据，并向用户确认应采用的版本规则。

## 发现技能

安装 GPUMDkit 并设置 `GPUMDkit_path` 后，运行：

```bash
gpumdkit.sh -skill
```

规范技能目录：

```bash
${GPUMDkit_path}/skills/gpumdkit-skill
${GPUMDkit_path}/skills/gpumdkit-skill-zh
```

## 为 Agent Skills 兼容客户端安装

跨客户端约定是 `.agents/skills/`。在创建链接之前，询问用户安装应该是全局的还是项目本地的。

全局安装：

```bash
target_dir="${HOME}/.agents/skills"
mkdir -p "${target_dir}"
ln -s "${GPUMDkit_path}/skills/gpumdkit-skill" "${target_dir}/gpumdkit-skill"
ln -s "${GPUMDkit_path}/skills/gpumdkit-skill-zh" "${target_dir}/gpumdkit-skill-zh"
```

仅当前项目：

```bash
target_dir=".agents/skills"
mkdir -p "${target_dir}"
ln -s "${GPUMDkit_path}/skills/gpumdkit-skill" "${target_dir}/gpumdkit-skill"
ln -s "${GPUMDkit_path}/skills/gpumdkit-skill-zh" "${target_dir}/gpumdkit-skill-zh"
```

常见的客户端特定替代方案：

| 客户端模式 | 全局目录 | 项目目录 |
|---|---|---|
| Agent Skills 标准 | `~/.agents/skills` | `.agents/skills` |
| OpenCode | `~/.config/opencode/skills` | `.opencode/skills` |
| Claude 兼容 | `~/.claude/skills` | `.claude/skills` |

使用相同的 `ln -s` 命令和适当的目标目录。支持自定义技能路径的客户端可以直接指向任一语言目录。

## 便携性说明

- `SKILL.md` 使用标准 YAML frontmatter；英文名称为 `gpumdkit-skill`，中文名称为 `gpumdkit-skill-zh`，避免安装时发生名称冲突。
- 核心指令保持简洁；模块详情使用一级 `references/` 进行渐进式披露。
- 该技能不依赖实验性的 `allowed-tools` 元数据或产品特定的激活函数。
- `agents/openai.yaml` 是可选的 UI 元数据。不识别它的客户端可以忽略它。
- 命令和资源使用仓库相对路径或 `${GPUMDkit_path}`。
- 任何未解决的科学或执行选择必须询问用户，而不是猜测。

## 参考资料

- [Agent Skills 规范](https://agentskills.io/specification)
- [Agent Skills 客户端实现指南](https://agentskills.io/client-implementation/adding-skills-support)
- [GPUMDkit 文档](https://gpumdkit.cn/)
