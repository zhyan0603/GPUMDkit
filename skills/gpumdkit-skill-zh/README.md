# GPUMDkit 智能体技能

GPUMDkit 提供两个便携式的 Agent Skills 兼容入口点（英文和中文）：

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

单一的 `gpumdkit-skill`（及其对应的中文版本 `gpumdkit-skill-zh`）处理 GPUMDkit 命令、GPUMD 模拟、NEP 训练、验证、后处理和 GPUMDkit 开发。其 `SKILL.md` 将智能体路由到 `references/` 下的专注文件，因此详细的模块内容仅在需要时加载。

GPUMD 和 NEP 参考文档是自包含的摘要，按输入、设置/动作、系综、计算、输出、训练数据和模型参数组织。它们不指向或需要本地 GPUMD 文档源码树。当可执行文件与捆绑快照不同时，智能体记录版本和解析器证据，并询问版本特定信息。

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

使用相同的 `ln -s` 命令和适当的目标目录。支持自定义技能路径的客户端可以直接指向 `${GPUMDkit_path}/skills/gpumdkit-skill`。

## 便携性说明

- `SKILL.md` 使用标准 YAML frontmatter，包含目录匹配的 `name` 和触发器聚焦的 `description`。
- 核心指令保持简洁；模块详情使用一级 `references/` 进行渐进式披露。
- 该技能不依赖实验性的 `allowed-tools` 元数据或产品特定的激活函数。
- `agents/openai.yaml` 是可选的 UI 元数据。不识别它的客户端可以忽略它。
- 命令和资源使用仓库相对路径或 `${GPUMDkit_path}`。
- 任何未解决的科学或执行选择必须询问用户，而不是猜测。

## 参考资料

- [Agent Skills 规范](https://agentskills.io/specification)
- [Agent Skills 客户端实现指南](https://agentskills.io/client-implementation/adding-skills-support)
- [GPUMDkit 文档](https://gpumdkit.cn/)
