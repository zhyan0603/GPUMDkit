---
name: gpumdkit-skill
description: 端到端使用 GPUMDkit、GPUMD 和 NEP。适用于安装或导航 GPUMDkit；转换、采样、过滤、分析、计算或绘制原子数据；准备批量或主动学习工作流；创建 model.xyz、run.in 或 nep.in；设置、运行、验证或后处理 GPUMD 模拟；训练、微调或使用 NEP 进行预测；执行扩散、离子电导率、阿伦尼乌斯、热传输或相关工作流；或修改 GPUMDkit 代码和文档。
---

# GPUMDkit

使用一个入口点处理 GPUMDkit 工具、GPUMD 模拟、NEP 训练、验证和后处理。仅加载当前任务所需的参考文档。

## 强制规则

- 在给出命令、编辑输入、运行工具或解释结果之前，阅读相关的参考文件。当本地参考或帮助命令可用时，不要依赖记忆。
- 当未解决的选择可能影响科学意义、数值稳定性、数据选择、单位、成本、文件或解释时，询问用户。永远不要编造看似合理的参数。
- 在没有用户输入或明确识别的现有项目设置的情况下，不要选择势函数、物种/类型映射、电荷、温度、时间步长、系综、压力、约束、运行长度、采样间隔、拟合窗口、截断半径、种子或收敛标准。
- 区分文档默认值、项目约定、建议和用户决定。不要默默地将默认值变成生产选择。
- 除非用户明确授权，否则不要覆盖输入、丢弃结构、启动 GPUMD/NEP/DFT、提交调度器作业或开始长时间/昂贵的计算。
- 遇到解析器错误、NaN/Inf 值、缺失输出、不稳定行为或无法解释的警告时停止。报告证据并在更改科学参数或重试之前询问。
- 保留原始数据和无关的工作树更改。记录转换、过滤器、排除项、命令、版本和假设。

## 参考路由

根据任务读取一个或多个参考文件。不要默认加载每个文件。

| 任务 | 所需参考 |
|---|---|
| 安装、命令发现、模块选择、常见 CLI 用法 | [overview.md](references/overview.md) |
| 结构/数据转换、标签、权重、复制、帧提取 | [format-conversion.md](references/format-conversion.md) |
| 均匀/随机/FPS 采样、扰动、力偏差选择 | [sampling.md](references/sampling.md) |
| 批量 SCF/MD 准备和主动学习工作流 | [workflows.md](references/workflows.md) |
| MSD、离子电导率、描述符、NEB、最小化、极化 | [calculators.md](references/calculators.md) |
| 组成、范围、距离、过滤、异常值、概率密度 | [analyzers.md](references/analyzers.md) |
| NEP/MD/传输/结构图和拟合输出 | [visualization.md](references/visualization.md) |
| 新建或调整 GPUMDkit 绘图；项目绘图风格与审查清单 | [plotting-style.md](references/plotting-style.md) |
| 任何 GPUMD 模拟任务；工作流和捆绑参考路由 | [gpumd.md](references/gpumd.md) |
| GPUMD `model.xyz`、`run.in`、组、单位、辅助输入 | [gpumd-inputs.md](references/gpumd-inputs.md) |
| GPUMD 势函数/设置、盒子控制、力、约束、最小化 | [gpumd-setup.md](references/gpumd-setup.md) |
| GPUMD NVE/NVT/NPT、MTTK、QTB、热、PIMD、TI、冲击系综 | [gpumd-ensembles.md](references/gpumd-ensembles.md) |
| GPUMD 计算关键字、签名、输出、兼容性 | [gpumd-computes.md](references/gpumd-computes.md) |
| GPUMD 转储、活动/观察者模式、输出写入行为 | [gpumd-outputs.md](references/gpumd-outputs.md) |
| 任何 NEP 训练/预测任务；工作流和捆绑参考路由 | [nep.md](references/nep.md) |
| NEP `train.xyz`/`test.xyz` 模式、目标、单位、审计 | [nep-data.md](references/nep-data.md) |
| 所有当前 `nep.in` 参数、默认值、范围和模式 | [nep-parameters.md](references/nep-parameters.md) |
| NEP 损失/模型/重启/奇偶输出和验证 | [nep-outputs.md](references/nep-outputs.md) |
| 扩散或电导率温度系列和阿伦尼乌斯活化能 | [arrhenius.md](references/arrhenius.md) |
| GPUMDkit 代码、CLI、脚本、文档、测试或技能维护 | [contributing.md](references/contributing.md) |

对于跨模块工作，加载每个相关参考。示例：

- 阿伦尼乌斯研究：读取 `gpumd.md`、`gpumd-inputs.md`、`gpumd-ensembles.md`、`gpumd-computes.md`、`gpumd-outputs.md`、`arrhenius.md`、`calculators.md` 和 `visualization.md`；如果准备 `model.xyz`，还读取 `format-conversion.md`。
- NEP 训练管道：读取 `nep.md`、`nep-data.md`、`nep-parameters.md`、`nep-outputs.md`、`format-conversion.md`、`analyzers.md`、`sampling.md` 和 `visualization.md`。
- 用于主动学习的 GPUMD 批量采样：读取 `gpumd.md`、协议使用的 GPUMD 参数参考、`workflows.md`、`sampling.md` 和 `analyzers.md`。
- 新的 GPUMDkit 命令：读取 `contributing.md` 加上受影响模块的参考。
- 新建或调整绘图：读取 `visualization.md`、`plotting-style.md` 和 `contributing.md`；还要读取所绘物理量对应的科学参考。

## 来源优先级

按此顺序使用证据：

1. 与此技能捆绑的自包含 GPUMD 和 NEP 参考。
2. 本地可执行文件帮助、解析器消息和版本输出。
3. GPUMDkit `gpumdkit.sh -h`、模块帮助和目标 Python 脚本 `-h` 输出。
4. 已确认与相同可执行文件版本一起工作的现有项目输入。
5. 当本地行为与此捆绑快照不同时，用户提供的版本特定文档。

如果来源不一致，显示确切的冲突、可执行文件版本和解析器证据，然后询问哪个软件版本控制任务。不要猜测版本特定的语法。将旧的独立 GPUMD/NEP 技能视为非权威的。

## 操作工作流

1. 分类请求并加载最小充分的参考集。
2. 检查可用文件、可执行文件版本、本地帮助和现有项目约定。
3. 将缺失的决定列为简洁的问题。继续安全的只读检查，但不要自己填补科学空白。
4. 在昂贵的执行之前，呈现解析的文件计划、命令、预期输出和验证标准。
5. 当文档记录时，优先使用 `gpumdkit.sh` CLI 命令。仅当参考识别出仅菜单/调试路径或不存在 CLI 路由时，才使用 `${GPUMDkit_path}/Scripts/` 下的直接脚本。
6. 当执行被授权时，捕获确切的命令、工作目录、版本、退出状态和警告。适当时首先使用用户批准的冒烟测试。
7. 在下游分析之前验证文件结构和物理行为。将平衡与生产分开并报告排除项。
8. 提供带有单位、来源、收敛/拟合证据、不确定性和限制的结果。

## 便携性

- 解析所有相对于此技能目录的链接。该技能不需要产品特定的激活工具。
- 使用标准 Markdown、YAML frontmatter、shell 命令和仓库相对资源，以便任何 Agent Skills 兼容客户端可以加载它。
- 将 `agents/openai.yaml` 视为可选的接口元数据；永远不要依赖它来获取指令或行为。
- 当可用时，使用 `${GPUMDkit_path}` 访问仓库资源。如果未设置，从用户或当前工作区定位仓库，而不是猜测路径。
