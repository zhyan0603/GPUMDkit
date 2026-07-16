---
name: gpumdkit-skill-zh
description: 端到端使用 GPUMDkit、GPUMD 与 NEP。适用于安装和了解 GPUMDkit，转换、采样、划分、筛选、分析、计算或绘制原子尺度数据集，准备批处理或主动学习流程，创建 model.xyz、run.in 或 nep.in，设置、运行、验证和后处理 GPUMD 模拟，训练、微调或评估 NEP 模型，开展扩散、离子电导率、Arrhenius、热输运等研究，以及维护 GPUMDkit 的代码、文档和技能文件。
---

# GPUMDkit

将本技能作为 GPUMDkit 工具、GPUMD 模拟、NEP 训练与评估、结果验证和后处理的统一入口。每次只加载当前任务真正需要的参考文件。

## 必须遵守的规则

- 在给出命令、修改输入、运行工具或解释结果之前，先阅读相应的参考文件。本地参考或帮助信息可用时，不要凭记忆作答。
- 只要尚未确定的选择可能影响科学含义、数值稳定性、数据筛选、单位、计算成本、文件内容或结果解释，就必须询问用户。不得自行编造“看起来合理”的参数。
- 未获得用户输入或未找到明确的项目既有设置时，不得擅自决定势函数、元素/类型映射、电荷、温度、时间步长、系综、压力、约束、运行步数、采样间隔、拟合区间、截断半径、随机种子或收敛标准。
- 明确区分文档默认值、项目约定、一般建议和用户决定。不得把默认值悄悄当成正式计算参数。
- 未经用户明确授权，不得覆盖输入文件、丢弃结构、启动 GPUMD/NEP/DFT、提交作业或开始耗时、昂贵的计算。
- 遇到解析错误、NaN/Inf、输出缺失、模拟不稳定或无法解释的警告时，应立即停止。先报告证据，再询问是否调整科学参数或重试。
- 保留原始数据和工作区中与当前任务无关的改动。记录数据转换、筛选条件、排除项、执行命令、软件版本和所作假设。

## 参考文件路由

根据任务读取一个或多个参考文件，不要默认加载全部内容。

| 任务 | 必读参考 |
|---|---|
| 安装、查找命令、选择模块、常见 CLI 用法 | [overview.md](references/overview.md) |
| 结构/数据格式转换、标签、权重、扩胞、帧提取 | [format-conversion.md](references/format-conversion.md) |
| 均匀/随机/FPS 采样、训练集/测试集划分、结构扰动、力偏差筛选 | [sampling.md](references/sampling.md) |
| 批量 SCF/MD 准备和主动学习流程 | [workflows.md](references/workflows.md) |
| MSD、离子电导率、描述符、NEB、结构优化、极化 | [calculators.md](references/calculators.md) |
| 组成、数值范围、原子间距、结构筛选、异常值、概率密度 | [analyzers.md](references/analyzers.md) |
| NEP/MD/输运/结构绘图与拟合输出 | [visualization.md](references/visualization.md) |
| 新增或调整 GPUMDkit 图形；绘图风格与审查清单 | [plotting-style.md](references/plotting-style.md) |
| 任意 GPUMD 模拟任务；模拟流程和内置参考路由 | [gpumd.md](references/gpumd.md) |
| GPUMD `model.xyz`、`run.in`、分组、单位和辅助输入 | [gpumd-inputs.md](references/gpumd-inputs.md) |
| GPUMD 势函数与初始化、模拟盒控制、受力、约束、结构优化 | [gpumd-setup.md](references/gpumd-setup.md) |
| GPUMD NVE/NVT/NPT、MTTK、QTB、热输运、PIMD、TI、冲击系综 | [gpumd-ensembles.md](references/gpumd-ensembles.md) |
| GPUMD 计算关键字、参数格式、输出和兼容性 | [gpumd-computes.md](references/gpumd-computes.md) |
| GPUMD 轨迹输出、`active`/`observer` 模式和写出行为 | [gpumd-outputs.md](references/gpumd-outputs.md) |
| 任意 NEP 训练或预测任务；训练流程和内置参考路由 | [nep.md](references/nep.md) |
| NEP `train.xyz`/`test.xyz` 数据格式、目标量、单位和审查 | [nep-data.md](references/nep-data.md) |
| 当前 `nep.in` 参数、默认值、取值范围和模式 | [nep-parameters.md](references/nep-parameters.md) |
| NEP 损失、模型、断点续训、parity 图和验证 | [nep-outputs.md](references/nep-outputs.md) |
| 扩散或电导率的变温序列与 Arrhenius 活化能 | [arrhenius.md](references/arrhenius.md) |
| GPUMDkit 代码、CLI、脚本、文档、测试或技能维护 | [contributing.md](references/contributing.md) |

跨模块任务必须加载所有相关参考。例如：

- Arrhenius 研究：读取 `gpumd.md`、`gpumd-inputs.md`、`gpumd-ensembles.md`、`gpumd-computes.md`、`gpumd-outputs.md`、`arrhenius.md`、`calculators.md` 和 `visualization.md`；如需准备 `model.xyz`，再读取 `format-conversion.md`。
- NEP 训练流程：读取 `nep.md`、`nep-data.md`、`nep-parameters.md`、`nep-outputs.md`、`format-conversion.md`、`analyzers.md`、`sampling.md` 和 `visualization.md`。
- 面向主动学习的 GPUMD 批量采样：读取 `gpumd.md`、采样方案涉及的 GPUMD 参数参考、`workflows.md`、`sampling.md` 和 `analyzers.md`。
- 新增 GPUMDkit 命令：读取 `contributing.md` 和受影响模块的参考文件。
- 新增或调整绘图：读取 `visualization.md`、`plotting-style.md`、`contributing.md`，以及该物理量对应的科学参考。

## 信息来源优先级

按以下顺序采用证据：

1. 本技能内置且自成体系的 GPUMD 与 NEP 参考。
2. 本地可执行文件的帮助、解析器报错和版本输出。
3. `gpumdkit.sh -h`、模块帮助和目标 Python 脚本的 `-h` 输出。
4. 已确认可在同一可执行文件版本下正常工作的现有项目输入。
5. 当本地行为与技能内置快照不一致时，由用户提供的特定版本文档。

若不同来源互相矛盾，应列出具体冲突、可执行文件版本和解析器证据，再询问本任务应以哪个版本为准。不要猜测特定版本的语法，也不要把旧的独立 GPUMD/NEP 技能当作权威来源。

## 工作流程

1. 判断任务类型，并加载满足任务所需的最小参考集合。
2. 检查现有文件、可执行文件版本、本地帮助和项目约定。
3. 用简洁问题列出仍缺少的决定。可以继续安全的只读检查，但不得自行补上科学参数。
4. 在执行昂贵计算前，说明已经确认的文件方案、命令、预期输出和验证标准。
5. 有文档支持时，优先使用 `gpumdkit.sh` 命令。只有参考文件明确指出该功能仅能从菜单/调试路径调用，或确实没有 CLI 入口时，才直接运行 `${GPUMDkit_path}/Scripts/` 下的脚本。
6. 获得执行授权后，记录完整命令、工作目录、版本、退出状态和警告；合适时先运行经用户同意的小规模测试。
7. 开始后续分析前，先验证文件结构和物理行为。将平衡阶段与正式采样阶段分开，并说明排除的数据。
8. 交付结果时给出单位、数据来源、收敛或拟合证据、不确定性和适用限制。

## 可移植性

- 所有链接均相对于本技能目录解析，不依赖特定产品的激活工具。
- 使用标准 Markdown、YAML frontmatter、Shell 命令和仓库相对资源，使其可被兼容 Agent Skills 的客户端读取。
- `agents/openai.yaml` 只是可选的界面元数据，不得依赖其中内容决定行为。
- 优先通过 `${GPUMDkit_path}` 访问仓库资源；若变量未设置，应从当前工作区定位仓库或询问用户，不要猜测路径。
