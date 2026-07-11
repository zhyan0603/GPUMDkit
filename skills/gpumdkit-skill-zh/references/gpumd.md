# GPUMD 模拟路由与工作流

对于每个 GPUMD 模拟任务，请先阅读本文件，然后加载覆盖所请求命令的参数参考。不要凭记忆构造 `model.xyz` 或 `run.in`。

## 目录

- 附带参考范围
- 任务路由
- 安全构建模拟
- 执行前后验证
- 已知版本敏感点

## 附带参考范围

本技能中的 GPUMD 文件构成了 GPUMDkit 附带的 GPUMD 输入输出规范的自包含、分类摘要。它们定义了支持的关键字族、当前签名、默认值、约束、生命周期和预期输出，无需依赖 GPUMD 源码树。

使用本地可执行文件的解析器和版本输出作为运行时证据。如果它拒绝这些参考中的语法，请停止并报告命令、输入行、确切错误和可执行文件版本。向用户询问版本特定的行为或文档；不要静默切换到记忆中的或较旧的语法。

## 任务路由

从主 `SKILL.md` 路由器直接加载以下文件：

| 需求 | 参考文件 |
|---|---|
| `model.xyz`、`run.in`、分组、单位、辅助输入 | `references/gpumd-inputs.md` |
| 势函数加载、初始化、盒子变更、力、约束、最小化、运行 | `references/gpumd-setup.md` |
| NVE/NVT/NPT、MTTK、QTB、热传导、PIMD、TI、冲击和 TTM 积分器 | `references/gpumd-ensembles.md` |
| 静态计算、MSD、RDF、输运、声子、粘度、LSQT | `references/gpumd-computes.md` |
| dump/active 关键字、输出文件名、追加/覆盖行为 | `references/gpumd-outputs.md` |
| NEP 势函数创建或预测 | `references/nep.md` 及其参数参考 |
| 扩散/电导率温度序列 | `references/arrhenius.md` |

## 安全构建模拟

1. 确认可观测量、结构、PBC、势函数文件、元素顺序、GPUMD 构建版本和交付物。
2. 加载 `gpumd-inputs.md`；验证 `model.xyz` 的原子数、晶格、属性、单位以及 `run.in` 引用的每个分组。
3. 为 `run.in` 中将出现的每个关键字加载参数参考。
4. 向用户询问所有未解决的科学选择：时间步长、温度/压力路径、系综、耦合周期、约束、阶段长度、种子、采样间隔、相关长度以及拟合/收敛准则。
5. 将最小化、平衡和生产分为明确的阶段。在需要时为每个 `run` 块重新发出非传播系综、计算、dump 和控制命令。
6. 显式计算物理时长和输出频率：

```text
duration_ps = run_steps * time_step_fs / 1000
frame_dt_fs = dump_interval_steps * time_step_fs
```

7. 在运行前列出预期输出、写入模式、估计大小和下游 GPUMDkit 命令。
8. 在用户明确授权执行之前，不要启动 GPUMD。

## 执行前后验证

执行前：

- 确保所有引用的文件存在且势函数支持每种元素；
- 保留旧的追加模式输出或使用干净目录；
- 检查盒子/PBC 与势函数和压力控制的兼容性；
- 检查相关窗口是否在生产时长范围内；
- 确认 GPU/运行时和确切的可执行文件路径。

执行后：

- 检查退出状态和日志中的解析器错误、警告、NaN/Inf 和不稳定性；
- 验证预期文件、行数、列数和输出频率；
- 将平衡阶段与生产平均值分开；
- 检查温度、能量、压力、盒子和轨迹行为；
- 在拟合或报告前验证可观测量的收敛性。

## 已知版本敏感点

- 当前 `dump_xyz` 以分组方法和分组 ID 开头；不要使用过时的 interval/frame-count 签名。
- `compute_phonon` 当前仅接受 `<displacement>` 并需要 `kpoints.in`；较旧的指南可能显示 cutoff 参数。
- HNEMD 需要温度控制；使用 Nose-Hoover 链，不要为此目的使用 Langevin。
- `compute_hac` 属于 EMD 生产运行；它不是用于先前 `compute ... jp jk` 数据的单步后处理命令。
- 许多关键字是运行作用域的，不会传播。除非其附带条目明确说明，否则为每个 `run` 块重新发出系综、计算、dump、约束、变形和驱动控制。
