# Arrhenius 扩散与离子电导率工作流

当用户询问 Arrhenius 任务、活化能、温度依赖的扩散系数或离子电导率时，使用本参考。本工作流将 GPUMD 模拟设置与 GPUMDkit 分析结合；它不授权选择科学参数或启动模拟。

## 目录

- 明确目标
- 规划温度系列
- 准备 GPUMD 输入
- 生成并验证 MSD
- 运行 GPUMDkit 分析
- 解读与报告
- 当前工具限制

## 明确目标

请用户确认：

- 扩散系数 `D`、离子电导率 `sigma`，或两者都要；
- 移动物种、分组定义和形式/有效电荷；
- 各向同性总量还是方向性结果；
- 温度列表以及是否预期单一 Arrhenius 机制；
- 结构、势函数、压力/体积条件和相稳定性；
- 系综序列、时间步长、平衡时长、生产时长和重复次数；
- MSD 采样间隔/相关长度或轨迹帧间隔；
- 拟合接受标准、不确定度方法，以及是否允许外推（例如到 300 K）。

不要从文件夹名、化学式或常见做法中推断这些信息。如果用户只说"做 Arrhenius 计算"，在创建生产输入之前先问这些问题。

## 规划温度系列

每个整数温度使用一个独立目录，因为当前绘图脚本会发现匹配 `<integer>K` 的名称，例如：

```text
project/
  500K/
    model.xyz
    nep.txt
    run.in
  600K/
  700K/
```

除非研究设计另有说明，否则在各温度之间使用相同的组成和受控的、用户批准的协议。记录种子和重复次数。两个有效温度是直线的数学最低要求，但科学证据较弱；请用户确认充分的设计，而不是自己选择数量。

不要假设所有温度属于同一线性机制。在拟合之前检查相变、机制变化、非线性和不稳定的轨迹。

## 准备 GPUMD 输入

对于直接 GPUMD MSD，当只有一种物种/子集是移动的时，在 `model.xyz` 中添加分组映射，然后使用经过验证的命令，例如：

```text
compute_msd <sample_interval> <Nc> group <group_method> <group_id> [save_every <interval>]
```

确切的系综和运行时长必须来自用户批准的设计。在正确的生产运行中重新发出非传播的 compute/dump 命令。

对于基于轨迹的 GPUMDkit MSD，输出带有足够元数据的 extxyz 帧用于周期性展开。将存储帧之间的时间传递给 `-calc msd`；按以下公式计算：

```text
frame_dt_fs = integration_time_step_fs * dump_interval_steps
```

从 `run.in` 和实际轨迹节奏确认此值。当帧的转储频率较低时，不要传递积分时间步长。

## 生成并验证 MSD

两条支持的路径是：

### GPUMD 原生 MSD

使用 `compute_msd` 运行 GPUMD。原生 `msd.out` 包含时间、三个 MSD 列和三个 SDC 列（一个分组）。GPUMDkit Arrhenius 脚本使用前四列。

### GPUMDkit 轨迹 MSD

从每个温度目录：

```bash
gpumdkit.sh -calc msd <trajectory.xyz> <element> <frame_dt_fs> [max_corr_steps]
```

这会在当前目录写入四列的 `msd.out`。

在拟合任何点之前：

- 检查轨迹和热力学稳定性；
- 验证选定的原子和周期性展开；
- 确认物理上有意义的扩散线性机制，而非弹道、笼式或饱和行为；
- 在适用时检查 `gpumdkit.sh -plt msd`、`-plt sdc` 和 `-plt msd_conv`；
- 根据用户的不确定度计划比较独立重复或块估计；
- 仅在有明确记录决定的情况下保留或拒绝可疑点。

## 运行 GPUMDkit 分析

### 扩散系数 Arrhenius 拟合

从包含 `<integer>K/msd.out` 的父目录：

```bash
gpumdkit.sh -plt arrhenius_d save
gpumdkit.sh -plt D_xyz save
```

当前脚本拟合每个 MSD 系列的中间 40%-80%，将斜率从 Angstrom^2/ps 转换为 cm^2/s，拟合 `log10(D)` 对 `1000/T`，并以 eV 报告活化能。`D_xyz` 报告总量和方向性值。

40%-80% 范围是 GPUMDkit 项目的既定选择。不要随意更改。如果它与观察到的扩散机制不匹配，请展示证据并询问用户如何处理。

### 单温度离子电导率

在包含 `msd.out`、`thermo.out`、`model.xyz` 和可选 `run.in` 的温度目录中：

```bash
gpumdkit.sh -calc ionic-cond <element> <charge>
```

这使用 Nernst-Einstein 关系，读取温度/体积，计数请求的物种，在找到 `replicate` 时考虑它，并打印方向性/总扩散系数和电导率。确认独立离子 Nernst-Einstein 电导率是预期的可观测量；除非方法被扩展，否则它不包含关联效应。

### 电导率 Arrhenius 拟合

对于满足当前脚本假设的一价 Li/Na 系统：

```bash
gpumdkit.sh -plt arrhenius_sigma save
gpumdkit.sh -plt sigma_xyz save
```

这些命令拟合 `ln(sigma*T)` 对 `1000/T` 并报告活化能。使用前请检查下方的当前限制。

## 解读与报告

至少报告：

- 包含/排除的温度及原因；
- 系综、生产时长、采样、物种/分组和重复次数；
- MSD 拟合区间（当前 GPUMDkit 脚本为 40%-80%）及其为扩散机制的证据；
- `D` 或 `sigma` 值，含单位和不确定度；
- Arrhenius 方程、变换轴、活化能、拟合质量和不确定度；
- 是否有任何值是插值或外推的；
- 电导率的 Nernst-Einstein 假设；
- 相/机制变化或非 Arrhenius 行为。

不要仅因为绘图命令完成了就从未验证的拟合中报告精确的活化能。

## 当前工具限制

- `plt_arrhenius_d.py` 和方向性变体仅发现整数 `<temperature>K` 目录，并需要可用的 `msd.out` 文件。
- Arrhenius 脚本使用固定的 40%-80% MSD 拟合区间。
- `plt_arrhenius_sigma.py` 和 `plt_arrhenius_sigma_xyz.py` 仅从第一个温度的 `model.xyz` 计数 Li 和 Na，并硬编码电荷量 `z = 1`。
- 电导率绘图从第一个温度的 `run.in` 读取复制信息，并假设整个系列的离子数相同。
- 电导率脚本外推到 300 K。在将其作为结果呈现之前，询问用户该外推在科学上是否可接受。
- 对于其他物种、多价离子、混合移动物种、变化的组成或关联电导率，不要将自动电导率 Arrhenius 脚本当作通用工具使用。使用单温度验证计算和用户批准的拟合方法，或请求代码更改。
