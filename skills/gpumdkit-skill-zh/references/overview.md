# GPUMDkit 概览与命令查询

## 目录

- 快速开始
- 模块概览
- 常用工作流
- 输出文件
- 实用命令
- 故障排除与文档

GPUMDkit 是面向 GPUMD（Graphics Processing Units Molecular Dynamics）和
NEP（Neuroevolution Potential）的命令行工具包，用于简化计算材料科学中的常见数据处理、模拟和分析任务。

## 快速开始

### 安装
```bash
conda create -n gpumdkit -c gpumdkit -c conda-forge gpumdkit
conda activate gpumdkit
```

部分功能需要额外安装可选依赖：

```bash
pip install neptrain calorine
```

如需从源码安装：

```bash
git clone https://github.com/zhyan0603/GPUMDkit.git
cd GPUMDkit && source ./install.sh
```

### 基本用法
```bash
# 交互模式（菜单驱动）
gpumdkit.sh

# 命令行模式
gpumdkit.sh -<option> [args...]

# 获取帮助
gpumdkit.sh -h                    # 通用帮助
gpumdkit.sh -<option> -h          # Python 支持的命令帮助（已实现时）
gpumdkit.sh -plt -h               # 绘图帮助
gpumdkit.sh -calc -h              # 计算器帮助
```

对于由 Python 脚本实现的命令，详细用法、类型转换、文件检查和错误信息通常由目标脚本负责，可通过 `-h` 查看。`gpumdkit.sh` 主要转发参数；`-plt`、`-calc` 等大型分发入口提供模块级帮助。`-pynep`、`-time` 等旧式 Shell/菜单入口不一定支持同样的逐命令帮助方式，应查阅模块参考或实际交互提示，不要默认它们支持 `-h`。

## 模块概览

| 模块 | 描述 | 参考 |
|--------|-------------|-----------|
| 格式转换 | 在 VASP、LAMMPS、CP2K、ABACUS、CIF、extxyz 之间转换 | `references/format-conversion.md` |
| 计算器 | 计算离子电导率、描述符、MSD、NEB 等 | `references/calculators.md` |
| 分析器 | 结构验证、过滤、成分分析 | `references/analyzers.md` |
| 可视化 | 绘制训练结果、输运性质、结构数据 | `references/visualization.md` |
| 工作流 | DFT 和 MD 模拟的批处理 | `references/workflows.md` |
| 采样 | 使用均匀、随机或 FPS 方法选取结构并划分训练集/测试集 | `references/sampling.md` |
| GPUMD/NEP | 规划模拟、训练、验证和后处理 | `references/gpumd.md`、`references/nep.md` |

## 实用命令

| 命令 | 行为与安全边界 |
|---|---|
| `gpumdkit.sh -doctor` | 检查配置路径、Python/Bash 版本以及常用或特定功能所需的 Python 软件包；缺少可选包不影响无关功能 |
| `gpumdkit.sh -skill` | 打印规范的 Skill 路径和跨客户端安装提示 |
| `gpumdkit.sh -time <gpumd\|nep>` | 旧版耗时分析器；仅使用支持的 `gpumd` 或 `nep` 选择器 |
| `gpumdkit.sh -nep_modifier [nep.txt] [nep.restart|-] [nep.in|-]` | 通过 calorine 检查和修改 NEP4 模型；保留源文件，扩展、缩减和添加元素需要匹配的 restart 数据及用户明确的科学选择 |
| `gpumdkit.sh -clean` | 从当前目录中删除生成的/多余的文件；运行前预览清理实现并获得明确删除批准 |
| `gpumdkit.sh -update` | 运行 GPUMDkit 的网络自更新；先检查工作区变更并获得明确更新授权 |

自定义命令教程中显示的命令（如 `-greet`、`-batch_plot` 或 `-prep_training`）是扩展 GPUMDkit 的示例，除非用户已安装这些自定义项，否则不是内置命令。

## 常用工作流

### 工作流 1：NEP 模型训练流程
```bash
# 1. 将 DFT 数据转换为 extxyz
gpumdkit.sh -out2xyz ./vasp_results/

# 2. 仅在下游工作流使用 GPUMD 分组时添加分组标签
gpumdkit.sh -addgroup POSCAR Li Y Cl

# 3. 采样多样化结构
gpumdkit.sh  # 选择：2) Sample Structures -> 203) FPS by NepTrain

# 4. 训练 NEP 模型（外部）
# ... 训练 nep.txt ...

# 5. 验证训练
gpumdkit.sh -plt train
gpumdkit.sh -plt prediction

# 6. 检查数据质量
gpumdkit.sh -range train.xyz force
gpumdkit.sh -min_dist_pbc train.xyz
```

### 工作流 2：离子电导率计算
```bash
# 路径 A：从存储的 extxyz 轨迹计算 MSD
gpumdkit.sh -calc msd trajectory.xyz Li 10

# 路径 B：使用 GPUMD compute_msd 直接生成 msd.out
# 除非比较实现，否则不要同时运行两条路径。

# 从验证过的 msd.out/model.xyz/run.in 输入计算离子电导率
gpumdkit.sh -calc ionic-cond Li 1

# 绘制结果
gpumdkit.sh -plt msd
gpumdkit.sh -plt sdc

# Arrhenius 分析（多温度）
gpumdkit.sh -plt sigma
gpumdkit.sh -plt D
```

### 工作流 3：结构分析流程
```bash
# 1. 检查成分
gpumdkit.sh -analyze_comp train.xyz

# 2. 检查最小距离
gpumdkit.sh -min_dist_pbc train.xyz

# 3. 过滤结构
gpumdkit.sh  # 选择：5) Analyzer -> 506) Filter by distance

# 4. 检查能量/力范围
gpumdkit.sh -range train.xyz force
gpumdkit.sh -range train.xyz energy
```

## 输出文件约定

| 文件 | 内容 |
|------|---------|
| `model.xyz` | extxyz 格式的结构文件 |
| `train.xyz` | 训练数据集 |
| `nep.txt` | NEP 模型文件 |
| `msd.out` | 均方位移数据 |
| `thermo.out` | 热力学性质 |
| `loss.out` | 训练损失历史 |
| `rdf.out` | 径向分布函数 |
| `sdc.out` | 自扩散系数数据 |

## 故障排除

**问题**：命令未找到
**解决方案**：确保 GPUMDkit 已安装并在 PATH 中：`source ./install.sh`

**问题**：缺少 Python 包
**解决方案**：激活 conda 环境：`conda activate gpumdkit`

**问题**：脚本错误
**解决方案**：使用 `-h` 标志检查输入文件格式和必需参数

## 详细文档

需要面向用户的教程详情时，请阅读 `${GPUMDkit_path}/docs/tutorials/en/` 或 `/zh/` 下的对应页面。
