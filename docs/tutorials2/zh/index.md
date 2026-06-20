# GPUMDkit 教程

<div align="center">
  <p>
    <strong>中文</strong> | <a href="../en/index.md">English</a>
  </p>
</div>

欢迎使用 **GPUMDkit** - 一个强大的 GPUMD 和 NEP 命令行工具包！

## 什么是 GPUMDkit？

GPUMDkit 可以简化您的分子动力学工作流程：

- **格式转换**：在 VASP、LAMMPS、CP2K、ABACUS、CIF 和 extxyz 格式之间转换
- **计算器**：计算离子电导率、描述符、MSD、NEB 等
- **分析工具**：结构验证、过滤、成分分析
- **可视化**：31+ 种绘图类型，用于训练、传输和结构分析
- **工作流**：DFT 和 MD 模拟的批处理
- **结构采样**：使用均匀、随机或 FPS 方法选择结构

## 快速开始

### 安装

```bash
# 克隆仓库
git clone https://github.com/zhyan0603/GPUMDkit.git

# 安装
cd GPUMDkit
source ./install.sh
```

### 依赖

```bash
# 创建 conda 环境
conda create -n gpumdkit python=3.12
conda activate gpumdkit

# 安装必需的包
pip install neptrain ase pymatgen dpdata numpy scipy matplotlib
```

### 基本用法

```bash
# 交互模式（菜单驱动）
gpumdkit.sh

# 命令行模式
gpumdkit.sh -<选项> [参数...]

# 获取帮助
gpumdkit.sh -h                    # 通用帮助
gpumdkit.sh -<选项> -h            # 特定命令帮助
gpumdkit.sh -plt -h               # 绘图帮助
gpumdkit.sh -calc -h              # 计算器帮助
```

## 教程列表

| 教程 | 描述 |
|------|------|
| [快速入门](quickstart.md) | 安装、配置和第一步 |
| [格式转换](format_conversion.md) | 在结构文件格式之间转换 |
| [计算器](calculators.md) | 计算材料属性 |
| [分析工具](analyzers.md) | 结构分析和验证 |
| [可视化](visualization.md) | 绘图和数据可视化 |
| [工作流](workflows.md) | 批处理和自动化 |
| [结构采样](sampling.md) | 结构选择方法 |
| [NEP 训练指南](nep_training.md) | 完整的 NEP 模型训练工作流 |

## 常见工作流程

### 工作流 1：NEP 模型训练

```bash
# 1. 转换 DFT 数据
gpumdkit.sh -out2xyz ./vasp_results/

# 2. 添加组标签
gpumdkit.sh -addgroup POSCAR Li Y Cl

# 3. 采样结构
gpumdkit.sh  # 选择: 2) Sample Structures -> 203

# 4. 训练 NEP（外部）
# ...

# 5. 验证
gpumdkit.sh -plt train
gpumdkit.sh -plt prediction
```

### 工作流 2：离子电导率计算

```bash
# 1. 计算 MSD
gpumdkit.sh -calc msd trajectory.xyz Li 10

# 2. 计算电导率
gpumdkit.sh -calc ionic-cond Li 1

# 3. 可视化
gpumdkit.sh -plt msd
gpumdkit.sh -plt sdc
```

### 工作流 3：结构分析

```bash
# 1. 检查成分
gpumdkit.sh -analyze_comp train.xyz

# 2. 检查距离
gpumdkit.sh -min_dist_pbc train.xyz

# 3. 检查力范围
gpumdkit.sh -range train.xyz force
```

## 模块概览

| 模块 | 菜单 | CLI 标志 | 描述 |
|------|------|----------|------|
| 格式转换 | 1 | `-out2xyz`, `-pos2exyz` 等 | 结构文件转换 |
| 结构采样 | 2 | `-pynep` | 结构采样 |
| 工作流 | 3 | - | 批处理 |
| 计算器 | 4 | `-calc <类型>` | 属性计算 |
| 分析工具 | 5 | `-range`, `-min_dist` 等 | 结构分析 |
| 可视化 | 6 | `-plt <类型>` | 绘图工具 |
| 实用工具 | 7 | `-time`, `-clean` | 实用功能 |

## 输出文件约定

| 文件 | 内容 |
|------|------|
| `model.xyz` | extxyz 格式的结构文件 |
| `train.xyz` | 训练数据集 |
| `nep.txt` | NEP 模型文件 |
| `msd.out` | 均方位移数据 |
| `thermo.out` | 热力学属性 |
| `loss.out` | 训练损失历史 |
| `rdf.out` | 径向分布函数 |
| `sdc.out` | 自扩散系数数据 |

## 获取帮助

- **通用帮助**：`gpumdkit.sh -h`
- **绘图帮助**：`gpumdkit.sh -plt -h`
- **计算器帮助**：`gpumdkit.sh -calc -h`
- **特定命令**：`gpumdkit.sh -<选项> -h`

## 联系方式

- **开发者**：严子涵 (yanzihan@westlake.edu.cn)
- **GitHub**：https://github.com/zhyan0603/GPUMDkit
- **文档**：https://zhyan0603.github.io/GPUMDkit/
