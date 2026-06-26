<div align="center">
  <h1>GPUMDkit 教程</h1>
  <p>
    <strong>简体中文</strong> | <a href="../en/README.md">English</a>
  </p>
</div>

欢迎使用 **GPUMDkit** - GPUMD 和 NEP 的命令行工具包。

## 什么是 GPUMDkit？

GPUMDkit 简化了计算材料科学中的常见任务：

- **格式转换**：在 VASP、LAMMPS、CP2K、ABACUS、CIF 和 extxyz 格式之间转换
- **可视化**：用于 NEP 训练、MD 模拟和分析的绘图工具
- **计算器**：计算离子电导率、描述符、MSD、NEB 等
- **分析工具**：结构验证、过滤、成分分析
- **工作流**：DFT 和 MD 模拟的批处理
- **结构采样**：使用均匀、随机或 FPS 方法选择结构

## 快速开始

### 安装

```bash
git clone https://github.com/zhyan0603/GPUMDkit.git
cd GPUMDkit
source ./install.sh
```

### 依赖

```bash
conda create -n gpumdkit python=3.12
conda activate gpumdkit
pip install neptrain ase pymatgen dpdata numpy scipy matplotlib
```

### 使用方法

```bash
# 交互模式
gpumdkit.sh

# 命令行模式
gpumdkit.sh -<选项> [参数...]

# 获取帮助
gpumdkit.sh -h
```

## 教程列表

| 教程 | 描述 |
|------|------|
| [快速入门](quickstart.md) | 安装和第一步 |
| [格式转换](format_conversion.md) | 在文件格式之间转换 |
| [计算器](calculators.md) | 计算材料属性 |
| [分析工具](analyzers.md) | 结构分析和验证 |
| [可视化](visualization.md) | 绘图和可视化 |
| [工作流](workflows.md) | 批处理和自动化 |
| [结构采样](sampling.md) | 结构选择方法 |
| [NEP 训练](nep_training.md) | 完整的 NEP 训练工作流 |

## 示例：NEP 训练流程

```bash
# 转换 DFT 数据
gpumdkit.sh -out2xyz ./vasp_results/

# 添加组标签
gpumdkit.sh -addgroup POSCAR Li Y Cl

# 采样结构
gpumdkit.sh  # 选择: 2) Sample Structures -> 203

# 验证训练
gpumdkit.sh -plt train
```

## 示例：离子电导率

```bash
gpumdkit.sh -calc msd trajectory.xyz Li 10
gpumdkit.sh -calc ionic-cond Li 1
gpumdkit.sh -plt msd
```

## 模块概览

| 模块 | 菜单 | CLI | 描述 |
|------|------|-----|------|
| 格式转换 | 1 | `-out2xyz`, `-pos2exyz` | 文件转换 |
| 结构采样 | 2 | `-pynep` | 结构采样 |
| 工作流 | 3 | - | 批处理 |
| 计算器 | 4 | `-calc <类型>` | 属性计算 |
| 分析工具 | 5 | `-range`, `-min_dist` | 结构分析 |
| 可视化 | 6 | `-plt <类型>` | 绘图工具 |
| 实用工具 | 7 | `-time`, `-clean` | 实用功能 |

## 链接

- GitHub: https://github.com/zhyan0603/GPUMDkit
- 文档: https://zhyan0603.github.io/GPUMDkit/
