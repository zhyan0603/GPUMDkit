<div align="center">
  <h1>📚 GPUMDkit 教程</h1>
  <p style="text-align: justify;">欢迎使用 <strong>GPUMDkit</strong>，这是一个面向 GPUMD 和 NEP 的命令行工具包。</p>
</div>

## 功能简介

GPUMDkit 帮助你在计算材料科学研究中完成常见任务，无需编写自定义脚本。它支持格式转换、结构分析、性质计算、可视化和批量工作流——全部通过一个命令行工具完成。

## 常见任务

| 我想要... | 教程 |
|-----------|------|
| 安装 GPUMDkit 并运行第一个命令 | [快速入门](快速入门.md) |
| 将 VASP、LAMMPS、CP2K 或 CIF 文件转换为 extxyz | [格式转换](格式转换.md) |
| 检查结构距离、过滤数据集或查找异常值 | [分析工具](分析工具.md) |
| 计算 MSD、离子电导率或描述符 | [计算器脚本](计算器脚本.md) |
| 分析极性材料、铁电体或 ABO3 体系 | [极性材料分析](极性材料分析.md) |
| 绘制 NEP 训练结果或热力学数据 | [绘图脚本](绘图脚本.md) |
| 运行批量 DFT 或 MD 模拟 | [工作流脚本](工作流脚本.md) |
| 从数据集中选择代表性结构 | [结构采样](结构采样.md) |
| 添加自己的 GPUMDkit 快捷命令 | [自定义命令](自定义命令.md) |

## 准备工作

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
pip install neptrain dpdata calorine
```

## 交互模式

```bash
gpumdkit.sh
```

打开菜单后，按数字选择模块。

## 命令行模式

```bash
gpumdkit.sh -<选项> [参数...]
```

例如：

```bash
gpumdkit.sh -pos2exyz POSCAR model.xyz
gpumdkit.sh -plt train
gpumdkit.sh -calc msd trajectory.xyz Li 10
```

使用 `gpumdkit.sh -h` 查看所有可用选项。

## 模块概览

| 模块 | 菜单 | CLI | 描述 |
|------|------|-----|------|
| 格式转换 | 1 | `-out2xyz`, `-pos2exyz` | 文件转换 |
| 结构采样 | 2 | `-neptrain`、交互菜单 | 结构采样 |
| 工作流 | 3 | - | 批处理 |
| 计算器 | 4 | `-calc <类型>` | 属性计算 |
| 分析工具 | 5 | `-range`, `-min_dist` | 结构分析 |
| 可视化 | 6 | `-plt <类型>` | 绘图工具 |
| 实用工具 | 7 | `-time`, `-clean` | 实用功能 |

## 全部教程

| 教程 | 描述 |
|------|------|
| [快速入门](快速入门.md) | 安装和第一步 |
| [命令参考](命令参考.md) | 常用 CLI 和菜单入口速查 |
| [格式转换](格式转换.md) | 在文件格式之间转换 |
| [计算器脚本](计算器脚本.md) | 计算材料属性 |
| [分析工具](分析工具.md) | 结构分析和验证 |
| [绘图脚本](绘图脚本.md) | 绘图和可视化 |
| [工作流脚本](工作流脚本.md) | 批处理和自动化 |
| [结构采样](结构采样.md) | 结构选择方法 |
| [自定义命令](自定义命令.md) | 用户自定义 GPUMDkit 快捷命令 |
| [主动学习工作流](主动学习工作流.md) | 批量主动学习流程说明 |
| [极性材料分析](极性材料分析.md) | 铁电和极化相关工具 |
| [贡献指南](贡献指南.md) | 开发和贡献说明 |

## 链接

- GitHub: https://github.com/zhyan0603/GPUMDkit
- 文档: https://zhyan0603.github.io/GPUMDkit/
