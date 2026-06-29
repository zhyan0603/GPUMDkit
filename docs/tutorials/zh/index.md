<div align="center">
  <h1>📚 GPUMDkit 教程</h1>
  <p style="text-align: justify;">欢迎使用 <strong>GPUMDkit</strong>，这是一个面向 GPUMD 和 NEP 的命令行工具包。</p>
</div>

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
pip install neptrain dpdata calorine
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

## 示例：转换 VASP OUTCAR

```bash
gpumdkit.sh -out2xyz ./vasp_results/
```

## 示例：绘制 NEP 结果

```bash
gpumdkit.sh -plt train
gpumdkit.sh -plt prediction
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
| 结构采样 | 2 | `-neptrain`、交互菜单 | 结构采样 |
| 工作流 | 3 | - | 批处理 |
| 计算器 | 4 | `-calc <类型>` | 属性计算 |
| 分析工具 | 5 | `-range`, `-min_dist` | 结构分析 |
| 可视化 | 6 | `-plt <类型>` | 绘图工具 |
| 实用工具 | 7 | `-time`, `-clean` | 实用功能 |

## 链接

- GitHub: https://github.com/zhyan0603/GPUMDkit
- 文档: https://zhyan0603.github.io/GPUMDkit/
