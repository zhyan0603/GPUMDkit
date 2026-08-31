<p align="center">
  <img src="./Gallery/gpumdkit_logo_lateral.png" width="55%" alt="GPUMDkit Logo">
</p>
<p align="center">
  <a href="../README.md">English</a>
  &nbsp;·&nbsp;
  <strong>简体中文</strong>
  &nbsp;·&nbsp;
  <a href="https://gpumdkit.cn/">Website</a> &nbsp;·&nbsp;
  <a href="https://gpumdkit.cn/htmls/index.html">Documentation</a>
  &nbsp;·&nbsp;
  <a href="https://gpumdkit.cn/gallery.html">Gallery</a>
  &nbsp;
</p>
<p align="center">
  <a href="https://github.com/zhyan0603/GPUMDkit/releases"><img src="https://img.shields.io/github/v/tag/zhyan0603/GPUMDkit?label=version&style=flat-square&color=brightgreen" alt="Version"></a>
  <a href="https://github.com/zhyan0603/GPUMDkit/blob/main/LICENCE"><img src="https://img.shields.io/badge/license-GPL--3.0-blue" alt="License"></a>
  <a href="https://github.com/zhyan0603/GPUMDkit/stargazers"><img src="https://img.shields.io/github/stars/zhyan0603/GPUMDkit?style=social" alt="Stars"></a>
  <img src="https://img.shields.io/github/languages/code-size/zhyan0603/GPUMDkit" alt="Code Size">
  <a href="https://github.com/zhyan0603/GPUMDkit/graphs/contributors"><img src="https://img.shields.io/github/contributors/zhyan0603/GPUMDkit?style=flat-square&color=brightgreen" alt="Contributors"></a>
</p>
<p style="text-align: justify;"><strong>GPUMDkit</strong> 是面向 GPUMD（<em>Graphics Processing Units Molecular Dynamics</em>）和 NEP（<em>neuroevolution potential</em>）程序的工具包。它提供统一的命令行入口，用于调用常用脚本、完成格式转换、结构采样、NEP 数据准备、分析和可视化。</p>


## 功能特点
- **数据准备**：转换、标记、采样、划分、筛选和检查原子尺度数据集。
- **工作流自动化**：准备批量 DFT/MD 计算和主动学习工作流。
- **计算与分析**：计算和分析结构、输运及 NEP 相关性质。
- **可视化与后处理**：可视化 NEP 训练、分子动力学、扩散和热输运结果。
- **灵活的使用方式**：支持交互菜单和直接命令行选项。

## 安装

### Conda（推荐）

```bash
conda create -n gpumdkit -c gpumdkit -c conda-forge gpumdkit
conda activate gpumdkit
```

部分功能需要额外安装可选依赖：

```bash
pip install neptrain calorine
```

### 从源码安装

```bash
git clone https://github.com/zhyan0603/GPUMDkit.git
cd GPUMDkit
source ./install.sh
```

### GPUMDkit Agent Skill

GPUMDkit 包含用于 AI 协助处理 GPUMDkit、GPUMD 和 NEP 工作流的中英文
Agent Skill。安装后运行：

```bash
gpumdkit.sh -skill
```

然后让 Agent 按照命令输出的提示全局安装相关 Skill。正常使用时建议全局安装；
如果没有说明应安装到全局还是当前项目，Agent 应先询问确认，再创建链接。

## 更新

### Conda 安装

如果通过 Conda 安装了 `GPUMDkit`，请使用以下命令更新：

```bash
conda activate gpumdkit
conda update -c gpumdkit -c conda-forge gpumdkit
```

通过 pip 安装的可选依赖项可根据需要单独更新：

```bash
pip install --upgrade neptrain calorine
```

### 源码安装

如果从源码仓库安装了 `GPUMDkit`，请运行：

```bash
gpumdkit.sh -update
```

该命令会检查当前安装的 Git 分支，并从同一分支拉取最新更新。

也可以手动下载最新的源码压缩包：

```bash
wget https://github.com/zhyan0603/GPUMDkit/archive/refs/heads/main.zip
```

## 使用方法

如果需要按步骤完成首次安装，请阅读[快速入门](./tutorials/zh/快速入门.md)；
如果已经知道任务，只需查阅[命令参考](./tutorials/zh/命令参考.md)中的准确语法。

### 交互模式

运行菜单并按数字选择模块：

```bash
gpumdkit.sh
```

### 命令行模式

对于可重复执行的任务，直接使用命令行选项：

```bash
gpumdkit.sh -h
gpumdkit.sh -doctor
gpumdkit.sh -<选项> [参数...]
```

常用示例：

```bash
gpumdkit.sh -pos2exyz POSCAR model.xyz
gpumdkit.sh -plt train
gpumdkit.sh -plt thermo save
```

使用 `gpumdkit.sh -plt -h` 列出绘图类型，使用 `gpumdkit.sh -calc -h`
列出计算器命令。对于提供选项帮助的命令，可以在命令末尾加上 `-h`；各模块
所需输入、输出和科学假设请参阅教程。

#### 自定义命令

`GPUMDkit` 支持通过 `~/.gpumdkit.in` 文件自定义命令。

你可以在该文件中定义函数来添加自己的快捷命令（例如 `gpumdkit.sh -yourcommand`），从而扩展 `GPUMDkit` 的功能。详细用法请参见[自定义命令文档](https://gpumdkit.cn/)。

#### Tab 补全支持

`gpumdkit.sh` 提供了可选的 Bash `Tab` 补全功能，增强命令行使用体验。按 `Tab` 键即可自动补全主选项（如 `-h`、`-plt`、`-calc`）及其二级参数（如 `thermo`、`train`）。

##### 使用示例

- 输入 `gpumdkit.sh -<Tab>` 查看所有可用选项。
- 输入 `gpumdkit.sh -plt <Tab>` 列出绘图子选项，如 `thermo`、`train` 等。
- 输入 `gpumdkit.sh -time <Tab>` 查看计算器选项，如 `gpumd`、`nep`。

## 加入我们

欢迎为 **GPUMDkit** 贡献力量！参与方式：

- 通过 [Pull Requests](https://github.com/zhyan0603/GPUMDkit/pulls) 贡献 Python/Shell 脚本。
- 通过 [issues](https://github.com/zhyan0603/GPUMDkit/issues) 报告问题或提出功能建议。
- 通过邮箱 [yanzihan@westlake.edu.cn](mailto:yanzihan@westlake.edu.cn) 联系我。

也欢迎加入我们的 QQ 群（[825696376](https://qun.qq.com/universal-share/share?ac=1&authKey=buBNi1ADDzIFF2oZ1yA5FywG3LA9EL9yKZmb%2BN2MMz7nNuuxTas54wH7BgPEqP0s&busi_data=eyJncm91cENvZGUiOiI4MjU2OTYzNzYiLCJ0b2tlbiI6IlRxL1RLTDlOK3U2ekRSUXJ1TkNTUWd3ODNVV3BrdG9HN2lWWmJKMHAraGlDNzBZWFFyRUY2dUlSaW8rbUd4MisiLCJ1aW4iOiIxNDg5NjQ3MTc5In0%3D&data=fa4zSsT_IdI4ftCT_wwpytYHf--TaTB35lH0Jac5JHVpYoyXw3_3bZ1l1NZejsOZnGJku5u3BCbf5_bgrCkhZg&svctype=4&tempid=h5_group_info)）。一起创造有用的工具！🌟

## 引用

**GPUMDkit** 是一款面向所有人的开源工具。如果在你的研究或工作中它有所帮助，欢迎 ⭐ [在 GitHub 上给我们点亮 Star](https://github.com/zhyan0603/GPUMDkit)。此外，如果 GPUMDkit 对你的发表工作有贡献，请引用我们的论文：

> Z. Yan\*, D. Li, X. Wu, Z. Liu, C. Hua, B. Situ, H. Yang, S. Tang, B. Tang, Z. Wang, S. Yi, H. Wang, D. Huang, K. Li, Q. Guo, Z. Chen, K. Xu, Y. Wang, Z. Wang, G. Tang, S. Liu, Z. Fan, and Y. Zhu\*. **GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP**. [MGE Advances, 2026, 4, e70074](https://doi.org/10.1002/mgea.70074).

在论文中可以参考类似表述：

> Data processing and figure generation were performed using GPUMDkit [x].
