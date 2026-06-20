# 快速入门指南

<div align="center">
  <p>
    <strong>中文</strong> | <a href="../en/quickstart.md">English</a>
  </p>
</div>

本指南将帮助您快速开始使用 GPUMDkit。

## 前提条件

在安装 GPUMDkit 之前，请确保您拥有：

- **Linux 或 macOS**（Windows 通过 WSL）
- **Bash** shell
- **Python 3.8+**（推荐 3.12）
- **Git**（用于克隆仓库）

## 安装

### 步骤 1：克隆仓库

```bash
git clone https://github.com/zhyan0603/GPUMDkit.git
cd GPUMDkit
```

### 步骤 2：运行安装脚本

```bash
source ./install.sh
```

此脚本将：
1. 将 GPUMDkit 添加到您的 PATH
2. 设置 bash 自动补全
3. 配置环境变量

### 步骤 3：重新加载 Shell

```bash
source ~/.bashrc
# 或
source ~/.zshrc
```

### 步骤 4：验证安装

```bash
gpumdkit.sh -h
```

您应该看到包含可用命令的帮助信息。

## Python 依赖

某些功能需要 Python 包。创建 conda 环境：

```bash
# 创建环境
conda create -n gpumdkit python=3.12
conda activate gpumdkit

# 安装核心包
pip install numpy ase scipy matplotlib

# 安装可选包（完整功能）
pip install neptrain pymatgen dpdata

# 用于钙钛矿分析（可选）
pip3 install git+https://github.com/MoseyQAQ/ferrodispcalc.git
```

## 两种使用模式

### 交互模式

启动交互式菜单：

```bash
gpumdkit.sh
```

您将看到：

```
       ____ ____  _   _ __  __ ____  _    _ _
      / ___|  _ \| | | |  \/  |  _ \| | _(_) |_
     | |  _| |_) | | | | |\/| | | | | |/ / | __|
     | |_| |  __/| |_| | |  | | |_| |   <| | |_
      \____|_|    \___/|_|  |_|____/|_|\_\_|\__|

      GPUMDkit Version 1.5.5 (dev) (2026-05-10)
Core Developer: Zihan YAN (yanzihan@westlake.edu.cn)

 ---------------------- GPUMD ------------------------
 1) Format Conversion          2) Sample Structures
 3) Workflow                   4) Calculators
 5) Analyzer                   6) Visualization
 7) Utilities                  8) Developing...
 0) Exit
 ------------>>
 Input the function number:
```

输入所需功能的编号进行导航。

### 命令行模式

直接执行命令：

```bash
# 显示帮助
gpumdkit.sh -h

# 将 POSCAR 转换为 extxyz
gpumdkit.sh -pos2exyz POSCAR model.xyz

# 绘制训练结果
gpumdkit.sh -plt train

# 计算离子电导率
gpumdkit.sh -calc ionic-cond Li 1
```

## 第一步示例

### 示例 1：转换结构

```bash
# 将 POSCAR 转换为 extxyz 格式
gpumdkit.sh -pos2exyz POSCAR model.xyz

# 添加组标签（NEP 训练必需）
gpumdkit.sh -addgroup POSCAR Li Y Cl
```

### 示例 2：分析训练数据

```bash
# 检查训练数据中的力范围
gpumdkit.sh -range train.xyz force

# 检查最小原子间距离
gpumdkit.sh -min_dist_pbc train.xyz

# 分析成分
gpumdkit.sh -analyze_comp train.xyz
```

### 示例 3：绘制结果

```bash
# 绘制 NEP 训练结果
gpumdkit.sh -plt train

# 绘制均方位移
gpumdkit.sh -plt msd

# 绘制径向分布函数
gpumdkit.sh -plt rdf
```

## 获取帮助

### 通用帮助

```bash
gpumdkit.sh -h
```

### 模块特定帮助

```bash
# 绘图帮助
gpumdkit.sh -plt -h

# 计算器帮助
gpumdkit.sh -calc -h

# 特定命令帮助
gpumdkit.sh -pos2exyz -h
gpumdkit.sh -calc ionic-cond -h
```

### 交互模式帮助

在交互模式中，每个子菜单都会显示可用选项及其描述。

## 常见问题

### 问题：命令未找到

**解决方案**：确保 GPUMDkit 在您的 PATH 中：

```bash
# 检查 gpumdkit.sh 是否在 PATH 中
which gpumdkit.sh

# 如果未找到，手动添加到 PATH
export GPUMDkit_path=/path/to/GPUMDkit
export PATH=${GPUMDkit_path}:${PATH}
```

### 问题：缺少 Python 包

**解决方案**：激活 conda 环境：

```bash
conda activate gpumdkit
```

### 问题：权限被拒绝

**解决方案**：使脚本可执行：

```bash
chmod +x /path/to/GPUMDkit/gpumdkit.sh
```

## 下一步

现在您已经安装了 GPUMDkit，请探索以下教程：

- [格式转换](format_conversion.md) - 学习在文件格式之间转换
- [计算器](calculators.md) - 计算材料属性
- [可视化](visualization.md) - 创建出版级图表
- [NEP 训练指南](nep_training.md) - 完整的训练工作流

## 更新 GPUMDkit

更新到最新版本：

```bash
gpumdkit.sh -update
```

或手动更新：

```bash
cd /path/to/GPUMDkit
git pull origin main
```
