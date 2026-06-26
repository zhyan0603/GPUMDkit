<div align="center">
  <h1>快速入门</h1>
  <p>
    <strong>简体中文</strong> | <a href="../en/quickstart.md">English</a>
  </p>
</div>

5 分钟内开始使用 GPUMDkit。

## 前提条件

- Linux 或 macOS（Windows 通过 WSL）
- Python 3.8+（推荐 3.12）
- Git

## 安装

### 步骤 1：克隆

```bash
git clone https://github.com/zhyan0603/GPUMDkit.git
cd GPUMDkit
```

### 步骤 2：安装

```bash
source ./install.sh
```

### 步骤 3：重新加载 Shell

```bash
source ~/.bashrc  # 或 source ~/.zshrc
```

### 步骤 4：验证

```bash
gpumdkit.sh -h
```

## Python 依赖

```bash
conda create -n gpumdkit python=3.12
conda activate gpumdkit
pip install neptrain ase pymatgen dpdata numpy scipy matplotlib
```

> 使用 GPUMDkit 前请确保已激活 `gpumdkit` 环境。

## 使用模式

### 交互模式

```bash
gpumdkit.sh
```

```
           ____ ____  _   _ __  __ ____  _    _ _
          / ___|  _ \| | | |  \/  |  _ \| | _(_) |_
         | |  _| |_) | | | | |\/| | | | | |/ / | __|
         | |_| |  __/| |_| | |  | | |_| |   <| | |_
          \____|_|    \___/|_|  |_|____/|_|\_\_|\__|

          GPUMDkit Version 1.5.5 (dev) (2026-05-10)

 ---------------------- GPUMD ------------------------
  1) Format Conversion          2) Sample Structures
  3) Workflow                   4) Calculators
  5) Analyzer                   6) Visualization
  7) Utilities                  0) Exit
 ------------>>
```

### 命令行模式

```bash
gpumdkit.sh -h                    # 显示帮助
gpumdkit.sh -pos2exyz POSCAR model.xyz  # 直接命令
```

## 第一个例子

### 转换结构

```bash
gpumdkit.sh -pos2exyz POSCAR model.xyz
gpumdkit.sh -addgroup POSCAR Li Y Cl
```

### 分析数据

```bash
gpumdkit.sh -range train.xyz force
gpumdkit.sh -min_dist_pbc train.xyz
gpumdkit.sh -analyze_comp train.xyz
```

### 绘制结果

```bash
gpumdkit.sh -plt train
gpumdkit.sh -plt msd
gpumdkit.sh -plt rdf
```

## 获取帮助

```bash
gpumdkit.sh -h              # 通用帮助
gpumdkit.sh -plt -h         # 绘图帮助
gpumdkit.sh -calc -h        # 计算器帮助
gpumdkit.sh -pos2exyz -h    # 特定命令帮助
```

## 故障排除

**命令未找到**
```bash
export GPUMDkit_path=/path/to/GPUMDkit
export PATH=${GPUMDkit_path}:${PATH}
```

**缺少包**
```bash
conda activate gpumdkit
```

## 下一步

- [格式转换](format_conversion.md)
- [计算器](calculators.md)
- [可视化](visualization.md)
- [NEP 训练](nep_training.md)
