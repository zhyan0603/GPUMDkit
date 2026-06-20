# 分析工具指南

<div align="center">
  <p>
    <strong>中文</strong> | <a href="../en/analyzers.md">English</a>
  </p>
</div>

本指南介绍 GPUMDkit 中用于结构验证、过滤和质量控制的所有分析工具。

## 可用工具

| 工具 | 命令 | 描述 |
|------|------|------|
| 成分分析 | `-analyze_comp` | 按成分分组结构 |
| 异常值检测 | 菜单 502 | 查找高误差结构 |
| 化学物种 | 菜单 503 | 列出唯一元素 |
| 电荷平衡 | `-cbc` | 检查氧化态平衡 |
| 属性范围 | `-range` | 能量/力/应力统计 |
| 距离过滤 | 菜单 506 | 按最小距离过滤（无 PBC） |
| 距离过滤 (PBC) | 菜单 506b | 按最小距离过滤（带 PBC） |
| 最小距离 | `-min_dist` | 计算最小距离（无 PBC） |
| 最小距离 (PBC) | `-min_dist_pbc` | 计算最小距离（带 PBC） |
| 概率密度 | 菜单 508 | 3D 扩散通道分析 |

## 交互模式

```bash
gpumdkit.sh
# 选择: 5) Analyzer
```

您将看到：

```
+------------------------------------------------------+
|                    ANALYZER TOOLS                    |
+------------------------------------------------------+
| 501) Analyze composition of extxyz                   |
| 502) Find outliers of extxyz                         |
| 503) Analyze chemical species of extxyz              |
| 504) Check charge balance of extxyz                  |
| 505) Analyze energy/force/virial range               |
| 506) Filter structures by minimum distance           |
| 507) Get minimum interatomic distance                |
| 508) Probability density analysis                    |
+------------------------------------------------------+
```

## 命令参考

### 成分分析

按化学成分分析和分组结构。

```bash
gpumdkit.sh -analyze_comp train.xyz
```

**输出：** 唯一成分的表格，包含原子数和结构数。可交互选择导出特定成分。

### 属性范围分析

计算 extxyz 文件中能量、力或应力的范围（最小值、最大值）。

```bash
gpumdkit.sh -range <文件名> <属性> [hist]

# 示例
gpumdkit.sh -range train.xyz energy
gpumdkit.sh -range train.xyz force
gpumdkit.sh -range train.xyz virial
gpumdkit.sh -range train.xyz force hist    # 带直方图
```

**输出：** 最小/最大值。使用 `hist` 参数可选择直方图。

### 最小距离

计算每个元素对的最小原子间距离。

```bash
# 快速计算（无 PBC）
gpumdkit.sh -min_dist dump.xyz

# 精确计算（带 PBC）
gpumdkit.sh -min_dist_pbc dump.xyz
```

**输出：** 所有元素对的最小距离表

### 电荷平衡检查

检查所有结构的氧化态平衡。

```bash
gpumdkit.sh -cbc train.xyz
```

**输出：**
- `balanced.xyz` - 电荷平衡的结构
- `unbalanced.xyz` - 不平衡的结构
- `indices.txt` - 摘要

### 化学物种识别

列出文件中所有唯一的化学元素。

```bash
# 交互模式
gpumdkit.sh  # 选择: 5) Analyzer -> 503
```

**输出：** 排序后的所有元素列表

### 结构过滤

#### 按最小距离（无 PBC）

```bash
gpumdkit.sh  # 选择: 5) Analyzer -> 506
```

**输入：** extxyz 文件，距离阈值
**输出：** `filtered_<file>.xyz`，`filtered_out_<file>.xyz`

#### 按最小距离（带 PBC）

```bash
gpumdkit.sh  # 选择: 5) Analyzer -> 506b
```

对周期性系统更精确。

#### 按元素对距离范围

```bash
gpumdkit.sh -filter_range <文件> <元素1> <元素2> <最小距离> <最大距离>

# 示例：过滤 Li-Li 距离在 1.9 到 2.0 Å 之间
gpumdkit.sh -filter_range dump.xyz Li Li 1.9 2.0
```

**输出：** `filtered_<elem1>_<elem2>_<min>_<max>.xyz`

#### 按盒子大小

```bash
gpumdkit.sh -filter_box <文件> <边缘限制>

# 示例：过滤盒子边缘 > 20 Å 的结构
gpumdkit.sh -filter_box dump.xyz 20
```

**输出：** `filtered_by_box.xyz`

#### 按属性值

```bash
gpumdkit.sh -filter_value <文件> <属性> <阈值>

# 示例：保留力 < 20 eV/Å 的结构
gpumdkit.sh -filter_value train.xyz force 20
```

**输出：** `filtered.xyz`

### 异常值检测

基于 RMSE 阈值查找异常值结构。

```bash
gpumdkit.sh  # 选择: 5) Analyzer -> 502
```

**必需文件：**
- `train.xyz`
- `energy_train.out`
- `force_train.out`
- `stress_train.out`

**输出：** `selected.xyz`（高误差），`remained.xyz`（低误差）

### 概率密度分析

从 AIMD 轨迹计算移动离子的 3D 概率密度。

```bash
gpumdkit.sh  # 选择: 5) Analyzer -> 508
```

**参数：**
- 参考结构（POSCAR）
- 轨迹文件（extxyz）
- 移动物种（如 Li）
- 网格间距（如 0.25 Å）

**输出：** `probability_density_<interval>.vasp`（CHGCAR 格式）

## 常见工作流程

### 数据质量检查

```bash
# 1. 检查成分
gpumdkit.sh -analyze_comp train.xyz

# 2. 检查最小距离
gpumdkit.sh -min_dist_pbc train.xyz

# 3. 检查力范围
gpumdkit.sh -range train.xyz force

# 4. 检查异常值
gpumdkit.sh  # 选择: 5) Analyzer -> 502
```

### 结构过滤管道

```bash
# 1. 按距离过滤
gpumdkit.sh  # 选择: 5) Analyzer -> 506

# 2. 按盒子大小过滤
gpumdkit.sh -filter_box filtered.xyz 20

# 3. 按力值过滤
gpumdkit.sh -filter_value filtered_by_box.xyz force 15
```

### 扩散通道分析

```bash
# 1. 计算概率密度
gpumdkit.sh  # 选择: 5) Analyzer -> 508

# 2. 使用 VESTA 或类似软件可视化
# 打开 probability_density_0.25.vasp
```

## 附加工具

### 时间监控

```bash
# 监控 GPUMD 进度
gpumdkit.sh -time gpumd

# 监控 NEP 训练进度
gpumdkit.sh -time nep
```

### 体积分析

```bash
# 获取每个温度的平均体积
python Scripts/analyzer/get_volume.py
# 需要：*/K/thermo.out 目录
```

## CLI 参考

| 标志 | 描述 | 语法 |
|------|------|------|
| `-analyze_comp` | 成分分析 | `gpumdkit.sh -analyze_comp <file>` |
| `-range` | 属性范围 | `gpumdkit.sh -range <file> <prop>` |
| `-min_dist` | 最小距离（无 PBC） | `gpumdkit.sh -min_dist <file>` |
| `-min_dist_pbc` | 最小距离（PBC） | `gpumdkit.sh -min_dist_pbc <file>` |
| `-cbc` | 电荷平衡 | `gpumdkit.sh -cbc <file>` |
| `-filter_range` | 按距离范围过滤 | `gpumdkit.sh -filter_range <file> <e1> <e2> <min> <max>` |
| `-filter_box` | 按盒子大小过滤 | `gpumdkit.sh -filter_box <file> <limit>` |
| `-filter_value` | 按属性过滤 | `gpumdkit.sh -filter_value <file> <prop> <thresh>` |
| `-time` | 时间监控 | `gpumdkit.sh -time <gpumd\|nep>` |

## 依赖

| 工具 | 必需的包 |
|------|----------|
| 所有 Python 脚本 | `ase`, `numpy` |
| 距离计算 | `scipy` |
| 电荷平衡 | `pymatgen`, `tqdm` |
| 属性范围 | `matplotlib` |
| 概率密度 | `pymatgen` |

## 参见

- [格式转换](format_conversion.md) - 转换文件格式
- [计算器](calculators.md) - 计算属性
- [可视化](visualization.md) - 绘制结果
