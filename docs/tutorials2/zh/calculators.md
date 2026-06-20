# 计算器指南

<div align="center">
  <p>
    <strong>中文</strong> | <a href="../en/calculators.md">English</a>
  </p>
</div>

本指南介绍 GPUMDkit 中用于计算材料属性的所有计算器工具。

## 可用计算器

| 计算器 | 命令 | 描述 |
|--------|------|------|
| 离子电导率 | `-calc ionic-cond` | 从 MSD 计算离子电导率 |
| NEP 属性 | `-calc nep` | 使用 NEP 计算能量/力/应力 |
| 描述符 | `-calc des` | 计算 NEP 描述符用于分析 |
| DOAS | `-calc doas` | 原子态密度 |
| NEB | 直接 Python | 推移弹性带计算 |
| 邻居列表 | `-calc nlist` | 构建用于分析的邻居列表 |
| 位移 | `-calc disp` | 计算原子位移 |
| 平均结构 | `-calc avg-struct` | 从轨迹计算时间平均结构 |
| 八面体倾斜 | `-calc oct-tilt` | 钙钛矿八面体倾斜角 |
| 极化 | `-calc pol-abo3` | ABO3 局部极化 |
| 最小化 | `-calc minimize` | 使用 NEP 进行结构弛豫 |
| MSD | `-calc msd` | 从轨迹计算均方位移 |

## 交互模式

```bash
gpumdkit.sh
# 选择: 4) Calculators
```

您将看到：

```
+----------------------------------------------------------+
|                     CALCULATOR TOOLS                     |
+----------------------------------------------------------+
| 401) Calc ionic conductivity                             |
| 402) Calc properties by nep                              |
| 403) Calc descriptors of specific elements               |
| 404) Calc density of atomistic states (DOAS)             |
| 405) Calc nudged elastic band (NEB) by nep               |
| 406) Build neighbor list                                 |
| 407) Calc displacement from trajectory                   |
| 408) Calc averaged structure                             |
| 409) Calc octahedral tilt                                |
| 410) Calc polarization for ABO3                          |
| 411) Minimize structure by nep                           |
| 412) Calc mean square displacement (MSD) from trajectory |
+----------------------------------------------------------+
```

## 命令参考

### 离子电导率

使用 Nernst-Einstein 方程从 MSD 数据计算离子电导率。

```bash
gpumdkit.sh -calc ionic-cond <元素> <电荷>

# 示例
gpumdkit.sh -calc ionic-cond Li 1    # 锂离子 (Li+)
gpumdkit.sh -calc ionic-cond Na 1    # 钠离子 (Na+)
gpumdkit.sh -calc ionic-cond O -2    # 氧离子 (O2-)
```

**当前目录中需要的文件：**
- `msd.out` - 来自 GPUMD 的 `compute_msd`
- `thermo.out` - 用于获取温度
- `model.xyz` - 用于获取体积
- `run.in` - 用于获取模拟参数

**输出：** 离子扩散率和电导率值

### NEP 属性预测

使用 NEP 模型作为 DFT 代理计算能量、力和应力。

```bash
gpumdkit.sh -calc nep <input.xyz> <output.xyz> <nep.txt>

# 示例
gpumdkit.sh -calc nep structures.xyz predictions.xyz nep.txt
```

**注意：** 先使用 `gpumdkit.sh -clean_xyz` 清理输入以删除现有属性。

### 描述符

计算 NEP 描述符用于降维和结构分析。

```bash
gpumdkit.sh -calc des <input.xyz> <output.npy> <nep.txt> <元素>

# 示例
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li
```

**可视化：**
```bash
gpumdkit.sh -plt des pca    # PCA 可视化
gpumdkit.sh -plt des umap   # UMAP 可视化
```

### 原子态密度 (DOAS)

按元素类型分组计算每原子能量分布。

```bash
gpumdkit.sh -calc doas <input.xyz> <nep.txt> <output.txt>

# 示例
gpumdkit.sh -calc doas structures.xyz nep.txt doas.out
```

**可视化：**
```bash
gpumdkit.sh -plt doas doas.out Li
```

### 均方位移 (MSD)

从 extxyz 轨迹计算方向性 MSD。

```bash
gpumdkit.sh -calc msd <trajectory.xyz> <元素> <dt_fs> [max_corr_steps]

# 示例：Li，时间步长 10 fs
gpumdkit.sh -calc msd dump.xyz Li 10
```

**输出：** `msd.out`（时间/ps, MSD_x, MSD_y, MSD_z）

### 邻居列表

为钙钛矿分析构建邻居列表。

```bash
gpumdkit.sh -calc nlist -i <input> -c <cutoff> -n <num_neighbors> -C <中心元素> -E <邻居元素>

# 示例
gpumdkit.sh -calc nlist -i model.xyz -c 4 -n 12 -C Ti -E O
gpumdkit.sh -calc nlist -i model.xyz -c 4 -n 12 -C Pb Sr -E O
```

**输出：** `nl-<center>-<neighbor>.dat`

### 位移

使用邻居列表从轨迹计算原子位移。

```bash
gpumdkit.sh -calc disp -i <trajectory.xyz> -n <邻居列表> -o <输出>

# 示例
gpumdkit.sh -calc disp -i movie.xyz -n nl-Pb-O.dat -o displacements.dat
```

**可选帧切片：** `-s <start> -t <stop> -p <step> -l <last_fraction>`

### 平均结构

从轨迹生成时间平均结构。

```bash
gpumdkit.sh -calc avg-struct -i <trajectory.xyz> -l <比例> -o <输出>

# 示例：平均最后 20% 的帧
gpumdkit.sh -calc avg-struct -i movie.xyz -l 0.2 -o averaged.xyz
```

### 八面体倾斜

计算钙钛矿体系的八面体倾斜角。

```bash
gpumdkit.sh -calc oct-tilt -i <input.xyz> -n <B-O 邻居列表> -o <输出>

# 示例
gpumdkit.sh -calc oct-tilt -i model.xyz -n nl-Ti-O.dat -o octahedral_tilt.dat
```

### 极化 (ABO3)

计算 ABO3 钙钛矿的局部极化。

```bash
gpumdkit.sh -calc pol-abo3 -i <input.xyz> \
  --nl-ba <B-A 邻居列表> \
  --nl-bo <B-O 邻居列表> \
  --bec <元素=电荷 ...>

# 示例
gpumdkit.sh -calc pol-abo3 -i model.xyz \
  --nl-ba nl-Ti-Pb.dat \
  --nl-bo nl-Ti-O.dat \
  --bec Pb=2.5 Sr=2.0 Ti=4.0 O=-2.0
```

### 结构最小化

使用 BFGS 优化器和 NEP 模型弛豫结构。

```bash
gpumdkit.sh -calc minimize <结构> <nep.txt> [fmax] [max_steps]

# 示例
gpumdkit.sh -calc minimize POSCAR nep.txt 0.01 1000
```

**参数：**
- `fmax`：力收敛阈值（默认：0.01 eV/Å）
- `max_steps`：最大优化步数（默认：1000）

**输出：** `minimize.xyz`（优化轨迹）

### NEB 计算

进行推移弹性带计算以获得迁移势垒。

```bash
# 直接 Python 执行
python Scripts/calculators/neb_calculation.py <initial.xyz> <final.xyz> <n_images> <nep.txt>

# 示例
python Scripts/calculators/neb_calculation.py init.xyz fin.xyz 9 nep.txt
```

## 常见工作流程

### 离子传输分析

```bash
# 1. 在 run.in 中使用 compute_msd 运行 MD
# 2. 计算 MSD
gpumdkit.sh -calc msd dump.xyz Li 10

# 3. 计算电导率
gpumdkit.sh -calc ionic-cond Li 1

# 4. 可视化
gpumdkit.sh -plt msd
gpumdkit.sh -plt sdc
```

### 描述符分析

```bash
# 1. 计算描述符
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li

# 2. 使用 PCA 可视化
gpumdkit.sh -plt des pca

# 3. 或使用 UMAP 可视化
gpumdkit.sh -plt des umap
```

### 钙钛矿分析

```bash
# 1. 构建邻居列表
gpumdkit.sh -calc nlist -i model.xyz -c 4 -n 12 -C Ti -E O
gpumdkit.sh -calc nlist -i model.xyz -c 4 -n 12 -C Ti -E Pb

# 2. 计算位移
gpumdkit.sh -calc disp -i movie.xyz -n nl-Ti-O.dat -o disp.dat

# 3. 计算倾斜
gpumdkit.sh -calc oct-tilt -i movie.xyz -n nl-Ti-O.dat -o tilt.dat

# 4. 计算极化
gpumdkit.sh -calc pol-abo3 -i movie.xyz \
  --nl-ba nl-Ti-Pb.dat --nl-bo nl-Ti-O.dat \
  --bec Pb=2.5 Ti=4.0 O=-2.0
```

### 结构弛豫

```bash
# 1. 最小化结构
gpumdkit.sh -calc minimize POSCAR nep.txt 0.01 1000

# 2. 使用最小化后的结构进行 MD
cp minimize.xyz relaxed_model.xyz
```

## 依赖

| 计算器 | 必需的包 |
|--------|----------|
| 所有 | `numpy`, `ase` |
| ionic-cond | `scipy` |
| nep, des, doas, minimize | `calorine` |
| nlist, disp, oct-tilt, pol-abo3 | `ferrodispcalc` |
| neb | `calorine`, `matplotlib` |

安装依赖：
```bash
pip install numpy ase scipy calorine matplotlib tqdm
pip3 install git+https://github.com/MoseyQAQ/ferrodispcalc.git
```

## 参见

- [可视化](visualization.md) - 绘制计算结果
- [分析工具](analyzers.md) - 结构分析
- [NEP 训练指南](nep_training.md) - 完整的训练工作流
