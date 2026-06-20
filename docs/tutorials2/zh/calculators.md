<div align="center">
  <h1>计算器</h1>
  <p>
    <strong>简体中文</strong> | <a href="../en/calculators.md">English</a>
  </p>
</div>

从分子动力学数据计算材料属性。

## 可用计算器

| 计算器 | 命令 | 描述 |
|--------|------|------|
| 离子电导率 | `-calc ionic-cond` | 从 MSD 计算离子电导率 |
| NEP 属性 | `-calc nep` | 使用 NEP 计算能量/力/应力 |
| 描述符 | `-calc des` | 计算 NEP 描述符用于分析 |
| DOAS | `-calc doas` | 原子态密度 |
| MSD | `-calc msd` | 从轨迹计算均方位移 |
| 邻居列表 | `-calc nlist` | 构建用于分析的邻居列表 |
| 位移 | `-calc disp` | 计算原子位移 |
| 平均结构 | `-calc avg-struct` | 从轨迹计算时间平均结构 |
| 八面体倾斜 | `-calc oct-tilt` | 钙钛矿八面体倾斜角 |
| 极化 | `-calc pol-abo3` | ABO3 局部极化 |
| 最小化 | `-calc minimize` | 使用 NEP 进行结构弛豫 |
| NEB | 直接 Python | 推移弹性带计算 |

## 命令参考

### 离子电导率

```bash
gpumdkit.sh -calc ionic-cond <元素> <电荷>

# 示例
gpumdkit.sh -calc ionic-cond Li 1    # Li+
gpumdkit.sh -calc ionic-cond O -2    # O2-
```

需要文件：`msd.out`, `thermo.out`, `model.xyz`, `run.in`

### NEP 属性

```bash
gpumdkit.sh -calc nep <input.xyz> <output.xyz> <nep.txt>

# 示例
gpumdkit.sh -calc nep structures.xyz predictions.xyz nep.txt
```

### 描述符

```bash
gpumdkit.sh -calc des <input.xyz> <output.npy> <nep.txt> <元素>

# 示例
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li
```

可视化：`gpumdkit.sh -plt des pca` 或 `gpumdkit.sh -plt des umap`

### DOAS

```bash
gpumdkit.sh -calc doas <input.xyz> <nep.txt> <output.txt>

# 示例
gpumdkit.sh -calc doas structures.xyz nep.txt doas.out
```

### MSD

```bash
gpumdkit.sh -calc msd <trajectory.xyz> <元素> <dt_fs> [max_corr_steps]

# 示例
gpumdkit.sh -calc msd dump.xyz Li 10
```

输出：`msd.out`（时间/ps, MSD_x, MSD_y, MSD_z）

### 邻居列表

```bash
gpumdkit.sh -calc nlist -i <input> -c <cutoff> -n <num> -C <中心元素> -E <邻居元素>

# 示例
gpumdkit.sh -calc nlist -i model.xyz -c 4 -n 12 -C Ti -E O
```

### 位移

```bash
gpumdkit.sh -calc disp -i <trajectory.xyz> -n <邻居列表> -o <输出>

# 示例
gpumdkit.sh -calc disp -i movie.xyz -n nl-Pb-O.dat -o displacements.dat
```

### 平均结构

```bash
gpumdkit.sh -calc avg-struct -i <trajectory.xyz> -l <比例> -o <输出>

# 示例：平均最后 20% 的帧
gpumdkit.sh -calc avg-struct -i movie.xyz -l 0.2 -o averaged.xyz
```

### 八面体倾斜

```bash
gpumdkit.sh -calc oct-tilt -i <input.xyz> -n <B-O 列表> -o <输出>

# 示例
gpumdkit.sh -calc oct-tilt -i model.xyz -n nl-Ti-O.dat -o tilt.dat
```

### 极化 (ABO3)

```bash
gpumdkit.sh -calc pol-abo3 -i <input.xyz> --nl-ba <B-A 列表> --nl-bo <B-O 列表> --bec <元素=电荷 ...>

# 示例
gpumdkit.sh -calc pol-abo3 -i model.xyz --nl-ba nl-Ti-Pb.dat --nl-bo nl-Ti-O.dat --bec Pb=2.5 Ti=4.0 O=-2.0
```

### 最小化

```bash
gpumdkit.sh -calc minimize <结构> <nep.txt> [fmax] [max_steps]

# 示例
gpumdkit.sh -calc minimize POSCAR nep.txt 0.01 1000
```

输出：`minimize.xyz`

### NEB

```bash
python Scripts/calculators/neb_calculation.py <initial.xyz> <final.xyz> <n_images> <nep.txt>

# 示例
python Scripts/calculators/neb_calculation.py init.xyz fin.xyz 9 nep.txt
```

## 依赖

| 计算器 | 包 |
|--------|-----|
| 所有 | numpy, ase |
| ionic-cond | scipy |
| nep, des, doas, minimize | calorine |
| nlist, disp, oct-tilt, pol-abo3 | ferrodispcalc |
