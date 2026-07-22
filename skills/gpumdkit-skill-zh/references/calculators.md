# 计算器

## 目录

- 可用计算器
- 命令参考
- 常用工作流
- 依赖

## 可用计算器

| 计算器 | 命令 | 描述 |
|------------|---------|-------------|
| 离子电导率 | `-calc ionic-cond` | 从 MSD 计算离子电导率 |
| NEP 性质 | `-calc nep` | 使用 NEP 计算能量/力/应力 |
| 描述符 | `-calc des` | 计算 NEP 描述符用于分析 |
| DOAS | `-calc doas` | 原子态密度 |
| NEB | `-calc neb` | 弹性带路径计算 |
| 邻居列表 | `-calc nlist` | 构建邻居列表用于分析 |
| 位移 | `-calc disp` | 计算原子位移 |
| 平均结构 | `-calc avg-struct` | 轨迹的时间平均结构 |
| 八面体倾斜 | `-calc oct-tilt` | 钙钛矿八面体倾斜角 |
| 极化 | `-calc pol-abo3` | ABO3 局部极化 |
| 最小化 | `-calc minimize` | 使用 NEP 进行结构弛豫 |
| MSD | `-calc msd` | 从轨迹计算均方位移 |

## 命令参考

### 离子电导率
```bash
# 从 MSD 数据计算离子电导率
gpumdkit.sh -calc ionic-cond <element> <charge>

# 示例：锂离子（Li+）
gpumdkit.sh -calc ionic-cond Li 1

# 示例：氧离子（O2-）
gpumdkit.sh -calc ionic-cond O -2

# 推荐的非交互式输入（当前目录）：
# - msd.out（来自 GPUMD compute_msd）
# - thermo.out（用于温度）
# - model.xyz（用于体积）
# - run.in（可选；用于检测复制因子）
```

### NEP 性质预测
```bash
# 使用 NEP 模型计算性质
gpumdkit.sh -calc nep <input.xyz> <output.xyz> <nep.txt>

# 示例
gpumdkit.sh -calc nep structures.xyz predictions.xyz nep.txt

# 注意：先用 gpumdkit.sh -clean_xyz 清理输入以移除已有性质
```

### 描述符
```bash
# 计算特定元素的 NEP 描述符
gpumdkit.sh -calc des <input.xyz> <output.npy> <nep.txt> <element>

# 示例
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li

# 可视化：gpumdkit.sh -plt des pca
# 或：gpumdkit.sh -plt des umap
```

### 原子态密度（DOAS）
```bash
# 计算 DOAS
gpumdkit.sh -calc doas <input.xyz> <nep.txt> <output.txt>

# 示例
gpumdkit.sh -calc doas structures.xyz nep.txt doas.out

# 可视化：gpumdkit.sh -plt doas doas.out Li
```

### 均方位移
```bash
# 从轨迹计算 MSD
gpumdkit.sh -calc msd <trajectory.xyz> <element> <dt_fs> [max_corr_steps]

# 示例：Li，时间步长 10 fs
gpumdkit.sh -calc msd dump.xyz Li 10

# 输出：msd.out（Time/ps、MSD_x、MSD_y、MSD_z）
```

### 邻居列表
```bash
# 为钙钛矿分析构建邻居列表
gpumdkit.sh -calc nlist -i <input> -c <cutoff> -n <num_neighbors> -C <center_elements> -E <neighbor_elements>

# 示例：BaTiO3 中的 Ti-O 邻居
gpumdkit.sh -calc nlist -i model.xyz -c 4 -n 12 -C Ti -E O

# 示例：PZT 中的 Pb/Sr-O 邻居
gpumdkit.sh -calc nlist -i model.xyz -c 4 -n 12 -C Pb Sr -E O

# 输出：nl-<center>-<neighbor>.dat
```

### 位移
```bash
# 从轨迹计算位移
gpumdkit.sh -calc disp -i <trajectory.xyz> -n <neighbor_list> -o <output>

# 示例
gpumdkit.sh -calc disp -i movie.xyz -n nl-Pb-O.dat -o displacements.dat

# 可选帧切片：-s <start> -t <stop> -p <step> -l <last_fraction>
```

### 平均结构
```bash
# 计算时间平均结构
gpumdkit.sh -calc avg-struct -i <trajectory.xyz> -l <fraction> -o <output>

# 示例：平均最后 20% 的帧
gpumdkit.sh -calc avg-struct -i movie.xyz -l 0.2 -o averaged.xyz
```

### 八面体倾斜
```bash
# 计算八面体倾斜角
gpumdkit.sh -calc oct-tilt -i <input.xyz> -n <B-O neighbor list> -o <output>

# 示例
gpumdkit.sh -calc oct-tilt -i model.xyz -n nl-Ti-O.dat -o octahedral_tilt.dat
```

### 极化（ABO3）
```bash
# 计算钙钛矿的局部极化
gpumdkit.sh -calc pol-abo3 -i <input.xyz> \
  --nl-ba <B-A neighbor list> \
  --nl-bo <B-O neighbor list> \
  --bec <Element=charge ...>

# 示例
gpumdkit.sh -calc pol-abo3 -i model.xyz \
  --nl-ba nl-Ti-Pb.dat \
  --nl-bo nl-Ti-O.dat \
  --bec Pb=2.5 Sr=2.0 Ti=4.0 O=-2.0
```

### 结构最小化
```bash
# 使用 NEP 最小化结构
gpumdkit.sh -calc minimize <structure> <nep.txt> [fmax] [max_steps]

# 示例
gpumdkit.sh -calc minimize POSCAR nep.txt 0.01 1000

# 输出：minimize.xyz（优化轨迹）
# 默认 fmax：0.01 eV/A，默认 max_steps：1000
```

### NEB 计算
```bash
# CLI 快捷方式
gpumdkit.sh -calc neb <initial.xyz> <final.xyz> <n_images> <nep.txt>

# 示例
gpumdkit.sh -calc neb init.xyz fin.xyz 9 nep.txt

# 直接 Python 执行（也可以）
python Scripts/calculators/neb_calculation.py init.xyz fin.xyz 9 nep.txt

# 使用 NepTrainKit 的替代方案：
python Scripts/calculators/neb_calculation_neptrain.py init.xyz fin.xyz 9 nep.txt
```

## 常用工作流

### 离子输运分析
```bash
# 路径 A：从 extxyz 轨迹推导 msd.out
gpumdkit.sh -calc msd dump.xyz Li 10

# 路径 B：使用 GPUMD compute_msd 直接生成 msd.out。
# 除非比较实现，否则不要同时运行两条路径。

# 电导率需要验证过的 msd.out、thermo.out、model.xyz 和 run.in
gpumdkit.sh -calc ionic-cond Li 1

# 可视化
gpumdkit.sh -plt msd
gpumdkit.sh -plt sdc
```

### 描述符分析
```bash
# 1. 计算描述符
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li
# 2. PCA 可视化
gpumdkit.sh -plt des pca
# 3. 或 UMAP 可视化
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
gpumdkit.sh -calc pol-abo3 -i movie.xyz --nl-ba nl-Ti-Pb.dat --nl-bo nl-Ti-O.dat --bec Pb=2.5 Ti=4.0 O=-2.0
```

### 结构弛豫
```bash
# 1. 最小化结构
gpumdkit.sh -calc minimize POSCAR nep.txt 0.01 1000
# 2. 检查结果
gpumdkit.sh -plt thermo  # 如果弛豫后运行 MD
```

## 依赖

| 计算器 | 所需包 |
|------------|------------------|
| 全部 | `numpy`、`ase` |
| ionic-cond | `scipy` |
| nep、des、doas、minimize | `calorine` |
| nlist、disp、oct-tilt、pol-abo3 | `ferrodispcalc` |
| neb | `calorine`、`matplotlib` |

安装依赖：
```bash
pip install numpy ase scipy calorine matplotlib tqdm
pip3 install git+https://github.com/MoseyQAQ/ferrodispcalc.git
```

## 详细文档

用户指南请参见 `${GPUMDkit_path}/docs/tutorials/en/calculator_scripts.md` 或中文对应版本。
