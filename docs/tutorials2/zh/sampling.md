# 结构采样指南

<div align="center">
  <p>
    <strong>中文</strong> | <a href="../en/sampling.md">English</a>
  </p>
</div>

本指南介绍 GPUMDkit 中的结构采样和选择工具。

## 可用方法

| 方法 | 菜单 | 描述 |
|------|------|------|
| 均匀/随机 | 201 | 选择均匀间隔或随机帧 |
| FPS by PyNEP | 202 | 最远点采样（已弃用） |
| FPS by NepTrain | 203 | 最远点采样（推荐） |
| 扰动 | 204 | 生成扰动结构 |
| 力偏差 | 205 | 选择高力偏差结构 |

## 交互模式

```bash
gpumdkit.sh
# 选择: 2) Sample Structures
```

您将看到：

```
+------------------------------------------------------+
|                 SAMPLE STRUCTURE TOOLS               |
+------------------------------------------------------+
| 201) Sample structures from extxyz                   |
| 202) FPS sampling by PyNEP [deprecated]              |
| 203) FPS sampling by NepTrain [preferred]            |
| 204) Perturb structure                               |
| 205) Select max force deviation structs              |
+------------------------------------------------------+
```

## 命令参考

### 均匀/随机采样

使用均匀或随机采样选择结构。

```bash
python Scripts/sample_structures/sample_structures.py <input.xyz> <method> <num_samples> [skip_initial]

# 示例
python Scripts/sample_structures/sample_structures.py dump.xyz uniform 50
python Scripts/sample_structures/sample_structures.py dump.xyz random 100 500
```

**参数：**
- `input.xyz`：输入轨迹文件
- `method`：`uniform` 或 `random`
- `num_samples`：要选择的帧数
- `skip_initial`：（可选）跳过前 N 帧

**工作原理：**
- **均匀**：使用 `numpy.linspace` 选择均匀间隔的索引
- **随机**：使用 `numpy.random.choice` 无放回地选择随机索引

**输出：** `sampled_structures.xyz`

### 最远点采样 (FPS) - NepTrain（推荐）

使用 NEP 描述符上的 FPS 选择多样化结构。

```bash
python Scripts/sample_structures/neptrain_select_structs.py <sample.xyz> <train.xyz> <nep.txt>

# 示例
python Scripts/sample_structures/neptrain_select_structs.py dump.xyz train.xyz nep.txt
```

**选择方法（交互式）：**
1. **最小距离**：选择直到与现有训练集的最大距离低于阈值
2. **结构数量**：选择 `min_select` 到 `max_select` 个结构

**输出：**
- `selected.xyz`
- `select.png`（PCA 可视化）
- `pca_sample.txt`, `pca_train.txt`, `pca_selected.txt`

**引用**：Chen et al., Comput. Phys. Commun., 2025, 317, 109859

### 最远点采样 (FPS) - PyNEP（已弃用）

```bash
python Scripts/sample_structures/pynep_select_structs.py <sample.xyz> <train.xyz> <nep.txt>
```

**注意**：PyNEP 包不再积极维护。请使用 NepTrain（方法 203）。

### 结构扰动

为训练数据生成扰动结构。

```bash
python Scripts/sample_structures/perturb_structure.py <input.vasp> <pert_num> <cell_pert> <atom_pert> <style>

# 示例
python Scripts/sample_structures/perturb_structure.py POSCAR 20 0.03 0.2 uniform
```

**参数：**
- `input.vasp`：POSCAR/CONTCAR 文件
- `pert_num`：扰动结构数量（如 20）
- `cell_pert`：晶胞扰动比例（如 0.03 = 3%）
- `atom_pert`：原子扰动距离（Å）（如 0.2）
- `style`：`normal`、`uniform` 或 `const`

**输出：** `POSCAR_01.vasp`, `POSCAR_02.vasp`, ..., `POSCAR_<N>.vasp`

### 力偏差选择

从主动学习中选择具有高力偏差的结构。

```bash
python Scripts/sample_structures/select_max_modev.py <top_n> <min_deviation>

# 示例
python Scripts/sample_structures/select_max_modev.py 200 0.15
```

**必需文件：**
- `active.out` - 每个结构的最大力偏差（来自 GPUMD 的 `active` 命令）
- `active.xyz` - 对应的结构

**输出：** `selected.xyz`

### 帧范围提取

从轨迹中提取帧范围。

```bash
python Scripts/sample_structures/frame_range.py <input.xyz> <start_fraction> <end_fraction>

# 示例：提取前 80% 的帧
python Scripts/sample_structures/frame_range.py dump.xyz 0 0.8
```

**输出：** `dump_0.00_0.80.xyz`

## 常见工作流程

### 训练数据准备

```bash
# 1. 运行 MD 模拟
# ... 生成轨迹 ...

# 2. 按距离过滤
gpumdkit.sh -min_dist_pbc dump.xyz
gpumdkit.sh -filter_dist_pbc dump.xyz 1.0

# 3. 采样多样化结构
gpumdkit.sh  # 选择: 2) Sample Structures -> 203

# 4. 添加到训练集
cat selected.xyz >> train.xyz
```

### 主动学习数据选择

```bash
# 1. 使用当前 NEP 运行 MD
# ... 生成 active.out 和 active.xyz ...

# 2. 选择高偏差结构
python Scripts/sample_structures/select_max_modev.py 100 0.1

# 3. 或使用 FPS 获取多样性
python Scripts/sample_structures/neptrain_select_structs.py active.xyz train.xyz nep.txt
```

### 初始数据的结构扰动

```bash
# 1. 从平衡结构开始
# 2. 生成扰动
python Scripts/sample_structures/perturb_structure.py POSCAR 50 0.05 0.3 uniform

# 3. 对扰动结构运行 DFT
# ... DFT 计算 ...

# 4. 添加到训练集
```

## 依赖

| 方法 | 必需的包 |
|------|----------|
| 均匀/随机 | `numpy`, `ase` |
| FPS (NepTrain) | `numpy`, `ase`, `matplotlib`, `scikit-learn`, `scipy`, `NepTrain` |
| FPS (PyNEP) | `numpy`, `ase`, `matplotlib`, `scikit-learn`, `pynep` |
| 扰动 | `dpdata` |
| 力偏差 | `numpy`, `ase` |

安装依赖：
```bash
pip install numpy ase matplotlib scikit-learn scipy dpdata
pip install neptrain  # 用于 FPS by NepTrain
```

## 重要提示

1. **FPS 是首选**：最远点采样提供比均匀/随机更好的多样性
2. **NepTrain 优于 PyNEP**：NepTrain 积极维护；PyNEP 已弃用
3. **采样前先过滤**：首先删除不合理的结构（太近、太大）
4. **帧索引从 0 开始**：第一帧索引为 0
5. **输出是 extxyz**：所有采样方法输出 extxyz 格式

## 参见

- [格式转换](format_conversion.md) - 转换文件格式
- [分析工具](analyzers.md) - 过滤和分析结构
- [NEP 训练指南](nep_training.md) - 完整的训练工作流
