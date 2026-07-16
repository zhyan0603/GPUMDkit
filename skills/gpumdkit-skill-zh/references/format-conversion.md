# 格式转换

## 目录

- 支持的格式
- 命令参考
- 详细 CLI 标志
- 示例
- 注意事项与依赖

## 支持的格式

| 格式 | 扩展名 | 描述 |
|--------|-----------|-------------|
| VASP | POSCAR、OUTCAR、XDATCAR | Vienna Ab initio Simulation Package |
| LAMMPS | .data、dump.* | Large-scale Atomic/Molecular Massively Parallel Simulator |
| CP2K | .log、pos.xyz、frc.xyz、cell.cell | 量子化学与固态物理 |
| ABACUS | running_scf.log、running_md.log | 基于原子轨道的从头计算 |
| CIF | .cif | 晶体学信息文件 |
| MTP | .cfg | Moment Tensor Potential 格式 |
| ASE | .traj | Atomic Simulation Environment 轨迹 |
| extxyz | .xyz | 扩展 XYZ（主要工作格式） |

## 命令参考

### VASP 转换

```bash
# OUTCAR 转 extxyz（目录，Shell 版本）
gpumdkit.sh -out2xyz <directory>

# OUTCAR 转 extxyz（Python 版本）
gpumdkit.sh -out2exyz <directory>

# XDATCAR 转 extxyz
gpumdkit.sh -xdat2exyz XDATCAR output.xyz

# POSCAR 转 extxyz
gpumdkit.sh -pos2exyz POSCAR model.xyz

# extxyz 转 POSCAR（所有帧）
gpumdkit.sh -exyz2pos structures.xyz
```

### LAMMPS 转换

```bash
# LAMMPS dump 转 extxyz
# 重要：元素符号必须与 dump 文件中的原子类型 ID 匹配
gpumdkit.sh -lmp2exyz dump.lammpstrj Li Y Cl

# POSCAR 转 LAMMPS data
gpumdkit.sh -pos2lmp POSCAR lammps.data
```

### CIF 转换

```bash
# CIF 转 POSCAR
gpumdkit.sh -cif2pos input.cif POSCAR.vasp

# CIF 转 extxyz
gpumdkit.sh -cif2exyz input.cif model.xyz
```

### 其他转换

```bash
# ASE 轨迹转 extxyz
gpumdkit.sh -traj2exyz input.traj output.xyz
```

### 仅菜单或旧版转换

```bash
# MTP cfg 转 extxyz
python3 ${GPUMDkit_path}/Scripts/format_conversion/mtp2xyz.py train.cfg Pd Ag

# 通过 CLI 菜单助手进行 CP2K 转换
gpumdkit.sh -cp2k2xyz

# ABACUS 转换可通过交互菜单使用：
# gpumdkit.sh  # 选择：1) Format Conversion -> 104
```

### dp2xyz

将 DeepMD npy 数据集转换为 extxyz 格式。递归扫描目录中包含 `type.raw`、`type_map.raw` 和 `set.000/` 的数据集。

**用法：**
```
gpumdkit.sh -dp2xyz database train.xyz
```
依赖：`dpdata`、`ase`。
作者：Denan LI (lidenan@westlake.edu.cn)

### 结构操作

```bash
# 当下游 GPUMD 工作流使用分组感知命令时添加分组标签
gpumdkit.sh -addgroup POSCAR Li Y Cl

# 为结构添加权重
gpumdkit.sh -addweight input.xyz output.xyz 2.0

# 复制结构（按倍数）
gpumdkit.sh -replicate POSCAR supercell.vasp 2 2 2

# 复制结构（按目标原子数）
gpumdkit.sh -replicate POSCAR supercell.vasp 256

# 提取特定帧（从 1 开始的索引）
gpumdkit.sh -get_frame trajectory.xyz 1000

# 清理 extxyz 元数据
gpumdkit.sh -clean_xyz input.xyz clean.xyz
```

## 详细 CLI 标志参考

| CLI 标志 | 转换 | 语法 |
|----------|-----------|--------|
| `-out2xyz` | OUTCAR -> extxyz（Shell） | `gpumdkit.sh -out2xyz <dir>` |
| `-out2exyz` | OUTCAR -> extxyz（Python） | `gpumdkit.sh -out2exyz <dir>` |
| `-pos2exyz` | POSCAR -> extxyz | `gpumdkit.sh -pos2exyz <poscar> <xyz>` |
| `-exyz2pos` | extxyz -> POSCAR | `gpumdkit.sh -exyz2pos <xyz>` |
| `-pos2lmp` | POSCAR -> LAMMPS data | `gpumdkit.sh -pos2lmp <poscar> <lmp>` |
| `-lmp2exyz` | LAMMPS dump -> extxyz | `gpumdkit.sh -lmp2exyz <dump> <elem...>` |
| `-cif2pos` | CIF -> POSCAR | `gpumdkit.sh -cif2pos <cif> <output>` |
| `-cif2exyz` | CIF -> extxyz | `gpumdkit.sh -cif2exyz <cif> <output>` |
| `-xdat2exyz` | XDATCAR -> extxyz | `gpumdkit.sh -xdat2exyz XDATCAR dump.xyz` |
| `-traj2exyz` | ASE traj -> extxyz | `gpumdkit.sh -traj2exyz <traj> <xyz>` |
| `-dp2xyz` | DeepMD npy -> extxyz（通过 dpdata） | `gpumdkit.sh -dp2xyz <input_dir/> [output.xyz]` |
| `-addgroup` | 添加分组标签 | `gpumdkit.sh -addgroup <poscar> <elem...>` |
| `-addweight` | 添加权重 | `gpumdkit.sh -addweight <in> <out> <weight>` |
| `-replicate` | 复制结构 | `gpumdkit.sh -replicate <in> <out> a b c` |
| `-get_frame` | 提取帧（从 1 开始的索引） | `gpumdkit.sh -get_frame <xyz> <index>` |
| `-clean_xyz` | 清理 extxyz 信息 | `gpumdkit.sh -clean_xyz <in> <out>` |

## 示例

### 示例 1：转换 VASP MD 输出
```bash
# 转换当前目录中的所有 OUTCAR 文件
gpumdkit.sh -out2xyz .

# 仅在后续工作流需要时添加分组标签
gpumdkit.sh -addgroup POSCAR Pb Ti O

# 结果：model.xyz 可用于 NEP 训练
```

### 示例 2：准备 LAMMPS 模拟
```bash
# 将 POSCAR 转换为 LAMMPS data
gpumdkit.sh -pos2lmp POSCAR system.data

# LAMMPS 模拟后，将 dump 转换回来
gpumdkit.sh -lmp2exyz dump.lammpstrj Li P S
```

### 示例 3：批量转换
```bash
# 转换多个 OUTCAR 文件
for dir in run_*; do
    gpumdkit.sh -out2xyz "$dir"
    mv "$dir"/model.xyz "$dir"/trajectory.xyz
done
```

### 示例 4：结构复制
```bash
# 复制为 2x2x2 超胞
gpumdkit.sh -replicate POSCAR supercell_222.vasp 2 2 2

# 复制到约 256 个原子
gpumdkit.sh -replicate POSCAR supercell_256.vasp 256
```

## 注意事项

1. **extxyz 是主要格式**：大多数 GPUMDkit 工具使用 extxyz 文件
2. **分组标签是可选元数据**：仅在分组感知的 GPUMD 计算或下游工具中添加；化学种类决定原子类型
3. **`-get_frame` 的帧索引从 1 开始**：第一帧索引为 1
4. **LAMMPS 的元素排序很重要**：必须与 dump 文件中的原子类型 ID 匹配
5. **exyz2pos 导出所有帧**：为每帧创建单独的 POSCAR

## 依赖

大多数 Python 脚本需要：
- `ase`（Atomic Simulation Environment）
- `numpy`

## 详细文档

用户指南请参见 `${GPUMDkit_path}/docs/tutorials/en/format_conversion.md` 或中文对应版本。
