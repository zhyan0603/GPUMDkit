<div align="center">
  <h1>格式转换</h1>
  <p>
    <strong>简体中文</strong> | <a href="../en/format_conversion.md">English</a>
  </p>
</div>

在计算材料科学文件格式之间转换。

## 支持的格式

| 格式 | 扩展名 | 软件 |
|------|--------|------|
| VASP | POSCAR, OUTCAR, XDATCAR, vasprun.xml | VASP |
| LAMMPS | .data, dump.* | LAMMPS |
| CP2K | .log, pos.xyz, frc.xyz, cell.cell | CP2K |
| ABACUS | running_scf.log, running_md.log | ABACUS |
| CIF | .cif | 各种软件 |
| MTP | .cfg | MTP |
| ASE | .traj | ASE |
| extxyz | .xyz | 各种软件（主要格式） |

## 快速参考

```bash
# VASP
gpumdkit.sh -out2xyz <目录>              # OUTCAR -> extxyz
gpumdkit.sh -pos2exyz POSCAR model.xyz   # POSCAR -> extxyz
gpumdkit.sh -exyz2pos structures.xyz     # extxyz -> POSCAR

# LAMMPS
gpumdkit.sh -lmp2exyz dump.lammpstrj Li Y Cl  # LAMMPS -> extxyz
gpumdkit.sh -pos2lmp POSCAR lammps.data        # POSCAR -> LAMMPS

# CIF
gpumdkit.sh -cif2pos input.cif POSCAR.vasp     # CIF -> POSCAR
gpumdkit.sh -cif2exyz input.cif model.xyz      # CIF -> extxyz

# 其他
gpumdkit.sh -xdat2exyz XDATCAR output.xyz      # XDATCAR -> extxyz
gpumdkit.sh -traj2exyz input.traj output.xyz   # ASE 轨迹 -> extxyz
```

## 结构操作

```bash
# 添加组标签（NEP 必需）
gpumdkit.sh -addgroup POSCAR Li Y Cl

# 添加权重
gpumdkit.sh -addweight input.xyz output.xyz 2.0

# 复制结构
gpumdkit.sh -replicate POSCAR supercell.vasp 2 2 2
gpumdkit.sh -replicate POSCAR supercell.vasp 256  # 按原子数

# 提取帧（从 0 开始）
gpumdkit.sh -get_frame trajectory.xyz 1000

# 清理元数据
gpumdkit.sh -clean_xyz input.xyz clean.xyz
```

## 示例

### 转换 VASP 输出

```bash
gpumdkit.sh -out2xyz .
gpumdkit.sh -addgroup POSCAR Pb Ti O
```

### 准备 LAMMPS 模拟

```bash
gpumdkit.sh -pos2lmp POSCAR system.data
# 模拟后
gpumdkit.sh -lmp2exyz dump.lammpstrj Li P S
```

### 批量转换

```bash
for dir in run_*; do
    gpumdkit.sh -out2xyz "$dir"
done
```

## 注意事项

- extxyz 是 GPUMDkit 的主要格式
- GPUMD/NEP 需要组标签
- 帧索引从 0 开始
- LAMMPS 元素顺序必须与原子类型 ID 匹配

## 依赖

```bash
pip install ase numpy
```
