# 格式转换指南

<div align="center">
  <p>
    <strong>中文</strong> | <a href="../en/format_conversion.md">English</a>
  </p>
</div>

本指南介绍 GPUMDkit 中的所有格式转换工具。

## 支持的格式

| 格式 | 扩展名 | 描述 |
|------|--------|------|
| VASP | POSCAR, OUTCAR, XDATCAR, vasprun.xml | 维也纳从头算模拟包 |
| LAMMPS | .data, dump.* | 大规模原子/分子大规模并行模拟器 |
| CP2K | .log, pos.xyz, frc.xyz, cell.cell | 量子化学和固态物理 |
| ABACUS | running_scf.log, running_md.log | 基于原子轨道的从头算计算 |
| CIF | .cif | 晶体学信息文件 |
| MTP | .cfg | 矩量张量势格式 |
| ASE | .traj | 原子模拟环境轨迹 |
| extxyz | .xyz | 扩展 XYZ（主要工作格式） |

## 交互模式

启动交互式菜单并选择选项 1：

```bash
gpumdkit.sh
# 选择: 1) Format Conversion
```

您将看到：

```
+-------------------------------------------------------------+
|                   FORMAT CONVERSION TOOLS                   |
+-------------------------------------------------------------+
| 101) VASP to extxyz            106) Add group labels        |
| 102) MTP to extxyz             107) Add weight to extxyz    |
| 103) CP2K to extxyz            108) Extract frame extxyz    |
| 104) ABACUS to extxyz          109) Clean XYZ info          |
| 105) extxyz to POSCAR          110) Replicate structure     |
+-------------------------------------------------------------+
| out2exyz) OUTCAR to extxyz     xdat2exyz) XDATCAR to extxyz |
| pos2exyz) POSCAR to extxyz     pos2lmp)   POSCAR to LAMMPS  |
| cif2pos)  CIF to POSCAR        lmp2exyz)  LAMMPS to extxyz  |
| cif2exyz) CIF to extxyz        traj2exyz) ASE traj to extxyz|
+-------------------------------------------------------------+
| 000) Return to main menu                                    |
+-------------------------------------------------------------+
```

## 命令行模式

### VASP 转换

#### OUTCAR 转 extxyz

```bash
# Shell 版本（支持 OUTCAR 和 vasprun.xml）
gpumdkit.sh -out2xyz <目录>

# Python 版本（仅 OUTCAR）
gpumdkit.sh -out2exyz <目录>

# 示例
gpumdkit.sh -out2xyz ./vasp_results/
```

#### XDATCAR 转 extxyz

```bash
gpumdkit.sh -xdat2exyz XDATCAR output.xyz
```

#### POSCAR 转 extxyz

```bash
gpumdkit.sh -pos2exyz POSCAR model.xyz
```

#### extxyz 转 POSCAR

```bash
# 将所有帧转换为单独的 POSCAR 文件
gpumdkit.sh -exyz2pos structures.xyz
```

### LAMMPS 转换

#### LAMMPS dump 转 extxyz

```bash
# 重要：元素符号必须与 dump 文件中的原子类型 ID 匹配
gpumdkit.sh -lmp2exyz dump.lammpstrj Li Y Cl

# 示例
gpumdkit.sh -lmp2exyz trajectory.dump Na Cl
```

#### POSCAR 转 LAMMPS data

```bash
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

## 结构操作

### 添加组标签

组标签是 GPUMD/NEP 识别原子类型所必需的：

```bash
gpumdkit.sh -addgroup POSCAR <元素1> <元素2> ...

# 示例
gpumdkit.sh -addgroup POSCAR Li Y Cl
```

这将创建带有组标签的 `model.xyz` 文件。

### 添加权重

为结构添加权重标签用于加权训练：

```bash
gpumdkit.sh -addweight input.xyz output.xyz <权重>

# 示例
gpumdkit.sh -addweight train.xyz weighted.xyz 2.0
```

### 复制结构

从单胞创建超胞：

```bash
# 按复制因子
gpumdkit.sh -replicate POSCAR supercell.vasp 2 2 2

# 按目标原子数
gpumdkit.sh -replicate POSCAR supercell.vasp 256
```

### 提取帧

从轨迹中提取特定帧：

```bash
# 从 0 开始的索引
gpumdkit.sh -get_frame trajectory.xyz 1000
```

### 清理 XYZ 元数据

从 extxyz 文件中删除额外的元数据：

```bash
gpumdkit.sh -clean_xyz input.xyz clean.xyz
```

## 完整 CLI 参考

| 标志 | 描述 | 语法 |
|------|------|------|
| `-out2xyz` | OUTCAR 转 extxyz (shell) | `gpumdkit.sh -out2xyz <dir>` |
| `-out2exyz` | OUTCAR 转 extxyz (python) | `gpumdkit.sh -out2exyz <dir>` |
| `-pos2exyz` | POSCAR 转 extxyz | `gpumdkit.sh -pos2exyz <poscar> <xyz>` |
| `-exyz2pos` | extxyz 转 POSCAR | `gpumdkit.sh -exyz2pos <xyz>` |
| `-pos2lmp` | POSCAR 转 LAMMPS data | `gpumdkit.sh -pos2lmp <poscar> <lmp>` |
| `-lmp2exyz` | LAMMPS dump 转 extxyz | `gpumdkit.sh -lmp2exyz <dump> <elem...>` |
| `-cif2pos` | CIF 转 POSCAR | `gpumdkit.sh -cif2pos <cif> <output>` |
| `-cif2exyz` | CIF 转 extxyz | `gpumdkit.sh -cif2exyz <cif> <output>` |
| `-xdat2exyz` | XDATCAR 转 extxyz | `gpumdkit.sh -xdat2exyz XDATCAR dump.xyz` |
| `-traj2exyz` | ASE 轨迹转 extxyz | `gpumdkit.sh -traj2exyz <traj> <xyz>` |
| `-addgroup` | 添加组标签 | `gpumdkit.sh -addgroup <poscar> <elem...>` |
| `-addweight` | 添加权重 | `gpumdkit.sh -addweight <in> <out> <weight>` |
| `-replicate` | 复制结构 | `gpumdkit.sh -replicate <in> <out> a b c` |
| `-get_frame` | 提取帧 | `gpumdkit.sh -get_frame <xyz> <index>` |
| `-clean_xyz` | 清理 extxyz 信息 | `gpumdkit.sh -clean_xyz <in> <out>` |

## 示例

### 示例 1：转换 VASP MD 输出

```bash
# 转换当前目录中的所有 OUTCAR 文件
gpumdkit.sh -out2xyz .

# 为 NEP 训练添加组标签
gpumdkit.sh -addgroup POSCAR Pb Ti O

# 结果：准备好用于 NEP 训练的 model.xyz
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
# 创建 2x2x2 超胞
gpumdkit.sh -replicate POSCAR supercell_222.vasp 2 2 2

# 创建约 256 个原子的超胞
gpumdkit.sh -replicate POSCAR supercell_256.vasp 256
```

## 重要提示

1. **extxyz 是主要格式**：大多数 GPUMDkit 工具使用 extxyz 文件
2. **组标签是必需的**：GPUMD/NEP 需要组标签来识别原子类型
3. **帧索引从 0 开始**：第一帧的索引为 0
4. **LAMMPS 的元素顺序很重要**：必须与 dump 文件中的原子类型 ID 匹配
5. **exyz2pos 导出所有帧**：为每帧创建单独的 POSCAR 文件

## 依赖

大多数 Python 脚本需要：
- `ase`（原子模拟环境）
- `numpy`

## 故障排除

### 问题：缺少元素信息
**解决方案**：确保 POSCAR 头部有正确的元素名称

### 问题：原子数不正确
**解决方案**：检查 POSCAR 是否启用了选择性动力学

### 问题：LAMMPS 转换物种错误
**解决方案**：验证元素顺序是否与 dump 文件中的原子类型 ID 匹配

## 参见

- [计算器](calculators.md) - 计算材料属性
- [分析工具](analyzers.md) - 结构分析
- [NEP 训练指南](nep_training.md) - 完整的训练工作流
