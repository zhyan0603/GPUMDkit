# 可视化指南

<div align="center">
  <p>
    <strong>中文</strong> | <a href="../en/visualization.md">English</a>
  </p>
</div>

本指南介绍 GPUMDkit 中的所有绘图和可视化工具。

## 快速参考

```bash
gpumdkit.sh -plt <类型>              # 交互式显示图表
gpumdkit.sh -plt <类型> save         # 保存为 PNG
gpumdkit.sh -plt <类型> -h           # 获取特定图表的帮助
```

## 图表类别

### NEP 训练和评估（12 种图表类型）

| 命令 | 输入文件 | 描述 |
|------|----------|------|
| `train` | `loss.out`, `*_train.out` | 训练损失曲线和对比图 |
| `prediction` / `test` | `*_test.out` | 测试集对比图 |
| `train_test` | `*_train.out`, `*_test.out` | 组合训练/测试对比图 |
| `parity_density` | `*_train.out` | 基于密度的对比图（大数据集） |
| `train_density` | `loss.out`, `*_train.out` | 带密度可视化的训练 |
| `force_errors` | `force_train.out` | 力误差分析 |
| `restart` | `nep.restart` | NEP 重启文件可视化 |
| `charge` | `charge_train.out` | 电荷分布 (qNEP) |
| `born_charge` / `bec` | `bec_train.out`, `bec_test.out` | Born 有效电荷 |
| `dimer` | NEP 模型 | 二聚体相互作用曲线 |
| `des` | `descriptors.npy` | 描述符 PCA/UMAP 可视化 |
| `lr` | `loss.out` (gnep) | 学习率衰减 |

```bash
# 训练结果
gpumdkit.sh -plt train
gpumdkit.sh -plt prediction
gpumdkit.sh -plt force_errors

# 描述符可视化（需要预先计算）
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li
gpumdkit.sh -plt des pca
gpumdkit.sh -plt des umap

# 二聚体图
gpumdkit.sh -plt dimer Li Li nep.txt
```

### 传输属性（9 种图表类型）

| 命令 | 输入文件 | 描述 |
|------|----------|------|
| `msd` | `msd.out` | 均方位移 |
| `msd_conv` | `msd_step*.out` | MSD 收敛检查 |
| `msd_all` | `msd.out` (all_groups) | 每个物种的 MSD |
| `sdc` | `msd.out` | 自扩散系数 |
| `msd_sdc` | `msd.out` | MSD 和 SDC 组合 |
| `sigma` / `arrhenius_sigma` | `*K/` 目录 | Arrhenius 离子电导率 |
| `D` / `arrhenius_d` | `*K/` 目录 | Arrhenius 扩散率 |
| `sigma_xyz` | `*K/` 目录 | 方向性 Arrhenius 电导率 |
| `D_xyz` | `*K/` 目录 | 方向性 Arrhenius 扩散率 |

```bash
# MSD 和扩散
gpumdkit.sh -plt msd
gpumdkit.sh -plt sdc
gpumdkit.sh -plt msd_sdc

# 所有物种的 MSD（需要在 run.in 中设置 all_groups）
gpumdkit.sh -plt msd_all Li P S

# Arrhenius 图（需要按温度组织的目录）
# 目录结构：300K/, 350K/, 400K/, ... 每个包含 msd.out
gpumdkit.sh -plt arrhenius_sigma
gpumdkit.sh -plt arrhenius_d
```

### 结构分析（10 种图表类型）

| 命令 | 输入文件 | 描述 |
|------|----------|------|
| `thermo` | `thermo.out` | 热力学属性 |
| `thermo2` | `thermo.out` | 替代热力学样式 |
| `thermo3` | `thermo.out` | 第三种热力学样式 |
| `rdf` | `rdf.out` | 径向分布函数 |
| `rdf_pmf` | `rdf.out` | RDF + 平均力势 |
| `vac` | `sdc.out` | 速度自相关 |
| `cohesive` | `cohesive.out` | 内聚能曲线 |
| `net_force` | extxyz 文件 | 净力分布 |
| `doas` | `doas.out` | 原子态密度 |
| `plane-grid` | `model.xyz`, `displacements.dat` | 位移网格可视化 |

```bash
# 热力学属性
gpumdkit.sh -plt thermo

# RDF 分析
gpumdkit.sh -plt rdf
gpumdkit.sh -plt rdf 2              # 特定列
gpumdkit.sh -plt rdf_pmf 300        # 300K 时的 PMF

# DOAS 可视化（需要预先计算）
gpumdkit.sh -plt doas doas.out Li

# 平面网格位移
gpumdkit.sh -plt plane-grid -i model.xyz -d displacements.dat -e Pb Sr
```

### 热传输（4 种图表类型）

| 命令 | 输入文件 | 描述 |
|------|----------|------|
| `emd` | EMD 输出 | EMD 热导率 |
| `nemd` | NEMD 输出 | NEMD 热传输 |
| `hnemd` | HNEMD 输出 | HNEMD 热传输 |
| `viscosity` | `viscosity.out` | 粘度分量 |

```bash
# EMD 热导率
gpumdkit.sh -plt emd x

# NEMD 热传输
# 参数：real_length scale_eff_size cutoff_freq
gpumdkit.sh -plt nemd <real_length> <scale_eff_size> <cutoff_freq> save

# HNEMD 热传输
gpumdkit.sh -plt hnemd <scale_eff_size> <cutoff_freq> save

# 粘度
gpumdkit.sh -plt viscosity save
```

### 声子（1 种图表类型）

| 命令 | 输入文件 | 描述 |
|------|----------|------|
| `pdos` | `model.xyz`, `run.in`, `dos.out`, `mvac.out` | 声子态密度和热容 |

```bash
gpumdkit.sh -plt pdos save
```

## 常见工作流程

### NEP 训练验证

```bash
# 1. 绘制训练损失
gpumdkit.sh -plt train

# 2. 检查测试预测
gpumdkit.sh -plt prediction

# 3. 分析力误差
gpumdkit.sh -plt force_errors

# 4. 可视化描述符
gpumdkit.sh -plt des pca
```

### 扩散分析

```bash
# 1. 绘制 MSD
gpumdkit.sh -plt msd

# 2. 绘制自扩散系数
gpumdkit.sh -plt sdc

# 3. 组合 MSD-SDC 图
gpumdkit.sh -plt msd_sdc

# 4. Arrhenius 分析（多温度）
gpumdkit.sh -plt arrhenius_d
```

### 热传输

```bash
# 1. 绘制热力学属性
gpumdkit.sh -plt thermo

# 2. 绘制热导率
gpumdkit.sh -plt emd x

# 3. 或 NEMD/HNEMD
gpumdkit.sh -plt nemd 10 1 60 save
```

## 输出文件

| 图表类型 | PNG 文件名（使用 `save`） |
|----------|---------------------------|
| `train` | `train.png` |
| `prediction` | `prediction.png` |
| `train_test` | `train_test.png` |
| `force_errors` | `force_errors.png` |
| `des` | `descriptors.png` |
| `msd` | `msd.png` |
| `sdc` | `sdc.png` |
| `msd_sdc` | `MSD_SDC.png` |
| `thermo` | `thermo.png` |
| `rdf` | `rdf.png` |
| `arrhenius_sigma` | `Arrhenius_sigma.png` |
| `arrhenius_d` | `Arrhenius_D.png` |
| `emd` | `emd.png` |
| `nemd` | `nemd.png` |
| `hnemd` | `hnemd.png` |
| `viscosity` | `viscosity.png` |
| `cohesive` | `Cohesive.png` |

## 依赖

所有绘图脚本需要：
- `matplotlib`
- `numpy`

## 参见

- [计算器](calculators.md) - 计算要绘制的属性
- [分析工具](analyzers.md) - 在绘图前分析数据
- [NEP 训练指南](nep_training.md) - 完整的训练工作流
