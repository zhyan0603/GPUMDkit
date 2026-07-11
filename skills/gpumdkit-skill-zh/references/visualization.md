# 可视化

## 目录

- 快速参考
- 绘图类别
- 常用工作流
- 输出文件与依赖

## 快速参考

```bash
gpumdkit.sh -plt <type>              # 交互式显示绘图
gpumdkit.sh -plt <type> save         # 常用保存形式；请验证下方绘图条目
gpumdkit.sh -plt -h                  # 列出绘图类型
```

分发器未提供统一的 `-plt <type> -h` 约定，参数位置因绘图脚本而异。请使用本参考中的签名，并在传递额外参数前检查目标脚本。不要假设 `save` 在每个绘图中都以相同位置接受。

## 绘图类别

### NEP 训练与评估（13 种绘图类型）

| 命令 | 输入文件 | 描述 |
|---------|-------------|-------------|
| `train` | `loss.out`、`*_train.out` | 训练损失曲线和对比图 |
| `prediction` / `test` | `*_test.out` | 测试集对比图 |
| `train_test` | `*_train.out`、`*_test.out` | 合并训练/测试对比图 |
| `parity_density` | `*_train.out` | 大数据集的密度对比图 |
| `train_density` | `loss.out`、`*_train.out` | 带密度可视化的训练 |
| `force_errors` | `force_train.out` | 力误差分析 |
| `restart` | `nep.restart` | NEP 重启文件可视化 |
| `charge` | `charge_train.out` | 电荷分布（qNEP） |
| `born_charge` / `bec` | `bec_train.out`、`bec_test.out` | Born 有效电荷 |
| `dimer` | NEP 模型 | 二聚体相互作用曲线 |
| `des` | `descriptors.npy` | 描述符 PCA/UMAP 可视化 |
| `lr` | `loss.out`（gnep） | 学习率衰减 |

```bash
# 训练结果
gpumdkit.sh -plt train
gpumdkit.sh -plt prediction
gpumdkit.sh -plt force_errors

# 描述符可视化（需先计算）
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li
gpumdkit.sh -plt des pca
gpumdkit.sh -plt des umap

# 二聚体绘图
gpumdkit.sh -plt dimer Li Li nep.txt
```

### 输运性质（10 种绘图类型）

| 命令 | 输入文件 | 描述 |
|---------|-------------|-------------|
| `msd` | `msd.out` | 均方位移 |
| `msd_conv` | `msd_step*.out` | MSD 收敛检查 |
| `msd_all` | `msd.out`（all_groups） | 各物种的 MSD |
| `sdc` | `msd.out` | 自扩散系数 |
| `msd_sdc` | `msd.out` | MSD 和 SDC 合并 |
| `sigma` / `arrhenius_sigma` | `*K/` 目录 | Arrhenius 离子电导率 |
| `D` / `arrhenius_d` | `*K/` 目录 | Arrhenius 扩散系数 |
| `sigma_xyz` | `*K/` 目录 | 方向性 Arrhenius 电导率 |
| `D_xyz` | `*K/` 目录 | 方向性 Arrhenius 扩散系数 |
| `doas` | `doas.out` | 原子态密度 |

```bash
# MSD 和扩散
gpumdkit.sh -plt msd
gpumdkit.sh -plt sdc
gpumdkit.sh -plt msd_sdc

# 所有物种的 MSD（需要 run.in 中的 all_groups）
gpumdkit.sh -plt msd_all msd.out Li P S

# Arrhenius 绘图（需要按温度组织的目录）
# 目录结构：300K/、350K/、400K/、...，每个包含 msd.out
gpumdkit.sh -plt arrhenius_sigma
gpumdkit.sh -plt arrhenius_d

# DOAS 可视化（需先计算）
gpumdkit.sh -plt doas doas.out Li
```

### 结构分析（9 种绘图类型）

| 命令 | 输入文件 | 描述 |
|---------|-------------|-------------|
| `thermo` | `thermo.out` | 热力学性质 |
| `thermo2` | `thermo.out` | 替代热力学样式 |
| `thermo3` | `thermo.out` | 第三种热力学样式 |
| `rdf` | `rdf.out` | 径向分布函数 |
| `rdf_pmf` | `rdf.out` | RDF + 平均力势 |
| `vac` | `sdc.out` | 速度自相关 |
| `cohesive` | `cohesive.out` | 内聚能曲线 |
| `net_force` | extxyz 文件 | 净力分布 |
| `plane-grid` | `model.xyz`、`displacements.dat` | 位移网格可视化 |

```bash
# 热力学性质
gpumdkit.sh -plt thermo

# RDF 分析
gpumdkit.sh -plt rdf
gpumdkit.sh -plt rdf 2              # 特定列
gpumdkit.sh -plt rdf_pmf 300        # 300K 下的 PMF

# 平面网格位移
gpumdkit.sh -plt plane-grid -i model.xyz -d displacements.dat -e Pb Sr
```

### 热输运（4 种绘图类型）

| 命令 | 输入文件 | 描述 |
|---------|-------------|-------------|
| `emd` | EMD 输出 | EMD 热导率 |
| `nemd` | NEMD 输出 | NEMD 热输运 |
| `hnemd` | HNEMD 输出 | HNEMD 热输运 |
| `viscosity` | `viscosity.out` | 粘度分量 |

```bash
# EMD 热导率
gpumdkit.sh -plt emd x

# NEMD 热输运
# 参数：real_length scale_eff_size cutoff_freq
gpumdkit.sh -plt nemd <real_length> <scale_eff_size> <cutoff_freq> save

# HNEMD 热输运
gpumdkit.sh -plt hnemd <scale_eff_size> <cutoff_freq> save

# 粘度
gpumdkit.sh -plt viscosity save
```

### 声子（1 种绘图类型）

| 命令 | 输入文件 | 描述 |
|---------|-------------|-------------|
| `pdos` | `model.xyz`、`run.in`、`dos.out`、`mvac.out` | 声子态密度和热容 |

```bash
gpumdkit.sh -plt pdos save
```

## 常用工作流

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
# 3. 合并 MSD-SDC 绘图
gpumdkit.sh -plt msd_sdc
# 4. Arrhenius 分析（多温度）
gpumdkit.sh -plt arrhenius_d
```

### 热输运
```bash
# 1. 绘制热力学性质
gpumdkit.sh -plt thermo
# 2. 绘制热导率
gpumdkit.sh -plt emd x
# 3. 或 NEMD/HNEMD
gpumdkit.sh -plt nemd 10 1 60 save
```

## 输出文件

| 绘图类型 | PNG 文件名（使用 `save`） |
|-----------|---------------------------|
| `train` | `train.png` |
| `prediction` | `prediction.png` |
| `train_test` | `train_test.png` |
| `force_errors` | `force_errors.png` |
| `des` | `descriptors.png` |
| `msd` | `msd.png` |
| `sdc` | `sdc.png` |
| `msd_sdc` | `msd_sdc.png` |
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

并非所有绘图共用一套依赖。多数脚本使用 `matplotlib` 和 `numpy`；单个
绘图还可能需要 `pandas`、`scipy`、`seaborn`、`ase`、`scikit-learn`、
`umap-learn`、`calorine` 或 `ferrodispcalc`。执行前检查所选绘图的导入和
帮助。新建或调整绘图风格时参见 `plotting-style.md`。

## 详细文档

用户指南请参见 `${GPUMDkit_path}/docs/tutorials/en/plot_scripts.md` 或中文对应版本。
