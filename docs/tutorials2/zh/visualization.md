<div align="center">
  <h1>可视化</h1>
  <p>
    <strong>简体中文</strong> | <a href="../en/visualization.md">English</a>
  </p>
</div>

模拟结果的绘图和可视化工具。

## 快速参考

```bash
gpumdkit.sh -plt <类型>              # 显示图表
gpumdkit.sh -plt <类型> save         # 保存为 PNG
gpumdkit.sh -plt <类型> -h           # 获取帮助
```

## NEP 训练图表

| 命令 | 描述 |
|------|------|
| `train` | 训练损失曲线和对比图 |
| `prediction` / `test` | 测试集对比图 |
| `train_test` | 组合训练/测试对比图 |
| `force_errors` | 力误差分析 |
| `des` | 描述符 PCA/UMAP 可视化 |
| `dimer` | 二聚体相互作用曲线 |
| `charge` | 电荷分布 (qNEP) |
| `born_charge` / `bec` | Born 有效电荷 |

```bash
gpumdkit.sh -plt train
gpumdkit.sh -plt prediction
gpumdkit.sh -plt force_errors
gpumdkit.sh -plt des pca
gpumdkit.sh -plt dimer Li Li nep.txt
```

## 传输属性图表

| 命令 | 描述 |
|------|------|
| `msd` | 均方位移 |
| `sdc` | 自扩散系数 |
| `msd_sdc` | MSD 和 SDC 组合 |
| `arrhenius_sigma` | Arrhenius 离子电导率 |
| `arrhenius_d` | Arrhenius 扩散率 |

```bash
gpumdkit.sh -plt msd
gpumdkit.sh -plt sdc
gpumdkit.sh -plt arrhenius_sigma
gpumdkit.sh -plt arrhenius_d
```

## 结构图表

| 命令 | 描述 |
|------|------|
| `thermo` | 热力学属性 |
| `rdf` | 径向分布函数 |
| `rdf_pmf` | RDF + 平均力势 |
| `vac` | 速度自相关 |
| `doas` | 原子态密度 |
| `net_force` | 净力分布 |

```bash
gpumdkit.sh -plt thermo
gpumdkit.sh -plt rdf
gpumdkit.sh -plt doas doas.out Li
gpumdkit.sh -plt net_force train.xyz
```

## 热传输图表

| 命令 | 描述 |
|------|------|
| `emd` | EMD 热导率 |
| `nemd` | NEMD 热传输 |
| `hnemd` | HNEMD 热传输 |
| `viscosity` | 粘度分量 |

```bash
gpumdkit.sh -plt emd x
gpumdkit.sh -plt nemd <real_length> <scale_eff_size> <cutoff_freq> save
gpumdkit.sh -plt hnemd <scale_eff_size> <cutoff_freq> save
```

## 输出文件

| 图表 | PNG 文件名 |
|------|-----------|
| `train` | `train.png` |
| `prediction` | `prediction.png` |
| `msd` | `msd.png` |
| `sdc` | `sdc.png` |
| `thermo` | `thermo.png` |
| `rdf` | `rdf.png` |
| `arrhenius_sigma` | `Arrhenius_sigma.png` |
| `arrhenius_d` | `Arrhenius_D.png` |

## 依赖

```bash
pip install matplotlib numpy
```
