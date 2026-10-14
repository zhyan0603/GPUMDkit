<div align="center">
  <h1>🌡️ MD 与结构</h1>
  <p style="text-align: justify;">MD 热力学输出和结构分析的绘图，包括极性材料的平面分布图。</p>
</div>

← 返回[绘图脚本总览](绘图脚本.md)

## 热力学量

### thermo (`plt_thermo.py`)

综合可视化热力学性质的主脚本。

**输入文件：** `thermo.out`

```bash
gpumdkit.sh -plt thermo
```

<div align="center">
  <img src="../../Gallery/thermo.png" alt="Thermo plot" width="75%" />
</div>

---

### thermo2 / thermo3 (`plt_thermo2.py`, `plt_thermo3.py`)

不同风格的热力学量可视化。

**输入文件：** `thermo.out`

```bash
gpumdkit.sh -plt thermo2
gpumdkit.sh -plt thermo3
```


## 结构分析

### rdf (`plt_rdf.py`)

绘制 `rdf.out` 中的所有 RDF 数值列（第 1 列半径用作横坐标）。

**输入文件：** `rdf.out`

```bash
gpumdkit.sh -plt rdf               # 显示所有 RDF 对
gpumdkit.sh -plt rdf save          # 将图保存为 rdf.png
```

`rdf` 绘图不提供单列选择。若需要选定一列进行 PMF 分析，请使用
`rdf_pmf`。

**RDF 输出：**

<div align="center">
  <img src="../../Gallery/rdf1.png" alt="Complete RDF" width="90%" />
</div>

---

### xrd (`plt_xrd.py`)

绘制计算器 413 生成的 X 射线衍射（XRD）输出。

**输入文件：** 默认读取 `xrd.out`；也可以将其他 XRD 输出文件路径作为第一个参数传入。

```bash
gpumdkit.sh -plt xrd
gpumdkit.sh -plt xrd path/to/xrd.out save
```

使用 `save` 时，图片保存为当前工作目录下的 `xrd.png`。

**示例：**

<div align="center">
  <img src="../../Gallery/xrd.png" alt="XRD 绘图示例" width="75%" />
</div>

---

### xrd_comp (`plt_xrd_comp.py`)

比较当前工作目录下多个温度目录中的 XRD 曲线。每个目录必须命名为
`<temperature>K`，并包含由计算器 413 写出的 `xrd.out`。图像按照高温在上、
低温在下的顺序堆叠显示。

请在包含多个 `*K` 子目录的目录中执行以下命令；该命令不接受目录参数：

```bash
cd /path/to/xrd_series
gpumdkit.sh -plt xrd_comp
gpumdkit.sh -plt xrd_comp save
```

使用 `save` 时，图片保存为当前工作目录下的 `xrd_comp.png`。

---

### rdf_pmf (`plt_rdf_pmf.py`)

绘制 RDF 与平均力势 (PMF) 的组合图。

**输入文件：** `rdf.out`

```bash
gpumdkit.sh -plt rdf_pmf
gpumdkit.sh -plt rdf_pmf 300 2 save
```

两个参数都可以省略：第一个是把 RDF 换算为 PMF 所用的温度（K），第二个 `column_index` 指定使用 `rdf.out` 的哪一列（第 2 列是总 RDF，第 3 列及以后是各原子对的 RDF）。

---

### vac (`plt_vac.py`)

绘制速度自相关函数 (VAC)。可用于分析声子性质和原子动力学。

**输入文件：** GPUMD `compute_sdc` 写入的 `sdc.out`（绘图脚本会从该文件中读取
VAC 相关列）。该文件与 `-plt sdc`、`-plt msd_sdc` 读取的 `msd.out` 不同。

**输出：** 交互式图表或 `vac.png`（使用 `save` 选项时）

```bash
gpumdkit.sh -plt vac
```

---

### cohesive (`plt_cohesive.py`)

从 `cohesive.out` 绘制内聚能曲线。可用于分析晶格稳定性和平衡晶格常数。

**输入文件：** `cohesive.out`（各向同性缩放因子 vs 内聚能）

**输出：** 交互式图表或 `cohesive.png`（使用 `save` 选项时）

```bash
gpumdkit.sh -plt cohesive
```

---

### net_force (`plt_net_force.py`)

绘制结构上的净力分布，可用于识别问题构型。

**输入文件：** `train.xyz`（extxyz 格式）

```bash
gpumdkit.sh -plt net_force train.xyz
```

**参考文献：** [arXiv:2510.19774](https://arxiv.org/abs/2510.19774)

<div align="center">
  <img src="../../Gallery/net_force_distribution.png" alt="Net force distribution" width="58%" />
</div>


## 极性材料的 Plane-Grid 绘图

此工作流将位移或极化数据映射到网格上，并绘制选定平面的轮廓。详细用法和实际案例请参见 **[极性材料分析](极性材料分析.md)**。

**依赖：**

```bash
pip3 install git+https://github.com/MoseyQAQ/ferrodispcalc.git
```

**典型前置步骤：**

```bash
gpumdkit.sh -calc nlist -i model.xyz -c 4 -n 12 -C Pb Sr -E O
gpumdkit.sh -calc disp -i movie.xyz -n nl-Pb_Sr-O.dat -o displacements.dat
gpumdkit.sh -calc avg-struct -i movie.xyz -l 0.2 -o averaged_structure.xyz
```

**使用方法：**

```bash
gpumdkit.sh -plt plane-grid -i averaged_structure.xyz -d displacements.dat -e Pb Sr
gpumdkit.sh -plt plane-grid -i averaged_structure.xyz -d displacements.dat -e Pb Sr --select-xy 0 1
```
