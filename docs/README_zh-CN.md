<p align="center">
  <img src="./Gallery/gpumdkit_logo_lateral.png" width="55%" alt="GPUMDkit Logo">
</p>
<p align="center">
  <a href="../README.md">English</a>
  &nbsp;·&nbsp;
  <strong>简体中文</strong>
  &nbsp;·&nbsp;
  <a href="https://zhyan0603.github.io/GPUMDkit/">Website</a> &nbsp;·&nbsp;
  <a href="https://zhyan0603.github.io/GPUMDkit/htmls/index.html">Documentation</a>
  &nbsp;·&nbsp;
  <a href="https://zhyan0603.github.io/GPUMDkit/gallery.html">Gallery</a>
  &nbsp;
</p>
<p align="center">
  <a href="https://github.com/zhyan0603/GPUMDkit/releases"><img src="https://img.shields.io/github/v/tag/zhyan0603/GPUMDkit?label=version&style=flat-square&color=brightgreen" alt="Version"></a>
  <a href="https://github.com/zhyan0603/GPUMDkit/blob/main/LICENCE"><img src="https://img.shields.io/badge/license-GPL--3.0-blue" alt="License"></a>
  <a href="https://github.com/zhyan0603/GPUMDkit/stargazers"><img src="https://img.shields.io/github/stars/zhyan0603/GPUMDkit?style=social" alt="Stars"></a>
  <img src="https://img.shields.io/github/languages/code-size/zhyan0603/GPUMDkit" alt="Code Size">
  <a href="https://github.com/zhyan0603/GPUMDkit/graphs/contributors"><img src="https://img.shields.io/github/contributors/zhyan0603/GPUMDkit?style=flat-square&color=brightgreen" alt="Contributors"></a>
</p>
<p style="text-align: justify;"><strong>GPUMDkit</strong> 是面向 GPUMD（<em>Graphics Processing Units Molecular Dynamics</em>）和 NEP（<em>neuroevolution potential</em>）程序的工具包。它提供用户友好的命令行界面，简化常见脚本和工作流程，涵盖脚本调用、格式转换、结构采样、NEP 构建流程及各类分析，旨在提升用户工作效率。</p>


## 功能特点
- **简化脚本调用**：轻松运行 GPUMD 和 NEP 相关脚本。
- **工作流自动化**：自动化常见任务，节省时间，减少人工干预。
- **用户友好界面**：直观的 shell 命令，提升使用体验。

## 安装
按以下步骤安装 `GPUMDkit`：

1. 克隆仓库或下载整个项目。

    ```
    git clone https://github.com/zhyan0603/GPUMDkit.git
    ```

    如需下载指定分支，可使用 `-b` 参数，例如：

    ```
    git clone -b dev https://github.com/zhyan0603/GPUMDkit.git
    ```

2. 执行以下命令：

    ```
    cd GPUMDkit; source ./install.sh
    ```


## 依赖项

`GPUMDkit` 的部分高级功能需要以下 Python 包：

```bash
# 创建干净的 conda 环境
conda create -n gpumdkit python=3.12
conda activate gpumdkit

# 安装所需包
pip install neptrain ase pymatgen dpdata
```

提示：使用 `GPUMDkit` 功能前，请确保已激活 `gpumdkit` 环境。

## 更新

如果设备可以访问 `github`，直接运行：

```
gpumdkit.sh -update
```

否则需要手动下载新版本：

```
wget https://github.com/zhyan0603/GPUMDkit/archive/refs/heads/main.zip
```

## 使用方法

提供两种模式：<u>*交互模式*</u> 和 <u>*命令行模式*</u>

#### 交互模式

---

1. 打开终端。

2. 执行 `gpumdkit.sh`：

   ```
   gpumdkit.sh
   ```

3. 根据屏幕提示交互式选择并运行所需功能。

    ```
               ____ ____  _   _ __  __ ____  _    _ _
              / ___|  _ \| | | |  \/  |  _ \| | _(_) |_
             | |  _| |_) | | | | |\/| | | | | |/ / | __|
             | |_| |  __/| |_| | |  | | |_| |   <| | |_
              \____|_|    \___/|_|  |_|____/|_|\_\_|\__|
    
              GPUMDkit Version 1.5.6 (dev) (2026-06-17)
        Core Developer: Zihan YAN (yanzihan@westlake.edu.cn)
     Main Contributors: Denan LI, Xin WU, Zhoulin LIU & Chen HUA
    
      ---------------------- GPUMD ------------------------
      1) Format Conversion          2) Sample Structures
      3) Workflow                   4) Calculators
      5) Analyzer                   6) Visualization
      7) Utilities                  8) Developing...
      0) Exit
      ------------>>
      Input the function number:
    ```

#### 命令行模式

----

对于熟悉 `GPUMDkit` 的用户，命令行模式可直接向 `gpumdkit.sh` 传递参数，执行更快捷。以下是一些示例：

##### 示例 1：查看帮助信息

```
gpumdkit.sh -h
```

帮助信息如下：

```
+-------------------------------------------------------------------------------------------------------+
|                          GPUMDkit 1.5.6 (dev) (2026-06-17) Command Help                               |
+-------------------------------------------------------------------------------------------------------+
|                                          MAIN FUNCTIONS                                               |
+-------------------------------------------------------------------------------------------------------+
| -h            Show this help table            | -plt <type>        Plot and visualization tools       |
| -calc <type>  Calculator tools                | -time <gpumd|nep>  Time-consuming analyzer            |
| -update       Update GPUMDkit                 | -clean             Clean extra files in current dir   |
+-------------------------------------------------------------------------------------------------------+
|                                         FORMAT CONVERSION                                             |
+-------------------------------------------------------------------------------------------------------+
| -out2xyz      OUTCAR -> extxyz (shell)        | -out2exyz          OUTCAR -> extxyz (python)          |
| -cp2k2xyz     CP2K log -> xyz                 | -xdat2exyz         XDATCAR -> extxyz                  |
| -cif2pos      cif -> POSCAR                   | -cif2exyz          cif -> extxyz                      |
| -pos2exyz     POSCAR -> extxyz                | -exyz2pos          extxyz -> POSCAR                   |
| -pos2lmp      POSCAR -> LAMMPS data           | -lmp2exyz          LAMMPS dump -> extxyz              |
| -traj2exyz    ASE traj -> extxyz              | -replicate         Replicate structure                |
| -addgroup     Add group labels                | -addweight         Add structure weight in extxyz     |
| -clean_xyz    Clean extra info in extxyz      | -get_frame         Extract specific frame             |
| -frame_range  Extract frames by range         |                                                       |
+-------------------------------------------------------------------------------------------------------+
|                                            ANALYSIS                                                   |
+-------------------------------------------------------------------------------------------------------+
| -range        Energy/force/virial statistics  | -analyze_comp      Analyze composition                |
| -chem_species Analyze chemical species        | -cbc               Charge balance check               |
| -min_dist     Min distance (no PBC)           | -min_dist_pbc      Min distance with PBC              |
| -filter_dist  Filter by min_dist (no PBC)     | -filter_dist_pbc   Filter by min_dist (PBC)           |
| -pda          Probability density analysis    | -filter_box        Filter by box-edge length          |
| -pynep        Deprecated PyNEP sampling       |                                                       |
+-------------------------------------------------------------------------------------------------------+
| Detailed usage: gpumdkit.sh -<option> -h    Plot details: gpumdkit.sh -plt <type> -h                  |
+-------------------------------------------------------------------------------------------------------+
```

##### 示例 2：查看 -plt 帮助信息

```
gpumdkit.sh -plt -h
```

帮助信息如下：

```
 +-----------------------------------------------------------------------------------------------+
 |                     GPUMDkit 1.5.6 (dev) (2026-06-17) PLOT & VISUALIZATION TOOLS              |
 +-----------------------------------------------------------------------------------------------+
 |  Usage: gpumdkit.sh -plt <type>                        Help: gpumdkit.sh -plt <type> -h       |
 +-----------------------------------------------------------------------------------------------+
 |                                    NEP Training & Evaluation                                  |
 +-----------------------------------------------------------------------------------------------+
 |  train          - NEP training results           prediction     - NEP prediction results      |
 |  train_test     - NEP train and test results     parity_density - Parity density plot         |
 |  train_density  - Training results density plot  restart        - Parameters in nep.restart   |
 |  charge         - Charge distribution            born_charge    - Born effective charges      |
 |  dimer          - Dimer energy/force curve       force_errors   - Force errors                |
 |  des            - Descriptors                    lr             - Learning rate for gnep      |
 +-----------------------------------------------------------------------------------------------+
 |                                     Diffusion & Transport                                     |
 +-----------------------------------------------------------------------------------------------+
 |  msd            - Mean square displacement       msd_conv       - MSD convergence             |
 |  msd_all        - MSD for all species            sdc            - Self diffusion coefficient  |
 |  msd_sdc        - MSD and SDC together           sigma          - Arrhenius ionic conductivity|
 |  D              - Arrhenius diffusivity          sigma_xyz      - Directional Arrhenius sigma |
 |  D_xyz          - Directional Arrhenius D                                                     |
 +-----------------------------------------------------------------------------------------------+
 |                                    MD & Structural Analysis                                   |
 +-----------------------------------------------------------------------------------------------+
 |  thermo         - thermo info in thermo.out      thermo2/3      - Thermo in different styles  |
 |  rdf            - Radial distribution function   rdf_pmf        - Potential of mean force     |
 |  vac            - Velocity autocorrelation       cohesive       - Cohesive energy curve       |
 |  net_force      - Net force distribution         plane-grid     - Displacement plane grid     |
 |  doas           - Density of atomistic states                                                 |
 +-----------------------------------------------------------------------------------------------+
 |                                        Heat Transport                                         |
 +-----------------------------------------------------------------------------------------------+
 |  emd            - EMD results                    nemd           - NEMD results                |
 |  hnemd          - HNEMD results                  viscosity      - Viscosity                   |
 +-----------------------------------------------------------------------------------------------+
 |                                          Phonons                                              |
 +-----------------------------------------------------------------------------------------------+
 |  pdos           - VAC and PDOS                                                                |
 +-----------------------------------------------------------------------------------------------+
```

##### 示例 3：转换 VASP OUTCAR 为 extxyz

将 `VASP` `OUTCAR` 文件转换为 extxyz 格式：

```
gpumdkit.sh -out2xyz <OUTCAR所在目录>

示例: gpumdkit.sh -out2xyz .
```

##### 示例 4：绘制损失曲线和 parity 图

可视化训练过程中各项损失函数的变化和 parity 图：

```
gpumdkit.sh -plt train
```

<div align="center">
    <img src="./Gallery/train.png" alt="msd" width="75%" />
</div>

##### 示例 5：绘制 parity 图

```
gpumdkit.sh -plt test
```

<div align="center">
    <img src="./Gallery/prediction.png" alt="msd" width="95%" />
</div>

##### 示例 6：绘制 thermo 演化图

可视化 `thermo.out` 中的热力学量演化：

```
gpumdkit.sh -plt thermo
```

![](./Gallery/thermo.png)

如果当前设备不支持显示图形，也可以将图片保存为 PNG：

```
gpumdkit.sh -plt thermo save
```

更多详细示例和命令选项请参考我们的[文档](https://zhyan0603.github.io/GPUMDkit/)。

#### 自定义命令

`GPUMDkit` 支持通过 `~/.gpumdkit.in` 文件自定义命令。

你可以在该文件中定义函数来添加自己的快捷命令（例如 `gpumdkit.sh -yourcommand`），从而扩展 `GPUMDkit` 的功能。详细用法请参见[自定义命令文档](https://zhyan0603.github.io/GPUMDkit/)。

#### Tab 补全支持

`gpumdkit.sh` 提供了可选的 Bash `Tab` 补全功能，增强命令行使用体验。按 `Tab` 键即可自动补全主选项（如 `-h`、`-plt`、`-calc`）及其二级参数（如 `thermo`、`train`）。

##### 使用示例

- 输入 `gpumdkit.sh -<Tab>` 查看所有可用选项。
- 输入 `gpumdkit.sh -plt <Tab>` 列出绘图子选项，如 `thermo`、`train` 等。
- 输入 `gpumdkit.sh -time <Tab>` 查看计算器选项，如 `gpumd`、`nep`。

## 加入我们

欢迎为 **GPUMDkit** 贡献力量！参与方式：

- 通过 [Pull Requests](https://github.com/zhyan0603/GPUMDkit/pulls) 贡献 Python/Shell 脚本。
- 通过 [issues](https://github.com/zhyan0603/GPUMDkit/issues) 报告问题或提出功能建议。
- 通过邮箱 [yanzihan@westlake.edu.cn](mailto:yanzihan@westlake.edu.cn) 联系我。

也欢迎加入我们的 QQ 群（[825696376](https://qun.qq.com/universal-share/share?ac=1&authKey=buBNi1ADDzIFF2oZ1yA5FywG3LA9EL9yKZmb%2BN2MMz7nNuuxTas54wH7BgPEqP0s&busi_data=eyJncm91cENvZGUiOiI4MjU2OTYzNzYiLCJ0b2tlbiI6IlRxL1RLTDlOK3U2ekRSUXJ1TkNTUWd3ODNVV3BrdG9HN2lWWmJKMHAraGlDNzBZWFFyRUY2dUlSaW8rbUd4MisiLCJ1aW4iOiIxNDg5NjQ3MTc5In0%3D&data=fa4zSsT_IdI4ftCT_wwpytYHf--TaTB35lH0Jac5JHVpYoyXw3_3bZ1l1NZejsOZnGJku5u3BCbf5_bgrCkhZg&svctype=4&tempid=h5_group_info)）。一起创造有用的工具！🌟

## 引用

**GPUMDkit** 是一款面向所有人的开源工具。如果在你的研究或工作中它有所帮助，欢迎 ⭐ [在 GitHub 上给我们点亮 Star](https://github.com/zhyan0603/GPUMDkit)。此外，如果 GPUMDkit 对你的发表工作有贡献，请引用我们的论文：

> Z. Yan\*, D. Li, X. Wu, Z. Liu, C. Hua, B. Situ, H. Yang, S. Tang, B. Tang, Z. Wang, S. Yi, H. Wang, D. Huang, K. Li, Q. Guo, Z. Chen, K. Xu, Y. Wang, Z. Wang, G. Tang, S. Liu, Z. Fan, and Y. Zhu\*. **GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP**. [MGE Advances, 2026, e70074](https://doi.org/10.1002/mgea.70074).
