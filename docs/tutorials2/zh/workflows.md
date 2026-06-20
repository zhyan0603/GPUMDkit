# 工作流指南

<div align="center">
  <p>
    <strong>中文</strong> | <a href="../en/workflows.md">English</a>
  </p>
</div>

本指南介绍 GPUMDkit 中的批处理和工作流自动化工具。

## 可用工作流

| 工作流 | 菜单 | 描述 |
|--------|------|------|
| SCF 批处理 (VASP) | 301.1 | VASP 单点批处理设置 |
| SCF 批处理 (CP2K) | 301.2 | CP2K 单点批处理设置 |
| MD 批处理 (GPUMD) | 302 | GPUMD MD 采样批处理设置 |
| MD 批处理 (LAMMPS) | 303 | LAMMPS MD 采样批处理设置 |

## 交互模式

```bash
gpumdkit.sh
# 选择: 3) Workflow
```

您将看到：

```
+---------------------------------------------------------+
|                      WORKFLOW TOOLS                     |
+---------------------------------------------------------+
| 301) SCF batch pretreatment                             |
| 302) MD sample batch pretreatment (gpumd)               |
| 303) MD sample batch pretreatment (lmp)                 |
+---------------------------------------------------------+
```

## SCF 批处理

### VASP SCF 批处理

为 VASP 准备批量单点能量/力计算。

```bash
gpumdkit.sh  # 选择: 3) Workflow -> 301 -> VASP
```

**前提条件：**

1. 将 POSCAR 文件（或单个 `.xyz`）放在当前目录
2. 准备 `fp/` 目录，包含：
   - `INCAR`
   - `POTCAR`
   - `KPOINTS`

**输出结构：**

```
struct_fp/
├── POSCAR_1.vasp
├── POSCAR_2.vasp
└── ...
fp/
├── INCAR
├── POTCAR
└── KPOINTS
<prefix>_1/
├── POSCAR -> ../struct_fp/POSCAR_1.vasp
├── POTCAR -> ../fp/POTCAR
├── INCAR -> ../fp/INCAR
└── KPOINTS -> ../fp/KPOINTS
<prefix>_2/
└── ...
presub.sh
```

### CP2K SCF 批处理

为 CP2K 准备批量单点计算。

```bash
python Scripts/workflow/scf_batch_pretreatment_cp2k.py <extxyz_file> <template.inp> <prefix_name>

# 示例
python Scripts/workflow/scf_batch_pretreatment_cp2k.py structures.xyz cp2k_template.inp H2O_batch
```

**前提条件：**

1. 准备 CP2K 输入模板（`cp2k_template.inp`）
2. 模板应从 `pos.xyz` 读取坐标

**输出结构：**

```
<prefix>_1/
├── input.inp
└── pos.xyz
<prefix>_2/
└── ...
```

## MD 采样批处理

### GPUMD MD 批处理

为 GPUMD 设置批量 MD 模拟。

```bash
gpumdkit.sh  # 选择: 3) Workflow -> 302
```

**前提条件：**

1. 将 POSCAR 文件（或单个 `.xyz`）放在当前目录
2. 准备 `md/` 目录，包含：
   - `nep.txt`（NEP 模型）
   - `run_1.in`, `run_2.in`, ...（GPUMD 输入文件）

**输出结构：**

```
struct_md/
├── model_1.xyz
├── model_2.xyz
└── ...
md/
├── nep.txt
├── run_1.in
├── run_2.in
└── ...
sample_1/
├── model.xyz -> ../struct_md/model_1.xyz
├── nep.txt -> ../md/nep.txt
└── run.in -> ../md/run_1.in
sample_2/
└── ...
presub.sh
```

### LAMMPS MD 批处理

为 LAMMPS 设置批量 MD 模拟。

```bash
gpumdkit.sh  # 选择: 3) Workflow -> 303
```

**前提条件：**

1. 将 POSCAR 文件（或单个 `.xyz`）放在当前目录
2. 准备 `md/` 目录，包含：
   - `nep.txt`（NEP 模型）
   - `lmprun.in`（LAMMPS 输入文件）

**输出结构：**

```
struct_md/
├── lammps_1.data
├── lammps_2.data
└── ...
md/
├── nep.txt
└── lmprun.in
sample_1/
├── lammps.data -> ../struct_md/lammps_1.data
├── nep.txt -> ../md/nep.txt
└── lmprun.in -> ../md/lmprun.in
sample_2/
└── ...
presub.sh
```

## 主动学习工作流

### 概述

主动学习通过以下步骤迭代改进 NEP 模型：

1. 使用当前 NEP 模型运行 MD
2. 选择多样化的结构
3. 计算 DFT 参考数据
4. 重新训练 NEP 模型

### 手动主动学习循环

```bash
# 步骤 1：使用当前 NEP 进行 MD 采样
gpumdkit.sh  # 选择: 3) Workflow -> 302

# 步骤 2：过滤和采样结构
gpumdkit.sh -min_dist_pbc dump.xyz
gpumdkit.sh -filter_box dump.xyz 13
gpumdkit.sh  # 选择: 2) Sample Structures -> 203

# 步骤 3：准备 SCF 计算
gpumdkit.sh  # 选择: 3) Workflow -> 301

# 步骤 4：运行 DFT 计算（外部）

# 步骤 5：转换和验证
gpumdkit.sh -out2xyz ./scf_results/
gpumdkit.sh -range new_data.xyz force

# 步骤 6：重新训练 NEP 模型（外部）
```

### 可配置参数

| 参数 | 描述 | 示例 |
|------|------|------|
| `prefix_name` | 目录前缀 | `LiF_iter01` |
| `min_dist` | 最小距离过滤 | `1.4` |
| `box_limit` | 最大盒子大小过滤 | `13` |
| `max_fp_num` | 最大 DFT 计算数 | `50` |
| `sample_method` | 采样方法 | `uniform`, `random` 或 `pynep` |
| `pynep_sample_dist` | FPS 距离阈值 | `0.01` |

## 示例

### 示例 1：温度依赖的 MD

```bash
# 为多个温度设置目录
for temp in 300 500 700 900; do
    mkdir -p "${temp}K"
    cp model.xyz "${temp}K/"
    cp nep.txt "${temp}K/"
    sed "s/TEMPERATURE/$temp/g" run_template.in > "${temp}K/run.in"
done
```

### 示例 2：高通量筛选

```bash
# 为多个计算准备结构
for struct in structures/*.vasp; do
    name=$(basename "$struct" .vasp)
    mkdir -p "calc_$name"
    cp "$struct" "calc_$name/POSCAR"
    ln -s ../POTCAR "calc_$name/"
    ln -s ../INCAR "calc_$name/"
    ln -s ../KPOINTS "calc_$name/"
done
```

### 示例 3：NEP 训练数据生成

```bash
# 1. 从初始结构开始
gpumdkit.sh -replicate POSCAR supercell.vasp 3 3 3

# 2. 扰动结构
python Scripts/sample_structures/perturb_structure.py supercell.vasp 20 0.03 0.2 uniform

# 3. 运行 MD 采样
gpumdkit.sh  # 选择: 3) Workflow -> 302

# 4. 选择多样化结构
gpumdkit.sh  # 选择: 2) Sample Structures -> 203

# 5. 准备 DFT 计算
gpumdkit.sh  # 选择: 3) Workflow -> 301
```

## 重要提示

1. **目录结构很重要**：工作流期望特定的目录布局
2. **使用符号链接**：输出目录使用符号链接以节省磁盘空间
3. **presub.sh**：生成的提交脚本可能需要为您的集群进行修改
4. **模板文件**：在批处理之前确保输入模板正确
5. **备份重要文件**：在运行工作流之前始终备份数据

## 参见

- [格式转换](format_conversion.md) - 转换文件格式
- [结构采样](sampling.md) - 为训练采样结构
- [NEP 训练指南](nep_training.md) - 完整的训练工作流
