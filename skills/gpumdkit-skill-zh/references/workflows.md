# 批处理与主动学习工作流

## 目录

- 可用工作流
- SCF 批处理准备
- MD 采样批处理准备
- 主动学习
- 示例与注意事项

## 可用工作流

| 工作流 | 菜单 | 描述 |
|----------|------|-------------|
| SCF 批处理（VASP） | 301 -> 1 | VASP 单点批处理设置 |
| SCF 批处理（CP2K） | 301 -> 2 | CP2K 单点批处理设置 |
| MD 批处理（GPUMD） | 302 | GPUMD MD 采样批处理设置 |
| MD 批处理（LAMMPS） | 303 | LAMMPS MD 采样批处理设置 |

## SCF 批处理预处理

### VASP SCF 批处理

```bash
# 交互模式
gpumdkit.sh  # 选择：3) Workflow -> 301 -> VASP
```

**前置条件：**
1. 将 POSCAR 文件（或单个 `.xyz`）放在当前目录
2. 准备 `fp/` 目录，包含：
   - `INCAR`
   - `POTCAR`
   - `KPOINTS`

**输出结构：**
```
struct_fp/POSCAR_1.vasp、POSCAR_2.vasp、...
fp/POTCAR、INCAR、KPOINTS          （用户提供）
<prefix>_1/POSCAR -> struct_fp/POSCAR_1.vasp
<prefix>_1/POTCAR -> fp/POTCAR
<prefix>_2/...
presub.sh
```

### CP2K SCF 批处理

```bash
# 直接 Python 执行
python Scripts/workflow/scf_batch_pretreatment_cp2k.py <extxyz_file> <template.inp> <prefix_name>

# 示例
python Scripts/workflow/scf_batch_pretreatment_cp2k.py structures.xyz cp2k_template.inp H2O_batch
```

**前置条件：**
1. 准备 CP2K 输入模板（`cp2k_template.inp`）
2. 模板应从 `pos.xyz` 读取坐标

**输出结构：**
```
<prefix>_1/input.inp、pos.xyz
<prefix>_2/input.inp、pos.xyz
...
```

## MD 采样批处理预处理

### GPUMD MD 批处理

```bash
# 交互模式
gpumdkit.sh  # 选择：3) Workflow -> 302
```

**前置条件：**
1. 将 POSCAR 文件（或单个 `.xyz`）放在当前目录
2. 准备 `md/` 目录，包含：
   - `nep.txt`（NEP 模型）
   - `run_1.in`、`run_2.in`、...（GPUMD 输入文件）

**输出结构：**
```
struct_md/model_1.xyz、model_2.xyz、...
md/nep.txt、run_1.in、run_2.in、...   （用户提供）
sample_1/model.xyz -> struct_md/model_1.xyz
sample_1/nep.txt -> md/nep.txt
sample_1/run.in -> md/run_1.in
presub.sh
```

### LAMMPS MD 批处理

```bash
# 交互模式
gpumdkit.sh  # 选择：3) Workflow -> 303
```

**前置条件：**
1. 将 POSCAR 文件（或单个 `.xyz`）放在当前目录
2. 准备 `md/` 目录，包含：
   - `nep.txt`（NEP 模型）
   - `lmprun.in`（LAMMPS 输入文件）

**输出结构：**
```
struct_md/lammps_1.data、lammps_2.data、...
md/lmprun.in、nep.txt                  （用户提供）
sample_1/lammps.data -> struct_md/lammps_1.data
sample_1/lmprun.in -> md/lmprun.in
sample_1/nep.txt -> md/nep.txt
presub.sh
```

## 主动学习工作流

### 概述

主动学习通过以下步骤迭代改进 NEP 模型：
1. 使用当前 NEP 模型运行 MD
2. 选择多样化结构
3. 计算 DFT 参考数据
4. 重新训练 NEP 模型

### 手动主动学习循环

```bash
# 步骤 1：使用当前 NEP 进行 MD 采样
gpumdkit.sh  # 选择：3) Workflow -> 302

# 步骤 2：过滤和采样结构
gpumdkit.sh -min_dist_pbc dump.xyz
gpumdkit.sh -filter_box dump.xyz 13
gpumdkit.sh  # 选择：2) Sample Structures -> 203

# 步骤 3：准备 SCF 计算
gpumdkit.sh  # 选择：3) Workflow -> 301

# 步骤 4：运行 DFT 计算（外部）

# 步骤 5：转换和验证
gpumdkit.sh -out2xyz ./scf_results/
gpumdkit.sh -range new_data.xyz force

# 步骤 6：重新训练 NEP 模型（外部）
```

### 自动化主动学习（开发中）

```bash
# 脚本可用但仍在开发中：
# Scripts/workflow/workflow_active_learning_dev.sh
# Scripts/workflow/workflow_active_learning_dev_multielement.sh
```

除非用户明确要求实验性主动学习自动化，否则不要运行这些开发脚本。

**可配置参数：**

| 参数 | 描述 | 示例 |
|-----------|-------------|---------|
| `prefix_name` | 目录前缀 | `LiF_iter01` |
| `min_dist` | 最小距离过滤 | `1.4` |
| `box_limit` | 最大盒子尺寸过滤 | `13` |
| `max_fp_num` | 最大 DFT 计算数 | `50` |
| `sample_method` | 采样方法 | `uniform`、`random` 或 `pynep` |
| `pynep_sample_dist` | FPS 距离阈值 | `0.01` |

## 示例

### 示例 1：VASP 高通量筛选
```bash
# 准备结构
for struct in structures/*.vasp; do
    name=$(basename "$struct" .vasp)
    mkdir -p "calc_$name"
    cp "$struct" "calc_$name/POSCAR"
    ln -s ../POTCAR "calc_$name/"
    ln -s ../INCAR "calc_$name/"
    ln -s ../KPOINTS "calc_$name/"
done
```

### 示例 2：温度依赖的 MD
```bash
# 为多个温度设置目录
for temp in 300 500 700 900; do
    mkdir -p "${temp}K"
    cp model.xyz "${temp}K/"
    cp nep.txt "${temp}K/"
    sed "s/TEMPERATURE/$temp/g" run_template.in > "${temp}K/run.in"
done
```

### 示例 3：NEP 训练数据生成
```bash
# 1. 从初始结构开始
gpumdkit.sh -replicate POSCAR supercell.vasp 3 3 3

# 2. 微扰结构
python Scripts/sample_structures/perturb_structure.py supercell.vasp 20 0.03 0.2 uniform

# 3. 运行 MD 采样
gpumdkit.sh  # 选择：3) Workflow -> 302

# 4. 选择多样化结构
gpumdkit.sh  # 选择：2) Sample Structures -> 203

# 5. 准备 DFT 计算
gpumdkit.sh  # 选择：3) Workflow -> 301
```

## 注意事项

1. **目录结构很重要**：工作流期望特定的目录布局
2. **使用符号链接**：输出目录使用符号链接节省磁盘空间
3. **presub.sh**：生成的提交脚本可能需要根据集群进行修改
4. **模板文件**：批处理前确保输入模板正确

## 详细文档

- `${GPUMDkit_path}/docs/tutorials/en/workflow_scripts.md`
- `${GPUMDkit_path}/docs/tutorials/en/active_learning_workflow.md`
