<div align="center">
  <h1>工作流</h1>
  <p>
    <strong>简体中文</strong> | <a href="../en/workflows.md">English</a>
  </p>
</div>

批处理和工作流自动化工具。

## 可用工作流

| 工作流 | 菜单 | 描述 |
|--------|------|------|
| SCF 批处理 (VASP) | 301.1 | VASP 单点批处理设置 |
| SCF 批处理 (CP2K) | 301.2 | CP2K 单点批处理设置 |
| MD 批处理 (GPUMD) | 302 | GPUMD MD 采样批处理设置 |
| MD 批处理 (LAMMPS) | 303 | LAMMPS MD 采样批处理设置 |

## SCF 批处理

### VASP

```bash
gpumdkit.sh  # 选择: 3) Workflow -> 301 -> VASP
```

前提条件：
1. POSCAR 文件在当前目录
2. `fp/` 目录包含 INCAR, POTCAR, KPOINTS

输出：`<prefix>_1/`, `<prefix>_2/`, ..., `presub.sh`

### CP2K

```bash
python Scripts/workflow/scf_batch_pretreatment_cp2k.py <extxyz> <template.inp> <prefix>
```

## MD 批处理

### GPUMD

```bash
gpumdkit.sh  # 选择: 3) Workflow -> 302
```

前提条件：
1. POSCAR 文件在当前目录
2. `md/` 目录包含 `nep.txt`, `run_*.in`

输出：`sample_1/`, `sample_2/`, ..., `presub.sh`

### LAMMPS

```bash
gpumdkit.sh  # 选择: 3) Workflow -> 303
```

前提条件：
1. POSCAR 文件在当前目录
2. `md/` 目录包含 `nep.txt`, `lmprun.in`

## 主动学习

迭代改进 NEP 模型：

```bash
# 1. MD 采样
gpumdkit.sh  # 选择: 3) Workflow -> 302

# 2. 过滤和采样
gpumdkit.sh -min_dist_pbc dump.xyz
gpumdkit.sh -filter_box dump.xyz 13

# 3. 准备 SCF
gpumdkit.sh  # 选择: 3) Workflow -> 301

# 4. 运行 DFT（外部）

# 5. 重新训练 NEP
```

## 示例

### 温度依赖的 MD

```bash
for temp in 300 500 700 900; do
    mkdir -p "${temp}K"
    cp model.xyz "${temp}K/"
    cp nep.txt "${temp}K/"
    sed "s/TEMPERATURE/$temp/g" run_template.in > "${temp}K/run.in"
done
```

### 批量结构准备

```bash
for struct in structures/*.vasp; do
    name=$(basename "$struct" .vasp)
    mkdir -p "calc_$name"
    cp "$struct" "calc_$name/POSCAR"
    ln -s ../POTCAR "calc_$name/"
    ln -s ../INCAR "calc_$name/"
done
```
