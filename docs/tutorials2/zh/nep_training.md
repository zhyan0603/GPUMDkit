# NEP 训练指南

<div align="center">
  <p>
    <strong>中文</strong> | <a href="../en/nep_training.md">English</a>
  </p>
</div>

本指南提供使用 GPUMDkit 训练 NEP（神经演化势）模型的完整工作流。

## 概述

NEP 训练包括：

1. **数据准备** - 转换和组织训练数据
2. **结构采样** - 选择多样化结构
3. **模型训练** - 训练 NEP 模型（外部）
4. **验证** - 检查训练质量
5. **主动学习** - 迭代改进模型

## 步骤 1：数据准备

### 将 DFT 数据转换为 extxyz

```bash
# 从 VASP OUTCAR
gpumdkit.sh -out2xyz ./vasp_results/

# 从 VASP XDATCAR
gpumdkit.sh -xdat2exyz XDATCAR trajectory.xyz

# 从 CP2K
gpumdkit.sh  # 选择: 1) Format Conversion -> 103

# 从 LAMMPS
gpumdkit.sh -lmp2exyz dump.lammpstrj Li Y Cl
```

### 添加组标签

组标签是 NEP 识别原子类型所必需的：

```bash
gpumdkit.sh -addgroup POSCAR Li Y Cl
```

这将创建带有组标签的 `model.xyz`。

### 合并训练数据

```bash
# 合并多个 extxyz 文件
cat file1.xyz file2.xyz file3.xyz > train.xyz

# 或使用成分分析选择特定成分
gpumdkit.sh -analyze_comp all_structures.xyz
```

## 步骤 2：数据质量检查

### 检查成分

```bash
gpumdkit.sh -analyze_comp train.xyz
```

### 检查最小距离

```bash
# 快速检查（无 PBC）
gpumdkit.sh -min_dist train.xyz

# 精确检查（带 PBC）
gpumdkit.sh -min_dist_pbc train.xyz
```

### 检查属性范围

```bash
# 检查力范围
gpumdkit.sh -range train.xyz force

# 检查能量范围
gpumdkit.sh -range train.xyz energy

# 带直方图
gpumdkit.sh -range train.xyz force hist
```

### 过滤结构

```bash
# 按最小距离过滤
gpumdkit.sh -filter_dist_pbc train.xyz 1.0

# 按盒子大小过滤
gpumdkit.sh -filter_box train.xyz 20

# 按力值过滤
gpumdkit.sh -filter_value train.xyz force 20
```

## 步骤 3：结构采样

### 采样多样化结构

```bash
# FPS 采样（首选）
gpumdkit.sh  # 选择: 2) Sample Structures -> 203

# 或均匀采样
python Scripts/sample_structures/sample_structures.py train.xyz uniform 100
```

### 扰动结构

从现有结构生成额外的训练数据：

```bash
python Scripts/sample_structures/perturb_structure.py POSCAR 20 0.03 0.2 uniform
```

## 步骤 4：NEP 模型训练

### 准备 nep.in

创建包含训练参数的 `nep.in`：

```
type 3 Li Y Cl
cutoff 6 4
n_max 4
l_max 4
neuron 100
batch 1000
generation 100000
```

### 运行训练

```bash
# GPUMDkit 外部
gpumd
```

## 步骤 5：验证

### 绘制训练结果

```bash
# 绘制训练损失和对比图
gpumdkit.sh -plt train

# 绘制测试预测
gpumdkit.sh -plt prediction

# 组合训练/测试图
gpumdkit.sh -plt train_test

# 力误差分析
gpumdkit.sh -plt force_errors
```

### 检查收敛性

```bash
# 绘制损失曲线
gpumdkit.sh -plt train
```

观察：
- 损失平滑下降
- 无振荡
- 收敛到低值

### 分析描述符

```bash
# 计算描述符
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li

# 使用 PCA 可视化
gpumdkit.sh -plt des pca

# 或 UMAP
gpumdkit.sh -plt des umap
```

## 步骤 6：主动学习

### 使用当前 NEP 运行 MD

```bash
# 设置 MD 模拟
gpumdkit.sh  # 选择: 3) Workflow -> 302
```

### 选择多样化结构

```bash
# 按距离过滤
gpumdkit.sh -min_dist_pbc dump.xyz
gpumdkit.sh -filter_box dump.xyz 13

# 采样多样化结构
gpumdkit.sh  # 选择: 2) Sample Structures -> 203
```

### 计算 DFT 参考

```bash
# 准备 SCF 计算
gpumdkit.sh  # 选择: 3) Workflow -> 301

# 运行 DFT（外部）
```

### 重新训练模型

```bash
# 将新数据添加到训练集
cat new_structures.xyz >> train.xyz

# 重新训练 NEP（外部）
```

## 完整工作流示例

```bash
# 1. 转换 DFT 数据
gpumdkit.sh -out2xyz ./dft_results/

# 2. 添加组标签
gpumdkit.sh -addgroup POSCAR Li Y Cl

# 3. 合并训练数据
cat */model.xyz > train.xyz

# 4. 检查数据质量
gpumdkit.sh -analyze_comp train.xyz
gpumdkit.sh -min_dist_pbc train.xyz
gpumdkit.sh -range train.xyz force

# 5. 过滤结构
gpumdkit.sh -filter_dist_pbc train.xyz 1.0

# 6. 采样多样化结构
gpumdkit.sh  # 选择: 2) Sample Structures -> 203

# 7. 训练 NEP（外部）
# ... 创建 nep.in 并运行 gpumd ...

# 8. 验证训练
gpumdkit.sh -plt train
gpumdkit.sh -plt prediction
gpumdkit.sh -plt force_errors

# 9. 分析描述符
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li
gpumdkit.sh -plt des pca
```

## 更好训练的技巧

1. **多样化数据**：包含各种结构、温度和成分
2. **平衡数据集**：确保所有元素类型都有良好表示
3. **质量过滤**：删除原子太近或力太大的结构
4. **足够数据**：至少 1000-5000 个结构以获得良好精度
5. **验证集**：保留 10-20% 的数据用于测试

## 故障排除

### 问题：高 RMSE

**解决方案**：
- 添加更多多样化训练数据
- 增加 `neuron` 和 `n_max` 参数
- 检查训练数据中的异常值

### 问题：预测不准确

**解决方案**：
- 验证训练数据质量
- 添加来自不同相/温度的结构
- 增加训练 `generation`

### 问题：过拟合

**解决方案**：
- 使用更多训练数据
- 降低模型复杂度
- 使用正则化

## 参见

- [格式转换](format_conversion.md) - 转换文件格式
- [计算器](calculators.md) - 计算属性
- [可视化](visualization.md) - 绘制结果
- [结构采样](sampling.md) - 采样结构
