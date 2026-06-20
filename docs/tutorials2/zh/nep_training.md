<div align="center">
  <h1>NEP 训练指南</h1>
  <p>
    <strong>简体中文</strong> | <a href="../en/nep_training.md">English</a>
  </p>
</div>

训练 NEP（神经演化势）模型的完整工作流。

## 概述

1. **数据准备** - 转换和组织训练数据
2. **质量检查** - 验证和过滤结构
3. **采样** - 选择多样化结构
4. **训练** - 训练 NEP 模型（外部）
5. **验证** - 检查训练质量
6. **主动学习** - 迭代改进模型

## 步骤 1：数据准备

```bash
# 转换 DFT 数据
gpumdkit.sh -out2xyz ./vasp_results/

# 添加组标签
gpumdkit.sh -addgroup POSCAR Li Y Cl

# 合并训练数据
cat */model.xyz > train.xyz
```

## 步骤 2：质量检查

```bash
# 检查成分
gpumdkit.sh -analyze_comp train.xyz

# 检查距离
gpumdkit.sh -min_dist_pbc train.xyz

# 检查力范围
gpumdkit.sh -range train.xyz force

# 过滤结构
gpumdkit.sh -filter_dist_pbc train.xyz 1.0
```

## 步骤 3：采样

```bash
# FPS 采样（推荐）
gpumdkit.sh  # 选择: 2) Sample Structures -> 203

# 或均匀采样
python Scripts/sample_structures/sample_structures.py train.xyz uniform 100
```

## 步骤 4：训练

创建 `nep.in`：

```
type 3 Li Y Cl
cutoff 6 4
n_max 4
l_max 4
neuron 100
batch 1000
generation 100000
```

运行训练：

```bash
gpumd
```

## 步骤 5：验证

```bash
# 绘制训练结果
gpumdkit.sh -plt train
gpumdkit.sh -plt prediction
gpumdkit.sh -plt force_errors

# 分析描述符
gpumdkit.sh -calc des train.xyz descriptors.npy nep.txt Li
gpumdkit.sh -plt des pca
```

## 步骤 6：主动学习

```bash
# 使用当前 NEP 运行 MD
gpumdkit.sh  # 选择: 3) Workflow -> 302

# 过滤和采样
gpumdkit.sh -min_dist_pbc dump.xyz
gpumdkit.sh -filter_box dump.xyz 13

# 准备 DFT 计算
gpumdkit.sh  # 选择: 3) Workflow -> 301

# 用新数据重新训练
cat new_structures.xyz >> train.xyz
```

## 建议

- 包含多样化的结构、温度和成分
- 至少 1000-5000 个结构以获得良好精度
- 保留 10-20% 的数据用于测试
- 训练前检查异常值
