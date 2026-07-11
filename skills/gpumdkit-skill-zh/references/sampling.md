# 结构采样

## 目录

- 可用方法
- 命令参考
- 常用工作流
- 依赖与注意事项

## 可用方法

| 方法 | 菜单 | 描述 |
|--------|------|-------------|
| 均匀/随机 | 201 | 选择等间距或随机帧 |
| PyNEP 提示 | 202 | 仅打印提示；通过 `gpumdkit.sh -pynep` 运行 PyNEP |
| FPS（PyNEP） | `-pynep` | 最远点采样（已弃用的兼容入口） |
| FPS（NepTrain） | 203 | 最远点采样（推荐） |
| 微扰 | 204 | 生成微扰结构 |
| 力偏差 | 205 | 选择高力偏差结构 |

## 命令参考

### 均匀/随机采样

```bash
# 直接 Python 执行
python3 ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py <input.xyz> <method> <num_samples> [skip_initial]

# 示例
python3 ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py dump.xyz uniform 50
python3 ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py dump.xyz random 100 500

# 参数：
# - input.xyz：输入轨迹文件
# - method：'uniform'（均匀）或 'random'（随机）
# - num_samples：要选择的帧数
# - skip_initial：（可选）跳过前 N 帧

# 输出：sampled_structures.xyz
```

**工作原理：**
- **均匀**：使用 `numpy.linspace` 在轨迹中选择等间距索引
- **随机**：使用 `numpy.random.choice` 无放回地选择随机索引

### 最远点采样（FPS）

#### 使用 NepTrain（推荐）

```bash
# 交互模式
gpumdkit.sh  # 选择：2) Sample Structures -> 203

# 半交互式直接 Python 执行
python3 ${GPUMDkit_path}/Scripts/sample_structures/neptrain_select_structs.py <sample.xyz> <train.xyz> <nep.txt>

# 示例
python3 ${GPUMDkit_path}/Scripts/sample_structures/neptrain_select_structs.py dump.xyz train.xyz nep.txt

# 选择方法（交互式）：
# 1. 最小距离：选择直到最大距离 < 阈值
# 2. 结构数量：选择 N 个结构

# 输出：
# - selected.xyz
# - select.png（PCA 可视化）
# - pca_sample.txt、pca_train.txt、pca_selected.txt
```

**算法**：使用 NEP 描述符（每结构的平均描述符）的最远点采样

**引用**：Chen et al., Comput. Phys. Commun., 2025, 317, 109859

#### 使用 PyNEP（已弃用）

```bash
# GPUMDkit 入口
gpumdkit.sh -pynep

# 注意：PyNEP 包已不再积极维护
# 请改用 NepTrain（方法 203）
```

### 结构微扰

```bash
# 交互模式
gpumdkit.sh  # 选择：2) Sample Structures -> 204

# 直接 Python 执行
python3 ${GPUMDkit_path}/Scripts/sample_structures/perturb_structure.py <input.vasp> <pert_num> <cell_pert> <atom_pert> <style>

# 示例
python3 ${GPUMDkit_path}/Scripts/sample_structures/perturb_structure.py POSCAR 20 0.03 0.2 uniform

# 参数：
# - input.vasp：POSCAR/CONTCAR 文件
# - pert_num：微扰结构数量（例如 20）
# - cell_pert：晶胞微扰比例（例如 0.03 = 3%）
# - atom_pert：原子微扰距离（埃，例如 0.2）
# - style：'normal'、'uniform' 或 'const'

# 输出：POSCAR_01.vasp、POSCAR_02.vasp、...、POSCAR_<N>.vasp
```

**工作原理**：使用 `dpdata.System.perturb()` 生成输入结构的微扰变体，同时改变晶胞向量和原子位置。

### 力偏差选择

```bash
# 交互模式
gpumdkit.sh  # 选择：2) Sample Structures -> 205

# 直接 Python 执行
python3 ${GPUMDkit_path}/Scripts/sample_structures/select_max_modev.py <top_n> <min_deviation>

# 示例
python3 ${GPUMDkit_path}/Scripts/sample_structures/select_max_modev.py 200 0.15

# 当前目录中需要的文件：
# - active.out（来自 GPUMD active 命令）
# - active.xyz（对应结构）

# 输出：selected.xyz
```

**工作原理**：
1. 读取 `active.out` 并过滤力偏差高于 `min_deviation` 的结构
2. 按偏差降序排序并提取前 N 个结构
3. 从 `active.xyz` 读取对应帧并写入输出

### 帧范围提取

```bash
# GPUMDkit CLI
gpumdkit.sh -frame_range <input.xyz> <start_fraction> <end_fraction>

# 直接 Python 执行也可以
python3 ${GPUMDkit_path}/Scripts/sample_structures/frame_range.py <input.xyz> <start_fraction> <end_fraction>

# 示例：提取前 80% 的帧
gpumdkit.sh -frame_range dump.xyz 0 0.8

# 输出：dump_0.00_0.80.xyz
```

## 常用工作流

### 训练数据准备

```bash
# 1. 运行 MD 模拟
# ... 生成轨迹 ...

# 2. 按距离过滤
gpumdkit.sh -min_dist_pbc dump.xyz
gpumdkit.sh -filter_dist_pbc dump.xyz 1.0

# 3. 采样多样化结构
gpumdkit.sh  # 选择：2) Sample Structures -> 203

# 4. 添加到训练集
cat selected.xyz >> train.xyz
```

### 主动学习数据选择

```bash
# 1. 使用当前 NEP 运行 MD
# ... 生成 active.out 和 active.xyz ...

# 2. 选择高偏差结构
python3 ${GPUMDkit_path}/Scripts/sample_structures/select_max_modev.py 100 0.1

# 3. 或使用 FPS 获取多样性
python3 ${GPUMDkit_path}/Scripts/sample_structures/neptrain_select_structs.py active.xyz train.xyz nep.txt
```

### 初始数据的结构微扰

```bash
# 1. 从平衡结构开始
# 2. 生成微扰
python3 ${GPUMDkit_path}/Scripts/sample_structures/perturb_structure.py POSCAR 50 0.05 0.3 uniform

# 3. 对微扰结构运行 DFT
# ... DFT 计算 ...

# 4. 添加到训练集
```

### PyNEP 兼容路径

```bash
# gpumdkit.sh -pynep 调用并行 PyNEP 兼容脚本。
gpumdkit.sh -pynep

# 启动后的示例输入：
dump.xyz train.xyz nep.txt 8
```

## 依赖

| 方法 | 所需包 |
|--------|------------------|
| 均匀/随机 | `numpy`、`ase` |
| FPS（NepTrain） | `numpy`、`ase`、`matplotlib`、`scikit-learn`、`scipy`、`NepTrain` |
| FPS（PyNEP） | `numpy`、`ase`、`matplotlib`、`scikit-learn`、`pynep` |
| 微扰 | `dpdata` |
| 力偏差 | `numpy`、`ase` |

安装依赖：
```bash
pip install numpy ase matplotlib scikit-learn scipy dpdata
pip install neptrain  # 用于 NepTrain FPS
```

## 注意事项

1. **FPS 是首选**：最远点采样比均匀/随机提供更好的多样性
2. **NepTrain 优于 PyNEP**：NepTrain 积极维护；PyNEP 已弃用
3. **采样前先过滤**：先移除非物理结构（太近、太大）
4. **帧范围使用分数**：`-frame_range` 使用 0.0 到 1.0 之间的起止分数
5. **输出为 extxyz**：所有采样方法输出 extxyz 格式

## 详细文档

用户指南请参见 `${GPUMDkit_path}/docs/tutorials/en/structure_sampling.md` 或中文对应版本。
