# 结构采样与数据集划分

## 本页内容

- 可用方法
- 命令参考
- 常见工作流程
- 依赖与注意事项

## 可用方法

| 方法 | 菜单/入口 | 用途 |
|---|---|---|
| 均匀/随机采样 | 201 | 从轨迹中等间隔或随机抽取结构 |
| PyNEP 提示 | 202 | 仅显示迁移提示；如需兼容入口，使用 `gpumdkit.sh -pynep` |
| PyNEP FPS | `-pynep` | 最远点采样的旧兼容入口，已弃用 |
| NepTrain FPS | 203 | 基于描述符的最远点采样，推荐使用 |
| 结构扰动 | 204 | 从已有结构生成扰动构型 |
| 力偏差筛选 | 205 | 从主动学习输出中挑选高力偏差结构 |
| 训练集/测试集划分 | 206 | 使用均匀、随机或 NepTrain FPS 方法划分 extxyz 数据集 |

## 命令参考

### 均匀与随机采样

```bash
# 直接运行 Python 脚本
python3 ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py <input.xyz> <method> <num_samples> [skip_initial]

# 示例
python3 ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py dump.xyz uniform 50
python3 ${GPUMDkit_path}/Scripts/sample_structures/sample_structures.py dump.xyz random 100 500
```

参数：

- `input.xyz`：输入的 extxyz 轨迹；
- `method`：`uniform` 或 `random`；
- `num_samples`：需要抽取的帧数；
- `skip_initial`：可选，跳过开头的 N 帧，默认值为 0。

输出文件为 `sampled_structures.xyz`。

- **Uniform**：通过 `numpy.linspace` 在可用轨迹范围内生成等间隔索引。
- **Random**：通过 `numpy.random.choice` 无放回随机抽样。

### 最远点采样（FPS）

#### NepTrain 版本（推荐）

```bash
# 交互模式：2) Sample Structures -> 203
gpumdkit.sh

# 半交互式直接调用
python3 ${GPUMDkit_path}/Scripts/sample_structures/parallel_neptrain_select_structs.py <sample.xyz> <train.xyz> <nep.txt> [threads]

# 示例
python3 ${GPUMDkit_path}/Scripts/sample_structures/parallel_neptrain_select_structs.py dump.xyz train.xyz nep.txt 4
```

`threads` 为可选参数，默认值为 `1`。设为正整数即可并行计算描述符。
每个工作进程会独立加载 NEP 模型并使用 1 个原生 OpenMP 线程，结果保持输入帧顺序。

运行后可选择：

1. 按最小描述符距离阈值选择，直到最大距离小于阈值；
2. 按结构数量选择，指定最少和最多结构数。

输出：

- `selected.xyz`；
- `select.png`（PCA 可视化）；
- `pca_sample.txt`、`pca_train.txt`、`pca_selected.txt`。

算法使用 NepTrain 计算逐原子 NEP 描述符，并以每个结构的平均描述符进行 FPS。

使用此功能时建议引用：Chen et al., *Comput. Phys. Commun.* 317, 109859 (2025)。

#### PyNEP 版本（已弃用）

```bash
gpumdkit.sh -pynep
```

PyNEP 已不再积极维护。新任务应优先使用功能 203；`-pynep` 仅用于兼容既有流程。

### 结构扰动

```bash
# 交互模式：2) Sample Structures -> 204
gpumdkit.sh

# 直接调用
python3 ${GPUMDkit_path}/Scripts/sample_structures/perturb_structure.py <input.vasp> <pert_num> <cell_pert> <atom_pert> <style>

# 示例
python3 ${GPUMDkit_path}/Scripts/sample_structures/perturb_structure.py POSCAR 20 0.03 0.2 uniform
```

参数：

- `input.vasp`：POSCAR/CONTCAR 文件；
- `pert_num`：生成的扰动结构数；
- `cell_pert`：晶胞扰动比例，例如 `0.03` 表示 3%；
- `atom_pert`：原子位置扰动距离，单位为 Å；
- `style`：`normal`、`uniform` 或 `const`。

脚本使用 `dpdata.System.perturb()` 同时扰动晶胞矢量和原子位置，输出
`POSCAR_01.vasp` 至 `POSCAR_<N>.vasp`。

### 力偏差筛选

```bash
# 交互模式：2) Sample Structures -> 205
gpumdkit.sh

# 直接调用
python3 ${GPUMDkit_path}/Scripts/sample_structures/select_max_modev.py <top_n> <min_deviation>

# 示例
python3 ${GPUMDkit_path}/Scripts/sample_structures/select_max_modev.py 200 0.15
```

当前目录必须包含来自同一次 GPUMD `active` 运行的：

- `active.out`；
- `active.xyz`。

脚本先筛选力偏差高于 `min_deviation` 的结构，再按偏差从高到低排序，最多保留 `top_n` 个结构，并从 `active.xyz` 提取对应帧写入 `selected.xyz`。

### 训练集和测试集划分

```bash
# 交互模式：2) Sample Structures -> 206
gpumdkit.sh

# 直接进入脚本；其余参数仍由交互提示输入
python3 ${GPUMDkit_path}/Scripts/sample_structures/split_train_test.py <input.xyz>

# 示例
python3 ${GPUMDkit_path}/Scripts/sample_structures/split_train_test.py data.xyz
```

脚本先输出数据集中的元素、总帧数以及每帧原子数范围，然后询问测试集大小：

- 严格位于 `0` 和 `1` 之间的数按比例解释，例如 `0.1` 表示 10%；
- 正整数按精确帧数解释，例如 `100` 表示 100 帧；
- 比例换算后的帧数按四舍五入（半数向上）处理，并至少选择 1 帧；
- 测试集必须小于完整数据集，确保训练集不为空。

可选方法：

- **Uniform**：在完整数据集索引范围内等间隔选择测试结构；
- **Random**：无放回随机选择；可提供整数随机种子以复现结果；
- **FPS**：使用 NepTrain 计算逐结构平均的 NEP 原子描述符，从第 1 帧开始，逐次加入与已选结构距离最远的候选结构；需要兼容的 `nep.txt`。

若输入为 `data.xyz`，输出为 `data_train.xyz` 和 `data_test.xyz`。两个文件中的结构都保持其在原始数据集中的先后顺序。选中的结构写入测试集，其余结构全部写入训练集。

### 帧区间提取

```bash
gpumdkit.sh -frame_range <input.xyz> <start_fraction> <end_fraction>
python3 ${GPUMDkit_path}/Scripts/sample_structures/frame_range.py <input.xyz> <start_fraction> <end_fraction>

# 示例：保留轨迹前 80%
gpumdkit.sh -frame_range dump.xyz 0 0.8
```

输出文件名为 `dump_0.00_0.80.xyz`。起止参数都是 0 到 1 之间的轨迹比例。

## 常见工作流程

### 准备训练数据

```bash
# 1. 对候选轨迹检查并筛除不合理结构
gpumdkit.sh -min_dist_pbc dump.xyz
gpumdkit.sh -filter_dist_pbc dump.xyz 1.0

# 2. 使用 FPS 选取多样化结构
gpumdkit.sh  # 选择 2 -> 203

# 3. 经确认后加入训练集
cat selected.xyz >> train.xyz
```

### 主动学习数据筛选

```bash
# 1. 从 active.out/active.xyz 中挑选高偏差结构
python3 ${GPUMDkit_path}/Scripts/sample_structures/select_max_modev.py 100 0.1

# 2. 如需进一步提高多样性，再使用 NepTrain FPS
python3 ${GPUMDkit_path}/Scripts/sample_structures/parallel_neptrain_select_structs.py active.xyz train.xyz nep.txt 4
```

### 生成初始扰动结构

```bash
python3 ${GPUMDkit_path}/Scripts/sample_structures/perturb_structure.py POSCAR 50 0.05 0.3 uniform
```

在进行 DFT 之前，应检查扰动结构的最小原子间距和晶胞是否合理。不得擅自采用示例中的扰动幅度作为正式参数。

## 依赖

| 方法 | Python 依赖 |
|---|---|
| 均匀/随机采样 | `numpy`、`ase` |
| NepTrain FPS | `numpy`、`ase`、`matplotlib`、`scikit-learn`、`scipy`、`NepTrain` |
| PyNEP FPS | `numpy`、`ase`、`matplotlib`、`scikit-learn`、`pynep` |
| 结构扰动 | `dpdata` |
| 力偏差筛选 | `numpy`、`ase` |
| 训练集/测试集划分（Uniform/Random） | `numpy`、`ase` |
| 训练集/测试集划分（FPS） | `numpy`、`ase`、`scipy`、`NepTrain` |

```bash
pip install numpy ase matplotlib scikit-learn scipy dpdata
pip install neptrain
```

## 注意事项

1. FPS 通常更适合追求结构多样性的任务，但是否采用仍取决于数据用途。
2. 新流程优先使用 NepTrain；PyNEP 入口仅为兼容旧流程保留。
3. 采样前先检查并筛除原子距离过短、晶胞异常等非物理结构。
4. `-frame_range` 的参数是轨迹比例，不是帧编号。
5. 本模块的采样输出均采用 extxyz 格式。
6. 功能 206 选择的是测试集结构；所有未被选中的结构写入训练集。

## 用户教程

需要面向用户的完整说明时，请参阅
`${GPUMDkit_path}/docs/tutorials/en/structure_sampling.md` 或对应的中文教程。
