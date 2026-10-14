# 批处理与手动主动学习工作流

本参考用于准备 SCF/MD 批处理目录，以及组织经过人工审查的 NEP 主动学习
迭代。GPUMDkit 负责准备目录和数据，不替用户选择科学设置，也不会自动提交
完整的主动学习循环。

## 目录

- 可用工作流
- SCF 批处理准备
- MD 采样批处理准备
- 手动主动学习协议
- 验证与交接记录

## 可用工作流

| 工作流 | 菜单 | 说明 |
|---|---|---|
| SCF 批处理（VASP） | 301 -> 1 | 准备 VASP 单点计算目录 |
| SCF 批处理（CP2K） | 301 -> 2 | 准备 CP2K 单点计算目录 |
| MD 批处理（GPUMD） | 302 | 准备 GPUMD MD 采样目录 |
| MD 批处理（LAMMPS） | 303 | 准备 LAMMPS MD 采样目录 |

应在专用工作副本中运行这些工具。VASP/XYZ 分支可能会把输入结构移动到
生成的 `struct_fp/` 或 `struct_md/` 目录；继续前应在其他位置保留源结构。

## SCF 批处理准备

### VASP SCF 批处理

1. 把待处理结构放入专用工作目录。工具接受 `.vasp` 文件或一个选定的
   命名为 `*.xyz` 的多帧 extxyz 文件；本流程禁止混放格式或放入多个结构
   XYZ 文件。当前脚本混放时优先处理 `.vasp`，但不得依赖该回退行为。
2. 运行菜单，并按提示输入目录前缀：

```bash
gpumdkit.sh  # 选择：3) Workflow -> 301 -> 1) VASP
```

3. 准备完成后，把用户批准的 VASP 输入文件放入 `fp/`。extxyz 中只有一种
   元素集合时，沿用原来的平铺结构并共用 `fp/POTCAR`。检测到多个元素集合时，
   POSCAR 放入 `struct_fp/<元素顺序>/`，各计算目录链接对应的
   `fp/POTCAR_<元素顺序>`。元素比例不同但元素集合相同的结构仍共用一组。
   分组依据是源 extxyz 中每帧的元素符号。转换脚本沿用轨迹中元素首次出现的
   顺序，并通过读回 POSCAR 核对元素顺序和数量。所有组共用 `fp/INCAR`，
   POTCAR 元素顺序要与 POSCAR 头部一致。每个计算目录都会链接到
   `../fp/KPOINTS`。显式指定 k 点时，请在 `fp/` 中提供该文件。如果目标文件
   不存在，链接会悬空，VASP 会按 KPOINTS 文件缺失处理并使用 `KSPACING`。

```text
struct_fp/Al_O/POSCAR_1.vasp、...           （仅在元素集合不同时分组）
fp/INCAR、POTCAR_Al_O、...                 （准备后由用户提供）
<prefix>_1/POSCAR -> ../struct_fp/Al_O/POSCAR_1.vasp
<prefix>_1/POTCAR -> ../fp/POTCAR_Al_O
<prefix>_1/KPOINTS -> ../fp/KPOINTS
<prefix>_1/INCAR  -> ../fp/INCAR
presub.sh
```

`presub.sh` 只是模板，不代表已经获准运行 VASP。使用前应审查可执行文件、
资源、调度器集成和全部 DFT 输入。

### CP2K SCF 批处理

优先使用正式菜单入口：

```bash
gpumdkit.sh  # 选择：3) Workflow -> 301 -> 2) CP2K
```

交互提示需要：

```text
<extxyz_file> <template.inp> <prefix_name>
```

CP2K 模板必须从 `pos.xyz` 读取坐标。工具写出
`<prefix>_<index>/input.inp` 和 `pos.xyz`，只准备文件，不提交 CP2K。

## MD 采样批处理准备

### GPUMD MD 批处理

1. 把 `.vasp` 结构或一个选定的 extxyz 文件放入专用工作目录。
2. 运行：

```bash
gpumdkit.sh  # 选择：3) Workflow -> 302
```

3. 准备完成后，把批准的 `nep.txt` 以及与每个样本对应的
   `run_<index>.in` 放入生成的 `md/` 目录。

```text
struct_md/model_1.xyz、model_2.xyz、...
md/nep.txt、run_1.in、run_2.in、...       （准备后由用户提供）
sample_1/model.xyz -> ../struct_md/model_1.xyz
sample_1/nep.txt   -> ../md/nep.txt
sample_1/run.in    -> ../md/run_1.in
presub.sh
```

### LAMMPS MD 批处理

运行：

```bash
gpumdkit.sh  # 选择：3) Workflow -> 303
```

准备完成后，把批准的 `lmprun.in` 和 `nep.txt` 放入生成的 `md/` 目录。
工具写出 `struct_md/lammps_<index>.data`、`sample_<index>/` 和
`presub.sh`。

对于两个 MD 工作流，执行前都要检查每个符号链接和输入映射。只有得到用户
明确授权后，才能运行生成的模板并启动 MD。

## 手动主动学习协议

### 目录约定与阶段衔接

统一使用 `iter00`、`iter01` 等，不使用 `round_00`。每个阶段的输入、输出、
日志、构型映射和检查记录都放在该阶段内。下列结构是流程约定，并非 GPUMDkit
自动生成的完整目录树：

```text
active_learning/
+-- iter00/
|   +-- perturbation/
|   +-- nep89md/
|   +-- 01.fps/
|   +-- 02.scf/
|   +-- 03.neptrain/
|   +-- 04.diffusion/
+-- iter01/
|   +-- 00.md/
|   +-- 01.fps/
|   +-- 02.scf/
|   +-- 03.neptrain/
|   +-- 04.diffusion/
+-- iter02/
    +-- ...
```

`iter00` 是初始数据准备阶段，可按用户选择使用微扰、NEP89 MD 或两者结合。
后续迭代通常从 `00.md` 开始，使用上一轮 `03.neptrain/nep.txt` 中通过验收的模型。
`04.xx` 是占位名，必须将 `xx` 替换成实际验证项目，如 `diffusion`；更多性质
阶段继续编号，如 `05.elastic`。只创建实际需要的阶段。工具生成的 `md/`、
`fp/`、`struct_md/`、`struct_fp/` 和计算目录保留在所属阶段内，避免在平级
另建用途重叠的候选、转换或训练目录。未经授权，不迁移已有项目目录。

下表文件名用于说明衔接。保留工具实际输出名，或记录复制/改名关系；不要假定
工具会自动生成这些交接文件名。

| 阶段 | 输入与模型来源 | 输出及下游 |
|---|---|---|
| `iter00/perturbation` | 批准的初始结构和微扰参数；不需要模型 | `POSCAR_*.vasp`、`perturb.txt`；候选几何交给 `01.fps` |
| `iter00/nep89md` | 批准的初始结构、NEP89 模型/版本和 MD 输入 | 原始轨迹与 MD 日志；候选几何交给 `01.fps` |
| `iterNN/00.md`（`NN >= 01`） | 起始结构、上一轮 `iter(NN-1)/03.neptrain/nep.txt` 中验收的模型、批准的 MD 输入 | 原始轨迹与 MD 日志；候选几何交给 `01.fps` |
| `01.fps` | 上游几何、批准的描述符模型与基准构型集 | 原始合并候选、过滤后候选、`selected.xyz`、FPS 输出和源帧映射；入选几何交给 `02.scf` |
| `02.scf` | 入选几何的工作副本、批准的 DFT 设置 | 计算输入/结果、完成与收敛检查、转换后的 DFT 参考数据（如 `new_reference.xyz`）；合格标签交给 `03.neptrain` |
| `03.neptrain` | 合格新 DFT 标签、已有数据集（如有）、合并决定与 `nep.in` | 最终累计 `train.xyz`、用户批准时才有的 `test.xyz`、`nep.txt`、训练日志和数据成员记录；模型交给验证和下一轮 |
| `04.<property>` | 本轮模型、批准的性质计算方案与参考 | 性质结果、检查和验收决定；为下一轮提供反馈 |

### 参数决定与执行授权

在已授权任务内，目录命名、不覆盖原文件的复制和常规整理可以直接执行。
科学参数沿用用户已确认的项目设置；缺失时先提出具体建议及理由，供审核后
再写入可执行输入。应明确：

- 应用范围的组成/相、结构、NEP89 适用性或上一轮模型；
- 微扰数量、原子位移（Angstrom）、晶胞微扰比例和分布（`normal`、`uniform`、`const`）；
- MD 温度、时间步长、时长、采样间隔和系综；NPT 还需明确压力、`iso/aniso/tri`、
  恒温/恒压耦合周期及单位；
- 几何筛选、描述符模型、FPS 基准集与选样数量/停止条件；
- DFT 泛函、赝势、截断能、k 点处理、电子收敛、标签和任务预算；
- 数据纳入、权重、是否划分；NEP 架构（版本、截断、基组、阶数和神经元）、
  损失权重、ZBL 和训练代数；
- 验证方案、验收条件及迭代停止标准。

用户修改参数后，同步更新受影响的输入与阶段记录，标明旧参数结果已被替代，
重新检查依赖阶段；不得继续默用旧方案。材料特定参数和文献设置只能作为示例，
不能当作通用默认值。

准备文件不等于授权运行 MD、DFT、NEP 或提交作业。只要求准备时，完成输入检查
后停止；已授权执行某阶段时，在批准范围内继续完成，不反复询问同一项授权。
只有缺失决定或超出原范围的变化才需再询问。从最后一个完成且检查通过的阶段
继续，不为重建目录而覆盖结果或重跑计算。

### 阶段 0：初始候选或后续 MD

微扰使用菜单 `2 -> 204`，读取 `sampling.md`，明确记录
`<input.vasp> <pert_num> <cell_pert> <atom_pert> <style>`。
`POSCAR_*.vasp` 和 `perturb.txt` 保留在 `perturbation/`。报告只描述变化，
不自动判定结构是否可接受。按照 `format-conversion.md` 支持的方式，在
`01.fps` 内转换几何以供 FPS 使用。

NEP89 或后续 MD 需读取相应 GPUMD 输入、系综和输出参考，检查元素覆盖、
模型来源及安装版本，在 `nep89md/` 或 `00.md/` 准备 `model.xyz`、`run.in`
和模型。批处理可使用菜单 302，生成目录留在该阶段内。只有获授权才运行。
保留原始轨迹和日志；检查稳定性，区分平衡与正式采样，再交接候选几何。

NEP89 和后续 NEP 轨迹中的能量/力是预测值，不是 DFT 标签。不得混入参考数据
或最终训练输入；替换标签必须来自对应构型实际完成且检查通过的 DFT 结果。

### 阶段 1：几何检查与 NepTrain FPS

在 `01.fps` 操作。合并选定的上游候选，保留源文件、运行、帧号和迭代编号，
保留原始候选，单独记录过滤与排除。只使用批准且有物理意义的阈值：

```bash
gpumdkit.sh -min_dist_pbc <candidates.xyz>
gpumdkit.sh -filter_dist_pbc <candidates.xyz> <minimum_distance>
# Optional: only with an approved box criterion, using the actual filtered file
gpumdkit.sh -filter_box <filtered_candidates.xyz> <maximum_box_edge>
gpumdkit.sh  # Select: 2) Sample Structures -> 203) FPS by NepTrain
```

`-min_dist_pbc` 只报告距离，`-filter_dist_pbc` 才实际剔除结构。盒子筛选写出
`filtered_by_box.xyz`。将实际保留的文件交给 FPS。读取 `sampling.md` 的提示：
候选文件、基准构型文件、`nep.txt`，以及距离阈值或结构数量上下限。优先 NepTrain。
明确描述符模型来自哪轮及其版本；候选与基准必须由同一个模型计算描述符。
后续基准通常是上一轮累计训练几何；`iter00` 需明确选择种子/基准文件，因为
当前工具要求该文件。不得为基准构型编造 DFT 标签。

当前实现使用每个结构的原子 NEP 描述符均值。记录阈值/数量范围和实际入选数；
描述符距离依赖模型，不是经过标定的能量/力误差。将 `selected.xyz`、`select.png`、
`pca_sample.txt`、`pca_train.txt`、`pca_selected.txt` 留在 `01.fps`。
另行保留入选帧映射，PCA 文件本身不构成来源记录。更换描述符模型后，两组描述符
都需重算并重新检查距离尺度。PCA 图只是诊断，不能证明覆盖充分。

### 阶段 2：DFT 预处理、标注与转换

1. 只复制入选几何到 `02.scf`。执行 **301** 前，工作目录中的结构输入必须
   是一批 `.vasp`，或单个命名为 `*.xyz` 的多帧 extxyz；不能混放，也不能有
   多个结构 XYZ 文件。本约定比当前 VASP 脚本更严格：脚本混放时优先处理
   `.vasp`，流程不得依赖这一回退行为。
2. 使用 `gpumdkit.sh` 菜单 `3 -> 301`，选择相应 DFT 分支。VASP 先生成计算
   目录、`struct_fp/` 和 `fp/`，再将批准的 `INCAR`、`POTCAR`（或分组 POTCAR）
   及其他输入放入 `fp/`。`INCAR` 已批准使用 `KSPACING` 时，不要求 `KPOINTS`；
   检查没有意外残留的 KPOINTS 覆盖它。CP2K 则在 301 提示时提供要求的已批准
   模板，见前述 CP2K 部分。
3. 检查输入/元素顺序、链接、入选帧到计算目录的映射、资源和任务数。在用户
   执行授权范围内才运行 DFT。
4. 逐个检查入选构型对应任务、完成状态、电子收敛、有限标签和几何对应关系。
   明确报告缺失、失败或未收敛任务及其编号，不能悄悄将部分结果当作最终版本；
   在完成重算或用户批准排除后，才定稿。
5. 在 `02.scf` 内通过相应入口转换合格结果：

```bash
gpumdkit.sh -out2xyz <accepted_scf_results_directory>
gpumdkit.sh -cp2k2xyz  # CP2K alternative: interactive converter
# Inspect the actual extxyz output from the chosen converter
gpumdkit.sh -range <new_reference.xyz> energy
gpumdkit.sh -range <new_reference.xyz> force
gpumdkit.sh -min_dist_pbc <new_reference.xyz>
gpumdkit.sh -analyze_comp <new_reference.xyz>
```

读取 `nep-data.md` 和 `format-conversion.md`。确认能量为 eV（每个构型总能量）、
力为 eV/Angstrom；存在 virial 时检查其约定和单位。核对晶胞、PBC、元素/数量、
有限标签及一一对应的构型映射。文件存在或转换成功不代表电子收敛。
原始候选、入选几何、合格 DFT 参考标签是不同产物，只有最后一项提供新训练标签。

### 阶段 3：累计数据与 NEP 训练

在 `03.neptrain` 内将合格新标签与上一版检查过的数据集合并，不修改上一版。
检查重复结构与参考计算设置的一致性，记录最终每个构型的源帧、DFT 任务、
纳入/权重/划分。`train.xyz` 是最终累计训练输入，不是未检查的轨迹或转换结果拼接。

不要默认创建测试集。由用户决定是否划分及具体方式；若划分，避免相关轨迹帧/
母结构跨集合泄漏，并记录成员关系。若不需要测试集，省略 `test.xyz`，不得复制
训练数据充当测试数据。核实安装 NEP 版本对可选测试集的支持，并说明尚未测量
独立留出数据误差。

读取 `nep.md`、`nep-parameters.md` 和 `nep-outputs.md`，只写入已确认的
`nep.in` 参数。获授权训练后，检查损失、适用的能量/力/virial 误差、异常点及
各元素/环境误差。模型来源和训练输出留在本阶段；将模型交给性质验证，不仅凭
总体 RMSE 验收。

### 阶段 4 及以后：模型验证与下一轮

按实际性质命名目录，如 `04.diffusion`。依据批准的方案和参考检查目标应用中的
动力学、几何、稳定性及性质；可考虑状态方程、弹性、相/缺陷相对能、迁移势垒
或扩散。只选择与用途有关的检查；文献不提供通用阈值，也不授权昂贵验证计算。

分别报告训练拟合、独立验证（如已进行）、性质一致性和未覆盖范围。训练拟合
或 MD 稳定本身不能证明可迁移性。记录验收的模型版本，将它和已检查的累计
数据交给下一轮。未通过时保留证据并建议有针对性的候选。停止迭代需满足用户
约定的覆盖、误差/性质及预算条件，不能只看 FPS 距离小。

### 可选选样路线：委员会不确定度

只在用户选择时使用。读取 `gpumd-outputs.md`：GPUMD `active` 使用兼容的
NEP 委员会写出 `active.out` 和 `active.xyz`。确认模型顺序、检查间隔、输出
开关、阈值、`top_n` 和 `min_deviation`；获授权运行后先检查日志和几何，再用
菜单 205。配对输出留在采样阶段。委员会分歧是选样信号，不是保证误差上界；
相关模型可能有共同盲区。

### 方法建议与文献

以下文献支持方法选择，项目是否采用仍需用户批准。

- **覆盖：** 描述符多样性补充 MD 采样；除轨迹数量外，检查组成、相和局部环境
  覆盖。参见
  [Training data selection for accuracy and transferability of interatomic potentials](https://doi.org/10.1038/s41524-022-00872-x)。
- **罕见环境：** 普通 MD 可能不足以采到缺陷和过渡态；必要时建议有针对性的
  起始结构或增强采样，并审核具体方案。参见
  [Hyperactive learning for data-driven interatomic potentials](https://doi.org/10.1038/s41524-023-01104-6)。
- **验证：** 平均能量/力误差低，仍可能存在动力学、缺陷和罕见事件误差；检查
  相关局部误差及目标性质。参见
  [Discrepancies and error evaluation metrics for machine learning interatomic potentials](https://doi.org/10.1038/s41524-023-01123-3)。

对于均值描述符 FPS，另行检查少数元素、缺陷及异常局部环境，因为结构平均
可能掩盖它们。这是当前表示方式的实际局限，不意味着每个此类候选都必须标注。
选样、参考标注、再训练和验证始终是不同阶段。

## 验证与交接记录

每个阶段内保留简洁记录，不另建平行 QA 目录。记录输入路径/模型版本、确认
参数及修改、授权范围、命令/工作目录/版本/状态、输出路径和帧映射、检查与
验收、失败/排除项及剩余限制。继续迭代时根据实际文件核对记录，指出下一个
未完成阶段。必要输出缺失或错误未解决时，停止依赖它的工作，明确报告，不得
宣称本轮已完成。

## 详细文档

- `${GPUMDkit_path}/docs/tutorials/zh/工作流脚本.md`
- `${GPUMDkit_path}/docs/tutorials/zh/主动学习工作流.md`
