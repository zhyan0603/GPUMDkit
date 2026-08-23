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
   extxyz 文件；两种格式同时存在时优先处理 `.vasp`。
2. 运行菜单，并按提示输入目录前缀：

```bash
gpumdkit.sh  # 选择：3) Workflow -> 301 -> 1) VASP
```

3. 准备完成后，把用户批准的 `INCAR`、`POTCAR` 和 `KPOINTS` 放入生成的
   `fp/` 目录。

```text
struct_fp/POSCAR_1.vasp、POSCAR_2.vasp、...
fp/INCAR、POTCAR、KPOINTS                 （准备后由用户提供）
<prefix>_1/POSCAR -> ../struct_fp/POSCAR_1.vasp
<prefix>_1/INCAR  -> ../fp/INCAR
<prefix>_2/...
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

每轮迭代应由经过审查的 GPUMDkit 阶段组成，使体系相关的科学选择和调度器
设置保持明确。推荐路线是普通 MD 采样，然后进行几何筛选和 NepTrain FPS
选取。

### 执行前确认本轮方案

要求用户提供或批准：

- 模型修订、目标组成/相和预期应用范围；
- 候选结构来源和 GPUMD 条件；
- MD 轨迹、距离/盒子筛选、NepTrain FPS 规则和结构数量条件；
- DFT 程序、设置、赝势、收敛标准和任务预算；
- 数据集合并/划分规则和数据泄漏控制；
- NEP 训练配置和模型验收标准。

不得用其他体系的示例或默认值替代任何尚未确定的项目。

### 推荐路径：MD 轨迹与 NepTrain FPS

针对用户批准的起始结构和条件运行已获授权的普通 MD。先检查稳定性，再应用
用户批准的筛选条件：

```bash
gpumdkit.sh -min_dist_pbc dump.xyz
gpumdkit.sh -filter_dist_pbc dump.xyz <minimum_distance>
gpumdkit.sh -filter_box filtered_dump.xyz <maximum_box_edge>
gpumdkit.sh  # 选择：2) Sample Structures -> 203) FPS by NepTrain
```

`-min_dist_pbc` 只报告距离，不会筛选结构；需要实际剔除结构时应使用
`-filter_dist_pbc`。盒子筛选是可选步骤，只有盒子边长条件对该数据集有
意义且获得批准时才使用。记录每个被剔除的结构和筛选阈值。NepTrain FPS
随后相对于当前训练集选取能增加描述符空间多样性的结构。记录 FPS 停止
规则和入选帧映射。

### 可选替代：委员会不确定度

只有项目明确选择时才使用委员会不确定度。读取 `gpumd-outputs.md`：GPUMD
`active` 需要兼容的 NEP 委员会模型，并写出 `active.out` 和 `active.xyz`。
确认模型顺序、检查间隔、输出开关、阈值、`top_n` 和 `min_deviation`；得到
授权并完成运行后，先验证日志和结构，再使用菜单 205。确保 `active.out`
和 `active.xyz` 来自同一次运行。

### 参考计算与数据集更新

1. 使用菜单 301 从批准的 `selected.xyz` 准备 DFT 目录。
2. 审查生成的结构、DFT 输入、资源请求和任务数量。
3. 通过用户自己的调度流程，并在另行获得明确授权后运行 DFT。
4. 确认每个计算完整结束；记录原因后排除失败或不完整输出。
5. 根据 DFT 程序选择对应入口转换参考计算输出，再检查转换器生成的数据集：

```bash
# VASP 结果
gpumdkit.sh -out2xyz <scf_results_directory>

# CP2K 结果（交互式转换器）
gpumdkit.sh -cp2k2xyz

# 对所选转换器生成的 extxyz 文件继续检查
gpumdkit.sh -range <new_reference.xyz> energy
gpumdkit.sh -range <new_reference.xyz> force
gpumdkit.sh -min_dist_pbc <new_reference.xyz>
gpumdkit.sh -analyze_comp <new_reference.xyz>
```

6. 读取 `nep-data.md`；合并前检查标签、单位、晶胞、元素/类型顺序、重复
   结构、训练/测试泄漏和构型覆盖。
7. 保留上一版数据集，并准确记录加入、剔除、加权或重新划分的结构。
8. 读取 `nep.md`，修改配置时读取 `nep-parameters.md`，并读取
   `nep-outputs.md`。只有获得授权后才训练；随后评估训练/测试行为、异常点、
   目标范围稳定性和用户批准的验收标准。不得只依据总体 RMSE 批准模型。

## 验证与交接记录

每轮迭代至少报告：

| 字段 | 记录内容 |
|---|---|
| 输入 | 模型哈希/修订、结构、GPUMD 和 DFT 输入 |
| 决策 | 候选路径和每个已批准的科学阈值 |
| 授权 | 哪些 MD、DFT 和 NEP 运行得到明确批准 |
| 命令 | 完整命令、工作目录、版本和退出状态 |
| 输出 | 候选、入选、参考数据、数据集和模型文件 |
| 验证 | 解析器/日志检查、结构/数据审查、误差和稳定性测试 |
| 限制 | 失败任务、排除项、未覆盖范围和未解决风险 |

缺少必要输出、日志存在未解决错误、结构不合理或验证关卡尚未批准时，不得
进入下一阶段。

## 详细文档

- `${GPUMDkit_path}/docs/tutorials/zh/工作流脚本.md`
- `${GPUMDkit_path}/docs/tutorials/zh/主动学习工作流.md`
