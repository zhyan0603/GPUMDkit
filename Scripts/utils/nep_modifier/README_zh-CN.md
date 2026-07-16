# NEP 模型修改器

> [English](README.md) | **简体中文**

NEP 模型修改器是对
[calorine](https://calorine.materialsmodeling.org/dev/get_started/modifying_models.html)
模型修改功能的 GPUMDkit 交互式封装。它将 `augment()`、`prune()`、
`add_species()`、`remove_species()` 和 `keep_species()` 转换成带输入检查、架构查看、修改审阅和文件包导出的终端流程。

该功能仍处于测试阶段。在开始后续训练前，应将生成的模型、restart 和 `nep.in` 作为一个整体进行检查。

## 适用场景

当模型的适用范围发生变化时，并非总要从头训练。该工具适合用于：

- 原模型容量不足时增加神经元或描述符；
- 引入更高体阶描述符或电荷输出头；
- 构建更小的架构，作为面向推理速度的再训练起点；
- 从多元素模型中提取单元素或少量元素子模型；
- 在保留已有元素已学习参数的基础上加入新元素。

模型修改只是新训练流程的起点。增减参数会改变优化问题，新增元素及其描述符对的参数也尚未训练。必须使用能代表目标应用的数据重新训练并验证生成的模型。

## 使用要求

- `calorine >= 3.4`
- `nep.txt` 格式的 NEP4 模型
- 扩展模型、缩减模型或添加元素时，需要匹配的 `nep.restart`
- 可选的源 `nep.in`，用于保留训练相关设置

安装 calorine：

```bash
pip install 'calorine>=3.4'
```

## 使用方法

使用当前目录中的默认文件启动：

```bash
gpumdkit.sh -nep_modifier
```

显式指定模型文件包：

```bash
gpumdkit.sh -nep_modifier path/nep.txt path/nep.restart path/nep.in
```

使用 `-` 跳过 restart 或训练输入：

```bash
gpumdkit.sh -nep_modifier path/nep.txt - path/nep.in
gpumdkit.sh -nep_modifier path/nep.txt path/nep.restart -
```

无需导入 calorine 即可查看帮助：

```bash
gpumdkit.sh -nep_modifier -h
```

未提供路径时，程序会分别询问 `nep.txt`、`nep.restart` 和 `nep.in`。默认路径从所选模型的同一目录解析，不会读取无关工作目录中的文件。

## 主菜单

```text
+-------------------------------------------------------------------------------+
|                              NEP MODEL MODIFIER                               |
+-------------------------------------------------------------------------------+
| 1) Expand model capacity             4) Remove chemical species               |
| 2) Reduce model capacity             5) Keep selected species                 |
| 3) Add chemical species              6) Inspect current model                 |
| 7) Review pending changes            8) Export model package                  |
| 0) Exit                                                                       |
+-------------------------------------------------------------------------------+
```

如果存在尚未导出的修改，退出前程序会进行提醒。

## 功能说明

| 操作 | Calorine 方法 | 需要 restart | 后续训练 |
|---|---|---:|---:|
| 扩展模型容量 | `augment()` | 是 | 必须 |
| 缩减模型容量 | `prune()` | 是 | 必须 |
| 添加元素 | `add_species()` | 是 | 必须 |
| 移除元素 | `remove_species()` | 否 | 可选 |
| 仅保留所选元素 | `keep_species()` | 否 | 可选 |

### 扩展模型容量

一次操作可以选择多项修改，并通过一个原子性的 `augment()` 调用完成：

- 增加 `n_neuron`；
- 增加 `l_max_4b` 或 `l_max_5b`；
- 启用 `q_112`、`q_123`、`q_233` 或 `q_134` 描述符项；
- 添加电荷输出头。

如果当前 calorine 提供 `charge_mode`，用户必须明确选择模式 1 或 2。程序会显示 calorine 记录的 SNES 默认值，并询问是接受默认值，还是手动输入 `sigma_new`、`sigma_factor` 和 `sigma_floor`。

已有参数的训练均值 `mu` 会被保留。新参数按照 calorine 的规则保持静默或完成初始化，相应的 SNES 搜索宽度由所选 sigma 设置生成。因此，后续训练可以调整扩展后的架构，同时保留模型已经学到的部分。

### 缩减模型容量

一次 `prune()` 操作可以组合：

- 减少 `n_neuron`；
- 减少或移除 `l_max_4b`、`l_max_5b` 项；
- 禁用已经启用的高体阶描述符项；
- 从带电荷模型中移除电荷输出头。

减少神经元时，calorine 会根据模型权重对神经元排序。缩减后的模型在用于生产前必须重新训练和验证。
完全关闭某一描述符族会删除对应的描述符列；但只是把非零 `l_max_4b` 调低时，可能仅改变文件头而不会减少已存储的描述符维数。

### 添加元素

`add_species()` 仅适用于带 restart 数据的 NEP4 模型。程序会：

- 拒绝重复元素和模型中已存在的元素；
- 要求用户明确提供随机 seed；
- 对 typewise cutoff 模型逐元素询问径向和角向 cutoff；
- 记录 seed、cutoff 和已接受的 SNES 设置。

新增元素参数尚未经过训练，因此导出的模型在完成适当参考数据训练和验证前不能用于生产计算。

### 移除或保留元素

这两项操作可以在没有 restart 数据时执行。如果已加载 restart，程序会在操作前显示自适应 sigma 的文档默认值。模型中必须至少保留一个元素。

待删除元素较少时适合使用 `remove_species()`；从大型基础模型中提取少量元素时，使用反向选择的 `keep_species()` 通常更方便。两者都会移除被舍弃元素对应的 ANN 子网络和描述符权重对。

## 交互示例

以下示例只展示关键输入，实际界面中的当前值取决于载入的模型。

### 增加神经元并启用描述符项

```text
Input the function number:
------------>>
1
Input one or more choices, separated by spaces:
------------>>
1 4
Input target neuron count (current: 50):
------------>>
60
Use these SNES initialization defaults? (Y/n)
------------>>
y
Apply these changes? (y/N)
------------>>
y
```

这会执行一次组合的 `augment()`：神经元数变为 60，同时启用 `q_112`。导出后，应配套使用生成的 `.in`、`.txt` 和 `.restart` 继续训练。

### 提取 Li-O 子模型

```text
Input the function number:
------------>>
5
Input species to keep, separated by spaces:
------------>>
Li O
Keep only these species? (y/N)
------------>>
y
```

对于大型多元素模型，这比逐一列出所有待删除元素更简洁。元素顺序以 calorine 返回的保留顺序为准；准备新训练数据或 GPUMD `model.xyz` 前应再次检查。

### 使用固定 seed 添加碳元素

```text
Input the function number:
------------>>
3
Input species to add, separated by spaces:
------------>>
C
Input random seed for reproducible initialization:
------------>>
20260712
Use these SNES initialization defaults? (Y/n)
------------>>
y
Add these species? (y/N)
------------>>
y
```

如果原模型使用 typewise cutoff，程序会在 seed 之前询问碳元素的径向和角向 cutoff。这些是具有科学含义的模型选择，应依据目标训练设置确定，不能直接照抄本示例。所选 seed 和 cutoff 会写入修改摘要。

### 审阅并导出

选择 `7` 查看全部已确认操作，然后选择 `8` 并输入输出目录。成功后会打印所有生成文件路径。如果仍有未导出的修改，选择 `0` 时程序会询问是否直接退出。

## 模型检查和修改历史

模型检查页面会显示：

- NEP 版本和模型类型；
- 元素顺序；
- 径向、角向 cutoff；
- `n_max`、`basis_size` 和 `l_max`；
- 高体阶描述符开关；
- 描述符、ANN、神经元和总参数数量；
- restart 状态、ZBL，以及可用时的 charge mode。

`Review pending changes` 会显示已确认操作、参数、发生变化的模型字段，以及每项操作是否已经导出。

## 安全导出

工具不会覆盖源文件包。它会为全部输出选择同一个无冲突后缀，先在临时目录中写出并检查所有文件，然后再移动到用户指定的输出目录。

对于名为 `nep.txt` 的输入，第一组输出为：

| 文件 | 用途 |
|---|---|
| `nep_modified.txt` | 修改后的模型 |
| `nep_modified.restart` | 已加载 restart 时生成的匹配重启文件 |
| `nep_modified.in` | 更新后的架构和保留的训练设置 |
| `nep_modified.changes.txt` | 源路径、calorine 版本、参数、seed 和操作历史 |

只要其中任意文件已经存在，整组输出都会使用下一个公共后缀，例如 `nep_modified_1.*`。

提供源 `nep.in` 时，工具会保留其中的训练设置，并使用 `model.training_parameters` 替换架构字段。未提供源 `nep.in` 时，生成文件只包含架构字段。

必须使用后续训练所采用的 NEP 可执行文件版本检查生成的 `nep.in`。不要把新模型或 restart 与旧的、不兼容的 `nep.in` 混合使用。

## 错误处理

- 缺少 restart 时，需要 restart 的操作会被拒绝，但不会终止整个会话。
- 非法数字、选项、元素和确认输入会给出明确错误，不显示 Python traceback。
- EOF 和 `Ctrl-C` 会干净退出。
- 加载和导出失败会被明确报告，不会显示为成功的部分导出。

## 引用

本功能直接调用 calorine 的 NEP 模型修改实现。用于科研工作时，应按照 calorine 项目的要求引用其论文；同时在适当情况下引用 GPUMDkit：

- E. Lindgren, J. M. Rahm, E. Fransson, F. Eriksson, N. Österbacka, Z. Fan
  和 P. Erhart，“calorine: A Python package for constructing and sampling
  neuroevolution potential models”，*Journal of Open Source Software*
  **9**(95)，6264 (2024)，<https://doi.org/10.21105/joss.06264>。
- Z. Yan et al., *MGE Advances* **4**, e70074 (2026)，
  <https://doi.org/10.1002/mgea.70074>。
