# NEP 输入参数

使用本自成体系目录查阅当前 `nep.in` 语法、约束和文档默认值。在更改影响模型或拟合目标的参数之前，请先确认。

## 目录

- 文件规则与最小输入
- 模式与元素种类参数
- 描述符参数
- 损失参数
- 优化与输出参数
- 微调与冲突检查

## 文件规则与最小输入

每行格式为 `keyword parameter...`；空行和 `#` 后的文字会被忽略。关键字可以按任意顺序出现，但 `type_weight` 必须跟在 `type` 之后。`type` 是必需的；其余关键字均有文档默认值。

最小的解析器级训练输入可以只包含 `type`，但正式使用的模型需要用户明确审查架构、损失、优化器、检查点和特殊物理设置。

## 模式与元素种类参数

| 关键字 | 语法/值 | 文档默认值/约束 |
|---|---|---|
| `version` | `version <integer>` | 本内置版本快照仅支持 `4` |
| `prediction` | `prediction 0\|1` | `0` 优化；`1` 预测 `train.xyz`；默认 0 |
| `model_type` | `model_type 0\|1\|2` | `0` 势函数，`1` 偶极矩，`2` 极化率；默认 0 |
| `type` | `type <number_of_species> <species...>` | 必需；区分大小写的周期表符号 |
| `type_weight` | `type_weight <weight_for_each_species...>` | 按 `type` 顺序为每个元素提供一个非负力损失权重；默认各为 1.0；必须跟在 `type` 之后 |
| `charge_mode` | `charge_mode 0\|1\|2` | 0 原始；1 实空间+倒空间 qNEP；2 仅倒空间 qNEP；默认 0 |
| `atomic_v` | `atomic_v 0\|1` | 0 全局默认；1 偶极矩/极化率使用原子级；势函数维里仅支持全局 |

## 描述符参数

| 关键字 | 语法 | 文档默认值/约束 |
|---|---|---|
| `zbl` | `zbl <outer_cutoff>` | 未设置表示关闭；截断距离 1-3 Angstrom；内截断距离为外截断的一半 |
| `use_typewise_cutoff_zbl` | `use_typewise_cutoff_zbl [<factor>]` | 未设置表示关闭；无量纲因子默认 0.7；启用模式将内 ZBL 截断设为零 |
| `cutoff` | `cutoff <radial> <angular>` 或每元素一对径向/角向值 | 默认 8,4 Angstrom；`3 <= angular <= radial <= 10`；跨元素截断为算术平均值 |
| `n_max` | `n_max <radial> <angular>` | 默认 6,6；径向 0-12，角向 0-8 |
| `basis_size` | `basis_size <radial> <angular>` | 默认 6,6；各为 0-8 |
| `l_max` | `l_max <l3> [q222 q1111 q112 q123 q233 q134]` | 默认 `4 1` = `4 1 0 0 0 0 0`；`l3` 2-8；开关为 0/1 |
| `neuron` | `neuron <number>` | 默认 30；范围 1-120 |

ZBL 说明：

- `zbl.in` 支持灵活的特定对参数；需要 `n(n+1)/2` 行按对排序的数据，每对 10 个值。
- 即使 `zbl.in` 提供了自己的截断距离，`zbl` 的截断参数在语法上仍然是必需的。
- `q1111` 仅为向后兼容保留，强烈不推荐使用；较新的可选四体开关尚未经过全面测试。

## 损失参数

| 关键字 | 语法 | 文档默认值/约束 |
|---|---|---|
| `lambda_1` | `lambda_1 <weight>` | 非负；默认 `sqrt(N_parameters)/1000` |
| `lambda_2` | `lambda_2 <weight>` | 非负；默认 `sqrt(N_parameters)/1000` |
| `lambda_e` | `lambda_e <weight>` | 非负；默认 1.0 |
| `lambda_f` | `lambda_f <weight>` | 非负；默认 1.0 |
| `lambda_v` | `lambda_v <weight>` | 非负；默认 0.1 |
| `lambda_q` | `lambda_q <weight>` | 仅 qNEP；非负；默认 0.1 |
| `lambda_z` | `lambda_z <weight>` | 仅 qNEP；非负；默认 0.5 |
| `lambda_shear` | `lambda_shear <weight>` | 额外剪切维里乘数；有效剪切权重为 `lambda_v*lambda_shear`；非负；默认 1 |
| `force_delta` | `force_delta <delta>` | 非负，单位为 eV/Angstrom；默认 0；正值会强调较小目标力上的误差 |

损失权重改变拟合目标，不是通用的调节旋钮。在更改之前，确认可用目标、预期可观测量和归一化方式。

## 优化与输出参数

| 关键字 | 语法 | 文档默认值/约束 |
|---|---|---|
| `batch` | `batch <batch_size>` | 整数 >=1；默认 1000 |
| `population` | `population <size>` | 整数 10-100；默认 50 |
| `generation` | `generation <count>` | 整数 0 到 10^7；默认 100000 |
| `save_potential` | `save_potential <interval> <format> <save_restart>` | 间隔默认 100000；格式 0 使用代数命名，1 使用时间戳/扩展命名（默认）；重启标志为 0/1 |
| `output_descriptor` | `output_descriptor 0\|1\|2` | 仅预测模式；0 关闭，1 逐结构，2 逐原子；默认 0 |

对于 `save_potential`，`save_restart=1` 表示请求保存匹配的重启检查点。已保存的 `nep.restart` 是重启/基础模型工作流所必需的。如果本地检查点名称不同，请在自动化收集之前确认可执行文件版本。

## 微调与冲突检查

```text
fine_tune <nep_model_file> <nep_restart_file>
type <number_of_types> <species...>
```

微调会将 `version`、`zbl`、`cutoff`、`n_max`、`basis_size`、`l_max` 和 `neuron` 固定为源模型架构。检查模型/重启文件来源；不要替换来自其他模型的值。

执行前：

- 比较所有元素种类和架构值与模型/重启文件头部；
- 确保预测/微调/重启模式不会意外混用；
- 确保理解 qNEP 设置和 k 空间部署需求；
- 记录所有显式值和隐式默认值；
- 在正式计算中使用文档默认值之前请先确认。
