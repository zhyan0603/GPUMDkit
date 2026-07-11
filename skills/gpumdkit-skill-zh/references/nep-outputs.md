# NEP 输出、重启与模型验证

使用本参考文件解读训练/预测文件，并判断模型是否可以部署。

## 目录

- 写入行为与频率
- 损失与模型文件
- 预测/拟合图文件
- 特殊模型输出
- 验证清单

## 写入行为与频率

附带版本的写入行为：

- `loss.out` 为追加模式；
- 其他输出文件会被持续覆盖；
- 能量/力/维里/应力/偶极矩/极化率的训练/测试输出每 1000 步更新一次；
- 其他输出文件每 100 步更新一次。

在新一轮训练之前，使用干净的目录或归档旧的 `loss.out`。

## 损失与模型文件

| 文件 | 含义 |
|---|---|
| `loss.out` | 代数、总损失/正则化损失以及训练/测试 RMSE 各项 |
| `nep.txt` | 当前训练好的模型参数；可在 GPUMD 中部署 |
| `nep.restart` | 持续更新的优化状态 |
| 检查点 `nep_*.txt`/restart | 由 `save_potential` 控制的周期性模型/状态文件 |

势函数模型 `loss.out` 列：

```text
gen L_t L_1 L_2 L_e_train L_f_train L_v_train L_e_test L_f_test L_v_test
```

偶极矩模型：

```text
gen L_t L_1 L_2 L_mu_train L_mu_test
```

极化率模型：

```text
gen L_t L_1 L_2 L_alpha_train L_alpha_test
```

能量/维里 RMSE 值记录为每原子；力 RMSE 单位为 eV/Angstrom。偶极矩/极化率 RMSE 请按数据集记录的单位解读。

如果存在 `nep.restart`，优化将从中断处继续。描述符相关超参数必须与产生该重启文件的状态一致。在没有来源检查的情况下，不要重命名或复用重启文件。

## 预测与拟合图文件

| 文件 | 列/含义 |
|---|---|
| `energy_train.out`, `energy_test.out` | 2 列：预测值然后目标值，单位为 eV/atom |
| `force_train.out`, `force_test.out` | 6 列：预测 xyz 然后目标 xyz，单位为 eV/Angstrom；每个原子一行 |
| `virial_train.out`, `virial_test.out` | 12 列：预测值然后目标值 xx yy zz xy yz zx，单位为 eV/atom |
| `stress_train.out`, `stress_test.out` | 12 列：预测值然后目标值 xx yy zz xy yz zx，单位为 GPa |
| `descriptor.out` | 预测模式下请求时输出归一化描述符 |

对于没有目标维里/应力的结构，输出目标哨兵值为 `-1e6`。请显式排除哨兵值并报告排除操作；不要将它们视为物理标签。

预测模式仅评估 `train.xyz`。当以完整批次评估时，文件行顺序与源构型/原子顺序一致；保留到源帧的映射并验证行数。

## 特殊模型输出

| 文件 | 含义 |
|---|---|
| `dipole_train.out`, `dipole_test.out` | 6 列：预测 xyz 然后目标 xyz；值为逐原子归一化 |
| `polarizability_train.out`, `polarizability_test.out` | 12 列：预测值然后目标值 xx yy zz xy yz zx；逐原子归一化 |
| `charge_train.out`, `charge_test.out` | 预测的 qNEP 电荷 |
| `bec_train.out`, `bec_test.out` | 18 列：预测值然后目标值 11 12 13 21 22 23 31 32 33；每个原子一行 |

分量顺序必须与训练数据约定保持一致。如果本地可执行文件输出了不同的布局，请停止解读并先确认其版本。

## 验证清单

- 确认训练已完成或被有意中止；保留日志和重启文件。
- 使用 `references/visualization.md` 绘制损失和拟合图输出。
- 比较训练/测试误差，按元素种类、组成、相态和构型类别检查异常值。
- 检查能量偏移、力的尾部行为、维里/应力哨兵值处理和标签单位。
- 针对预期的下游可观测量评估误差，而不仅看全局 RMSE。
- 在目标 MD 工况下测试外推/稳定性，如果相关则检查短程行为/ZBL。
- 当不确定性或委员会方法需要时，运行不同的随机种子/模型。
- 记录模型文件哈希、训练数据版本、`nep.in`、可执行文件/源代码版本和验收决定。
- 在没有用户批准标准的情况下，不要声称存在通用的 RMSE 阈值或生产就绪状态。
