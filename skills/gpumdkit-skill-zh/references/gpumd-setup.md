# GPUMD 设置、控制和操作

使用此自成体系目录了解非系综设置命令、外部控制、最小化和 `run`。高级物理参数仍需要用户明确决定。

## 目录

- 核心设置
- 模拟盒和外部控制
- 约束和混合 MD/MC
- 立即执行的操作
- 版本处理

## 核心设置

| 关键字 | 当前参数格式 | 重要含义/约束 |
|---|---|---|
| `replicate` | `replicate <n_a> <n_b> <n_c>` | 正整数复制计数；在加载势函数之前应用 |
| `potential` | `potential <potential_filename>` | 相对/绝对文件；格式来自其内容；验证元素/类型顺序 |
| `velocity` | `velocity <initial_temperature> [seed <seed_number>]` | 温度单位为 K；种子可选。如果速度不存在且此命令被省略，GPUMD 在 300 K 初始化；建议使用用户批准的明确值 |
| `correct_velocity` | `correct_velocity <interval> [<group_method>]` | 周期性去除平移/旋转漂移；间隔必须至少为 10 步 |
| `time_step` | `time_step <dt_in_fs> [<max_distance_per_step>]` | 步长单位为 fs；文档默认值为 1 fs；可选的正位移限制器单位为 Angstrom |
| `compute_extrapolation` | `compute_extrapolation asi_file <file> gamma_low <v> gamma_high <v> check_interval <n> dump_interval <n>` | ASI 监控。默认值：`gamma_low=0`，实际上无上限的 `gamma_high`，两个间隔均为 1；在正式计算中依赖这些值之前需要用户批准 |
| `dftd3` | `dftd3 <functional> <potential_cutoff> <coordination_number_cutoff>` | 截断距离单位为 Angstrom；泛函和截断距离必须匹配预期的 D3 参数化 |
| `kspace` | `kspace <ewald|pppm>` | 倒空间静电方法；默认为 `pppm`；仅用于兼容的电荷模型 |

识别的势函数族包括 Tersoff-1988/1989/mini、EAM、FCP、LJ、NEP、NEP+ILP、SW+ILP、Tersoff+ILP 和 Deep Potential。势函数文件模式是特定于族的，不可互换。不要从此摘要创建或重写势函数文件；使用用户提供的已验证模型或获取其确切格式规范。

## 模拟盒和外部控制

| 关键字 | 当前参数格式 | 重要含义/约束 |
|---|---|---|
| `change_box` | `change_box <delta>`；`change_box <delta_xx> <delta_yy> <delta_zz>`；或 6 参数形式后跟 `<epsilon_yz> <epsilon_xz> <epsilon_xy>` | 即时模拟盒变更。长度增量单位为 Angstrom；剪切项无量纲；6 参数形式需要三斜模拟盒 |
| `deform` | `deform <A_per_step> <deform_x> <deform_y> <deform_z>` 或分量速率后跟三个标志 | 运行作用域变形；如果解析器拒绝所选形式，请确认可执行文件版本 |
| `add_force` | `add_force <group_method> <group_id> <Fx> <Fy> <Fz>` 或 `add_force <group_method> <group_id> <add_force_file>` | 向选定分组中的每个原子添加力，单位为 eV/Angstrom；文件形式提供周期性力序列 |
| `add_efield` | `add_efield <method> <group> <Ex> <Ey> <Ez> [charge|bec]` 或 `add_efield <method> <group> <field_file> [charge|bec]` | 电场单位为 V/Angstrom。默认力使用 qNEP BEC 或 `model.xyz` 电荷；显式 `charge` 使用预测/提供的电荷，而 `bec` 需要 qNEP |
| `add_spring` | `add_spring ghost_com|ghost_atom|com_com ... couple|decouple ...` | 使用以下六个完整参数格式之一；弹簧常数和参考值由用户提供 |
| `electron_stop` | `electron_stop <file>` | 加载电子阻止本领表；使用前验证表格单位/模式 |
| `plumed` | `plumed <plumed_file> <interval> <if_restart>` | 需要启用 PLUMED 的 GPUMD 构建；间隔单位为 MD 步，重启标志为 0/1 |

`add_spring` 参数格式：

```text
add_spring ghost_com  <method> <group> <gvx> <gvy> <gvz> couple   <k> <R0> <ox> <oy> <oz>
add_spring ghost_com  <method> <group> <gvx> <gvy> <gvz> decouple <kx> <ky> <kz> <ox> <oy> <oz>
add_spring ghost_atom <method> <group> <gvx> <gvy> <gvz> couple   <k> <R0> <ox> <oy> <oz>
add_spring ghost_atom <method> <group> <gvx> <gvy> <gvz> decouple <kx> <ky> <kz> <ox> <oy> <oz>
add_spring com_com <method> <group1> <group2> couple   <k> <R0>
add_spring com_com <method> <group1> <group2> decouple <kx> <ky> <kz>
```

`electron_stop` 表以 `N E_min E_max` 开头，后跟 `N` 行均匀间隔的能量行。每行包含每个势函数元素的一个阻止本领列，顺序与势函数文件中的元素顺序相同。

## 约束和混合 MD/MC

| 关键字 | 当前参数格式 | 重要含义/约束 |
|---|---|---|
| `fix` | `fix <group_label>` 或 `fix <grouping_method> <group_label>` | 冻结选定原子；单参数形式使用方法 0 |
| `move` | `move <group_id> <vx> <vy> <vz>` 或 `move <method> <group_id> <vx> <vy> <vz>` | 恒定速度，单位为 Angstrom/fs；支持 `nvt_ber`、`nvt_nhc`、`nvt_bdp` 和 `heat_lan`，不支持 NPT；与 `fix` 组合时，两者必须使用相同的分组方法 |
| `mc canonical` | `mc canonical <md_steps> <mc_trials> <T_i> <T_f> [group <method> <id>]` | 正则 MC/MD；温度单位为 K |
| `mc sgc` | `mc sgc <md_steps> <mc_trials> <T_i> <T_f> <num_species> {<species> <mu> ...} [group ...]` | 半巨正则；化学势是特定于元素的物理输入 |
| `mc vcsgc` | `mc vcsgc <md_steps> <mc_trials> <T_i> <T_f> <num_species> {<species> <phi> ...} kappa [group ...]` | 方差约束 SGC；`phi` 和 `kappa` 需要定义的热力学方案 |

## 立即执行的操作

| 关键字 | 当前参数格式 | 重要含义/约束 |
|---|---|---|
| `minimize` | `minimize <sd|fire> <force_tolerance> <max_steps> [<box_change> [<hydrostatic_strain>]]` | 容差单位为 eV/Angstrom；模拟盒优化目前需要 FIRE |
| `run` | `run <number_of_steps>` | 正步数；执行当前运行作用域设置 |

要输出最小化结构，使用零时间步、单步 NVE 块配合 `dump_xyz -1 0 1 relaxed.xyz`。使用前请加载 `gpumd-outputs.md`。

## 版本处理

这些参数格式来自本技能的内置参考。如果本地解析器拒绝某个参数格式，请捕获可执行文件版本和确切错误，然后询问匹配的版本行为。在选择任何物理值之前，无论解析器是否接受，都要询问用户。
