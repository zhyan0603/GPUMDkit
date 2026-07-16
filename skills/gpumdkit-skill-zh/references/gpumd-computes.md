# GPUMD 计算关键字

使用此自成体系参考查找计算的参数格式、输出、兼容性规则和参数角色。选择物理或收敛参数前请询问用户。

## 目录

- 通用和结构性计算
- 动力学和输运
- 声子和模态输运
- 电子和有序度计算
- 兼容性和收敛性

## 通用和结构性计算

| 关键字 | 当前参数格式 | 主要输出 |
|---|---|---|
| `compute` | `compute <grouping_method> <sample_interval> <output_interval> {<quantity>}` | `compute.out` |
| `compute_chunk` | `compute_chunk <sample_interval> <output_interval> bin/1d\|bin/2d\|bin/3d <bin parameters> {<quantity>}` | `compute_chunk.out` |
| `compute_adf` | `compute_adf <interval> <num_bins> <rc_min> <rc_max>` 或类型化多三元组形式 | `adf.out` |
| `compute_rdf` | `compute_rdf <cutoff> <num_bins> <interval>` | `rdf.out` |
| `compute_angular_rdf` | `compute_angular_rdf <cutoff> <r_bins> <angle_bins> <interval> [atom <i> <j> ...]` | `angular_rdf.out` |
| `compute_cohesive` | `compute_cohesive <scale_start> <scale_end> <direction>` | `cohesive.out` |
| `compute_elastic` | `compute_elastic <strain_value>` | `elastic.out` |
| `compute_phonon` | `compute_phonon <displacement>` | `D.out`、`omega2.out` |

对于 `compute_cohesive`，方向值 0-6 分别表示 x、y、z、xy、yz、zx 和 xyz 缩放；点数为 `(scale_end-scale_start)*1000+1`。`compute_elastic` 使用无量纲有限应变。两者都必须跟在 `potential` 之后。

`compute_phonon` 需要 `kpoints.in`；当需要复制时，将 `replicate` 放在 `run.in` 的开头，然后在此即时计算之前加载势函数。位移单位为 Angstrom。当前参数格式没有 cutoff 参数。

通用采样规则：

- `sample_interval` 是采样之间的 MD 步数。`output_interval` 是每次写入的平均采样数，因此每个 `sample_interval * output_interval` 步写入一个块。
- `compute` 接受一个或多个不同的量：`temperature`、`potential`、`force`、`virial`、`jp`（势热流）、`jk`（动能热流）和 `momentum`。输出顺序遵循命令顺序，按所选分组方法中的每个分组解析。
- `compute_chunk` 从原子当前位置进行分箱。对于每个轴使用 `<dim> lower <delta>`，其中维度为 `x`、`y` 或 `z`，原点目前仅为 `lower`，正 `delta` 是以 Angstrom 为单位的箱宽。在 2D/3D 中轴必须不同。
- `compute_chunk` 的量为 `temperature`（K）、`density/number`（Angstrom^-3）、`density/mass`（amu/Angstrom^3）、`vx vy vz`（Angstrom/fs）和 `fx fy fz`（eV/Angstrom）。每个输出块每个分箱一行：从零开始的 ID、分箱中心坐标、平均原子数，然后是请求的值。
- 结构截断距离、分箱数、应变振幅、声子位移和元素/类型选择会改变分辨率或物理内容。请询问用户而不是复制示例值。

## 动力学和输运

| 关键字 | 当前参数格式 | 主要输出 |
|---|---|---|
| `compute_msd` | `compute_msd <sample_interval> <Nc> [group <method> <id> \| all_groups <method>] [save_every <interval>]` | `msd.out`、快照 |
| `compute_sdc` | `compute_sdc <sample_interval> <Nc> [group <method> <id>]` | `sdc.out` |
| `compute_ic` | `compute_ic <sample_interval> <Nc> <type_index> <charge>` | `ic.out` |
| `compute_hac` | `compute_hac <sampling_interval> <correlation_steps> <output_interval>` | `hac.out` |
| `compute_hnemd` | `compute_hnemd <output_interval> <Fe_x> <Fe_y> <Fe_z>` | `kappa.out` |
| `compute_hnemdec` | `compute_hnemdec <drive_type> <output_interval> <Fe_x> <Fe_y> <Fe_z>` | `onsager.out` |
| `compute_shc` | `compute_shc <sample_interval> <Nc> <0\|1\|2> <num_omega> <max_omega> [group <method> <id>]` | `shc.out` |
| `compute_viscosity` | `compute_viscosity <sampling_interval> <correlation_steps>` | `viscosity.out` |

解释规则：

- `sample_interval` 单位为 MD 步；相关时间为 `sample_interval * Nc * time_step`。
- `Nc` 是最大相关采样数。第一个完整的相关输出至少需要 `sample_interval * Nc` 个 MD 步。
- `compute_msd` 输出包含每个选定分组的时间、方向性 MSD 和方向性 SDC。
- `compute_ic` 使用从零开始的势函数类型索引和用户提供的离子电荷；验证元素/类型顺序并说明分析中使用的离子电导率定义。
- `compute_hac` 在 EMD 正式采样轨迹中运行，通常是在平衡后进行 NVE。
- `compute_hnemd` 需要温度控制；使用一个非零驱动力分量，除非方法明确要求其他方式。
- HNEMD/HNEMDEC 驱动力分量单位为 Angstrom^-1。对于 `compute_hnemdec`，`drive_type=0` 选择热驱动力，力单位为 Angstrom^-1；正整数 `i` 为 `nep.txt` 头中的第 i 个元素选择扩散驱动力，力单位为 eV/Angstrom。Langevin 恒温器与 HNEMD 和 HNEMDEC 动力学都不兼容。
- `compute_shc` 要求 `1 <= sample_interval <= 10`、`100 <= Nc <= 1000`，方向 0/1/2 对应 x/y/z。`max_omega` 单位为 THz。`group <method> -1` 计算每个非零分组，可能计算量很大。
- 驱动力、相关长度和正式采样时长需要与用户一起进行收敛测试。

## 声子和模态输运

| 关键字 | 当前参数格式 | 主要输出 |
|---|---|---|
| `compute_dos` | `compute_dos <sample_interval> <Nc> <omega_max> [group <method> <id>] [num_dos_points <n>]` | `mvac.out`、`dos.out` |
| `compute_gkma` | `compute_gkma <sample_interval> <first_mode> <last_mode> <bin_size\|f_bin_size> <size>` | `heatmode.out` |
| `compute_hnema` | `compute_hnema <sample_interval> <output_interval> <Fe_x> <Fe_y> <Fe_z> <first_mode> <last_mode> <bin_option> <size>` | `kappamode.out` |

GKMA/HNEMA 需要匹配的 `eigenvector.in`；模式索引和分箱必须与该文件匹配。`compute_gkma` 和 `compute_hnema` 不能在同一运行中同时激活；后定义的生效。模态计算可能需要大量内存和磁盘空间。
对于 HNEMA，`sample_interval` 必须整除 `output_interval`。`bin_size` 按固定模式数分组；`f_bin_size` 按频率宽度分组模式。

## 电子、偶极和有序度计算

| 关键字 | 当前参数格式 | 主要输出 |
|---|---|---|
| `compute_dpdt` | `compute_dpdt <sampling_interval>` | `dpdt.out` |
| `compute_lsqt` | `compute_lsqt <direction> <num_moments> <num_energies> <E_1> <E_2> <E_max>` | `lsqt_dos.out`、`lsqt_velocity.out`、`lsqt_sigma.out` |
| `compute_orientorder` | `compute_orientorder <interval> <cutoff\|nnn> <mode_value> <ndegrees> <l...> [<average> <wl> <wlhat>]` | `orientorder.out` |

对于 `compute_orientorder`，`cutoff` 使用距离，`nnn` 使用邻居数。可选的 `average`、`wl` 和 `wlhat` 标志为 0/1，控制邻居平均和三阶不变量。

LSQT 接受输运方向 `x`、`y` 或 `z`；能量单位为 eV。它是初步的，目前仅支持硬编码的碳紧束缚模型。`E_max` 必须略大于模型的最大绝对能量，需要用户批准的收敛程序。

## 兼容性和收敛性

- `compute_dos` 和 `compute_sdc` 不能在同一运行中使用。
- 对于 `compute_dos`，角频率最大值单位为 THz，`num_dos_points` 默认为 `Nc`，`group <method> -1` 请求所有分组。
- 计算关键字通常只在当前 `run` 段生效，不会自动延续；每个正式采样段都应重新给出所需命令。
- 确保 `Nc * sample_interval` 明显短于正式采样阶段时长，并在适用时使用保存的窗口/副本验证收敛性。
- 重新运行前确认输出追加行为。
- 如果输出列或解析器行为与此快照不同，请停止并在解释数据前确认可执行文件版本。
