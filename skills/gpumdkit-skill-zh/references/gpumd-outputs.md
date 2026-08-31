# GPUMD 输出和 Dump 关键字

使用此自成体系参考选择输出、估算存储空间并防止意外的追加污染。

## 目录

- 轨迹和重启 dump
- 热力学和模型特定 dump
- Active/观察者模式
- 输出文件和写入模式
- 输出审核

## 当前轨迹和重启 dump（GPUMD 5.7 及更高版本）

| 关键字 | 当前参数格式 | 输出 | 行为 |
|---|---|---|---|
| `dump_xyz` | `dump_xyz <interval> <filename> [group <grouping_method> <group_id>] [precision <single\|double>] {<property> ...}` | 用户文件名 | 追加 extended XYZ 帧；以 `*` 结尾的文件名每帧写入一个文件；精度默认为 single |
| `dump_netcdf` | `dump_netcdf <interval> <filename> [{<optional_args> ...}]` | 用户指定的 NetCDF 文件 | 需要 NetCDF 包；坐标始终写出；可选参数包括属性、`group`、`precision <single\|double>` 和 `compression none\|deflate [<level>]` |
| `dump_restart` | `dump_restart <interval>` | `restart.xyz` | 覆盖最新的重启状态；interval 应接近总运行步数（例如只在末尾或少数几次输出） |
| `dump_beads` | 无参数 | `beads_dump_<k>.xyz` | 仅 PIMD 珠子数据 |

`dump_xyz` 始终包含包装坐标。可选属性目前包括 `mass`、`velocity`、`force`、`potential`、`virial`、`charge`、`bec`、`group_labels` 和 `unwrapped_position`。不使用 `group` 时输出整个系统；使用 `group` 时仅输出选定分组。以 `*` 结尾的文件名每帧生成一个文件。

## 旧版 dump 关键字（仅适用于 GPUMD 5.7 之前）

GPUMD 5.7 已删除 `dump_position`、`dump_velocity`、`dump_force` 和 `dump_exyz`。对于 GPUMD 5.7 及更高版本，不要生成这些关键字。它们可能出现在 5.7 之前版本的旧输入文件中；如果可执行文件版本未知，应先询问用户，再复用或翻译这些命令。

从 5.7 之前版本迁移到 5.7 及更高版本时，常见替换关系为：

```text
dump_position <i>                   -> dump_xyz <i> movie.xyz
dump_position <i> group <m> <g>     -> dump_xyz <i> movie.xyz group <m> <g>
dump_velocity <i>                   -> dump_xyz <i> dump.xyz velocity
dump_force <i>                      -> dump_xyz <i> dump.xyz force
dump_exyz <i> <v> <f> <p> <s>       -> dump_xyz <i> dump.xyz[*] velocity force potential
```

这些是迁移示例，不是当前命令。只有需要每帧独立文件时才在文件名中使用 `*`。extended XYZ 仍是支持的数据格式，`dump_observer` 仍有 `interval_exyz` 参数；不要删除这些格式/参数引用。

对于基于轨迹的 GPUMDkit MSD，请求非包装坐标或验证下游解包装对模拟盒和帧频率是否有效。

## 热力学和模型特定 dump

| 关键字 | 当前参数格式 | 输出 | 约束 |
|---|---|---|---|
| `dump_thermo` | `dump_thermo <interval>` | `thermo.out` | 正间隔，单位为 MD 步 |
| `dump_dipole` | `dump_dipole <interval>` | `dipole.out` | 需要支持偶极的 NEP 设置 |
| `dump_polarizability` | `dump_polarizability <interval>` | `polarizability.out` | 需要极化率模型 |
| `dump_shock_nemd` | `dump_shock_nemd interval <n> [bin_size <Angstrom>]` | 冲击空间热力学输出 | 分箱大小默认为 10 Angstrom |

在本技能内置的 GPUMD 版本快照中，`thermo.out` 有 18 列：`T K U Pxx Pyy Pzz Pyz Pxz Pxy ax ay az bx by bz cx cy cz`。温度单位为 K，能量单位为 eV，压力分量单位为 GPa，模拟盒矢量分量单位为 Angstrom。较旧的可执行文件可能不同；在解释不同布局之前确认其版本。

## Active 和观察者模式

```text
active <interval> <has_velocity> <has_force> <has_uncertainty> <threshold>
dump_observer <observe|average> <interval_thermo> <interval_exyz> <has_velocity> <has_force>
```

`active` 需要一个 NEP 势函数委员会，并使用第一个势函数传播 MD。它每检查一次写入 `active.out`，并将超过力不确定性阈值的结构追加到 `active.xyz`。检查选定结构是否有爆炸/非物理几何。

`dump_observer observe` 使用第一个势函数传播，并为每个模型写入 `observer0.out/.xyz`、`observer1.out/.xyz` 等。`average` 每步评估所有提供的 NEP 势函数，在其平均值上传播，并写入不带编号的 `observer.out` 和 `observer.xyz`。所有模型必须使用相同的元素顺序。

## 输出文件和写入模式

| 输出 | 生成器 | 典型写入模式 |
|---|---|---|
| `thermo.out` | `dump_thermo` | 追加 |
| 用户命名的 `.xyz` 文件 | `dump_xyz` | 追加 |
| `restart.xyz` | `dump_restart` | 覆盖 |
| `observer*.out`、`observer*.xyz` | `dump_observer` | 追加 |
| `active.out`、`active.xyz` | `active` | 追加 |
| `compute.out` | `compute` | 追加 |
| `ttm_electron_temperature.out` | `ensemble ttm` 或 `heat_ttm` | 覆盖 |
| `hac.out`、`kappa.out`、`shc.out` | 输运计算 | 追加 |
| `msd.out`、`sdc.out`、`ic.out` | 扩散计算 | 追加 |
| `rdf.out`、`adf.out`、`angular_rdf.out` | 结构计算 | 追加 |
| `mcmd.out` | `mc` | 追加 |
| `elastic.out`、`cohesive.out`、`D.out`、`omega2.out` | 即时/静态计算 | 查看匹配页面；通常覆盖 |

`ttm_electron_temperature.out` 以网格/活动范围/源元数据开头。每个快照以步长标记开头，后跟每个单元一个 `ix iy iz T_e` 行；索引从 1 开始，温度单位为 K，排序为 x 最快，然后 y，然后 z。`mcmd.out` 列为 MD 步、MC 接受率，然后是命令顺序的元素浓度。

## 输出审核

- 使用干净的运行目录或在重新运行前归档追加模式输出。
- 在正式计算前计算预期的帧/行数和估算存储空间。
- 确认 dump 属性、精度、分组和频率与下游分析匹配。
- 验证文件非空且具有文档中的列/单位。
- 在拟合前检测不完整的最终块、重复的追加数据和混合运行。
- 记录哪个 `run` 块生成了每个输出。
