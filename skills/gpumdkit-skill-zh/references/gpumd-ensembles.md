# GPUMD 系综和积分器

使用此自包含参考选择正确的命令族并解释其参数。系综选择和耦合值是科学决策；请询问用户。

## 目录

- 标准系综
- 压力控制形式
- MTTK 和 QTB
- 热输运和 TTM
- PIMD、热力学积分和冲击
- 选择和验证规则

## 标准系综

```text
ensemble nve
ensemble nvt_ber <T_1> <T_2> <T_coup>
ensemble nvt_nhc <T_1> <T_2> <T_coup>
ensemble nvt_bdp <T_1> <T_2> <T_coup>
ensemble nvt_lan <T_1> <T_2> <T_coup>
ensemble nvt_bao <T_1> <T_2> <T_coup>
ensemble npt_ber <T_1> <T_2> <T_coup> {<pressure_control_parameters>}
ensemble npt_scr <T_1> <T_2> <T_coup> {<pressure_control_parameters>}
```

`T_1` 和 `T_2` 是初始/最终目标温度，单位为 K，在运行期间线性渐变。`T_coup` 以时间步单位控制恒温器耦合。不要在没有用户协议或有来源支持的方法的情况下选择它。

## 压力控制形式

标准 `npt_ber` 和 `npt_scr` 族接受以下形式之一：

```text
# 各向同性；正交盒子；所有方向周期性
<p_hydro> <C_hydro> <p_coup>

# 正交各向异性；周期性方向独立控制
<p_xx> <p_yy> <p_zz> <C_xx> <C_yy> <C_zz> <p_coup>

# 三斜；所有方向周期性；六个盒子分量
<p_xx> <p_yy> <p_zz> <p_yz> <p_xz> <p_xy> <C_xx> <C_yy> <C_zz> <C_yz> <C_xz> <C_xy> <p_coup>
```

压力和弹性模量单位为 GPa。弹性常数是粗略的转换参数，不是代理可以编造的值。超过 2000 GPa 的分量会在文档中的标准压力控制器中禁用该盒子分量的耦合。

## MTTK 和 QTB

```text
ensemble nvt_mttk temp <T_1> <T_2> [tperiod <tau_temp>]
ensemble npt_mttk temp <T_1> <T_2> <direction> <p_1> <p_2> [tperiod <tau_temp>] [pperiod <tau_press>]
ensemble nph_mttk <direction> <p_1> <p_2> [pperiod <tau_press>]

ensemble nvt_qtb <T_1> <T_2> <T_coup> [f_max <value>] [N_f <value>]
ensemble npt_qtb <direction> <p_1> <p_2> temp <T_1> <T_2> tperiod <tau_T> pperiod <tau_p> [f_max <value>] [N_f <value>]
```

MTTK 方向包括 `iso`、`aniso`、`tri`，以及单独的 `x`、`y`、`z`、`xy`、`yz`、`xz` 分量。文档中 MTTK 的默认周期为温度 100 个时间步、压力 1000 个时间步；压力周期应至少为 200 个时间步。将这些视为文档默认值，而非自动的生产选择。QTB `f_max` 是最大色噪声频率，单位为 ps^-1（默认 200，应大于系统最高声子频率）；`N_f` 是频率点数（默认 100，滤波器点数为 `2*N_f`）。QTB 压力周期必须至少为 200 个时间步。

## 热输运和 TTM

```text
ensemble heat_nhc <T> <T_coup> <delta_T> <label_source> <label_sink>
ensemble heat_bdp <T> <T_coup> <delta_T> <label_source> <label_sink>
ensemble heat_lan <T> <T_coup> <delta_T> <label_source> <label_sink>

ensemble ttm <grouping_method> <group_id> <Ce> <rho_e> <kappa_e> <gamma_p> <gamma_s> <v_0> <nx> <ny> <nz> <T_e_init> [{optional_args}]
ensemble heat_ttm <T> <T_coup> <delta_T> <label_source> <label_sink> <grouping_method> <group_id> <Ce> <rho_e> <kappa_e> <gamma_p> <gamma_s> <v_0> <nx> <ny> <nz> <T_e_init> [{optional_args}]
```

对于热源/汇系综，标签引用分组方法 0，目标温度为 `T + delta_T` 和 `T - delta_T`。TTM 单位为：`Ce` 每电子 eV/K，`rho_e` Angstrom^-3，`kappa_e` eV/(ps K Angstrom)，`gamma_p` 和 `gamma_s` amu/ps，`v_0` Angstrom/ps，`T_e_init` K。`Ce*rho_e` 的乘积是体积热容，单位为 eV/(K Angstrom^3)；`nx ny nz` 是正电子网格计数。

TTM 可选参数：

| 键 | 含义/模式 |
|---|---|
| `ttm_out_interval <n>` | 电子温度快照间隔；默认 1 步 |
| `ttm_infile <file>` | 初始网格温度；每个单元一个 `ix iy iz T_e` 行，索引从 1 开始 |
| `ttm_properties_file <file>` | 每个单元一个 `ix iy iz C_vol kappa_e gamma_p eta` 行；覆盖均匀热容、导热系数和耦合 |
| `ttm_source <value>` | 体积热源，单位为 eV/(ps Angstrom^3)；提供时乘以每个单元的 `eta` |
| `ttm_active_x|y|z <range>` | `all`、单个从 1 开始的索引或包含端点的 `3:10`/`3-10`；非活动单元保持零电子温度 |

HNEMD 不是系综关键字。使用温度控制系综配合 `compute_hnemd`；推荐使用 Nose-Hoover 链，Langevin 不适用于此目的。

## PIMD、热力学积分和冲击

| 族 | 当前签名 | 参数含义 |
|---|---|---|
| PIMD | `ensemble pimd <num_beads> <T_1> <T_2> <T_coup> [{pressure parameters}]` | 珠子数必须是不超过 128 的正偶数；在第一个 PIMD 相关运行中设置，之后不要更改 |
| RPMD | `ensemble rpmd <num_beads>` | 使用请求的珠子数进行实时环聚合物动力学 |
| TRPMD | `ensemble trpmd <num_beads>` | 使用请求的珠子数进行恒温环聚合物动力学 |
| 液体 TI | `ensemble ti_liquid temp <T> [tperiod <tau_T>] [tequil <n>] [tswitch <n>] [press <p>] sigmasqrd <v> <v>` | 恒温器周期默认为 100；省略的平衡/切换长度按 1:4 比例分配；实现的最终 `p` 值为 1、25、50、75、100 |
| 弹簧 TI | `ensemble ti_spring temp <T> [tperiod <tau_T>] [tequil <n>] [tswitch <n>] [press <p>] [spring <element> <k> ...]` | 恒温器周期默认为 100；弹簧常数单位为 eV/Angstrom^2，必须放在最后；如果缺失则从原子 MSD 估算 |
| 绝热切换 | `ensemble ti_as temp <T> [tperiod <tau_T>] <pressure_control> <pmin> <pmax> [pperiod <tau_p>] [tswitch <n>] [tequil <n>]` | 等温压力路径；省略的切换/平衡长度使用 4:1 比例 |
| 可逆缩放 | `ensemble ti_rs temp <tmin> <tmax> [tperiod <tau_T>] <pressure_control> <p> [pperiod <tau_p>] [tswitch <n>] [tequil <n>]` | 等压温度缩放路径；省略的切换/平衡长度使用 4:1 比例 |
| 固定 lambda TI | `ensemble ti lambda <lambda> temp <T> [tperiod <tau_T>] spring <element> <k> ...` | 固定耦合参数的测试积分器；恒温器周期默认为 100 |
| 壁面活塞 | `ensemble wall_piston vp <vp> [thickness <thickness>]` | 活塞速度和可选壁厚，单位为 Angstrom；壁厚默认 20 |
| 镜像壁面 | `ensemble wall_mirror vp <vp>` | 运动反射壁面速度 |
| 谐波壁面 | `ensemble wall_harmonic vp <vp> [k <k>]` | 壁面速度和可选谐波力常数，单位为 eV/Angstrom^2；`k` 默认 10 |
| MSST | `ensemble msst <x|y|z> <shock_velocity> qmass <q> mu <mu> [tscale <v>] [p0 <p0>] [v0 <v0>] [e0 <e0>]` | 冲击速度单位为 km/s；`qmass` 和 `mu` 是必需的；`tscale` 默认为 0；省略的初始状态值在第一步计算 |
| NPHug | `ensemble nphug <direction> <p_1> <p_2> [tperiod <tau_T>] [pperiod <tau_p>] [p0 <p0>] [v0 <v0>] [e0 <e0>]` | 目标压力单位为 GPa，应相等；压力周期默认为 1000；省略的初始状态值在第一步计算 |

此表记录语法和角色，而非推荐设置。这些方法需要用户批准的物理参数，以及涉及辅助文件时版本特定的文件模式。

TI 输出是特定于方法的：液体和弹簧 TI 写入 CSV 数据加上 YAML 自由能摘要；绝热切换写入压力/体积 CSV 数据；可逆缩放写入 lambda/dlambda/焓 CSV 数据。评估收敛性时保留屏幕输出和正向/反向分支。

## 选择和验证规则

- 询问盒子形状/体积是否需要固定、各向同性、各向异性或完全三斜。
- 压力控制前检查 PBC 和盒子形状限制。
- 询问温度/压力渐变是否是有意的。
- 将系综命令放在它们控制的同一 `run` 块中；它们不传播。
- 验证达到的温度/压力和盒子行为，而不仅仅是解析器成功。
- 对于输运，保留所选方法要求的动力学；不要替代恒温器。
