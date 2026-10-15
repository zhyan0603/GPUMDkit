# 分析器

## 目录

- 可用工具
- 命令参考
- 常用工作流
- CLI 标志与依赖

## 可用工具

| 工具 | 命令 | 描述 |
|------|---------|-------------|
| 成分分析 | `-analyze_comp` | 按成分分组结构 |
| 离群值检测 | 菜单 502 | 查找高误差结构 |
| 化学种类 | `-chem_species` 或菜单 503 | 列出唯一元素 |
| 电荷平衡 | `-cbc` | 检查氧化态平衡 |
| 性质范围 | `-range` | 能量/力/维里统计 |
| 距离过滤 | `-filter_dist` 或菜单 506 | 按最小距离过滤（无 PBC） |
| 距离过滤（PBC） | `-filter_dist_pbc` | 按最小距离过滤（含 PBC） |
| 最小距离 | `-min_dist` | 计算最小距离（无 PBC） |
| 最小距离（PBC） | `-min_dist_pbc` | 计算最小距离（含 PBC） |
| 概率密度 | `-pda` 或菜单 508 | 3D 扩散通道分析 |
| 能量参考平移 | `-shift_energy` 或菜单 509 | 交互式平移 extxyz 能量参考 |

## 命令参考

### 成分分析
```bash
# 分析 extxyz 文件的成分
gpumdkit.sh -analyze_comp train.xyz

# 输出：唯一成分及计数表
# 交互式：选择要导出为单独文件的成分
```

### 性质范围分析
```bash
# 分析能量范围
gpumdkit.sh -range train.xyz energy

# 分析力范围
gpumdkit.sh -range train.xyz force

# 分析维里范围
gpumdkit.sh -range train.xyz virial

# 带直方图
gpumdkit.sh -range train.xyz force hist
```

能量和维里从 `energy=`/`virial=` 读取，键名不区分大小写，ASE 3.23+ 存入 calculator 的键也能读到。缺少该性质的帧会被跳过并给出提示；只有所有帧都缺少时才报错。能量为文件中的总能量（eV）；力的范围按每个原子的力的大小统计；维里范围包含所有存储的分量。

### 能量参考平移
```bash
# 打开交互式能量参考平移页面
gpumdkit.sh -shift_energy

# 也可以从菜单 5 -> 509 进入
```

`align_energy_reference.py` 会保留原子行和其他 header 字段，只更新
`energy=` 值。支持参考组对齐、零原子基线对齐以及 DFT 到 NEP 对齐；后一
模式需要 `calorine` 和兼容的 NEP 模型文件。

### 最小距离
```bash
# 快速计算（无 PBC）
gpumdkit.sh -min_dist dump.xyz

# 精确计算（含 PBC）
gpumdkit.sh -min_dist_pbc dump.xyz

# 输出：所有元素对的最小距离表
```

两个版本都按元素在所有帧中首次出现的顺序收集元素，并跳过某帧中不存在的元素对。

### 电荷平衡检查
```bash
# 检查电荷平衡
gpumdkit.sh -cbc train.xyz

# 输出：
# - balanced.xyz（电荷平衡的结构）
# - unbalanced.xyz（不平衡的结构）
# - indices.txt（摘要）
```

### 化学种类识别
```bash
# 列出唯一元素
gpumdkit.sh -chem_species train.xyz

# 交互模式也可以：
# gpumdkit.sh  # 选择：5) Analyzer -> 503

# 输出：文件中所有元素的排序列表
```

### 结构过滤

#### 按最小距离（无 PBC）
```bash
gpumdkit.sh -filter_dist dump.xyz <min_distance>
# 输出：filtered_<输入文件名>、filtered_out_<输入文件名>
# 例如 filtered_dump.xyz 和 filtered_out_dump.xyz
```

#### 按最小距离（含 PBC）
```bash
gpumdkit.sh -filter_dist_pbc dump.xyz <min_distance>
# 对周期性系统更精确；输出文件名与 -filter_dist 相同
```

#### 按元素对距离范围
```bash
# 过滤 Li-Li 距离在 1.9 到 2.0 埃之间的结构
gpumdkit.sh -filter_range dump.xyz Li Li 1.9 2.0

# 输出：filtered_<elem1>_<elem2>_<min>_<max>.xyz
```

#### 按模拟盒尺寸
```bash
# 丢弃任一模拟盒边长超过 <max_edge> 埃的结构
gpumdkit.sh -filter_box dump.xyz <max_edge>

# 输出：filtered_by_box.xyz
```

#### 按性质值
```bash
# 保留所有值都不超过阈值的结构
gpumdkit.sh -filter_value train.xyz force <max_force_component>

# 输出：filtered.xyz
```

所有值都不超过阈值的结构才会保留：`force` 比较每个力分量的绝对值（eV/埃），`energy` 比较文件中的总能量（eV，不是每原子能量），`virial` 比较每个维里分量。能量和维里的键名不区分大小写，ASE 3.23+ 存入 calculator 的键也能读到。没有存储该性质的帧无法检查，会直接保留，命令会报告其数量。没有任何帧含该性质、性质名无法识别或阈值不是数字时，命令报错且不写出文件。

### 离群值检测
```bash
# 基于 RMSE 查找离群结构
gpumdkit.sh  # 选择：5) Analyzer -> 502

# 非交互式阈值模式也可以
python3 ${GPUMDkit_path}/Scripts/analyzer/find_outliers.py <energy_rmse_meV_atom> <force_rmse_meV_A> <stress_rmse_GPa>

# 当前目录中需要的文件：
# - train.xyz
# - energy_train.out
# - force_train.out
# - stress_train.out

# 输出：selected.xyz（高误差）、remained.xyz（低误差）
```

三个可选阈值分别是能量（meV/atom）、力（meV/埃）和应力（GPa）的目标 RMSE。对每种性质，脚本逐个剔除误差最大的结构，直到剩余结构的 RMSE 低于阈值；被剔除的结构写入 `selected.xyz`。阈值需与用户确认。

含有 `-1e6` 应力缺失标记，或任一应力分量绝对值超过 `1e6` 的结构，会从 RMSE 计算和两个输出数据集中排除。若没有剩余的有效结构，脚本会报错且不写出结果。

### 概率密度分析
```bash
# 计算移动离子的 3D 概率密度
gpumdkit.sh -pda LLZO.vasp dump.xyz Li 0.25

# 参数：
# - 参考结构（POSCAR）
# - 轨迹文件（extxyz）
# - 移动物种（例如 Li）
# - 网格间距（例如 0.25 埃）

# 输出：probability_density_<interval>.vasp（CHGCAR 格式）
# 可用 VESTA 或类似软件可视化
```

## 常用工作流

### 数据质量检查
```bash
# 1. 检查成分
gpumdkit.sh -analyze_comp train.xyz

# 2. 检查最小距离
gpumdkit.sh -min_dist_pbc train.xyz

# 3. 检查力范围
gpumdkit.sh -range train.xyz force

# 4. 检查离群值
gpumdkit.sh  # 选择：5) Analyzer -> 502
```

### 结构过滤流程

各阈值与体系相关，需由用户提供。

```bash
# 1. 按距离过滤
gpumdkit.sh -filter_dist_pbc dump.xyz <min_distance>

# 2. 按模拟盒尺寸过滤
gpumdkit.sh -filter_box filtered_dump.xyz <max_edge>

# 3. 按力值过滤
gpumdkit.sh -filter_value filtered_by_box.xyz force <max_force_component>
```

### 扩散通道分析
```bash
# 1. 计算概率密度
gpumdkit.sh -pda LLZO.vasp dump.xyz Li 0.25

# 2. 用 VESTA 或类似软件可视化
# 打开 probability_density_0.25.vasp
```

## 附加工具

### 时间监控
```bash
# 监控 GPUMD 进度
gpumdkit.sh -time gpumd

# 监控启用 WPE 的 GPUMD 运行
gpumdkit.sh -time wpe

# 监控 NEP 训练进度
gpumdkit.sh -time nep
```

## CLI 标志参考

| 标志 | 描述 | 语法 |
|------|-------------|--------|
| `-analyze_comp` | 成分分析 | `gpumdkit.sh -analyze_comp <file>` |
| `-range` | 性质范围 | `gpumdkit.sh -range <file> <prop>` |
| `-shift_energy` | 交互式能量参考平移 | `gpumdkit.sh -shift_energy` |
| `-min_dist` | 最小距离（无 PBC） | `gpumdkit.sh -min_dist <file>` |
| `-min_dist_pbc` | 最小距离（PBC） | `gpumdkit.sh -min_dist_pbc <file>` |
| `-chem_species` | 化学种类列表 | `gpumdkit.sh -chem_species <file>` |
| `-cbc` | 电荷平衡 | `gpumdkit.sh -cbc <file>` |
| `-filter_dist` | 按最小距离过滤 | `gpumdkit.sh -filter_dist <file> <min_dist>` |
| `-filter_dist_pbc` | 按最小距离过滤（PBC） | `gpumdkit.sh -filter_dist_pbc <file> <min_dist>` |
| `-filter_range` | 按距离范围过滤 | `gpumdkit.sh -filter_range <file> <e1> <e2> <min> <max>` |
| `-filter_box` | 按模拟盒尺寸过滤 | `gpumdkit.sh -filter_box <file> <limit>` |
| `-filter_value` | 按性质过滤 | `gpumdkit.sh -filter_value <file> <prop> <thresh>` |
| `-time` | 时间监控 | `gpumdkit.sh -time <gpumd\|wpe\|nep\|gnep>` |

## 依赖

| 工具 | 所需包 |
|------|------------------|
| 所有 Python 脚本 | `ase`、`numpy` |
| 距离计算 | `scipy` |
| 电荷平衡 | `pymatgen`、`tqdm` |
| 性质范围 | `matplotlib` |
| 概率密度 | `pymatgen` |
| DFT 到 NEP 能量平移 | `calorine` |

## 详细文档

用户指南请参见 `${GPUMDkit_path}/docs/tutorials/zh/分析工具.md`。
