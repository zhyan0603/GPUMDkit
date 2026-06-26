<div align="center">
  <h1>结构采样</h1>
  <p>
    <strong>简体中文</strong> | <a href="../en/sampling.md">English</a>
  </p>
</div>

结构选择和采样工具。

## 可用方法

| 方法 | 菜单 | 描述 |
|------|------|------|
| 均匀/随机 | 201 | 选择均匀间隔或随机帧 |
| FPS by NepTrain | 203 | 最远点采样（推荐） |
| 扰动 | 204 | 生成扰动结构 |
| 力偏差 | 205 | 选择高力偏差结构 |

## 均匀/随机采样

```bash
python Scripts/sample_structures/sample_structures.py <input.xyz> <method> <num> [skip]

# 示例
python Scripts/sample_structures/sample_structures.py dump.xyz uniform 50
python Scripts/sample_structures/sample_structures.py dump.xyz random 100 500
```

输出：`sampled_structures.xyz`

## 最远点采样 (FPS)

```bash
python Scripts/sample_structures/neptrain_select_structs.py <sample.xyz> <train.xyz> <nep.txt>

# 示例
python Scripts/sample_structures/neptrain_select_structs.py dump.xyz train.xyz nep.txt
```

输出：`selected.xyz`, `select.png`

## 结构扰动

```bash
python Scripts/sample_structures/perturb_structure.py <input.vasp> <num> <cell_pert> <atom_pert> <style>

# 示例
python Scripts/sample_structures/perturb_structure.py POSCAR 20 0.03 0.2 uniform
```

参数：
- `num`：扰动结构数量
- `cell_pert`：晶胞扰动比例（如 0.03 = 3%）
- `atom_pert`：原子扰动距离（埃）
- `style`：`normal`, `uniform` 或 `const`

输出：`POSCAR_01.vasp`, `POSCAR_02.vasp`, ...

## 力偏差选择

```bash
python Scripts/sample_structures/select_max_modev.py <top_n> <min_deviation>

# 示例
python Scripts/sample_structures/select_max_modev.py 200 0.15
```

需要文件：`active.out`, `active.xyz`

输出：`selected.xyz`

## 帧范围提取

```bash
python Scripts/sample_structures/frame_range.py <input.xyz> <start> <end>

# 示例：前 80% 的帧
python Scripts/sample_structures/frame_range.py dump.xyz 0 0.8
```

## 依赖

| 方法 | 包 |
|------|-----|
| 均匀/随机 | numpy, ase |
| FPS | numpy, ase, matplotlib, scikit-learn, scipy, NepTrain |
| 扰动 | dpdata |
| 力偏差 | numpy, ase |
