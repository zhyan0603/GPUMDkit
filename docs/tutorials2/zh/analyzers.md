<div align="center">
  <h1>分析工具</h1>
  <p>
    <strong>简体中文</strong> | <a href="../en/analyzers.md">English</a>
  </p>
</div>

结构分析、验证和过滤工具。

## 可用工具

| 工具 | 命令 | 描述 |
|------|------|------|
| 成分分析 | `-analyze_comp` | 按成分分组结构 |
| 属性范围 | `-range` | 能量/力/应力统计 |
| 最小距离 | `-min_dist` | 计算最小距离（无 PBC） |
| 最小距离 (PBC) | `-min_dist_pbc` | 计算最小距离（带 PBC） |
| 电荷平衡 | `-cbc` | 检查氧化态平衡 |
| 距离过滤 | `-filter_dist_pbc` | 按最小距离过滤 |
| 盒子过滤 | `-filter_box` | 按盒子大小过滤 |
| 值过滤 | `-filter_value` | 按属性值过滤 |
| 距离范围 | `-filter_range` | 按元素对距离过滤 |
| 时间监控 | `-time` | 监控 GPUMD/NEP 进度 |

## 命令参考

### 成分分析

```bash
gpumdkit.sh -analyze_comp train.xyz
```

### 属性范围

```bash
gpumdkit.sh -range <文件> <属性> [hist]

# 示例
gpumdkit.sh -range train.xyz force
gpumdkit.sh -range train.xyz energy
gpumdkit.sh -range train.xyz force hist
```

### 最小距离

```bash
# 快速（无 PBC）
gpumdkit.sh -min_dist dump.xyz

# 精确（带 PBC）
gpumdkit.sh -min_dist_pbc dump.xyz
```

### 电荷平衡

```bash
gpumdkit.sh -cbc train.xyz
```

输出：`balanced.xyz`, `unbalanced.xyz`

### 过滤

```bash
# 按最小距离
gpumdkit.sh -filter_dist_pbc train.xyz 1.0

# 按盒子大小
gpumdkit.sh -filter_box dump.xyz 20

# 按属性值
gpumdkit.sh -filter_value train.xyz force 20

# 按元素对距离范围
gpumdkit.sh -filter_range dump.xyz Li Li 1.9 2.0
```

### 时间监控

```bash
gpumdkit.sh -time gpumd
gpumdkit.sh -time nep
```

## 依赖

| 工具 | 包 |
|------|-----|
| 所有 | ase, numpy |
| 距离 | scipy |
| 电荷平衡 | pymatgen, tqdm |
| 属性范围 | matplotlib |
