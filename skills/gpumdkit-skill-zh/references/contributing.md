# 为 GPUMDkit 做贡献

在更改 GPUMDkit 代码、脚本、CLI 路由、文档或 skill 文件之前，请阅读本参考。保持编辑最小化，保留无关更改，在更改科学默认值或用户可见行为之前询问维护者。

## 目录

- 项目结构与路由
- 添加交互和 CLI 功能
- Python 和 Shell 约定
- 维护者决定与被拒绝的更改
- 验证清单

## 项目结构

```
GPUMDkit/
├── gpumdkit.sh              # 主入口：CLI 路由 + 交互菜单
├── install.sh               # 安装脚本
├── src/                     # Shell 菜单模块（交互模式包装器）
│   ├── f1_format_conversions.sh
│   ├── f2_sample_structures.sh
│   ├── f3_workflows.sh
│   ├── f4_calculators.sh
│   ├── f5_analyzers.sh
│   ├── f6_plots.sh
│   └── f7_utilities.sh
├── Scripts/                 # 实现脚本（Python + Bash）
│   ├── format_conversion/
│   ├── sample_structures/
│   ├── workflow/
│   ├── calculators/
│   ├── analyzer/
│   ├── plt_scripts/
│   └── utils/
├── skills/                  # AI agent skill 定义
├── docs/                    # 文档
│   ├── tutorials/en/        # 英文教程
│   ├── tutorials/zh/        # 中文教程
│   ├── mkdocs.yml           # MkDocs 配置
│   └── htmls/               # 生成的 HTML
```

## 路由工作原理

`gpumdkit.sh` 有两种模式：

1. **交互模式**：`gpumdkit.sh`（无参数）-> 显示编号菜单
   - 用户选择类别（1-7）-> 加载对应的 `src/fN_*.sh`
   - 用户选择功能（例如 403）-> 调用函数（例如 `f403_calc_descriptors`）
   - 函数提示输入，调用 `Scripts/` 中的 Python/Bash 脚本，打印代码路径

2. **CLI 模式**：`gpumdkit.sh -flag [args]` -> 通过 `case $1 in` 路由
   - 每个标志映射到 `Scripts/` 中的一个脚本
   - 帮助信息放在 `help_info_table()` 和 `calculator_help_table()` 中

## 添加新的交互模式功能

### 步骤 1：创建实现脚本

将其放在适当的 `Scripts/` 子目录中：

```bash
# Python 脚本示例：Scripts/calculators/calc_new_feature.py
```

### 步骤 2：在 `src/fN_*.sh` 中添加包装函数

遵循以下精确模式：

```bash
function f413_new_feature(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/calculators       |"
echo " | Script: calc_new_feature.py                     |"
echo " | Developer: Your Name (your@email.com)           |"
echo " >-------------------------------------------------<"
echo " Input <param1> <param2> [optional_param3]"
echo " Example: input.xyz nep.txt"
echo " ------------>>"
read_menu_choice feature_args || return 1
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/calc_new_feature.py "${feature_args[@]}"
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_new_feature.py"
echo " ---------------------------------------------------"
}
```

关键约定：
- 横幅使用 `>-----<` 和 `| ... |` 格式（内部 49 个破折号，总计 53 个字符）
- `echo " Input <...>"` 描述参数
- `echo " Example: ..."` 展示具体示例
- 单值输入使用 `read_menu_choice var || return 1`，数组输入使用 `read_menu_array arr || return 1`
- **不要使用裸 `read -p` 或 `read -r -a`** — 这些在 stdin 关闭时会挂起。辅助函数定义在 `gpumdkit.sh` 中
- `"${varname[@]}"` 正确引用传递所有参数
- 函数名格式：`f<Category><Number>_<descriptive_name>`
- **脚本横幅必须命名精确文件名**（不是缩写或错误后缀）
- 所有横幅使用统一宽度：内部 49 个破折号（第 82-86 行样式）。不一致的宽度（51/52/53/56）在修改文件时应规范化。

### 步骤 3：在 `gpumdkit.sh` 中注册

1. 将选择编号添加到 `array_choice`（约第 64 行）
2. 在 `main()` 函数的嵌套 case 语句中添加 case：
   ```bash
   "413") f413_new_feature ;;
   ```

### 步骤 4：更新菜单显示

如果是计算器，更新 `f4_calculators.sh` 中的计算器菜单显示。

## 添加新的 CLI 标志

### 步骤 1：创建实现脚本

同上，放在 `Scripts/` 中。

### 步骤 2：在 `gpumdkit.sh` 中添加 CLI 处理器

在主 `case $1 in` 块中添加新 case。保持处理器轻量：
`gpumdkit.sh` 应路由到实现脚本，并仅为大型分发器（如 `-plt` 和 `-calc`）打印广泛的模块帮助。详细的参数检查、用法文本、类型转换、文件存在检查和面向用户的错误信息属于 Python 实现脚本。

```bash
-new_flag)
    run_python_script "Your Name (your@email.com)" "${analyzer_path}/new_script.py" "${@:2}" ;;
```

关键约定：
- Python CLI 条目优先使用 `run_python_script "Author" "${path}/script.py" "${@:2}"`，以确保参数一致转发。
- 不要在 `gpumdkit.sh` 中重复详细的参数验证；在 Python 脚本中实现 `-h`/`--help`、缺少参数处理和验证。
- 仅在顶级分发决策（`-plt <type>`、`-calc <type>`）或无法自行验证的旧版 shell 脚本中保留 shell 侧检查。
- 为可读性按字母顺序放置在现有标志中

> 当脚本同时需要交互菜单条目时，还需在 `src/f1_format_conversions.sh` 中添加关键字转换包装函数，遵循相同的横幅模式。参考现有的 `dp2extxyz` 函数作为模板。

### 步骤 3：添加到帮助表

更新 `gpumdkit.sh` 中的 `help_info_table()` 或 `calculator_help_table()`：

```bash
echo "| -new_flag  Description of new flag           | -existing_flag    Description                   |"
```

### 步骤 4：更新 Tab 补全

编辑 `Scripts/utils/completion.sh`：

1. 将 `-new_flag` 添加到 `opts` 字符串（第 16 行）
2. 如果需要文件补全，将其添加到 `|` 分隔的 case（第 38 行）：
   ```bash
   -out2xyz|-...|-new_flag)
       COMPREPLY=($(compgen -f -- "$cur")) ;;
   ```
3. 如果需要二级选项，添加新的 case 块

## Python 脚本约定

### 文件头（完整文档字符串）

超过 50 行的脚本必须有模块级文档字符串，包含所有部分：

```python
"""
=============================================================================
GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP
Repository: https://github.com/zhyan0603/GPUMDkit
Citation: Z. Yan et al., GPUMDkit: A User-Friendly Toolkit for GPUMD and NEP,
          MGE Advances, 2026, 4, e70074 (https://doi.org/10.1002/mgea.70074)
=============================================================================
Script:     script_name.py
Category:   Category Name
Purpose:    One-line description of what this script does.
            Multi-line elaboration if needed.
Usage:      gpumdkit.sh -flag <param1> <param2> [optional]
            python3 script_name.py <param1> <param2> [optional]
Arguments:
  param1    Description of param1
  param2    Description of param2
  optional  Description of optional (default: value)
Output:
  <file>     (description of output)
Author:     Name (email)
Last-modified: YYYY-MM-DD
=============================================================================
"""
```

### 参数解析 — 标准模式

简单脚本使用 `sys.argv` 位置参数。使用以下精确模式（与 `traj2exyz.py`、`dp2xyz.py`、`pos2exyz.py` 等一致）：

```python
import sys

args = sys.argv[1:]
if len(args) < N or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -flag <param1> <param2>")
    print("    or: python3 script_name.py <param1> <param2>")
    print("")
    print(" Arguments:")
    print("   param1    Description")
    print("   param2    Description")
    print("")
    print(" Example: gpumdkit.sh -flag input.xyz output.xyz")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)
```

此模式的关键要点：
- `-h`/`--help` 使用 `sys.exit(0 ...)`（成功退出），缺少参数使用 `sys.exit(1)`（错误退出）
- 使用 `args = sys.argv[1:]` 切片模式，然后 `args[0]`、`args[1]` 而非 `sys.argv[1]`、`sys.argv[2]`
- `Usage:` 行同时显示 gpumdkit.sh 调用和直接 python 调用
- **关键：将此块放在任何重量级/可选导入之前**（如 `import dpdata`、`import calorine`、`import ovito`）。普通导入（`os`、`sys`、`numpy`）先放，然后是帮助检查，再是可选导入。这确保即使可选包未安装，`-h` 也能工作。
- 检查后读取位置参数：`input_file = args[0]`、`output_file = args[1] if len(args) >= 2 else "default"`

仅在脚本有许多带复杂默认值的可选标志时使用 `argparse`（例如 `calc_neighbor_list.py`、`calc_displacement.py`）。

### 主守卫

`if __name__ == "__main__":` 对于 GPUMDkit 脚本是可选的（它们独立运行）。但是，**鼓励新脚本使用它**以提高模块级可读性。没有它的现有脚本不应被改造（维护者明确拒绝）。

### 错误处理

Python 脚本使用 `print()` 配合 `sys.exit(1)`，不要使用原始异常：

```python
# ✓ 正确 — 前导空格，f-string 提高清晰度，之后 sys.exit
if not os.path.isfile(input_file):
    print(f" Error: file '{input_file}' does not exist.")
    sys.exit(1)

# ✓ 正确 — 带前导空格的 f-string
print(f" Error: failed to load {dataset_dir}: {e}")

# ✗ 错误 — 缺少前导空格
print(f"Error: file '{input_file}' does not exist.")

# ✗ 错误 — raise Exception 而非 print + exit
raise FileNotFoundError(input_file)
```

对于捕获预期失败的 try/except 块：

```python
try:
    result = do_something()
except Exception as e:
    print(f" Error: operation failed: {e}")
    sys.exit(1)
```

### 输入验证

读取前始终验证文件/目录存在性：

```python
if not os.path.isdir(input_dir):
    print(f" Error: input directory '{input_dir}' does not exist.")
    sys.exit(1)

if not os.path.isfile(input_file):
    print(f" Error: file '{input_file}' does not exist.")
    sys.exit(1)
```

### 输出文件命名

- 输出文件名应包含足够上下文以避免冲突（例如 `filtered_Li_Li_1.8_2.0.xyz`，而非仅 `filtered.xyz`）
- 覆写现有文件无需警告（项目约定 — 这是 CLI 工具包，不是交互式软件）

### 依赖提示

如果你的脚本需要重量级/特殊包（`NepTrain`、`calorine`、`dpdata`、`ovito`），添加提示：

```python
def print_dependency_notice():
    print(" This function requires the 'calorine' package.")
    print(" If you use this function, please cite:")
    print("   Bochkov et al., npj Comput Mater 6, 170 (2020)")

# 在脚本早期调用
print_dependency_notice()
```

不要为常见包（如 `ase`、`numpy`、`pymatgen`）添加提示。

## Shell 脚本约定

### 文件头（在 `src/` 中）

```bash
# ============================================================
# GPUMDkit <module name> module
# Repository: https://github.com/zhyan0603/GPUMDkit
# Author: Your Name (your@email.com)
# ============================================================
```

### 变量检查

仅在 shell 控制流和广泛分发决策中使用 shell 变量检查。新的 Python 支持的 CLI 命令应将参数转发给 Python，让 Python 脚本处理详细验证和错误信息。

```bash
# 非空检查（项目约定）
if [ ! -z "$var" ]; then
    ...
fi

# 非 -h 标志检查
if [ ! -z "$2" ] && [ "$2" != "-h" ]; then
    ...
fi
```

### 交互菜单模式

```bash
function f413_new_feature(){
echo " >-------------------------------------------------<"
echo " | Calling the script in Scripts/calculators       |"
echo " | Script: calc_new_feature.py                     |"
echo " | Developer: Your Name (your@email.com)           |"
echo " >-------------------------------------------------<"
echo " Input <param1> <param2>"
echo " Example: input.xyz nep.txt"
echo " ------------>>"
read_menu_choice feature_args || return 1   # EOF 安全；见项目约定
echo " ---------------------------------------------------"
python ${GPUMDkit_path}/Scripts/calculators/calc_new_feature.py "${feature_args[@]}"
echo " Code path: ${GPUMDkit_path}/Scripts/calculators/calc_new_feature.py"
echo " ---------------------------------------------------"
}
```

### 新 Python 脚本的文件结构

使用 `dp2xyz.py` 或 `traj2exyz.py` 作为参考。导入顺序模式为：

```python
import os
import sys
# 重量级/可选导入在帮助检查之后

args = sys.argv[1:]
if len(args) < N or args[0] in ("-h", "--help"):
    print(" Usage: gpumdkit.sh -flag <param1> <param2>")
    print("    or: python3 script_name.py <param1> <param2>")
    print("")
    print(" Arguments:")
    print("   param1    Description")
    print("   param2    Description")
    print("")
    print(" Example: gpumdkit.sh -flag input.xyz output.xyz")
    print("")
    sys.exit(0 if args and args[0] in ("-h", "--help") else 1)

import heavy_dependency  # 仅在帮助检查之后
from ase.io import write  # 仅在帮助检查之后

def helper_functions():
    """所有辅助函数的文档字符串。"""
    ...

if __name__ == "__main__":
    input_file = args[0]
    output_file = args[1] if len(args) == 2 else "default.xyz"
    ...
```

关键规则：
- 普通导入（`os`、`sys`、`numpy`）放在帮助检查之前
- 重量级/可选导入（`dpdata`、`calorine`、`ovito`、`ase`）放在帮助检查**之后**，这样即使可选包缺失，`-h` 也能工作
- 所有辅助函数都有文档字符串
- 鼓励新脚本使用 `if __name__ == "__main__":` 守卫
- 使用 `args = sys.argv[1:]` 模式，然后访问 `args[0]`、`args[1]` 等

## 项目约定与个人偏好

这些是维护者的约定。严格遵守。

### Echo / Print 前导空格

**每个面向用户的 `echo`/`print` 必须以空格开头。**

```bash
# ✓ 正确
echo " Error: File not found."
echo " Deleting the files..."
echo " Operation canceled."

# ✗ 错误
echo "Error: File not found."
```

```python
# ✓ 正确
print(" Error: 'input.xyz' not found.")
print(" Deleting the files...")

# ✗ 错误
print("Error: 'input.xyz' not found.")
```

这适用于所有脚本类型：交互包装器、工作流脚本、实用消息、教程记录和 Python print() 输出。唯一的例外是 `" ------------>>"` 提示标记，它按设计以空格开头。

### 纯 ASCII 框线绘制

所有边框字符必须是纯 ASCII：`|`、`+`、`-`、`/`、`\`、`>`、`<`。**不要使用** Unicode 框线字符（如 `│`（U+2502））— 它们在非 UTF-8 终端上会显示为乱码。

### EOF 安全输入

所有交互输入必须通过 `gpumdkit.sh` 中定义的辅助函数：

- `read_menu_choice var || return 1` — 读取单行到变量
- `read_menu_array arr || return 1` — 读取空格分隔的标记到数组

**不要使用裸 `read -p` 或 `read -r -a`。** 当 stdin 关闭（管道、重定向、/dev/null）时，裸 read 会以空/未设置的值静默成功或永远挂起。辅助函数打印 `" Input closed. Exiting."` 并返回 1。

### 菜单验证

- 验证数组：始终命名为 `valid_menu_choices`
- 验证循环：`while ! echo "${valid_menu_choices[@]}" | grep -wq "$num_choice"; do ...; done`
- 返回主菜单：`case ... "000") menu; main ;;`

### 脚本横幅

运行脚本前显示的信息横幅遵循以下精确格式：

```bash
echo " >-------------------------------------------------<"   # 恰好 49 个破折号
echo " | Calling the script in Scripts/<category>       |"
echo " | Script: <exact_filename>                         |"   # 必须匹配实际文件名
echo " | Developer: Name (email)                          |"
echo " >-------------------------------------------------<"
```

- 横幅中的脚本名称**必须是精确文件名**（例如 `scf_batch_pretreatment_vasp.sh`，而非 `scf_batch_pretreatment.sh`）
- 宽度：内部 49 个破折号（`>` 和 `<` 之间总计 53 个字符）。在项目中保持统一。

### 提示箭头

用户输入遵循标记模式：读取输入前 `echo " ------------>>"`。这在所有交互函数中保持一致。

### 错误消息格式

Shell：`echo " Error: <message>."` — 前导空格，大写 E，冒号，句号结尾。
Python：`print(" Error: <message>."); sys.exit(1)` — 相同格式，之后始终调用 sys.exit。

### 不要提议或执行的操作

这些更改已被项目维护者**明确拒绝**。不要建议它们。

| 被拒绝 | 原因 |
|---|---|
| 将 `if __name__ == "__main__":` 改造到现有 Python 脚本中 | 脚本独立运行；不要仅为添加守卫而搅动工作脚本 |
| 将 `from pylab import *` 替换为显式导入 | 保持现状 |
| 统一 `-plt` "save" 参数位置 | 不同脚本有不同参数数量；保持灵活 |
| 统一 DPI（150 -> 300）到所有绘图脚本 | 保持每脚本的 DPI 设置 |
| 将工作流 sourced 脚本中的 `exit` 改为 `return` | 保持现状 |
| 将 `-filter_value` / `-filter_range` / `-get_volume` / `-re_atoms` 添加到帮助表 | 刻意不记录 |
| 使交互模式在一个功能后循环回菜单 | 保持单次执行；按此设计 |
| 添加颜色/ANSI 转义码 | 用户明确拒绝；保持单色 |
| 在每个教程页面添加调试示例 | 不必要；保持教程精简 |
| 随意更改 MSD 拟合范围 | 科学选择；先问维护者 |
| 移除故障排除页面 | 由维护者明确要求移除；不要重新添加 |
| 将 `des_compare` 添加到 CLI | 脚本存在但维护者选择不接入 |
| 调试后保留 `__pycache__` | 必须始终 `find . -type d -name __pycache__ -exec rm -rf {} +` |
| 将 `-get_volume` / `-re_atoms` 添加到 completion.sh | 刻意从补全和帮助表中排除 |
| 修改 plt_scripts 中的逻辑 | 除非明确请求，否则仅对绘图脚本进行外观/格式更改 |

### macOS / 跨平台注意事项

- macOS 自带 **bash 3.2**（无 `local -n` nameref）。使用 `read -a "$varname"` 配合 eval 或 herestring，而非 namerefs。
- macOS **sed 是 BSD** 版本，非 GNU。跨平台 sed：改用 `perl -i -pe`。
- **zsh 不支持 `read -a`**（仅 bash）。所有脚本在 bash 下运行。
- `timeout` 命令在 macOS 上不可用。使用后台作业或依赖 EOF 安全读取。

### 测试新脚本

对于任何新脚本，提交前运行以下两项检查：

```bash
python3 Scripts/path/to/script.py -h    # 帮助必须在无依赖时工作
python3 Scripts/path/to/script.py        # 缺少参数必须打印用法 + exit 1
```

### 特定设计决策

- **calc_ion_conductivity.py**：MSD 拟合范围为数据的 40%-80%（第 107-108 行）。Arrhenius 绘图脚本使用相同范围。未经维护者批准不要更改。
- **clean_extra_files.sh**：保留子串删除使用**全词匹配**，而非子串移除。不要退化到 `${var/pattern/}` 子串移除（它会损坏共享子串的文件名）。
- **troubleshooting.md**：仅 5 个 FAQ 项。不询问不要扩展。
- **select_max_modev.py 第 111 行**：索引映射 `atoms_list[filtered_indices[i]]` 最近从 bug（`atoms_list[i]`）修复。不要撤销。
- **charge_balance_check.py 第 118-131 行**：现在使用手动循环而非 `dict()` 来处理 3 元组错误返回。不要恢复为 `dict()`。

---

## 验证清单

提交前，运行以下检查：

```bash
# Shell 语法
bash -n gpumdkit.sh
find src Scripts -name '*.sh' -exec bash -n {} +

# Python 语法
python3 -m py_compile Scripts/path/to/your_script.py

# py_compile 后清理 __pycache__
find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null

# MkDocs 构建（如果文档有更改）
mkdocs build -f docs/mkdocs.yml

# 空白检查
git diff --check
```
