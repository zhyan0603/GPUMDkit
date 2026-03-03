## 主要目的

将你在 `Example-Scripts-tmp` 中已经完成且**不能修改**的正式脚本，按 GPUMDkit 现有框架进行“无侵入接入”：只做入口集成、菜单与补全、文档同步，不改脚本算法与参数行为。

# GPUMDkit 扩展计划（基于现有正式脚本）

## 1. 输入脚本现状（以 `Example-Scripts-tmp` 为准）

已确认你提供的 4 个正式脚本如下：

- `calc_neighbor_list.py`
  - 功能：从结构构建近邻并保存；
  - 关键参数：`-i/--input`、`-x/--index`、`-c/--cutoff`、`-n/--neighbor-num`、`-C/--center-elements`、`-E/--neighbor-elements`、`-d/--defect`、`-o/--output`。
- `calc_displacement.py`
  - 功能：读取轨迹 + 近邻列表，计算位移并输出文本；
  - 关键参数：`-i/--input`、`-n/--neighbor-list`、`-o/--output`、`-s/-t/-p`、`-l/--last`（与切片互斥）。
- `calc_averaged_structure.py`
  - 功能：按切片/最后若干帧求平均结构；
  - 关键参数：`-i/--input`、`-o/--output`、`-s/-t/-p`、`-l/--last`（与切片互斥）。
- `plt_plane_grid.py`
  - 功能：读取结构 + 位移文件，映射到网格并绘制平面图；
  - 关键参数：`-i/--input`、`-d/--disp`、`-e/--elements`、`-m/--tol`、`-g/--target-size`、`-o/--save-dir`、`--select-xy/xz/yz`、`-f/--frame`。

依赖统一为 `ferrodispcalc`（脚本内已有导入失败提示）。

## 2. 仓库结构与接入位点

- 主入口：`gpumdkit.sh`
  - 命令行二级分发：`-calc ...` 与 `-plt ...`
  - help 输出：`help_info_table()` 与 `plot_info_table()`
- 交互层：`src/f4_calculators.sh`
  - 适合承接三个计算类脚本（neighbor/displacement/averaged_structure）
- 实现层：
  - 计算类脚本放 `Scripts/calculators/`
  - 绘图脚本放 `Scripts/plt_scripts/`
- 补全：`Scripts/utils/completion.sh`

## 3. 需要修改/新增的文件（按“脚本不可改”原则）

## 3.1 新增文件（直接复制，不改内容）

- `Scripts/calculators/calc_neighbor_list.py`
  - 来源：`Example-Scripts-tmp/calc_neighbor_list.py`
- `Scripts/calculators/calc_displacement.py`
  - 来源：`Example-Scripts-tmp/calc_displacement.py`
- `Scripts/calculators/calc_averaged_structure.py`
  - 来源：`Example-Scripts-tmp/calc_averaged_structure.py`
- `Scripts/plt_scripts/plt_plane_grid.py`
  - 来源：`Example-Scripts-tmp/plt_plane_grid.py`

## 3.2 修改文件（仅接入，不改算法）

- `gpumdkit.sh`
  - `-calc` 新增子命令建议：
    - `nlist` -> `calc_neighbor_list.py`
    - `disp` -> `calc_displacement.py`
    - `avg-struct` -> `calc_averaged_structure.py`
  - `-plt` 新增子命令建议：
    - `plane-grid` -> `plt_plane_grid.py`
  - **参数透传策略（必须）**：
    - `gpumdkit.sh -calc nlist ...` 使用 `python .../calc_neighbor_list.py ${@:3}`
    - `gpumdkit.sh -calc disp ...` 使用 `python .../calc_displacement.py ${@:3}`
    - `gpumdkit.sh -calc avg-struct ...` 使用 `python .../calc_averaged_structure.py ${@:3}`
    - `gpumdkit.sh -plt plane-grid ...` 使用 `python .../plt_plane_grid.py ${@:3}`
    - 不做固定位置参数截断，不在外层重复实现参数校验；校验统一交给原脚本 `argparse`。
  - **`-h` 透传策略（必须）**：
    - 支持 `gpumdkit.sh -calc nlist -h` / `disp -h` / `avg-struct -h`
    - 支持 `gpumdkit.sh -plt plane-grid -h`
  - 更新 `help_info_table()`、`plot_info_table()`。

- `src/f4_calculators.sh`
  - 新增交互函数建议：
    - `f406_calc_neighbor_list`
    - `f407_calc_displacement`
    - `f408_calc_averaged_structure`
  - 更新菜单与合法输入数组。
  - 函数内部调用同样走参数透传思路，不重写算法逻辑。
  - 若保留“主菜单直达编号”风格，则同步更新 `gpumdkit.sh`：
    - `array_choice` 增加 `406/407/408`
    - `case` 分发增加 `406/407/408` 对应函数入口

- `Scripts/utils/completion.sh`
  - `-calc` 二级补全增加：`nlist disp avg-struct`
  - `-plt` 二级补全增加：`plane-grid`

## 3.3 推荐文档同步

- `Scripts/calculators/README.md`
- `Scripts/plt_scripts/README.md`
- `Scripts/README.md`

## 4. 功能流程（对应你的目标）

1. 构建近邻列表  
`gpumdkit.sh -calc nlist ...`

2. 轨迹 + 近邻 -> 位移  
`gpumdkit.sh -calc disp ...`

3. 轨迹 -> 平均结构（用于更稳定的结构参考，可选但推荐）  
`gpumdkit.sh -calc avg-struct ...`

4. 结构 + 位移 -> 平面格点图  
`gpumdkit.sh -plt plane-grid ...`

## 5. 分阶段执行计划

1. 复制 4 个正式脚本到目标目录（原样复制）。
2. 修改 `gpumdkit.sh` 完成 CLI 二级命令接入，并落实“`${@:3}` 参数透传 + 子命令 `-h` 透传”。
3. 修改 `src/f4_calculators.sh` 完成交互菜单接入；若启用主菜单直达，则同步修改 `gpumdkit.sh` 的 `array_choice/case`（`406/407/408`）。
4. 修改 `completion.sh` 完成补全。
5. 补 README 示例。
6. 冒烟测试（只测接入链路，不测算法本体），最小测试矩阵如下：
   - `gpumdkit.sh -calc nlist -h`
   - `gpumdkit.sh -calc disp -h`
   - `gpumdkit.sh -calc avg-struct -h`
   - `gpumdkit.sh -plt plane-grid -h`
   - 各子命令至少 1 条真实参数调用（确保可变参数如 `-C Pb Sr`、`--select-xy 0 1` 可透传）
   - `gpumdkit.sh` 交互模式：`4` 菜单可见并可调用新函数
   - 若启用直达编号：`406/407/408` 可直接触发目标函数
   - `completion.sh` 中 `-calc/-plt` 二级补全出现新子命令

## 6. 注意事项（基于现有正式脚本行为）

- **严格不改脚本内部逻辑**
  - 仅允许复制与入口调用；参数语义按脚本当前定义执行。
- **外层入口只做分发，不做二次解析**
  - 对 4 个新脚本统一使用 `${@:3}` 透传参数，避免破坏 `argparse` 的多参数与可选参数行为。
- **依赖准备**
  - 使用前需确保 `ferrodispcalc` 可导入。
- **`plt_plane_grid.py` 的结构帧行为**
  - 当前脚本对结构使用 `read(args.input)`（单帧读取）；`-f` 仅用于位移帧选择。
  - 多帧结构可先通过 `calc_averaged_structure.py` 生成单帧平均结构再绘图。
- **位移文件格式约束**
  - `plt_plane_grid.py` 假设位移文件可按目标元素数量整除并可重排为 `(nframe, n_ele, 3)`。
- **互斥参数**
  - `calc_displacement.py` 与 `calc_averaged_structure.py` 中，`--last` 与 `-s/-t/-p` 互斥。

## 7. 代码风格总结（仓库 + 你的新脚本）

## 7.1 仓库现有风格

- Shell 入口统一 `case` 分发，提示信息含 `Usage/Examp/Code path`。
- 交互函数以 `f<编号>_<name>` 命名，位于 `src/f*.sh`。
- Python 旧脚本多为 `sys.argv` 风格，简单直跑。

## 7.2 你提供脚本的风格

- 使用 `argparse` + 自定义 `HelpFormatter`，帮助信息更完整。
- 大多采用 `main()` 结构，参数校验更严格（互斥条件、格式检查）。
- 对外部依赖有明确 ImportError 指引（安装命令直接给出）。

## 7.3 集成策略结论

- 以“仓库入口风格”包裹“新脚本 argparse 风格”是最稳妥方案：
  - 外层保持 GPUMDkit 一致的命令体系；
  - 内层保持你脚本原始行为，不做逻辑改写。
