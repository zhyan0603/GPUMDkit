# 为 GPUMDkit 贡献代码

<p align="center">
  <a href="../CONTRIBUTING.md">English</a>
  &nbsp;·&nbsp;
  <strong>简体中文</strong>
</p>

感谢你有兴趣为 `GPUMDkit` 做出贡献！我们感谢你为改善这个工具包所付出的时间和努力。`GPUMDkit` 是一个开源项目，我们欢迎社区的各种贡献，无论是修复 bug、添加新功能、改进文档还是提出建议。

---

## 目录

- [开发规范与准则](#开发规范与准则)
- [报告问题](#报告问题)
- [功能建议](#功能建议)
- [贡献代码](#贡献代码)
  - [Fork 并克隆仓库](#fork-并克隆仓库)
  - [创建功能分支](#创建功能分支)
  - [进行修改](#进行修改)
    - [交互模式贡献](#交互模式贡献)
    - [命令行模式贡献](#命令行模式贡献)
    - [更新 Tab 补全](#更新-tab-补全)
  - [代码风格与最佳实践](#代码风格与最佳实践)
  - [测试你的修改](#测试你的修改)
  - [提交信息](#提交信息)
  - [推送并创建 Pull Request](#推送并创建-pull-request)

---

## 开发规范与准则

为保持代码质量和项目一致性，请遵循以下准则：

### 语言要求

- **所有代码应使用英文编写**：包括变量名、函数名、注释、文档字符串、提交信息和文档。
- 虽然我们理解贡献者来自不同背景，但使用英文能确保代码库对最广泛的用户群体可访问。

### 模块化与可复用性

- **编写模块化和可复用的代码**：函数和脚本的设计应具有灵活性。
- **优先传递参数而非硬编码值**：使用函数参数、命令行参数或配置文件来替代硬编码值。但某些文件名由 GPUMD 和 NEP 程序固定要求（如 `train.xyz`、`thermo.out`），这些可以按需硬编码。
- **示例**：
  ```bash
  # 好的做法：接受参数
  function convert_file() {
      local input_file=$1
      local output_file=$2
      # 处理逻辑
  }
  
  # 避免：硬编码值（除非 GPUMD/NEP 要求）
  function convert_file() {
      local input_file="hardcoded_input.xyz"
      # 处理逻辑
  }
  ```

### 代码风格

- 编写清晰、易于他人维护和修改的代码。
- 使用有意义的变量和函数名。
- 在必要处添加注释以解释复杂逻辑。

### 文档

- 确保帮助信息（如 `-h` 标志）清晰准确。
- 教程文档位于 `docs/tutorials/en/`（英文）和 `docs/tutorials/zh/`（中文）。如果你添加了新功能，请考虑更新相关教程页面。
- 编辑教程 Markdown 文件后，使用 `mkdocs build -f docs/mkdocs.yml` 重新构建 HTML。

---

## 报告问题

如果你遇到 bug，请通过 GitHub Issues 帮助我们修复：

1. **搜索现有 issues** 查看该 bug 是否已被报告。
2. **创建新 issue**（如果尚未被报告）：[创建 Issue](https://github.com/zhyan0603/GPUMDkit/issues/new)
3. **详细描述问题**，如果可能请提供测试文件。

你也可以通过以下方式联系我们：
- **QQ 群**：825696376
- **邮箱**：yanzihan@westlake.edu.cn
- **脚本开发者**：联系脚本头部列出的开发者

---

## 功能建议

我们欢迎功能建议！提出新功能：

1. **搜索现有 issues** 查看是否已有人提出。
2. **创建新功能请求**：[创建功能请求](https://github.com/zhyan0603/GPUMDkit/issues/new) 并 mention @zhyan0603
3. **清楚描述你的需求** 并解释其用途。

---

## 贡献代码

### Fork 并克隆仓库

1. 在 GitHub 上 **Fork 仓库**，点击 [仓库页面](https://github.com/zhyan0603/GPUMDkit) 右上角的 "Fork" 按钮。

2. **克隆你的 fork** 到本地：

   ```bash
   git clone https://github.com/YOUR_USERNAME/GPUMDkit.git
   cd GPUMDkit
   ```

3. **设置 upstream remote** 以保持 fork 同步：
   ```bash
   git remote add upstream https://github.com/zhyan0603/GPUMDkit.git
   ```

### 创建功能分支

为你的修改创建新分支：

```bash
# 从 main 创建并切换到新分支
git checkout -b your-branch-name
```

分支名称随意。

### 进行修改

`GPUMDkit` 的结构包括：
- **`gpumdkit.sh`**：主入口（Bash 脚本），处理交互菜单模式和命令行模式
- **`install.sh`**：安装脚本，设置环境变量和 shell 配置
- **`Scripts/`**：按功能组织的 Python 和 Bash 实现脚本
  - `plt_scripts/`：绘图脚本
  - `calculators/`：计算工具
  - `format_conversion/`：格式转换工具
  - `workflow/`：工作流自动化
  - `sample_structures/`：结构采样工具
  - `analyzer/`：分析工具
  - `utils/`：实用函数（包括 `completion.sh` Tab 补全）
- **`src/`**：交互模式菜单函数的 Shell 脚本
  - `f1_format_conversions.sh`
  - `f2_sample_structures.sh`
  - `f3_workflows.sh`
  - `f4_calculators.sh`
  - `f5_analyzers.sh`
  - `f6_plots.sh`
  - `f7_utilities.sh`
- **`skills/`**：AI agent 技能定义（用于 opencode/Claude Code 集成的 SKILL.md 文件）
- **`docs/`**：文档文件
  - `tutorials/en/` 和 `tutorials/zh/`：双语教程页面
  - `mkdocs.yml`：构建教程 HTML 的 MkDocs 配置
  - `command_reference.tsv`：机器可读的命令参考
  - `htmls/`：MkDocs 生成的 HTML 输出

#### 交互模式贡献

要添加通过交互菜单访问的新功能：

1. **在相应的 `Scripts/` 子目录中创建实现脚本**：

   ```bash
   # 示例：添加新的格式转换脚本
   touch Scripts/format_conversion/new_converter.py
   ```

2. **实现你的功能**，遵循模块化准则（接受参数，避免硬编码）。

3. **在相应的 `src/` 文件中添加包装函数**：
   ```bash
   # 示例：编辑 src/f1_format_conversions.sh
   vim src/f1_format_conversions.sh
   ```

   添加新函数：
   ```bash
   function f111_new_converter(){
   echo " >-------------------------------------------------<"
   echo " | Calling the script in Scripts/format_conversion |"
   echo " | Script: new_converter.py                        |"
   echo " | Developer: Your Name (your@email.com)           |"
   echo " >-------------------------------------------------<"
   echo " Input <required_param1> <required_param2>"
   echo " Example: input.xyz output.lmp"
   echo " ------------>>"
   read -r -a converter_args
   echo " ---------------------------------------------------"
   python ${GPUMDkit_path}/Scripts/format_conversion/new_converter.py "${converter_args[@]}"
   echo " Code path: ${GPUMDkit_path}/Scripts/format_conversion/new_converter.py"
   echo " ---------------------------------------------------"
   }
   ```

4. **在 `gpumdkit.sh` 中更新菜单显示**：
   ```bash
   # 找到 menu() 函数并按需更新
   # 找到 array_choice 数组并添加你的新选择编号
   array_choice=(
       "0" "1" "101" "102" "103" "104" "105" "106" "107" "108" "109" "110" "111"  # 添加 "111"
       # ... 其余选择
   )
   ```

5. **在 `gpumdkit.sh` 中添加 case 语句**：
   ```bash
   # 在 main() 函数中找到相应的 case 语句
   case "${choice:0:1}" in
       # ...
       "1")
           source ${GPUMDkit_path}/src/f1_format_conversions.sh
           case $choice in
               "1") f1_format_conversion ;;
               "101") f101_out2xyz ;;
               # ... 现有 case ...
               "111") f111_new_converter ;;  # 添加你的 case
           esac ;;
       # ...
   esac
   ```

#### 命令行模式贡献

要添加新的命令行标志或子命令：

1. **在相应的 `Scripts/` 子目录中创建实现脚本**：
   ```bash
   # 示例：添加新的分析器脚本
   touch Scripts/analyzer/analyze_bonds.py
   ```

2. **实现你的功能**，明确参数要求：
   ```python
   # 示例：Scripts/analyzer/analyze_bonds.py
   import sys
   
   if len(sys.argv) != 3:
       print("Usage: gpumdkit.sh -analyze_bonds <input.xyz> <cutoff_distance>")
       print("Example: gpumdkit.sh -analyze_bonds structure.xyz 3.0")
       sys.exit(1)
   
   input_file = sys.argv[1]
   cutoff = float(sys.argv[2])
   
   # 你的实现
   ```

   > **注意**：GPUMDkit Python 脚本设计为直接运行（不作为模块导入）。
   > 如果你愿意，可以使用 `main()` 函数配合 `if __name__ == "__main__":` 保护，
   > 但这不是必需的。大多数现有脚本在模块级别执行。

3. **在 `gpumdkit.sh` 中添加命令行标志处理**：
   ```bash
   # 找到命令行解析部分（大的 "case $1 in" 块）
   # 在适当位置添加你的新标志
   
   case $1 in
       # ... 现有 case ...
       -analyze_bonds)
           if [ ! -z "$2" ] && [ "$2" != "-h" ] && [ ! -z "$3" ]; then
               python ${analyzer_path}/analyze_bonds.py $2 $3
           else
               echo " Usage: gpumdkit.sh -analyze_bonds <input.xyz> <cutoff_distance>"
               echo " Example: gpumdkit.sh -analyze_bonds structure.xyz 3.0"
               echo " Code path: ${analyzer_path}/analyze_bonds.py"
           fi ;;
       # ... 其余 case ...
   esac
   ```

4. **在 `gpumdkit.sh` 中更新帮助信息**：
   ```bash
   # 找到 help_info_table() 函数并添加你的命令
   function help_info_table(){
       echo "+====================================== Analysis ===============================================+"
       echo "| -analyze_bonds Analyze bond lengths in struct | -analyze_comp Analyze composition of extxyz      |"
       # ... 其余帮助表 ...
   }
   ```

#### 更新 Tab 补全

添加新的命令行标志时，更新 Tab 补全脚本：

```bash
# 编辑 Scripts/utils/completion.sh
vim Scripts/utils/completion.sh
```

将你的新标志添加到 `opts` 变量：
```bash
# 找到 local opts=... 行
local opts="-h -help -update -U -clean -time -plt -calc ... -your_new_flag ..."
#                                                        在此添加你的标志
```

如果你的标志需要文件参数，将其添加到现有的文件补全 case 中：
```bash
# 找到需要文件参数的 case 并用 | 添加你的标志
-out2xyz|-out2exyz|-...|-your_new_flag)
    COMPREPLY=($(compgen -f -- "$cur")) ;;
```

如果你的标志接受二级选项（如 `-plt` 或 `-calc`），添加新 case：
```bash
case "$prev" in
    # ... 现有 case ...
    -your_new_flag)
        COMPREPLY=($(compgen -W "option1 option2 option3" -- "$cur")) ;;
    # ... 其余 case ...
esac
```

### 代码风格与最佳实践

- **Shell 脚本**：
  - 使用 `${variable}` 进行变量展开
  - 使用 `[ ! -z "$var" ]` 进行非空检查（项目约定）
  - 交互函数应遵循 banner 格式：
    ```bash
    echo " >-------------------------------------------------<"
    echo " | Calling the script in Scripts/<category>        |"
    echo " | Script: <script_name>.py                        |"
    echo " | Developer: <Name> (<email>)                     |"
    echo " >-------------------------------------------------<"
    ```
  - 使用 `read -r -a varname` 读取多词输入，然后用 `"${varname[@]}"` 传递
  - `src/` 中的 Shell 脚本应以文件头注释块开始：
    ```bash
    # ============================================================
    # GPUMDkit <module name> module
    # Repository: https://github.com/zhyan0603/GPUMDkit
    # Author: <Name> (<email>)
    # ============================================================
    ```

- **Python 脚本**：
  - 编写清晰、可维护的代码
  - 使用有意义的变量名
  - 使用 `sys.argv` 进行参数解析（与大多数现有脚本一致）
  - 提供清晰的使用信息：`print("Usage: ...")` 和 `print("Example: ...")`
  - 如果你的脚本使用重型/特殊包（`NepTrain`、`calorine`、`dpdata`），请添加
    `print_dependency_notice()` 函数以通知用户引用建议

### 测试你的修改

提交贡献前，运行相关验证命令：

```bash
# Shell 语法检查（始终运行）
bash -n gpumdkit.sh
find src Scripts -name '*.sh' -exec bash -n {} +

# Python 语法检查（针对修改的 Python 文件）
python3 -m py_compile path/to/modified_script.py

# MkDocs 构建（如果修改了文档）
mkdocs build -f docs/mkdocs.yml

# 检查尾随空格问题
git diff --check
```

然后测试功能：

1. **测试交互模式**（如适用）：
   ```bash
   gpumdkit.sh
   # 导航到你的新功能并充分测试
   ```

2. **测试命令行模式**（如适用）：
   ```bash
   gpumdkit.sh -your_new_flag [arguments]
   gpumdkit.sh -your_new_flag -h  # 测试帮助信息
   ```

3. **测试 Tab 补全**（如果你修改了 `completion.sh`）：
   ```bash
   source Scripts/utils/completion.sh
   gpumdkit.sh -<TAB>  # 应显示你的新标志
   ```

4. **用各种输入测试**：
   - 用典型输入测试
   - 用边缘情况测试（空文件、大文件等）
   - 测试错误处理（缺少文件、无效参数）

5. **验证无回归**：确保现有功能仍然正常工作。

### 提交信息

编写清晰、描述性的提交信息来解释你做了什么更改。没有严格的格式要求——只要确保你的信息是可理解的。

### 推送并创建 Pull Request

1. **提交你的修改**：
   ```bash
   git add .
   git commit -m "你的描述性信息"
   ```

2. **保持你的分支与最新的 `dev` 分支同步**：
   ```bash
   git fetch origin
   git rebase origin/dev
   ```

3. **推送你的分支**到你的 fork：
   ```bash
   git push origin your-branch-name
   ```

4. **创建 Pull Request**：
   - 前往 [GPUMDkit 仓库](https://github.com/zhyan0603/GPUMDkit)
   - 点击 "Pull Requests" → "New Pull Request"
   - 设置目标分支为 `dev` 进行代码审查
   - 设置比较分支为你的功能分支
   - 填写 PR 描述：
     - **标题**：简要描述更改
     - **描述**：详细解释更改内容和原因
     - **测试**：描述你如何测试这些更改

5. **回应审查反馈**：
   - 对建议和建设性批评保持开放
   - 及时进行请求的修改
   - 推送额外提交到你的分支（它们会自动出现在 PR 中）

6. **批准后**：维护者将合并你的 PR。

---

## 有问题或需要帮助？

如果你对贡献有疑问或需要帮助：

- **发起讨论**：使用 [GitHub Discussions](https://github.com/zhyan0603/GPUMDkit/discussions) 提出一般性问题
- **在 issue 中提问**：在相关 issue 中评论以获取具体问题的帮助
- **联系维护者**：联系 README 中列出的核心开发者

---

再次感谢你为 `GPUMDkit` 做出贡献！你的努力帮助这个工具包为整个 GPUMD 和 NEP 社区变得更好。🚀
