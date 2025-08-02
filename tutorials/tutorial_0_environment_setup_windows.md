# 🎯 教程目标

这是一个专门为**完全零基础的初学者**准备的Windows系统环境安装教程。我们将从最基础的软件下载开始，一步一步地搭建完整的ICellbioRpy使用环境。

**完成本教程后，您将拥有：**
- ✅ R和RStudio工作环境
- ✅ Anaconda Python环境
- ✅ 专用的conda环境
- ✅ ICellbioRpy包及其依赖
- ✅ 完全配置好的Python-R交互环境

**预计时间：** 1-2小时（取决于网络速度）

# 📋 开始前的准备

## 检查您的系统

1. **确认您使用的是Windows系统**（本教程适用于Windows 10/11）
2. **确保您有管理员权限**（安装软件时可能需要）
3. **确保网络连接稳定**（需要下载较大的安装文件）
4. **准备至少5GB的磁盘空间**

## 本教程的特点

- 🔍 **超详细步骤**：每一步都有截图说明
- 🚨 **错误预防**：提前说明可能遇到的问题
- 💡 **小贴士**：实用的操作建议
- 🔧 **故障排除**：常见问题的解决方法

# 第一步：安装R语言

## 1.1 下载R

1. **打开浏览器**（Chrome、Edge、Firefox都可以）

2. **访问R官方网站**：
   - 在地址栏输入：`https://www.r-project.org/`
   - 按回车键

3. **进入下载页面**：
   - 点击页面上的 **"CRAN"** 链接
   - 选择一个中国的镜像站点（建议选择清华大学或中科院的镜像）

4. **选择Windows版本**：
   - 点击 **"Download R for Windows"**
   - 点击 **"base"**
   - 点击 **"Download R 4.x.x for Windows"**（数字可能不同，选最新版本）

💡 **小贴士**：下载文件大约70-100MB，请耐心等待。

## 1.2 安装R

1. **找到下载的文件**：
   - 通常在"下载"文件夹中
   - 文件名类似：`R-4.4.3-win.exe`

2. **右键点击安装文件**：
   - 选择 **"以管理员身份运行"**

3. **按照安装向导操作**：
   - 选择语言：**中文（简体）** 或 **English**
   - 点击 **"下一步"**
   - 阅读许可协议，点击 **"我接受"**
   - 选择安装位置：**保持默认** `C:\Program Files\R\R-4.x.x\`
   - 选择组件：**保持默认选择**
   - 启动选项：**保持默认**
   - 点击 **"安装"**

4. **等待安装完成**：
   - 安装过程大约3-5分钟
   - 看到"安装完成"提示，点击 **"完成"**

## 1.3 验证R安装

1. **按 `Win + R` 键**打开运行对话框

2. **输入 `cmd`** 并按回车，打开命令提示符

3. **输入以下命令**：
   ```cmd
   R --version
   ```

4. **如果看到类似输出，说明安装成功**：
   ```
   R version 4.4.3 (2025-02-28) -- "Trophy Case"
   Copyright (C) 2025 The R Foundation for Statistical Computing
   ```

🚨 **如果出现错误**："R不是内部或外部命令"，说明R没有添加到系统路径中。请参考本教程末尾的故障排除部分。

# 第二步：安装RStudio

## 2.1 下载RStudio

1. **访问RStudio官网**：
   - 在浏览器中输入：`https://rstudio.com/products/rstudio/download/`

2. **选择免费版本**：
   - 找到 **"RStudio Desktop"**
   - 点击 **"Download"** 按钮

3. **下载Windows版本**：
   - 系统会自动检测您的操作系统
   - 点击 **"Download RStudio Desktop for Windows"**

💡 **小贴士**：RStudio安装文件约200MB，比R大一些。

## 2.2 安装RStudio

1. **运行安装文件**：
   - 找到下载的文件（类似`RStudio-2024.xx.x-xxx.exe`）
   - 双击运行（如果提示需要管理员权限，选择"是"）

2. **按照向导安装**：
   - 点击 **"Next"**
   - 选择安装位置：**保持默认**
   - 点击 **"Install"**
   - 等待安装完成，点击 **"Finish"**

## 2.3 首次启动RStudio

1. **在桌面或开始菜单找到RStudio图标**，双击启动

2. **RStudio界面说明**：
   - **左下角**：Console（控制台）- 在这里运行R代码
   - **左上角**：Source（源代码编辑器）
   - **右上角**：Environment（环境变量）
   - **右下角**：Files/Plots/Packages（文件/图形/包管理）

3. **测试RStudio**：
   - 在Console中输入：`1 + 1`
   - 按回车，应该看到结果：`[1] 2`

✅ **如果看到正确结果，说明R和RStudio都安装成功了！**

# 第三步：安装Anaconda

## 3.1 为什么需要Anaconda？

Anaconda是一个Python发行版，它：
- 包含Python解释器
- 包含科学计算常用的库
- 提供conda包管理器
- 可以创建独立的虚拟环境

## 3.2 下载Anaconda

1. **访问Anaconda官网**：
   - 在浏览器中输入：`https://www.anaconda.com/products/distribution`

2. **选择Windows版本**：
   - 点击 **"Download"** 按钮
   - 选择 **"Windows"** 标签
   - 下载 **64-Bit Graphical Installer**

💡 **小贴士**：Anaconda安装文件很大（约500MB），请确保网络稳定。

## 3.3 安装Anaconda

1. **运行安装文件**：
   - 找到下载的文件（类似`Anaconda3-2024.xx-Windows-x86_64.exe`）
   - 右键选择 **"以管理员身份运行"**

2. **安装步骤**：
   - 点击 **"Next"**
   - 点击 **"I Agree"** 接受许可协议
   - 选择 **"Just Me (recommended)"**
   - 选择安装位置：**建议保持默认** `C:\Users\YourName\anaconda3\`
   - **重要**：在"Advanced Options"页面：
     - ✅ **勾选** "Add Anaconda3 to my PATH environment variable"
     - ✅ **勾选** "Register Anaconda3 as my default Python 3.x"
   - 点击 **"Install"**

3. **等待安装完成**：
   - 安装过程可能需要10-20分钟
   - 完成后点击 **"Next"**，然后点击 **"Finish"**

## 3.4 验证Anaconda安装

1. **打开Anaconda Prompt**：
   - 按 `Win` 键，搜索 **"Anaconda Prompt"**
   - 点击打开

2. **测试conda命令**：
   ```bash
   conda --version
   ```
   应该看到类似：`conda 24.x.x`

3. **测试Python**：
   ```bash
   python --version
   ```
   应该看到类似：`Python 3.11.x`

✅ **如果看到版本信息，说明Anaconda安装成功！**

# 第四步：创建专用conda环境

## 4.1 为什么需要专用环境？

创建专用环境的好处：
- 避免包冲突
- 更好的版本管理
- 出问题时可以重新创建
- 不影响系统默认Python

## 4.2 创建ICellbioRpy专用环境

1. **打开Anaconda Prompt**（如果还没打开的话）

2. **创建新环境**：
   ```bash
   conda create -n 1cellbio python=3.11
   ```
   - `-n 1cellbio`：环境名称为"1cellbio"
   - `python=3.11`：指定Python版本

3. **确认创建**：
   - 当提示 "Proceed ([y]/n)?" 时，输入 `y` 并按回车
   - 等待环境创建完成

4. **激活环境**：
   ```bash
   conda activate 1cellbio
   ```

5. **验证环境**：
   - 命令提示符应该变成：`(1cellbio) C:\Users\YourName>`
   - 这表示您现在在icellbio环境中

## 4.3 在新环境中安装anndata

1. **确保在icellbio环境中**（命令提示符前有`(1cellbio)`）

2. **安装anndata**：
   ```bash
   conda install -c conda-forge anndata
   ```

3. **确认安装**：
   - 当提示时输入 `y`
   - 等待安装完成

4. **验证anndata安装**：
   ```bash
   python -c "import anndata; print('anndata version:', anndata.__version__)"
   ```
   应该看到：`anndata version: 0.11.x`

✅ **如果看到版本信息，说明conda环境和anndata都安装成功！**

# 第五步：安装R包和依赖

## 5.1 在RStudio中安装基础R包

1. **打开RStudio**

2. **安装基础包**：
   在Console中逐行运行以下代码：
   ```r
   # 安装基础包
   install.packages(c("devtools", "reticulate", "Matrix", "data.table"))
   ```

3. **等待安装完成**：
   - 这可能需要几分钟
   - 如果提示选择镜像，选择中国的镜像

4. **安装Bioconductor相关包**：
   ```r
   # 安装BiocManager
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   
   # 安装生物信息学相关包
   BiocManager::install(c("SingleCellExperiment", "SummarizedExperiment", "S4Vectors"))
   ```

## 5.2 安装ICellbioRpy包

**注意**：这里我们假设ICellbioRpy包已经发布到GitHub。如果还没有，请跳到第六步先配置Python环境。

```r
# 从GitHub安装ICellbioRpy包
devtools::install_github("your_username/ICellbioRpy")
```

如果遇到网络问题，可以尝试：
```r
# 设置代理（如果需要）
# options(download.file.method = "wininet")

# 或者从本地安装（如果您有源码包）
# devtools::install_local("path/to/ICellbioRpy")
```

# 第六步：配置Python-R交互环境

## 6.1 在RStudio中配置reticulate

1. **加载reticulate包**：
   ```r
   library(reticulate)
   ```

2. **查看可用的conda环境**：
   ```r
   conda_list()
   ```
   您应该能看到刚才创建的"1cellbio"环境

3. **指定使用icellbio环境**：
   ```r
   use_condaenv("1cellbio", required = TRUE)
   ```

4. **验证Python配置**：
   ```r
   py_config()
   ```
   检查输出中的Python路径是否指向icellbio环境

## 6.2 测试ICellbioRpy的Python环境配置

1. **加载ICellbioRpy包**：
   ```r
   library(ICellbioRpy)
   ```

2. **配置Python环境**：
   ```r
   configure_python_env(conda_env = "1cellbio", verbose = TRUE)
   ```

3. **检查anndata可用性**：
   ```r
   check_anndata_available()
   ```

✅ **如果看到 "✓ anndata is available (version: X.X.X)"，说明配置成功！**

# 第七步：完整功能测试

## 7.1 运行测试脚本

在RStudio中运行以下测试代码：

```r
# 测试脚本
cat("=== ICellbioRpy 完整功能测试 ===\n")

# 1. 加载包
library(ICellbioRpy)
cat("✓ ICellbioRpy包加载成功\n")

# 2. 配置Python环境
configure_python_env(conda_env = "1cellbio", verbose = TRUE)
cat("✓ Python环境配置成功\n")

# 3. 检查anndata
check_anndata_available()
cat("✓ anndata检查通过\n")

# 4. 检查主要函数
functions_to_check <- c("read1Cellbio", "iCellbio2H5ad", "seurat_to_h5ad", "read_10x_mtx_to_h5ad")
for (func in functions_to_check) {
  if (exists(func)) {
    cat("✓", func, "函数可用\n")
  } else {
    cat("✗", func, "函数不可用\n")
  }
}

cat("\n=== 测试完成 ===\n")
cat("如果所有项目都显示 ✓，说明环境配置完全成功！\n")
```

## 7.2 创建简单示例

测试一个简单的Python-R交互：

```r
# 简单的Python-R交互测试
library(reticulate)

# 在Python中创建数据
py_run_string("
import numpy as np
import pandas as pd
test_data = np.random.rand(5, 3)
df = pd.DataFrame(test_data, columns=['A', 'B', 'C'])
")

# 在R中访问Python数据
py_df <- py$df
print(py_df)

cat("✓ Python-R数据交互测试成功！\n")
```

# 🔧 常见问题排查

## 问题1：R命令不被识别

**症状**：在命令提示符中输入`R --version`时提示"R不是内部或外部命令"

**解决方法**：
1. 按 `Win + X`，选择"系统"
2. 点击"高级系统设置"
3. 点击"环境变量"
4. 在"系统变量"中找到"Path"，点击"编辑"
5. 点击"新建"，添加R的安装路径（通常是`C:\Program Files\R\R-4.x.x\bin`）
6. 点击"确定"保存
7. 重新打开命令提示符测试

## 问题2：Anaconda安装后conda命令无效

**症状**：在命令提示符中输入`conda --version`提示命令不存在

**解决方法**：
1. 重新运行Anaconda安装程序
2. 确保勾选"Add Anaconda3 to my PATH environment variable"
3. 或者手动添加conda路径到环境变量：
   - 添加：`C:\Users\YourName\anaconda3\Scripts`
   - 添加：`C:\Users\YourName\anaconda3\Library\bin`

## 问题3：reticulate找不到conda环境

**症状**：`conda_list()`返回空结果或找不到icellbio环境

**解决方法**：
```r
# 手动指定conda路径
Sys.setenv(RETICULATE_CONDA = "C:/Users/YourName/anaconda3/Scripts/conda.exe")

# 重新检查
conda_list()

# 或者直接指定Python路径
use_python("C:/Users/YourName/anaconda3/envs/1cellbio/python.exe")
```

## 问题4：anndata导入失败

**症状**：Python可以导入anndata但R中报错

**解决方法**：
```r
# 强制重新配置Python环境
library(reticulate)
use_condaenv("1cellbio", required = TRUE)

# 手动测试Python导入
py_run_string("import anndata; print(anndata.__version__)")

# 如果还是失败，重新安装anndata
system("conda activate 1cellbio && conda install -c conda-forge anndata")
```

## 问题5：网络连接问题

**症状**：下载速度慢或连接失败

**解决方法**：
1. **使用国内镜像**：
   ```bash
   # 配置conda使用清华镜像
   conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
   conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
   conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
   ```

2. **在R中设置镜像**：
   ```r
   # 设置CRAN镜像
   options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
   ```

# 💡 成功标志

当您完成所有步骤后，应该能够：

1. ✅ **在RStudio中运行**：`library(ICellbioRpy)`
2. ✅ **成功配置Python环境**：`configure_python_env(conda_env = "1cellbio")`
3. ✅ **验证anndata可用**：`check_anndata_available()`
4. ✅ **查看所有函数**：`ls("package:ICellbioRpy")`

# 🎉 恭喜！

如果您完成了所有步骤，恭喜您已经成功搭建了完整的ICellbioRpy使用环境！

## 下一步

现在您可以：
1. 学习[教程1：1CellBio转H5AD](tutorial_1_1cellbio_to_h5ad.md)
2. 学习[教程2：1CellBio转Seurat](tutorial_2_1cellbio_to_seurat.md)
3. 根据您的数据类型选择合适的教程

## 保存您的环境配置

建议在RStudio中创建一个启动脚本：

```r
# 保存为 setup_environment.R
library(ICellbioRpy)
configure_python_env(conda_env = "1cellbio", verbose = FALSE)

cat("✓ ICellbioRpy环境已就绪！\n")
cat("Python环境：1cellbio\n")
cat("可用函数：", paste(ls("package:ICellbioRpy"), collapse = ", "), "\n")
```

每次使用ICellbioRpy时，只需运行：`source("setup_environment.R")`

---

**需要帮助？**
- 查看其他教程获取使用指导
- 访问GitHub项目页面寻求支持
- 加入用户交流群获取帮助

**祝您单细胞数据分析顺利！** 🧬✨