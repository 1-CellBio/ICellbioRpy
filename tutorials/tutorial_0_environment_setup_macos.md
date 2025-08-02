# 🎯 教程目标

这是一个专门为**完全零基础的初学者**准备的macOS系统环境安装教程。我们将从最基础的软件下载开始，一步一步地搭建完整的ICellbioRpy使用环境。

**完成本教程后，您将拥有：**
- ✅ R和RStudio工作环境
- ✅ Anaconda Python环境
- ✅ 专用的conda环境
- ✅ ICellbioRpy包及其依赖
- ✅ 完全配置好的Python-R交互环境

**预计时间：** 1-2小时（取决于网络速度）

# 📋 开始前的准备

## 检查您的系统

1. **确认您使用的是macOS系统**（本教程适用于macOS 10.15+）
2. **检查您的Mac芯片类型**：
   - 点击屏幕左上角的🍎图标
   - 选择"关于本机"
   - 查看"芯片"或"处理器"信息：
     - **Apple Silicon (M1/M2/M3)**：显示"Apple M1"、"Apple M2"等
     - **Intel**：显示"Intel Core i5"、"Intel Core i7"等
3. **确保您有管理员权限**（安装软件时可能需要）
4. **确保网络连接稳定**（需要下载较大的安装文件）
5. **准备至少5GB的磁盘空间**

## 本教程的特点

- 🔍 **超详细步骤**：每一步都有说明
- 🍎 **Apple Silicon优化**：针对M系列芯片的特殊说明
- 🚨 **错误预防**：提前说明可能遇到的问题
- 💡 **小贴士**：实用的操作建议
- 🔧 **故障排除**：常见问题的解决方法

# 第一步：安装Xcode Command Line Tools

在macOS上，很多软件编译需要Xcode Command Line Tools。

## 1.1 安装Command Line Tools

1. **打开终端（Terminal）**：
   - 按 `Cmd + 空格键` 打开Spotlight搜索
   - 输入"Terminal"或"终端"
   - 按回车键打开

2. **安装Command Line Tools**：
   ```bash
   xcode-select --install
   ```

3. **在弹出的对话框中**：
   - 点击 **"Install"（安装）**
   - 同意许可协议
   - 等待下载和安装完成（可能需要10-30分钟）

4. **验证安装**：
   ```bash
   xcode-select -p
   ```
   应该看到类似：`/Applications/Xcode.app/Contents/Developer` 或 `/Library/Developer/CommandLineTools`

✅ **看到路径信息说明Command Line Tools安装成功！**

# 第二步：安装R语言

## 2.1 下载R

1. **打开Safari或Chrome浏览器**

2. **访问R官方网站**：
   - 在地址栏输入：`https://www.r-project.org/`
   - 按回车键

3. **进入下载页面**：
   - 点击页面上的 **"CRAN"** 链接
   - 选择一个镜像站点（建议选择日本或香港的镜像，速度较快）

4. **选择macOS版本**：
   - 点击 **"Download R for macOS"**
   - **重要**：根据您的Mac芯片选择合适版本：
     - **Apple Silicon (M1/M2/M3)**：选择 **"R-4.x.x-arm64.pkg"**
     - **Intel Mac**：选择 **"R-4.x.x.pkg"**

💡 **小贴士**：如果不确定芯片类型，选择Intel版本也可以在Apple Silicon上运行（通过Rosetta）。

## 2.2 安装R

1. **找到下载的文件**：
   - 通常在"下载"文件夹中
   - 文件名类似：`R-4.4.3-arm64.pkg` 或 `R-4.4.3.pkg`

2. **双击安装文件**

3. **按照安装向导操作**：
   - 点击 **"继续"**
   - 阅读许可协议，点击 **"同意"**
   - 选择安装位置：**保持默认**（通常是Macintosh HD）
   - 点击 **"安装"**
   - 输入您的密码（管理员密码）
   - 等待安装完成

4. **完成安装**：
   - 看到"安装成功"提示
   - 点击 **"关闭"**

## 2.3 验证R安装

1. **打开终端**

2. **输入以下命令**：
   ```bash
   R --version
   ```

3. **如果看到类似输出，说明安装成功**：
   ```
   R version 4.4.3 (2025-02-28) -- "Trophy Case"
   Copyright (C) 2025 The R Foundation for Statistical Computing
   Platform: aarch64-apple-darwin20 (64-bit)
   ```

✅ **看到版本信息说明R安装成功！**

# 第三步：安装RStudio

## 3.1 下载RStudio

1. **访问RStudio官网**：
   - 在浏览器中输入：`https://rstudio.com/products/rstudio/download/`

2. **选择免费版本**：
   - 找到 **"RStudio Desktop"**
   - 点击 **"Download"** 按钮

3. **下载macOS版本**：
   - 系统会自动检测您的操作系统
   - 点击 **"Download RStudio Desktop for macOS"**
   - 下载文件名类似：`RStudio-2024.xx.x-xxx.dmg`

## 3.2 安装RStudio

1. **打开下载的DMG文件**：
   - 双击下载的`.dmg`文件
   - 系统会挂载磁盘镜像

2. **拖动安装**：
   - 在打开的窗口中，将RStudio图标拖动到Applications文件夹
   - 等待复制完成

3. **启动RStudio**：
   - 打开Finder，进入"应用程序"文件夹
   - 找到RStudio，双击启动
   - 如果出现安全提示，点击 **"打开"**

## 3.3 首次启动RStudio

1. **RStudio界面说明**：
   - **左下角**：Console（控制台）- 在这里运行R代码
   - **左上角**：Source（源代码编辑器）
   - **右上角**：Environment（环境变量）
   - **右下角**：Files/Plots/Packages（文件/图形/包管理）

2. **测试RStudio**：
   - 在Console中输入：`1 + 1`
   - 按回车，应该看到结果：`[1] 2`

✅ **如果看到正确结果，说明R和RStudio都安装成功了！**

# 第四步：安装Anaconda

## 4.1 为什么需要Anaconda？

Anaconda是一个Python发行版，它：
- 包含Python解释器
- 包含科学计算常用的库
- 提供conda包管理器
- 可以创建独立的虚拟环境
- 在macOS上兼容性更好

## 4.2 下载Anaconda

1. **访问Anaconda官网**：
   - 在浏览器中输入：`https://www.anaconda.com/products/distribution`

2. **选择macOS版本**：
   - 点击 **"Download"** 按钮
   - **重要**：根据您的Mac芯片选择：
     - **Apple Silicon (M1/M2/M3)**：选择 **"64-Bit (M1) Graphical Installer"**
     - **Intel Mac**：选择 **"64-Bit Graphical Installer"**

💡 **小贴士**：Anaconda安装文件很大（约500MB），请确保网络稳定。

## 4.3 安装Anaconda

1. **运行安装文件**：
   - 找到下载的文件（类似`Anaconda3-2024.xx-MacOSX-arm64.pkg`）
   - 双击运行

2. **安装步骤**：
   - 点击 **"继续"**
   - 点击 **"同意"** 接受许可协议
   - 选择安装位置：**保持默认**（推荐）
   - 点击 **"安装"**
   - 输入管理员密码
   - 等待安装完成（可能需要10-20分钟）

3. **完成安装**：
   - 看到安装成功提示
   - 点击 **"关闭"**

## 4.4 配置Shell环境

Anaconda安装后需要配置Shell环境：

1. **打开终端**

2. **初始化conda**：
   ```bash
   /opt/anaconda3/bin/conda init zsh
   ```
   
   如果您使用bash shell，运行：
   ```bash
   /opt/anaconda3/bin/conda init bash
   ```

3. **重新启动终端**：
   - 关闭当前终端窗口
   - 重新打开终端

4. **验证conda安装**：
   ```bash
   conda --version
   ```
   应该看到类似：`conda 24.x.x`

5. **验证Python**：
   ```bash
   python --version
   ```
   应该看到类似：`Python 3.11.x`

✅ **如果看到版本信息，说明Anaconda安装成功！**

# 第五步：创建专用conda环境

## 5.1 为什么需要专用环境？

创建专用环境的好处：
- 避免包冲突
- 更好的版本管理
- 出问题时可以重新创建
- 不影响系统默认Python

## 5.2 创建ICellbioRpy专用环境

1. **确保终端已经配置好conda**（应该看到命令提示符前有`(base)`）

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
   - 命令提示符应该变成：`(1cellbio) username@MacBook-Pro ~ %`
   - 这表示您现在在icellbio环境中

## 5.3 在新环境中安装anndata

1. **确保在icellbio环境中**（命令提示符前有`(1cellbio)`）

2. **安装anndata和相关包**：
   ```bash
   conda install -c conda-forge anndata pandas numpy scipy matplotlib
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

# 第六步：安装R包和依赖

## 6.1 在RStudio中安装基础R包

1. **打开RStudio**

2. **安装基础包**：
   在Console中逐行运行以下代码：
   ```r
   # 安装基础包
   install.packages(c("devtools", "reticulate", "Matrix", "data.table", "hdf5r", "jsonlite"))
   ```

3. **处理编译问题**（如果遇到）：
   如果遇到编译错误，可以尝试安装二进制版本：
   ```r
   # 强制安装二进制版本
   install.packages(c("devtools", "reticulate", "Matrix", "data.table"), 
                    type = "binary")
   ```

4. **安装Bioconductor相关包**：
   ```r
   # 安装BiocManager
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   
   # 安装生物信息学相关包
   BiocManager::install(c("SingleCellExperiment", "SummarizedExperiment", "S4Vectors"))
   ```

## 6.2 处理特殊的macOS问题

### 6.2.1 解决hdf5r安装问题

hdf5r包在macOS上可能需要特殊处理：

1. **首先尝试brew安装hdf5**（在终端中）：
   ```bash
   # 如果没有brew，先安装Homebrew
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
   
   # 安装hdf5
   brew install hdf5
   ```

2. **然后在RStudio中安装hdf5r**：
   ```r
   # 设置环境变量
   Sys.setenv(HDF5_DIR = "/opt/homebrew/opt/hdf5")  # Apple Silicon
   # 或者对于Intel Mac:
   # Sys.setenv(HDF5_DIR = "/usr/local/opt/hdf5")
   
   # 安装hdf5r
   install.packages("hdf5r")
   ```

### 6.2.2 Apple Silicon特殊处理

如果您使用Apple Silicon Mac，可能需要：

```r
# 设置编译器路径
Sys.setenv(CC = "/opt/homebrew/bin/gcc-13")  # 如果安装了gcc
Sys.setenv(CXX = "/opt/homebrew/bin/g++-13")

# 或者使用系统默认编译器
Sys.setenv(CC = "clang")
Sys.setenv(CXX = "clang++")
```

## 6.3 安装ICellbioRpy包

**注意**：这里我们假设ICellbioRpy包已经发布到GitHub。

```r
# 从GitHub安装ICellbioRpy包
devtools::install_github("your_username/ICellbioRpy")

# 如果遇到网络问题，可以设置代理或使用镜像
# options(download.file.method = "libcurl")
```

如果从本地安装：
```r
# 从本地安装（如果您有源码包）
devtools::install_local("path/to/ICellbioRpy")
```

# 第七步：配置Python-R交互环境

## 7.1 在RStudio中配置reticulate

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

## 7.2 处理macOS特殊的路径问题

在macOS上，conda环境路径可能比较复杂：

```r
# 如果自动检测有问题，手动指定路径
# 检查conda环境位置
system("conda info --envs")

# Apple Silicon典型路径
# /opt/anaconda3/envs/1cellbio/bin/python

# Intel Mac典型路径  
# /Users/YourName/anaconda3/envs/1cellbio/bin/python

# 手动指定（根据实际路径调整）
use_python("/opt/anaconda3/envs/1cellbio/bin/python")
```

## 7.3 测试ICellbioRpy的Python环境配置

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

# 第八步：完整功能测试

## 8.1 运行测试脚本

在RStudio中运行以下测试代码：

```r
# 测试脚本
cat("=== ICellbioRpy macOS 完整功能测试 ===\n")

# 1. 检查系统信息
cat("系统信息:", Sys.info()["sysname"], Sys.info()["machine"], "\n")

# 2. 加载包
library(ICellbioRpy)
cat("✓ ICellbioRpy包加载成功\n")

# 3. 配置Python环境
configure_python_env(conda_env = "1cellbio", verbose = TRUE)
cat("✓ Python环境配置成功\n")

# 4. 检查anndata
check_anndata_available()
cat("✓ anndata检查通过\n")

# 5. 检查主要函数
functions_to_check <- c("read1Cellbio", "iCellbio2H5ad", "seurat_to_h5ad", "read_10x_mtx_to_h5ad")
for (func in functions_to_check) {
  if (exists(func)) {
    cat("✓", func, "函数可用\n")
  } else {
    cat("✗", func, "函数不可用\n")
  }
}

# 6. 检查编译工具
cat("\nXcode Command Line Tools:", 
    ifelse(system("xcode-select -p", ignore.stdout = TRUE) == 0, "✓ 已安装", "✗ 未安装"), "\n")

cat("\n=== 测试完成 ===\n")
cat("如果所有项目都显示 ✓，说明环境配置完全成功！\n")
```

## 8.2 创建简单示例

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
print('Python DataFrame created successfully!')
")

# 在R中访问Python数据
py_df <- py$df
print(py_df)

cat("✓ Python-R数据交互测试成功！\n")
```

# 🔧 常见问题排查

## 问题1：Command Line Tools安装失败

**症状**：`xcode-select --install`命令没有弹出安装对话框

**解决方法**：
```bash
# 方法1：直接从App Store安装Xcode
# 打开App Store，搜索"Xcode"，点击安装

# 方法2：重置xcode-select
sudo xcode-select --reset
xcode-select --install

# 方法3：手动下载
# 访问 https://developer.apple.com/download/more/
# 下载适合您系统版本的Command Line Tools
```

## 问题2：conda命令无效

**症状**：在终端中输入`conda --version`提示命令不存在

**解决方法**：
```bash
# 1. 检查Anaconda是否安装在默认位置
ls /opt/anaconda3/bin/conda

# 2. 手动添加到PATH
echo 'export PATH="/opt/anaconda3/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc

# 3. 重新初始化conda
/opt/anaconda3/bin/conda init zsh

# 4. 重启终端
```

## 问题3：Apple Silicon兼容性问题

**症状**：某些包在Apple Silicon上安装失败

**解决方法**：
```bash
# 1. 使用Rosetta运行终端
# 右键点击终端应用 → 显示简介 → 勾选"使用Rosetta打开"

# 2. 安装x86版本的conda
arch -x86_64 conda install package_name

# 3. 在R中设置特定的repository
options(repos = c(CRAN = "https://cloud.r-project.org"))
```

## 问题4：reticulate找不到Python环境

**症状**：`conda_list()`返回空结果或报错

**解决方法**：
```r
# 1. 检查conda路径
Sys.which("conda")

# 2. 手动设置conda路径
Sys.setenv(RETICULATE_CONDA = "/opt/anaconda3/bin/conda")

# 3. 直接指定Python路径
use_python("/opt/anaconda3/envs/1cellbio/bin/python")

# 4. 重启R会话
.rs.restartR()
```

## 问题5：hdf5r安装失败

**症状**：hdf5r包编译失败

**解决方法**：
```bash
# 在终端中安装HDF5库
brew install hdf5

# 如果是Apple Silicon
brew install --cask homebrew/cask-versions/adoptopenjdk8
```

然后在R中：
```r
# 设置环境变量
Sys.setenv(HDF5_DIR = "/opt/homebrew/opt/hdf5")  # Apple Silicon
# 或者
Sys.setenv(HDF5_DIR = "/usr/local/opt/hdf5")     # Intel

# 重新安装
install.packages("hdf5r", type = "source")
```

## 问题6：权限问题

**症状**：安装时提示权限不足

**解决方法**：
```bash
# 修复用户目录权限
sudo chown -R $(whoami):staff ~/anaconda3

# 修复R包库权限
sudo chown -R $(whoami):staff /Library/Frameworks/R.framework/Versions/*/Resources/library
```

# 💡 Apple Silicon (M1/M2/M3) 特殊注意事项

## 优化建议

1. **使用原生ARM64版本**：
   - R选择arm64版本
   - Anaconda选择M1版本
   - 避免通过Rosetta运行（除非必要）

2. **编译器设置**：
   ```r
   # 在R中设置编译器
   Sys.setenv(CC = "clang")
   Sys.setenv(CXX = "clang++")
   ```

3. **包安装策略**：
   ```r
   # 优先使用二进制包
   options(pkgType = "binary")
   
   # 如果必须编译，设置编译参数
   options(configure.args = "--disable-openmp")
   ```

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
cat("系统:", Sys.info()["sysname"], Sys.info()["machine"], "\n")
cat("Python环境：icellbio\n")
cat("可用函数：", paste(ls("package:ICellbioRpy"), collapse = ", "), "\n")
```

每次使用ICellbioRpy时，只需运行：`source("setup_environment.R")`

## macOS特殊的性能优化

```r
# 创建 ~/.Rprofile 文件进行全局配置
options(Ncpus = parallel::detectCores())  # 使用所有CPU核心
options(repos = c(CRAN = "https://cloud.r-project.org"))  # 设置快速镜像

# 为Apple Silicon优化
if (Sys.info()["machine"] == "arm64") {
  options(pkgType = "binary")  # 优先使用二进制包
}
```

---

**需要帮助？**
- 查看其他教程获取使用指导
- 访问GitHub项目页面寻求支持
- 加入用户交流群获取帮助

**祝您单细胞数据分析顺利！** 🧬✨