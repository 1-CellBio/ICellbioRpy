# 🎯 教程目标

在这个教程中，我们将学习如何使用`ICellbioRpy`包将1CellBio分析结果的ZIP文件直接转换为H5AD格式，方便在Python的scanpy包中进行下游分析。

# 📋 前提条件

## 必需的软件环境

1. **R (≥ 4.0.0)**
2. **Python (≥ 3.7)** 
3. **anndata Python包** - 用于处理H5AD文件

## 检查Python环境

在开始之前，请确保您的系统中安装了Python和anndata包：

```bash
# 在终端中检查Python版本
python --version

# 检查是否安装了anndata
python -c "import anndata; print('anndata version:', anndata.__version__)"
```

如果没有安装anndata，请运行：

```bash
# 使用pip安装
pip install anndata

# 或者使用conda安装
conda install -c conda-forge anndata
```

# 🚀 开始教程

## 第1步: 安装和加载包

```r
# 如果还没有安装ICellbioRpy包，请先安装
# install.packages("devtools")
# devtools::install_github("your_username/ICellbioRpy")

# 加载必需的包
library(ICellbioRpy)
library(reticulate)  # 用于Python交互
```

## 第2步: 配置Python环境

这是**最重要的一步**！必须正确配置Python环境才能使用H5AD相关功能。

```r
# 方法1: 使用当前默认Python环境（推荐新手使用）
configure_python_env(verbose = TRUE)

# 方法2: 如果您使用conda环境，指定环境名称
# configure_python_env(conda_env = "your_env_name", verbose = TRUE)

# 方法3: 如果您想使用特定的Python路径
# configure_python_env(python_path = "/usr/local/bin/python3", verbose = TRUE)

# 验证配置是否成功
check_anndata_available()
```

**💡 配置说明：**

- `verbose = TRUE`: 显示详细的配置信息，帮助排查问题
- 如果您看到 "✓ anndata is available" 消息，说明配置成功
- 如果配置失败，请检查Python环境和anndata安装

## 第3步: 准备数据文件

您需要有一个1CellBio分析结果的ZIP文件。这个文件通常包含以下内容：

- **元数据文件** (JSON格式)
- **表达矩阵** (HDF5格式)  
- **细胞注释** (HDF5格式)
- **基因注释** (HDF5格式)
- **降维结果** (如PCA, UMAP, t-SNE等)

```r
# 示例：设置您的ZIP文件路径
zip_file_path <- "path/to/your/1cellbio_results.zip"

# 检查文件是否存在
if (file.exists(zip_file_path)) {
  cat("✓ 找到ZIP文件:", zip_file_path, "\n")
  cat("文件大小:", round(file.size(zip_file_path) / 1024^2, 2), "MB\n")
} else {
  cat("✗ 找不到ZIP文件，请检查路径是否正确\n")
}
```

## 第4步: 执行转换

现在我们可以将ZIP文件转换为H5AD格式：

```r
# 设置输出文件路径
output_h5ad_path <- "1cellbio_data.h5ad"

# 执行转换
# 这个函数会：
# 1. 解压ZIP文件
# 2. 读取所有数据组件
# 3. 组装成标准的单细胞数据格式  
# 4. 保存为H5AD文件
iCellbio2H5ad(
  zip_path = zip_file_path,
  h5ad_path = output_h5ad_path,
  overwrite = FALSE,              # 新增：默认不覆盖已存在文件
  name_conflict = "make_unique",  # 新增：命名冲突策略（或设为 "error"）
  verbose = TRUE                  # 显示详细进度信息
)
```

**🔍 转换过程说明：**

转换过程中您会看到以下信息：

1. **解压文件**: 提取ZIP内容到临时目录
2. **读取元数据**: 解析实验配置和参数
3. **加载表达矩阵**: 读取基因表达数据（通常是稀疏矩阵）
4. **加载注释**: 读取细胞和基因的注释信息
5. **加载降维结果**: 读取PCA、UMAP等降维数据
6. **创建H5AD**: 组装所有数据并保存

## 第5步: 验证转换结果

```r
# 检查输出文件是否成功创建
if (file.exists(output_h5ad_path)) {
  cat("✓ H5AD文件创建成功!\n")
  cat("文件路径:", output_h5ad_path, "\n")
  cat("文件大小:", round(file.size(output_h5ad_path) / 1024^2, 2), "MB\n")
  
  # 可以尝试在R中加载验证（需要anndata包）
  tryCatch({
    ad <- reticulate::import("anndata")
    adata <- ad$read_h5ad(output_h5ad_path)
    cat("数据维度:", adata$shape[1], "细胞 ×", adata$shape[2], "基因\n")
    cat("可用的降维结果:", names(adata$obsm), "\n")
  }, error = function(e) {
    cat("注意: 无法在R中预览H5AD文件，但文件已成功创建\n")
  })
} else {
  cat("✗ H5AD文件创建失败，请检查错误信息\n")
}
```

# 🐍 在Python中使用转换后的数据

转换完成后，您可以在Python中使用scanpy来分析数据：

```python
import scanpy as sc
import pandas as pd

# 读取转换后的H5AD文件
adata = sc.read_h5ad('1cellbio_data.h5ad')

# 查看数据基本信息
print(f"数据维度: {adata.shape[0]} 细胞 × {adata.shape[1]} 基因")
print(f"细胞注释列: {list(adata.obs.columns)}")
print(f"基因注释列: {list(adata.var.columns)}")
print(f"降维结果: {list(adata.obsm.keys())}")

# 可视化UMAP（如果有的话）
if 'X_umap' in adata.obsm:
    sc.pl.umap(adata, color='cell_type')  # 根据实际的注释列名调整
```

# 🔧 常见问题排查

## 问题1: Python环境配置失败

**错误信息**: "Failed to configure Python environment"

**解决方法**:
```r
# 1. 检查reticulate配置
reticulate::py_config()

# 2. 手动指定Python路径
# 首先在终端运行: which python 或 which python3
# 然后使用返回的路径
configure_python_env(python_path = "/usr/bin/python3")

# 3. 如果使用conda
# reticulate::conda_list()  # 查看可用环境
# configure_python_env(conda_env = "base")
```

## 问题2: anndata包找不到

**错误信息**: "Package 'anndata' is not available"

**解决方法**:
```bash
# 在终端中安装anndata
pip install anndata

# 或者使用conda
conda install -c conda-forge anndata
```

## 问题3: ZIP文件损坏或格式不正确

**错误信息**: "Error reading ZIP file" 或 "Invalid file structure"

**解决方法**:
- 确保ZIP文件是完整的1CellBio输出
- 检查文件是否损坏（重新下载）
- 确认ZIP文件包含必需的组件

## 问题4: 内存不足

**错误信息**: "Memory allocation failed"

**解决方法**:
```r
# 1. 增加R的内存限制（仅限Windows）
# memory.limit(size = 8000)  # 设置为8GB

# 2. 清理环境
gc()  # 垃圾回收

# 3. 分批处理大文件（高级用法）
```

# 📚 下一步

恭喜！您已经成功将1CellBio数据转换为H5AD格式。现在您可以：

1. **在Python中分析**: 使用scanpy进行质控、聚类、差异分析等
2. **转换为其他格式**: 继续学习其他教程
3. **整合多个数据集**: 使用scanpy的整合功能

建议继续学习：
- 教程2: 1CellBio转Seurat格式
- 教程3: Seurat转H5AD格式  
- 教程4: 10X MTX数据读取和整合

# 💡 小贴士

1. **备份原始数据**: 转换前请备份原始ZIP文件
2. **检查磁盘空间**: H5AD文件可能比ZIP文件大
3. **版本兼容**: 确保anndata版本 ≥ 0.7.0
4. **性能优化**: 大数据集转换可能需要较长时间，请耐心等待

---

**需要帮助？** 
- 查看包文档: `?iCellbio2H5ad`
- 检查GitHub Issues: [项目地址]
- 联系维护者: [邮箱地址]