# ICellbioRpy 📊🧬

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **一站式单细胞数据格式转换解决方案**  
> 支持1CellBio、Seurat、SingleCellExperiment、10X MTX与H5AD格式间的灵活转换

## 🎯 概述

ICellbioRpy是一个专业的R包，用于单细胞RNA测序数据格式转换。它提供了完整的工具链，支持从1CellBio分析结果到主流单细胞分析格式的无缝转换，特别针对初学者优化了使用体验。

### ✨ 主要特性

- 🚀 **一键转换**：支持1CellBio ZIP → H5AD/Seurat/SingleCellExperiment
- 🗂️ **Stereo-seq支持**：原生读取GEF文件，完整保留细胞边界信息
- 📊 **多样本整合**：从多个10X MTX文件整合为单一H5AD文件
- 🎨 **空间可视化**：带细胞边界的高质量空间转录组可视化
- 💾 **内存高效**：自动保留稀疏矩阵格式，节省70-90%存储空间
- 🔄 **双向转换**：支持R对象与Python格式间的相互转换
- 🎨 **初学者友好**：提供超详细的环境安装教程
- 🐍 **智能Python配置**：自动检测和配置Python环境

---

## 🚀 快速开始

### 📋 环境配置（初学者必读）

**如果您是完全初学者，请先完成环境安装：**

| 操作系统 | 安装教程 | 说明 |
|---------|---------|------|
| **Windows** | [📘 Windows环境安装指南](tutorials/tutorial_0_environment_setup_windows.md) | 从R安装到conda环境的完整指南 |
| **macOS** | [📘 macOS环境安装指南](tutorials/tutorial_0_environment_setup_macos.md) | 支持Apple Silicon优化 |

### ⚡ 环境检测

```r
# 快速检测环境是否配置完成
source("tutorials/check_environment.R")
```

### 📦 安装

```r
# 安装依赖包
install.packages(c("devtools", "reticulate", "Matrix", "data.table"))

# 从GitHub安装ICellbioRpy
devtools::install_github("1-Cellbio/ICellbioRpy")
```

---

## 📚 完整教程

我们提供了面向初学者的详细教程，每个教程都可以在RStudio中直接运行：

### 🔧 环境安装教程
- [📘 Windows环境安装](tutorials/tutorial_0_environment_setup_windows.md) - Windows系统完整配置指南
- [📘 macOS环境安装](tutorials/tutorial_0_environment_setup_macos.md) - macOS系统完整配置指南

### 📊 数据转换教程

| 教程 | 输入格式 | 输出格式 | 难度 | 用时 |
|------|---------|---------|------|------|
| [教程1](tutorials/tutorial_1_1cellbio_to_h5ad.md) | 1CellBio ZIP | H5AD | ⭐⭐ | 30-45min |
| [教程2](tutorials/tutorial_2_1cellbio_to_seurat.md) | 1CellBio ZIP | Seurat | ⭐⭐ | 45-60min |
| [教程3](tutorials/tutorial_3_seurat_to_h5ad.md) | Seurat对象 | H5AD | ⭐⭐⭐ | 40-50min |
| [教程4](tutorials/tutorial_4_10x_mtx_to_h5ad.md) | 10X MTX文件 | H5AD | ⭐⭐⭐⭐ | 50-70min |

### 📖 其他资源
- [📍 新手导航](tutorials/START_HERE.md) - 从这里开始您的学习之旅
- [⚡ 快速指南](tutorials/QUICK_GUIDE.md) - 代码片段和快速参考
- [📚 教程总览](tutorials/README.md) - 完整教程列表和说明

---

## 💻 核心功能

### 🔄 Python环境配置

**重要**：所有H5AD相关功能都需要Python环境支持，ICellbioRpy提供智能配置：

```r
library(ICellbioRpy)

# 方法1：自动检测并使用当前环境（推荐）
configure_python_env(verbose = TRUE)

# 方法2：指定conda环境
configure_python_env(conda_env = "1cellbio", verbose = TRUE)

# 方法3：指定Python路径
configure_python_env(python_path = "/path/to/python", verbose = TRUE)

# 验证配置
check_anndata_available()
```

### 📈 主要使用场景

#### 1. 1CellBio → H5AD（用于Python分析）

```r
# 直接转换（推荐）
iCellbio2H5ad("path/to/1cellbio_results.zip", "output.h5ad")

# 或两步转换
data <- read1Cellbio("path/to/1cellbio_results.zip")
as.h5ad(data, "output.h5ad")
```

#### 2. 1CellBio → Seurat（用于R分析）

```r
# 读取数据
data <- read1Cellbio("path/to/1cellbio_results.zip")

# 转换为Seurat对象（需要指定基因名和细胞名列）
# 函数会自动显示所有可用的列名供参考
seurat_obj <- as.Seurat.1CB(data, 
                           rownames = "id",        # 基因名列
                           colnames = "cell_id")   # 细胞名列

# 可视化
library(Seurat)
DimPlot(seurat_obj, reduction = "umap", group.by = "level1class")
```

**💡 重要提示：**
- `rownames` 和 `colnames` 参数是必填的
- 如果不确定列名，调用函数时会自动显示所有可用选项
- 如遇到重复名称，默认会自动添加后缀（如 Gene-1, Gene-2）。若希望在名称冲突时直接报错，可在相关函数中设置 `name_conflict = "error"`。

#### 3. Seurat → H5AD（跨语言转换）

```r
# 将R中的Seurat对象转换为Python可用的H5AD格式
seurat_to_h5ad(seurat_obj, "output.h5ad")

# 高级选项
seurat_to_h5ad(seurat_obj, "output.h5ad",
               default_assay = "RNA",
               layer = "data",
               include_reductions = TRUE,
               overwrite = FALSE,                  # 若目标存在，默认拒绝覆盖（可设TRUE允许）
               name_conflict = "make_unique")     # 基因/细胞命名冲突策略: make_unique|error
```

#### 4. 10X MTX → H5AD（多样本整合）

```r
# 准备样本信息CSV文件
sample_info <- data.frame(
  Sample_id = c("sample1", "sample2"),
  mtx_fns = c("path/to/sample1/matrix.mtx.gz", "path/to/sample2/matrix.mtx.gz"),
  features_fns = c("path/to/sample1/features.tsv.gz", "path/to/sample2/features.tsv.gz"),
  barcodes_fns = c("path/to/sample1/barcodes.tsv.gz", "path/to/sample2/barcodes.tsv.gz")
)
write.csv(sample_info, "samples.csv", row.names = FALSE)

# 读取和整合
read_10x_mtx_to_h5ad(
  csv_file = "samples.csv",
  output_h5ad = "integrated.h5ad",
  min_counts_per_cell = 200
)
```

#### 5. Stereo-seq GEF → R对象（空间转录组）

```r
# 读取Stereo-seq GEF文件（带细胞边界）
stereo_data <- read_gef(
  file_path = "sample.cellbin.gef",
  bin_type = "cell_bins",
  include_cellborder = TRUE
)

# 转换为Seurat对象（自动保留细胞边界信息）
seurat_obj <- as.Seurat(stereo_data)

# 转换为SingleCellExperiment对象
sce_obj <- as.SingleCellExperiment(stereo_data)

# 空间可视化（带细胞边界）
plot_cells_with_borders(
  seurat_obj,
  color_by = "area",           # 按细胞面积着色
  show_borders = TRUE,         # 显示细胞边界
  border_color = "white"
)

# 内存优化读取（大文件处理）
stereo_subset <- read_gef(
  "large.cellbin.gef",
  max_cells = 5000,                           # 限制细胞数
  region = c(1000, 5000, 1000, 5000)        # 空间区域过滤
)
```

#### 6. Stereo-seq GEF → H5AD（Python分析）

```r
# 方法1：两步转换
stereo_data <- read_gef("sample.cellbin.gef", include_cellborder = TRUE)
stereo_to_h5ad(stereo_data, "output.h5ad", include_spatial = TRUE)

# 方法2：直接转换（推荐，内存高效）
gef_to_h5ad(
  gef_file = "sample.cellbin.gef",
  h5ad_file = "output.h5ad",
  bin_type = "cell_bins",
  include_cellborder = TRUE,    # cell_borders存储在adata.uns中
  include_spatial = TRUE,       # 空间坐标存储在adata.obsm中
  max_cells = 10000            # 内存管理
)

# 大文件处理：区域过滤 + 基因筛选
gef_to_h5ad(
  gef_file = "large_sample.cellbin.gef",
  h5ad_file = "filtered.h5ad",
  region = c(1000, 5000, 1000, 5000),      # 空间区域过滤
  gene_list = c("GAPDH", "ACTB", "CD45"),  # 基因过滤
  overwrite = TRUE
)
```

#### 7. H5AD → R对象（反向转换）

```r
# H5AD转SingleCellExperiment
sce <- h5ad_to_sce("data.h5ad")

# H5AD转Seurat
seurat_obj <- h5ad_to_seurat("data.h5ad")

# 指定数据层
sce <- h5ad_to_sce("data.h5ad", use_x_as = "counts")
```

---

## 🛠️ 主要函数

### 核心转换函数
- `read1Cellbio()` - 从ZIP文件读取1CellBio结果
- `read_gef()` - **新功能**：读取Stereo-seq GEF文件（支持细胞边界）
- `gef_to_h5ad()` - **新功能**：GEF文件直接转H5AD（细胞边界存储在uns中）
- `stereo_to_h5ad()` - **新功能**：StereoData对象转H5AD格式
- `iCellbio2H5ad()` - **直接转换**：ZIP → H5AD（内存高效）
- `as.Seurat.1CB()` - 1CellbioData → Seurat对象
- `as.SingleCellExperiment.1CB()` - 1CellbioData → SingleCellExperiment对象
- `plot_cells_with_borders()` - **新功能**：空间转录组可视化（细胞边界）
- `seurat_to_h5ad()` - Seurat对象 → H5AD文件
- `read_10x_mtx_to_h5ad()` - 多样本10X数据整合
- `h5ad_to_sce()` - H5AD → SingleCellExperiment
- `h5ad_to_seurat()` - H5AD → Seurat对象

### 环境配置函数
- `configure_python_env()` - 智能Python环境配置
- `check_anndata_available()` - 检查anndata可用性

---

## 🔬 数据结构支持

### 输入格式
- **1CellBio ZIP文件**：包含HDF5格式的表达矩阵、元数据和降维结果
- **Seurat对象**：标准Seurat对象，支持多个assays和降维结果
- **10X MTX文件**：Cell Ranger输出的稀疏矩阵格式
- **H5AD文件**：Python scanpy格式的单细胞数据

### 输出格式
- **H5AD**：保留稀疏矩阵，支持Python/scanpy分析
- **Seurat**：包含counts、data和降维结果
- **SingleCellExperiment**：Bioconductor标准格式
- **整合数据**：多样本合并的统一格式

---

## 💡 性能优势

### 内存效率
- **稀疏矩阵保留**：自动维持稀疏格式，节省70-90%存储空间
- **直接转换**：避免创建中间对象，减少内存占用
- **流式处理**：支持大数据集的高效处理

### 用户体验
- **智能检测**：自动检测和配置Python环境
- **详细文档**：每个函数都有完整的帮助文档和示例
- **错误处理**：友好的错误提示和解决建议
- **进度显示**：长时间操作显示详细进度

---

## 🎯 典型工作流程

### 场景1：Python用户接收1CellBio数据

```r
# 直接转换为H5AD格式
iCellbio2H5ad("1cellbio_results.zip", "analysis_data.h5ad")
```

```python
# 在Python中继续分析
import scanpy as sc
adata = sc.read_h5ad("analysis_data.h5ad")
sc.pl.umap(adata, color='level1class')
```

### 场景2：R用户进行单细胞分析

```r
# 读取数据
data <- read1Cellbio("1cellbio_results.zip")

# 转换为Seurat对象
seurat_obj <- as.Seurat.1CB(data)

# Seurat分析流程
library(Seurat)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindClusters(seurat_obj)
DimPlot(seurat_obj, reduction = "umap")
```

### 场景3：多样本10X数据整合

```r
# 创建样本信息文件
# 然后整合多个样本
read_10x_mtx_to_h5ad("samples.csv", "integrated.h5ad")
```

```python
# Python中的下游分析
import scanpy as sc
adata = sc.read_h5ad("integrated.h5ad")

# 批次校正和整合分析
sc.pp.combat(adata, key='sample_id')
sc.tl.leiden(adata)
sc.pl.umap(adata, color=['sample_id', 'leiden'])
```

---

## 🆘 技术支持

### 自助解决
1. **运行环境检测**：`source("tutorials/check_environment.R")`
2. **查看详细教程**：根据您的数据类型选择对应教程
3. **查看函数帮助**：`?function_name`

### 获取帮助
- 📧 **提交Issue**：[GitHub Issues](https://github.com/1-Cellbio/ICellbioRpy/issues)
- 📚 **查看教程**：[完整教程集合](tutorials/)
- 💬 **用户社区**：加入用户交流群

---

## 📋 系统要求

### R环境
- R ≥ 4.0.0
- 推荐包：Seurat, SingleCellExperiment, Matrix, data.table, hdf5r, jsonlite

### Python环境  
- Python ≥ 3.7
- 必需包：anndata ≥ 0.7.0
- 推荐包：scanpy, pandas, numpy

### 硬件建议
- **内存**：8GB+（大数据集需要更多）
- **存储**：确保足够的磁盘空间
- **CPU**：多核CPU可提高处理速度

---

## 📄 许可证

本项目采用 [MIT许可证](LICENSE)

---

**开始您的单细胞数据转换之旅！** 🚀

[![开始使用](https://img.shields.io/badge/开始使用-教程导航-blue)](tutorials/START_HERE.md)
[![快速指南](https://img.shields.io/badge/快速指南-代码示例-green)](tutorials/QUICK_GUIDE.md)
[![环境安装](https://img.shields.io/badge/环境安装-Windows|macOS-orange)](tutorials/)