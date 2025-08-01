
# ICellbioRpy

用于直接从 zip 文件读取 1Cellbio 管道结果，并将其转换为 Seurat、SingleCellExperiment 对象或 h5ad 格式。

---

## 概述

此包提供函数，可直接读取压缩的 1Cellbio 单细胞 RNA 测序分析结果，并将其转换为常用的 Bioconductor/Seurat 对象或 h5ad 格式，以便进行下游分析。

---

## 安装

```r
# 如果您尚未安装 devtools，请先安装
install.packages("devtools")

# 从 GitHub 安装 ICellbioRpy
devtools::install_github("1-Cellbio/ICellbioRpy")
````

-----

## Python 环境配置

**重要提示**：所有需要 Python/anndata 的函数都会自动使用您当前的环境配置 Python 环境！

### 自动配置

所有与 h5ad 相关的函数现在都会自动使用您当前的 Python 环境：

```r
library(ICellbioRpy)

# 这些函数会自动使用当前环境配置 Python 环境
iCellbio2H5ad("data.zip", "output.h5ad")  # 自动配置 Python
sce <- h5ad_to_sce("data.h5ad")             # 自动配置 Python  
seurat_obj <- h5ad_to_seurat("data.h5ad")   # 自动配置 Python
seurat_to_h5ad(seurat_obj, "output.h5ad")   # 自动配置 Python
```

### 手动配置（需要时）

如果您的当前环境中没有 anndata，您可以指定其他环境：

```r
# 配置使用特定的 conda 环境
configure_python_env(conda_env = "your_env_name")

# 或者配置使用特定的 Python 路径  
configure_python_env(python_path = "/path/to/python")

# 验证 anndata 是否可用
check_anndata_available()
```

**注意**：如果找不到 anndata，函数将提供有用的错误消息和安装说明。

有关详细的故障排除，请参阅 `anndata_installation_guide.md`。

-----

## 用法

### 方法 1：两步转换（通过 1CellbioData 对象）

```r
library(ICellbioRpy)

# 从 zip 文件读取 1Cellbio 结果
data <- read1Cellbio("path/to/1cellbio_results.zip")

# 转换为 SingleCellExperiment
sce <- as.SingleCellExperiment(data)

# 转换为 Seurat 对象
seurat <- as.Seurat(data)

# 转换为 h5ad 格式，用于 Scanpy
as.h5ad(data, "output.h5ad")
```

### 方法 2：直接转换为 h5ad（内存高效 + 稀疏矩阵保留）

```r
library(ICellbioRpy)

# 直接从 zip 文件转换为 h5ad 格式
# 这会保留稀疏矩阵格式，并且内存效率高
# 通常比密集存储减少 70-90% 的文件大小
iCellbio2H5ad("path/to/1cellbio_results.zip", "output.h5ad")
```

### 方法 3：将 h5ad 文件转换为 R 对象（反向转换）

```r
library(ICellbioRpy)

# 将 h5ad 文件转换为 SingleCellExperiment
sce <- h5ad_to_sce("data.h5ad")

# 将 h5ad 文件转换为 Seurat 对象
seurat_obj <- h5ad_to_seurat("data.h5ad")

# 指定哪个矩阵用作主要表达数据
# 使用 X 矩阵作为 counts 而不是 logcounts
sce <- h5ad_to_sce("data.h5ad", use_x_as = "counts")
seurat_obj <- h5ad_to_seurat("data.h5ad", use_x_as = "counts")
```

### 方法 4：将 Seurat 对象转换为 h5ad 文件

```r
library(ICellbioRpy)

# 将 Seurat 对象转换为 h5ad 文件
seurat_to_h5ad(seurat_obj, "output.h5ad")

# 使用 counts 而不是 data 层
seurat_to_h5ad(seurat_obj, "output.h5ad", layer = "counts")

# 指定不同的默认 assay (例如, SCT)
seurat_to_h5ad(seurat_obj, "output.h5ad", default_assay = "SCT")

# 如果不需要，排除降维结果
seurat_to_h5ad(seurat_obj, "output.h5ad", include_reductions = FALSE)
```

-----

## 函数

### 核心转换函数

  - `read1Cellbio()` - 从 zip 文件读取 1Cellbio 结果并创建 1CellbioData 对象
  - `as.SingleCellExperiment()` - 将 1CellbioData 对象转换为 SingleCellExperiment 对象
  - `as.Seurat()` - 将 1CellbioData 对象转换为 Seurat 对象
  - `iCellbio2H5ad()` - **新增**：直接将 zip 文件转换为 h5ad 格式，无需创建中间对象。**保留稀疏矩阵格式**，以实现最佳内存效率和存储（通常减少 70-90% 的文件大小）
  - `h5ad_to_sce()` - **新增**：将 h5ad 文件转换为 SingleCellExperiment 对象，并保留稀疏矩阵
  - `h5ad_to_seurat()` - **新增**：将 h5ad 文件转换为 Seurat 对象，并保留稀疏矩阵
  - `seurat_to_h5ad()` - **新增**：将 Seurat 对象转换为 h5ad 文件，并保留稀疏矩阵和元数据

### Python 环境配置

  - `configure_python_env()` - **新增**：配置 Python 环境并防止自动安装 anndata
  - `check_anndata_available()` - **新增**：检查当前 Python 环境中 anndata 是否可用

-----

## 详情

该包适用于标准的 1Cellbio 结果结构，包括：

  - 基因表达矩阵（counts 和 logcounts）
  - 细胞元数据
  - 基因元数据
  - 降维结果（PCA、tSNE、UMAP）

### 数据结构

1Cellbio 结果包含：

1.  **计数数据**：以 HDF5 稀疏矩阵形式存储的原始计数矩阵
2.  **Logcounts 数据**：以 HDF5 密集矩阵形式存储的对数转换归一化数据
3.  **细胞元数据**：包括组织、聚类分配、质量控制指标等每个细胞的信息
4.  **基因元数据**：每个基因的信息
5.  **降维结果**：PCA、tSNE 和 UMAP 嵌入

### 支持的输出格式

#### SingleCellExperiment

SingleCellExperiment 对象包含：

  - Assays：`counts`（稀疏矩阵）和 `logcounts`（密集矩阵）
  - ColData：细胞元数据
  - RowData：基因元数据
  - ReducedDims：PCA、tSNE 和 UMAP 嵌入

#### Seurat

Seurat 对象包含：

  - Assay：RNA assay，包含 `counts` 和 `data` 槽位
  - MetaData：细胞元数据
  - Reductions：pca、tsne 和 umap

#### h5ad (用于 Scanpy)

h5ad 文件包含：

  - X：原始计数矩阵（保留为稀疏矩阵以提高内存效率）
  - layers：'logcounts' 层中的对数转换数据（保留为稀疏矩阵）
  - obs：细胞元数据
  - var：基因元数据
  - obsm：降维结果（PCA、tSNE、UMAP）

**注意**：`iCellbio2H5ad()` 函数会自动保留稀疏矩阵格式，从而显著节省存储空间（通常比密集矩阵存储减少 70-90% 的文件大小）。

-----

## 示例工作流程

### 标准工作流程（两步转换）

```r
# 加载包
library(ICellbioRpy)

# 读取数据
data <- read1Cellbio("path/to/1cellbio_results.zip")

# 转换为 SingleCellExperiment，用于 Bioconductor 包
sce <- as.SingleCellExperiment(data)

# 或转换为 Seurat 对象，用于 Seurat 函数
seurat <- as.Seurat(data)

# 或转换为 h5ad，用于 Python 中的 Scanpy
as.h5ad(data, "results.h5ad")

# 执行下游分析
# 使用 SingleCellExperiment:
library(scater)
sce <- logNormCounts(sce)
plotPCA(sce, colour_by = "level1class")

# 使用 Seurat:
library(Seurat)
DimPlot(seurat, reduction = "umap", group.by = "level1class")

# 使用 Scanpy (在 Python 中):
# import scanpy as sc
# adata = sc.read_h5ad("results.h5ad")
# sc.pl.umap(adata, color='level1class')
```

### 内存高效工作流程（直接转换并保留稀疏矩阵）

```r
# 加载包
library(ICellbioRpy)

# 对于大型数据集，直接转换为 h5ad 格式
# 这会保留稀疏矩阵格式，并避免创建中间 R 对象
# 实现显著的内存节省和存储效率
iCellbio2H5ad("path/to/1cellbio_results.zip", "results.h5ad")

# 可选：指定自定义临时目录并保留临时文件以进行调试
iCellbio2H5ad("path/to/1cellbio_results.zip", "results.h5ad", 
            temp_dir = "/custom/temp/dir", cleanup = FALSE)

# 生成的 h5ad 文件将比密集存储小得多（减少 70-90%）
# 示例：一个作为密集矩阵约 470 MB 的数据集，作为稀疏矩阵变为约 99 MB

# 然后在 Python 中使用 Scanpy 处理 h5ad 文件：
# import scanpy as sc
# adata = sc.read_h5ad("results.h5ad")
# print(f"Matrix sparsity: {(1 - adata.X.nnz / (adata.shape[0] * adata.shape[1])) * 100:.1f}%")
# sc.pl.umap(adata, color='level1class')
```

-----

## 技术实现

该包使用 `hdf5r` 包处理 HDF5 文件读取，并正确管理：

  - 从 10x 格式转换稀疏矩阵
  - 正确处理所有数据类型的维度
  - 从 HDF5 数据集解析元数据
  - 自动处理 zip 文件解压
  - 使用 `anndata` 包转换为 h5ad 格式
  - 针对大型数据集的内存高效直接转换

-----

## 要求

  - R \>= 4.0
  - Bioconductor 包：SingleCellExperiment
  - CRAN 包：Seurat, hdf5r, jsonlite, Matrix, anndata, reticulate (\>= 1.39.0)

-----

## 许可证

MIT
