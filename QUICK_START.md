# 1CellbioRpy 快速开始指南

## 📦 安装

```r
# 安装devtools（如果尚未安装）
install.packages("devtools")

# 从GitHub安装1CellbioRpy
devtools::install_github("your_username/1CellbioRpy")

# 加载包
library(1CellbioRpy)
```

## 🚀 核心功能概览

1CellbioRpy提供完整的单细胞数据格式转换生态系统：

- **读取1Cellbio结果** → `read1Cellbio()`
- **转换为h5ad格式** → `iCellbio2H5ad()`
- **h5ad转R对象** → `h5ad_to_sce()`, `h5ad_to_seurat()`
- **R对象转h5ad** → `seurat_to_h5ad()`
- **Python环境配置** → `configure_python_env()`

## 🔧 Python环境配置

### 自动配置（推荐）

```r
# 包会自动检测并配置Python环境
library(1CellbioRpy)

# 验证配置是否成功
check_anndata_available()
```

### 手动配置

```r
# 使用特定conda环境
configure_python_env(conda_env = "scanpy")

# 使用特定Python路径
configure_python_env(python_path = "/usr/local/bin/python3")

# 详细输出（用于调试）
configure_python_env(verbose = TRUE)
```

### 避免自动安装提示

如果遇到anndata自动安装提示：

```r
# 方法1：在R会话开始时设置
Sys.setenv(RETICULATE_AUTOCONFIGURE = "FALSE")
library(1CellbioRpy)
configure_python_env(conda_env = "your_env")

# 方法2：直接指定环境
configure_python_env(conda_env = "atlas")
check_anndata_available()
```

## 📊 基本使用流程

### 1. 读取1Cellbio数据

```r
# 从zip文件读取1Cellbio结果
data <- read1Cellbio("1cellbio_results.zip")

# 查看数据结构
class(data)
#> [1] "1CellbioData"

# 查看数据信息
print(data)
```

### 2. 转换为不同格式

#### 转换为SingleCellExperiment

```r
# 转换为SCE对象
sce <- as.SingleCellExperiment(data)

# 查看数据
sce
#> class: SingleCellExperiment 
#> dim: 20000 3000 
#> assays(2): counts logcounts
#> reducedDims(3): PCA TSNE UMAP

# 访问数据
counts_matrix <- counts(sce)
cell_metadata <- colData(sce)
gene_metadata <- rowData(sce)
pca_coords <- reducedDim(sce, "PCA")
```

#### 转换为Seurat对象

```r
# 转换为Seurat对象
seurat <- as.Seurat(data)

# 查看数据
seurat
#> An object of class Seurat 
#> 20000 features across 3000 samples

# 访问数据
counts_matrix <- GetAssayData(seurat, slot = "counts")
cell_metadata <- seurat[[]]
pca_coords <- Embeddings(seurat, reduction = "pca")
```

### 3. 直接转换为h5ad格式

```r
# 一步转换：zip → h5ad
iCellbio2H5ad("1cellbio_results.zip", "output.h5ad")

# 检查输出文件
file.exists("output.h5ad")
#> [1] TRUE

# 查看文件大小
file.info("output.h5ad")$size / (1024^2)  # MB
```

## 🔄 双向h5ad转换

### 从h5ad到R对象

```r
# h5ad → SingleCellExperiment
sce <- h5ad_to_sce("data.h5ad")

# h5ad → Seurat
seurat <- h5ad_to_seurat("data.h5ad")

# 查看转换结果
assayNames(sce)
#> [1] "X" "raw"

names(seurat@reductions)
#> [1] "X_pca" "X_umap"
```

### 从R对象到h5ad

```r
# Seurat → h5ad
seurat_to_h5ad(seurat_object, "seurat_output.h5ad")

# 验证转换
file.exists("seurat_output.h5ad")
#> [1] TRUE
```

## 🔬 与分析工具集成

### 使用Bioconductor工具

```r
library(scater)
library(scran)

# 质量控制
sce <- addPerCellQC(sce)
sce <- addPerFeatureQC(sce)

# 可视化
plotPCA(sce, colour_by = "level1class")
plotUMAP(sce, colour_by = "total_counts")
```

### 使用Seurat工具

```r
library(Seurat)

# 标准Seurat流程
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)

# 可视化
DimPlot(seurat, reduction = "umap", group.by = "level1class")
FeaturePlot(seurat, features = c("gene1", "gene2"))
```

### 使用Python/Scanpy

转换为h5ad后在Python中使用：

```python
import scanpy as sc
import pandas as pd

# 读取h5ad文件
adata = sc.read_h5ad("output.h5ad")

# 基本信息
print(adata)

# 可视化
sc.pl.umap(adata, color='level1class')
sc.pl.violin(adata, keys=['gene1'], groupby='level1class')

# 差异表达分析
sc.tl.rank_genes_groups(adata, 'level1class')
sc.pl.rank_genes_groups(adata)
```

## 📁 批量处理

### 批量转换多个文件

```r
# 获取所有zip文件
zip_files <- list.files(pattern = "*.zip", full.names = TRUE)

# 批量转换为h5ad
for (zip_file in zip_files) {
  output_name <- gsub(".zip", ".h5ad", basename(zip_file))
  cat("转换:", zip_file, "→", output_name, "\n")
  iCellbio2H5ad(zip_file, output_name)
}
```

### 批量转换h5ad文件

```r
# 获取所有h5ad文件
h5ad_files <- list.files(pattern = "*.h5ad", full.names = TRUE)

# 批量转换为SCE对象
sce_list <- lapply(h5ad_files, function(file) {
  cat("读取:", file, "\n")
  h5ad_to_sce(file)
})

names(sce_list) <- gsub(".h5ad", "", basename(h5ad_files))
```

## 🛠️ 故障排除

### Python环境问题

```r
# 检查Python配置
reticulate::py_config()

# 检查anndata可用性
check_anndata_available()

# 重新配置环境
configure_python_env(conda_env = "base", verbose = TRUE)
```

### 内存问题

```r
# 检查对象大小
object.size(sce)

# 对于大数据集，考虑分批处理
# 或直接使用h5ad格式在Python中处理
```

### 文件路径问题

```r
# 使用绝对路径
file_path <- file.path(getwd(), "data.zip")
iCellbio2H5ad(file_path, "output.h5ad")

# 检查文件是否存在
if (!file.exists("data.zip")) {
  stop("文件不存在: data.zip")
}
```

## 💡 最佳实践

### 1. 工作流程建议

```r
# 推荐的分析流程
library(1CellbioRpy)

# 1. 配置环境
configure_python_env()

# 2. 读取数据
data <- read1Cellbio("results.zip")

# 3. 根据需要选择格式
if (use_bioconductor) {
  sce <- as.SingleCellExperiment(data)
  # Bioconductor分析...
} else if (use_seurat) {
  seurat <- as.Seurat(data)
  # Seurat分析...
} else if (use_python) {
  iCellbio2H5ad("results.zip", "analysis.h5ad")
  # Python/Scanpy分析...
}
```

### 2. 性能优化

```r
# 对于大数据集
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB

# 使用稀疏矩阵
library(Matrix)
# 包会自动保持稀疏性
```

### 3. 数据验证

```r
# 转换后验证数据完整性
original_dims <- dim(counts(sce))
converted_sce <- h5ad_to_sce("temp.h5ad")
new_dims <- dim(counts(converted_sce))

identical(original_dims, new_dims)
#> [1] TRUE
```

## 📚 更多资源

- **详细文档**: `vignette("introduction", package = "1CellbioRpy")`
- **函数帮助**: `?iCellbio2H5ad`, `?h5ad_to_sce`
- **安装指南**: `anndata_installation_guide.md`
- **完整示例**: `README.md`

## 🆘 获取帮助

```r
# 查看包信息
packageVersion("1CellbioRpy")

# 查看会话信息
sessionInfo()

# 报告问题时请包含上述信息
```

---

**提示**: 如果您是第一次使用，建议先阅读 `vignette("introduction")` 获取完整的使用指南。