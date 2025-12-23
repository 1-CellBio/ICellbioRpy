# ICellbioRpy 快速开始指南

## 📦 安装

```r
# 安装devtools（如果尚未安装）
install.packages("devtools")

# 从GitHub安装ICellbioRpy
devtools::install_github("1-Cellbio/ICellbioRpy")

# 加载包
library(ICellbioRpy)
```

## 🚀 核心功能概览

ICellbioRpy提供完整的单细胞和空间转录组数据格式转换生态系统：

- **读取1Cellbio结果** → `read1Cellbio()`
- **读取Stereo-seq GEF文件** → `read_gef()`
- **转换为h5ad格式** → `iCellbio2H5ad()`, `gef_to_h5ad()`
- **h5ad转R对象** → `h5ad_to_sce()`, `h5ad_to_seurat()`
- **R对象转h5ad** → `seurat_to_h5ad()`
- **空间转录组可视化** → `plot_cells_with_borders()`
- **GMT基因集预处理** → `preprocess_gmt_custom()`
- **Python环境配置** → `configure_python_env()`, `smart_python_config()`

## 🔧 Python环境配置

### 自动配置（推荐）

```r
# 包会自动检测并配置Python环境
library(ICellbioRpy)

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

### 智能交互式配置（推荐）

```r
# 自动检测所有conda环境，并列出包含anndata的环境供选择
smart_python_config(verbose = TRUE, interactive = TRUE)

# 输出示例：
# 📋 Found multiple environments with anndata:
#   1. 1cellbio (anndata 0.12.0)
#   2. atlas (anndata 0.11.3)
#   3. scanpy (anndata 0.10.9)
# Please select environment to use (1-3):
```

如果只有一个环境包含 anndata，会自动选择；如果有多个，会提示用户选择。

### 避免自动安装提示

如果遇到anndata自动安装提示：

```r
# 方法1：在R会话开始时设置
Sys.setenv(RETICULATE_AUTOCONFIGURE = "FALSE")
library(ICellbioRpy)
configure_python_env(conda_env = "your_env")

# 方法2：直接指定环境
configure_python_env(conda_env = "atlas")
check_anndata_available()
```

## 📊 基本使用流程

### 1. 读取1Cellbio数据

```r
# 从zip文件读取1Cellbio结果
data <- read1Cellbio("1Cellbio_results.zip")

# 查看数据结构
class(data)
#> [1] "1CellbioData"

# 查看数据信息
print(data)
```

### 2. 转换为不同格式

#### 转换为SingleCellExperiment

```r
# 方法1：使用新的.1CB函数（推荐，自动检测列名）
sce <- as.SingleCellExperiment.1CB(data)

# 方法2：手动指定列名
sce <- as.SingleCellExperiment.1CB(data,
                                   rownames = "id",        # 基因名列
                                   colnames = "cell_id")   # 细胞名列

# 方法3：查看可用选项后再转换
show_column_options(data)
# 这会显示：
# === Column Detection Results ===
#
# Available gene identifier columns (rownames):
#   id, gene_symbol, gene_name
#   → Detected: id
#
# Available cell identifier columns (colnames):
#   cell_id, barcode, sample_id
#   → Detected: cell_id

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
# 方法1：使用新的.1CB函数（推荐，自动检测列名）
seurat <- as.Seurat.1CB(data)

# 方法2：手动指定列名
seurat <- as.Seurat.1CB(data,
                        rownames = "id",        # 基因名列
                        colnames = "cell_id")   # 细胞名列

# 方法3：关闭自动检测并指定列名
seurat <- as.Seurat.1CB(data,
                        rownames = "specific_gene_col",
                        colnames = "specific_cell_col",
                        auto_detect = FALSE)

# 查看数据
seurat
#> An object of class Seurat
#> 20000 features across 3000 samples

# 访问数据
counts_matrix <- GetAssayData(seurat, layer = "counts")
cell_metadata <- seurat[[]]
pca_coords <- Embeddings(seurat, reduction = "pca")
```

### 3. 直接转换为h5ad格式

```r
# 一步转换：zip → h5ad
iCellbio2H5ad("1cellbio_results.zip", "output.h5ad")

# 带参数控制
iCellbio2H5ad(
  "1cellbio_results.zip",
  "output.h5ad",
  overwrite = FALSE,              # 不覆盖已存在的文件
  name_conflict = "make_unique"   # 命名冲突时自动重命名 ("make_unique" 或 "error")
)

# 检查输出文件
file.exists("output.h5ad")
#> [1] TRUE

# 查看文件大小
file.info("output.h5ad")$size / (1024^2)  # MB
```

**参数说明：**
- `overwrite`: 是否覆盖已存在的输出文件（默认 FALSE）
- `name_conflict`: 命名冲突处理策略
  - `"make_unique"`: 自动重命名冲突的行名列名
  - `"error"`: 在冲突时报错

## 🔄 双向h5ad转换

### 从h5ad到R对象

```r
# h5ad → SingleCellExperiment
sce <- h5ad_to_sce("data.h5ad")

# h5ad → Seurat
seurat <- h5ad_to_seurat("data.h5ad")

# 带参数控制
sce <- h5ad_to_sce(
  "data.h5ad",
  use_x_as = "auto",            # 自动检测 X 层类型 ("auto"/"logcounts"/"counts")
  name_conflict = "make_unique" # 命名冲突处理策略
)

# 查看转换结果
assayNames(sce)
#> [1] "X" "raw"

names(seurat@reductions)
#> [1] "X_pca" "X_umap"
```

**参数说明：**
- `use_x_as`: X 矩阵的解析方式
  - `"auto"`: 自动检测（默认）
  - `"logcounts"`: 作为标准化数据
  - `"counts"`: 作为原始计数
- `name_conflict`: 命名冲突处理策略（"make_unique" 或 "error"）

### 从R对象到h5ad

```r
# Seurat → h5ad（新增覆盖与命名冲突控制）
seurat_to_h5ad(
  seurat_object,
  "seurat_output.h5ad",
  overwrite = FALSE,              # 默认不覆盖已存在文件
  name_conflict = "make_unique"   # 或设置为 "error" 以在命名冲突时报错
)

# 验证转换
file.exists("seurat_output.h5ad")
#> [1] TRUE
```

## 🧬 Stereo-seq 空间转录组支持

ICellbioRpy 支持 Stereo-seq GEF 文件格式的读取和转换，包括细胞分割数据和细胞边界信息。

### 读取 GEF 文件

```r
# 读取 GEF 文件（包含细胞边界）
stereo_data <- read_gef(
  "sample1.gef",
  bin_type = "cell_bins",      # 使用细胞分割数据
  include_cellborder = TRUE    # 包含细胞边界信息
)
```

### 转换为 R 对象

```r
# 转换为 Seurat 对象
seurat <- as.Seurat(stereo_data)

# 空间坐标存储在 reductions 中
seurat@reductions$spatial
#> Coordinate system: spatial

# 细胞边界存储在 @misc 中
seurat@misc$cell_borders

# 转换为 SingleCellExperiment
sce <- as.SingleCellExperiment(stereo_data)

# 空间坐标在 reducedDims 中
reducedDim(sce, "spatial")
```

### 空间可视化

```r
# 绘制带细胞边界的空间细胞图
plot_cells_with_borders(
  stereo_data,
  color_by = "cluster",        # 按聚类着色
  show_borders = TRUE,         # 显示细胞边界
  border_color = "gray",       # 边界颜色
  border_size = 0.5,           # 边界线宽
  point_size = 1               # 细胞点大小
)

# 按基因表达着色
plot_cells_with_borders(
  stereo_data,
  color_by = "EPCAM",
  show_borders = TRUE
)
```

### 直接转换为 H5AD

```r
# GEF → H5AD 直接转换（内存高效）
gef_to_h5ad(
  "../C04042E3.cellbin.gef",
  "output.h5ad",
  bin_type = "cell_bins",
  include_cellborder = TRUE,    # 保留细胞边界
  include_spatial = TRUE,       # 保留空间坐标
  overwrite = FALSE
)

# 或者先读取再转换
stereo_to_h5ad(
  stereo_data,
  "output.h5ad",
  layer = "counts"
)
```

### 高级选项

```r
# 读取特定空间区域
stereo_data <- read_gef(
  "sample.gef",
  region = c(1000, 3000, 1000, 3000),  # minX, maxX, minY, maxY
  max_cells = 10000,                    # 限制细胞数量
  gene_list = c("EPCAM", "KRT8", "VIM") # 仅读取指定基因
)

# 使用方形 bins（不使用细胞分割）
stereo_bins <- read_gef(
  "sample.gef",
  bin_type = "bins",
  bin_size = 50  # 50x50 像素的 bin
)
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

## 🧫 GMT 基因集预处理

ICellbioRpy 提供 GMT（Gene Matrix Transposed）文件预处理功能，支持多种基因 ID 类型映射。

### 基本用法

```r
# 预处理 GMT 文件
preprocess_gmt_custom(
  gmt_file = "pathways.gmt",
  species = "9606",          # 物种 ID（9606 = 人类）
  output_dir = "gesel_output"
)
```

### 支持的物种

| 物种 ID | 物种名称 |
|---------|----------|
| 9606 | 人类 (Homo sapiens) |
| 10090 | 小鼠 (Mus musculus) |
| 10116 | 大鼠 (Rattus norvegicus) |
| 7227 | 果蝇 (Drosophila melanogaster) |
| 6239 | 线虫 (Caenorhabditis elegans) |
| 7955 | 斑马鱼 (Danio rerio) |
| 9598 | 猩猩 (Pan troglodytes) |

### 性能优化：预构建查找表

对于大量 GMT 文件处理，可以预构建基因查找表以提升性能：

```r
# 第一步：预构建所有物种的查找表
prebuild_gene_lookup_tables(
  data_dir = "~/gene_mapping_data",
  output_file = "gene_lookup_tables.rdata"
)

# 第二步：加载到全局环境（创建 master_lookup_tables 变量）
load("gene_lookup_tables.rdata")

# 第三步：调用预处理函数，会自动检测并使用预构建的表
preprocess_gmt_custom("pathways.gmt", species = "9606")

# 批量处理多个 GMT 文件（性能优势明显）
for (gmt in list.files(pattern = "\\.gmt$")) {
  preprocess_gmt_custom(gmt, species = "9606")
}
```

**说明：** 预构建的查找表通过全局环境变量 `master_lookup_tables` 隐式传递给预处理函数。加载后，`preprocess_gmt_custom()` 会自动检测并使用它，无需显式传递参数。

### 输出格式

预处理后会在输出目录生成以下文件：

```
gesel_output/
├── collections.tsv       # 集合元数据
├── sets.tsv              # 基因集定义
├── set2gene.tsv          # 集合到基因映射（差分编码）
└── gene2gene.tsv         # 基因 ID 映射（差分编码）
```

### 高级用法

```r
# 使用自定义映射文件
preprocess_gmt_with_custom_mapping(
  gmt_file = "custom_pathways.gmt",
  species = "10090",  # 小鼠
  gene_mapping_files = list(
    symbol = "mouse_symbols.tsv",
    entrez = "mouse_entrez.tsv",
    ensembl = "mouse_ensembl.tsv"
  ),
  collection_name = "mouse_pathways",
  collection_desc = "Custom mouse pathways",
  output_dir = "mouse_gesel",
  auto_download_missing = TRUE  # 自动下载缺失的映射文件
)
```

## 🛠️ 故障排除

### 列名检测问题

```r
# 如果自动检测失败，首先查看可用选项
show_column_options(data)

# 常见错误和解决方法：
# 错误："Cannot detect gene identifier column"
# 解决：手动指定基因名列
sce <- as.SingleCellExperiment.1CB(data, rownames = "gene_name", colnames = "cell_id")

# 错误："Column 'xxx' not found"
# 解决：查看正确的列名
show_column_options(data)
# 或手动检查
# coldata <- data$experiment$summarized_experiment$column_data$resource
# rowdata <- data$experiment$summarized_experiment$row_data$resource

# 关闭自动检测，完全手动控制
sce <- as.SingleCellExperiment.1CB(data,
                                    rownames = "your_gene_col",
                                    colnames = "your_cell_col",
                                    auto_detect = FALSE)
```

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

```
# 推荐的分析流程
library(ICellbioRpy)

# 1. 配置环境
configure_python_env()

# 2. 读取数据
data <- read1Cellbio("results.zip")

# 3. 根据需要选择格式
if (use_bioconductor) {
  # 新版本：自动检测列名（推荐）
  sce <- as.SingleCellExperiment.1CB(data)

  # 旧版本：手动指定列名
  # sce <- as.SingleCellExperiment(data, rownames = "id", colnames = "cell_id")

  # Bioconductor分析...
} else if (use_seurat) {
  # 新版本：自动检测列名（推荐）
  seurat <- as.Seurat.1CB(data)

  # 旧版本：手动指定列名
  # seurat <- as.Seurat(data, rownames = "id", colnames = "cell_id")

  # Seurat分析...
} else if (use_python) {
  iCellbio2H5ad("results.zip", "analysis.h5ad")
  # Python/Scanpy分析...
}

# 4. 查看可用列选项（如果自动检测失败）
# show_column_options(data)
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

- **详细文档**: `vignette("introduction", package = "ICellbioRpy")`
- **函数帮助**: `?iCellbio2H5ad`, `?h5ad_to_sce`
- **安装指南**: `anndata_installation_guide.md`
- **完整示例**: `README.md`

## 🆘 获取帮助

```r
# 查看包信息
packageVersion("ICellbioRpy")

# 查看会话信息
sessionInfo()

# 报告问题时请包含上述信息
```

## 6. 读取和整合10X MTX格式数据 🧬

`read_10x_mtx_to_h5ad()` 函数可以直接读取多个10X Cell Ranger输出的MTX格式数据，进行基础QC过滤，并整合为h5ad文件：

### 特性
- ✅ 支持压缩(.gz)和非压缩的MTX文件
- ✅ 使用data.table快速读取，依赖少
- ✅ 自动细胞ID重命名避免重复
- ✅ 简单QC过滤（可设置最小counts阈值）
- ✅ 直接输出h5ad格式

### 使用方法

```r
# 1. 准备样本信息CSV文件
sample_info <- data.frame(
  Sample_id = c("sample1", "sample2"),
  mtx_fns = c(
    "/path/to/sample1/matrix.mtx.gz",
    "/path/to/sample2/matrix.mtx.gz"
  ),
  features_fns = c(
    "/path/to/sample1/features.tsv.gz",
    "/path/to/sample2/features.tsv.gz"
  ),
  barcodes_fns = c(
    "/path/to/sample1/barcodes.tsv.gz",
    "/path/to/sample2/barcodes.tsv.gz"
  )
)
write.csv(sample_info, "samples.csv", row.names = FALSE)

# 2. 读取和整合数据
read_10x_mtx_to_h5ad(
  csv_file = "samples.csv",
  output_h5ad = "integrated_data.h5ad",
  min_counts_per_cell = 200,  # QC过滤阈值
  verbose = TRUE
)
```

### CSV文件格式要求
- `Sample_id`: 样本标识符
- `mtx_fns`: matrix.mtx文件路径
- `features_fns`: features.tsv或genes.tsv文件路径  
- `barcodes_fns`: barcodes.tsv文件路径

### 输出
- 细胞ID格式: `{Sample_id}_{原始barcode}`
- 矩阵格式: 基因 × 细胞 (输出时转换为细胞 × 基因)
- 自动过滤低质量细胞(总counts < 阈值)

---

**提示**: 如果您是第一次使用，建议先阅读 `vignette("introduction")` 获取完整的使用指南。