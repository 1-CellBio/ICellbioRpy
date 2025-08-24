
在这个教程中，我们将学习如何使用`ICellbioRpy`包将Seurat对象转换为H5AD格式，以便在Python的scanpy生态系统中进行分析。这在需要跨语言分析或使用Python特有工具时非常有用。

# 📋 前提条件

## 必需的软件和包

1. **R环境**: R (≥ 4.0.0)
2. **Python环境**: Python (≥ 3.7) + anndata包
3. **R包**: Seurat, ICellbioRpy, reticulate, Matrix

```r
# 检查R包安装状态
required_r_packages <- c("Seurat", "Matrix", "reticulate")
missing_r_packages <- required_r_packages[!sapply(required_r_packages, requireNamespace, quietly = TRUE)]

if (length(missing_r_packages) > 0) {
  cat("需要安装以下R包:", paste(missing_r_packages, collapse = ", "), "\n")
} else {
  cat("✓ 所有必需的R包已安装\n")
}
```

## Python环境检查

```{bash check-python, eval=FALSE}
# 在终端中检查Python和anndata
python --version
python -c "import anndata; print('anndata version:', anndata.__version__)"
```

如果anndata未安装：

```{bash install-anndata, eval=FALSE}
# 使用pip安装
pip install anndata

# 或使用conda安装  
conda install -c conda-forge anndata
```

# 🚀 开始教程

## 第1步: 环境设置

```r
# 加载必需的包
library(ICellbioRpy)
library(Seurat)
library(Matrix)
library(dplyr)

# 显示版本信息
cat("包版本信息:\n")
cat("ICellbioRpy:", as.character(packageVersion("ICellbioRpy")), "\n")
cat("Seurat:", as.character(packageVersion("Seurat")), "\n")
cat("Matrix:", as.character(packageVersion("Matrix")), "\n")
```

## 第2步: 配置Python环境（关键步骤！）

这是最重要的步骤，必须正确配置才能成功转换。

```r
# 方法1: 使用默认Python环境（推荐初学者）
configure_python_env(verbose = TRUE)

# 方法2: 指定conda环境（如果您使用conda）
# configure_python_env(conda_env = "your_env_name", verbose = TRUE)

# 方法3: 指定Python路径（如果需要特定版本）
# configure_python_env(python_path = "/usr/local/bin/python3", verbose = TRUE)

# 验证配置成功
check_anndata_available()
```

**💡 配置说明:**

- `verbose = TRUE`: 显示详细配置信息，便于诊断问题
- 看到"✓ anndata is available"表示配置成功
- 如果失败，请检查Python环境和anndata安装

**🔍 配置排查:**

```r
# 如果配置失败，运行以下诊断命令
if (FALSE) {  # 设为TRUE来运行诊断
  # 查看当前Python配置
  reticulate::py_config()
  
  # 列出可用的conda环境
  reticulate::conda_list()
  
  # 检查特定环境中的包
  # reticulate::py_list_packages()
}
```

## 第3步: 准备Seurat对象

您可以使用现有的Seurat对象，或创建一个示例对象进行测试。

### 选项A: 加载现有的Seurat对象

```r
# 如果您有保存的Seurat对象
seurat_file <- "path/to/your/seurat_object.rds"

if (file.exists(seurat_file)) {
  seurat_obj <- readRDS(seurat_file)
  cat("✓ 成功加载Seurat对象\n")
} else {
  cat("✗ 找不到Seurat文件，请检查路径或使用下面的示例数据\n")
}
```

### 选项B: 创建示例Seurat对象（用于测试）

```r
# 如果没有现成的Seurat对象，创建一个示例对象
if (!exists("seurat_obj")) {
  # 使用Seurat内置的示例数据
  data("pbmc_small")
  seurat_obj <- pbmc_small
  
  # 确保对象包含必要的信息
  if (!"nFeature_RNA" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj[["nFeature_RNA"]] <- Matrix::colSums(seurat_obj@assays$RNA@counts > 0)
    seurat_obj[["nCount_RNA"]] <- Matrix::colSums(seurat_obj@assays$RNA@counts)
  }
  
  cat("✓ 使用示例Seurat对象 (pbmc_small)\n")
}
```

## 第4步: 检查Seurat对象

在转换前，让我们检查Seurat对象的结构：

```r
cat("=== Seurat对象信息 ===\n")
print(seurat_obj)

cat("\n=== 基本统计 ===\n")
cat("细胞数量:", ncol(seurat_obj), "\n")
cat("基因数量:", nrow(seurat_obj), "\n")
cat("默认assay:", DefaultAssay(seurat_obj), "\n")

cat("\n=== 可用的Assays ===\n")
cat("Assays:", names(seurat_obj@assays), "\n")

# 检查每个assay的slots
for (assay_name in names(seurat_obj@assays)) {
  assay_obj <- seurat_obj@assays[[assay_name]]
  available_slots <- c()
  if (length(assay_obj@counts) > 0) available_slots <- c(available_slots, "counts")
  if (length(assay_obj@data) > 0) available_slots <- c(available_slots, "data")
  if (length(assay_obj@scale.data) > 0) available_slots <- c(available_slots, "scale.data")
  
  cat(sprintf("%s assay可用slots: %s\n", assay_name, paste(available_slots, collapse = ", ")))
}

cat("\n=== 元数据列 ===\n")
cat("元数据列数:", ncol(seurat_obj@meta.data), "\n")
cat("列名:", paste(colnames(seurat_obj@meta.data), collapse = ", "), "\n")

cat("\n=== 降维结果 ===\n")
if (length(seurat_obj@reductions) > 0) {
  for (reduction_name in names(seurat_obj@reductions)) {
    dims <- dim(seurat_obj@reductions[[reduction_name]]@cell.embeddings)
    cat(sprintf("%s: %d维度\n", reduction_name, dims[2]))
  }
} else {
  cat("无降维结果\n")
}
```

## 第5步: 执行转换

现在进行实际的转换操作：

```r
# 设置输出文件路径
output_h5ad <- "seurat_converted.h5ad"

cat("开始转换Seurat对象到H5AD格式...\n")

# 执行转换
# 重要参数说明：
# - default_assay: 指定主要的assay（通常是"RNA"）
# - layer: 指定使用哪个数据层（"data"表示标准化数据，"counts"表示原始计数）
# - include_reductions: 是否包含降维结果
seurat_to_h5ad(
  seurat_obj = seurat_obj,
  output_file = output_h5ad,
  default_assay = "RNA",           # 使用RNA assay作为主要数据
  layer = "data",                  # 使用标准化数据作为主矩阵
  include_reductions = TRUE,       # 包含降维结果
  overwrite = FALSE,               # 新增：默认不覆盖已存在文件
  name_conflict = "make_unique",  # 新增：命名冲突策略（或设为 "error"）
  verbose = TRUE                   # 显示详细进度
)

cat("✓ 转换完成！\n")
```

**🔍 参数详解（更新）:**

- **default_assay**: 指定哪个assay作为主要数据源
- **layer**: 
  - `"data"`: 标准化后的数据（推荐）
  - `"counts"`: 原始计数数据
- **include_reductions**: 是否包含PCA、UMAP等降维结果
- **overwrite**: 若输出文件已存在，`FALSE`（默认）时报错并拒绝覆盖；设为 `TRUE` 允许覆盖
- **name_conflict**: 当细胞/基因名称有重复时的处理策略：
  - `"make_unique"`（默认）：自动添加后缀保证名称唯一
  - `"error"`：发现冲突即报错，提醒用户显式处理
- **verbose**: 显示转换过程的详细信息

## 第6步: 验证转换结果

检查生成的H5AD文件：

```r
# 检查文件是否成功创建
if (file.exists(output_h5ad)) {
  cat("✓ H5AD文件创建成功!\n")
  cat("文件路径:", output_h5ad, "\n")
  
  # 显示文件大小
  file_size_mb <- round(file.size(output_h5ad) / 1024^2, 2)
  cat("文件大小:", file_size_mb, "MB\n")
  
  # 尝试在R中读取验证（需要reticulate和anndata）
  tryCatch({
    # 导入anndata模块
    ad <- reticulate::import("anndata")
    
    # 读取H5AD文件
    adata <- ad$read_h5ad(output_h5ad)
    
    cat("\n=== H5AD文件信息 ===\n")
    cat("数据维度:", adata$shape[1], "细胞 ×", adata$shape[2], "基因\n")
    
    # 检查观察值（细胞）注释
    cat("细胞注释列:", paste(names(adata$obs), collapse = ", "), "\n")
    
    # 检查变量（基因）注释
    cat("基因注释列:", paste(names(adata$var), collapse = ", "), "\n")
    
    # 检查降维结果
    if (length(names(adata$obsm)) > 0) {
      cat("降维结果:", paste(names(adata$obsm), collapse = ", "), "\n")
    } else {
      cat("无降维结果\n")
    }
    
    # 检查层（layers）
    if (length(names(adata$layers)) > 0) {
      cat("数据层:", paste(names(adata$layers), collapse = ", "), "\n")
    }
    
  }, error = function(e) {
    cat("无法在R中读取H5AD文件进行验证，但文件已创建\n")
    cat("错误信息:", e$message, "\n")
  })
  
} else {
  cat("✗ H5AD文件创建失败\n")
}
```

## 第7步: 高级转换选项

### 使用不同的数据层

```r
# 选项1: 使用原始计数数据
seurat_to_h5ad(
  seurat_obj = seurat_obj,
  output_file = "seurat_counts.h5ad",
  default_assay = "RNA",
  layer = "counts",              # 使用原始计数
  include_reductions = TRUE,
  verbose = TRUE
)

# 选项2: 不包含降维结果（节省空间）
seurat_to_h5ad(
  seurat_obj = seurat_obj,
  output_file = "seurat_no_reductions.h5ad",
  default_assay = "RNA",
  layer = "data",
  include_reductions = FALSE,    # 不包含降维结果
  verbose = TRUE
)

# 选项3: 如果有多个assays，指定不同的assay
if ("SCT" %in% names(seurat_obj@assays)) {
  seurat_to_h5ad(
    seurat_obj = seurat_obj,
    output_file = "seurat_sct.h5ad",
    default_assay = "SCT",       # 使用SCT assay
    layer = "data",
    include_reductions = TRUE,
    verbose = TRUE
  )
}
```

### 转换特定的细胞子集

```r
# 如果您只想转换特定的细胞群体
if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  # 选择特定聚类（例如聚类0和1）
  subset_obj <- subset(seurat_obj, subset = seurat_clusters %in% c("0", "1"))
  
  seurat_to_h5ad(
    seurat_obj = subset_obj,
    output_file = "seurat_subset.h5ad",
    default_assay = "RNA",
    layer = "data",
    include_reductions = TRUE,
    verbose = TRUE
  )
  
  cat("✓ 子集转换完成，包含", ncol(subset_obj), "个细胞\n")
}
```

# 🐍 在Python中使用转换后的数据

转换完成后，您可以在Python中加载和分析数据：

```{python python-analysis, eval=FALSE}
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 设置scanpy参数
sc.settings.verbosity = 3  # 显示详细信息
sc.settings.set_figure_params(dpi=80, facecolor='white')

# 读取转换后的H5AD文件
adata = sc.read_h5ad('seurat_converted.h5ad')

# 查看数据基本信息
print(f"数据维度: {adata.shape}")
print(f"细胞注释: {list(adata.obs.columns)}")
print(f"基因注释: {list(adata.var.columns)}")
print(f"降维结果: {list(adata.obsm.keys())}")

# 如果有UMAP结果，直接可视化
if 'X_umap' in adata.obsm:
    sc.pl.umap(adata, color='seurat_clusters', legend_loc='on data')

# 进行scanpy特有的分析
# 例如：计算邻居图（如果还没有）
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# 运行leiden聚类
sc.tl.leiden(adata, resolution=0.5)

# 可视化leiden聚类结果
sc.pl.umap(adata, color='leiden')

# 寻找标记基因
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False)
```

# 🔧 常见问题排查

## 问题1: Python环境问题

**错误信息**: "Python environment configuration failed"

```r
# 诊断步骤1: 检查Python配置
cat("当前Python配置:\n")
reticulate::py_config()

# 诊断步骤2: 检查anndata安装
tryCatch({
  ad <- reticulate::import("anndata")
  cat("✓ anndata可用，版本:", ad$`__version__`, "\n")
}, error = function(e) {
  cat("✗ anndata不可用:", e$message, "\n")
  cat("请在终端运行: pip install anndata\n")
})

# 诊断步骤3: 重新配置Python环境
# configure_python_env(conda_env = "base", verbose = TRUE)
```

## 问题2: 内存不足

**错误信息**: "Memory allocation failed"

```r
# 检查对象大小
cat("Seurat对象大小:", format(object.size(seurat_obj), units = "MB"), "\n")

# 内存清理
gc()

# 如果内存不足，可以尝试：
# 1. 转换子集数据
# 2. 不包含降维结果
# 3. 只使用必要的assay
```

## 问题3: 文件写入权限问题

**错误信息**: "Permission denied" 或 "Cannot write file"

```r
# 检查当前工作目录
cat("当前工作目录:", getwd(), "\n")

# 检查写入权限
test_file <- "test_write.txt"
tryCatch({
  writeLines("test", test_file)
  unlink(test_file)
  cat("✓ 当前目录有写入权限\n")
}, error = function(e) {
  cat("✗ 当前目录无写入权限\n")
  cat("建议使用绝对路径或更改工作目录\n")
})
```

## 问题4: Seurat对象格式问题

**错误信息**: "Invalid Seurat object" 或数据缺失

```r
# 检查Seurat对象完整性
cat("=== Seurat对象诊断 ===\n")

# 检查assays
if (length(seurat_obj@assays) == 0) {
  cat("✗ 没有可用的assays\n")
} else {
  cat("✓ 可用assays:", names(seurat_obj@assays), "\n")
}

# 检查默认assay的数据
default_assay <- DefaultAssay(seurat_obj)
assay_obj <- seurat_obj@assays[[default_assay]]

if (length(assay_obj@data) == 0 && length(assay_obj@counts) == 0) {
  cat("✗ 默认assay中没有数据\n")
} else {
  cat("✓ 默认assay包含数据\n")
}

# 检查元数据
if (nrow(seurat_obj@meta.data) == 0) {
  cat("✗ 没有元数据\n")
} else {
  cat("✓ 元数据包含", nrow(seurat_obj@meta.data), "个细胞\n")
}
```

# 📊 数据质量对比

让我们比较转换前后的数据一致性：

```r
# 这个代码块需要成功转换后运行
if (file.exists(output_h5ad)) {
  tryCatch({
    # 读取转换后的H5AD文件
    ad <- reticulate::import("anndata")
    adata <- ad$read_h5ad(output_h5ad)
    
    cat("=== 数据一致性检查 ===\n")
    
    # 检查维度
    seurat_dims <- c(ncol(seurat_obj), nrow(seurat_obj))
    h5ad_dims <- c(adata$shape[1], adata$shape[2])
    
    cat("Seurat维度 (细胞, 基因):", seurat_dims[1], ",", seurat_dims[2], "\n")
    cat("H5AD维度 (细胞, 基因):", h5ad_dims[1], ",", h5ad_dims[2], "\n")
    
    if (identical(seurat_dims, h5ad_dims)) {
      cat("✓ 维度一致\n")
    } else {
      cat("⚠ 维度不一致\n")
    }
    
    # 检查细胞名称
    seurat_cells <- colnames(seurat_obj)
    h5ad_cells <- adata$obs_names$to_list()
    
    if (length(intersect(seurat_cells, h5ad_cells)) == length(seurat_cells)) {
      cat("✓ 细胞名称一致\n")
    } else {
      cat("⚠ 细胞名称不完全一致\n")
    }
    
    # 检查基因名称
    seurat_genes <- rownames(seurat_obj)
    h5ad_genes <- adata$var_names$to_list()
    
    if (length(intersect(seurat_genes, h5ad_genes)) == length(seurat_genes)) {
      cat("✓ 基因名称一致\n")
    } else {
      cat("⚠ 基因名称不完全一致\n")
    }
    
  }, error = function(e) {
    cat("无法进行一致性检查:", e$message, "\n")
  })
}
```

# 📚 下一步分析

成功转换为H5AD格式后，您可以：

## 在Python中进行的分析

1. **质量控制**: 使用scanpy进行更精细的QC
2. **批次校正**: 使用scanorama、harmony等
3. **轨迹分析**: 使用scvelo进行RNA velocity分析
4. **空间分析**: 使用squidpy进行空间转录组分析
5. **多组学整合**: 整合ATAC-seq、蛋白质组学数据

## 推荐工具和包

- **scanpy**: 核心单细胞分析包
- **scvelo**: RNA velocity分析
- **cellrank**: 细胞命运预测
- **squidpy**: 空间转录组学
- **anndata**: 数据结构核心

# 💡 最佳实践建议

1. **数据备份**: 转换前备份原始Seurat对象
2. **版本记录**: 记录使用的软件版本
3. **参数记录**: 记录转换参数以便重现
4. **质量检查**: 转换后验证数据完整性
5. **文档记录**: 详细记录分析流程

## 性能优化技巧

```r
# 1. 对于大数据集，考虑分批处理
# 2. 使用SSD存储以提高I/O速度
# 3. 确保足够的RAM
# 4. 使用压缩格式减少文件大小

# 检查系统资源
cat("可用内存:", round(as.numeric(system("free -m | grep '^Mem:' | awk '{print $7}'", intern = TRUE)) / 1024, 1), "GB\n")
```

---

**需要帮助？**

- 查看函数文档: `?seurat_to_h5ad`
- 查看scanpy文档: https://scanpy.readthedocs.io/
- GitHub Issues: [项目地址]
- 联系维护者: [邮箱地址]

**相关教程:**
- 教程1: 1CellBio转H5AD
- 教程2: 1CellBio转Seurat  
- 教程4: 10X MTX数据整合