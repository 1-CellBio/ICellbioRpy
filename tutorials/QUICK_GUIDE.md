# 🚀 ICellbioRpy 教程快速使用指南

本指南帮助您快速开始使用ICellbioRpy的各种教程。

## 📋 开始前的检查清单

### 🚨 重要提醒
**如果您是完全初学者，请先完成环境安装教程：**
- Windows用户：[教程0a: Windows环境完整安装指南](tutorial_0_environment_setup_windows.html)
- macOS用户：[教程0b: macOS环境完整安装指南](tutorial_0_environment_setup_macos.html)

### ✅ 必需软件（有经验用户检查清单）
- [ ] R (≥ 4.0.0) 已安装
- [ ] Python (≥ 3.7) 已安装  
- [ ] anndata Python包已安装 (`pip install anndata`)
- [ ] RStudio 或其他R IDE（推荐）
- [ ] ICellbioRpy包已安装

### ✅ 快速环境检查
```r
# 1. 加载ICellbioRpy
library(ICellbioRpy)

# 2. 运行快速测试
source("tutorials/test_tutorials.R")

# 3. 如果所有项目显示 ✓，您就可以开始了！
```

**如果测试失败，请参考对应的环境安装教程进行配置。**

## 🎯 选择合适的教程

### 根据您的数据类型选择：

| 您有什么数据？ | 想要什么格式？ | 推荐教程 | 难度 |
|---------------|---------------|----------|------|
| 1CellBio ZIP文件 | H5AD (Python分析) | [教程1](tutorial_1_1cellbio_to_h5ad.html) | ⭐⭐ |
| 1CellBio ZIP文件 | Seurat (R分析) | [教程2](tutorial_2_1cellbio_to_seurat.html) | ⭐⭐ |
| Seurat对象 | H5AD (转到Python) | [教程3](tutorial_3_seurat_to_h5ad.html) | ⭐⭐⭐ |
| 10X MTX文件 | H5AD (多样本整合) | [教程4](tutorial_4_10x_mtx_to_h5ad.html) | ⭐⭐⭐⭐ |

## ⚡ 快速开始代码片段

### 基本Python环境设置（所有教程必需）
```r
library(ICellbioRpy)

# 配置Python环境
configure_python_env(verbose = TRUE)

# 验证环境
check_anndata_available()
# 应该看到: ✓ anndata available (version: X.X.X)
```

### 1. 1CellBio → H5AD
```r
# 最简单的转换
iCellbio2H5ad(
  zip_path = "path/to/your/1cellbio_results.zip",
  h5ad_path = "output.h5ad"
)
```

### 2. 1CellBio → Seurat
```r
# 读取数据
data <- read1Cellbio("path/to/your/1cellbio_results.zip")

# 转换为Seurat
seurat_obj <- as.Seurat.1CB(data)
```

### 3. Seurat → H5AD
```r
# 转换Seurat对象
seurat_to_h5ad(
  seurat_obj = your_seurat_object,
  output_file = "seurat_data.h5ad"
)
```

### 4. 10X MTX → H5AD
```r
# 准备样本信息CSV
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

## 🔧 常见问题快速解决

### Python环境问题
```r
# 检查Python配置
reticulate::py_config()

# 指定conda环境
configure_python_env(conda_env = "1cellbio")

# 或指定Python路径
configure_python_env(python_path = "/usr/local/bin/python3")
```

### 内存不足
```r
# 清理内存
gc()

# 检查对象大小
object.size(your_object)

# 对于大数据集，考虑：
# 1. 提高QC阈值
# 2. 分批处理
# 3. 使用云计算资源
```

### 文件路径问题
```r
# 使用绝对路径
file.path("/full/path/to/your/file")

# 检查文件是否存在
file.exists("your/file/path")

# 检查当前工作目录
getwd()
```

## 🐍 在Python中继续分析

教程转换完成后，在Python中使用：

```python
import scanpy as sc
import pandas as pd

# 读取H5AD文件
adata = sc.read_h5ad('your_file.h5ad')

# 基本信息
print(f"数据维度: {adata.shape}")
print(f"样本信息: {adata.obs.columns.tolist()}")

# 标准分析流程
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.pp.scale(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)

# 可视化
sc.pl.umap(adata, color='leiden')
```

## 📖 学习路径建议

### 初学者路径
1. **先学教程1** - 掌握基本的格式转换概念
2. **再学教程2** - 了解R中的单细胞分析
3. **然后学教程3** - 理解跨语言数据交换
4. **最后学教程4** - 掌握复杂的多样本整合

### 高级用户路径
- 直接跳到需要的教程
- 重点关注参数优化部分
- 学习性能优化技巧

## 🆘 获取帮助

1. **查看详细教程**: 每个.Rmd文件都有完整说明
2. **运行测试脚本**: `source("tutorials/test_tutorials.R")`
3. **查看函数帮助**: `?function_name`
4. **查看GitHub Issues**: [项目地址]

## 💡 专业提示

1. **数据备份**: 转换前备份原始数据
2. **参数记录**: 记录转换参数以便重现
3. **版本控制**: 记录软件版本信息
4. **批量处理**: 大项目使用脚本自动化
5. **资源监控**: 关注内存和存储使用

---

**准备好了吗？选择一个教程开始您的单细胞数据转换之旅！** 🧬✨