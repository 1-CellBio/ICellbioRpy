# ICellbioRpy Quick Start Guide

## 📦 Installation

```r
# Install devtools if not already installed
install.packages("devtools")

# Install ICellbioRpy from GitHub
devtools::install_github("1-Cellbio/ICellbioRpy")

# Load the package
library(ICellbioRpy)
```

## 🚀 Core Features Overview

ICellbioRpy provides a complete single-cell data format conversion ecosystem:

- **Read 1Cellbio Results** → `read1Cellbio()`
- **Convert to h5ad Format** → `iCellbio2H5ad()`
- **h5ad to R Objects** → `h5ad_to_sce()`, `h5ad_to_seurat()`
- **R Objects to h5ad** → `seurat_to_h5ad()`
- **Python Environment Configuration** → `configure_python_env()`

## 🔧 Python Environment Configuration

### Automatic Configuration (Recommended)

```r
# The package will automatically detect and configure Python environment
library(ICellbioRpy)

# Verify configuration success
check_anndata_available()
```

### Manual Configuration

```r
# Use specific conda environment
configure_python_env(conda_env = "scanpy")

# Use specific Python path
configure_python_env(python_path = "/usr/local/bin/python3")

# Verbose output (for debugging)
configure_python_env(verbose = TRUE)
```

### Avoid Automatic Installation Prompts

If you encounter anndata automatic installation prompts:

```r
# Method 1: Set at the beginning of R session
Sys.setenv(RETICULATE_AUTOCONFIGURE = "FALSE")
library(ICellbioRpy)
configure_python_env(conda_env = "your_env")

# Method 2: Directly specify environment
configure_python_env(conda_env = "atlas")
check_anndata_available()
```

## 📊 Basic Usage Workflow

### 1. Read 1Cellbio Data

```r
# Read 1Cellbio results from zip file
data <- read1Cellbio("1cellbio_results.zip")

# Check data structure
class(data)
#> [1] "1CellbioData"

# View data information
print(data)
```

### 2. Convert to Different Formats

#### Convert to SingleCellExperiment

```r
# Convert to SCE object
sce <- as.SingleCellExperiment(data)

# View data
sce
#> class: SingleCellExperiment 
#> dim: 20000 3000 
#> assays(2): counts logcounts
#> reducedDims(3): PCA TSNE UMAP

# Access data
counts_matrix <- counts(sce)
cell_metadata <- colData(sce)
gene_metadata <- rowData(sce)
pca_coords <- reducedDim(sce, "PCA")
```

#### Convert to Seurat Object

```r
# Convert to Seurat object
seurat <- as.Seurat(data)

# View data
seurat
#> An object of class Seurat 
#> 20000 features across 3000 samples

# Access data
counts_matrix <- GetAssayData(seurat, slot = "counts")
cell_metadata <- seurat[[]]
pca_coords <- Embeddings(seurat, reduction = "pca")
```

### 3. Direct Conversion to h5ad Format

```r
# One-step conversion: zip → h5ad
iCellbio2H5ad("1cellbio_results.zip", "output.h5ad")

# Check output file
file.exists("output.h5ad")
#> [1] TRUE

# Check file size
file.info("output.h5ad")$size / (1024^2)  # MB
```

## 🔄 Bidirectional h5ad Conversion

### From h5ad to R Objects

```r
# h5ad → SingleCellExperiment
sce <- h5ad_to_sce("data.h5ad")

# h5ad → Seurat
seurat <- h5ad_to_seurat("data.h5ad")

# Check conversion results
assayNames(sce)
#> [1] "X" "raw"

names(seurat@reductions)
#> [1] "X_pca" "X_umap"
```

### From R Objects to h5ad

```r
# Seurat → h5ad
seurat_to_h5ad(seurat_object, "seurat_output.h5ad")

# Verify conversion
file.exists("seurat_output.h5ad")
#> [1] TRUE
```

## 🔬 Integration with Analysis Tools

### Using Bioconductor Tools

```r
library(scater)
library(scran)

# Quality control
sce <- addPerCellQC(sce)
sce <- addPerFeatureQC(sce)

# Visualization
plotPCA(sce, colour_by = "level1class")
plotUMAP(sce, colour_by = "total_counts")
```

### Using Seurat Tools

```r
library(Seurat)

# Standard Seurat workflow
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)

# Visualization
DimPlot(seurat, reduction = "umap", group.by = "level1class")
FeaturePlot(seurat, features = c("gene1", "gene2"))
```

### Using Python/Scanpy

After converting to h5ad, use in Python:

```python
import scanpy as sc
import pandas as pd

# Read h5ad file
adata = sc.read_h5ad("output.h5ad")

# Basic information
print(adata)

# Visualization
sc.pl.umap(adata, color='level1class')
sc.pl.violin(adata, keys=['gene1'], groupby='level1class')

# Differential expression analysis
sc.tl.rank_genes_groups(adata, 'level1class')
sc.pl.rank_genes_groups(adata)
```

## 📁 Batch Processing

### Batch Convert Multiple Files

```r
# Get all zip files
zip_files <- list.files(pattern = "*.zip", full.names = TRUE)

# Batch convert to h5ad
for (zip_file in zip_files) {
  output_name <- gsub(".zip", ".h5ad", basename(zip_file))
  cat("Converting:", zip_file, "→", output_name, "\n")
  iCellbio2H5ad(zip_file, output_name)
}
```

### Batch Convert h5ad Files

```r
# Get all h5ad files
h5ad_files <- list.files(pattern = "*.h5ad", full.names = TRUE)

# Batch convert to SCE objects
sce_list <- lapply(h5ad_files, function(file) {
  cat("Reading:", file, "\n")
  h5ad_to_sce(file)
})

names(sce_list) <- gsub(".h5ad", "", basename(h5ad_files))
```

## 🛠️ Troubleshooting

### Python Environment Issues

```r
# Check Python configuration
reticulate::py_config()

# Check anndata availability
check_anndata_available()

# Reconfigure environment
configure_python_env(conda_env = "base", verbose = TRUE)
```

### Memory Issues

```r
# Check object size
object.size(sce)

# For large datasets, consider batch processing
# or use h5ad format directly in Python
```

### File Path Issues

```r
# Use absolute paths
file_path <- file.path(getwd(), "data.zip")
iCellbio2H5ad(file_path, "output.h5ad")

# Check if file exists
if (!file.exists("data.zip")) {
  stop("File does not exist: data.zip")
}
```

## 💡 Best Practices

### 1. Recommended Workflow

```r
# Recommended analysis workflow
library(ICellbioRpy)

# 1. Configure environment
configure_python_env()

# 2. Read data
data <- read1Cellbio("results.zip")

# 3. Choose format based on needs
if (use_bioconductor) {
  sce <- as.SingleCellExperiment(data)
  # Bioconductor analysis...
} else if (use_seurat) {
  seurat <- as.Seurat(data)
  # Seurat analysis...
} else if (use_python) {
  iCellbio2H5ad("results.zip", "analysis.h5ad")
  # Python/Scanpy analysis...
}
```

### 2. Performance Optimization

```r
# For large datasets
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB

# Use sparse matrices
library(Matrix)
# Package automatically maintains sparsity
```

### 3. Data Validation

```r
# Verify data integrity after conversion
original_dims <- dim(counts(sce))
converted_sce <- h5ad_to_sce("temp.h5ad")
new_dims <- dim(counts(converted_sce))

identical(original_dims, new_dims)
#> [1] TRUE
```

## 📚 Additional Resources

- **Detailed Documentation**: `vignette("introduction", package = "ICellbioRpy")`
- **Function Help**: `?iCellbio2H5ad`, `?h5ad_to_sce`
- **Installation Guide**: `anndata_installation_guide.md`
- **Complete Examples**: `README.md`

## 🆘 Getting Help

```r
# Check package information
packageVersion("ICellbioRpy")

# Check session information
sessionInfo()

# Please include the above information when reporting issues
```

---

**Tip**: If you're using this for the first time, we recommend reading `vignette("introduction")` for a complete usage guide.