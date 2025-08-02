# ICellbioRpy 📊🧬

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/1-Cellbio/ICellbioRpy/workflows/R-CMD-check/badge.svg)](https://github.com/1-Cellbio/ICellbioRpy/actions)

> **One-Stop Single-Cell Data Format Conversion Solution**  
> Seamless conversion between 1CellBio, Seurat, SingleCellExperiment, 10X MTX and H5AD formats

## 🎯 Overview

ICellbioRpy is a professional R package for single-cell RNA-seq data format conversion. It provides a complete toolkit for seamless conversion from 1CellBio analysis results to mainstream single-cell analysis formats, with optimized user experience especially for beginners.

### ✨ Key Features

- 🚀 **One-Click Conversion**: 1CellBio ZIP → H5AD/Seurat/SingleCellExperiment
- 📊 **Multi-Sample Integration**: Integrate multiple 10X MTX files into single H5AD
- 💾 **Memory Efficient**: Automatically preserves sparse matrix format, saving 70-90% storage
- 🔄 **Bidirectional Conversion**: R objects ↔ Python formats
- 🎨 **Beginner-Friendly**: Super detailed environment setup tutorials
- 🐍 **Smart Python Configuration**: Automatic Python environment detection and configuration

---

## 🚀 Quick Start

### 📋 Environment Setup (Must Read for Beginners)

**If you are a complete beginner, please complete the environment setup first:**

| OS | Installation Tutorial | Description |
|----|--------------------|-------------|
| **Windows** | [📘 Windows Setup Guide](tutorials/tutorial_0_environment_setup_windows.html) | Complete guide from R to conda environment |
| **macOS** | [📘 macOS Setup Guide](tutorials/tutorial_0_environment_setup_macos.html) | Optimized for Apple Silicon |

### ⚡ Environment Check

```r
# Quick environment check
source("tutorials/check_environment.R")
```

### 📦 Installation

```r
# Install dependencies
install.packages(c("devtools", "reticulate", "Matrix", "data.table"))

# Install ICellbioRpy from GitHub
devtools::install_github("1-Cellbio/ICellbioRpy")
```

---

## 📚 Complete Tutorials

We provide detailed beginner-friendly tutorials that can be run directly in RStudio:

### 🔧 Environment Setup Tutorials
- [📘 Windows Setup](tutorials/tutorial_0_environment_setup_windows.html) - Complete Windows configuration guide
- [📘 macOS Setup](tutorials/tutorial_0_environment_setup_macos.html) - Complete macOS configuration guide

### 📊 Data Conversion Tutorials

| Tutorial | Input | Output | Difficulty | Duration |
|----------|-------|--------|------------|----------|
| [Tutorial 1](tutorials/tutorial_1_1cellbio_to_h5ad.html) | 1CellBio ZIP | H5AD | ⭐⭐ | 30-45min |
| [Tutorial 2](tutorials/tutorial_2_1cellbio_to_seurat.html) | 1CellBio ZIP | Seurat | ⭐⭐ | 45-60min |
| [Tutorial 3](tutorials/tutorial_3_seurat_to_h5ad.html) | Seurat Object | H5AD | ⭐⭐⭐ | 40-50min |
| [Tutorial 4](tutorials/tutorial_4_10x_mtx_to_h5ad.html) | 10X MTX Files | H5AD | ⭐⭐⭐⭐ | 50-70min |

### 📖 Additional Resources
- [📍 Start Here](tutorials/START_HERE.md) - Begin your learning journey here
- [⚡ Quick Guide](tutorials/QUICK_GUIDE.md) - Code snippets and quick reference
- [📚 Tutorial Overview](tutorials/README.md) - Complete tutorial list and descriptions

---

## 💻 Core Functions

### 🔄 Python Environment Configuration

**Important**: All H5AD-related functions require Python environment support. ICellbioRpy provides smart configuration:

```r
library(ICellbioRpy)

# Method 1: Auto-detect and use current environment (recommended)
configure_python_env(verbose = TRUE)

# Method 2: Specify conda environment
configure_python_env(conda_env = "1cellbio", verbose = TRUE)

# Method 3: Specify Python path
configure_python_env(python_path = "/path/to/python", verbose = TRUE)

# Verify configuration
check_anndata_available()
```

### 📈 Main Use Cases

#### 1. 1CellBio → H5AD (for Python analysis)

```r
# Direct conversion (recommended)
iCellbio2H5ad("path/to/1cellbio_results.zip", "output.h5ad")

# Or two-step conversion
data <- read1Cellbio("path/to/1cellbio_results.zip")
as.h5ad(data, "output.h5ad")
```

#### 2. 1CellBio → Seurat (for R analysis)

```r
# Read and convert
data <- read1Cellbio("path/to/1cellbio_results.zip")
seurat_obj <- as.Seurat.1CB(data)

# Visualization
library(Seurat)
DimPlot(seurat_obj, reduction = "umap", group.by = "level1class")
```

#### 3. Seurat → H5AD (cross-language conversion)

```r
# Convert R Seurat object to Python-compatible H5AD format
seurat_to_h5ad(seurat_obj, "output.h5ad")

# Advanced options
seurat_to_h5ad(seurat_obj, "output.h5ad",
               default_assay = "RNA",
               layer = "data",
               include_reductions = TRUE)
```

#### 4. 10X MTX → H5AD (multi-sample integration)

```r
# Prepare sample information CSV file
sample_info <- data.frame(
  Sample_id = c("sample1", "sample2"),
  mtx_fns = c("path/to/sample1/matrix.mtx.gz", "path/to/sample2/matrix.mtx.gz"),
  features_fns = c("path/to/sample1/features.tsv.gz", "path/to/sample2/features.tsv.gz"),
  barcodes_fns = c("path/to/sample1/barcodes.tsv.gz", "path/to/sample2/barcodes.tsv.gz")
)
write.csv(sample_info, "samples.csv", row.names = FALSE)

# Read and integrate
read_10x_mtx_to_h5ad(
  csv_file = "samples.csv",
  output_h5ad = "integrated.h5ad",
  min_counts_per_cell = 200
)
```

#### 5. H5AD → R Objects (reverse conversion)

```r
# H5AD to SingleCellExperiment
sce <- h5ad_to_sce("data.h5ad")

# H5AD to Seurat
seurat_obj <- h5ad_to_seurat("data.h5ad")

# Specify data layer
sce <- h5ad_to_sce("data.h5ad", use_x_as = "counts")
```

---

## 🛠️ Main Functions

### Core Conversion Functions
- `read1Cellbio()` - Read 1CellBio results from ZIP file
- `iCellbio2H5ad()` - **Direct conversion**: ZIP → H5AD (memory efficient)
- `as.Seurat.1CB()` - 1CellbioData → Seurat object
- `as.SingleCellExperiment.1CB()` - 1CellbioData → SingleCellExperiment object
- `seurat_to_h5ad()` - Seurat object → H5AD file
- `read_10x_mtx_to_h5ad()` - **New feature**: Multi-sample 10X data integration
- `h5ad_to_sce()` - H5AD → SingleCellExperiment
- `h5ad_to_seurat()` - H5AD → Seurat object

### Environment Configuration Functions
- `configure_python_env()` - Smart Python environment configuration
- `check_anndata_available()` - Check anndata availability

---

## 🔬 Supported Data Structures

### Input Formats
- **1CellBio ZIP files**: Contains HDF5-format expression matrices, metadata, and dimensionality reductions
- **Seurat objects**: Standard Seurat objects with multiple assays and reductions
- **10X MTX files**: Cell Ranger output in sparse matrix format
- **H5AD files**: Python scanpy format single-cell data

### Output Formats
- **H5AD**: Preserves sparse matrices, supports Python/scanpy analysis
- **Seurat**: Contains counts, data, and dimensionality reductions
- **SingleCellExperiment**: Bioconductor standard format
- **Integrated data**: Unified format for multi-sample merged data

---

## 💡 Performance Advantages

### Memory Efficiency
- **Sparse Matrix Preservation**: Automatically maintains sparse format, saving 70-90% storage space
- **Direct Conversion**: Avoids intermediate objects, reducing memory usage
- **Streaming Processing**: Supports efficient processing of large datasets

### User Experience
- **Smart Detection**: Automatic Python environment detection and configuration
- **Comprehensive Documentation**: Complete help documentation and examples for every function
- **Error Handling**: Friendly error messages and solution suggestions
- **Progress Display**: Detailed progress for long-running operations

---

## 🎯 Typical Workflows

### Scenario 1: Python users receiving 1CellBio data

```r
# Direct conversion to H5AD format
iCellbio2H5ad("1cellbio_results.zip", "analysis_data.h5ad")
```

```python
# Continue analysis in Python
import scanpy as sc
adata = sc.read_h5ad("analysis_data.h5ad")
sc.pl.umap(adata, color='level1class')
```

### Scenario 2: R users performing single-cell analysis

```r
# Read data
data <- read1Cellbio("1cellbio_results.zip")

# Convert to Seurat object
seurat_obj <- as.Seurat.1CB(data)

# Seurat analysis workflow
library(Seurat)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindClusters(seurat_obj)
DimPlot(seurat_obj, reduction = "umap")
```

### Scenario 3: Multi-sample 10X data integration

```r
# Create sample information file
# Then integrate multiple samples
read_10x_mtx_to_h5ad("samples.csv", "integrated.h5ad")
```

```python
# Downstream analysis in Python
import scanpy as sc
adata = sc.read_h5ad("integrated.h5ad")

# Batch correction and integration analysis
sc.pp.combat(adata, key='sample_id')
sc.tl.leiden(adata)
sc.pl.umap(adata, color=['sample_id', 'leiden'])
```

---

## 🆘 Technical Support

### Self-Help Solutions
1. **Run Environment Check**: `source("tutorials/check_environment.R")`
2. **Check Detailed Tutorials**: Choose tutorials based on your data type
3. **View Function Help**: `?function_name`

### Get Help
- 📧 **Submit Issue**: [GitHub Issues](https://github.com/1-Cellbio/ICellbioRpy/issues)
- 📚 **View Tutorials**: [Complete Tutorial Collection](tutorials/)
- 💬 **User Community**: Join user discussion groups

---

## 📋 System Requirements

### R Environment
- R ≥ 4.0.0
- Recommended packages: Seurat, SingleCellExperiment, Matrix, data.table, hdf5r, jsonlite

### Python Environment  
- Python ≥ 3.7
- Required packages: anndata ≥ 0.7.0
- Recommended packages: scanpy, pandas, numpy

### Hardware Recommendations
- **Memory**: 8GB+ (more needed for large datasets)
- **Storage**: Ensure sufficient disk space
- **CPU**: Multi-core CPU improves processing speed

---

## 📄 License

This project is licensed under the [MIT License](LICENSE)

---

**Start your single-cell data conversion journey!** 🚀

[![Get Started](https://img.shields.io/badge/Get_Started-Tutorial_Navigation-blue)](tutorials/START_HERE.md)
[![Quick Guide](https://img.shields.io/badge/Quick_Guide-Code_Examples-green)](tutorials/QUICK_GUIDE.md)
[![Environment Setup](https://img.shields.io/badge/Environment_Setup-Windows|macOS-orange)](tutorials/)