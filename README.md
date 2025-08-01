# 1CellbioRpy

An R package for reading 1Cellbio pipeline results directly from zip files and converting them to Seurat, SingleCellExperiment objects, or h5ad format.

## Overview

This package provides functions to directly read compressed 1Cellbio single-cell RNA-seq analysis results and convert them to commonly used Bioconductor/Seurat objects or h5ad format for downstream analysis.

## Installation

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install 1CellbioRpy from GitHub
devtools::install_github("your_username/1CellbioRpy")
```

## Python Environment Configuration

**Important**: Starting from version 1.1.0, all functions that require Python/anndata automatically configure the Python environment using your current environment!

### Automatic Configuration

All h5ad-related functions now automatically use your current Python environment:

```r
library(1CellbioRpy)

# These functions automatically configure Python environment using current environment
iCellbio2H5ad("data.zip", "output.h5ad")  # Auto-configures Python
sce <- h5ad_to_sce("data.h5ad")             # Auto-configures Python  
seurat_obj <- h5ad_to_seurat("data.h5ad")   # Auto-configures Python
seurat_to_h5ad(seurat_obj, "output.h5ad")   # Auto-configures Python
```

### Manual Configuration (When Needed)

If anndata is not available in your current environment, you can specify a different one:

```r
# Configure to use a specific conda environment
configure_python_env(conda_env = "your_env_name")

# Or configure to use a specific Python path  
configure_python_env(python_path = "/path/to/python")

# Verify anndata is available
check_anndata_available()
```

**Note**: If anndata is not found, the functions will provide helpful error messages with installation instructions.

For detailed troubleshooting, see `anndata_installation_guide.md`.

## Usage

### Method 1: Two-step conversion (via 1CellbioData object)

```r
library(1CellbioRpy)

# Read 1Cellbio results from zip file
data <- read1Cellbio("path/to/1cellbio_results.zip")

# Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(data)

# Convert to Seurat object
seurat <- as.Seurat(data)

# Convert to h5ad format for Scanpy
as.h5ad(data, "output.h5ad")
```

### Method 2: Direct conversion to h5ad (Memory Efficient + Sparse Matrix Preservation)

```r
library(1CellbioRpy)

# Convert directly from zip file to h5ad format
# This preserves sparse matrix format and is highly memory efficient
# Typically achieves 70-90% reduction in file size compared to dense storage
iCellbio2H5ad("path/to/1cellbio_results.zip", "output.h5ad")
```

### Method 3: Convert h5ad files to R objects (Reverse conversion)

```r
library(1CellbioRpy)

# Convert h5ad file to SingleCellExperiment
sce <- h5ad_to_sce("data.h5ad")

# Convert h5ad file to Seurat object
seurat_obj <- h5ad_to_seurat("data.h5ad")

# Specify which matrix to use as main expression data
# Use X matrix as counts instead of logcounts
sce <- h5ad_to_sce("data.h5ad", use_x_as = "counts")
seurat_obj <- h5ad_to_seurat("data.h5ad", use_x_as = "counts")
```

### Method 4: Convert Seurat objects to h5ad files

```r
library(1CellbioRpy)

# Convert Seurat object to h5ad file
seurat_to_h5ad(seurat_obj, "output.h5ad")

# Use counts instead of data layer
seurat_to_h5ad(seurat_obj, "output.h5ad", layer = "counts")

# Specify different default assay (e.g., SCT)
seurat_to_h5ad(seurat_obj, "output.h5ad", default_assay = "SCT")

# Exclude dimensionality reductions if not needed
seurat_to_h5ad(seurat_obj, "output.h5ad", include_reductions = FALSE)
```

## Functions

### Core Conversion Functions
- `read1Cellbio()` - Reads 1Cellbio results from a zip file and creates a 1CellbioData object
- `as.SingleCellExperiment()` - Converts a 1CellbioData object to a SingleCellExperiment object
- `as.Seurat()` - Converts a 1CellbioData object to a Seurat object
- `iCellbio2H5ad()` - **NEW**: Directly converts a zip file to h5ad format without creating intermediate objects. **Preserves sparse matrix format** for optimal memory efficiency and storage (typically 70-90% file size reduction)
- `h5ad_to_sce()` - **NEW**: Converts h5ad files to SingleCellExperiment objects with sparse matrix preservation
- `h5ad_to_seurat()` - **NEW**: Converts h5ad files to Seurat objects with sparse matrix preservation
- `seurat_to_h5ad()` - **NEW**: Converts Seurat objects to h5ad files with sparse matrix preservation and metadata retention

### Python Environment Configuration
- `configure_python_env()` - **NEW**: Configures Python environment and prevents automatic anndata installation
- `check_anndata_available()` - **NEW**: Checks if anndata is available in the current Python environment

## Details

The package works with the standard 1Cellbio results structure which includes:

- Gene expression matrices (counts and logcounts)
- Cell metadata
- Gene metadata
- Dimensionality reductions (PCA, tSNE, UMAP)

### Data Structure

The 1Cellbio results contain:

1. **Count Data**: Raw count matrix stored as HDF5 sparse matrix
2. **Logcounts Data**: Log-transformed normalized data stored as HDF5 dense matrix
3. **Cell Metadata**: Per-cell information including tissue, cluster assignments, quality control metrics
4. **Gene Metadata**: Per-gene information
5. **Dimensionality Reductions**: PCA, tSNE, and UMAP embeddings

### Supported Output Formats

#### SingleCellExperiment

The SingleCellExperiment object contains:

- Assays: `counts` (sparse matrix) and `logcounts` (dense matrix)
- ColData: Cell metadata
- RowData: Gene metadata
- ReducedDims: PCA, tSNE, and UMAP embeddings

#### Seurat

The Seurat object contains:

- Assay: RNA assay with `counts` and `data` slots
- MetaData: Cell metadata
- Reductions: pca, tsne, and umap

#### h5ad (for Scanpy)

The h5ad file contains:

- X: Raw count matrix (preserved as sparse matrix for memory efficiency)
- layers: Log-transformed data in the 'logcounts' layer (preserved as sparse matrix)
- obs: Cell metadata
- var: Gene metadata
- obsm: Dimensionality reductions (PCA, tSNE, UMAP)

**Note**: The `iCellbio2H5ad()` function automatically preserves sparse matrix format, resulting in significant storage savings (typically 70-90% reduction in file size compared to dense matrix storage).

## Example Workflows

### Standard Workflow (Two-step conversion)

```r
# Load the package
library(1CellbioRpy)

# Read the data
data <- read1Cellbio("path/to/1cellbio_results.zip")

# Convert to SingleCellExperiment for use with Bioconductor packages
sce <- as.SingleCellExperiment(data)

# Or convert to Seurat object for use with Seurat functions
seurat <- as.Seurat(data)

# Or convert to h5ad for use with Scanpy in Python
as.h5ad(data, "results.h5ad")

# Perform downstream analysis
# With SingleCellExperiment:
library(scater)
sce <- logNormCounts(sce)
plotPCA(sce, colour_by = "level1class")

# With Seurat:
library(Seurat)
DimPlot(seurat, reduction = "umap", group.by = "level1class")

# With Scanpy (in Python):
# import scanpy as sc
# adata = sc.read_h5ad("results.h5ad")
# sc.pl.umap(adata, color='level1class')
```

### Memory-Efficient Workflow (Direct conversion with Sparse Matrix Preservation)

```r
# Load the package
library(1CellbioRpy)

# For large datasets, convert directly to h5ad format
# This preserves sparse matrix format and avoids creating intermediate R objects
# Achieves significant memory savings and storage efficiency
iCellbio2H5ad("path/to/1cellbio_results.zip", "results.h5ad")

# Optional: specify custom temporary directory and keep temp files for debugging
iCellbio2H5ad("path/to/1cellbio_results.zip", "results.h5ad", 
            temp_dir = "/custom/temp/dir", cleanup = FALSE)

# The resulting h5ad file will be much smaller (70-90% reduction) compared to dense storage
# Example: A dataset that would be ~470 MB as dense matrices becomes ~99 MB as sparse

# Then use the h5ad file in Python with Scanpy:
# import scanpy as sc
# adata = sc.read_h5ad("results.h5ad")
# print(f"Matrix sparsity: {(1 - adata.X.nnz / (adata.shape[0] * adata.shape[1])) * 100:.1f}%")
# sc.pl.umap(adata, color='level1class')
```

## Technical Implementation

The package handles HDF5 file reading with the `hdf5r` package and properly manages:

- Sparse matrix conversion from 10x format
- Proper dimension handling for all data types
- Metadata parsing from HDF5 datasets
- Automatic handling of zip file extraction
- Conversion to h5ad format using the `anndata` package
- Memory-efficient direct conversion for large datasets

## Requirements

- R >= 4.0
- Bioconductor packages: SingleCellExperiment
- CRAN packages: Seurat, hdf5r, jsonlite, Matrix, anndata, reticulate (>= 1.39.0)

## License

MIT