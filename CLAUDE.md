# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

ICellbioRpy is an R package for converting single-cell RNA sequencing data between formats, specifically:
- 1CellBio ZIP files → Seurat/SingleCellExperiment/H5AD
- Stereo-seq GEF files → Seurat/SingleCellExperiment/H5AD
- 10X MTX files → H5AD
- Seurat objects → H5AD
- H5AD → R objects

The package is designed for memory efficiency with sparse matrix preservation and includes spatial transcriptomics visualization capabilities.

## Development Environment

### Required Setup
1. Activate conda environment before running R: `source ~/anaconda3/etc/profile.d/conda.sh && conda activate 1cellbio`
2. R version: 4.3.3 (installed in conda env)
3. Python: 3.13 (in conda env, required for anndata/scanpy integration)

### Key Dependencies
- **Seurat** (>= 4.3.0) - Single-cell analysis in R
- **SingleCellExperiment** - Bioconductor S4 class for single-cell data
- **anndata** - Python package for annotated data (via reticulate)
- **hdf5r** - HDF5 file interface
- **Matrix** - Sparse matrix support
- **reticulate** - R-Python interface

## Common Development Commands

```bash
# Activate conda environment
source ~/anaconda3/etc/profile.d/conda.sh && conda activate 1cellbio

# Install package locally
R CMD INSTALL --no-multiarch --with-keep.source .

# Run tests
Rscript -e "devtools::test()"

# Check package (lightweight)
R CMD check --no-manual --no-build-vignettes .

# Generate documentation
Rscript -e "devtools::document()"

# Build source package
R CMD build .
```

## Architecture

### Core Data Structures
- **1CellbioData**: S3 class for 1CellBio zip file contents
- **StereoData**: S3 class for Stereo-seq GEF file contents
- Conversion functions use S3 methods for `as.Seurat()` and `as.SingleCellExperiment()`

### Key Modules

1. **Data Readers** (`R/read*.R`):
   - `read1Cellbio()` - Extracts and parses 1CellBio zip files
   - `read_gef()` - Reads Stereo-seq GEF files with cell border support
   - `read_10x_mtx()` - Handles 10X MTX format

2. **Format Converters** (`R/*_to_*.R`):
   - `iCellbio2H5ad()` - Direct zip to H5AD conversion
   - `seurat_to_h5ad()` - Seurat to H5AD with assay/reduction support
   - `gef_to_h5ad()` - GEF to H5AD with spatial data

3. **Python Integration** (`R/python_config.R`):
   - `configure_python_env()` - Auto-detects/configures Python environment
   - `check_anndata_available()` - Validates anndata installation
   - Uses reticulate for Python-R bridge

4. **Utilities** (`R/utils.R`):
   - `read_hdf5_sparse_matrix()` - HDF5 to dgCMatrix conversion
   - `read_hdf5_matrix()` - General HDF5 matrix reader
   - `construct_sparse_matrix()` - Builds dgCMatrix from 10X format

### Memory Efficiency
- Preserves sparse matrices throughout the pipeline
- Uses lazy evaluation where possible
- Direct file-to-file conversion to avoid intermediate objects

## Testing

Test structure follows testthat convention:
- `tests/testthat/test-*.R` files for each module
- Tests cover both success and error scenarios
- Include Python integration tests with anndata availability checks

## Known Issues

1. NAMESPACE warning: `as.Seurat.1CellbioData_impl` declared but not found
2. Generic function tests failing for `as.Seurat.1CB` and `as.SingleCellExperiment.1CB`
3. Description file uses `Authors@R` field instead of separate Author/Maintainer fields

## Spatial Transcriptomics Support

- Cell border information preserved in `uns` field of H5AD files
- `plot_cells_with_borders()` for spatial visualization
- Supports GEF file cell border extraction

## Python Environment Handling

The package provides intelligent Python environment configuration:
- Auto-detects conda environments
- Validates anndata availability
- Falls back to system Python if needed
- All H5AD operations require proper Python setup