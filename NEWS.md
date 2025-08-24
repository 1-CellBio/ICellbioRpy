# ICellbioRpy 0.2.0

## Major Improvements

### API Consistency and Cleanup
* **BREAKING**: Deprecated `as.Seurat.1CellbioData()` and `as.SingleCellExperiment.1CellbioData()` 
* **NEW**: Unified API with `as.Seurat.1CB()` and `as.SingleCellExperiment.1CB()`
* Cleaned up NAMESPACE exports to reduce confusion
* Added deprecation warnings for old functions

### Enhanced Python/AnnData Integration
* **NEW**: Robust scipy.sparse matrix detection and conversion (`icb_is_py_sparse()`, `icb_py_sparse_to_dgc()`)
* **NEW**: Unified AnnData matrix conversion interface (`icb_anndata_to_csparse()`)
* Improved memory efficiency for large sparse matrices
* Better error handling for Python environment issues

### Multi-language Support
* **NEW**: Complete internationalization framework (`icb_i18n()`, `icb_get_lang()`)
* Support for Chinese and English interfaces
* Language control via `options(ICellbioRpy.lang = "zh"|"en")`
* Multilingual error messages and documentation

### Name Conflict Resolution
* **NEW**: Standardized name uniqueness handling (`icb_make_unique()`)
* Added `name_conflict` parameter to all conversion functions
* Options: `"make_unique"` (default) or `"error"`
* Consistent duplicate name handling across all functions

### File Overwrite Control
* **NEW**: `overwrite` parameter in `seurat_to_h5ad()` and related functions
* Default behavior: prevent accidental overwrites
* Explicit control over file replacement

### Testing and Quality Assurance
* **NEW**: Comprehensive test suite covering core functionality
* Unit tests for utility functions, conversions, and error handling
* Mock data generation for reproducible testing
* Improved code coverage and reliability

### Documentation and User Experience
* **NEW**: pkgdown website configuration with modern Bootstrap 5 theme
* **NEW**: Comprehensive vignettes for quick start and advanced usage
* **NEW**: Structured function reference with logical groupings
* Improved error messages with clear guidance
* Enhanced progress reporting and diagnostic information

## Technical Improvements

### Memory and Performance
* Optimized sparse matrix handling throughout the conversion pipeline
* Reduced memory footprint for large datasets
* Improved performance for multi-sample 10X integration
* Better handling of empty or zero matrices

### HDF5 Function Consolidation
* Unified HDF5 utility functions in `R/utils.R`
* Eliminated duplicate function definitions
* Improved maintainability and consistency
* Better error handling for malformed HDF5 files

### Dependency Management
* Updated DESCRIPTION with comprehensive dependency specifications
* Clear minimum version requirements for all dependencies
* Proper Bioconductor package version constraints
* SystemRequirements specification for Python dependencies

## Bug Fixes
* Fixed edge cases in matrix transposition during format conversion
* Improved handling of missing or empty data layers
* Better validation of input file formats
* Enhanced error recovery and user guidance

## Backwards Compatibility
* Deprecated functions remain available with warning messages
* Existing user code continues to work with deprecation notices
* Clear migration path to new API documented
* Gradual transition period for API changes

---

# ICellbioRpy 0.1.0

## Initial Release

### Core Features
* Read and convert 1CellBio ZIP results to R objects
* Convert between Seurat, SingleCellExperiment, and H5AD formats
* Basic 10X MTX file integration
* Python environment configuration helpers

### Basic Functionality
* `read1Cellbio()` for reading 1CellBio ZIP files
* `as.Seurat.1CellbioData()` and `as.SingleCellExperiment.1CellbioData()` for format conversion
* `seurat_to_h5ad()` for exporting to Python-compatible format
* `configure_python_env()` for Python environment setup

### Foundation
* Established basic package structure
* Initial documentation and examples
* Basic error handling and validation
* Support for sparse matrix preservation
