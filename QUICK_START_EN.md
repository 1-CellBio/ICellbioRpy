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

ICellbioRpy provides a complete single-cell and spatial transcriptomics data format conversion ecosystem:

- **Read 1Cellbio Results** → `read1Cellbio()`
- **Read Stereo-seq GEF Files** → `read_gef()`
- **Convert to h5ad Format** → `iCellbio2H5ad()`, `gef_to_h5ad()`
- **h5ad to R Objects** → `h5ad_to_sce()`, `h5ad_to_seurat()`
- **R Objects to h5ad** → `seurat_to_h5ad()`
- **Spatial Transcriptomics Visualization** → `plot_cells_with_borders()`
- **GMT Gene Set Preprocessing** → `preprocess_gmt_custom()`
- **Python Environment Configuration** → `configure_python_env()`, `smart_python_config()`

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

### Smart Interactive Configuration (Recommended)

```r
# Automatically detect all conda environments and list those with anndata
smart_python_config(verbose = TRUE, interactive = TRUE)

# Output example:
# 📋 Found multiple environments with anndata:
#   1. 1cellbio (anndata 0.12.0)
#   2. atlas (anndata 0.11.3)
#   3. scanpy (anndata 0.10.9)
# Please select environment to use (1-3):
```

If only one environment has anndata, it will be auto-selected. If multiple exist, you'll be prompted to choose.

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
# Method 1: Using new .1CB functions (recommended, auto-detects column names)
sce <- as.SingleCellExperiment.1CB(data)

# Method 2: Manually specify column names
sce <- as.SingleCellExperiment(data,
                              rownames = "id",        # gene name column
                              colnames = "cell_id")   # cell name column

# Method 3: View available options before converting
show_column_options(data)
# This will display:
# === Column Detection Results ===
#
# Available gene identifier columns (rownames):
#   id, gene_symbol, gene_name
#   → Detected: id
#
# Available cell identifier columns (colnames):
#   cell_id, barcode, sample_id
#   → Detected: cell_id

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
# Method 1: Using new .1CB functions (recommended, auto-detects column names)
seurat <- as.Seurat.1CB(data)

# Method 2: Manually specify column names
seurat <- as.Seurat(data,
                   rownames = "id",        # gene name column
                   colnames = "cell_id")   # cell name column

# Method 3: Disable auto-detection and specify columns
seurat <- as.Seurat.1CB(data,
                        rownames = "specific_gene_col",
                        colnames = "specific_cell_col",
                        auto_detect = FALSE)

# View data
seurat
#> An object of class Seurat
#> 20000 features across 3000 samples

# Access data
counts_matrix <- GetAssayData(seurat, layer = "counts")
cell_metadata <- seurat[[]]
pca_coords <- Embeddings(seurat, reduction = "pca")
```

### 3. Direct Conversion to h5ad Format

```r
# One-step conversion: zip → h5ad
iCellbio2H5ad("1cellbio_results.zip", "output.h5ad")

# With parameter controls
iCellbio2H5ad(
  "1cellbio_results.zip",
  "output.h5ad",
  overwrite = FALSE,              # do not overwrite existing files
  name_conflict = "make_unique"   # auto-rename on conflicts ("make_unique" or "error")
)

# Check output file
file.exists("output.h5ad")
#> [1] TRUE

# Check file size
file.info("output.h5ad")$size / (1024^2)  # MB
```

**Parameter Notes:**
- `overwrite`: Whether to overwrite existing output file (default FALSE)
- `name_conflict`: Naming conflict handling strategy
  - `"make_unique"`: Automatically rename conflicting row/column names
  - `"error"`: Error on conflicts

## 🔄 Bidirectional h5ad Conversion

### From h5ad to R Objects

```r
# h5ad → SingleCellExperiment
sce <- h5ad_to_sce("data.h5ad")

# h5ad → Seurat
seurat <- h5ad_to_seurat("data.h5ad")

# With parameter controls
sce <- h5ad_to_sce(
  "data.h5ad",
  use_x_as = "auto",            # auto-detect X layer type ("auto"/"logcounts"/"counts")
  name_conflict = "make_unique" # naming conflict handling
)

# Check conversion results
assayNames(sce)
#> [1] "X" "raw"

names(seurat@reductions)
#> [1] "X_pca" "X_umap"
```

**Parameter Notes:**
- `use_x_as`: X matrix parsing method
  - `"auto"`: Auto-detect (default)
  - `"logcounts"`: Use as normalized data
  - `"counts"`: Use as raw counts
- `name_conflict`: Naming conflict handling strategy ("make_unique" or "error")

### From R Objects to h5ad

```r
# Seurat → h5ad (with overwrite and naming conflict control)
seurat_to_h5ad(
  seurat_object,
  "seurat_output.h5ad",
  overwrite = FALSE,              # do not overwrite by default
  name_conflict = "make_unique"   # or "error" to fail on conflicts
)

# Verify conversion
file.exists("seurat_output.h5ad")
#> [1] TRUE
```

## 🧬 Stereo-seq Spatial Transcriptomics Support

ICellbioRpy supports Stereo-seq GEF file format reading and conversion, including cell segmentation data and cell boundary information.

### Read GEF Files

```r
# Read GEF file (including cell boundaries)
stereo_data <- read_gef(
  "sample1.gef",
  bin_type = "cell_bins",      # Use cell segmentation data
  include_cellborder = TRUE    # Include cell boundary information
)

# View data summary
summary(stereo_data)
#> StereoData Object
#> Genes: 30000
#> Cells: 5000
#> Spatial coordinates available
#> Cell borders available

# View spatial coordinate range
head(stereo_data$spatial_coords)
```

### Convert to R Objects

```r
# Convert to Seurat object
seurat <- as.Seurat(stereo_data)

# Spatial coordinates stored in reductions
seurat@reductions$spatial
#> Coordinate system: spatial

# Cell boundaries stored in @misc
seurat@misc$cell_borders

# Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(stereo_data)

# Spatial coordinates in reducedDims
reducedDim(sce, "spatial")
```

### Spatial Visualization

```r
# Plot spatial cells with borders
plot_cells_with_borders(
  stereo_data,
  color_by = "cluster",        # Color by cluster
  show_borders = TRUE,         # Show cell boundaries
  border_color = "gray",       # Border color
  border_size = 0.5,           # Border line width
  point_size = 1               # Cell point size
)

# Color by gene expression
plot_cells_with_borders(
  stereo_data,
  color_by = "EPCAM",
  show_borders = TRUE
)
```

### Direct Conversion to H5AD

```r
# GEF → H5AD direct conversion (memory efficient)
gef_to_h5ad(
  "sample1.gef",
  "output.h5ad",
  bin_type = "cell_bins",
  include_cellborder = TRUE,    # Preserve cell boundaries
  include_spatial = TRUE,       # Preserve spatial coordinates
  overwrite = FALSE
)

# Or read first then convert
stereo_to_h5ad(
  stereo_data,
  "output.h5ad",
  layer = "counts"
)
```

### Advanced Options

```r
# Read specific spatial region
stereo_data <- read_gef(
  "sample.gef",
  region = c(1000, 3000, 1000, 3000),  # minX, maxX, minY, maxY
  max_cells = 10000,                    # Limit cell count
  gene_list = c("EPCAM", "KRT8", "VIM") # Only read specified genes
)

# Use square bins (without cell segmentation)
stereo_bins <- read_gef(
  "sample.gef",
  bin_type = "bins",
  bin_size = 50  # 50x50 pixel bins
)
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

## 🧫 GMT Gene Set Preprocessing

ICellbioRpy provides GMT (Gene Matrix Transposed) file preprocessing functionality with support for multiple gene ID type mappings.

### Basic Usage

```r
# Preprocess GMT file
preprocess_gmt_custom(
  gmt_file = "pathways.gmt",
  species = "9606",          # Species ID (9606 = human)
  output_dir = "gesel_output"
)
```

### Prepare Gene Mapping Files

Gene mapping files should be in TSV format with two columns: gene ID and gene symbol.

```r
# symbol.tsv format example
# gene_id	symbol
# 1	A1BG
# 2	A1BG-AS1

# entrez.tsv format example
# gene_id	entrez_id
# A1BG	1
# A1BG-AS1	503538

# ensembl.tsv format example
# gene_id	ensembl_id
# A1BG	ENSG00000121410
# A1BG-AS1	ENSG00000272398
```

### Supported Species

| Species ID | Species Name |
|------------|--------------|
| 9606 | Human (Homo sapiens) |
| 10090 | Mouse (Mus musculus) |
| 10116 | Rat (Rattus norvegicus) |
| 7227 | Fruit fly (Drosophila melanogaster) |
| 6239 | Nematode (Caenorhabditis elegans) |
| 7955 | Zebrafish (Danio rerio) |
| 9598 | Chimpanzee (Pan troglodytes) |

### Performance Optimization: Pre-build Lookup Tables

For large-scale GMT file processing, pre-build gene lookup tables for better performance:

```r
# Step 1: Pre-build lookup tables for all species
prebuild_gene_lookup_tables(
  data_dir = "~/gene_mapping_data",
  output_file = "gene_lookup_tables.rdata"
)

# Step 2: Load into global environment (creates master_lookup_tables variable)
load("gene_lookup_tables.rdata")

# Step 3: Call preprocessing function, will auto-detect and use pre-built tables
preprocess_gmt_custom("pathways.gmt", species = "9606")

# Batch process multiple GMT files (significant performance advantage)
for (gmt in list.files(pattern = "\\.gmt$")) {
  preprocess_gmt_custom(gmt, species = "9606")
}
```

**Note:** Pre-built lookup tables are passed implicitly via the global environment variable `master_lookup_tables`. After loading, `preprocess_gmt_custom()` will auto-detect and use them without requiring explicit parameters.

### Output Format

After preprocessing, the following files are generated in the output directory:

```
gesel_output/
├── collections.tsv       # Collection metadata
├── sets.tsv              # Gene set definitions
├── set2gene.tsv          # Set to gene mapping (delta-encoded)
└── gene2gene.tsv         # Gene ID mapping (delta-encoded)
```

### Advanced Usage

```r
# Use with custom mapping files
preprocess_gmt_with_custom_mapping(
  gmt_file = "custom_pathways.gmt",
  species = "10090",  # Mouse
  gene_mapping_files = list(
    symbol = "mouse_symbols.tsv",
    entrez = "mouse_entrez.tsv",
    ensembl = "mouse_ensembl.tsv"
  ),
  collection_name = "mouse_pathways",
  collection_desc = "Custom mouse pathways",
  output_dir = "mouse_gesel",
  auto_download_missing = TRUE  # Auto-download missing mapping files
)
```

## 🛠️ Troubleshooting

### Column Detection Issues

```r
# If auto-detection fails, first view available options
show_column_options(data)

# Common errors and solutions:
# Error: "Cannot detect gene identifier column"
# Solution: Manually specify gene column
sce <- as.SingleCellExperiment.1CB(data, rownames = "gene_name", colnames = "cell_id")

# Error: "Column 'xxx' not found"
# Solution: View correct column names
show_column_options(data)
# Or check manually
# coldata <- data$experiment$summarized_experiment$column_data$resource
# rowdata <- data$experiment$summarized_experiment$row_data$resource

# Disable auto-detection for full manual control
sce <- as.SingleCellExperiment.1CB(data,
                                    rownames = "your_gene_col",
                                    colnames = "your_cell_col",
                                    auto_detect = FALSE)
```

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
  # New version: auto-detect column names (recommended)
  sce <- as.SingleCellExperiment.1CB(data)

  # Old version: manually specify column names
  # sce <- as.SingleCellExperiment(data, rownames = "id", colnames = "cell_id")

  # Bioconductor analysis...
} else if (use_seurat) {
  # New version: auto-detect column names (recommended)
  seurat <- as.Seurat.1CB(data)

  # Old version: manually specify column names
  # seurat <- as.Seurat(data, rownames = "id", colnames = "cell_id")

  # Seurat analysis...
} else if (use_python) {
  iCellbio2H5ad("results.zip", "analysis.h5ad")
  # Python/Scanpy analysis...
}

# 4. View available column options (if auto-detection fails)
# show_column_options(data)
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