# Example: Using the improved .1CB functions with automatic column detection
#
# This example demonstrates how the new column detection features work
# to make data conversion easier for users.

library(ICellbioRpy)

# Load your 1Cellbio data
# data <- read1Cellbio("path/to/your/results.zip")

# Option 1: Automatic column detection (easiest)
# The function will automatically detect:
# - Gene ID column (looks for 'id', 'gene_id', 'gene_symbol', etc.)
# - Cell ID column (looks for 'cell_id', 'barcode', etc.)
# seurat_obj <- as.Seurat.1CB(data)
# sce_obj <- as.SingleCellExperiment.1CB(data)

# Option 2: Show available column options first
# This helps you understand what columns are available
# show_column_options(data)

# Option 3: Manual specification (if auto-detection doesn't work)
# seurat_obj <- as.Seurat.1CB(data,
#                            rownames = "id",        # Gene identifier column
#                            colnames = "cell_id")   # Cell identifier column

# Option 4: Turn off auto-detection for more control
# seurat_obj <- as.Seurat.1CB(data,
#                            rownames = "specific_gene_col",
#                            colnames = "specific_cell_col",
#                            auto_detect = FALSE)

# The new functions will provide helpful error messages if columns are not found:
# "Cannot detect gene identifier column. Available options: id, gene_name, symbol.
#  Use: as.Seurat.1CB(data, rownames = 'id') or see show_column_options(data)"

# Additional helper functions available:
# - detect_gene_id_column(colnames) - Detect gene ID column from a vector of column names
# - detect_cell_id_column(colnames) - Detect cell ID column from a vector of column names
# - validate_column_names(object, rownames, colnames) - Validate specified column names