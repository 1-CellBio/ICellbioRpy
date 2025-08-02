# Example: Reading and integrating 10X MTX format data
# This example demonstrates how to use read_10x_mtx_to_h5ad function

library(ICellbioRpy)

# Step 1: Create a CSV file with sample information
# The CSV should have 4 columns: Sample_id, mtx_fns, features_fns, barcodes_fns

sample_info <- data.frame(
  Sample_id = c("sample1", "sample2", "sample3"),
  mtx_fns = c(
    "/path/to/sample1/filtered_feature_bc_matrix/matrix.mtx.gz",
    "/path/to/sample2/filtered_feature_bc_matrix/matrix.mtx.gz",
    "/path/to/sample3/filtered_feature_bc_matrix/matrix.mtx.gz"
  ),
  features_fns = c(
    "/path/to/sample1/filtered_feature_bc_matrix/features.tsv.gz",
    "/path/to/sample2/filtered_feature_bc_matrix/features.tsv.gz",
    "/path/to/sample3/filtered_feature_bc_matrix/features.tsv.gz"
  ),
  barcodes_fns = c(
    "/path/to/sample1/filtered_feature_bc_matrix/barcodes.tsv.gz",
    "/path/to/sample2/filtered_feature_bc_matrix/barcodes.tsv.gz",
    "/path/to/sample3/filtered_feature_bc_matrix/barcodes.tsv.gz"
  ),
  stringsAsFactors = FALSE
)

# Save the sample information to a CSV file
write.csv(sample_info, "samples.csv", row.names = FALSE)

# Step 2: Read and integrate the 10X data
# This will:
# - Read each sample's MTX, features, and barcodes files
# - Filter cells with total counts < 200 (default)
# - Rename cells to avoid duplicates (sample_id + "_" + barcode)
# - Integrate all samples into a single matrix
# - Export as h5ad file

read_10x_mtx_to_h5ad(
  csv_file = "samples.csv",
  output_h5ad = "integrated_data.h5ad",
  min_counts_per_cell = 200,  # Remove cells with < 200 total counts
  verbose = TRUE
)

# The resulting h5ad file can be loaded in Python/scanpy for downstream analysis:
# import scanpy as sc
# adata = sc.read_h5ad("integrated_data.h5ad")

# Alternative: Use different QC thresholds
read_10x_mtx_to_h5ad(
  csv_file = "samples.csv",
  output_h5ad = "integrated_data_strict.h5ad",
  min_counts_per_cell = 500,  # More stringent filtering
  verbose = TRUE
)

# Example with mixed file formats (some .gz, some not)
sample_info_mixed <- data.frame(
  Sample_id = c("sample1", "sample2"),
  mtx_fns = c(
    "/path/to/sample1/matrix.mtx.gz",      # compressed
    "/path/to/sample2/matrix.mtx"          # uncompressed
  ),
  features_fns = c(
    "/path/to/sample1/features.tsv.gz",    # compressed
    "/path/to/sample2/genes.tsv"           # uncompressed, different name
  ),
  barcodes_fns = c(
    "/path/to/sample1/barcodes.tsv.gz",    # compressed
    "/path/to/sample2/barcodes.tsv"        # uncompressed
  ),
  stringsAsFactors = FALSE
)

write.csv(sample_info_mixed, "samples_mixed.csv", row.names = FALSE)

read_10x_mtx_to_h5ad(
  csv_file = "samples_mixed.csv",
  output_h5ad = "integrated_mixed.h5ad"
)