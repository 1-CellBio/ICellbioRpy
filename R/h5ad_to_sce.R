#' Convert h5ad file to SingleCellExperiment object
#'
#' This function reads an h5ad file and converts it to a SingleCellExperiment object,
#' preserving sparse matrix format for memory efficiency.
#'
#' @param h5ad_file Path to the h5ad file
#' @param use_x_as Character string specifying which data to use as the main expression matrix.
#'   Options: "logcounts" (default) or "counts". If "logcounts", X matrix becomes logcounts
#'   and counts layer becomes counts. If "counts", X matrix becomes counts.
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#'
#' @return A SingleCellExperiment object with:
#'   \itemize{
#'     \item assays: counts and logcounts (preserved as sparse matrices)
#'     \item colData: cell metadata from obs
#'     \item rowData: gene metadata from var
#'     \item reducedDims: dimensionality reductions from obsm
#'   }
#'
#' @examples
#' \dontrun{
#' # Convert h5ad to SingleCellExperiment
#' sce <- h5ad_to_sce("data.h5ad")
#' 
#' # Use X matrix as counts instead of logcounts
#' sce <- h5ad_to_sce("data.h5ad", use_x_as = "counts")
#' }
#'
#' @export
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom Matrix t
h5ad_to_sce <- function(h5ad_file, use_x_as = "logcounts", verbose = TRUE) {
  
  # Check if required packages are available
  if (!requireNamespace("anndata", quietly = TRUE)) {
    stop("Package 'anndata' is required but not installed. Please install it with: install.packages('anndata')")
  }
  
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required but not installed. Please install it from Bioconductor.")
  }
  
  if (!requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("Package 'S4Vectors' is required but not installed. Please install it from Bioconductor.")
  }
  
  # Check if file exists
  if (!file.exists(h5ad_file)) {
    stop("File not found: ", h5ad_file)
  }
  
  if (verbose) cat("Reading h5ad file:", h5ad_file, "\n")
  
  # Configure Python environment using unified configuration
  if (verbose) cat("Configuring Python environment...\n")
  tryCatch({
    configure_python_env(verbose = FALSE)
  }, error = function(e) {
    warning("Failed to configure Python environment: ", as.character(e$message))
    warning("Please ensure anndata is installed in your Python environment or specify a conda environment:")
    warning("  configure_python_env(conda_env = \"your_env_name\")")
    stop("Python environment configuration failed")
  })
  
  # Read h5ad file using anndata package
  adata <- anndata::read_h5ad(h5ad_file)
  
  if (verbose) {
    cat("Dataset dimensions:", adata$n_obs, "cells x", adata$n_vars, "genes\n")
  }
  
  # Prepare assays
  assays_list <- list()
  
  if (use_x_as == "logcounts") {
    # X matrix as logcounts, counts layer as counts
    if (verbose) cat("Using X matrix as logcounts...\n")
    assays_list[["logcounts"]] <- as(Matrix::t(adata$X), 'CsparseMatrix')
    
    # Check if counts layer exists
    if ("counts" %in% names(adata$layers)) {
      if (verbose) cat("Adding counts layer as counts assay...\n")
      assays_list[["counts"]] <- as(Matrix::t(adata$layers[['counts']]), 'CsparseMatrix')
    } else {
      if (verbose) cat("Warning: No 'counts' layer found, only logcounts will be available\n")
    }
    
  } else if (use_x_as == "counts") {
    # X matrix as counts
    if (verbose) cat("Using X matrix as counts...\n")
    assays_list[["counts"]] <- as(Matrix::t(adata$X), 'CsparseMatrix')
    
    # Check if logcounts layer exists
    if ("logcounts" %in% names(adata$layers)) {
      if (verbose) cat("Adding logcounts layer as logcounts assay...\n")
      assays_list[["logcounts"]] <- as(Matrix::t(adata$layers[['logcounts']]), 'CsparseMatrix')
    }
    
  } else {
    stop("use_x_as must be either 'logcounts' or 'counts'")
  }
  
  # Prepare cell metadata (colData)
  if (verbose) cat("Processing cell metadata...\n")
  col_data <- S4Vectors::DataFrame(adata$obs)
  
  # Prepare gene metadata (rowData)
  if (verbose) cat("Processing gene metadata...\n")
  row_data <- S4Vectors::DataFrame(adata$var)
  
  # Prepare dimensionality reductions
  reduced_dims <- list()
  if (length(adata$obsm) > 0) {
    if (verbose) cat("Processing dimensionality reductions...\n")
    for (reduction_name in names(adata$obsm)) {
      reduced_dims[[reduction_name]] <- adata$obsm[[reduction_name]]
    }
  }
  
  # Create SingleCellExperiment object
  if (verbose) cat("Creating SingleCellExperiment object...\n")
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assays_list,
    colData = col_data,
    rowData = row_data,
    reducedDims = reduced_dims
  )
  
  if (verbose) {
    cat("✓ Conversion completed successfully!\n")
    cat("SingleCellExperiment object created with:\n")
    cat("  - Assays:", paste(names(assays_list), collapse = ", "), "\n")
    cat("  - Cell metadata columns:", ncol(col_data), "\n")
    cat("  - Gene metadata columns:", ncol(row_data), "\n")
    cat("  - Dimensionality reductions:", paste(names(reduced_dims), collapse = ", "), "\n")
  }
  
  return(sce)
}