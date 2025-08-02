#' Read and integrate 10X MTX format single-cell data
#'
#' This function reads 10X MTX format data (matrix.mtx, features/genes.tsv, barcodes.tsv)
#' from multiple samples defined in a CSV file, performs basic QC filtering,
#' and exports the integrated data as an h5ad file.
#'
#' @param csv_file Character string specifying the path to the CSV file containing sample information.
#'   The CSV file should have 4 columns: Sample_id, mtx_fns, features_fns, barcodes_fns
#' @param output_h5ad Character string specifying the output h5ad file path
#' @param min_counts_per_cell Numeric value specifying the minimum total counts per cell for filtering (default: 200)
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#'
#' @return Invisibly returns the output file path
#'
#' @details
#' This function performs the following steps:
#' \itemize{
#'   \item Reads the CSV file containing sample information
#'   \item For each sample, reads the MTX, features, and barcodes files (supports .gz compression)
#'   \item Creates a features x cells sparse matrix for each sample
#'   \item Performs basic QC filtering (removes cells with total counts < min_counts_per_cell)
#'   \item Renames cell IDs to "Sample_id_cellID" to avoid duplicates
#'   \item Integrates all samples into a single matrix
#'   \item Exports the integrated data as an h5ad file
#' }
#'
#' The CSV file format should be:
#' Sample_id,mtx_fns,features_fns,barcodes_fns
#' sample1,/path/to/matrix.mtx.gz,/path/to/features.tsv.gz,/path/to/barcodes.tsv.gz
#' sample2,/path/to/matrix.mtx,/path/to/genes.tsv,/path/to/barcodes.tsv
#'
#' @examples
#' \dontrun{
#' # Create a CSV file with sample information
#' sample_info <- data.frame(
#'   Sample_id = c("sample1", "sample2"),
#'   mtx_fns = c("/path/to/sample1/matrix.mtx.gz", "/path/to/sample2/matrix.mtx"),
#'   features_fns = c("/path/to/sample1/features.tsv.gz", "/path/to/sample2/genes.tsv"),
#'   barcodes_fns = c("/path/to/sample1/barcodes.tsv.gz", "/path/to/sample2/barcodes.tsv")
#' )
#' write.csv(sample_info, "samples.csv", row.names = FALSE)
#'
#' # Read and integrate the data
#' read_10x_mtx_to_h5ad("samples.csv", "integrated_data.h5ad")
#' }
#'
#' @export
read_10x_mtx_to_h5ad <- function(csv_file, 
                                  output_h5ad, 
                                  min_counts_per_cell = 200, 
                                  verbose = TRUE) {
  
  # Check if required packages are available
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed")
  }
  
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required but not installed")
  }
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required but not installed")
  }
  
  if (!requireNamespace("anndata", quietly = TRUE)) {
    stop("Package 'anndata' is required but not installed")
  }
  
  # Configure Python environment
  if (verbose) cat("Configuring Python environment...\n")
  tryCatch({
    configure_python_env(verbose = FALSE)
  }, error = function(e) {
    warning("Failed to configure Python environment: ", as.character(e$message))
    stop("Python environment configuration failed")
  })
  
  if (verbose) cat("Reading sample information from CSV file:", csv_file, "\n")
  
  # Read CSV file with sample information
  sample_info <- read_csv_sample_info(csv_file)
  
  if (verbose) {
    cat("Found", nrow(sample_info), "samples to process:\n")
    for (i in 1:nrow(sample_info)) {
      cat("  -", sample_info$Sample_id[i], "\n")
    }
  }
  
  # Initialize lists to store matrices and metadata
  all_matrices <- list()
  all_cell_metadata <- list()
  all_gene_names <- NULL
  
  # Process each sample
  for (i in 1:nrow(sample_info)) {
    sample_id <- sample_info$Sample_id[i]
    mtx_file <- sample_info$mtx_fns[i]
    features_file <- sample_info$features_fns[i]
    barcodes_file <- sample_info$barcodes_fns[i]
    
    if (verbose) cat("\nProcessing sample:", sample_id, "\n")
    
    # Read 10X data for this sample
    sample_data <- read_single_10x_sample(
      mtx_file = mtx_file,
      features_file = features_file,
      barcodes_file = barcodes_file,
      sample_id = sample_id,
      min_counts_per_cell = min_counts_per_cell,
      verbose = verbose
    )
    
    # Store the matrix and metadata
    all_matrices[[sample_id]] <- sample_data$matrix
    all_cell_metadata[[sample_id]] <- sample_data$cell_metadata
    
    # Check gene consistency
    current_gene_names <- rownames(sample_data$matrix)
    if (is.null(all_gene_names)) {
      all_gene_names <- current_gene_names
    } else {
      if (!identical(all_gene_names, current_gene_names)) {
        warning("Gene names differ between samples. Using intersection of genes.")
        all_gene_names <- intersect(all_gene_names, current_gene_names)
      }
    }
  }
  
  if (verbose) cat("\nIntegrating data from all samples...\n")
  
  # Integrate all samples
  integrated_data <- integrate_samples(
    matrices = all_matrices,
    cell_metadata_list = all_cell_metadata,
    gene_names = all_gene_names,
    verbose = verbose
  )
  
  if (verbose) cat("Creating h5ad file...\n")
  
  # Create h5ad file
  create_h5ad_from_integrated_data(
    matrix = integrated_data$matrix,
    cell_metadata = integrated_data$cell_metadata,
    gene_metadata = integrated_data$gene_metadata,
    output_file = output_h5ad,
    verbose = verbose
  )
  
  if (verbose) {
    cat("✓ Integration completed successfully!\n")
    cat("H5AD file created:", output_h5ad, "\n")
    file_size <- file.size(output_h5ad)
    cat("File size:", round(file_size / 1024^2, 2), "MB\n")
  }
  
  invisible(output_h5ad)
}


#' Read CSV file with sample information
#' 
#' @param csv_file Path to the CSV file
#' @return A data.frame with sample information
#' @keywords internal
read_csv_sample_info <- function(csv_file) {
  if (!file.exists(csv_file)) {
    stop("CSV file not found: ", csv_file)
  }
  
  # Read CSV file
  sample_info <- data.table::fread(csv_file, data.table = FALSE)
  
  # Check required columns
  required_cols <- c("Sample_id", "mtx_fns", "features_fns", "barcodes_fns")
  missing_cols <- setdiff(required_cols, colnames(sample_info))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in CSV file: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check if all files exist
  all_files <- c(sample_info$mtx_fns, sample_info$features_fns, sample_info$barcodes_fns)
  missing_files <- all_files[!file.exists(all_files)]
  if (length(missing_files) > 0) {
    stop("The following files do not exist:\n", paste(missing_files, collapse = "\n"))
  }
  
  return(sample_info)
}


#' Read a single 10X sample
#' 
#' @param mtx_file Path to the MTX file
#' @param features_file Path to the features/genes file
#' @param barcodes_file Path to the barcodes file
#' @param sample_id Sample identifier
#' @param min_counts_per_cell Minimum counts per cell for filtering
#' @param verbose Whether to print verbose messages
#' @return A list containing matrix and cell metadata
#' @keywords internal
read_single_10x_sample <- function(mtx_file, features_file, barcodes_file, 
                                   sample_id, min_counts_per_cell, verbose) {
  
  if (verbose) cat("  Reading MTX file:", basename(mtx_file), "\n")
  
  # Read MTX file using data.table for speed
  # MTX format: first line is header, second line is dimensions, then data
  mtx_lines <- readLines(mtx_file)
  
  # Skip comment lines (start with %)
  data_start <- which(!grepl("^%", mtx_lines))[1]
  
  # Read dimensions from second non-comment line
  dims_line <- mtx_lines[data_start]
  dims <- as.numeric(strsplit(dims_line, "\\s+")[[1]])
  n_genes <- dims[1]
  n_cells <- dims[2]
  n_entries <- dims[3]
  
  if (verbose) cat("    Matrix dimensions:", n_genes, "genes x", n_cells, "cells\n")
  
  # Read the actual data (skip header and dimension lines)
  if (n_entries > 0) {
    # Use data.table::fread to read the MTX data quickly
    temp_file <- tempfile()
    writeLines(mtx_lines[(data_start + 1):length(mtx_lines)], temp_file)
    
    mtx_data <- data.table::fread(temp_file, header = FALSE, col.names = c("gene", "cell", "count"))
    unlink(temp_file)
    
    # Convert to 1-based indexing (MTX is 1-based already, but just to be sure)
    i <- mtx_data$gene
    j <- mtx_data$cell
    x <- mtx_data$count
  } else {
    # Empty matrix
    i <- integer(0)
    j <- integer(0)
    x <- numeric(0)
  }
  
  if (verbose) cat("  Reading features file:", basename(features_file), "\n")
  
  # Read features/genes file
  features <- data.table::fread(features_file, header = FALSE, data.table = FALSE)
  if (ncol(features) >= 2) {
    gene_ids <- features[, 1]
    gene_symbols <- features[, 2]
  } else {
    gene_ids <- features[, 1]
    gene_symbols <- gene_ids
  }
  
  if (verbose) cat("  Reading barcodes file:", basename(barcodes_file), "\n")
  
  # Read barcodes file
  barcodes <- data.table::fread(barcodes_file, header = FALSE, data.table = FALSE)[, 1]
  
  # Create sparse matrix (genes x cells)
  if (verbose) cat("  Creating sparse matrix...\n")
  count_matrix <- Matrix::sparseMatrix(
    i = i,
    j = j,
    x = x,
    dims = c(n_genes, n_cells),
    dimnames = list(gene_symbols, barcodes)
  )
  
  # Perform QC filtering
  if (verbose) cat("  Performing QC filtering (min counts per cell:", min_counts_per_cell, ")...\n")
  
  # Calculate total counts per cell
  cell_counts <- Matrix::colSums(count_matrix)
  
  # Filter cells
  cells_to_keep <- cell_counts >= min_counts_per_cell
  n_cells_before <- ncol(count_matrix)
  n_cells_after <- sum(cells_to_keep)
  
  if (verbose) {
    cat("    Cells before filtering:", n_cells_before, "\n")
    cat("    Cells after filtering:", n_cells_after, "\n")
    cat("    Cells removed:", n_cells_before - n_cells_after, "\n")
  }
  
  # Apply filtering
  count_matrix <- count_matrix[, cells_to_keep, drop = FALSE]
  barcodes <- barcodes[cells_to_keep]
  cell_counts <- cell_counts[cells_to_keep]
  
  # Rename cells to include sample ID
  new_cell_names <- paste0(sample_id, "_", barcodes)
  colnames(count_matrix) <- new_cell_names
  
  # Create cell metadata
  cell_metadata <- data.frame(
    cell_id = new_cell_names,
    original_barcode = barcodes,
    sample_id = sample_id,
    total_counts = cell_counts,
    stringsAsFactors = FALSE,
    row.names = new_cell_names
  )
  
  return(list(
    matrix = count_matrix,
    cell_metadata = cell_metadata
  ))
}


#' Integrate matrices from multiple samples
#' 
#' @param matrices List of sparse matrices (one per sample)
#' @param cell_metadata_list List of cell metadata data.frames
#' @param gene_names Vector of gene names to use for integration
#' @param verbose Whether to print verbose messages
#' @return A list containing integrated matrix and metadata
#' @keywords internal
integrate_samples <- function(matrices, cell_metadata_list, gene_names, verbose) {
  
  # Subset all matrices to common genes
  for (i in seq_along(matrices)) {
    matrices[[i]] <- matrices[[i]][gene_names, , drop = FALSE]
  }
  
  if (verbose) {
    cat("  Common genes across all samples:", length(gene_names), "\n")
    total_cells <- sum(sapply(matrices, ncol))
    cat("  Total cells after filtering:", total_cells, "\n")
  }
  
  # Combine matrices horizontally (cbind)
  integrated_matrix <- do.call(cbind, matrices)
  
  # Combine cell metadata
  integrated_cell_metadata <- do.call(rbind, cell_metadata_list)
  
  # Create gene metadata
  gene_metadata <- data.frame(
    gene_id = gene_names,
    gene_symbol = gene_names,
    stringsAsFactors = FALSE,
    row.names = gene_names
  )
  
  return(list(
    matrix = integrated_matrix,
    cell_metadata = integrated_cell_metadata,
    gene_metadata = gene_metadata
  ))
}


#' Create h5ad file from integrated data
#' 
#' @param matrix Integrated sparse matrix (genes x cells)
#' @param cell_metadata Cell metadata data.frame
#' @param gene_metadata Gene metadata data.frame
#' @param output_file Output h5ad file path
#' @param verbose Whether to print verbose messages
#' @keywords internal
create_h5ad_from_integrated_data <- function(matrix, cell_metadata, gene_metadata, 
                                              output_file, verbose) {
  
  # Transpose matrix for h5ad format (cells x genes)
  matrix_t <- Matrix::t(matrix)
  
  if (verbose) {
    cat("  Final matrix dimensions:", nrow(matrix_t), "cells x", ncol(matrix_t), "genes\n")
    cat("  Matrix sparsity:", round((1 - Matrix::nnzero(matrix_t) / length(matrix_t)) * 100, 2), "%\n")
  }
  
  # Create AnnData object
  adata <- anndata::AnnData(
    X = matrix_t,
    obs = cell_metadata,
    var = gene_metadata
  )
  
  # Write to h5ad file
  adata$write_h5ad(output_file)
  
  invisible(NULL)
}