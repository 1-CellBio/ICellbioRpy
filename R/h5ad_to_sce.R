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
h5ad_to_sce <- function(h5ad_file, use_x_as = "auto", name_conflict = c("make_unique", "error"), verbose = TRUE) {
  name_conflict <- match.arg(name_conflict)
  
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
  
  if (verbose) cat(icb_i18n("读取 h5ad 文件:", "Reading h5ad file:"), h5ad_file, "\n")
  
  # Configure Python environment (REQUIRED)
  tryCatch({
    smart_python_config(verbose = verbose, interactive = FALSE)
  }, error = function(e) {
    stop(icb_i18n(
      zh = paste0("Python环境配置失败: ", e$message, "\n",
                 "请手动配置或安装所需环境:\n",
                 "  smart_python_config(verbose = TRUE, interactive = TRUE)"),
      en = paste0("Python environment configuration failed: ", e$message, "\n",
                 "Please configure manually or install required environment:\n",
                 "  smart_python_config(verbose = TRUE, interactive = TRUE)")
    ))
  })
  
  # Read h5ad file using anndata package
  adata <- anndata::read_h5ad(h5ad_file)
  
  if (verbose) {
    cat(icb_i18n("数据维度:", "Dataset dimensions:"), adata$n_obs, icb_i18n("细胞 x", "cells x"), adata$n_vars, icb_i18n("基因\n", "genes\n"))
  }
  
  # Check if use_x_as is set to auto, then determine based on data characteristics
  if (use_x_as == "auto") {
    # Check if X matrix values are integers (typical for counts)
    X_data <- adata$X
    if (inherits(X_data, "sparseMatrix") || inherits(X_data, "dgCMatrix")) {
      # For sparse matrices, randomly sample values to determine if they are integers
      x_vals <- X_data@x
      if (length(x_vals) > 0) {
        # Randomly sample up to 1000 values to determine if they are integers
        sample_size <- min(1000, length(x_vals))
        if (sample_size < length(x_vals)) {
          sampled_indices <- sample(length(x_vals), sample_size)
          sampled_vals <- x_vals[sampled_indices]
        } else {
          sampled_vals <- x_vals
        }
        
        are_integers <- all(sampled_vals == floor(sampled_vals))
        if (are_integers) {
          use_x_as <- "counts"
          if (verbose) cat("Detected X matrix as counts (integer values)\n")
        } else {
          use_x_as <- "logcounts"
          if (verbose) cat("Detected X matrix as logcounts (non-integer values)\n")
        }
      } else {
        # All values are zero, default to counts
        use_x_as <- "counts"
        if (verbose) cat("Detected X matrix as counts (all zero values)\n")
      }
    } else {
      # For dense matrices, check a sample of values
      x_vals <- X_data[1:min(1000, length(X_data))]  # Check first 1000 values or all if fewer
      x_vals <- x_vals[x_vals != 0]  # Remove zeros
      if (length(x_vals) > 0) {
        are_integers <- all(x_vals == floor(x_vals))
        if (are_integers) {
          use_x_as <- "counts"
          if (verbose) cat("Detected X matrix as counts (integer values)\n")
        } else {
          use_x_as <- "logcounts"
          if (verbose) cat("Detected X matrix as logcounts (non-integer values)\n")
        }
      } else {
        # All values are zero, default to counts
        use_x_as <- "counts"
        if (verbose) cat("Detected X matrix as counts (all zero values)\n")
      }
    }
  }
  
  # Prepare assays
  assays_list <- list()
  
  if (use_x_as == "logcounts") {
    # X matrix as logcounts, counts layer as counts
    if (verbose) cat(icb_i18n("使用 X 作为 logcounts...\n", "Using X matrix as logcounts...\n"))
    assays_list[["logcounts"]] <- icb_anndata_to_csparse(adata$X, transpose = TRUE, verbose = verbose)
    
    # Check if counts layer exists
    if ("counts" %in% names(adata$layers)) {
      if (verbose) cat(icb_i18n("加入 counts 层作为 counts assay...\n", "Adding counts layer as counts assay...\n"))
      assays_list[["counts"]] <- icb_anndata_to_csparse(adata$layers[['counts']], transpose = TRUE, verbose = verbose)
    } else {
      if (verbose) cat(icb_i18n("警告: 未找到 'counts' 层，仅提供 logcounts\n", "Warning: No 'counts' layer found, only logcounts will be available\n"))
    }
    
  } else if (use_x_as == "counts") {
    # X matrix as counts
    if (verbose) cat(icb_i18n("使用 X 作为 counts...\n", "Using X matrix as counts...\n"))
    assays_list[["counts"]] <- icb_anndata_to_csparse(adata$X, transpose = TRUE, verbose = verbose)
    
    # Check if logcounts layer exists
    if ("logcounts" %in% names(adata$layers)) {
      if (verbose) cat(icb_i18n("加入 logcounts 层作为 logcounts assay...\n", "Adding logcounts layer as logcounts assay...\n"))
      assays_list[["logcounts"]] <- icb_anndata_to_csparse(adata$layers[['logcounts']], transpose = TRUE, verbose = verbose)
    }
    
  } else {
    stop("use_x_as must be either 'logcounts', 'counts', or 'auto'")
  }
  
  # Prepare cell metadata (colData)
  if (verbose) cat("Processing cell metadata...\n")
  col_data <- S4Vectors::DataFrame(adata$obs)
  if (!is.null(rownames(col_data))) {
    rownames(col_data) <- icb_make_unique(rownames(col_data), strategy = name_conflict, sep = "-")
  }
  
  # Prepare gene metadata (rowData)
  if (verbose) cat("Processing gene metadata...\n")
  row_data <- S4Vectors::DataFrame(adata$var)
  if (!is.null(rownames(row_data))) {
    rownames(row_data) <- icb_make_unique(rownames(row_data), strategy = name_conflict, sep = "-")
  }
  
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