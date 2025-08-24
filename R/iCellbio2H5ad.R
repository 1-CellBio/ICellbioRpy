#' Convert 1Cellbio results zip file directly to h5ad format
#'
#' This function reads a 1Cellbio pipeline results zip file and directly converts it
#' to h5ad format without creating intermediate objects. This is more memory efficient
#' for large datasets and preserves sparse matrix format for optimal storage.
#'
#' The function automatically preserves the sparse matrix format of counts and logcounts
#' data, resulting in significant memory savings (typically 70-90% reduction in file size
#' compared to dense matrix storage).
#'
#' @param zip_path Path to the zip file containing 1Cellbio results
#' @param h5ad_path Path where the h5ad file should be saved
#' @param temp_dir Optional temporary directory for extraction. If NULL, uses system temp directory
#' @param cleanup Whether to clean up temporary files after conversion (default: TRUE)
#' @return NULL (file is saved to disk)
#' @export
#' @import anndata
#' @import reticulate
#' @import jsonlite
#' @import hdf5r
#' @import Matrix
#' @examples
#' \dontrun{
#' # Basic usage - preserves sparse matrix format automatically
#' iCellbio2H5ad("path/to/1Cellbio_results.zip", "output.h5ad")
#' 
#' # With custom temporary directory
#' iCellbio2H5ad("path/to/1Cellbio_results.zip", "output.h5ad", 
#'               temp_dir = "/custom/temp/dir")
#' }
iCellbio2H5ad <- function(zip_path, h5ad_path, temp_dir = NULL, cleanup = TRUE,
                          overwrite = FALSE, 
                          name_conflict = c("make_unique", "error"),
                          verbose = TRUE) {
  # Check if required packages are available
  required_packages <- c("anndata", "reticulate", "jsonlite", "hdf5r", "Matrix")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("The '", pkg, "' package is required for this function but is not installed.")
    }
  }
  
  name_conflict <- match.arg(name_conflict)
  # Check if zip file exists
  if (!file.exists(zip_path)) {
    stop("Zip file not found: ", zip_path)
  }
  # Overwrite behavior
  if (file.exists(h5ad_path) && !isTRUE(overwrite)) {
    stop(icb_i18n(
      zh = paste0("输出目标已存在：", h5ad_path, "。设置 overwrite=TRUE 以允许覆盖，或更换文件名。"),
      en = paste0("Output file exists: ", h5ad_path, ". Set overwrite=TRUE to allow overwrite, or choose a new path.")
    ))
  }
  
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
  
  # Create temporary directory for extraction
  if (is.null(temp_dir)) {
    temp_dir <- tempdir()
  }
  extract_dir <- file.path(temp_dir, paste0("1cellbio_", basename(tempfile())))
  dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Ensure cleanup happens even if function fails
  if (cleanup) {
    on.exit({
      if (dir.exists(extract_dir)) {
        unlink(extract_dir, recursive = TRUE)
      }
    })
  }
  
  message("Extracting zip file...")
  # Extract zip file
  utils::unzip(zip_path, exdir = extract_dir)
  
  # Find the sce.json file to locate the experiment data
  sce_json_path <- list.files(extract_dir, pattern = "sce\\.json$", recursive = TRUE, full.names = TRUE)
  
  if (length(sce_json_path) == 0) {
    stop("Could not find sce.json in the provided zip file")
  }
  
  message("Reading experiment metadata...")
  # Read the sce.json file
  sce_json <- jsonlite::read_json(sce_json_path)
  experiment_path <- file.path(extract_dir, sce_json$redirection$targets[[1]]$location)
  
  # Read the experiment.json file
  experiment <- jsonlite::read_json(experiment_path)
  
  message("Reading count matrix...")
  # Read count data
  counts_path <- file.path(extract_dir, experiment$summarized_experiment$assays[[1]]$resource$path)
  counts_matrix <- read_hdf5_sparse_matrix(counts_path)
  
  message("Reading logcounts matrix...")
  # Read logcounts data
  logcounts_path <- file.path(extract_dir, experiment$summarized_experiment$assays[[2]]$resource$path)
  logcounts_matrix <- read_hdf5_matrix(logcounts_path)
  
  message("Reading cell metadata...")
  # Read column data (cell metadata)
  coldata_path <- file.path(extract_dir, experiment$summarized_experiment$column_data$resource$path)
  coldata_df <- read_hdf5_dataframe(coldata_path)
  
  message("Reading gene metadata...")
  # Read row data (gene metadata)
  rowdata_path <- file.path(extract_dir, experiment$summarized_experiment$row_data$resource$path)
  rowdata_df <- read_hdf5_dataframe(rowdata_path)
  
  # Set row and column names
  if ("cell_id" %in% colnames(coldata_df)) {
    rownames(coldata_df) <- icb_make_unique(coldata_df$cell_id, strategy = name_conflict, sep = "-")
  } else {
    rownames(coldata_df) <- seq_len(nrow(coldata_df))
  }
  
  if ("id" %in% colnames(rowdata_df)) {
    rownames(rowdata_df) <- icb_make_unique(rowdata_df$id, strategy = name_conflict, sep = "-")
  } else {
    rownames(rowdata_df) <- seq_len(nrow(rowdata_df))
  }
  
  colnames(counts_matrix) <- rownames(coldata_df)
  rownames(counts_matrix) <- rownames(rowdata_df)
  colnames(logcounts_matrix) <- rownames(coldata_df)
  rownames(logcounts_matrix) <- rownames(rowdata_df)
  
  message("Reading dimensionality reduction data...")
  # Read dimensionality reduction data
  reduced_dims <- list()
  if (!is.null(experiment$single_cell_experiment$reduced_dimensions)) {
    for (reddim in experiment$single_cell_experiment$reduced_dimensions) {
      reddim_name <- reddim$name
      reddim_path <- file.path(extract_dir, reddim$resource$path)
      reduced_dims[[reddim_name]] <- read_hdf5_matrix(reddim_path)
      if (nrow(reduced_dims[[reddim_name]]) == nrow(coldata_df)) {
        rownames(reduced_dims[[reddim_name]]) <- rownames(coldata_df)
      }
    }
  }
  
  message("Creating AnnData object...")
  # Create AnnData object preserving sparse matrix format
  # Keep counts as sparse matrix for memory efficiency
  # Check if logcounts is also sparse and preserve it
  
  # Check matrix types
  counts_is_sparse <- inherits(counts_matrix, "sparseMatrix")
  logcounts_is_sparse <- inherits(logcounts_matrix, "sparseMatrix")
  
  message("Counts matrix type: ", if(counts_is_sparse) "sparse" else "dense")
  message("Logcounts matrix type: ", if(logcounts_is_sparse) "sparse" else "dense")
  
  # Check dimensions match
  message("Counts matrix dimensions: ", nrow(counts_matrix), " genes x ", ncol(counts_matrix), " cells")
  message("ColData dimensions: ", nrow(coldata_df), " cells x ", ncol(coldata_df), " metadata columns")
  message("RowData dimensions: ", nrow(rowdata_df), " genes x ", ncol(rowdata_df), " metadata columns")
  
  # Create AnnData object with sparse matrices
  # For AnnData, X should typically be the raw counts (sparse)
  # and layers can contain processed data
  ad <- anndata::AnnData(
    X = t(counts_matrix),  # Keep sparse format, transpose to have cells as rows
    layers = list(logcounts = t(logcounts_matrix)),  # Keep original format, transpose to have cells as rows
    obs = coldata_df,
    var = rowdata_df
  )
  
  # Add dimensionality reductions
  if (length(reduced_dims) > 0) {
    message("Adding dimensionality reductions...")
    for (reddim_name in names(reduced_dims)) {
      # Convert to lowercase for consistency with Scanpy conventions
      reduc_name <- tolower(reddim_name)
      
      # Special handling for common dimensionality reductions
      if (reddim_name == "TSNE") {
        reduc_name <- "tsne"
      } else if (reddim_name == "UMAP") {
        reduc_name <- "umap"
      } else if (reddim_name == "PCA") {
        reduc_name <- "pca"
      }
      
      # Assign the reduced dimensions
      ad$obsm[[paste0("X_", reduc_name)]] <- reduced_dims[[reddim_name]]
    }
  }
  
  message("Writing h5ad file...")
  # Write to h5ad file
  ad$write_h5ad(h5ad_path)
  
  message("Successfully converted ", zip_path, " to ", h5ad_path)
  
  return(invisible(NULL))
}