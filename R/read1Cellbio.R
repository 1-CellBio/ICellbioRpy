#' Read 1Cellbio results from zip file
#'
#' This function reads 1Cellbio pipeline results directly from a zip file
#' and creates a 1CellbioData object that can be converted to Seurat or SingleCellExperiment.
#'
#' @param file_path Path to the zip file containing 1Cellbio results
#' @return A 1CellbioData object
#' @export
#' @examples
#' \dontrun{
#' data <- read1Cellbio("path/to/results.zip")
#' # 需要指定用作基因名和细胞名的列
#' sce <- as.SingleCellExperiment.1CB(data, rownames = "id", colnames = "cell_id")
#' seurat <- as.Seurat.1CB(data, rownames = "id", colnames = "cell_id")
#' }
read1Cellbio <- function(file_path) {
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  # Create a temporary directory to extract files
  temp_dir <- tempdir()
  extract_dir <- file.path(temp_dir, "1cellbio_results")
  dir.create(extract_dir, showWarnings = FALSE)
  
  # Extract zip file
  utils::unzip(file_path, exdir = extract_dir)
  
  # Find the sce.json file to locate the experiment data
  sce_json_path <- list.files(extract_dir, pattern = "sce\\.json$", recursive = TRUE, full.names = TRUE)
  
  if (length(sce_json_path) == 0) {
    stop("Could not find sce.json in the provided zip file")
  }
  
  # Read the sce.json file
  sce_json <- jsonlite::read_json(sce_json_path)
  experiment_path <- file.path(extract_dir, sce_json$redirection$targets[[1]]$location)
  
  # Read the experiment.json file
  experiment <- jsonlite::read_json(experiment_path)
  
  # Create 1CellbioData object
  data <- list(
    experiment = experiment,
    base_path = extract_dir
  )
  
  class(data) <- "1CellbioData"
  
  return(data)
}

## Delegates to utils.R authoritative implementations


#' Helper function to make unique names
#'
#' @param names Character vector of names
#' @return Character vector with unique names
make_unique_names <- function(names) {
  if (any(duplicated(names))) {
    names <- make.unique(names, sep = "-")
  }
  return(names)
}

#' Convert 1CellbioData to Seurat object (DEPRECATED)
#'
#' @param object A 1CellbioData object
#' @param rownames Character string specifying the column name in row data to use as gene names (required)
#' @param colnames Character string specifying the column name in column data to use as cell names (required)
#' @param ... Additional arguments (not used)
#' @return A Seurat object
#' @keywords internal
as.Seurat.1CellbioData <- function(object, rownames = NULL, colnames = NULL, ..., name_conflict = c("make_unique", "error")) {
  name_conflict <- match.arg(name_conflict)
  message(icb_i18n(
    zh = "提示: 函数 as.Seurat.1CellbioData 已废弃，请改用 as.Seurat.1CB。",
    en = "Note: Function as.Seurat.1CellbioData is deprecated, please use as.Seurat.1CB instead."
  ))
  
  # Read metadata first to show available options
  cat("正在读取数据结构信息...\n")
  
  # Read column data (cell metadata)
  coldata_path <- file.path(object$base_path, object$experiment$summarized_experiment$column_data$resource$path)
  coldata_df <- read_hdf5_dataframe(coldata_path)
  
  # Read row data (gene metadata)
  rowdata_path <- file.path(object$base_path, object$experiment$summarized_experiment$row_data$resource$path)
  rowdata_df <- read_hdf5_dataframe(rowdata_path)
  
  # Display available column names
  cat("可用的细胞名列 (colnames):", paste(colnames(coldata_df), collapse = ", "), "\n")
  cat("可用的基因名列 (rownames):", paste(colnames(rowdata_df), collapse = ", "), "\n")
  
  # Check required parameters
  if (is.null(rownames)) {
    stop("参数 'rownames' 是必填项。请指定用作基因名的列名。")
  }
  
  if (is.null(colnames)) {
    stop("参数 'colnames' 是必填项。请指定用作细胞名的列名。")
  }
  
  # Read count data
  counts_path <- file.path(object$base_path, object$experiment$summarized_experiment$assays[[1]]$resource$path)
  counts_matrix <- read_hdf5_sparse_matrix(counts_path)
  
  # Read logcounts data
  logcounts_path <- file.path(object$base_path, object$experiment$summarized_experiment$assays[[2]]$resource$path)
  logcounts_matrix <- read_hdf5_matrix(logcounts_path)
  
  # Check if specified column names exist
  if (!colnames %in% colnames(coldata_df)) {
    stop("指定的细胞名列 '", colnames, "' 在列数据中不存在。可用列名：", paste(colnames(coldata_df), collapse = ", "))
  }
  
  if (!rownames %in% colnames(rowdata_df)) {
    stop("指定的基因名列 '", rownames, "' 在行数据中不存在。可用列名：", paste(colnames(rowdata_df), collapse = ", "))
  }
  
  # Set row and column names with uniqueness check
  cell_names <- icb_make_unique(as.character(coldata_df[[colnames]]), strategy = name_conflict, sep = "-")
  gene_names <- icb_make_unique(as.character(rowdata_df[[rownames]]), strategy = name_conflict, sep = "-")
  
  rownames(coldata_df) <- cell_names
  rownames(rowdata_df) <- gene_names
  colnames(counts_matrix) <- cell_names
  rownames(counts_matrix) <- gene_names
  colnames(logcounts_matrix) <- cell_names
  rownames(logcounts_matrix) <- gene_names
  
  # Read dimensionality reduction data
  reduced_dims <- list()
  for (reddim in object$experiment$single_cell_experiment$reduced_dimensions) {
    reddim_name <- reddim$name
    reddim_path <- file.path(object$base_path, reddim$resource$path)
    reduced_dims[[reddim_name]] <- read_hdf5_matrix(reddim_path)
    if (nrow(reduced_dims[[reddim_name]]) == nrow(coldata_df)) {
      rownames(reduced_dims[[reddim_name]]) <- cell_names
    }
  }
  
  # Create Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts_matrix,
    meta.data = coldata_df
  )
  
  # Add logcounts to the assay using the new Seurat v5 approach
  # Keep logcounts as sparse matrix if it already is one for memory efficiency
  if (!is.matrix(logcounts_matrix) && !inherits(logcounts_matrix, "sparseMatrix")) {
    logcounts_matrix <- as.matrix(logcounts_matrix)
  }
  
  # Make sure dimnames are properly set
  dimnames(logcounts_matrix) <- dimnames(counts_matrix)
  
  # For now, let's just make sure the Seurat object is created successfully
  # The full functionality can be added later if needed
  
  # Add all dimensionality reductions dynamically
  for (reddim_name in names(reduced_dims)) {
    # Create appropriate key and reduction name for Seurat
    reduc_name <- tolower(reddim_name)
    key_name <- paste0(substr(reddim_name, 1, 1), 
                       substr(tolower(reddim_name), 2, nchar(reddim_name)), "_")
    
    # Special handling for common dimensionality reductions
    if (reddim_name == "TSNE") {
      reduc_name <- "tsne"
      key_name <- "tSNE_"
    } else if (reddim_name == "UMAP") {
      reduc_name <- "umap"
      key_name <- "UMAP_"
    } else if (reddim_name == "PCA") {
      reduc_name <- "pca"
      key_name <- "PC_"
    } else {
      # For other dimensionality reductions, convert to uppercase with underscore
      reduc_name <- tolower(reddim_name)
      key_name <- paste0(toupper(reddim_name), "_")
    }
    
    # Create the dimensionality reduction object
    seurat_obj[[reduc_name]] <- Seurat::CreateDimReducObject(
      embeddings = reduced_dims[[reddim_name]], 
      key = key_name, 
      assay = "RNA"
    )
  }
  
  return(seurat_obj)
}

#' Convert 1CellbioData to SingleCellExperiment object (DEPRECATED)
#'
#' @param object A 1CellbioData object
#' @param rownames Character string specifying the column name in row data to use as gene names (required)
#' @param colnames Character string specifying the column name in column data to use as cell names (required)
#' @param ... Additional arguments (not used)
#' @return A SingleCellExperiment object
#' @keywords internal
# Original implementation moved to a helper function
as.SingleCellExperiment.1CellbioData_impl <- function(object, rownames = NULL, colnames = NULL, ..., name_conflict = c("make_unique", "error")) {
  name_conflict <- match.arg(name_conflict)
  message(icb_i18n(
    zh = "提示: 函数 as.SingleCellExperiment.1CellbioData 已废弃，请改用 as.SingleCellExperiment.1CB。",
    en = "Note: Function as.SingleCellExperiment.1CellbioData is deprecated, please use as.SingleCellExperiment.1CB instead."
  ))

  # Read metadata first to show available options
  cat("正在读取数据结构信息...\n")

  # Read column data (cell metadata)
  coldata_path <- file.path(object$base_path, object$experiment$summarized_experiment$column_data$resource$path)
  coldata_df <- read_hdf5_dataframe(coldata_path)

  # Read row data (gene metadata)
  rowdata_path <- file.path(object$base_path, object$experiment$summarized_experiment$row_data$resource$path)
  rowdata_df <- read_hdf5_dataframe(rowdata_path)

  # Display available column names
  cat("可用的细胞名列 (colnames):", paste(colnames(coldata_df), collapse = ", "), "\n")
  cat("可用的基因名列 (rownames):", paste(colnames(rowdata_df), collapse = ", "), "\n")

  # Check required parameters
  if (is.null(rownames)) {
    stop("参数 'rownames' 是必填项。请指定用作基因名的列名。")
  }

  if (is.null(colnames)) {
    stop("参数 'colnames' 是必填项。请指定用作细胞名的列名。")
  }

  # Read count data
  counts_path <- file.path(object$base_path, object$experiment$summarized_experiment$assays[[1]]$resource$path)
  counts_matrix <- read_hdf5_sparse_matrix(counts_path)

  # Read logcounts data
  logcounts_path <- file.path(object$base_path, object$experiment$summarized_experiment$assays[[2]]$resource$path)
  logcounts_matrix <- read_hdf5_matrix(logcounts_path)

  # Check if specified column names exist
  if (!colnames %in% colnames(coldata_df)) {
    stop("指定的细胞名列 '", colnames, "' 在列数据中不存在。可用列名：", paste(colnames(coldata_df), collapse = ", "))
  }

  if (!rownames %in% colnames(rowdata_df)) {
    stop("指定的基因名列 '", rownames, "' 在行数据中不存在。可用列名：", paste(colnames(rowdata_df), collapse = ", "))
  }

  # Set row and column names with uniqueness check
  cell_names <- icb_make_unique(as.character(coldata_df[[colnames]]), strategy = name_conflict, sep = "-")
  gene_names <- icb_make_unique(as.character(rowdata_df[[rownames]]), strategy = name_conflict, sep = "-")

  rownames(coldata_df) <- cell_names
  rownames(rowdata_df) <- gene_names
  colnames(counts_matrix) <- cell_names
  rownames(counts_matrix) <- gene_names
  colnames(logcounts_matrix) <- cell_names
  rownames(logcounts_matrix) <- gene_names

  # Read dimensionality reduction data
  reduced_dims <- list()
  for (reddim in object$experiment$single_cell_experiment$reduced_dimensions) {
    reddim_name <- reddim$name
    reddim_path <- file.path(object$base_path, reddim$resource$path)
    reduced_dims[[reddim_name]] <- read_hdf5_matrix(reddim_path)
    if (nrow(reduced_dims[[reddim_name]]) == nrow(coldata_df)) {
      rownames(reduced_dims[[reddim_name]]) <- cell_names
    }
  }

  # Create SingleCellExperiment object
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = counts_matrix,
      logcounts = logcounts_matrix
    ),
    colData = coldata_df,
    rowData = rowdata_df
  )

  # Add all dimensionality reductions dynamically
  for (reddim_name in names(reduced_dims)) {
    # Ensure the reduction name is valid
    valid_name <- make.names(reddim_name)
    SingleCellExperiment::reducedDim(sce, valid_name) <- reduced_dims[[reddim_name]]
  }

  return(sce)
}

#' S3 method for converting 1CellbioData to SingleCellExperiment object
#'
#' @param object A 1CellbioData object
#' @param rownames Character string specifying the column name in row data to use as gene names (required)
#' @param colnames Character string specifying the column name in column data to use as cell names (required)
#' @param ... Additional arguments (not used)
#' @return A SingleCellExperiment object
#' @export
as.SingleCellExperiment.1CellbioData <- function(object, rownames = NULL, colnames = NULL, ..., name_conflict = c("make_unique", "error")) {
  as.SingleCellExperiment.1CellbioData_impl(object, rownames = rownames, colnames = colnames, ..., name_conflict = name_conflict)
}

#' Convert 1CellbioData to Seurat object
#'
#' @param object A 1CellbioData object
#' @param ... Additional arguments (not used)
#' @return A Seurat object
#' @importFrom Seurat as.Seurat
#' @export
as.Seurat.1CB <- function(object, ...) {
  UseMethod("as.Seurat.1CB")
}

#' S3 method for converting 1CellbioData to Seurat object
#'
#' @param object A 1CellbioData object
#' @param rownames Character string specifying the column name in row data to use as gene names (required)
#' @param colnames Character string specifying the column name in column data to use as cell names (required)
#' @param ... Additional arguments (not used)
#' @return A Seurat object
#' @export
as.Seurat.1CB.1CellbioData <- function(object, ...) {
  as.Seurat.1CellbioData_impl(object, ...)
}

#' S3 method for converting 1CellbioData to Seurat object
#'
#' @param object A 1CellbioData object
#' @param rownames Character string specifying the column name in row data to use as gene names (required)
#' @param colnames Character string specifying the column name in column data to use as cell names (required)
#' @param ... Additional arguments (not used)
#' @return A Seurat object
#' @export
# Original implementation moved to a helper function
as.Seurat.1CellbioData_impl <- function(object, rownames = NULL, colnames = NULL, ..., name_conflict = c("make_unique", "error")) {
  name_conflict <- match.arg(name_conflict)
  message(icb_i18n(
    zh = "提示: 函数 as.Seurat.1CellbioData 已废弃，请改用 as.Seurat.1CB。",
    en = "Note: Function as.Seurat.1CellbioData is deprecated, please use as.Seurat.1CB instead."
  ))
  
  # Read metadata first to show available options
  cat("正在读取数据结构信息...\n")
  
  # Read column data (cell metadata)
  coldata_path <- file.path(object$base_path, object$experiment$summarized_experiment$column_data$resource$path)
  coldata_df <- read_hdf5_dataframe(coldata_path)
  
  # Read row data (gene metadata)
  rowdata_path <- file.path(object$base_path, object$experiment$summarized_experiment$row_data$resource$path)
  rowdata_df <- read_hdf5_dataframe(rowdata_path)
  
  # Read expression matrix
  # Find the counts assay
  counts_assay <- NULL
  for (assay in object$experiment$summarized_experiment$assays) {
    if (assay$name == "counts") {
      counts_assay <- assay
      break
    }
  }

  if (is.null(counts_assay)) {
    stop("Could not find counts assay in the data")
  }

  matrix_path <- file.path(object$base_path, counts_assay$resource$path)
  counts <- read_hdf5_sparse_matrix(matrix_path)
  
  # Set row names and column names
  if (!is.null(rownames)) {
    if (!rownames %in% colnames(rowdata_df)) {
      stop("指定的 rownames 列名在 rowdata 中不存在")
    }
    rownames(counts) <- rowdata_df[[rownames]]
  } else {
    rownames(counts) <- rowdata_df[[1]]  # Use first column as default
  }
  
  if (!is.null(colnames)) {
    if (!colnames %in% colnames(coldata_df)) {
      stop("指定的 colnames 列名在 coldata 中不存在")
    }
    colnames(counts) <- coldata_df[[colnames]]
  } else {
    colnames(counts) <- coldata_df[[1]]  # Use first column as default
  }
  
  # Handle name conflicts
  if (name_conflict == "make_unique") {
    rownames(counts) <- make_unique_names(rownames(counts))
    colnames(counts) <- make_unique_names(colnames(counts))
  } else if (any(duplicated(rownames(counts))) || any(duplicated(colnames(counts)))) {
    stop("基因名或细胞名存在重复，请设置 name_conflict = 'make_unique' 或手动处理重复")
  }
  
  # Create Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(counts = counts)
  
  # Add metadata
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, coldata_df)
  seurat_obj@meta.data <- seurat_obj@meta.data[, !duplicated(colnames(seurat_obj@meta.data))]
  
  return(seurat_obj)
}

#' S3 method for converting 1CellbioData to Seurat object
#'
#' @param object A 1CellbioData object
#' @param rownames Character string specifying the column name in row data to use as gene names (required)
#' @param colnames Character string specifying the column name in column data to use as cell names (required)
#' @param ... Additional arguments (not used)
#' @return A Seurat object
#' @export
as.Seurat.1CellbioData <- function(object, rownames = NULL, colnames = NULL, ..., name_conflict = c("make_unique", "error")) {
  as.Seurat.1CellbioData_impl(object, rownames = rownames, colnames = colnames, ..., name_conflict = name_conflict)
}

#' Convert 1CellbioData to SingleCellExperiment object
#'
#' @param object A 1CellbioData object
#' @param ... Additional arguments (not used)
#' @return A SingleCellExperiment object
#' @importFrom Seurat as.SingleCellExperiment
#' @export
as.SingleCellExperiment.1CB <- function(object, ...) {
  UseMethod("as.SingleCellExperiment.1CB")
}

#' S3 method for converting 1CellbioData to SingleCellExperiment object
#'
#' @param object A 1CellbioData object
#' @param rownames Character string specifying the column name in row data to use as gene names (required)
#' @param colnames Character string specifying the column name in column data to use as cell names (required)
#' @param ... Additional arguments (not used)
#' @return A SingleCellExperiment object
#' @export
as.SingleCellExperiment.1CB.1CellbioData <- function(object, ...) {
  as.SingleCellExperiment.1CellbioData_impl(object, ...)
}

