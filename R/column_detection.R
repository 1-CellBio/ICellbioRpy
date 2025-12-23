#' Column name detection utilities for 1CellBio data conversion
#'
#' These functions help automatically detect appropriate column names for gene
#' and cell identifiers in 1CellBio data.
#'
#' @param colnames Character vector of column names to search through
#' @return The first matching column name or NULL if no match found
#' @keywords internal
#'
#' @examples
#' # Detect gene ID column
#' gene_cols <- c("gene_id", "symbol", "id", "expression")
#' detect_gene_id_column(gene_cols)  # Returns "gene_id"
#'
#' # Detect cell ID column
#' cell_cols <- c("cell_id", "barcode", "sample_id")
#' detect_cell_id_column(cell_cols)  # Returns "cell_id"

#' Detect gene identifier column
#'
#' Attempts to automatically detect which column contains gene identifiers.
#' Searches for common patterns in order of preference.
#'
#' @param colnames Character vector of column names
#' @return The best matching column name or NULL
#' @export
detect_gene_id_column <- function(colnames) {
  # Common patterns for gene identifiers, in order of preference
  patterns <- c("id", "gene_id", "gene_symbol", "symbol", "gene_name", "gene")

  for (pattern in patterns) {
    if (pattern %in% colnames) {
      return(pattern)
    }
  }

  return(NULL)
}

#' Detect cell identifier column
#'
#' Attempts to automatically detect which column contains cell identifiers.
#' Searches for common patterns in order of preference.
#'
#' @param colnames Character vector of column names
#' @return The best matching column name or NULL
#' @export
detect_cell_id_column <- function(colnames) {
  # Common patterns for cell identifiers, in order of preference
  patterns <- c("cell_id", "cell_barcode", "barcode", "cellid", "cell", "sample_id")

  for (pattern in patterns) {
    if (pattern %in% colnames) {
      return(pattern)
    }
  }

  return(NULL)
}

#' Show available column options
#'
#' Display available options for gene and cell identifier columns
#' with detected suggestions.
#'
#' @param object A 1CellbioData object
#' @return Invisible list with suggested column names
#' @export
show_column_options <- function(object) {
  # Validate input
  if (!inherits(object, "1CellbioData")) {
    stop("Object must be of class 1CellbioData")
  }

  # Read metadata
  coldata_path <- file.path(object$base_path, object$experiment$summarized_experiment$column_data$resource$path)
  rowdata_path <- file.path(object$base_path, object$experiment$summarized_experiment$row_data$resource$path)

  coldata_df <- read_hdf5_dataframe(coldata_path)
  rowdata_df <- read_hdf5_dataframe(rowdata_path)

  # Detect suggestions
  gene_suggestion <- detect_gene_id_column(colnames(rowdata_df))
  cell_suggestion <- detect_cell_id_column(colnames(coldata_df))

  # Display results
  cat("\n=== Column Detection Results ===\n\n")

  cat("Available gene identifier columns (rownames):\n")
  cat("  ", paste(colnames(rowdata_df), collapse = ", "), "\n")
  if (!is.null(gene_suggestion)) {
    cat("  → Detected: ", gene_suggestion, "\n", sep = "")
  } else {
    cat("  → No standard gene ID column detected\n")
  }

  cat("\nAvailable cell identifier columns (colnames):\n")
  cat("  ", paste(colnames(coldata_df), collapse = ", "), "\n")
  if (!is.null(cell_suggestion)) {
    cat("  → Detected: ", cell_suggestion, "\n", sep = "")
  } else {
    cat("  → No standard cell ID column detected\n")
  }

  cat("\nSuggested usage:\n")
  if (!is.null(gene_suggestion) && !is.null(cell_suggestion)) {
    cat("  as.Seurat.1CB(data, rownames = '", gene_suggestion,
        "', colnames = '", cell_suggestion, "')\n", sep = "")
    cat("  # Or simply use auto-detection:\n")
    cat("  as.Seurat.1CB(data)\n")
  } else {
    cat("  Manually specify columns:\n")
    cat("  as.Seurat.1CB(data, rownames = 'your_gene_column', colnames = 'your_cell_column')\n")
  }

  invisible(list(
    gene_columns = colnames(rowdata_df),
    cell_columns = colnames(coldata_df),
    gene_suggestion = gene_suggestion,
    cell_suggestion = cell_suggestion
  ))
}

#' Validate column names
#'
#' Check if specified column names exist in the data and return
#' validated names or error messages.
#'
#' @param object A 1CellbioData object
#' @param rownames Column name for gene identifiers
#' @param colnames Column name for cell identifiers
#' @param auto_detect Whether to try auto-detection if names are NULL
#' @param verbose Whether to print messages when auto-generating column names
#' @return List with validated column names
#' @keywords internal
validate_column_names <- function(object, rownames = NULL, colnames = NULL, auto_detect = TRUE, verbose = TRUE) {
  # Read metadata
  coldata_path <- file.path(object$base_path, object$experiment$summarized_experiment$column_data$resource$path)
  rowdata_path <- file.path(object$base_path, object$experiment$summarized_experiment$row_data$resource$path)

  coldata_df <- read_hdf5_dataframe(coldata_path)
  rowdata_df <- read_hdf5_dataframe(rowdata_path)

  # Auto-detect if needed
  if (is.null(rownames) && auto_detect) {
    rownames <- detect_gene_id_column(colnames(rowdata_df))
  }

  if (is.null(colnames) && auto_detect) {
    colnames <- detect_cell_id_column(colnames(coldata_df))
  }

  # Handle undetected rownames - auto-select first column
  if (is.null(rownames)) {
    rownames <- colnames(rowdata_df)[1]
    if (verbose) {
      cat(icb_i18n(
        paste0("⚠️  无法检测到基因标识符列，自动使用第一列: ", rownames, "\n"),
        paste0("⚠️  Cannot detect gene identifier column, auto-using first column: ", rownames, "\n")
      ))
    }
  }

  if (!rownames %in% colnames(rowdata_df)) {
    stop("Gene identifier column '", rownames, "' not found. Available options: ",
         paste(colnames(rowdata_df), collapse = ", "))
  }

  # Handle undetected colnames - auto-select first column
  if (is.null(colnames)) {
    colnames <- colnames(coldata_df)[1]
    if (verbose) {
      cat(icb_i18n(
        paste0("⚠️  无法检测到细胞标识符列，自动使用第一列: ", colnames, "\n"),
        paste0("⚠️  Cannot detect cell identifier column, auto-using first column: ", colnames, "\n")
      ))
    }
  }

  if (!colnames %in% colnames(coldata_df)) {
    stop("Cell identifier column '", colnames, "' not found. Available options: ",
         paste(colnames(coldata_df), collapse = ", "))
  }

  return(list(rownames = rownames, colnames = colnames))
}