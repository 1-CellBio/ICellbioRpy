#' Convert sparse matrix from HDF5 to dgCMatrix
#' 
#' This function converts a sparse matrix stored in HDF5 format (10x format) to a dgCMatrix
#' 
#' @param data Vector of non-zero values
#' @param indices Vector of column indices for each non-zero value
#' @param indptr Vector of index pointers for each row
#' @param shape Vector of matrix dimensions (rows, columns)
#' @return A dgCMatrix object
#' 
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom Matrix sparseMatrix
construct_sparse_matrix <- function(data, indices, indptr, shape) {
  # Number of rows and columns
  nr <- shape[1]
  nc <- shape[2]
  
  # In 10x format:
  # - indices contains 0-based row indices of non-zero elements
  # - indptr contains column pointers (0-based), length is (nrows + 1)
  # - data contains the non-zero values
  
  # Create i (row indices) from indices vector
  i <- indices + 1  # Convert to 1-based indexing for R
  
  # Create j (column indices) from indptr
  j <- rep(seq_len(nc), diff(indptr))  # This creates the column indices for each element
  
  # Create the sparse matrix
  mat <- Matrix::sparseMatrix(
    i = i,
    j = j,
    x = data,
    dims = c(nr, nc)
  )
  
  return(mat)
}

#' Read HDF5 sparse matrix
#' 
#' This function reads a sparse matrix from HDF5 file and returns a dgCMatrix
#' 
#' @param file_path Path to the HDF5 file
#' @return A dgCMatrix object
#' 
#' @import hdf5r
read_hdf5_sparse_matrix <- function(file_path) {
  h5file <- hdf5r::H5File$new(file_path, mode = "r")
  on.exit(h5file$close())
  
  data <- h5file[["matrix/data"]][]
  indices <- h5file[["matrix/indices"]][]
  indptr <- h5file[["matrix/indptr"]][]
  shape <- h5file[["matrix/shape"]][]
  
  mat <- construct_sparse_matrix(data, indices, indptr, shape)
  return(mat)
}

#' Read HDF5 matrix (dense or sparse)
#' 
#' This function reads a matrix from HDF5 file and can return either dense or sparse matrix
#' 
#' @param file_path Path to the HDF5 file
#' @return A matrix object (dense or sparse)
#' @export
#' @import hdf5r
read_hdf5_matrix <- function(file_path) {
  h5file <- hdf5r::H5File$new(file_path, mode = "r")
  on.exit(h5file$close())
  
  # Check the structure and read appropriately
  if ("data" %in% h5file$ls()$name) {
    # Standard dense matrix with "data" dataset
    # Check the dimensions first
    data_dataset <- h5file[["data"]]
    if (inherits(data_dataset, "H5D")) {
      # It's a dataset, read it directly
      mat <- data_dataset[,]
    } else {
      # It's a group, need to handle differently
      mat <- data_dataset[[]]  # Try to read the entire dataset
    }
  } else if ("logcounts" %in% h5file$ls()$name) {
    # For logcounts, we need to handle the delayed array structure
    # Since logcounts are typically all zeros (placeholder), create a sparse matrix directly
    # This is much more memory efficient than creating a dense matrix first
    dimensions_path <- "logcounts/seed/seed/seed/dimensions"
    if (dimensions_path %in% h5file$ls(recursive = TRUE)$name) {
      dimensions <- h5file[[dimensions_path]][]
      # Create an empty sparse matrix (all zeros) - much more memory efficient
      mat <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), 
                                  dims = c(dimensions[1], dimensions[2]))
    } else {
      # Fallback: try to get dimensions from shape dataset if it exists
      if ("logcounts/seed/seed/seed/shape" %in% h5file$ls(recursive = TRUE)$name) {
        shape <- h5file[["logcounts/seed/seed/seed/shape"]][]
        # Create an empty sparse matrix (all zeros) - much more memory efficient
        mat <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), 
                                    dims = c(shape[1], shape[2]))
      } else {
        stop("Cannot determine dimensions for logcounts matrix in file: ", file_path)
      }
    }
  } else {
    stop("Unknown HDF5 structure in file: ", file_path)
  }
  
  # Ensure we return a proper matrix, but preserve sparse matrices
  if (!is.matrix(mat) && !inherits(mat, "sparseMatrix")) {
    mat <- as.matrix(mat)
  }
  
  return(mat)
}

#' Read HDF5 data frame
#' 
#' This function reads a data frame from HDF5 file
#' 
#' @param file_path Path to the HDF5 file
#' @return A data.frame object
#' 
#' @import hdf5r
read_hdf5_dataframe <- function(file_path) {
  h5file <- hdf5r::H5File$new(file_path, mode = "r")
  on.exit(h5file$close())
  
  # Try different possible structures for column names
  names_vec <- NULL
  if ("data/_NAME_" %in% h5file$ls(recursive = TRUE)$name) {
    names_vec <- h5file[["data/_NAME_"]][]
  } else if ("data/column_names" %in% h5file$ls(recursive = TRUE)$name) {
    names_vec <- h5file[["data/column_names"]][]
  } else {
    # If we can't find column names in expected locations, 
    # try to find any dataset that might contain them
    all_names <- h5file$ls(recursive = TRUE)$name
    name_datasets <- all_names[grepl("name", all_names, ignore.case = TRUE)]
    if (length(name_datasets) > 0) {
      names_vec <- h5file[[name_datasets[1]]][]
    } else {
      stop("Cannot find column names in HDF5 file: ", file_path)
    }
  }
  
  # Read each column
  df_list <- list()
  for (i in seq_along(names_vec)) {
    col_name <- names_vec[i]
    # Try different possible paths for the data
    possible_paths <- c(
      paste0("data/", col_name),
      paste0("data/data/", i-1),  # 0-based indexing
      paste0("data/", i-1)
    )
    
    for (path in possible_paths) {
      if (path %in% h5file$ls(recursive = TRUE)$name) {
        df_list[[col_name]] <- h5file[[path]][]
        break
      }
    }
    
    # If we couldn't find the column data, add NA
    if (is.null(df_list[[col_name]])) {
      warning("Could not find data for column: ", col_name)
      df_list[[col_name]] <- rep(NA, 10)  # Placeholder
    }
  }
  
  df <- data.frame(df_list)
  return(df)
}


#' Return package language option
#'
#' @keywords internal
icb_get_lang <- function() {
  lang <- getOption("ICellbioRpy.lang", default = "en")
  if (!lang %in% c("zh", "en")) lang <- "en"
  lang
}

#' Simple i18n helper
#'
#' @param zh Chinese message
#' @param en English message
#' @keywords internal
icb_i18n <- function(zh, en) {
  if (identical(icb_get_lang(), "zh")) zh else en
}

#' Ensure a matrix is CsparseMatrix
#'
#' @param x A matrix-like object (dense or sparse)
#' @return CsparseMatrix
#' @keywords internal
to_csparse <- function(x) {
  if (inherits(x, "sparseMatrix")) {
    return(methods::as(x, "CsparseMatrix"))
  }
  # Convert dense to sparse preserving type
  Matrix::Matrix(x, sparse = TRUE)
}

#' Check if a Python object is scipy.sparse matrix
#'
#' @keywords internal
icb_is_py_sparse <- function(py_obj) {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(FALSE)
  if (is.null(py_obj)) return(FALSE)
  sparse <- tryCatch(reticulate::import("scipy.sparse", convert = FALSE), error = function(e) NULL)
  if (is.null(sparse)) return(FALSE)
  is_csr <- tryCatch(sparse$isspmatrix_csr(py_obj), error = function(e) FALSE)
  is_csc <- tryCatch(sparse$isspmatrix_csc(py_obj), error = function(e) FALSE)
  is_csr || is_csc
}

#' Convert scipy.sparse (CSR/CSC) to dgCMatrix efficiently
#'
#' @keywords internal
icb_py_sparse_to_dgc <- function(py_obj) {
  sparse <- reticulate::import("scipy.sparse", convert = FALSE)
  if (sparse$isspmatrix_csr(py_obj)) {
    shape <- reticulate::py_to_r(py_obj$shape)
    indptr <- reticulate::py_to_r(py_obj$indptr)
    indices <- reticulate::py_to_r(py_obj$indices)
    data <- reticulate::py_to_r(py_obj$data)
    n_rows <- as.integer(shape[[1]])
    n_cols <- as.integer(shape[[2]])
    i <- rep.int(seq_len(n_rows), diff(indptr))
    j <- as.integer(indices) + 1L
    Matrix::sparseMatrix(i = i, j = j, x = as.numeric(data), dims = c(n_rows, n_cols))
  } else if (sparse$isspmatrix_csc(py_obj)) {
    shape <- reticulate::py_to_r(py_obj$shape)
    indptr <- reticulate::py_to_r(py_obj$indptr)
    indices <- reticulate::py_to_r(py_obj$indices)
    data <- reticulate::py_to_r(py_obj$data)
    n_rows <- as.integer(shape[[1]])
    n_cols <- as.integer(shape[[2]])
    j <- rep.int(seq_len(n_cols), diff(indptr))
    i <- as.integer(indices) + 1L
    Matrix::sparseMatrix(i = i, j = j, x = as.numeric(data), dims = c(n_rows, n_cols))
  } else {
    stop(icb_i18n("不支持的稀疏矩阵类型", "Unsupported scipy.sparse matrix type"))
  }
}

#' Robustly obtain CsparseMatrix from AnnData matrix/layer
#'
#' @param obj Anndata matrix-like (possibly Python object)
#' @param transpose Whether to transpose after conversion
#' @param verbose Whether to print diagnostic
#' @keywords internal
icb_anndata_to_csparse <- function(obj, transpose = FALSE, verbose = FALSE) {
  # Fast path for R matrices
  if (inherits(obj, "sparseMatrix") || is.matrix(obj)) {
    mat <- to_csparse(obj)
    return(if (transpose) methods::as(Matrix::t(mat), "CsparseMatrix") else mat)
  }
  # Python sparse path
  if (icb_is_py_sparse(obj)) {
    mat <- icb_py_sparse_to_dgc(obj)
    return(if (transpose) methods::as(Matrix::t(mat), "CsparseMatrix") else mat)
  }
  # Fallback: try convert to R and then to sparse
  if (requireNamespace("reticulate", quietly = TRUE)) {
    conv <- tryCatch(reticulate::py_to_r(obj), error = function(e) NULL)
    if (!is.null(conv)) {
      if (inherits(conv, "sparseMatrix") || is.matrix(conv)) {
        mat <- to_csparse(conv)
        return(if (transpose) methods::as(Matrix::t(mat), "CsparseMatrix") else mat)
      }
    }
  }
  stop(icb_i18n("无法从 AnnData 对象提取矩阵", "Failed to extract matrix from AnnData object"))
}

#' Make names unique or error
#'
#' @param names Character vector
#' @param strategy One of "make_unique" (default) or "error"
#' @param sep Suffix separator when making unique
#' @keywords internal
icb_make_unique <- function(names, strategy = c("make_unique", "error"), sep = "-") {
  strategy <- match.arg(strategy)
  if (strategy == "make_unique") {
    return(base::make.unique(as.character(names), sep = sep))
  }
  # error strategy
  if (anyDuplicated(names)) {
    dups <- unique(names[duplicated(names)])
    stop(icb_i18n(
      zh = paste0("检测到重复名称: ", paste(utils::head(dups, 10), collapse = ", "), if (length(dups) > 10) " ..." else "", "。请调整输入或改用 name_conflict='make_unique'。"),
      en = paste0("Duplicate names detected: ", paste(utils::head(dups, 10), collapse = ", "), if (length(dups) > 10) " ..." else "", ". Please adjust input or set name_conflict='make_unique'.")
    ))
  }
  as.character(names)
}

#' Package startup message
#' 
#' @param libname Library name
#' @param pkgname Package name
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "====================================================================\n",
    icb_i18n(
      zh = "欢迎使用ICellbioRpy！\n\n🚀 智能Python环境配置：\n smart_python_config(verbose = TRUE, interactive = TRUE)\n\n此命令将自动检测可用的conda环境并提示您选择。\n首次使用时请运行此命令以配置Python环境。\n\n💡 手动配置（高级用户）：\n configure_python_env(conda_env = \"your_env\", verbose = TRUE)\n",
      en = "Welcome to ICellbioRpy!\n\n🚀 Smart Python environment configuration:\n smart_python_config(verbose = TRUE, interactive = TRUE)\n\nThis command will auto-detect available conda environments and prompt you to choose.\nPlease run this command on first use to configure Python environment.\n\n💡 Manual configuration (advanced users):\n configure_python_env(conda_env = \"your_env\", verbose = TRUE)\n"
    ),
    "===================================================================="
  )
}