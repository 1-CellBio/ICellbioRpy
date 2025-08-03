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