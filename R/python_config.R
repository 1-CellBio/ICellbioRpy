#' Configure Python environment for anndata
#'
#' This function configures the Python environment to use existing anndata installation
#' and prevents automatic package installation/upgrade.
#'
#' @param python_path Optional path to Python executable. If NULL, uses current environment.
#' @param conda_env Optional conda environment name. If NULL, uses current environment.
#' @param verbose Logical, whether to print configuration messages (default: TRUE)
#'
#' @return Invisible TRUE if configuration successful
#'
#' @examples
#' \dontrun{
#' # Use current Python environment (default)
#' configure_python_env()
#' 
#' # Use specific conda environment
#' configure_python_env(conda_env = "your_env_name")
#' 
#' # Use specific Python path
#' configure_python_env(python_path = "/usr/local/bin/python3")
#' }
#'
#' @export
configure_python_env <- function(python_path = NULL, conda_env = NULL, verbose = TRUE) {
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required but not installed. Please install it with: install.packages('reticulate')")
  }
  
  # Disable automatic configuration and installation
  Sys.setenv(RETICULATE_AUTOCONFIGURE = "FALSE")
  Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
  
  if (verbose) cat("Configuring Python environment...\n")
  
  # Configure Python environment
  if (!is.null(conda_env)) {
    if (verbose) cat("Using specified conda environment:", conda_env, "\n")
    reticulate::use_condaenv(conda_env, required = FALSE)
  } else if (!is.null(python_path)) {
    if (verbose) cat("Using specified Python path:", python_path, "\n")
    reticulate::use_python(python_path, required = FALSE)
  } else {
    # Use current environment - no explicit configuration needed
    if (verbose) cat("Using current Python environment...\n")
  }
  
  # Check if anndata is available in the configured environment
  tryCatch({
    py_config <- reticulate::py_config()
    if (verbose) {
      cat("Python configuration:\n")
      cat("  - Python:", py_config$python, "\n")
      cat("  - Version:", py_config$version, "\n")
    }
    
    # Test anndata import
    anndata_available <- reticulate::py_module_available("anndata")
    if (anndata_available) {
      tryCatch({
        anndata <- reticulate::import("anndata")
        anndata_version <- anndata$`__version__`
        if (verbose) cat("  - anndata version:", anndata_version, "\n")
      }, error = function(e) {
        if (verbose) cat("  - anndata found but version check failed\n")
      })
    } else {
      # anndata not available - provide helpful message
      warning("anndata module not found in current Python environment.\n",
              "Please install anndata or specify a conda environment that has anndata installed:\n",
              "  configure_python_env(conda_env = \"your_env_name\")\n",
              "Or install anndata in your current environment:\n",
              "  pip install anndata")
      return(invisible(FALSE))
    }
    
  }, error = function(e) {
    error_msg <- if(is.character(e$message)) e$message else as.character(e$message)
    warning("Failed to configure Python environment: ", error_msg)
    return(invisible(FALSE))
  })
  
  if (verbose) cat("✓ Python environment configured successfully\n")
  return(invisible(TRUE))
}

#' Check if Python environment is properly configured for anndata
#'
#' @param verbose Logical, whether to print status messages (default: TRUE)
#' @return Logical, TRUE if anndata is available
#'
#' @export
check_anndata_available <- function(verbose = TRUE) {
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    if (verbose) cat("reticulate package not available\n")
    return(FALSE)
  }
  
  # Check if Python is configured
  if (!reticulate::py_available()) {
    if (verbose) cat("Python not available\n")
    return(FALSE)
  }
  
  # Check if anndata module is available
  anndata_available <- reticulate::py_module_available("anndata")
  
  if (verbose) {
    if (anndata_available) {
      tryCatch({
        anndata <- reticulate::import("anndata")
        anndata_version <- anndata$`__version__`
        cat("✓ anndata available (version:", anndata_version, ")\n")
      }, error = function(e) {
        cat("✓ anndata module found but version check failed\n")
      })
    } else {
      cat("✗ anndata module not available\n")
    }
  }
  
  return(anndata_available)
}