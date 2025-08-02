#' Convert h5ad file to Seurat object
#'
#' This function reads an h5ad file and converts it to a Seurat object,
#' preserving sparse matrix format for memory efficiency. The conversion
#' goes through SingleCellExperiment as an intermediate step.
#'
#' @param h5ad_file Path to the h5ad file
#' @param use_x_as Character string specifying which data to use as the main expression matrix.
#'   Options: "logcounts" (default), "counts", or "auto". If "auto", automatically detect based on data.
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#'
#' @return A Seurat object with:
#'   \itemize{
#'     \item RNA assay with counts and data slots (preserved as sparse matrices)
#'     \item Cell metadata in meta.data
#'     \item Dimensionality reductions (PCA, tSNE, UMAP if available)
#'   }
#'
#' @examples
#' \dontrun{
#' # Convert h5ad to Seurat
#' seurat_obj <- h5ad_to_seurat("data.h5ad")
#' 
#' # Use X matrix as counts instead of logcounts
#' seurat_obj <- h5ad_to_seurat("data.h5ad", use_x_as = "counts")
#' }
#'
#' @export
#' @importFrom Seurat as.Seurat
h5ad_to_seurat <- function(h5ad_file, use_x_as = "auto", verbose = TRUE) {
  
  # Check if required packages are available
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required but not installed. Please install it with: install.packages('Seurat')")
  }
  
  if (verbose) cat("Converting h5ad to Seurat via SingleCellExperiment...\n")
  
  # First convert to SingleCellExperiment
  sce <- h5ad_to_sce(h5ad_file, use_x_as = use_x_as, verbose = verbose)
  
  if (verbose) cat("Converting SingleCellExperiment to Seurat...\n")
  
  # Check available assays and use the first one as default
  available_assays <- names(SummarizedExperiment::assays(sce))
  if (length(available_assays) == 0) {
    stop("No assays found in SingleCellExperiment object")
  }
  
  # Use the first available assay for Seurat conversion
  default_assay <- available_assays[1]
  
  # Create Seurat object with the available assay
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = SummarizedExperiment::assay(sce, default_assay),
    meta.data = as.data.frame(SummarizedExperiment::colData(sce))
  )
  
  # Add other assays if available
  if (length(available_assays) > 1) {
    for (assay_name in available_assays[-1]) {
      seurat_obj[[assay_name]] <- Seurat::CreateAssayObject(
        data = SummarizedExperiment::assay(sce, assay_name)
      )
    }
  }
  
  # Add dimensionality reductions if available
  if (length(SingleCellExperiment::reducedDims(sce)) > 0) {
    for (reduction_name in names(SingleCellExperiment::reducedDims(sce))) {
      # Clean reduction name for Seurat (remove X_ prefix if present)
      clean_name <- gsub("^X_", "", reduction_name)
      
      reduction_data <- SingleCellExperiment::reducedDims(sce)[[reduction_name]]
      
      # Create DimReduc object
      seurat_obj[[clean_name]] <- Seurat::CreateDimReducObject(
        embeddings = reduction_data,
        key = paste0(toupper(substr(clean_name, 1, 1)), substr(clean_name, 2, nchar(clean_name)), "_"),
        assay = Seurat::DefaultAssay(seurat_obj)
      )
    }
  }
  
  if (verbose) {
    cat("✓ Conversion to Seurat completed successfully!\n")
    cat("Seurat object created with:\n")
    cat("  - Assays:", paste(names(seurat_obj@assays), collapse = ", "), "\n")
    cat("  - Metadata columns:", ncol(seurat_obj@meta.data), "\n")
    if (length(seurat_obj@reductions) > 0) {
      cat("  - Reductions:", paste(names(seurat_obj@reductions), collapse = ", "), "\n")
    }
  }
  
  return(seurat_obj)
}