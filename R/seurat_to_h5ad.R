#' Convert Seurat object to h5ad file
#'
#' This function converts a Seurat object to an h5ad file format using the anndata package.
#' It preserves sparse matrix format, metadata, and dimensionality reductions for memory efficiency.
#'
#' @param seurat_obj A Seurat object to convert
#' @param output_file Character string specifying the output h5ad file path
#' @param default_assay Character string specifying which assay to use as the main X matrix (default: "RNA")
#' @param layer Character string specifying which layer to use from the assay ("data" or "counts", default: "data")
#' @param include_reductions Logical indicating whether to include dimensionality reductions (default: TRUE)
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#'
#' @return Invisibly returns the output file path
#'
#' @details
#' This function converts a Seurat object to h5ad format while preserving:
#' \itemize{
#'   \item Sparse matrix format for memory efficiency
#'   \item Cell metadata from Seurat meta.data
#'   \item Gene metadata (feature names)
#'   \item Dimensionality reductions (PCA, UMAP, tSNE, etc.)
#'   \item Multiple assays as layers in the h5ad file
#' }
#'
#' The function uses the anndata package through reticulate to create the h5ad file.
#' The specified assay and layer will be used as the main X matrix, while other
#' assays will be stored as additional layers.
#'
#' @examples
#' \dontrun{
#' # Convert Seurat object to h5ad
#' seurat_to_h5ad(seurat_obj, "output.h5ad")
#'
#' # Use counts instead of data
#' seurat_to_h5ad(seurat_obj, "output.h5ad", layer = "counts")
#'
#' # Specify different default assay
#' seurat_to_h5ad(seurat_obj, "output.h5ad", default_assay = "SCT")
#' }
#'
#' @export
seurat_to_h5ad <- function(seurat_obj, 
                          output_file, 
                          default_assay = "RNA",
                          layer = "data",
                          include_reductions = TRUE,
                          overwrite = FALSE,
                          name_conflict = c("make_unique", "error"),
                          verbose = TRUE) {
  
  # Check if required packages are available
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required but not installed")
  }
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required but not installed")
  }
  
  if (!requireNamespace("anndata", quietly = TRUE)) {
    stop("Package 'anndata' is required but not installed")
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
  
  if (verbose) cat(icb_i18n("正在将 Seurat 对象转换为 h5ad 格式...\n", "Converting Seurat object to h5ad format...\n"))

  name_conflict <- match.arg(name_conflict)
  # Overwrite behavior
  if (file.exists(output_file) && !isTRUE(overwrite)) {
    stop(icb_i18n(
      zh = paste0("输出目标已存在：", output_file, "。设置 overwrite=TRUE 以允许覆盖，或更换文件名。"),
      en = paste0("Output file exists: ", output_file, ". Set overwrite=TRUE to allow overwrite, or choose a new path.")
    ))
  }
  
  # Check if the specified assay exists
  if (!default_assay %in% names(seurat_obj@assays)) {
    stop("Assay '", default_assay, "' not found in Seurat object. Available assays: ", 
         paste(names(seurat_obj@assays), collapse = ", "))
  }
  
  # Get the main matrix (X)
  if (verbose) cat("Extracting main matrix from assay:", default_assay, "layer:", layer, "\n")
  
  if (layer == "data") {
    main_matrix <- Seurat::GetAssayData(seurat_obj, assay = default_assay, layer = "data")
    # If data slot is empty, try counts slot
    if (is.null(main_matrix) || nrow(main_matrix) == 0 || all(main_matrix == 0)) {
      if (verbose) cat("Data slot is empty, using counts slot instead...\n")
      main_matrix <- Seurat::GetAssayData(seurat_obj, assay = default_assay, layer = "counts")
      layer <- "counts"  # Update layer info for later processing
    }
  } else if (layer == "counts") {
    main_matrix <- Seurat::GetAssayData(seurat_obj, assay = default_assay, layer = "counts")
  } else {
    stop("Layer must be either 'data' or 'counts'")
  }
  
  # Check if matrix is valid
  if (is.null(main_matrix) || nrow(main_matrix) == 0) {
    stop("No valid expression data found in the specified assay and layer")
  }
  
  # Transpose matrix (Seurat: genes x cells, h5ad: cells x genes)
  main_matrix <- Matrix::t(main_matrix)
  
  if (verbose) {
    cat("Main matrix dimensions:", dim(main_matrix)[1], "cells x", dim(main_matrix)[2], "genes\n")
    cat("Matrix is sparse:", is(main_matrix, "sparseMatrix"), "\n")
  }
  
  # Get cell metadata
  if (verbose) cat(icb_i18n("提取细胞元数据...\n", "Extracting cell metadata...\n"))
  cell_metadata <- seurat_obj@meta.data
  
  # Get gene names
  gene_names <- rownames(seurat_obj@assays[[default_assay]])
  cell_names <- colnames(seurat_obj@assays[[default_assay]])
  
  # Create gene metadata
  gene_metadata <- data.frame(
    gene_ids = gene_names,
    row.names = gene_names,
    stringsAsFactors = FALSE
  )
  
  if (verbose) cat(icb_i18n("创建 AnnData 对象...\n", "Creating AnnData object...\n"))
  
  # Normalize names (conflict strategy)
  if (!is.null(colnames(main_matrix))) {
    colnames(main_matrix) <- icb_make_unique(colnames(main_matrix), strategy = name_conflict, sep = "-")
  }
  if (!is.null(rownames(main_matrix))) {
    rownames(main_matrix) <- icb_make_unique(rownames(main_matrix), strategy = name_conflict, sep = "-")
  }

  # Create AnnData object using anndata package
  adata <- anndata::AnnData(
    X = main_matrix,
    obs = cell_metadata,
    var = gene_metadata
  )
  
  # Add other assays as layers
  if (length(seurat_obj@assays) > 1) {
    if (verbose) cat(icb_i18n("添加其它 assay 到 layers...\n", "Adding additional assays as layers...\n"))
    
    for (assay_name in names(seurat_obj@assays)) {
      if (assay_name != default_assay) {
        # Add both data and counts if available
        tryCatch({
          data_matrix <- Seurat::GetAssayData(seurat_obj, assay = assay_name, layer = "data")
          if (!is.null(data_matrix) && nrow(data_matrix) > 0) {
            data_matrix <- Matrix::t(data_matrix)
            adata$layers[[paste0(assay_name, "_data")]] <- data_matrix
            if (verbose) cat(icb_i18n("已添加层:", "Added layer:"), paste0(assay_name, "_data"), "\n")
          }
        }, error = function(e) {
          if (verbose) cat("Warning: Could not add", assay_name, "data layer:", as.character(e$message), "\n")
        })
        
        tryCatch({
          counts_matrix <- Seurat::GetAssayData(seurat_obj, assay = assay_name, layer = "counts")
          if (!is.null(counts_matrix) && nrow(counts_matrix) > 0) {
            counts_matrix <- Matrix::t(counts_matrix)
            adata$layers[[paste0(assay_name, "_counts")]] <- counts_matrix
            if (verbose) cat(icb_i18n("已添加层:", "Added layer:"), paste0(assay_name, "_counts"), "\n")
          }
        }, error = function(e) {
          if (verbose) cat("Warning: Could not add", assay_name, "counts layer:", as.character(e$message), "\n")
        })
      }
    }
  }
  
  # Add the alternative layer for the default assay
  if (layer == "data") {
    # Add counts as a layer if available
    tryCatch({
      counts_matrix <- Seurat::GetAssayData(seurat_obj, assay = default_assay, layer = "counts")
      if (!is.null(counts_matrix) && nrow(counts_matrix) > 0) {
        counts_matrix <- Matrix::t(counts_matrix)
        adata$layers[["counts"]] <- counts_matrix
        if (verbose) cat(icb_i18n("已从默认 assay 添加 counts 层\n", "Added counts layer from default assay\n"))
      }
    }, error = function(e) {
      if (verbose) cat("Warning: Could not add counts layer:", as.character(e$message), "\n")
    })
  } else {
    # Add data as a layer if available
    tryCatch({
      data_matrix <- Seurat::GetAssayData(seurat_obj, assay = default_assay, layer = "data")
      if (!is.null(data_matrix) && nrow(data_matrix) > 0) {
        data_matrix <- Matrix::t(data_matrix)
        adata$layers[["logcounts"]] <- data_matrix
        if (verbose) cat(icb_i18n("已从默认 assay 添加 logcounts 层\n", "Added logcounts layer from default assay\n"))
      }
    }, error = function(e) {
      if (verbose) cat("Warning: Could not add logcounts layer:", as.character(e$message), "\n")
    })
  }
  
  # Add dimensionality reductions
  if (include_reductions && length(seurat_obj@reductions) > 0) {
    if (verbose) cat(icb_i18n("添加降维结果...\n", "Adding dimensionality reductions...\n"))
    
    for (reduction_name in names(seurat_obj@reductions)) {
      tryCatch({
        reduction_data <- seurat_obj@reductions[[reduction_name]]@cell.embeddings
        
        # Convert reduction name to h5ad convention
        h5ad_name <- paste0("X_", tolower(reduction_name))
        
        adata$obsm[[h5ad_name]] <- reduction_data
        if (verbose) cat("Added reduction:", h5ad_name, "with dimensions", dim(reduction_data)[2], "\n")
      }, error = function(e) {
        if (verbose) cat("Warning: Could not add reduction", reduction_name, ":", as.character(e$message), "\n")
      })
    }
  }
  
  # Write to h5ad file
  if (verbose) cat(icb_i18n("写入 h5ad 文件:", "Writing to h5ad file:"), output_file, "\n")
  adata$write_h5ad(output_file)
  
  if (verbose) {
    cat("✓ Conversion completed successfully!\n")
    cat("H5AD file created:", output_file, "\n")
    file_size <- file.size(output_file)
    cat("File size:", round(file_size / 1024^2, 2), "MB\n")
  }
  
  invisible(output_file)
}