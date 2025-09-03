#' Read Stereo-seq GEF files with cellborder support
#'
#' This function reads Stereo-seq GEF (Gene Expression File) format and creates
#' a StereoData object that can be converted to Seurat or SingleCellExperiment.
#' Specifically designed to handle cellBin data with cell border information.
#'
#' @param file_path Path to the GEF file
#' @param bin_type Type of binning: "bins" for square bins or "cell_bins" for cell segmentation
#' @param bin_size Bin size for square bins (only used when bin_type = "bins")
#' @param gene_list Optional vector of genes to filter
#' @param region Optional vector c(minX, maxX, minY, maxY) to filter spatial region
#' @param max_cells Maximum number of cells to read (for memory management)
#' @param include_cellborder Whether to include cell border information (default TRUE)
#' @return A StereoData object
#' @export
#' @examples
#' \dontrun{
#' # Read cellbin GEF file with cell borders
#' stereo_data <- read_gef("data.cellbin.gef", bin_type = "cell_bins")
#' 
#' # Convert to Seurat
#' seurat_obj <- as.Seurat(stereo_data)
#' 
#' # Convert to SingleCellExperiment
#' sce_obj <- as.SingleCellExperiment.StereoData(stereo_data)
#' }
read_gef <- function(file_path, 
                     bin_type = c("cell_bins", "bins"),
                     bin_size = 50,
                     gene_list = NULL,
                     region = NULL,
                     max_cells = NULL,
                     include_cellborder = TRUE) {
  
  bin_type <- match.arg(bin_type)
  
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("GEF file not found: ", file_path)
  }
  
  # Load required libraries
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required for reading GEF files. Please install it with: install.packages('hdf5r')")
  }
  
  message("Opening GEF file: ", file_path)
  h5file <- hdf5r::H5File$new(file_path, mode = "r")
  on.exit(h5file$close())
  
  if (bin_type == "cell_bins") {
    return(read_gef_cellbins(h5file, gene_list, region, max_cells, include_cellborder))
  } else {
    return(read_gef_bins(h5file, bin_size, gene_list, region))
  }
}

#' Read GEF file with cell bin data
#' @keywords internal
read_gef_cellbins <- function(h5file, gene_list = NULL, region = NULL, max_cells = NULL, include_cellborder = TRUE) {
  
  message("Reading cell bin data...")
  
  # Check if cellBin group exists
  if (!"cellBin" %in% h5file$ls()$name) {
    stop("cellBin group not found in GEF file. This may not be a cellbin GEF file.")
  }
  
  cellbin_group <- h5file[["cellBin"]]
  
  # Read basic information
  message("Reading cell information...")
  cell_info <- read_gef_cell_info(cellbin_group, max_cells, region)
  
  message("Reading gene information...")
  gene_info <- read_gef_gene_info(cellbin_group, gene_list)
  
  message("Reading expression matrix...")
  expr_matrix <- read_gef_expression_matrix(cellbin_group, cell_info, gene_info)
  
  # Read cell borders if requested
  cell_borders <- NULL
  if (include_cellborder && "cellBorder" %in% cellbin_group$ls()$name) {
    message("Reading cell border information...")
    cell_borders <- read_gef_cell_borders(cellbin_group, cell_info)
  }
  
  # Create StereoData object
  stereo_data <- list(
    expression = expr_matrix,
    cell_info = cell_info,
    gene_info = gene_info,
    cell_borders = cell_borders,
    bin_type = "cell_bins",
    file_path = h5file$filename
  )
  
  class(stereo_data) <- "StereoData"
  
  message("Successfully loaded ", nrow(cell_info), " cells and ", nrow(gene_info), " genes")
  if (!is.null(cell_borders)) {
    message("Cell border information included for spatial visualization")
  }
  
  return(stereo_data)
}

#' Read GEF file with square bin data
#' @keywords internal
read_gef_bins <- function(h5file, bin_size, gene_list = NULL, region = NULL) {
  stop("Square bin reading not yet implemented. Please use bin_type = 'cell_bins'")
}

#' Read cell information from GEF cellBin group
#' @keywords internal
read_gef_cell_info <- function(cellbin_group, max_cells = NULL, region = NULL) {
  
  if (!"cell" %in% cellbin_group$ls()$name) {
    stop("cell dataset not found in cellBin group")
  }
  
  cell_dataset <- cellbin_group[["cell"]]
  cell_data <- cell_dataset$read()
  
  # Convert to data frame
  cell_df <- data.frame(
    cell_id = cell_data$id,
    x = cell_data$x,
    y = cell_data$y,
    offset = cell_data$offset,
    gene_count = cell_data$geneCount,
    exp_count = cell_data$expCount,
    dnb_count = cell_data$dnbCount,
    area = cell_data$area,
    cell_type_id = cell_data$cellTypeID,
    cluster_id = cell_data$clusterID,
    stringsAsFactors = FALSE
  )
  
  # Apply region filter if specified
  if (!is.null(region) && length(region) == 4) {
    min_x <- region[1]
    max_x <- region[2]
    min_y <- region[3]
    max_y <- region[4]
    
    region_filter <- cell_df$x >= min_x & cell_df$x <= max_x & 
                     cell_df$y >= min_y & cell_df$y <= max_y
    cell_df <- cell_df[region_filter, ]
    message("Filtered to ", nrow(cell_df), " cells in specified region")
  }
  
  # Apply max_cells limit if specified
  if (!is.null(max_cells) && nrow(cell_df) > max_cells) {
    cell_df <- cell_df[1:max_cells, ]
    message("Limited to first ", max_cells, " cells")
  }
  
  rownames(cell_df) <- paste0("Cell_", cell_df$cell_id)
  
  return(cell_df)
}

#' Read gene information from GEF cellBin group
#' @keywords internal
read_gef_gene_info <- function(cellbin_group, gene_list = NULL) {
  
  if (!"gene" %in% cellbin_group$ls()$name) {
    stop("gene dataset not found in cellBin group")
  }
  
  gene_dataset <- cellbin_group[["gene"]]
  gene_data <- gene_dataset$read()
  
  # Convert to data frame
  gene_df <- data.frame(
    gene_id = gene_data$geneID,
    gene_name = gene_data$geneName,
    offset = gene_data$offset,
    cell_count = gene_data$cellCount,
    exp_count = gene_data$expCount,
    max_mid_count = gene_data$maxMIDcount,
    stringsAsFactors = FALSE
  )
  
  # Apply gene filter if specified
  if (!is.null(gene_list)) {
    gene_filter <- gene_df$gene_name %in% gene_list | gene_df$gene_id %in% gene_list
    gene_df <- gene_df[gene_filter, ]
    message("Filtered to ", nrow(gene_df), " genes from gene list")
  }
  
  # Use gene names as row names, fallback to gene IDs if names are empty
  gene_names <- ifelse(nchar(gene_df$gene_name) > 0, gene_df$gene_name, gene_df$gene_id)
  rownames(gene_df) <- make.unique(gene_names)
  
  return(gene_df)
}

#' Read expression matrix from GEF cellBin group (fast vectorized construction)
#' @keywords internal
read_gef_expression_matrix <- function(cellbin_group, cell_info, gene_info) {
  
  # Read expression data
  if (!"cellExp" %in% cellbin_group$ls()$name) {
    stop("cellExp dataset not found in cellBin group")
  }
  
  cell_exp_dataset <- cellbin_group[["cellExp"]]
  cell_exp_data <- cell_exp_dataset$read()
  
  n_cells <- nrow(cell_info)
  n_genes <- nrow(gene_info)
  
  # Extract all data vectors at once for maximum efficiency
  all_gene_indices_0based <- cell_exp_data$geneID
  all_counts <- cell_exp_data$count
  
  # Convert to 1-based gene indices
  all_gene_indices <- all_gene_indices_0based + 1
  
  # Filter for valid genes
  valid_gene_mask <- all_gene_indices > 0 & all_gene_indices <= n_genes
  
  # Apply filter to reduce data size early
  filtered_gene_indices <- all_gene_indices[valid_gene_mask]
  filtered_counts <- all_counts[valid_gene_mask]
  
  # Build cell indices using vectorized approach with offsets
  # Create a lookup table for which cell each expression entry belongs to
  n_total_entries <- length(all_gene_indices_0based)
  cell_lookup <- integer(n_total_entries)
  
  # Use vectorized assignment based on cell offsets and counts
  for (i in seq_len(n_cells)) {
    if (cell_info$exp_count[i] > 0) {
      start_pos <- cell_info$offset[i] + 1  # Convert to 1-based
      end_pos <- min(cell_info$offset[i] + cell_info$exp_count[i], n_total_entries)
      
      if (start_pos <= end_pos) {
        cell_lookup[start_pos:end_pos] <- i
      }
    }
  }
  
  # Apply the same filtering to cell lookup as we did to gene data
  filtered_cell_indices <- cell_lookup[valid_gene_mask]
  
  # Remove any entries that couldn't be assigned to a cell (safety check)
  assigned_mask <- filtered_cell_indices > 0
  final_gene_indices <- filtered_gene_indices[assigned_mask]
  final_cell_indices <- filtered_cell_indices[assigned_mask]
  final_counts <- filtered_counts[assigned_mask]
  
  # Create sparse matrix directly from the final vectors
  if (length(final_gene_indices) > 0) {
    expr_matrix <- Matrix::sparseMatrix(
      i = final_gene_indices,
      j = final_cell_indices,
      x = final_counts,
      dims = c(n_genes, n_cells),
      dimnames = list(rownames(gene_info), rownames(cell_info))
    )
  } else {
    # Create empty sparse matrix
    expr_matrix <- Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      x = numeric(0),
      dims = c(n_genes, n_cells),
      dimnames = list(rownames(gene_info), rownames(cell_info))
    )
  }
  
  return(expr_matrix)
}

#' Read cell border information from GEF cellBin group
#' @keywords internal
read_gef_cell_borders <- function(cellbin_group, cell_info) {
  
  if (!"cellBorder" %in% cellbin_group$ls()$name) {
    warning("cellBorder dataset not found in cellBin group")
    return(NULL)
  }
  
  cellborder_dataset <- cellbin_group[["cellBorder"]]
  
  # Get dataset dimensions - format is (2, max_vertices, n_cells) in HDF5
  dims <- cellborder_dataset$dims
  cat(sprintf("Border dataset dimensions: %s\n", paste(dims, collapse=" x ")))
  
  # Read the full 3D array
  border_data <- cellborder_dataset[, , ]
  
  # The data format is (2, max_vertices, n_cells) where:
  # - first dimension: 0=x coordinates, 1=y coordinates  
  # - second dimension: vertex index (max 32)
  # - third dimension: cell index
  # 32767 is used as a sentinel value for missing vertices
  
  max_vertices <- dims[2]
  n_cells <- dims[3]
  
  # Create a list of border polygons for each cell
  cell_borders <- vector("list", nrow(cell_info))
  names(cell_borders) <- rownames(cell_info)
  
  for (i in seq_len(min(n_cells, nrow(cell_info)))) {
    # Extract coordinates for this cell - note the corrected indexing
    x_coords <- border_data[1, , i]  # x coordinates are in first dimension
    y_coords <- border_data[2, , i]  # y coordinates are in second dimension
    
    # Remove sentinel values (32767)
    valid_coords <- x_coords != 32767 & y_coords != 32767
    
    if (any(valid_coords)) {
      # Get cell center coordinates
      cell_center_x <- cell_info$x[i]
      cell_center_y <- cell_info$y[i]
      
      # Convert relative coordinates to absolute coordinates
      abs_x <- cell_center_x + x_coords[valid_coords]
      abs_y <- cell_center_y + y_coords[valid_coords]
      
      cell_borders[[i]] <- data.frame(
        x = abs_x,
        y = abs_y
      )
    } else {
      cell_borders[[i]] <- NULL  # Changed from empty data.frame to NULL
    }
  }
  
  return(cell_borders)
}

#' Convert StereoData to Seurat object
#' @param object A StereoData object
#' @param ... Additional arguments
#' @return A Seurat object
#' @method as.Seurat StereoData
#' @export
as.Seurat.StereoData <- function(object, ...) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required. Please install it.")
  }
  
  # Create Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = object$expression,
    meta.data = object$cell_info
  )
  
  # Add spatial coordinates
  spatial_coords <- object$cell_info[, c("x", "y")]
  colnames(spatial_coords) <- c("spatial_1", "spatial_2")
  rownames(spatial_coords) <- rownames(object$cell_info)
  
  # Add cell borders if available
  if (!is.null(object$cell_borders)) {
    # Store cell borders in misc slot
    seurat_obj@misc$cell_borders <- object$cell_borders
    message("Cell border information stored in @misc$cell_borders")
  }
  
  # Add spatial coordinates to reductions
  seurat_obj[["spatial"]] <- Seurat::CreateDimReducObject(
    embeddings = as.matrix(spatial_coords),
    key = "spatial_",
    assay = "RNA"
  )
  
  return(seurat_obj)
}

#' Convert StereoData to SingleCellExperiment object
#' @param object A StereoData object
#' @param ... Additional arguments
#' @return A SingleCellExperiment object
#' @method as.SingleCellExperiment StereoData
#' @export
as.SingleCellExperiment.StereoData <- function(object, ...) {
  
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required. Please install it.")
  }
  
  # Create SingleCellExperiment object
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = object$expression),
    colData = object$cell_info,
    rowData = object$gene_info
  )
  
  # Add spatial coordinates
  spatial_coords <- object$cell_info[, c("x", "y")]
  colnames(spatial_coords) <- c("x", "y")
  SingleCellExperiment::reducedDim(sce, "spatial") <- as.matrix(spatial_coords)
  
  # Add cell borders if available
  if (!is.null(object$cell_borders)) {
    # Store in metadata
    S4Vectors::metadata(sce)$cell_borders <- object$cell_borders
    message("Cell border information stored in metadata")
  }
  
  return(sce)
}

#' Plot cells with borders
#' @param object A StereoData object or Seurat object with cell borders
#' @param color_by Variable to color cells by (from cell metadata)
#' @param show_borders Whether to show cell borders
#' @param border_color Color for cell borders
#' @param border_size Line width for cell borders
#' @param point_size Size of cell center points
#' @return A ggplot object
#' @export
plot_cells_with_borders <- function(object, 
                                   color_by = NULL,
                                   show_borders = TRUE,
                                   border_color = "black",
                                   border_size = 0.1,
                                   point_size = 0.5) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }
  
  # Extract data based on object type
  if (inherits(object, "StereoData")) {
    cell_info <- object$cell_info
    cell_borders <- object$cell_borders
  } else if (inherits(object, "Seurat")) {
    cell_info <- object@meta.data
    cell_info$x <- object@reductions$spatial@cell.embeddings[, 1]
    cell_info$y <- object@reductions$spatial@cell.embeddings[, 2]
    cell_borders <- object@misc$cell_borders
  } else {
    stop("Object must be StereoData or Seurat with spatial information")
  }
  
  # Create base plot
  p <- ggplot2::ggplot(cell_info, ggplot2::aes(x = x, y = y))
  
  # Add cell borders if available and requested
  if (show_borders && !is.null(cell_borders)) {
    # Prepare border data for ggplot
    border_df_list <- list()
    border_count <- 0
    
    for (i in seq_along(cell_borders)) {
      if (!is.null(cell_borders[[i]]) && is.data.frame(cell_borders[[i]]) && nrow(cell_borders[[i]]) > 0) {
        border_count <- border_count + 1
        border_df_list[[border_count]] <- data.frame(
          x = cell_borders[[i]]$x,
          y = cell_borders[[i]]$y,
          cell_id = names(cell_borders)[i],
          stringsAsFactors = FALSE
        )
      }
    }
    
    if (length(border_df_list) > 0) {
      border_df <- do.call(rbind, border_df_list)
      
      # Use geom_polygon for filled borders or geom_path for outline only
      p <- p + ggplot2::geom_polygon(
        data = border_df,
        ggplot2::aes(x = x, y = y, group = cell_id),
        fill = NA,
        color = border_color,
        linewidth = border_size,  # Updated from size to linewidth
        alpha = 0.8
      )
      
      message(sprintf("Added borders for %d cells to plot", length(unique(border_df$cell_id))))
    } else {
      message("No valid cell borders found for plotting")
    }
  }
  
  # Add cell center points
  if (!is.null(color_by) && color_by %in% colnames(cell_info)) {
    p <- p + ggplot2::geom_point(ggplot2::aes_string(color = color_by), size = point_size)
  } else {
    p <- p + ggplot2::geom_point(size = point_size)
  }
  
  # Style the plot
  p <- p + 
    ggplot2::coord_equal() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Spatial Cell Plot",
      x = "X coordinate",
      y = "Y coordinate"
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
  
  return(p)
}

#' Get summary information about StereoData object
#' @param object A StereoData object
#' @export
summary.StereoData <- function(object) {
  cat("StereoData Object\n")
  cat("=================\n")
  cat("File:", object$file_path, "\n")
  cat("Bin type:", object$bin_type, "\n")
  cat("Number of cells:", nrow(object$cell_info), "\n")
  cat("Number of genes:", nrow(object$gene_info), "\n")
  cat("Cell borders included:", !is.null(object$cell_borders), "\n")
  
  if (!is.null(object$cell_borders)) {
    non_empty_borders <- sum(sapply(object$cell_borders, nrow) > 0)
    cat("Cells with border data:", non_empty_borders, "\n")
  }
  
  cat("\nExpression matrix:\n")
  cat("  Dimensions:", nrow(object$expression), "x", ncol(object$expression), "\n")
  cat("  Non-zero entries:", Matrix::nnzero(object$expression), "\n")
  cat("  Sparsity:", round((1 - Matrix::nnzero(object$expression) / length(object$expression)) * 100, 2), "%\n")
  
  cat("\nSpatial range:\n")
  cat("  X: [", min(object$cell_info$x), ",", max(object$cell_info$x), "]\n")
  cat("  Y: [", min(object$cell_info$y), ",", max(object$cell_info$y), "]\n")
}

#' Convert StereoData to H5AD format
#'
#' This function converts a StereoData object (from read_gef) to an H5AD file,
#' which can be used in Python with scanpy and other tools. Cell borders are
#' stored in the uns section of the AnnData object.
#'
#' @param object A StereoData object from read_gef()
#' @param output_file Output H5AD file path
#' @param layer Which layer to use as main X matrix ("counts" or "data")
#' @param include_spatial Whether to include spatial coordinates (default TRUE)
#' @param overwrite Logical, whether to overwrite existing file (default FALSE)
#' @param name_conflict Strategy for handling duplicate names ("make_unique" or "error")
#' @param verbose Whether to print progress messages (default TRUE)
#'
#' @return Invisibly returns the output file path
#'
#' @details
#' The function automatically configures a suitable Python environment with 
#' anndata support. The resulting H5AD file contains:
#' - Expression matrix in adata.X (scipy sparse format)
#' - Cell metadata in adata.obs
#' - Gene metadata in adata.var
#' - Spatial coordinates in adata.obsm['spatial'] (if include_spatial = TRUE)
#' - Cell borders in adata.uns['cell_borders'] as optimized deck.gl format (if available):
#'   - vertices: all vertices as (n_total_vertices, 2) matrix [[x,y], [x,y], ...]
#'   - polygon_starts: starting index for each polygon in vertices array
#'   - polygon_lengths: number of vertices per polygon
#' - Conversion metadata in adata.uns['stereo_metadata']
#'
#' @examples
#' \dontrun{
#' # Read GEF file
#' stereo_data <- read_gef("sample.cellbin.gef")
#' 
#' # Convert to H5AD
#' stereo_to_h5ad(stereo_data, "output.h5ad")
#' 
#' # With custom settings
#' stereo_to_h5ad(stereo_data, "output.h5ad", 
#'                layer = "counts", 
#'                overwrite = TRUE)
#' }
#'
#' @export
stereo_to_h5ad <- function(object, 
                           output_file,
                           layer = c("counts", "data"),
                           include_spatial = TRUE,
                           overwrite = FALSE,
                           name_conflict = c("make_unique", "error"),
                           verbose = TRUE) {
  
  # Check input object
  if (!inherits(object, "StereoData")) {
    stop("Object must be a StereoData object from read_gef()")
  }
  
  # Check required packages
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required for H5AD conversion. Please install it.")
  }
  
  if (!requireNamespace("anndata", quietly = TRUE)) {
    stop("Package 'anndata' is required for H5AD conversion. Please install it.")
  }
  
  layer <- match.arg(layer)
  name_conflict <- match.arg(name_conflict)
  
  # Check overwrite
  if (file.exists(output_file) && !overwrite) {
    stop(icb_i18n(
      zh = paste0("输出文件已存在: ", output_file, "。设置 overwrite=TRUE 以允许覆盖。"),
      en = paste0("Output file exists: ", output_file, ". Set overwrite=TRUE to allow overwrite.")
    ))
  }
  
  # Configure Python environment
  if (verbose) cat(icb_i18n("配置Python环境...\n", "Configuring Python environment...\n"))
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
  
  if (verbose) cat(icb_i18n("正在转换StereoData到H5AD格式...\n", "Converting StereoData to H5AD format...\n"))
  
  # Prepare expression matrix
  expr_matrix <- object$expression
  
  # Ensure matrix is sparse for memory efficiency
  if (!inherits(expr_matrix, "sparseMatrix")) {
    if (verbose) cat("Converting to sparse matrix...\n")
    expr_matrix <- Matrix::Matrix(expr_matrix, sparse = TRUE)
  }
  
  # Transpose matrix (StereoData: genes x cells, H5AD: cells x genes)
  expr_matrix <- Matrix::t(expr_matrix)
  
  if (verbose) {
    cat("Expression matrix dimensions:", dim(expr_matrix)[1], "cells x", dim(expr_matrix)[2], "genes\n")
    cat("Matrix sparsity:", round((1 - Matrix::nnzero(expr_matrix) / length(expr_matrix)) * 100, 2), "%\n")
  }
  
  # Prepare cell metadata
  cell_metadata <- object$cell_info
  
  # Prepare gene metadata
  gene_metadata <- object$gene_info
  
  # Ensure unique names
  if (verbose) cat("Ensuring unique cell and gene names...\n")
  cell_names <- icb_make_unique(rownames(cell_metadata), strategy = name_conflict, sep = "-")
  gene_names <- icb_make_unique(rownames(gene_metadata), strategy = name_conflict, sep = "-")
  
  # Update names
  rownames(cell_metadata) <- cell_names
  rownames(gene_metadata) <- gene_names
  rownames(expr_matrix) <- cell_names
  colnames(expr_matrix) <- gene_names
  
  if (verbose) cat(icb_i18n("创建AnnData对象...\n", "Creating AnnData object...\n"))
  
  # Create AnnData object
  adata <- anndata::AnnData(
    X = expr_matrix,
    obs = cell_metadata,
    var = gene_metadata
  )
  
  # Add spatial coordinates if requested
  if (include_spatial && all(c("x", "y") %in% colnames(cell_metadata))) {
    if (verbose) cat("Adding spatial coordinates to obsm...\n")
    spatial_coords <- as.matrix(cell_metadata[, c("x", "y")])
    rownames(spatial_coords) <- cell_names
    adata$obsm[["spatial"]] <- spatial_coords
  }
  
  # Add cell borders to uns if available
  if (!is.null(object$cell_borders)) {
    if (verbose) cat("Adding cell border information to uns...\n")
    
    # Optimize: use single arrays for all polygon data (efficient HDF5 storage)
    all_vertices <- list()  # All vertices as [x, y] pairs
    polygon_starts <- integer(0)  # Starting index for each polygon
    polygon_lengths <- integer(0)  # Number of vertices per polygon
    border_count <- 0
    current_vertex_count <- 0
    
    for (i in seq_along(object$cell_borders)) {
      border_data <- object$cell_borders[[i]]
      
      if (!is.null(border_data) && is.data.frame(border_data) && nrow(border_data) > 0) {
        # Record polygon start position (0-based for Python)
        polygon_starts <- c(polygon_starts, current_vertex_count)
        
        # Record number of vertices for this polygon
        n_vertices <- nrow(border_data)
        polygon_lengths <- c(polygon_lengths, n_vertices)
        
        # Add vertices to the global list
        for (j in seq_len(n_vertices)) {
          all_vertices[[length(all_vertices) + 1]] <- c(border_data$x[j], border_data$y[j])
        }
        
        current_vertex_count <- current_vertex_count + n_vertices
        border_count <- border_count + 1
      } else {
        # For cells without borders, record empty polygon
        polygon_starts <- c(polygon_starts, current_vertex_count)
        polygon_lengths <- c(polygon_lengths, 0L)
      }
    }
    
    if (border_count > 0) {
      # Convert vertices list to matrix for efficient storage
      vertices_matrix <- do.call(rbind, all_vertices)
      
      # Create optimized deck.gl format with single arrays
      cell_borders_data <- list(
        vertices = vertices_matrix,           # All vertices as (n_total_vertices, 2) matrix
        polygon_starts = as.integer(polygon_starts),  # Start index for each polygon
        polygon_lengths = as.integer(polygon_lengths) # Number of vertices per polygon
      )
      
      adata$uns[["cell_borders"]] <- cell_borders_data
      if (verbose) {
        cat(sprintf("Stored borders for %d cells in optimized format:\n", border_count))
        cat(sprintf("  - Total vertices: %d\n", nrow(vertices_matrix)))
        cat(sprintf("  - Polygons: %d\n", length(polygon_starts)))
        cat(sprintf("  - Storage: 3 arrays instead of %d datasets\n", length(object$cell_borders)))
      }
    }
  }
  
  # Add stereo metadata to uns
  stereo_metadata <- list(
    bin_type = object$bin_type,
    source_file = basename(object$file_path),
    layer_used = layer,
    n_cells = nrow(cell_metadata),
    n_genes = nrow(gene_metadata),
    n_cells_with_borders = if (!is.null(object$cell_borders)) {
      sum(sapply(object$cell_borders, function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0))
    } else {
      0L
    },
    spatial_range = if (include_spatial && all(c("x", "y") %in% colnames(cell_metadata))) {
      list(
        x_min = min(cell_metadata$x),
        x_max = max(cell_metadata$x),
        y_min = min(cell_metadata$y),
        y_max = max(cell_metadata$y)
      )
    } else {
      NULL
    },
    conversion_time = as.character(Sys.time()),
    software_version = "ICellbioRpy"
  )
  
  adata$uns[["stereo_metadata"]] <- stereo_metadata
  
  # Write to H5AD file
  if (verbose) cat(icb_i18n("写入H5AD文件: ", "Writing H5AD file: "), output_file, "\n")
  adata$write_h5ad(output_file)
  
  if (verbose) {
    cat("✓ Conversion completed successfully!\n")
    cat("H5AD file created:", output_file, "\n")
    file_size <- file.size(output_file)
    cat("File size:", round(file_size / 1024^2, 2), "MB\n")
    
    # Summary
    cat("\nH5AD file contents:\n")
    cat(sprintf("  - Expression matrix: %d cells × %d genes\n", 
                nrow(cell_metadata), nrow(gene_metadata)))
    cat("  - Cell metadata in adata.obs\n")
    cat("  - Gene metadata in adata.var\n")
    if (include_spatial && all(c("x", "y") %in% colnames(cell_metadata))) {
      cat("  - Spatial coordinates in adata.obsm['spatial']\n")
    }
    if (!is.null(object$cell_borders)) {
      border_count <- sum(sapply(object$cell_borders, function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0))
      cat(sprintf("  - Cell borders for %d cells in adata.uns['cell_borders']\n", border_count))
    }
    cat("  - Metadata in adata.uns['stereo_metadata']\n")
  }
  
  invisible(output_file)
}

#' Convert GEF file directly to H5AD format
#'
#' This function reads a Stereo-seq GEF file and directly converts it to H5AD format
#' without creating intermediate R objects. This is memory efficient for large datasets.
#'
#' @param gef_file Path to the GEF file
#' @param h5ad_file Output H5AD file path
#' @param bin_type Type of binning ("cell_bins" or "bins", default "cell_bins")
#' @param layer Which layer to use as main X matrix ("counts" or "data", default "counts")
#' @param max_cells Maximum number of cells to read (for memory management)
#' @param region Spatial region filter c(minX, maxX, minY, maxY)
#' @param gene_list Optional vector of genes to filter
#' @param include_cellborder Whether to include cell borders (default TRUE)
#' @param include_spatial Whether to include spatial coordinates (default TRUE)
#' @param overwrite Whether to overwrite existing H5AD file (default FALSE)
#' @param name_conflict Strategy for handling duplicate names ("make_unique" or "error")
#' @param verbose Whether to print progress messages (default TRUE)
#'
#' @return Invisibly returns the output file path
#'
#' @examples
#' \dontrun{
#' # Basic conversion
#' gef_to_h5ad("sample.cellbin.gef", "output.h5ad")
#' 
#' # With filtering and custom settings
#' gef_to_h5ad(
#'   gef_file = "large_sample.cellbin.gef",
#'   h5ad_file = "filtered.h5ad", 
#'   max_cells = 5000,
#'   region = c(1000, 5000, 1000, 5000),
#'   overwrite = TRUE
#' )
#' }
#'
#' @export
gef_to_h5ad <- function(gef_file, 
                        h5ad_file,
                        bin_type = c("cell_bins", "bins"),
                        layer = c("counts", "data"),
                        max_cells = NULL,
                        region = NULL,
                        gene_list = NULL,
                        include_cellborder = TRUE,
                        include_spatial = TRUE,
                        overwrite = FALSE,
                        name_conflict = c("make_unique", "error"),
                        verbose = TRUE) {
  
  bin_type <- match.arg(bin_type)
  layer <- match.arg(layer)
  name_conflict <- match.arg(name_conflict)
  
  if (verbose) cat(icb_i18n("开始GEF到H5AD直接转换...\n", "Starting GEF to H5AD direct conversion...\n"))
  
  # Check if GEF file exists
  if (!file.exists(gef_file)) {
    stop("GEF file not found: ", gef_file)
  }
  
  # Check if output file exists
  if (file.exists(h5ad_file) && !overwrite) {
    stop(icb_i18n(
      zh = paste0("输出文件已存在: ", h5ad_file, "。设置 overwrite=TRUE 以允许覆盖。"),
      en = paste0("Output file exists: ", h5ad_file, ". Set overwrite=TRUE to allow overwrite.")
    ))
  }
  
  # Configure Python environment first
  if (verbose) cat(icb_i18n("配置Python环境...\n", "Configuring Python environment...\n"))
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
  
  # Read GEF file
  if (verbose) cat(icb_i18n("读取GEF文件...\n", "Reading GEF file...\n"))
  stereo_data <- read_gef(
    file_path = gef_file,
    bin_type = bin_type,
    max_cells = max_cells,
    region = region,
    gene_list = gene_list,
    include_cellborder = include_cellborder
  )
  
  # Convert to H5AD
  if (verbose) cat(icb_i18n("转换为H5AD格式...\n", "Converting to H5AD format...\n"))
  stereo_to_h5ad(
    object = stereo_data,
    output_file = h5ad_file,
    layer = layer,
    include_spatial = include_spatial,
    overwrite = overwrite,
    name_conflict = name_conflict,
    verbose = verbose
  )
  
  invisible(h5ad_file)
}
