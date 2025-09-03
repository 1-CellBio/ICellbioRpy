# Stereo-seq GEF File Reading Example for ICellbioRpy
# 
# This example demonstrates how to read Stereo-seq GEF files with cellborder support
# and visualize the spatial data using the ICellbioRpy package.

library(ICellbioRpy)
library(ggplot2)

# Example 1: Basic GEF file reading
# Read cellbin GEF file with cell borders
cat("Reading GEF file...\n")
stereo_data <- read_gef(
  file_path = "C04042E3.cellbin.gef",
  bin_type = "cell_bins",
  max_cells = 1000,  # Limit for demonstration
  include_cellborder = TRUE
)

# Display summary
cat("Data summary:\n")
summary(stereo_data)

# Example 2: Convert to Seurat object
cat("\nConverting to Seurat object...\n")
seurat_obj <- as.Seurat(stereo_data)

# Check basic information
cat("Seurat object created with", ncol(seurat_obj), "cells and", nrow(seurat_obj), "genes\n")
cat("Cell border information available:", !is.null(seurat_obj@misc$cell_borders), "\n")

# Example 3: Convert to SingleCellExperiment object
cat("\nConverting to SingleCellExperiment object...\n")
sce_obj <- as.SingleCellExperiment(stereo_data)

cat("SCE object created with", ncol(sce_obj), "cells and", nrow(sce_obj), "genes\n")
cat("Spatial coordinates available:", "spatial" %in% names(SingleCellExperiment::reducedDims(sce_obj)), "\n")

# Example 4: Spatial visualization
cat("\nCreating spatial plots...\n")

# Plot cells without borders
p1 <- plot_cells_with_borders(
  object = seurat_obj,
  show_borders = FALSE,
  point_size = 1.0
) + 
  ggtitle("Cells without borders") +
  theme(plot.title = element_text(hjust = 0.5))

# Plot cells with borders
p2 <- plot_cells_with_borders(
  object = seurat_obj,
  show_borders = TRUE,
  border_color = "darkblue",
  border_size = 0.2,
  point_size = 0.8
) + 
  ggtitle("Cells with borders") +
  theme(plot.title = element_text(hjust = 0.5))

# Color by cell area if available
if ("area" %in% colnames(seurat_obj@meta.data)) {
  p3 <- plot_cells_with_borders(
    object = seurat_obj,
    color_by = "area",
    show_borders = TRUE,
    border_color = "white",
    border_size = 0.1,
    point_size = 0.7
  ) + 
    ggtitle("Cells colored by area") +
    scale_color_viridis_c() +
    theme(plot.title = element_text(hjust = 0.5))
  
  cat("Created plot colored by cell area\n")
}

# Example 5: Filter by spatial region
cat("\nReading specific spatial region...\n")
stereo_region <- read_gef(
  file_path = "C04042E3.cellbin.gef",
  bin_type = "cell_bins",
  region = c(8000, 12000, 2000, 6000),  # minX, maxX, minY, maxY
  include_cellborder = TRUE
)

cat("Region data contains", nrow(stereo_region$cell_info), "cells\n")

# Plot the region
seurat_region <- as.Seurat(stereo_region)
p4 <- plot_cells_with_borders(
  object = seurat_region,
  show_borders = TRUE,
  border_color = "red",
  border_size = 0.3,
  point_size = 0.1
) + 
  ggtitle("Filtered spatial region") +
  theme(plot.title = element_text(hjust = 0.5))

# Example 6: Gene filtering
cat("\nReading with gene filtering...\n")
# Note: This requires knowing gene names in the dataset
# stereo_genes <- read_gef(
#   file_path = "C04042E3.cellbin.gef",
#   bin_type = "cell_bins",
#   gene_list = c("GAPDH", "ACTB", "RPL13A"),  # Example housekeeping genes
#   max_cells = 500
# )

cat("\nExample completed successfully!\n")
cat("Available plots: p1 (no borders), p2 (with borders)")
if (exists("p3")) cat(", p3 (colored by area)")
cat(", p4 (spatial region)\n")

# Display plots (uncomment to show)
# print(p1)
# print(p2)
# if (exists("p3")) print(p3)
# print(p4)

# Example 7: Working with cell borders directly
cat("\nExamining cell border data structure...\n")
cell_borders <- seurat_obj@misc$cell_borders

if (!is.null(cell_borders)) {
  # Check the first few cells' border data
  for (i in 1:min(3, length(cell_borders))) {
    cell_name <- names(cell_borders)[i]
    border_points <- cell_borders[[i]]
    cat("Cell", cell_name, "has", nrow(border_points), "border points\n")
    
    if (nrow(border_points) > 0) {
      cat("  Border coordinates range: X[", min(border_points$x), ",", max(border_points$x), 
          "], Y[", min(border_points$y), ",", max(border_points$y), "]\n")
    }
  }
}

cat("\nAll examples completed!\n")





###########
devtools::load_all('1CellbioRpy')
library(data.table)
library(Seurat)
library(future.apply)
seu = readRDS('~/Downloads/sempreqc.RDS')
polygon = fread('~/Downloads/S0-polygons.csv.gz', data.table = F)

processed_data <- data.frame(
  cell_id = polygon$cell,  # 使用 'cell' 作为唯一ID
  x = polygon$x_global_px, # 全局x坐标
  y = polygon$y_global_px  # 全局y坐标
)
polygon_list <- split(processed_data[, c("x", "y")], processed_data$cell_id)
polygon_list <- lapply(polygon_list, function(x){rownames(x) = 1:nrow(x);return(x)})
spatial = lapply(polygon_list, colMeans)
spatial = do.call(rbind, spatial)
colnames(spatial) = c('spatial_1','spatial_2')

spatial = as.data.frame(spatial)
spatial_sub = spatial[spatial$spatial_2 > 119121.5/1.5,]
dim(spatial_sub)

seu_sub = seu[,rownames(spatial_sub)]


custom_reduc_obj <- CreateDimReducObject(
  embeddings = as.matrix(spatial_sub),
  key = "spatial_", # 务必以 '_' 结尾
  assay = "RNA" # 指定来源 Assay
)

seu_sub[['spatial']] = custom_reduc_obj


polygon_list_sub = future_lapply(rownames(spatial_sub), function(x){polygon_list[[x]]})
names(polygon_list_sub) = rownames(spatial_sub)

seu_sub@misc$cell_borders = polygon_list_sub


# Plot cells with borders
p2 <- plot_cells_with_borders(
  object = seu_sub,
  show_borders = TRUE,
  border_color = "darkblue",
  border_size = 0.2,
  point_size = 0.01
)

configure_python_env(conda_env = '1cellbio')
adata = anndata::AnnData(X = t(seu_sub@assays$RNA@counts),obs = seu_sub@meta.data, var = seu_sub@assays$RNA@meta.features)

adata$obsm[['spatial']] = as.matrix(seu_sub@reductions$spatial@cell.embeddings)

object = seu_sub@misc
# Optimize: use single arrays for all polygon data (efficient HDF5 storage)
all_vertices <- list()  # All vertices as [x, y] pairs
polygon_starts <- integer(0)  # Starting index for each polygon
polygon_lengths <- integer(0)  # Number of vertices per polygon
border_count <- 0
current_vertex_count <- 0
verbose = T

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

adata$write_h5ad('cosMX_colon_ST_sub.h5ad')
