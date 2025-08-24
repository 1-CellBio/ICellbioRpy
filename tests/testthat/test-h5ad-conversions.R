context("H5AD conversion functions")

# Helper function to create a minimal test h5ad file
create_test_h5ad <- function() {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("anndata")
  
  # Create minimal test data
  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 1, 3),
    j = c(1, 1, 2, 3),
    x = c(1, 2, 3, 4),
    dims = c(3, 3)
  )
  
  # Create cell metadata
  obs <- data.frame(
    cell_type = c("TypeA", "TypeB", "TypeA"),
    total_counts = c(4, 2, 7),
    row.names = c("Cell1", "Cell2", "Cell3")
  )
  
  # Create gene metadata
  var <- data.frame(
    gene_name = c("Gene1", "Gene2", "Gene3"),
    highly_variable = c(TRUE, FALSE, TRUE),
    row.names = c("Gene1", "Gene2", "Gene3")
  )
  
  # Try to configure Python environment
  tryCatch({
    configure_python_env(conda_env = "1cellbio", verbose = FALSE)
  }, error = function(e) {
    skip("Python environment not available")
  })
  
  # Create AnnData object
  adata <- anndata::AnnData(
    X = Matrix::t(counts),  # transpose for cells x genes
    obs = obs,
    var = var
  )
  
  # Create temporary file
  temp_file <- tempfile(fileext = ".h5ad")
  adata$write_h5ad(temp_file)
  
  return(temp_file)
}

test_that("h5ad_to_sce conversion works", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("S4Vectors")
  
  temp_h5ad <- create_test_h5ad()
  on.exit(unlink(temp_h5ad))
  
  # Test basic conversion
  sce <- h5ad_to_sce(temp_h5ad, verbose = FALSE)
  
  # Check object class
  expect_true(inherits(sce, "SingleCellExperiment"))
  
  # Check dimensions
  expect_equal(ncol(sce), 3)  # 3 cells
  expect_equal(nrow(sce), 3)  # 3 genes
  
  # Check metadata
  expect_true("cell_type" %in% colnames(SummarizedExperiment::colData(sce)))
  expect_true("gene_name" %in% colnames(SummarizedExperiment::rowData(sce)))
})

test_that("h5ad_to_seurat conversion works", {
  skip_if_not_installed("Seurat")
  
  temp_h5ad <- create_test_h5ad()
  on.exit(unlink(temp_h5ad))
  
  # Test basic conversion
  seurat_obj <- h5ad_to_seurat(temp_h5ad, verbose = FALSE)
  
  # Check object class
  expect_true(inherits(seurat_obj, "Seurat"))
  
  # Check dimensions
  expect_equal(ncol(seurat_obj), 3)  # 3 cells
  expect_equal(nrow(seurat_obj), 3)  # 3 genes
  
  # Check metadata
  expect_true("cell_type" %in% colnames(seurat_obj@meta.data))
})

test_that("seurat_to_h5ad conversion works", {
  skip_if_not_installed("Seurat")
  
  # Create minimal Seurat object
  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 1, 3),
    j = c(1, 1, 2, 3),
    x = c(1, 2, 3, 4),
    dims = c(3, 3),
    dimnames = list(
      c("Gene1", "Gene2", "Gene3"),
      c("Cell1", "Cell2", "Cell3")
    )
  )
  
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts,
    meta.data = data.frame(
      cell_type = c("TypeA", "TypeB", "TypeA"),
      row.names = c("Cell1", "Cell2", "Cell3")
    )
  )
  
  # Test conversion
  temp_h5ad <- tempfile(fileext = ".h5ad")
  on.exit(unlink(temp_h5ad))
  
  # Try to configure Python environment
  tryCatch({
    configure_python_env(conda_env = "1cellbio", verbose = FALSE)
  }, error = function(e) {
    skip("Python environment not available")
  })
  
  expect_no_error({
    seurat_to_h5ad(seurat_obj, temp_h5ad, verbose = FALSE)
  })
  
  # Check that file was created
  expect_true(file.exists(temp_h5ad))
  expect_true(file.size(temp_h5ad) > 0)
})

test_that("name_conflict parameter works", {
  skip_if_not_installed("SingleCellExperiment")
  
  temp_h5ad <- create_test_h5ad()
  on.exit(unlink(temp_h5ad))
  
  # Test make_unique strategy (should work without error)
  expect_no_error({
    sce1 <- h5ad_to_sce(temp_h5ad, name_conflict = "make_unique", verbose = FALSE)
  })
  
  # Test error strategy (should work with unique names)
  expect_no_error({
    sce2 <- h5ad_to_sce(temp_h5ad, name_conflict = "error", verbose = FALSE)
  })
})
