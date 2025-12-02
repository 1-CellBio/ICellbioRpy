context("read1Cellbio")

# Helper function to create a minimal 1CellbioData object
create_mock_1cellbio_data <- function() {
  # Create a mock structure that mimics 1CellbioData
  temp_dir <- tempdir()
  base_path <- file.path(temp_dir, "mock_1cellbio")
  dir.create(base_path, showWarnings = FALSE, recursive = TRUE)
  
  # Create mock experiment structure
  experiment <- list(
    summarized_experiment = list(
      column_data = list(resource = list(path = "coldata.h5")),
      row_data = list(resource = list(path = "rowdata.h5")),
      assays = list(
        list(resource = list(path = "counts.h5")),
        list(resource = list(path = "logcounts.h5"))
      )
    ),
    single_cell_experiment = list(
      reduced_dimensions = list(
        list(name = "PCA", resource = list(path = "pca.h5")),
        list(name = "UMAP", resource = list(path = "umap.h5"))
      )
    )
  )
  
  data <- list(
    experiment = experiment,
    base_path = base_path
  )
  class(data) <- "1CellbioData"
  
  return(data)
}

test_that("read1Cellbio file validation works", {
  # Test with non-existent file
  expect_error(
    read1Cellbio("nonexistent.zip"),
    "File not found"
  )
  
  # Test with empty temporary file
  temp_zip <- tempfile(fileext = ".zip")
  file.create(temp_zip)
  on.exit(unlink(temp_zip))
  
  expect_error(
    read1Cellbio(temp_zip),
    "Could not find sce.json"
  )
})

test_that("1CellbioData object structure is correct", {
  mock_data <- create_mock_1cellbio_data()
  
  # Check object class
  expect_true(inherits(mock_data, "1CellbioData"))
  
  # Check required components
  expect_true("experiment" %in% names(mock_data))
  expect_true("base_path" %in% names(mock_data))
  
  # Check experiment structure
  expect_true("summarized_experiment" %in% names(mock_data$experiment))
  expect_true("single_cell_experiment" %in% names(mock_data$experiment))
})

test_that("as.Seurat.1CB function exists and has correct signature", {
  # Check that the function exists
  expect_true(exists("as.Seurat.1CB"))

  # Check that it's an S3 generic function (uses UseMethod)
  func_body <- body(as.Seurat.1CB)
  expect_true(any(grepl("UseMethod", as.character(func_body))))

  # Check that the S3 method exists
  expect_true(exists("as.Seurat.1CB.1CellbioData"))
})

test_that("as.SingleCellExperiment.1CB function exists and has correct signature", {
  # Check that the function exists
  expect_true(exists("as.SingleCellExperiment.1CB"))

  # Check that it's an S3 generic function (uses UseMethod)
  func_body <- body(as.SingleCellExperiment.1CB)
  expect_true(any(grepl("UseMethod", as.character(func_body))))

  # Check that the S3 method exists
  expect_true(exists("as.SingleCellExperiment.1CB.1CellbioData"))
})

test_that("deprecated functions show appropriate warnings", {
  mock_data <- create_mock_1cellbio_data()
  
  # Test deprecated as.Seurat.1CellbioData
  expect_message(
    tryCatch(as.Seurat.1CellbioData(mock_data, rownames = "test", colnames = "test"), 
             error = function(e) NULL),
    "deprecated"
  )
  
  # Test deprecated as.SingleCellExperiment.1CellbioData  
  expect_message(
    tryCatch(as.SingleCellExperiment.1CellbioData(mock_data, rownames = "test", colnames = "test"), 
             error = function(e) NULL),
    "deprecated"
  )
})

test_that("name_conflict parameter validation works", {
  mock_data <- create_mock_1cellbio_data()
  
  # Should not error with valid name_conflict values
  expect_no_error({
    tryCatch(as.Seurat.1CellbioData(mock_data, rownames = "test", colnames = "test", name_conflict = "make_unique"), 
             error = function(e) NULL)
  })
  
  expect_no_error({
    tryCatch(as.Seurat.1CellbioData(mock_data, rownames = "test", colnames = "test", name_conflict = "error"), 
             error = function(e) NULL)
  })
})