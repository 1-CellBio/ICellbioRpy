test_that("read_csv_sample_info works correctly", {
  # Create a temporary CSV file for testing
  temp_csv <- tempfile(fileext = ".csv")
  
  # Create temporary files for MTX data
  temp_mtx <- tempfile(fileext = ".mtx")
  temp_features <- tempfile(fileext = ".tsv")
  temp_barcodes <- tempfile(fileext = ".tsv")
  
  # Create minimal test files
  writeLines(c(
    "%%MatrixMarket matrix coordinate integer general",
    "3 2 4",
    "1 1 1",
    "2 1 2", 
    "1 2 3",
    "3 2 4"
  ), temp_mtx)
  
  writeLines(c("Gene1\tGene1", "Gene2\tGene2", "Gene3\tGene3"), temp_features)
  writeLines(c("AAAA", "BBBB"), temp_barcodes)
  
  # Create CSV file
  csv_data <- data.frame(
    Sample_id = "test_sample",
    mtx_fns = temp_mtx,
    features_fns = temp_features,
    barcodes_fns = temp_barcodes,
    stringsAsFactors = FALSE
  )
  write.csv(csv_data, temp_csv, row.names = FALSE)
  
  # Test the function
  expect_no_error({
    sample_info <- read_csv_sample_info(temp_csv)
  })
  
  # Test that it returns correct structure
  sample_info <- read_csv_sample_info(temp_csv)
  expect_equal(nrow(sample_info), 1)
  expect_equal(sample_info$Sample_id, "test_sample")
  
  # Clean up
  unlink(c(temp_csv, temp_mtx, temp_features, temp_barcodes))
})

test_that("read_csv_sample_info handles missing files", {
  temp_csv <- tempfile(fileext = ".csv")
  
  # Create CSV with non-existent files
  csv_data <- data.frame(
    Sample_id = "test_sample",
    mtx_fns = "/nonexistent/path.mtx",
    features_fns = "/nonexistent/features.tsv",
    barcodes_fns = "/nonexistent/barcodes.tsv",
    stringsAsFactors = FALSE
  )
  write.csv(csv_data, temp_csv, row.names = FALSE)
  
  # Should throw an error
  expect_error(read_csv_sample_info(temp_csv), "do not exist")
  
  unlink(temp_csv)
})

test_that("read_csv_sample_info handles missing columns", {
  temp_csv <- tempfile(fileext = ".csv")
  
  # Create CSV with missing columns
  csv_data <- data.frame(
    Sample_id = "test_sample",
    mtx_fns = "/some/path.mtx",
    stringsAsFactors = FALSE
  )
  write.csv(csv_data, temp_csv, row.names = FALSE)
  
  # Should throw an error
  expect_error(read_csv_sample_info(temp_csv), "Missing required columns")
  
  unlink(temp_csv)
})