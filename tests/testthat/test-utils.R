context("Utility functions")

test_that("icb_i18n works correctly", {
  # Test default language (English)
  options(ICellbioRpy.lang = "en")
  expect_equal(icb_i18n("中文", "English"), "English")
  
  # Test Chinese language
  options(ICellbioRpy.lang = "zh")
  expect_equal(icb_i18n("中文", "English"), "中文")
  
  # Test invalid language defaults to English
  options(ICellbioRpy.lang = "invalid")
  expect_equal(icb_i18n("中文", "English"), "English")
  
  # Reset to default
  options(ICellbioRpy.lang = NULL)
})

test_that("icb_make_unique works correctly", {
  # Test make_unique strategy
  names1 <- c("A", "B", "A", "C")
  result1 <- icb_make_unique(names1, strategy = "make_unique")
  expect_equal(result1, c("A", "B", "A-1", "C"))
  
  # Test error strategy with duplicates
  expect_error(
    icb_make_unique(names1, strategy = "error"),
    "Duplicate names detected"
  )
  
  # Test no duplicates
  names2 <- c("A", "B", "C", "D")
  result2 <- icb_make_unique(names2, strategy = "error")
  expect_equal(result2, names2)
})

test_that("to_csparse works correctly", {
  skip_if_not_installed("Matrix")
  
  # Test dense matrix conversion
  dense_mat <- matrix(c(1, 0, 2, 0, 0, 3), nrow = 2)
  sparse_result <- to_csparse(dense_mat)
  expect_true(inherits(sparse_result, "CsparseMatrix"))
  
  # Test sparse matrix passthrough
  sparse_mat <- Matrix::Matrix(dense_mat, sparse = TRUE)
  sparse_result2 <- to_csparse(sparse_mat)
  expect_true(inherits(sparse_result2, "CsparseMatrix"))
})

test_that("construct_sparse_matrix works correctly", {
  skip_if_not_installed("Matrix")
  
  # Test basic sparse matrix construction
  data <- c(1, 2, 3)
  indices <- c(0, 1, 2)  # 0-based indexing
  indptr <- c(0, 1, 2, 3)
  shape <- c(3, 3)
  
  result <- construct_sparse_matrix(data, indices, indptr, shape)
  expect_true(inherits(result, "sparseMatrix"))
  expect_equal(dim(result), c(3, 3))
  expect_equal(Matrix::nnzero(result), 3)
})

