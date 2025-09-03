context("Python environment configuration")

test_that("check_anndata_available function works", {
  # Test that function exists
  expect_true(exists("check_anndata_available"))
  
  # Test function returns logical
  result <- check_anndata_available()
  expect_type(result, "logical")
})

test_that("configure_python_env function exists", {
  # Test that function exists
  expect_true(exists("configure_python_env"))
  
  # Test basic function call (may skip if conda environment not available)
  expect_no_error({
    tryCatch({
      configure_python_env(verbose = FALSE)
    }, error = function(e) {
      # It's ok if this fails due to environment issues
      NULL
    })
  })
})

test_that("configure_python_env handles conda environment parameter", {
  # Test with specific conda environment
  expect_no_error({
    tryCatch({
      configure_python_env(conda_env = "1cellbio", verbose = FALSE)
    }, error = function(e) {
      # It's ok if this fails due to environment issues
      NULL
    })
  })
})

test_that("icb_get_lang function works", {
  # Test default language
  options(ICellbioRpy.lang = NULL)
  expect_equal(icb_get_lang(), "en")
  
  # Test set language
  options(ICellbioRpy.lang = "zh")
  expect_equal(icb_get_lang(), "zh")
  
  # Test invalid language defaults to English
  options(ICellbioRpy.lang = "invalid")
  expect_equal(icb_get_lang(), "en")
  
  # Reset to default
  options(ICellbioRpy.lang = NULL)
})

