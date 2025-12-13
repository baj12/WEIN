library("testthat")

# Source the functions directly instead of loading the package
source("../../R/genesignatures.R")

context("Testing read_gmt function")

test_that("read_gmt function exists", {
  # Test the function exists
  expect_true(exists("read_gmt"))
})

test_that("read_gmt works with basic GMT file", {
  # Create a temporary GMT file for testing
  tmp_gmt <- tempfile(fileext = ".gmt")
  
  # Write test data to the file
  cat("PATHWAY1\tsource1\tGENE1\tGENE2\tGENE3\n", file = tmp_gmt)
  cat("PATHWAY2\tsource2\tGENE4\tGENE5\n", file = tmp_gmt, append = TRUE)
  
  # Test the function
  result <- read_gmt(tmp_gmt)
  
  # Check that it returns a list
  expect_type(result, "list")
  
  # Check that it has the right number of elements
  expect_equal(length(result), 2)
  
  # Check the names
  expect_equal(names(result), c("PATHWAY1", "PATHWAY2"))
  
  # Check the contents
  expect_equal(result$PATHWAY1, c("GENE1", "GENE2", "GENE3"))
  expect_equal(result$PATHWAY2, c("GENE4", "GENE5"))
  
  # Clean up
  unlink(tmp_gmt)
})

test_that("read_gmt handles empty file", {
  # Create an empty temporary file
  tmp_gmt <- tempfile(fileext = ".gmt")
  file.create(tmp_gmt)
  
  # Test the function with empty file
  expect_warning(result <- read_gmt(tmp_gmt))
  expect_equal(length(result), 0)
  
  # Clean up
  unlink(tmp_gmt)
})

test_that("read_gmt handles non-existent file", {
  # Test with non-existent file
  expect_error(read_gmt("non_existent_file.gmt"), "GMT file does not exist")
})