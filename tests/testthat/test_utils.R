library("testthat")
library("DESeq2")

# Source the functions directly instead of loading the package
source("../../R/WEIN.R")
source("../../R/helpers.R")

context("Testing Utility Functions")

test_that("multiAxPCA function works", {
  # Create a small test dataset
  dds <- DESeq2::makeExampleDESeqDataSet(n=100, m=8)
  
  # Test the function exists
  expect_true(exists("multiAxPCA"))
  
  # Test with default parameters
  p <- multiAxPCA(dds)
  expect_s3_class(p, "ggplot")
  
  # Test with custom parameters
  p2 <- multiAxPCA(dds, pc1=2, pc2=3)
  expect_s3_class(p2, "ggplot")
  
  # Test returnData parameter
  pct_var <- multiAxPCA(dds, returnData=TRUE)
  expect_type(pct_var, "double")
})

test_that("sepguesser function works", {
  # Test with a simple CSV-like text connection
  test_data <- "a,b,c\n1,2,3\n4,5,6"
  tc <- textConnection(test_data)
  writeLines(test_data, "test_file.txt")
  close(tc)
  
  # Test the function exists
  expect_true(exists("sepguesser"))
  
  # Test that it can guess comma separator
  sep <- sepguesser("test_file.txt")
  expect_equal(sep, ",")
  
  # Clean up
  unlink("test_file.txt")
})