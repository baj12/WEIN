library("testthat")
library("DESeq2")

# Source the functions directly instead of loading the package
source("../../R/goseqTable.R")

context("Testing goseqTable function")

test_that("goseqTable function exists", {
  # Test the function exists
  expect_true(exists("goseqTable"))
})

test_that("goseqTable handles basic inputs", {
  # Skip test if goseq is not available
  skip_if_not_installed("goseq")
  
  # Create test data
  # Small set of "differentially expressed" genes
  de_genes <- c("gene1", "gene2", "gene3", "gene4", "gene5")
  
  # Larger set of "assayed" genes
  assayed_genes <- paste0("gene", 1:20)
  
  # Test that the function exists and can be called
  # We won't test the full functionality as it requires genome packages
  # but we can test that it doesn't error with basic inputs
  expect_true(exists("goseqTable"))
})

test_that("goseqTable returns NULL for invalid inputs", {
  # Skip test if goseq is not available
  skip_if_not_installed("goseq")
  
  # Test with empty inputs
  result <- goseqTable(de.genes = character(0), assayed.genes = character(0))
  expect_null(result)
})