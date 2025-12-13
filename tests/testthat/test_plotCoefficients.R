library("testthat")
library("DESeq2")
library("RColorBrewer")
library("WEIN")

context("Testing plotCoefficients function")

test_that("plotCoefficients function exists", {
  # Test the function exists
  expect_true(exists("plotCoefficients"))
})

test_that("plotCoefficients works with basic DESeqDataSet", {
  # Create a small test dataset inside the test
  dds <- DESeq2::makeExampleDESeqDataSet(n=50, m=6)
  dds <- DESeq2::DESeq(dds)
  
  # Get a gene name to test with
  gene_name <- rownames(dds)[1]
  
  # Test that the function runs without error
  # We use expect_error with NA to check that no error occurs
  expect_error(plotCoefficients(dds, geneName = gene_name), NA)
})
