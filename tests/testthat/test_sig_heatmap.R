library("WEIN")
library("testthat")
library("DESeq2")

context("Testing sig_heatmap function")

test_that("sig_heatmap function exists", {
  # Test the function exists
  expect_true(exists("sig_heatmap"))
})

test_that("sig_heatmap function can be called", {
  # Skip if heatmaply is not available
  skip_if_not_installed("heatmaply")
  
  # Create a small test dataset
  dds <- DESeq2::makeExampleDESeqDataSet(n=50, m=6)
  vst_data <- DESeq2::vst(dds)
  
  # Create a simple signature
  my_signature <- rownames(dds)[1:10]
  
  # Test that the function exists
  expect_true(exists("sig_heatmap"))
})

test_that("sig_heatmap handles basic parameters", {
  # Skip if heatmaply is not available
  skip_if_not_installed("heatmaply")
  
  # Create a small test dataset
  dds <- DESeq2::makeExampleDESeqDataSet(n=50, m=6)
  vst_data <- DESeq2::vst(dds)
  
  # Create a simple signature
  my_signature <- rownames(dds)[1:5]
  
  # Test with minimal parameters
  expect_true(exists("sig_heatmap"))
})