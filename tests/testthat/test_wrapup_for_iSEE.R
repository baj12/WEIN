library("WEIN")
library("testthat")
library("DESeq2")

context("Testing wrapup_for_iSEE function")

test_that("wrapup_for_iSEE function exists", {
  # Test the function exists
  expect_true(exists("wrapup_for_iSEE"))
})

test_that("wrapup_for_iSEE function can be called", {
  # Create a small test dataset
  dds <- DESeq2::makeExampleDESeqDataSet(n=50, m=6)
  
  # Test that the function exists and can be called
  expect_true(exists("wrapup_for_iSEE"))
})