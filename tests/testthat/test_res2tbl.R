library("testthat")
library("DESeq2")

# Source the functions directly instead of loading the package
source("../../R/res2tbl.R")

context("Testing deseqresult2tbl and deseqresult2DEgenes functions")

test_that("deseqresult2tbl function exists", {
  # Test the function exists
  expect_true(exists("deseqresult2tbl"))
})

test_that("deseqresult2DEgenes function exists", {
  # Test the function exists
  expect_true(exists("deseqresult2DEgenes"))
})

test_that("deseqresult2tbl works with DESeqResults object", {
  # Create a small test dataset
  dds <- DESeq2::makeExampleDESeqDataSet(n=50, m=6)
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds)
  
  # Test the function
  result <- deseqresult2tbl(res)
  
  # Check that it returns a data frame
  expect_is(result, "data.frame")
  
  # Check that it has the expected columns
  expect_true("id" %in% colnames(result))
  
  # Check that it's arranged by padj
  if (nrow(result) > 1 && !all(is.na(result$padj))) {
    expect_true(all(diff(result$padj[!is.na(result$padj)]) >= 0))
  }
})

test_that("deseqresult2DEgenes works with DESeqResults object", {
  # Create a small test dataset
  dds <- DESeq2::makeExampleDESeqDataSet(n=50, m=6)
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds)
  
  # Test the function
  result <- deseqresult2DEgenes(res, FDR = 0.05)
  
  # Check that it returns a data frame
  expect_is(result, "data.frame")
  
  # Check that it has the expected columns
  expect_true("id" %in% colnames(result))
  
  # Check that all padj values are below the threshold (if any exist)
  if (nrow(result) > 0 && !all(is.na(result$padj))) {
    expect_true(all(result$padj[!is.na(result$padj)] <= 0.05))
  }
})

test_that("deseqresult2tbl handles invalid input", {
  # Test with invalid input
  expect_error(deseqresult2tbl("invalid_input"), "Not a DESeqResults object.")
})

test_that("deseqresult2DEgenes handles invalid input", {
  # Test with invalid input
  expect_error(deseqresult2DEgenes("invalid_input"), "Not a DESeqResults object.")
})