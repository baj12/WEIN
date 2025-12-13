library("testthat")
library("DESeq2")
library("ggplot2")
library("WEIN")

context("Integration Tests")

test_that("Complete workflow functions correctly", {
  # Create test data inside the test
  dds <- DESeq2::makeExampleDESeqDataSet(n=100, m=8)
  
  # Run DESeq analysis
  dds <- DESeq(dds)
  res <- results(dds)
  
  # Test that all major functions exist
  expect_true(exists("WEIN"))
  expect_true(exists("multiAxPCA"))
  expect_true(exists("plot_ma"))
  expect_true(exists("plot_volcano"))
  
  # Test that the main function creates an app
  app <- WEIN()
  expect_is(app, "shiny.appobj")
  
  # Test that the main function works with data
  app_with_data <- WEIN(dds_obj = dds, res_obj = res)
  expect_is(app_with_data, "shiny.appobj")
})

test_that("Plotting functions work with test data", {
  # Create test data inside the test
  dds <- DESeq2::makeExampleDESeqDataSet(n=100, m=8)
  dds <- DESeq(dds)
  res <- results(dds)
  
  # Test MA plot
  ma_plot <- plot_ma(res, FDR = 0.05)
  expect_is(ma_plot, "ggplot")
  
  # Test volcano plot
  volcano_plot <- plot_volcano(res)
  expect_is(volcano_plot, "ggplot")
})