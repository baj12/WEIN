library("WEIN")
library("testthat")
library("DESeq2")

context("Testing New Utility Functions")

test_that("read_countmatrix function works", {
  # Create a simple test file
  test_data <- "gene_id\tsample1\tsample2\nGene1\t10\t20\nGene2\t30\t40"
  writeLines(test_data, "test_countmatrix.txt")
  
  # Test the function exists
  expect_true(exists("read_countmatrix"))
  
  # Test reading the count matrix
  cm <- read_countmatrix("test_countmatrix.txt")
  expect_is(cm, "data.frame")
  expect_equal(nrow(cm), 2)
  expect_equal(ncol(cm), 2)
  expect_equal(rownames(cm), c("Gene1", "Gene2"))
  
  # Clean up
  unlink("test_countmatrix.txt")
})

test_that("read_metadata function works", {
  # Create a simple test file
  test_data <- "sample_id\tcondition\nsample1\tcontrol\nsample2\ttreated"
  writeLines(test_data, "test_metadata.txt")
  
  # Test the function exists
  expect_true(exists("read_metadata"))
  
  # Test reading the metadata
  md <- read_metadata("test_metadata.txt")
  expect_is(md, "data.frame")
  expect_equal(nrow(md), 2)
  expect_equal(ncol(md), 2)
  
  # Clean up
  unlink("test_metadata.txt")
})

test_that("create_dds function exists", {
  # Test the function exists
  expect_true(exists("create_dds"))
})

test_that("generate_cooccurrence_plots function works", {
  # Create test data
  metadata <- data.frame(
    sample_id = c("sample1", "sample2", "sample3", "sample4"),
    condition = c("control", "control", "treated", "treated"),
    time = c("early", "late", "early", "late")
  )
  rownames(metadata) <- c("sample1", "sample2", "sample3", "sample4")
  
  # Test the function exists
  expect_true(exists("generate_cooccurrence_plots"))
  
  # Test generating plots
  plots <- generate_cooccurrence_plots(metadata, ~ condition + time)
  expect_true(is.null(plots) || inherits(plots, "ggplot"))
  
  # Test with null inputs
  plots_null <- generate_cooccurrence_plots(NULL, ~ condition)
  expect_null(plots_null)
})

test_that("generate_gene_lists function exists", {
  # Test the function exists
  expect_true(exists("generate_gene_lists"))
})