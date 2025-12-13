library("testthat")
library("DESeq2")
library("ggplot2")

# Source the functions directly instead of loading the package
source("../../R/helpers.R")
source("../../R/utilities.R")

context("Testing New Utility Functions")

test_that("read_countmatrix function works", {
  # Create a simple test file
  test_data <- "gene_id\tsample1\tsample2\nGene1\t10\t20\nGene2\t30\t40"
  writeLines(test_data, "test_countmatrix.txt")
  
  # Test the function exists
  expect_true(exists("read_countmatrix"))
  
  # Test reading the count matrix
  cm <- read_countmatrix("test_countmatrix.txt")
  expect_s3_class(cm, "data.frame")
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
  expect_s3_class(md, "data.frame")
  expect_equal(nrow(md), 2)
  expect_equal(ncol(md), 1)
  expect_equal(rownames(md), c("sample1", "sample2"))
  
  # Clean up
  unlink("test_metadata.txt")
})

test_that("read_metadata handles typical metadata files correctly", {
  # Create a test file representing a typical metadata file
  # First column contains sample names, subsequent columns contain metadata
  test_data <- "sample\tcondition\ttreatment\nsample1\tcontrol\tuntreated\nsample2\ttreated\twith_drug"
  writeLines(test_data, "test_typical_metadata.txt")
  
  # Test reading the metadata
  md <- read_metadata("test_typical_metadata.txt")
  expect_s3_class(md, "data.frame")
  expect_equal(nrow(md), 2)
  # Should have 2 columns (condition, treatment) with row names from first column
  expect_equal(ncol(md), 2)
  expect_equal(rownames(md), c("sample1", "sample2"))
  expect_equal(colnames(md), c("condition", "treatment"))
  
  # Clean up
  unlink("test_typical_metadata.txt")
})

test_that("read_metadata handles header/data mismatch correctly", {
  # Create a test file where header has n-1 elements and data rows have n elements
  # This represents the case where first column contains row names but has no header
  test_data <- "condition\ttreatment\nsample1\tcontrol\tuntreated\nsample2\ttreated\twith_drug"
  writeLines(test_data, "test_header_mismatch.txt")
  
  # Test reading the metadata
  md <- read_metadata("test_header_mismatch.txt")
  expect_s3_class(md, "data.frame")
  expect_equal(nrow(md), 2)
  # Should have 2 columns (condition, treatment) with row names from first column
  expect_equal(ncol(md), 2)
  expect_equal(rownames(md), c("sample1", "sample2"))
  
  # Clean up
  unlink("test_header_mismatch.txt")
})

test_that("read_metadata handles files with n column headers and first column as sample names", {
  # Create a test file with n column headers where first column corresponds to sample names
  test_data <- "sample\tcondition\ttreatment\nsample1\tcontrol\tuntreated\nsample2\ttreated\twith_drug"
  writeLines(test_data, "test_metadata_n_columns.txt")
  
  # Test reading the metadata
  md <- read_metadata("test_metadata_n_columns.txt")
  expect_s3_class(md, "data.frame")
  expect_equal(nrow(md), 2)
  # Should have 2 columns (condition, treatment) with row names from first column
  expect_equal(ncol(md), 2)
  expect_equal(rownames(md), c("sample1", "sample2"))
  expect_equal(colnames(md), c("condition", "treatment"))
  
  # Clean up
  unlink("test_metadata_n_columns.txt")
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