library("testthat")
library("DESeq2")
library("ggplot2")
library("ggrepel")
library("WEIN")

context("Testing ggplotCounts function")

test_that("ggplotCounts function exists", {
  # Test the function exists
  expect_true(exists("ggplotCounts"))
})

test_that("ggplotCounts works with basic DESeqDataSet", {
  # Create a small test dataset inside the test
  dds <- DESeq2::makeExampleDESeqDataSet(n=50, m=6)
  
  # Get a gene name to test with
  gene_name <- rownames(dds)[1]
  
  # Test that the function runs without error
  p <- ggplotCounts(dds, gene = gene_name, intgroup = "condition")
  
  # Check that it returns a ggplot object
  expect_is(p, "ggplot")
})

test_that("ggplotCounts handles different parameters", {
  # Create a small test dataset inside the test
  dds <- DESeq2::makeExampleDESeqDataSet(n=50, m=6)
  
  # Get a gene name to test with
  gene_name <- rownames(dds)[1]
  
  # Test with transform = FALSE
  p1 <- ggplotCounts(dds, gene = gene_name, intgroup = "condition", transform = FALSE)
  expect_is(p1, "ggplot")
  
  # Test with labels_repel = FALSE
  p2 <- ggplotCounts(dds, gene = gene_name, intgroup = "condition", labels_repel = FALSE)
  expect_is(p2, "ggplot")
})

test_that("ggplotCounts works with annotation object", {
  # Create a small test dataset inside the test
  dds <- DESeq2::makeExampleDESeqDataSet(n=50, m=6)
  
  # Get a gene name to test with
  gene_name <- rownames(dds)[1]
  
  # Create a simple annotation object
  annotation_obj <- data.frame(
    gene_id = rownames(dds),
    gene_name = paste0("GENE_", 1:nrow(dds)),
    stringsAsFactors = FALSE
  )
  rownames(annotation_obj) <- annotation_obj$gene_id
  
  # Test that the function runs without error with annotation
  p <- ggplotCounts(dds, gene = gene_name, intgroup = "condition", annotation_obj = annotation_obj)
  
  # Check that it returns a ggplot object
  expect_is(p, "ggplot")
})