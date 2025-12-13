library("testthat")

# Source the functions directly instead of loading the package
source("../../R/utilities_internal.R")

context("Testing Internal Utility Functions")

test_that("mat_rowscale function works", {
  # Create a simple test matrix
  test_matrix <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
  
  # Test the function exists
  expect_true(exists("mat_rowscale"))
  
  # Test scaling the matrix
  scaled_matrix <- mat_rowscale(test_matrix)
  
  # Check that the result is a matrix
  expect_is(scaled_matrix, "matrix")
  
  # Check dimensions are preserved
  expect_equal(dim(scaled_matrix), dim(test_matrix))
  
  # Check that each row has mean approximately 0 and sd approximately 1
  row_means <- apply(scaled_matrix, 1, mean)
  row_sds <- apply(scaled_matrix, 1, sd)
  
  expect_equal(row_means, c(0, 0), tolerance = 1e-10)
  expect_equal(row_sds, c(1, 1), tolerance = 1e-10)
})

test_that("mat_rowscale handles NA values", {
  # Create a test matrix with NA values
  test_matrix <- matrix(c(1, 2, NA, 4, 5, 6), nrow = 2, ncol = 3)
  
  # Test scaling the matrix with NA values
  scaled_matrix <- mat_rowscale(test_matrix)
  
  # Check that the result is a matrix
  expect_is(scaled_matrix, "matrix")
  
  # Check dimensions are preserved
  expect_equal(dim(scaled_matrix), dim(test_matrix))
})

test_that("get_gene_list function works", {
  # Create mock values
  mock_values <- list(
    genelistUP = function() c("GENE1", "GENE2", "GENE3"),
    genelistDOWN = function() c("GENE4", "GENE5"),
    genelistUPDOWN = function() c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"),
    genelist1 = data.frame("Gene.Symbol" = c("GENE6", "GENE7")),
    genelist2 = data.frame("Gene.Symbol" = c("GENE8", "GENE9")),
    genelist3 = data.frame("Gene.Symbol" = c("GENE10")),
    genelist4 = data.frame("Gene.Symbol" = character(0))
  )
  
  # Test the function exists
  expect_true(exists("get_gene_list"))
  
  # Test retrieving UP gene list
  up_list <- get_gene_list("UP", mock_values)
  expect_equal(up_list, c("GENE1", "GENE2", "GENE3"))
  
  # Test retrieving DOWN gene list
  down_list <- get_gene_list("DOWN", mock_values)
  expect_equal(down_list, c("GENE4", "GENE5"))
  
  # Test retrieving UPDOWN gene list
  updown_list <- get_gene_list("UPDOWN", mock_values)
  expect_equal(updown_list, c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5"))
  
  # Test retrieving LIST1 gene list
  list1_result <- get_gene_list("LIST1", mock_values)
  expect_equal(list1_result, c("GENE6", "GENE7"))
  
  # Test retrieving LIST2 gene list
  list2_result <- get_gene_list("LIST2", mock_values)
  expect_equal(list2_result, c("GENE8", "GENE9"))
  
  # Test retrieving LIST3 gene list
  list3_result <- get_gene_list("LIST3", mock_values)
  expect_equal(list3_result, c("GENE10"))
  
  # Test retrieving LIST4 gene list (empty)
  list4_result <- get_gene_list("LIST4", mock_values)
  expect_equal(list4_result, character(0))
  
  # Test retrieving unknown gene list
  unknown <- get_gene_list("UNKNOWN", mock_values)
  expect_null(unknown)
})

test_that("create_goterm_heatmap function exists", {
  # Test the function exists
  expect_true(exists("create_goterm_heatmap"))
})

test_that("create_gene_list_handler function exists", {
  # Test the function exists
  expect_true(exists("create_gene_list_handler"))
})

test_that("create_ui_output function exists", {
  # Test the function exists
  expect_true(exists("create_ui_output"))
})

test_that("create_dt_output function exists", {
  # Test the function exists
  expect_true(exists("create_dt_output"))
})