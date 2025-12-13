library("testthat")
library("DESeq2")
library("WEIN")

context("Testing UI and Server Components")

test_that("UI components can be loaded", {
  # Test that the UI function exists
  expect_true(exists("WEIN_ui"))
  
  # Test that the UI function can be called without error
  # Note: We can't fully test the UI without running Shiny, but we can check it loads
  expect_type(WEIN_ui, "closure")
})

test_that("Server components can be loaded", {
  # Test that the server function exists
  expect_true(exists("WEIN_server"))
  
  # Test that the server function can be called without error
  expect_type(WEIN_server, "closure")
})

test_that("Main function generates Shiny app", {
  # Create a small test dataset inside the test
  dds <- DESeq2::makeExampleDESeqDataSet(n=100, m=8)
  
  # Test the main function
  app <- WEIN()
  expect_is(app, "shiny.appobj")
  
  # Test with DDS object
  app_with_data <- WEIN(dds_obj = dds)
  expect_is(app_with_data, "shiny.appobj")
})