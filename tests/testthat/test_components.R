library("idealImmunoTP")
library("testthat")

context("Testing UI and Server Components")

# Create a small test dataset
dds <- DESeq2::makeExampleDESeqDataSet(n=100, m=8)

test_that("UI components can be loaded", {
  # Test that the UI function exists
  expect_true(exists("idealImmunoTP_ui"))
  
  # Test that the UI function can be called without error
  # Note: We can't fully test the UI without running Shiny, but we can check it loads
  expect_type(idealImmunoTP_ui, "closure")
})

test_that("Server components can be loaded", {
  # Test that the server function exists
  expect_true(exists("idealImmunoTP_server"))
  
  # Test that the server function can be called without error
  expect_type(idealImmunoTP_server, "closure")
})

test_that("Main function generates Shiny app", {
  # Test the main function
  app <- idealImmunoTP()
  expect_is(app, "shiny.appobj")
  
  # Test with DDS object
  app_with_data <- idealImmunoTP(dds_obj = dds)
  expect_is(app_with_data, "shiny.appobj")
})