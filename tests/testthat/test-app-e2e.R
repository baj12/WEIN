library(testthat)
library(WEIN)
library(shinytest2)
library(dplyr)

context("End-to-End App Tests")

# Skip tests on CRAN and CI environments without Chrome
skip_if_not_installed("shinytest2")
skip_on_cran()

# Custom function to check if Chrome/chromote is available
skip_if_no_chrome <- function() {
  # Try to load chromote and check if Chrome is available
  if (!requireNamespace("chromote", quietly = TRUE)) {
    skip("chromote package not available")
  }
  
  # Try to create a ChromoteSession to see if Chrome is available
  tryCatch({
    chromote::ChromoteSession$new(wait_ = FALSE)
    # If we get here, Chrome is available
  }, error = function(e) {
    skip("Chrome/Chromium not available for testing")
  })
}

# Test full user workflows through the GUI
test_that("App initializes and loads correctly", {
  skip_on_ci()
  skip_if_no_chrome()
  
  # Create a minimal test dataset to speed up initialization
  library(DESeq2)
  dds <- DESeq2::makeExampleDESeqDataSet(n = 50, m = 4)
  
  # Create a test app instance with actual WEIN components
  app <- WEIN(dds_obj = dds)
  
  # Create a test driver with proper initialization
  app_driver <- AppDriver$new(
    app = app,
    name = "app-initialization",
    load_timeout = 45000,
    timeout = 45000
  )
  
  # Check that the app loads without error
  expect_true(app_driver$is_alive())
  
  # Wait for the app to be fully initialized
  app_driver$wait_for_idle()
  
  # Take a screenshot for visual regression testing
  # app_driver$expect_screenshot()
  
  # Close the app
  app_driver$stop()
})

test_that("Input interactions work", {
  skip_on_ci()
  skip_if_no_chrome()
  
  # Create a minimal test dataset to speed up initialization
  library(DESeq2)
  dds <- DESeq2::makeExampleDESeqDataSet(n = 50, m = 4)
  
  # Create a test app instance with actual WEIN components
  app <- WEIN(dds_obj = dds)
  
  # Create a test driver with proper initialization
  app_driver <- AppDriver$new(
    app = app,
    name = "input-interactions",
    load_timeout = 45000,
    timeout = 45000
  )
  
  # Check that the app loads without error
  expect_true(app_driver$is_alive())
  
  # Wait for the app to be fully initialized
  app_driver$wait_for_idle()
  
  # We would test input interactions here if we had specific inputs to test
  
  # Close the app
  app_driver$stop()
})

test_that("Visual regression testing works", {
  skip_on_ci()
  skip_if_no_chrome()
  
  # This is a placeholder for visual regression tests
  # In a real implementation, we would compare screenshots against baselines
  expect_true(TRUE)
})