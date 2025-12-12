library(testthat)
library(WEIN)
library(shiny)
library(shinytest2)

context("Shiny Modules")

# Test Shiny modules independently
test_that("ui_setup_server module functions correctly", {
  # Create a test server function that mimics the module
  test_server <- function() {
    # Create mock values
    values <- reactiveValues(
      dds_obj = NULL,
      annotation_obj = NULL,
      res_obj = NULL
    )
    
    # Test that the function exists
    expect_true(exists("ui_setup_server"))
    
    # We can't fully test the module without running the app,
    # but we can verify the function exists and is callable
    expect_type(ui_setup_server, "closure")
  }
  
  test_server()
})

test_that("module UI functions exist", {
  # Test that UI functions exist
  expect_true(exists("dash_Header_ui"))
  expect_true(exists("sideBar_ui"))
  expect_true(exists("welcome_panel_ui"))
  expect_true(exists("data_setup_ui"))
  expect_true(exists("panel_count_overview_ui"))
  expect_true(exists("extract_results_ui"))
  expect_true(exists("summaryPlots_ui"))
  expect_true(exists("geneFinder_ui"))
  expect_true(exists("panalAnalysis_ui"))
  expect_true(exists("signature_ui"))
  expect_true(exists("report_editor_ui"))
  expect_true(exists("about_ui"))
})

test_that("WEIN_ui function exists and is callable", {
  # Test that the main UI function exists
  expect_true(exists("WEIN_ui"))
  expect_type(WEIN_ui, "closure")
})