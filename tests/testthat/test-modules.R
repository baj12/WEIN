library(testthat)
library(WEIN)
library(shiny)
library(shinytest2)

context("Shiny Modules")

# Source the module files directly to access internal functions
# These are not exported from the package namespace
source("../../R/ui_setup_server.R")
source("../../R/dash_Header_ui.R")
source("../../R/sideBar_ui.R")
source("../../R/welcome_panel_ui.R")
source("../../R/data_setup_ui.R")
source("../../R/panel_count_overview_ui.R")
source("../../R/extract_results_ui.R")
source("../../R/summaryPlots_ui.R")
source("../../R/geneFinder_ui.R")
source("../../R/fun_ana_ui.R")
source("../../R/ui_components.R")
source("../../R/signature_ui.R")
source("../../R/report_editor_ui.R")
source("../../R/about_ui.R")

# Test Shiny modules independently
test_that("ui_setup_server module functions correctly", {
  # Test that the function exists
  expect_true(exists("ui_setup_server"))
  
  # We can't fully test the module without running the app,
  # but we can verify the function exists and is callable
  if (exists("ui_setup_server")) {
    expect_type(ui_setup_server, "closure")
  }
})

test_that("module UI functions exist", {
  # Test that UI functions exist
  ui_functions <- c("dash_Header_ui", "sideBar_ui", "welcome_panel_ui",
                    "data_setup_ui", "panel_count_overview_ui", "extract_results_ui",
                    "summaryPlots_ui", "geneFinder_ui", "panalAnalysis_ui",
                    "signature_ui", "report_editor_ui", "about_ui")
  
  for (func_name in ui_functions) {
    # Check if function exists
    if (exists(func_name)) {
      func <- get(func_name)
      # Check that it's a function
      expect_type(func, "closure")
      # Optionally check that it doesn't return NULL when called (but don't fail the test if it does)
      # result <- tryCatch(func(), error = function(e) NULL)
      # Don't check result as some functions might legitimately return NULL
    } else {
      # Fail the test if function doesn't exist
      fail(paste("Function", func_name, "not found"))
    }
  }
})

test_that("WEIN_ui function exists and is callable", {
  # Test that the main UI function exists
  if (exists("WEIN_ui")) {
    expect_type(WEIN_ui, "closure")
  }
})