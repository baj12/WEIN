library("testthat")
library("DESeq2")

# Source the functions directly instead of loading the package
source("../../R/WEIN.R")
source("../../R/ui_components.R")
source("../../R/server_components.R")

context("Check that shiny app is generated")

dds <- DESeq2::makeExampleDESeqDataSet(n=100, m=8)

test_that("Shiny app is generated", {
  expect_is(WEIN(), "shiny.appobj")
  expect_is(WEIN(dds_obj = dds), "shiny.appobj")
})
