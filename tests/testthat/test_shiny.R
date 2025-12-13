library("testthat")
library("DESeq2")
library("WEIN")

context("Check that shiny app is generated")

test_that("Shiny app is generated", {
  # Create test data inside the test block
  dds <- DESeq2::makeExampleDESeqDataSet(n=100, m=8)
  expect_s3_class(WEIN(), "shiny.appobj")
  expect_s3_class(WEIN(dds_obj = dds), "shiny.appobj")
})
