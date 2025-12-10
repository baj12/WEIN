library("idealImmunoTP")

context("Check that shiny app is generated")

dds <- DESeq2::makeExampleDESeqDataSet(n=100, m=8)

test_that("Shiny app is generated", {
  expect_is(idealImmunoTP(), "shiny.appobj")
  expect_is(idealImmunoTP(dds_obj = dds), "shiny.appobj")
})
