test_that("visualise_pipeline fails if sce has not been processed by run_pipeline()", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("scRNAseq")
  skip_if_not_installed("testthat")


  library(scRNAseq)
  library(SingleCellExperiment)


  sce <- ZeiselBrainData()
  sce <- sce[seq_len(500), seq_len(100)]

  expect_false("logcounts" %in% SummarizedExperiment::assayNames(sce))
  expect_false("PCA" %in% SingleCellExperiment::reducedDimNames(sce))
  expect_false("UMAP" %in% SingleCellExperiment::reducedDimNames(sce))
  expect_null(S4Vectors::metadata(sce)$hvg_data)


  # no preocessing done before applying visualise_pipeline function so fail
  expect_error(
    visualise_pipeline(sce))
})
