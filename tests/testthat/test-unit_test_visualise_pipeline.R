test_that("visualise_pipeline runs without error and returns a ggplot object (no colour_by)", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("scuttle")
  skip_if_not_installed("scran")
  skip_if_not_installed("scater")
  skip_if_not_installed("scRNAseq")
  skip_if_not_installed("ggplot2")

  library(scRNAseq)
  library(SingleCellExperiment)

  sce <- ZeiselBrainData()
  sce <- sce[seq_len(500), seq_len(min(200, ncol(sce)))]

  sce_proc <- run_pipeline(sce)

  # hvg_data is used to create the plot, if it doesn't exist then should expect error
  expect_true("hvg_data" %in% names(S4Vectors::metadata(sce_proc)))

  expect_no_error({
    p <- visualise_pipeline(sce_proc)
    expect_s3_class(p, "ggplot")
  })
})
