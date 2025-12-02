test_that("visualise_pipeline runs without error and returns list of ggplots (no colour_by)", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("scuttle")
  skip_if_not_installed("scran")
  skip_if_not_installed("scater")
  skip_if_not_installed("scRNAseq")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  library(scRNAseq)
  library(SingleCellExperiment)

  sce <- ZeiselBrainData()
  # subset for speed
  sce <- sce[seq_len(500), seq_len(min(200, ncol(sce)))]

  sce_proc <- run_pipeline(sce)

  # hvg_data is used to create the plot, if it doesn't exist then should expect error
  expect_true("hvg_data" %in% names(S4Vectors::metadata(sce_proc)))

  # Run visualise_pipeline, capture returned object, don't print plots, save to temp file
  tmp_file <- file.path(tempdir(), "scPipeVis_combined_test.png")

  plots <- NULL
  expect_no_error({
    plots <- visualise_pipeline(
      sce_proc,
      print_plots = FALSE,
      save_path   = tmp_file
    )
  })

  # Check that we got a list back
  expect_true(is.list(plots))

  # Check that some expected plots are present
  expect_true(all(c("mean_variance", "pca", "umap") %in% names(plots)))

  # Check that these are ggplot objects
  expect_s3_class(plots$mean_variance, "ggplot")
  expect_s3_class(plots$pca,           "ggplot")
  expect_s3_class(plots$umap,          "ggplot")

  # Optionally, check that the combined plot file was created
  expect_true(file.exists(tmp_file))
})
