test_that("run_pipeline adds logcounts, PCA, and UMAP using real scRNAseq data", {
  # Found this very useful function in R packages chapter 15
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("scuttle")
  skip_if_not_installed("scran")
  skip_if_not_installed("scater")
  skip_if_not_installed("scRNAseq")

  library(scRNAseq)
  library(SingleCellExperiment)

  sce <- ZeiselBrainData()  # Load dataset
  sce <- sce[seq_len(500), seq_len(200)]

  sce_proc <- run_pipeline(sce)

  expect_true("logcounts" %in% SummarizedExperiment::assayNames(sce_proc))
  expect_equal(dim(SummarizedExperiment::assay(sce_proc, "logcounts")),
               dim(SummarizedExperiment::assay(sce_proc, "counts")))

  rd_names <- SingleCellExperiment::reducedDimNames(sce_proc)
  expect_true("PCA" %in% rd_names)
  pca_mat <- SingleCellExperiment::reducedDim(sce_proc, "PCA")
  expect_gt(ncol(pca_mat), 1)

  expect_true("UMAP" %in% rd_names)
  umap_mat <- SingleCellExperiment::reducedDim(sce_proc, "UMAP")
  expect_equal(ncol(umap_mat), 2)

  expect_s4_class(sce_proc, "SingleCellExperiment")
})
