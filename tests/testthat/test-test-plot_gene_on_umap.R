test_that("plot_gene_on_umap runs without error and returns a ggplot", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("scuttle")
  skip_if_not_installed("scran")
  skip_if_not_installed("scater")
  skip_if_not_installed("scRNAseq")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("viridis")

  library(scRNAseq)
  library(SingleCellExperiment)

  # Use a small subset for speed
  sce <- ZeiselBrainData()
  sce <- sce[seq_len(500), seq_len(min(200, ncol(sce)))]

  # Run the pipeline to get logcounts + UMAP
  sce_proc <- run_pipeline(sce)

  # Pick a gene that exists in rownames
  gene_to_plot <- rownames(sce_proc)[1]

  p <- NULL
  expect_no_error({
    p <- plot_gene_on_umap(sce_proc, gene = gene_to_plot)
  })

  # Function should return a ggplot object (invisibly)
  expect_s3_class(p, "ggplot")
})

test_that("plot_gene_on_umap gives a clear error for missing gene", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("scuttle")
  skip_if_not_installed("scran")
  skip_if_not_installed("scater")
  skip_if_not_installed("scRNAseq")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("viridis")

  library(scRNAseq)
  library(SingleCellExperiment)

  sce <- ZeiselBrainData()
  sce <- sce[seq_len(500), seq_len(min(200, ncol(sce)))]

  sce_proc <- run_pipeline(sce)

  # Use a gene name that definitely doesn't exist
  expect_error(
    plot_gene_on_umap(sce_proc, gene = "THIS_GENE_DOES_NOT_EXIST"),
    "not found in rownames\\(sce\\)"
  )
})
