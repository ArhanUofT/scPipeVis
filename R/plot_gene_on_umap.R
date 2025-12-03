#' Plot Gene Expression on a UMAP Embedding
#'
#' This function visualizes the expression of a single gene on a UMAP embedding
#' by colouring each cell according to its log-normalized expression value. Such gene-level
#' visualizations are widely used in exploratory single-cell RNA-seq analysis to
#' assess marker gene localization across clusters and confirm biological
#' structure in the data .
#'
#' This function is intended to be applied **after** `run_pipeline()`, which
#' computes the the UMAP embedding.
#'
#' @param sce A processed `SingleCellExperiment` object containing both a
#'   `"logcounts"` assay and a `"UMAP"` matrix inside `reducedDims(sce)`. These
#'   components are generated automatically when running `run_pipeline()`.
#'
#' @param gene A character string giving the name of the gene to plot. The gene
#'   must be present in `rownames(sce)`.
#'
#' @details
#' The function extracts:
#' \enumerate{
#'   \item the two-dimensional UMAP coordinates from `reducedDims(sce)$UMAP`
#'   \item the log-normalized expression values for the specified gene from
#'         `assays(sce)$logcounts`
#' }
#' and constructs a scatter plot where each point represents a cell. The
#' expression values are mapped to colour using a **viridis** scale, which is
#' colourblind-friendly, making it suitable for
#' high-density scatterplots in scRNA-seq visualizations.
#'
#' UMAP-based gene-level visualization is a standard diagnostic tool used to:
#' \itemize{
#'   \item confirm cluster identity via known marker gene expression
#'   \item evaluate whether clusters align with expected biological signals
#'   \item detect potential technical artefacts such as gradients or batch
#'         effects
#' }
#'
#' @return Invisibly returns a `ggplot2` object. The function is called primarily
#'   for its side-effect of producing a UMAP visualization.
#'
#' @examples
#' \dontrun{
#' library(scRNAseq)
#' sce <- ZeiselBrainData()
#' sce_processed <- run_pipeline(sce)
#' plot_gene_on_umap(sce_processed, gene = "Tubb3")
#' }
#'
#' @references
#' Amezquita RA, Lun ATL, Becht E, Carey VJ, Carpp LN, Geistlinger L, Marini F,
#' Rue-Albrecht K, Risso D, Soneson C, et al. 2019. Orchestrating single-cell
#' analysis with Bioconductor. *Nature Methods*.
#'
#' Lun ATL, McCarthy DJ, Marioni JC. 2016. A step-by-step workflow for low-level
#' analysis of single-cell RNA-seq data with Bioconductor. *F1000Research*.
#'
#' McCarthy DJ, Campbell KR, Lun ATL, Wills QF. 2017. Scater: pre-processing,
#' quality control, normalization and visualization of single-cell RNA-seq data
#' in R. *Bioinformatics*.
#'
#' @export
#' @importFrom SummarizedExperiment assays
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom ggplot2 ggplot aes geom_point labs theme_minimal
#' @importFrom viridis scale_colour_viridis
plot_gene_on_umap <- function(sce, gene) {
  # Check UMAP embedding availability
  rd_names <- names(SingleCellExperiment::reducedDims(sce))
  if (!"UMAP" %in% rd_names) {
    stop(
      "Reduced dimension 'UMAP' not found in `reducedDims(sce)`. ",
      "Did you run `run_pipeline()` (with UMAP calculation) first?",
      call. = FALSE
    )
  }

  # Check normalized expression availability
  assay_names <- names(SummarizedExperiment::assays(sce))
  if (!"logcounts" %in% assay_names) {
    stop(
      "Assay 'logcounts' not found in assays(sce). ",
      "Did you run `run_pipeline()` (with normalization) first?",
      call. = FALSE
    )
  }

  # Check gene exists
  if (!gene %in% rownames(sce)) {
    stop(
      sprintf("Gene '%s' not found in rownames(sce).", gene),
      call. = FALSE
    )
  }

  # Extract coordinates and expression
  umap_mat <- SingleCellExperiment::reducedDims(sce)[["UMAP"]]
  expr <- SummarizedExperiment::assays(sce)[["logcounts"]][gene, ]
  expr <- as.numeric(expr)

  plot_df <- data.frame(
    UMAP1 = umap_mat[, 1],
    UMAP2 = umap_mat[, 2],
    expression = expr
  )

  # Plot with ggplot and viridis
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = UMAP1, y = UMAP2, colour = expression)) +
    ggplot2::geom_point(size = 0.7, alpha = 0.85) +
    viridis::scale_colour_viridis(
      option = "viridis",
      direction = 1,
      name = paste0(gene, " expression")
    ) +
    ggplot2::labs(
      title = paste("UMAP coloured by", gene, "expression"),
      x = "UMAP1",
      y = "UMAP2"
    ) +
    ggplot2::theme_minimal()

  print(p)
  invisible(p)
}
