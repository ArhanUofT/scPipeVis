#' Plot gene expression on UMAP using ggplot2 and viridis
#'
#' This function visualizes the expression of a single gene on the UMAP
#' embedding by colouring cells according to their (log-normalized) expression
#' values, using a viridis colour scale. It is intended to be used after
#' `run_pipeline()`, which computes log-normalized counts and UMAP coordinates.
#'
#' @param sce A processed `SingleCellExperiment` object containing
#'   a `"logcounts"` assay and a `"UMAP"` entry in `reducedDims(sce)`.
#' @param gene Character scalar giving the name of the gene to plot. Must exist
#'   in `rownames(sce)`.
#'
#' @return Invisibly returns the ggplot object. Called for its side-effect
#'   of producing a UMAP plot coloured by gene expression.
#'
#' @examples
#' \dontrun{
#' sce_processed <- run_pipeline(sce)
#' plot_gene_on_umap(sce_processed, gene = "Tubb3")
#' }
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
