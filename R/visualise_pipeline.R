#' Visualize key steps of the single-cell analysis pipeline
#'
#' This function produces quick diagnostic plots from a processed
#' `SingleCellExperiment` object. It first visualizes the mean–variance
#' relationship estimated during feature selection (i.e. the output of
#' `scran::modelGeneVar()` that was stored in `metadata(sce)$hvg_data`), and then
#' plots low-dimensional embeddings (PCA and UMAP). If `colour_by` is supplied,
#' the reduced-dimension plots will be coloured by the specified column in
#' `colData(sce)`. If `colour_by` is not supplied and a `cluster_id` column is
#' present in `colData(sce)`, cells will be coloured by `cluster_id` by default.
#'
#' @param sce A `SingleCellExperiment` object that has already been processed
#'   by the pipeline function `run_pipeline()`.
#' @param colour_by Optional character scalar naming a column in
#'   `colData(sce)` to use for colouring cells in PCA/UMAP plots. If `NULL`
#'   (default), the function will:
#'   \enumerate{
#'     \item Use `"cluster_id"` when present in `colData(sce)`, or
#'     \item Fall back to default colouring if `cluster_id` is not available.
#'   }
#'
#' @return Invisibly returns `NULL`. Called for its side-effect of producing plots.
#'
#' @export
#' @importFrom S4Vectors metadata
#' @importFrom scater plotReducedDim
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom SummarizedExperiment colData
visualise_pipeline <- function(sce, colour_by = NULL) {
  # 1) Check HVG statistics are available
  hvg_data <- S4Vectors::metadata(sce)$hvg_data
  if (is.null(hvg_data)) {
    stop(
      "No 'hvg_data' found in metadata(sce). ",
      "Did you run `run_pipeline()` (with feature selection) first?",
      call. = FALSE
    )
  }

  fit <- S4Vectors::metadata(hvg_data)
  if (is.null(fit) || is.null(fit$trend)) {
    stop(
      "No trend fit found in metadata(hvg_data). ",
      "Check that feature selection via scran::modelGeneVar() completed correctly.",
      call. = FALSE
    )
  }

  # 2) Check PCA and UMAP are present
  rd_names <- names(SingleCellExperiment::reducedDims(sce))
  if (!"PCA" %in% rd_names) {
    stop(
      "Reduced dimension 'PCA' not found in `reducedDims(sce)`. ",
      "Did you run `run_pipeline()` (with dimensionality reduction) first?",
      call. = FALSE
    )
  }
  if (!"UMAP" %in% rd_names) {
    stop(
      "Reduced dimension 'UMAP' not found in `reducedDims(sce)`. ",
      "Did you run `run_pipeline()` (with UMAP calculation) first?",
      call. = FALSE
    )
  }

  # 3) Decide on colouring: default to cluster_id if present
  cd <- SummarizedExperiment::colData(sce)
  if (is.null(colour_by)) {
    if ("cluster_id" %in% colnames(cd)) {
      colour_by <- "cluster_id"
    }
  } else {
    if (!colour_by %in% colnames(cd)) {
      stop(
        sprintf("Column '%s' not found in colData(sce). ", colour_by),
        "Available columns include: ",
        paste(colnames(cd), collapse = ", "),
        call. = FALSE
      )
    }
  }

  # 4) Mean–variance plot for HVGs (base R plot)
  plot(
    hvg_data$mean,
    hvg_data$var,
    xlab = "Mean of log-expression",
    ylab = "Variance of log-expression"
  )
  curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

  # 5) PCA and UMAP plots (ggplot objects -> must be printed)
  if (is.null(colour_by)) {
    p_pca  <- scater::plotReducedDim(sce, dimred = "PCA")
    p_umap <- scater::plotReducedDim(sce, dimred = "UMAP")
  } else {
    p_pca  <- scater::plotReducedDim(sce, dimred = "PCA",  colour_by = colour_by)
    p_umap <- scater::plotReducedDim(sce, dimred = "UMAP", colour_by = colour_by)
  }

  print(p_pca)
  print(p_umap)

  invisible(NULL)
}
