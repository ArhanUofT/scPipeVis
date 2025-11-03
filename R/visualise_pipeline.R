#' Visualize key steps of the single-cell analysis pipeline
#'
#' This function produces quick diagnostic plots from a processed
#' `SingleCellExperiment` object. It first visualizes the mean–variance
#' relationship estimated during feature selection (i.e. the output of
#' `scran::modelGeneVar()` that was stored in `metadata(sce)$hvg_data`), and then
#' plots low-dimensional embeddings (PCA and UMAP). If `colour_by` is supplied,
#' the reduced-dimension plots will be coloured by the specified column in
#' `colData(sce)`.
#'
#' @param sce A `SingleCellExperiment` object that has already been processed
#' by the pipeline function `run_pipeline()` and contains all the information required to create typical visualisations
#' that are part of analysing scRNA-seq data.
#' @param colour_by Optional. A character scalar naming a column in
#' `colData(sce)` to use for colouring cells in the PCA/UMAP plots. Defaults to
#' `NA`, in which case no colouring is applied.
#' Unsupervised clustering step is not yet implemented but extremely useful parameter after implementation of clustering
#'
#' @details
#' This function assumes that:
#' \itemize{
#'   \item `run_pipeline()` function was run on a `SingleCellExperiment` object and strictly that
#'    the output of that function is the input to this function.
#' }
#'
#' @return The function does not explicitely return anything but is called for its
#' side-effect of producing plots.
#'
#' @examples
#' \dontrun{
#' visualise_pipeline(sce)
#' visualise_pipeline(sce, colour_by = "cluster")
#' }
#'
#' @export
#' @importFrom S4Vectors metadata
#' @importFrom scater plotReducedDim
visualise_pipeline <- function(sce, colour_by = NA) {
  # Create few basic plots for now
  # Add more functionality after first submission
  # Right now the function is bare bones, but I will reimplement this using the tidyverse + ggplot2 framework

  # Visualising the mean variance relationship for the genes
  # maybe: visualise the elbow plot for PCA
  if (is.na(colour_by)) {

    dec <- metadata(sce)$hvg_data
    fit <- metadata(dec)
    plot(
      dec$mean,
      dec$var,
      xlab = "Mean of log-expression",
      ylab = "Variance of log-expression"
    )
    curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)
    # Add citation for code OCSA
    scater::plotReducedDim(sce, dimred="PCA")
    scater::plotReducedDim(sce, dimred="UMAP")
  } else {
    dec <- metadata(sce)$hvg_data
    fit <- metadata(dec)
    plot(
      dec$mean,
      dec$var,
      xlab = "Mean of log-expression",
      ylab = "Variance of log-expression"
    )
    curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)
    # Add citation for code OCSA

    # colour_by feature will be useful once unsupervised clustering is incorporated in the pipeline
    scater::plotReducedDim(sce, dimred="PCA", colour_by=colour_by)
    scater::plotReducedDim(sce, dimred="UMAP", colour_by=colour_by)
    }
}

# Break this functions into subfunctions
