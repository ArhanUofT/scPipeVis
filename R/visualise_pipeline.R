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
#' that are part of analysing scRNA-seq data (McCarthy et al. 2017).
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
#' @references
#' Amezquita RA, Lun ATL, Becht E, Carey VJ, Carpp LN, Geistlinger L, Marini F, Rue-Albrecht K, Risso D, Soneson C, et al. 2019 Dec 2. Orchestrating single-cell analysis with Bioconductor. Nature Methods. doi:https://doi.org/10.1038/s41592-019-0654-x.
#'
#' Lun ATL, McCarthy DJ, Marioni JC. 2016. A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor. F1000Research. 5:2122. doi:https://doi.org/10.12688/f1000research.9501.2.
#'
#' McCarthy DJ, Campbell KR, Lun ATL, Wills QF. 2017 Jan 14. Scater: pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R. Bioinformatics.:btw777. doi:https://doi.org/10.1093/bioinformatics/btw777.
#'
#' Wickham H. 2015. R Packages. “O’Reilly Media, Inc.”
#'
#' Zeisel A, Munoz-Manchado AB, Codeluppi S, Lonnerberg P, La Manno G, Jureus A, Marques S, Munguba H, He L, Betsholtz C, et al. 2015. Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. Science. 347(6226):1138–1142. doi:https://doi.org/10.1126/science.aaa1934.
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

    dec <- S4Vectors::metadata(sce)$hvg_data
    fit <- S4Vectors::metadata(dec)
    plot(
      dec$mean,
      dec$var,
      xlab = "Mean of log-expression",
      ylab = "Variance of log-expression"
    )
    curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

    scater::plotReducedDim(sce, dimred="PCA")
    scater::plotReducedDim(sce, dimred="UMAP")
  } else {
    dec <- S4Vectors::metadata(sce)$hvg_data
    fit <- S4Vectors::metadata(dec)
    plot(
      dec$mean,
      dec$var,
      xlab = "Mean of log-expression",
      ylab = "Variance of log-expression"
    )
    curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)


    # colour_by feature will be useful once unsupervised clustering is incorporated in the pipeline
    scater::plotReducedDim(sce, dimred="PCA", colour_by=colour_by)
    scater::plotReducedDim(sce, dimred="UMAP", colour_by=colour_by)

    # Reference for the code - in text citation as was told to be added in the class
    # (Amezquita et al. 2019)
    # Also added in the function documentation
    }
}

# Break this functions into subfunctions
