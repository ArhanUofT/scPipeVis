#' Run the Single-Cell RNA-Seq Processing Pipeline
#'
#' This function executes a standardized single-cell RNA-seq processing pipeline
#' using a provided `SingleCellExperiment` object. The pipeline includes normalization,
#' feature selection, dimensionality reduction, and calculating UMAP embedding. Methods and default values have been chosen for each step based on literature review so user can simply run the function for very preliminary results without having to worry about method selection and hyperparameter tuning at each step.
#' Each processing step is applied sequentially and the resulting object is stored in `sce` variable to conserve memory.
#'
#' @param sce A `SingleCellExperiment` object containing pre-quality-controlled
#' single-cell RNA-seq data.
#'
#' @return A processed `SingleCellExperiment` object containing normalized counts,
#' selected features, reduced dimensional representations (PCA), and UMAP coordinates.
#'
#' @details
#' The pipeline performs the following steps in order:
#' \enumerate{
#'   \item Normalization using `scuttle::logNormCounts()`
#'   \item Feature selection of highly variable genes using `scran::modelGeneVar()`
#'   \item Dimensionality reduction with `scran::fixedPCA()`
#'   \item UMAP embedding with `scater::runUMAP()`
#' }
#' Each step uses the corresponding helper function (`perform_normalization`,
#' `perform_feature_selection`, etc.).
#'
#' @note
#' Unsupervised clustering and automated annotations are to be added in to the pipeline to make it extensive.
#'
#' @examples
#' \dontrun{
#' # Load example data from the scRNAseq package
#' library(scRNAseq)
#' library(SingleCellExperiment)
#'
#' # Example dataset: Zeisel et al. 2015 mouse brain data
#' sce <- ZeiselBrainData()
#'
#' # Run the pipeline
#' sce_processed <- run_pipeline(sce)
#' }
#'
#' @export
#' @importFrom scuttle logNormCounts
#' @importFrom scran modelGeneVar getTopHVGs fixedPCA
#' @importFrom scater runUMAP addPerCellQC
#' @importFrom SingleCellExperiment rowData rowSubset
#' @importFrom S4Vectors metadata
run_pipeline <- function(sce) {
  # Assumption: Data is quality controlled and in a SingleCellExperiment object
  # We can design the pipeline to start from very scratch

  # we will keep updating the same sce object to save space
  sce <- perform_normalization(sce)
  sce <- perform_feature_selection(sce)
  sce <- perform_dimensionality_reduction(sce)
  sce <- calculate_umap(sce)

  # Clustering need to be added
  return(sce)
}

# Work in progress for qc function
perform_qc <- function(sce) {
  # Based on the literature perform best qc steps
  # Example:  removing mitochondrial genes, removing doublets, removing low gene count cells,
  # removing very lowly expressed genes?

  # performs qc and stores the results in sce
  # we will use ggplot2 to make some visulisation functions
  # There are a lot of things going on here, focus on only a few
  # We perform very basic qc filtering but rather assume that user has qc'd that data to a great extent
  # Extension of the pipeline can include full qc-ing

  # This adds basic qc information
  unfiltered_sce <- scater::addPerCellQC(sce)
  # Decide if we want to pass in the subsets argument

  # Do qc


  return(sce)
}

#####
# Helper functions

perform_normalization <- function(sce) {
  # Best kind of normalisation based on benchmarking of eltze and Huber, nature methods, 2023 is log(c+1) transforaation
  # can make this function more complicated by adding more arguments, like size factores true or false etc etc
  sce <- scuttle::logNormCounts(sce)
  return(sce)
}

perform_feature_selection <- function(sce, number=2000) {
  # Highly variable genes for feature slection step
  # Selecting highly variable genes provides a good way to do feature selection while preserving most of the biological variation.
  hvg_data <- scran::modelGeneVar(sce)
  chosen <- scran::getTopHVGs(hvg_data, n=number)
  SingleCellExperiment::rowSubset(sce) <- chosen

  metadata(sce)$hvg_data <- hvg_data
  # Instead of adding values
  return(sce)
}


perform_dimensionality_reduction <- function(sce) {
  # maybe let people calculate using ICA, NMF???
  set.seed(1008796812)
  sce <- scran::fixedPCA(sce, subset.row = SingleCellExperiment::rowData(sce)$subset)
  return(sce)
}

calculate_umap <- function(sce) {
  set.seed(1008796812) # My student number
  sce <- scater::runUMAP(sce, dimred="PCA")
  return(sce)
}
