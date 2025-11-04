#' Run the Single-Cell RNA-Seq Processing Pipeline
#'
#' This function executes a standardized single-cell RNA-seq processing pipeline
#' using a provided `SingleCellExperiment` object. The pipeline includes normalization,
#' feature selection, dimensionality reduction, and calculating UMAP embedding. Methods and default values have been chosen for each step based on literature review so user can simply run the function for very preliminary results without having to worry about method selection and hyperparameter tuning at each step (Amezquita et al. 2019, Lun et al. 2016).
#' Each processing step is applied sequentially and the resulting object is stored in `sce` variable to conserve memory.
#'
#' @param sce A `SingleCellExperiment` object containing pre-quality-controlled
#' single-cell RNA-seq data (McCarthy et al. 2017).
#'
#' @return A processed `SingleCellExperiment` object containing normalized counts,
#' selected features, reduced dimensional representations (PCA), and UMAP coordinates.
#'
#' @details
#' The pipeline performs the following steps in order:
#' \enumerate{
#'   \item Normalization using `scuttle::logNormCounts()` (McCarthy et al. 2017)
#'   \item Feature selection of highly variable genes using `scran::modelGeneVar()` (McCarthy et al. 2017)
#'   \item Dimensionality reduction with `scran::fixedPCA()` (McCarthy et al. 2017)
#'   \item UMAP embedding with `scater::runUMAP()` (McCarthy et al. 2017)
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
#'@references
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

  # Reference for the code - in text citation as was told to be added in the class
  # (Amezquita et al. 2019)
  # Also added in the function documentation
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
