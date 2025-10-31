run_pipeline <- function(sce) {
  # Assumption: Data is quality controlled and in a SingleCellExperiment object
  # We can design the pipeline to start from very scratch

}


perform_qc <- function(sce) {
  # performs qc and stores the results in sce
  # we will use ggplot2 to make some visulisation functions
  # There are a lot of things going on here, focus on only a few

  # This adds basic qc information
  sce <- scater::addPerCellQC(sce)
  # Decide if we want to pass in the subsets argument

  # Do qc


  return(sce)
}

perform_normalization <- function(sce) {
  # Based on the literature perform best qc steps
  # Example:  removing mitochondrial genes, removing doublets, removing low gene count cells,
  # removing very lowly expressed genes?
}

perform_feature_selection <- function(sce, number=2000) {

}


perform_dimensionality_reduction <- function(sce) {
  # maybe let people calculate using ICA, NMF???
}

