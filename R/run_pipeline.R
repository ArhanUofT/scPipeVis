run_pipeline <- function(sce) {
  # Assumption: Data is quality controlled and in a SingleCellExperiment object
  # We can design the pipeline to start from very scratch

}


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

perform_normalization <- function(sce) {
  # Best kind of normalisation based on benchmarking of eltze and Huber, nature methods, 2023 is log(c+1) transforaation
  # can make this function more complicated by adding more arguments, like size factores true or false etc etc
  sce <- scuttle::logNormCounts(sce)
  return(sce)
}

perform_feature_selection <- function(sce, number=2000) {
  # Highly variable genes for feature slection step

}


perform_dimensionality_reduction <- function(sce) {
  # maybe let people calculate using ICA, NMF???
}

