run_pipeline <- function(sce) {
  # Assumption: Data is quality controlled and in a SingleCellExperiment object
  # We can design the pipeline to start from very scratch

}


perform_qc <- function(sce) {

}

perform_normalization <- function(sce) {
  # Based on the literature perform best qc steps
  # Example:  removing mitochondrial genes, removing doublets, removing low gene count cells,
  # removing very lowly expressed genes?
}

perform_feature_selection <- function(sce, number=2000) {

}


perform_dimensionality_reduction <- function(sce) {

}

