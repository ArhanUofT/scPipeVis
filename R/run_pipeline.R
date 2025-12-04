#' Run the Single-Cell RNA-Seq Processing Pipeline
#'
#' This function executes a standardized single-cell RNA-seq processing pipeline
#' using a provided `SingleCellExperiment` object. The pipeline includes:
#' optional quality control filtering, normalization, feature selection,
#' dimensionality reduction, UMAP embedding, and unsupervised clustering.
#' Methods and default values have been chosen based on benchmarking and
#' workflow literature so users can run a preliminary analysis with a single
#' function call (Amezquita et al. 2019; Lun et al. 2016).
#'
#' @param sce A `SingleCellExperiment` object containing single-cell RNA-seq
#'   data with a `"counts"` assay (McCarthy et al. 2017).
#' @param run_qc Logical; if `TRUE` (default), basic per-cell QC filtering
#'   is performed using library size and number of detected genes.
#' @param min_genes Numeric; minimum number of detected genes required for a
#'   cell to pass QC. Default is `200`.
#' @param min_counts Numeric; minimum total counts required for a cell to pass
#'   QC. Default is `1000`.
#' @param n_hvg Numeric; number of highly variable genes to select for
#'   downstream analysis. Default is `2000`.
#' @param n_pcs Numeric; number of principal components to retain in PCA.
#'   Default is `20`.
#' @param umap_n_neighbors Numeric; number of neighbors used by UMAP.
#'   Default is `15`.
#' @param cluster_dimred Character; which reduced dimension to use for
#'   clustering (e.g., `"PCA"` or `"UMAP"`). Default is `"PCA"`.
#'
#' @return A processed `SingleCellExperiment` object containing:
#' \itemize{
#'   \item QC metrics and an optional `qc_pass` column in `colData(sce)`
#'   \item normalized counts (`logcounts` assay)
#'   \item highly variable gene information stored in `metadata(sce)$hvg_data`
#'   \item PCA coordinates (`reducedDims(sce)$PCA`)
#'   \item UMAP coordinates (`reducedDims(sce)$UMAP`)
#'   \item cluster labels in `colData(sce)$cluster_id`
#'   \item pipeline configuration stored in `metadata(sce)$pipeline`
#' }
#'
#' @details
#' The pipeline performs the following steps in order:
#' \enumerate{
#'   \item Input validation to ensure `sce` is a `SingleCellExperiment`
#'         with a valid `"counts"` assay.
#'   \item (Optional) Per-cell QC with `scater::addPerCellQC()`, followed by
#'         filtering based on `min_genes` and `min_counts` (McCarthy et al. 2017).
#'   \item Normalization of counts using `scuttle::logNormCounts()` (Amezquita et al. 2019; McCarthy et al. 2017).
#'   \item Feature selection: highly variable genes via `scran::modelGeneVar()`
#'         and `scran::getTopHVGs()` (Lun et al. 2016).
#'   \item Dimensionality reduction using `scran::fixedPCA()` (Lun et al. 2016).
#'   \item UMAP embedding using `scater::runUMAP()` on the PCA space (Lun et al. 2016).
#'   \item Unsupervised clustering using `scater::clusterCells()` on the
#'         reduced dimension specified by `cluster_dimred` (Lun et al. 2016).
#' }
#'
#' Each step is delegated to a dedicated helper function
#' (`validate_sce_input()`, `perform_qc()`, `perform_normalization()`,
#' `perform_feature_selection()`, `perform_dimensionality_reduction()`,
#' `calculate_umap()`, `perform_clustering()`), which improves modularity,
#' readability, and testability of the code. The overall configuration of the
#' pipeline run is stored in `metadata(sce)$pipeline` for transparency and
#' reproducibility.
#'
#' @note
#' This implementation is intended for preliminary exploratory analysis.
#' More advanced features (e.g., doublet detection, batch correction, automated
#' cell-type annotation) can be layered on top of the returned object in
#' downstream analyses.
#'
#' @examples
#' \dontrun{
#' library(scRNAseq)
#' library(SingleCellExperiment)
#'
#' sce <- ZeiselBrainData()
#' sce_processed <- run_pipeline(
#'   sce,
#'   run_qc = TRUE,
#'   n_hvg = 2000,
#'   n_pcs = 20,
#'   umap_n_neighbors = 15,
#'   cluster_dimred = "PCA"
#' )
#'
#' summarize_clusters(sce_processed)
#' }
#'
#' @references
#' Amezquita RA, Lun ATL, Becht E, Carey VJ, Carpp LN, Geistlinger L, Marini F,
#' Rue-Albrecht K, Risso D, Soneson C, et al. 2019. Orchestrating single-cell
#' analysis with Bioconductor. Nature Methods.
#'
#' Lun ATL, McCarthy DJ, Marioni JC. 2016. A step-by-step workflow for low-level
#' analysis of single-cell RNA-seq data with Bioconductor. F1000Research. 5:2122.
#'
#' McCarthy DJ, Campbell KR, Lun ATL, Wills QF. 2017. Scater: pre-processing,
#' quality control, normalization and visualization of single-cell RNA-seq data
#' in R. Bioinformatics.
#'
#' Wickham H. 2015. R Packages. “O’Reilly Media, Inc.”
#'
#' Zeisel A, Munoz-Manchado AB, Codeluppi S, Lonnerberg P, La Manno G, Jureus A,
#' Marques S, Munguba H, He L, Betsholtz C, et al. 2015. Cell types in the mouse
#' cortex and hippocampus revealed by single-cell RNA-seq. Science.
#'
#' @export
#' @importFrom scuttle logNormCounts
#' @importFrom scran clusterCells modelGeneVar getTopHVGs fixedPCA
#' @importFrom scater runUMAP addPerCellQC
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom SummarizedExperiment assays colData rowData
#' @importFrom S4Vectors metadata
run_pipeline <- function(
    sce,
    run_qc = TRUE,
    min_genes = 200,
    min_counts = 1000,
    n_hvg = 2000,
    n_pcs = 20,
    umap_n_neighbors = 15,
    cluster_dimred = "PCA"
) {
  # Design principle: validate inputs early ("fail fast") with informative
  # error messages.
  sce <- validate_sce_input(sce)

  # Optional QC step: filter out low-quality cells before downstream analyses.
  if (run_qc) {
    sce <- perform_qc(
      sce,
      min_genes = min_genes,
      min_counts = min_counts
    )
  }

  # Design principle: run_pipeline() orchestrates the workflow; each step is
  # handled by a dedicated helper function with a single responsibility.
  sce <- perform_normalization(sce)
  sce <- perform_feature_selection(sce, number = n_hvg)
  sce <- perform_dimensionality_reduction(sce, n_pcs = n_pcs)
  sce <- calculate_umap(sce, n_neighbors = umap_n_neighbors)
  sce <- perform_clustering(sce, use_dimred = cluster_dimred)

  # Record pipeline configuration in metadata for transparency.
  metadata(sce)$pipeline <- list(
    qc = list(
      run_qc = run_qc,
      min_genes = min_genes,
      min_counts = min_counts
    ),
    normalization = list(
      method = "logNormCounts"
    ),
    feature_selection = list(
      method = "modelGeneVar/getTopHVGs",
      n_hvg = n_hvg
    ),
    dimred = list(
      method = "fixedPCA",
      n_pcs = n_pcs
    ),
    umap = list(
      method = "runUMAP",
      n_neighbors = umap_n_neighbors,
      dimred = "PCA"
    ),
    clustering = list(
      method = "clusterCells",
      dimred = cluster_dimred
    )
  )

  return(sce)
}

# Input validation

validate_sce_input <- function(sce) {
  # Design principle: separate validation logic for clarity and easy testing.
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("`sce` must be a SingleCellExperiment object.", call. = FALSE)
  }

  counts_names <- names(SummarizedExperiment::assays(sce))
  if (!"counts" %in% counts_names) {
    stop(
      "`sce` must contain a 'counts' assay (raw count matrix). ",
      "Use `assays(sce)$counts <- your_counts_matrix` before calling run_pipeline().",
      call. = FALSE
    )
  }

  counts <- SummarizedExperiment::assays(sce)[["counts"]]

  if (nrow(counts) == 0L || ncol(counts) == 0L) {
    stop(
      "`sce` must contain at least one gene (row) and one cell (column) ",
      "in the 'counts' assay.",
      call. = FALSE
    )
  }

  return(sce)
}

# Quality control helper

perform_qc <- function(sce, min_genes, min_counts) {
  # Design principle: QC is optional and configurable, but defaults are
  # chosen to be reasonable for typical scRNA-seq data.

  sce <- scater::addPerCellQC(sce)

  cd <- SummarizedExperiment::colData(sce)
  qc_pass <- cd$detected >= min_genes & cd$sum >= min_counts

  # Safely add qc_pass column and assign colData back.
  cd$qc_pass <- qc_pass
  SummarizedExperiment::colData(sce) <- cd

  sce_filtered <- sce[, qc_pass, drop = FALSE]
  metadata(sce_filtered)$qc_thresholds <- list(
    min_genes = min_genes,
    min_counts = min_counts
  )

  return(sce_filtered)
}

# Normalization

perform_normalization <- function(sce) {
  # Use log-normalization (logNormCounts) as a robust default.
  sce <- scuttle::logNormCounts(sce)
  return(sce)
}

# Feature selection

perform_feature_selection <- function(sce, number = 2000) {
  # Select highly variable genes (HVGs) for downstream analyses.
  hvg_data <- scran::modelGeneVar(sce)
  chosen <- scran::getTopHVGs(hvg_data, n = number)

  rd <- SummarizedExperiment::rowData(sce)
  # create a logical indicator for HVGs
  rd$is_hvg <- rownames(sce) %in% chosen
  SummarizedExperiment::rowData(sce) <- rd

  metadata(sce)$hvg_data <- hvg_data

  return(sce)
}

# Dimensionality reduction

perform_dimensionality_reduction <- function(sce, n_pcs = 20) {
  set.seed(1008796812)  # reproducibility

  rd <- SummarizedExperiment::rowData(sce)
  if (!"is_hvg" %in% colnames(rd)) {
    stop(
      "`rowData(sce)$is_hvg` not found. ",
      "Run perform_feature_selection() before perform_dimensionality_reduction().",
      call. = FALSE
    )
  }

  hvgs <- which(rd$is_hvg)
  sce <- scran::fixedPCA(
    sce,
    subset.row = hvgs,
    rank = n_pcs
  )
  return(sce)
}

# UMAP

calculate_umap <- function(sce, n_neighbors = 15) {
  set.seed(1008796812)  # reproducibility
  sce <- scater::runUMAP(
    sce,
    dimred = "PCA",
    n_neighbors = n_neighbors
  )
  return(sce)
}

# Clustering

perform_clustering <- function(sce, use_dimred = "PCA") {

  # Check that the chosen reduced dimension exists
  if (!use_dimred %in% names(SingleCellExperiment::reducedDims(sce))) {
    stop(
      sprintf("Reduced dimension '%s' not found in object.", use_dimred),
      call. = FALSE
    )
  }

  # Perform clustering
  clusters <- scran::clusterCells(sce, use.dimred = use_dimred)

  # If clusters returns a vector, use it directly
  if (is.atomic(clusters)) {
    cluster_labels <- as.factor(clusters)
  } else if ("cluster" %in% names(clusters)) {
    cluster_labels <- as.factor(clusters$cluster)
  } else {
    stop("Unexpected return type from scran::clusterCells().", call. = FALSE)
  }

  # Assign cluster labels to colData
  cd <- SummarizedExperiment::colData(sce)
  cd$cluster_id <- cluster_labels
  SummarizedExperiment::colData(sce) <- cd

  # Store full clustering result, regardless of type
  metadata(sce)$cluster_info <- clusters

  return(sce)
}

# Cluster summary helper


#' Summarize cluster sizes
#'
#' This helper function summarizes the number of cells in each cluster after
#' running \code{run_pipeline()}.
#'
#' @param sce A processed `SingleCellExperiment` object that contains
#'   `colData(sce)$cluster_id`.
#'
#' @return A `data.frame` with cluster IDs and corresponding cell counts.
#' @examples
#' \dontrun{
#' sce_processed <- run_pipeline(sce)
#' summarize_clusters(sce_processed)
#' }
#' @export
summarize_clusters <- function(sce) {
  cd <- SummarizedExperiment::colData(sce)

  if (!"cluster_id" %in% colnames(cd)) {
    stop(
      "`cluster_id` not found in `colData(sce)`. ",
      "Run `run_pipeline()` with clustering first.",
      call. = FALSE
    )
  }

  cl <- cd$cluster_id
  tab <- table(cl)
  df <- data.frame(
    cluster_id = names(tab),
    n_cells = as.integer(tab),
    row.names = NULL
  )
  return(df)
}
