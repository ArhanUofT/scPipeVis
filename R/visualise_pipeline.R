#' Visualize key steps of the single-cell analysis pipeline
#'
#' This function produces quick diagnostic plots from a processed
#' `SingleCellExperiment` object. It visualizes the mean–variance relationship
#' estimated during feature selection (output of `scran::modelGeneVar()`
#' stored in `metadata(sce)$hvg_data`), low-dimensional embeddings (PCA and
#' UMAP), and, when available, basic QC and cluster-level summaries.
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
  # --- Step 1: Validate HVG data ---
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

  # --- Step 2: Confirm PCA and UMAP exist ---
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

  # --- Step 3: Decide colouring logic ---
  cd <- SummarizedExperiment::colData(sce)
  if (is.null(colour_by)) {
    if ("cluster_id" %in% colnames(cd)) {
      colour_by <- "cluster_id"
    }
  } else if (!colour_by %in% colnames(cd)) {
    stop(
      sprintf("Column '%s' not found in colData(sce). ", colour_by),
      "Available columns include: ",
      paste(colnames(cd), collapse = ", "),
      call. = FALSE
    )
  }

  # --- Step 4: Mean–variance plot for HVGs ---
  plot(
    hvg_data$mean,
    hvg_data$var,
    xlab = "Mean of log-expression",
    ylab = "Variance of log-expression",
    main = "Mean–Variance Trend of Genes"
  )
  curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

  # --- Step 5: PCA & UMAP plots ---
  if (is.null(colour_by)) {
    p_pca  <- scater::plotReducedDim(sce, dimred = "PCA")
    p_umap <- scater::plotReducedDim(sce, dimred = "UMAP")
  } else {
    p_pca  <- scater::plotReducedDim(sce, dimred = "PCA",  colour_by = colour_by)
    p_umap <- scater::plotReducedDim(sce, dimred = "UMAP", colour_by = colour_by)
  }

  print(p_pca)
  print(p_umap)

  # --- Step 6: QC plots (if QC metrics exist) ---
  if (all(c("sum", "detected") %in% colnames(cd))) {
    # Library size vs detected genes
    plot(
      cd$sum,
      cd$detected,
      xlab = "Library size (sum of counts)",
      ylab = "Number of detected genes",
      main = "QC: Library size vs detected genes"
    )

    # Histogram of library size
    hist(
      cd$sum,
      breaks = 50,
      xlab = "Library size (sum of counts)",
      main = "QC: Distribution of library sizes",
      col = "grey80",
      border = "grey40"
    )

    # Histogram of detected genes
    hist(
      cd$detected,
      breaks = 50,
      xlab = "Number of detected genes",
      main = "QC: Distribution of detected genes",
      col = "grey80",
      border = "grey40"
    )
  } else {
    message(
      "QC metrics 'sum' and/or 'detected' not found in colData(sce). ",
      "Run `run_pipeline()` with run_qc = TRUE to enable QC plots."
    )
  }

  # --- Step 7: Cluster-level summary plots (if clustering ran) ---
  if ("cluster_id" %in% colnames(cd)) {
    cl <- cd$cluster_id
    tab <- table(cl)

    # a) Cluster size barplot
    barplot(
      tab,
      xlab = "Cluster ID",
      ylab = "Number of cells",
      main = "Cluster sizes",
      col = "grey80",
      border = "grey50"
    )

    # b) Mean QC metric ('sum') per cluster – only if 'sum' exists
    if ("sum" %in% colnames(cd)) {
      mean_by_cluster <- tapply(cd$sum, cl, mean, na.rm = TRUE)
      barplot(
        mean_by_cluster,
        xlab = "Cluster ID",
        ylab = "Mean library size (sum)",
        main = "Mean library size by cluster",
        col = "grey85",
        border = "grey40"
      )
    }
  } else {
    message(
      "Clustering not detected (colData(sce)$cluster_id missing). ",
      "Skipping cluster-level plots."
    )
  }

  invisible(NULL)
}
