#' Visualize results of the single-cell analysis pipeline
#'
#' This function generates a set of diagnostic plots from a processed
#' `SingleCellExperiment` object. It is designed to be used after
#' `run_pipeline()`, which computes quality control metrics, normalization,
#' highly variable genes (HVGs), dimensionality reduction, UMAP, and
#' clustering using Bioconductor workflows for single-cell analysis
#' (Amezquita et al. 2019; Lun et al. 2016; McCarthy et al. 2017).
#' The function produces:
#' \itemize{
#'   \item A mean–variance plot for HVGs, overlaying the fitted trend used
#'         during feature selection with `scran::modelGeneVar()`
#'         (Lun et al. 2016; Amezquita et al. 2019).
#'   \item Low-dimensional embeddings (PCA and UMAP), optionally coloured by
#'         cluster labels or another column in `colData(sce)`, using
#'         visualisation tools from `scater` (McCarthy et al. 2017).
#'   \item Basic QC plots (library size vs detected genes, and distributions
#'         of library size and detected genes) when QC metrics are available,
#'         following standard scRNA-seq practice (McCarthy et al. 2017).
#'   \item Simple cluster-level summaries (cluster sizes and mean library size
#'         per cluster) when clustering has been performed.
#' }
#' All plots are returned as a named list of `ggplot` objects. Optionally, they
#' can be printed to the active device and combined into a single multi-panel
#' figure that is saved to disk, leveraging the `ggplot2` and `patchwork`
#' ecosystems (Wickham 2016; Wickham 2015).
#'
#' @param sce A processed `SingleCellExperiment` object that has been run
#'   through `run_pipeline()`
#' @param colour_by Optional character scalar naming a column in
#'   `colData(sce)` used to colour PCA and UMAP plots (e.g. `"cluster_id"` or
#'   a sample/covariate column). If `NULL` (default), the function will:
#'   \enumerate{
#'     \item Use `"cluster_id"` when present in `colData(sce)`, or
#'     \item Fall back to the default colouring used by
#'           `scater::plotReducedDim()` (McCarthy et al. 2017).
#'   }
#' @param print_plots Logical; if `TRUE` (default), each plot in the output
#'   list is printed to the active graphics device in sequence. Set to `FALSE`
#'   to suppress on-screen plotting (useful in automated scripts and tests).
#' @param save_path Character string giving the file path where a combined
#'   multi-panel figure should be saved. The combined plot is assembled with
#'   `patchwork::wrap_plots()` and saved via `ggplot2::ggsave()`. The default
#'   is `"scPipeVis_combined_plot.png"`. Set to a path inside `tempdir()` in
#'   testing contexts.
#'
#' @return
#' Invisibly returns a named `list` of `ggplot` objects. Depending on which
#' components are present in `sce`, this list may include some or all of:
#' \itemize{
#'   \item `mean_variance` – HVG mean–variance plot with fitted trend
#'         (Amezquita et al. 2019; Lun et al. 2016)
#'   \item `pca` – PCA embedding
#'   \item `umap` – UMAP embedding
#'   \item `qc_lib_vs_detected` – library size vs detected genes scatter plot
#'   \item `qc_lib_hist` – histogram of library sizes
#'   \item `qc_detected_hist` – histogram of detected genes
#'   \item `cluster_sizes` – bar plot of cluster sizes
#'   \item `cluster_mean_library` – mean library size per cluster
#' }
#' As a side-effect, the function can print these plots and saves a combined
#' multi-panel summary figure to `save_path` when at least one plot is
#' generated.
#'
#' @details
#' This function assumes that `sce` has been processed by `run_pipeline()`
#'
#' The visualisation logic is organised as follows:
#' \enumerate{
#'   \item Validate that `metadata(sce)$hvg_data` exists and contains a fitted
#'         variance trend from `scran::modelGeneVar()`; ensure PCA and UMAP
#'         exist in `reducedDims(sce)`.
#'   \item Build an HVG mean–variance plot using `ggplot2` (Wickham 2016),
#'         overlaying the fitted trend from the HVG model.
#'   \item Generate PCA and UMAP plots using `scater::plotReducedDim()`
#'         (McCarthy et al. 2017), with optional colouring by `colour_by`
#'         (or by `cluster_id` when present).
#'   \item If QC metrics (`sum`, `detected`) are available in `colData(sce)`,
#'         produce basic QC plots: library size vs detected genes and
#'         distributions of both.
#'   \item If a `cluster_id` column is present, create simple cluster-level
#'         summary plots: cluster sizes and mean library size per cluster.
#'   \item Optionally print all plots and assemble non-`NULL` plots into a
#'         single multi-panel figure with `patchwork::wrap_plots()`, saving it
#'         to `save_path` via `ggplot2::ggsave()`.
#' }
#'
#' Design-wise, `visualise_pipeline()` follows a “single responsibility”
#' principle: it does not modify `sce`, but focuses solely on providing a
#' systematic, reproducible way to visualise the main outputs of
#' `run_pipeline()`. This supports both interactive
#' exploration and scripted/automated workflows in downstream analysis.
#'
#' @note
#' This implementation is intended for preliminary exploratory visualisation.
#' If required components (e.g., HVG metadata, PCA, UMAP) are missing,
#' informative errors are raised to help users diagnose whether the pipeline
#' was run correctly. QC- and cluster-based plots are only generated if the
#' relevant columns are present in `colData(sce)`.
#'
#' @examples
#' \dontrun{
#' library(scRNAseq)
#' library(SingleCellExperiment)
#' library(scPipeVis)
#'
#' sce <- ZeiselBrainData()
#' sce_processed <- run_pipeline(sce)
#'
#' # Generate and print all diagnostic plots, and save a combined PNG
#' plots <- visualise_pipeline(sce_processed)
#'
#' # Access the PCA plot directly
#' plots$pca
#'
#' # Use a custom colouring variable (e.g., cluster_id or a sample covariate)
#' visualise_pipeline(sce_processed, colour_by = "cluster_id")
#'
#' # Suppress printing and save only the combined figure to a temporary file
#' tmp_file <- file.path(tempdir(), "scPipeVis_overview.png")
#' visualise_pipeline(
#'   sce_processed,
#'   print_plots = FALSE,
#'   save_path   = tmp_file
#' )
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
#' Wickham H. 2016. ggplot2: Elegant Graphics for Data Analysis. Springer.
#'
#' Zeisel A, Munoz-Manchado AB, Codeluppi S, Lonnerberg P, La Manno G, Jureus A,
#' Marques S, Munguba H, He L, Betsholtz C, et al. 2015. Cell types in the mouse
#' cortex and hippocampus revealed by single-cell RNA-seq. Science.
#'
#' @export
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment colData
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom scater plotReducedDim
#' @importFrom ggplot2 ggplot aes geom_point geom_histogram geom_col labs theme_minimal stat_function ggsave
#' @importFrom patchwork wrap_plots
visualise_pipeline <- function(
    sce,
    colour_by   = NULL,
    print_plots = TRUE,
    save_path   = "scPipeVis_combined_plot.png"
) {
  ## Validate HVG info
  hvg_data <- S4Vectors::metadata(sce)$hvg_data
  if (is.null(hvg_data)) {
    stop("No 'hvg_data' found in metadata(sce). Did you run run_pipeline()?", call. = FALSE)
  }

  fit <- S4Vectors::metadata(hvg_data)
  if (is.null(fit) || is.null(fit$trend)) {
    stop("No trend fit found in metadata(hvg_data). Check scran::modelGeneVar().", call. = FALSE)
  }

  ##Sanely pick variance column: var / total / bio
  hvg_cols <- colnames(hvg_data)
  var_candidate <- intersect(c("var", "total", "bio"), hvg_cols)
  if (length(var_candidate) == 0L) {
    stop(
      "Could not find a variance column in hvg_data. ",
      "Looked for one of: 'var', 'total', 'bio'. ",
      "Available columns: ", paste(hvg_cols, collapse = ", "),
      call. = FALSE
    )
  }
  var_col_name <- var_candidate[1]

  ## Check reduced dimansions
  rd_names <- names(SingleCellExperiment::reducedDims(sce))
  if (!"PCA" %in% rd_names) {
    stop("Reduced dimension 'PCA' not found in reducedDims(sce).", call. = FALSE)
  }
  if (!"UMAP" %in% rd_names) {
    stop("Reduced dimension 'UMAP' not found in reducedDims(sce).", call. = FALSE)
  }

  ## Decide colouring logic
  cd <- SummarizedExperiment::colData(sce)
  if (is.null(colour_by)) {
    if ("cluster_id" %in% colnames(cd)) {
      colour_by <- "cluster_id"
    }
  } else if (!colour_by %in% colnames(cd)) {
    stop(
      sprintf("Column '%s' not found in colData(sce). ", colour_by),
      "Available columns: ", paste(colnames(cd), collapse = ", "),
      call. = FALSE
    )
  }

  ##Plot container
  plt_list <- list()

  ## HVG mean–variance plot
  mv_df <- data.frame(
    mean = hvg_data$mean,
    var  = hvg_data[[var_col_name]]
  )

  p_mean_var <- ggplot2::ggplot(mv_df, ggplot2::aes(x = mean, y = var)) +
    ggplot2::geom_point(alpha = 0.4) +
    ggplot2::stat_function(fun = fit$trend, colour = "dodgerblue", linewidth = 1) +
    ggplot2::labs(
      x = "Mean log-expression",
      y = paste0("Variance (", var_col_name, ")"),
      title = "HVG mean–variance trend"
    ) +
    ggplot2::theme_minimal()

  plt_list$mean_variance <- p_mean_var

  ## PCA & UMAP
  if (is.null(colour_by)) {
    p_pca  <- scater::plotReducedDim(sce, dimred = "PCA")
    p_umap <- scater::plotReducedDim(sce, dimred = "UMAP")
  } else {
    p_pca  <- scater::plotReducedDim(sce, dimred = "PCA",  colour_by = colour_by)
    p_umap <- scater::plotReducedDim(sce, dimred = "UMAP", colour_by = colour_by)
  }
  plt_list$pca  <- p_pca
  plt_list$umap <- p_umap

  ## QC plots (if QC metrics present)
  if (all(c("sum", "detected") %in% colnames(cd))) {
    qc_df <- as.data.frame(cd)

    plt_list$qc_lib_vs_detected <- ggplot2::ggplot(
      qc_df,
      ggplot2::aes(x = sum, y = detected)
    ) +
      ggplot2::geom_point(alpha = 0.4) +
      ggplot2::labs(
        x = "Library size (sum of counts)",
        y = "Number of detected genes",
        title = "QC: Library size vs detected genes"
      ) +
      ggplot2::theme_minimal()

    plt_list$qc_lib_hist <- ggplot2::ggplot(qc_df, ggplot2::aes(x = sum)) +
      ggplot2::geom_histogram(bins = 50) +
      ggplot2::labs(
        x = "Library size (sum of counts)",
        y = "Number of cells",
        title = "QC: Distribution of library sizes"
      ) +
      ggplot2::theme_minimal()

    plt_list$qc_detected_hist <- ggplot2::ggplot(qc_df, ggplot2::aes(x = detected)) +
      ggplot2::geom_histogram(bins = 50) +
      ggplot2::labs(
        x = "Number of detected genes",
        y = "Number of cells",
        title = "QC: Distribution of detected genes"
      ) +
      ggplot2::theme_minimal()
  } else {
    message(
      "QC metrics 'sum' and/or 'detected' not found in colData(sce). ",
      "Run run_pipeline(..., run_qc = TRUE) to enable QC plots."
    )
  }

  ## Cluster-level summaries (fi clustering present)
  if ("cluster_id" %in% colnames(cd)) {
    cl <- cd$cluster_id
    tab <- table(cl)
    cl_df <- data.frame(
      cluster = names(tab),
      n_cells = as.numeric(tab)
    )

    plt_list$cluster_sizes <- ggplot2::ggplot(
      cl_df,
      ggplot2::aes(x = cluster, y = n_cells)
    ) +
      ggplot2::geom_col() +
      ggplot2::labs(
        x = "Cluster ID",
        y = "Number of cells",
        title = "Cluster sizes"
      ) +
      ggplot2::theme_minimal()

    if ("sum" %in% colnames(cd)) {
      mean_by_cluster <- tapply(cd$sum, cl, mean, na.rm = TRUE)
      mean_df <- data.frame(
        cluster  = names(mean_by_cluster),
        mean_sum = as.numeric(mean_by_cluster)
      )

      plt_list$cluster_mean_library <- ggplot2::ggplot(
        mean_df,
        ggplot2::aes(x = cluster, y = mean_sum)
      ) +
        ggplot2::geom_col() +
        ggplot2::labs(
          x = "Cluster ID",
          y = "Mean library size (sum of counts)",
          title = "Mean library size by cluster"
        ) +
        ggplot2::theme_minimal()
    }
  } else {
    message(
      "Clustering not detected (colData(sce)$cluster_id missing). ",
      "Skipping cluster-level plots."
    )
  }

  ## Optionally print all plots
  if (isTRUE(print_plots)) {
    for (nm in names(plt_list)) {
      if (!is.null(plt_list[[nm]])) print(plt_list[[nm]])
    }
  }

  ## Build and save combined plot
  non_null_plots <- plt_list[!vapply(plt_list, is.null, logical(1))]
  if (length(non_null_plots) > 0) {
    combined <- patchwork::wrap_plots(non_null_plots)
    ggplot2::ggsave(
      filename = save_path,
      plot     = combined,
      width    = 14,
      height   = 10,
      dpi      = 300
    )
    message("Combined pipeline visualisation saved to: ",
            normalizePath(save_path))
  } else {
    message("No plots generated; nothing to save.")
  }

  invisible(plt_list)
}
