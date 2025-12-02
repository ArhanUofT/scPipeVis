#' Visualize key steps of the single-cell analysis pipeline
#'
#' Produces diagnostic plots (HVG trend, PCA, UMAP, QC, clusters) and returns
#' them as a named list. Optionally prints them and also saves a *combined*
#' multi-panel figure as a PNG file using patchwork.
#'
#' @param sce A processed SingleCellExperiment object.
#' @param colour_by Optional column name in colData(sce) used to colour PCA/UMAP.
#' @param print_plots Logical; if TRUE prints all plots. Default TRUE.
#' @param save_path File path to save the combined plot. Default:
#'        "scPipeVis_combined_plot.png".
#'
#' @return A list of ggplot objects invisibly.
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
  ## --- Validate HVG info ---
  hvg_data <- S4Vectors::metadata(sce)$hvg_data
  if (is.null(hvg_data)) {
    stop("No 'hvg_data' found in metadata(sce). Did you run run_pipeline()?", call. = FALSE)
  }

  fit <- S4Vectors::metadata(hvg_data)
  if (is.null(fit) || is.null(fit$trend)) {
    stop("No trend fit found in metadata(hvg_data). Check scran::modelGeneVar().", call. = FALSE)
  }

  ## --- Sanely pick variance column: var / total / bio ---
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

  ## --- Check reduced dimensions ---
  rd_names <- names(SingleCellExperiment::reducedDims(sce))
  if (!"PCA" %in% rd_names) {
    stop("Reduced dimension 'PCA' not found in reducedDims(sce).", call. = FALSE)
  }
  if (!"UMAP" %in% rd_names) {
    stop("Reduced dimension 'UMAP' not found in reducedDims(sce).", call. = FALSE)
  }

  ## --- Decide colouring logic ---
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

  ## --- Plot container ---
  plt_list <- list()

  ## --- HVG mean–variance plot ---
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

  ## --- PCA & UMAP ---
  if (is.null(colour_by)) {
    p_pca  <- scater::plotReducedDim(sce, dimred = "PCA")
    p_umap <- scater::plotReducedDim(sce, dimred = "UMAP")
  } else {
    p_pca  <- scater::plotReducedDim(sce, dimred = "PCA",  colour_by = colour_by)
    p_umap <- scater::plotReducedDim(sce, dimred = "UMAP", colour_by = colour_by)
  }
  plt_list$pca  <- p_pca
  plt_list$umap <- p_umap

  ## --- QC plots (if QC metrics present) ---
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

  ## --- Cluster-level summaries (if clustering present) ---
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

  ## --- Optionally print all plots ---
  if (isTRUE(print_plots)) {
    for (nm in names(plt_list)) {
      if (!is.null(plt_list[[nm]])) print(plt_list[[nm]])
    }
  }

  ## --- Build and save combined plot ---
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
