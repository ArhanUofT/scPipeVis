# inst/shiny-scripts/app.R

library(shiny)
library(scPipeVis)             # your package
library(scRNAseq)              # for ZeiselBrainData()
library(SingleCellExperiment)
library(SummarizedExperiment)
library(ggplot2)
library(viridis)

# Simple local validator for SCE objects used in the app
validate_sce_input <- function(sce) {
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("Input object must be a SingleCellExperiment.", call. = FALSE)
  }
  assay_names <- SummarizedExperiment::assayNames(sce)
  if (!"counts" %in% assay_names) {
    stop("SingleCellExperiment must contain a 'counts' assay.", call. = FALSE)
  }
  sce
}

ui <- fluidPage(
  titlePanel("scPipeVis: Gene expression on UMAP"),

  sidebarLayout(
    sidebarPanel(
      h4("Pipeline settings"),

      checkboxInput("run_qc", "Run QC filtering", value = TRUE),

      numericInput("n_hvg", "Number of HVGs", value = 2000, min = 100, step = 100),
      numericInput("n_pcs", "Number of PCs", value = 20, min = 5, step = 1),
      numericInput("umap_nn", "UMAP neighbours", value = 15, min = 5, step = 1),

      selectInput(
        "cluster_dimred",
        "Clustering on:",
        choices = c("PCA", "UMAP"),
        selected = "PCA"
      ),

      actionButton("run_pipeline", "Run pipeline"),

      hr(),
      h4("Visualisation options"),

      selectizeInput(
        "gene_name",
        "Gene to plot on UMAP:",
        choices = NULL,       # populated after pipeline runs
        multiple = FALSE,
        options = list(placeholder = "Select a gene…")
      )
    ),

    mainPanel(
      plotOutput("gene_umap_plot", height = "650px"),
      verbatimTextOutput("gene_info")
    )
  )
)

server <- function(input, output, session) {

  # ------------------------------------------------------------------
  # Base SCE: either user-supplied (user_sce) or demo Zeisel dataset
  # ------------------------------------------------------------------

  sce_raw <- reactiveVal(NULL)

  observeEvent(1, {
    if (exists("user_sce", envir = .GlobalEnv)) {
      sce <- get("user_sce", envir = .GlobalEnv)
      sce_raw(validate_sce_input(sce))
    } else {
      demo <- scRNAseq::ZeiselBrainData()
      sce_raw(validate_sce_input(demo))
    }
  }, once = TRUE)

  # ------------------------------------------------------------------
  # Processed SCE after running scPipeVis pipeline
  # ------------------------------------------------------------------

  sce_processed <- reactiveVal(NULL)

  observeEvent(input$run_pipeline, {
    base_sce <- sce_raw()
    req(base_sce)

    # Run your package pipeline
    sce_proc <- run_pipeline(
      base_sce,
      run_qc = input$run_qc,
      n_hvg = input$n_hvg,
      n_pcs = input$n_pcs,
      umap_n_neighbors = input$umap_nn,
      cluster_dimred = input$cluster_dimred
    )

    sce_processed(sce_proc)

    # Update gene dropdown with available genes
    gene_choices <- rownames(sce_proc)
    updateSelectInput(session, "gene_name", choices = gene_choices)
  })

  # ------------------------------------------------------------------
  # Gene expression on UMAP (only view in the app)
  # ------------------------------------------------------------------

  output$gene_umap_plot <- renderPlot({
    sce_proc <- sce_processed()
    req(sce_proc)

    gene <- input$gene_name
    req(!is.null(gene), nchar(gene) > 0)

    if (!gene %in% rownames(sce_proc)) {
      plot.new()
      title(main = paste("Gene", shQuote(gene), "not found in rownames(sce)."))
      return(invisible(NULL))
    }

    plot_gene_on_umap(sce_proc, gene = gene)
  })

  output$gene_info <- renderPrint({
    sce_proc <- sce_processed()
    req(sce_proc)

    gene <- input$gene_name
    req(!is.null(gene), nchar(gene) > 0)

    if (!gene %in% rownames(sce_proc)) {
      cat("Gene", shQuote(gene), "not found in rownames(sce).")
    } else {
      expr <- SummarizedExperiment::assays(sce_proc)[["logcounts"]][gene, ]
      cat("Showing UMAP expression for gene:", gene, "\n")
      cat("Mean log-expression:", mean(expr), "\n")
      cat("Min log-expression:", min(expr), "\n")
      cat("Max log-expression:", max(expr), "\n")
    }
  })
}

shinyApp(ui = ui, server = server)
