library(shiny)
library(scPipeVis)
library(SingleCellExperiment)
library(scRNAseq)
library(SummarizedExperiment)
library(ggplot2)
library(viridis)

ui <- fluidPage(
  titlePanel("scPipeVis: scRNA-seq pipeline and visualisation"),

  sidebarLayout(
    sidebarPanel(
      h4("Pipeline settings"),

      checkboxInput("run_qc", "Run QC filtering", value = TRUE),

      numericInput("n_hvg", "Number of HVGs", value = 2000, min = 100),
      numericInput("n_pcs", "Number of PCs", value = 20, min = 5),
      numericInput("umap_nn", "UMAP neighbours", value = 15, min = 5),
      selectInput(
        "cluster_dimred",
        "Clustering on:",
        choices = c("PCA", "UMAP"),
        selected = "PCA"
      ),
      actionButton("run_pipeline", "Run pipeline"),

      hr(),
      h4("Visualisation options"),
      selectInput(
        "colour_by",
        "Colour PCA/UMAP by:",
        choices = c("Auto" = "auto"),
        selected = "auto"
      ),
      textInput("gene_name", "Gene to plot on UMAP:", value = "")
    ),

    mainPanel(
      tabsetPanel(
        id = "tabs",
        tabPanel("Pipeline overview", plotOutput("pipeline_plots", height = "900px")),
        tabPanel("Gene on UMAP",
                 plotOutput("gene_umap_plot", height = "650px"),
                 verbatimTextOutput("gene_info")
        )
      )
    )
  )
)

server <- function(input, output, session) {

  sce_raw <- reactiveVal(NULL)

  observeEvent(1, {
    if (exists("user_sce", envir = .GlobalEnv)) {
      sce_raw(validate_sce_input(get("user_sce", envir = .GlobalEnv)))
    } else {
      demo <- ZeiselBrainData()
      sce_raw(validate_sce_input(demo))
    }
  }, once = TRUE)

  sce_processed <- reactiveVal(NULL)

  observeEvent(input$run_pipeline, {
    base_sce <- sce_raw()
    req(base_sce)

    sce_proc <- run_pipeline(
      base_sce,
      run_qc = input$run_qc,
      n_hvg = input$n_hvg,
      n_pcs = input$n_pcs,
      umap_n_neighbors = input$umap_nn,
      cluster_dimred = input$cluster_dimred
    )

    sce_processed(sce_proc)

    cd <- colData(sce_proc)
    choices <- colnames(cd)
    choices <- c("Auto" = "auto", choices)
    updateSelectInput(session, "colour_by", choices = choices)
  })

  output$pipeline_plots <- renderPlot({
    sce_proc <- sce_processed()
    req(sce_proc)

    if (input$colour_by == "auto")
      visualise_pipeline(sce_proc, colour_by = NULL)
    else
      visualise_pipeline(sce_proc, colour_by = input$colour_by)
  })

  output$gene_umap_plot <- renderPlot({
    sce_proc <- sce_processed()
    req(sce_proc)

    gene <- input$gene_name
    req(nchar(gene) > 0)

    if (!gene %in% rownames(sce_proc)) {
      plot.new()
      title(main = paste("Gene", gene, "not found"))
      return(invisible(NULL))
    }
    plot_gene_on_umap(sce_proc, gene)
  })

  output$gene_info <- renderPrint({
    sce_proc <- sce_processed()
    req(sce_proc)

    gene <- input$gene_name
    req(nchar(gene) > 0)

    if (!gene %in% rownames(sce_proc)) {
      cat("Gene", gene, "not found in dataset.")
    } else {
      expr <- assays(sce_proc)[["logcounts"]][gene, ]
      cat("Summary for gene:", gene, "\n")
      print(summary(expr))
    }
  })
}

shinyApp(ui = ui, server = server)
