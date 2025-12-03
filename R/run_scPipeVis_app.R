#' Launch Shiny App for scPipeVis
#'
#' A function that launches the Shiny app for \code{scPipeVis}.
#' The purpose of this app is to demonstrate the single-cell RNA-seq
#' analysis pipeline implemented in the package. The app code has been
#' placed in \code{./inst/shiny-scripts} as required by the assignment guidelines.
#'
#' If a \code{SingleCellExperiment} object is supplied via the \code{sce}
#' argument, it will be made available to the app (as \code{user_sce})
#' and used as the input dataset. I am aware of the requirement of being able to upload or pass in users own data into the shiny app, however since the prerequisite of all the pipeline is to have the input in the form of an
#' \code{SingleCellExperiment} object and the user can't upload that to the shiny app like a csv or something so I decided that to pass in the data to the shiny app the user must provide the data in the form of a
#' \code{SingleCellExperiment} object directly as an argument to the function.
#'
#' Now as far as the demo dataset is concerened if \code{sce} argueent is \code{NULL}, the app
#' will instead load a built-in demo dataset automatically.
#'
#' The demoo/default dataset is the Zeisal et al 2015 mouse brain single cell dataset.
#'
#' @param sce Optional. A \code{SingleCellExperiment} object containing
#'   scRNA-seq data. If provided, this object is passed to the Shiny app
#'   via the global environment (as \code{user_sce}) and used instead of
#'   the built-in demo dataset.
#'
#' @return No meaningful return value; this function is called for its
#'   side-effect of launching a Shiny app in the browser.
#'
#' @examples
#' \dontrun{
#' # Launch app with demo dataset:
#' scPipeVis::run_scPipeVis_app()
#'
#' # Launch app with a user-supplied SingleCellExperiment:
#' # scPipeVis::run_scPipeVis_app(sce = my_sce)
#' }
#'
#' @references
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials.
#' \href{https://shiny.rstudio.com/tutorial/}{Link}
#'
#' Zeisel A, Munoz-Manchado AB, Codeluppi S, Lonnerberg P, La Manno G, Jureus A,
#' Marques S, Munguba H, He L, Betsholtz C, et al. 2015. Cell types in the mouse
#' cortex and hippocampus revealed by single-cell RNA-seq. Science.
#'
#' @export
#' @importFrom shiny runApp
run_scPipeVis_app <- function(sce = NULL) {
  # If the user supplies their own SingleCellExperiment, make it visible
  # to the Shiny app. The app code (in inst/shiny-scripts/app.R) can
  # check for `user_sce` and use it instead of the demo dataset.
  if (!is.null(sce)) {
    assign("user_sce", sce, envir = .GlobalEnv)
  }

  appDir <- system.file("shiny-scripts", package = "scPipeVis")

  if (appDir == "") {
    stop("Could not find Shiny app directory 'shiny-scripts' in the scPipeVis package.")
  }

  shiny::runApp(appDir, display.mode = "normal")
}
# [END]
