#' Launch Shiny App for scPipeVis
#'
#' A function that launches the Shiny app for \code{scPipeVis}.
#' The purpose of this app is to demonstrate the single-cell RNA-seq
#' analysis pipeline implemented in the package. The app code has been
#' placed in \code{./inst/shiny-scripts}.
#'
#' If a \code{SingleCellExperiment} object is supplied via the \code{sce}
#' argument, it will be made available to the app (as \code{user_sce})
#' and used as the input dataset. If \code{sce} is \code{NULL}, the app
#' will instead load a built-in demo dataset.
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
