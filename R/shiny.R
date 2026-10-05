#' Launch the interactive Monopoly Atlas Shiny app
#' @export
launch_shiny <- function() {
  if (!requireNamespace("shiny", quietly = TRUE))
    stop("shiny required: install.packages('shiny')")
  app_dir <- system.file("shiny", package = "monopolyAtlas")
  if (app_dir == "") stop("Shiny app not found in package installation")
  shiny::runApp(app_dir, display.mode = "normal")
}
