#' @export
launchGCIMSgui <- function() {
  gui_deps <- c("shiny", "plotly", "bslib", "knitr",
                "markdown", "BiocParallel", "ggplot2",
                "GCIMS", "shinyFiles", "shinycssloaders",
                "prompter", "yaml", "fs")
  require_pkgs(gui_deps)
  gui_file <- system.file("gui", "gui.R", package = "GCIMS")
  source(gui_file, local = TRUE)
  shiny::shinyApp(ui, server)
}

