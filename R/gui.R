#' Launch the GCIMS graphical interface
#'
#' Launches the graphical interface. You may need to install additional packages.
#' @export
launchGCIMSgui <- function() {
  gui_deps <- c("shiny", "plotly", "bslib", "knitr",
                "markdown", "BiocParallel", "ggplot2",
                "GCIMS", "shinyFiles", "shinycssloaders",
                "prompter", "yaml", "fs")
  require_pkgs(gui_deps)
  app_dir <- system.file("gui", package = "GCIMS")
  if (app_dir == "") {
    stop("Could not find gui/ directory. Try re-installing `GCIMS`.", call. = FALSE)
  }
  shiny::runApp(app_dir, display.mode = "normal")
}

