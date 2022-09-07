#' @export
launchGCIMSgui <- function() {
  gui_file <- system.file("gui", "gui.R", package = "GCIMS")
  source(gui_file, local = TRUE)
  shiny::shinyApp(ui, server)
}
