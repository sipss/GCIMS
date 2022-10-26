#' @export
plot.GCIMSChromatogram <- function(x, ...) {
  dt_ms <- dtime(x)
  if (length(dt_ms) == 1) {
    subtitle <- glue("Drift time {dt_ms} ms")
  } else if (length(dt_ms) >= 2) {
    subtitle <- glue("Drift time {min(dt_ms)} - {max(dt_ms)} ms")
  } else {
    subtitle <- NULL
  }
  ggplot2::ggplot(
    data.frame(
      x = rtime(x),
      y = intensity(x)
    )
  ) + ggplot2::geom_line(ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::labs(
      x = "Retention time (s)",
      y = "Intensity (a.u.)",
      title = description(x),
      subtitle = subtitle
    )
}
