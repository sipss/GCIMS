#' @export
plot.GCIMSSpectrum <- function(x, ...) {
  rts <- rtime(x)
  if (length(rts) == 1) {
    subtitle <- glue("Retention time {rts} s")
  } else if (length(rts) >= 2) {
    subtitle <- glue("Retention time {min(rts)} - {max(rts)} s")
  } else {
    subtitle <- NULL
  }

  ggplot2::ggplot(
    data.frame(
      x = dtime(x),
      y = intensity(x)
    )
  ) + ggplot2::geom_line(ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::labs(
      title = description(x),
      subtitle = subtitle,
      x = "Drift time (ms)",
      y = "Intensity (a.u.)"
    )
}
