#' @export
plot.GCIMSSpectrum <- function(x, ...) {
  rts <- x@retention_time_s
  if (length(rts) == 1) {
    subtitle <- glue("Retention time {rts} s")
  } else if (length(rts) >= 2) {
    subtitle <- glue("Retention time {min(rts)} - {max(rts)} s")
  } else {
    subtitle <- NULL
  }
  ggplot2::qplot(x = dtime(x), y = intensity(x), geom = "line") +
    ggplot2::labs(
      title = description(x),
      subtitle = subtitle,
      x = "Drift time (ms)",
      y = "Intensity (a.u.)"
    )
}
