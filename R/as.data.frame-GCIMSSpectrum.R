#' @inheritParams base::as.data.frame
#' @export
as.data.frame.GCIMSSpectrum <- function(x, row.names = NULL, optional = FALSE, ...) {
  # Why would anyone give row.names or optional is beyond my understanding...
  # but that's how the generic is defined
  out <- data.frame(
    drift_time_ms = dtime(x),
    intensity = intensity(x)
  )
  if (!is.null(x@baseline)) {
    basel <- baseline(x)
    out$baseline <- basel
  }
  if (!is.null(row.names)) {
    rownames(out) <- row.names
  }
  out
}
