#' @describeIn GCIMSChromatogram-class Coerce to data frame
#' @inheritParams base::as.data.frame
#' @family GCIMSChromatogram
#' @export
as.data.frame.GCIMSChromatogram <- function(x, row.names = NULL, optional = FALSE, ...) {
  # Why would anyone give row.names or optional is beyond my understanding...
  # but that's how the generic is defined
  out <- data.frame(
    retention_time_s = rtime(x),
    intensity = intensity(x)
  )
  if (!is.null(x@baseline)) {
    out$baseline <- baseline(x)
  }
  if (!is.null(row.names)) {
    rownames(out) <- row.names
  }
  out
}
