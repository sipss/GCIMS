#' GCIMSChromatogram class
#'
#' @description
#' GCIMSChromatogram is an S4 class to store a GCIMS Chromatogram It can be a single chromatogram
#' or the aggregation of several chromatograms.
#'
#' @slot retention_time A numeric vector with retention times
#' @slot intensity A numeric vector with the corresponding intensities
#' @slot drift_time_idx The index or indices used to get the intensity
#' @slot drift_time_s The drift times corresponding to the retention time indices.
#' @slot description A string with a description (used as plot title, useful e.g. to know the sample it came from)
#'
#' @export
methods::setClass(
  Class = "GCIMSChromatogram",
  slots = c(
    retention_time = "numeric",
    intensity = "numeric",
    drift_time_idx = "numeric",
    drift_time_ms = "numeric",
    description = "character"
  )
)

methods::setMethod(
  "initialize", "GCIMSChromatogram",
  function(.Object, retention_time, intensity, drift_time_idx, drift_time_ms, description) {
    stopifnot(length(retention_time) == length(intensity))
    .Object@retention_time <- retention_time
    .Object@drift_time_idx <- drift_time_idx
    .Object@drift_time_ms <- drift_time_ms
    .Object@intensity <- intensity
    .Object@description <- description
    .Object
  })


#' @describeIn GCIMSChromatogram class
#'
#' @param ... See the slots section in this page
#' @return A GCIMSChromatogram object
GCIMSChromatogram <- function(...) {
  methods::new("GCIMSChromatogram", ...)
}

methods::setMethod(
  "smooth", "GCIMSChromatogram",
  function(x, method = "savgol", rt_length_s = 3, rt_order = 2L) {
    rt <- rtime(x)
    rt_length_pts <- units_to_points(rt_length_s, rt[2] - rt[1], must_odd = TRUE)
    if (rt_length_pts >= 1L) {
      x@intensity <- signal::sgolayfilt(x@intensity, n = rt_length_pts, p = rt_order)
    }
    x
  }
)

#' @export
as.data.frame.GCIMSChromatogram <- function(x) {
  data.frame(
    retention_time_s = rtime(x),
    intensity = intensity(x)
  )
}


#' @export
plot.GCIMSChromatogram <- function(x) {
  dt_ms <- x@drift_time_ms
  if (length(dt_ms) == 1) {
    subtitle <- glue("Drift time {dt_ms} ms")
  } else if (length(dt_ms) >= 2) {
    subtitle <- glue("Drift time {min(dt_ms)} - {max(dt_ms)} ms")
  } else {
    subtitle <- NULL
  }
  ggplot2::qplot(x = rtime(x), y = intensity(x), geom = "line") +
    ggplot2::labs(
      title = x@description,
      subtitle = subtitle
    )
}
