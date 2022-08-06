#' GCIMSSpectrum class
#'
#' @description
#' GCIMSSpectrum is an S4 class to store a GCIMS Spectrum. It can be a single spectrum
#' or the aggregation of several spectra.
#'
#' @slot drift_time A numeric vector with drift times
#' @slot intensity A numeric vector witht the corresponding intensities
#' @slot retention_time_idx The index or indices used to get the intensity
#' @slot retention_time_s The retention times corresponding to the retention time indices.
#' @slot description A string with a description (used as plot title, useful e.g. to know the sample it came from)
#'
#' @export
methods::setClass(
  Class = "GCIMSSpectrum",
  slots = c(
    drift_time = "numeric",
    intensity = "numeric",
    retention_time_idx = "numeric",
    retention_time_s = "numeric",
    description = "character"
  )
)

methods::setMethod(
  "initialize", "GCIMSSpectrum",
  function(.Object, drift_time, intensity, retention_time_idx, retention_time_s, description) {
    stopifnot(length(drift_time) == length(intensity))
    .Object@drift_time <- drift_time
    .Object@retention_time_idx <- retention_time_idx
    .Object@retention_time_s <- retention_time_s
    .Object@intensity <- intensity
    .Object@description <- description
    .Object
  })


#' @describeIn GCIMSSpectrum class
#'
#' @param ... See the slots section in this page
#' @return A GCIMSSpectrum object
GCIMSSpectrum <- function(...) {
  methods::new("GCIMSSpectrum", ...)
}

methods::setMethod(
  "smooth", "GCIMSSpectrum",
  function(x, method = "savgol", dt_length_ms, dt_order = 2){
    dt <- dtime(x)
    dt_length_pts <- units_to_points(dt_length_ms, dt[2] - dt[1], must_odd = TRUE)
    x@intensity <- signal::sgolayfilt(x@intensity, n = dt_length_pts, p = dt_order)
    x
  })

#' @export
as.data.frame.GCIMSSpectrum <- function(x) {
  data.frame(
    drift_time_ms = dtime(x),
    intensity = intensity(x)
  )
}

#' @export
plot.GCIMSSpectrum <- function(x) {
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
      title = x@description,
      subtitle = subtitle,
      x = "Drift time (ms)",
      y = "Intensity (a.u.)"
    )
}
