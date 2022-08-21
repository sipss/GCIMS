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


#' @export
as.data.frame.GCIMSSpectrum <- function(x, ...) {
  data.frame(
    drift_time_ms = dtime(x),
    intensity = intensity(x)
  )
}


