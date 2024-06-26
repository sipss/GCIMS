#' GCIMSSpectrum class
#'
#' @description
#' GCIMSSpectrum is an S4 class to store a GCIMS Spectrum. It can be a single spectrum
#' or the aggregation of several spectra.
#'
#' @slot drift_time A numeric vector with drift times
#' @slot intensity A numeric vector with the corresponding intensities
#' @slot baseline A numeric vector of the same length as `intensity` with the corresponding baseline,
#' or `NULL` if not set. Use [estimateBaseline()] to estimate it, [baseline()] to directly access it.
#' @slot retention_time_idx The index or indices used to get the intensity
#' @slot retention_time_s The retention times corresponding to the retention time indices.
#' @slot description A string with a description (used as plot title, useful e.g. to know the sample it came from)
#' @slot peaks  A data frame with peaks, use [findPeaks()] to fill it, or [peaks()] to set/get it
#' @slot peaks_debug_info A list with arbitrary debug information from [findPeaks()]
#'
#' @export
methods::setClass(
  Class = "GCIMSSpectrum",
  slots = c(
    drift_time = "numeric",
    intensity = "numeric",
    baseline = "numericOrNULL",
    retention_time_idx = "numeric",
    retention_time_s = "numeric",
    description = "character",
    peaks = "DataFrameOrNULL",
    peaks_debug_info = "listOrNULL"
  )
)

methods::setMethod(
  "initialize", "GCIMSSpectrum",
  function(.Object, drift_time, intensity, retention_time_idx = NA_integer_, retention_time_s = NA_real_, description = "", ..., baseline = NULL) {
    dots <- list(...)
    if (length(dots) > 0) {
      if (!is.null(names(dots))) {
        invalid_names <- paste0(names(dots), collapse = ", ")
      }
      valid_names <- paste0(setdiff(names(formals()), c(".Object", "...")), collapse = ", ")
      cli_abort(
        message = c(
          "Invalid argument to GCIMSSpectrum constructor",
          "x" = "The following arguments are invalid: {invalid_names}",
          "i" = "The valid argument names are: {valid_names}"
          )
      )
    }
    stopifnot(length(drift_time) == length(intensity))
    stopifnot(is.null(baseline) || length(drift_time) == length(baseline))
    .Object@drift_time <- drift_time
    .Object@retention_time_idx <- retention_time_idx
    .Object@retention_time_s <- retention_time_s
    .Object@intensity <- intensity
    .Object@description <- description
    .Object@baseline <- baseline
    .Object
  })


#' @describeIn GCIMSSpectrum-class Friendly constructor
#' @param ... See the slots section
#' @return A GCIMSSpectrum object
#' @export
#' @examples
#' spec <- GCIMSSpectrum(drift_time = 1:10, intensity = c(1:5, 6:2))
GCIMSSpectrum <- function(...) {
  methods::new("GCIMSSpectrum", ...)
}
