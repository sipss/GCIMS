#' GCIMSChromatogram class
#'
#' @description
#' GCIMSChromatogram is an S4 class to store a GCIMS Chromatogram It can be a single chromatogram
#' or the aggregation of several chromatograms.
#'
#' @slot retention_time A numeric vector with retention times
#' @slot intensity A numeric vector with the corresponding intensities
#' @slot baseline A numeric vector of the same length as `intensity` with the corresponding baseline,
#' or `NULL` if not set. Use [estimateBaseline()] to estimate it, [baseline()] to directly access it.
#' @slot drift_time_idx The index or indices used to get the intensity
#' @slot drift_time_ms The drift times corresponding to `drift_time_idx`.
#' @slot description A string with a description (used as plot title, useful e.g. to know the sample it came from)
#' @slot peaks A data frame with peaks, use [findPeaks()] to fill it, or [peaks()] to set/get it
#' @slot peaks_debug_info A list with arbitrary debug infor from [findPeaks()]
#'
#' @export
methods::setClass(
  Class = "GCIMSChromatogram",
  slots = c(
    retention_time = "numeric",
    intensity = "numeric",
    baseline = "numericOrNULL",
    drift_time_idx = "numeric",
    drift_time_ms = "numeric",
    description = "character",
    peaks = "DataFrameOrNULL",
    peaks_debug_info = "listOrNULL" # arbitrary debug info from findPeaks()
  )
)

methods::setMethod(
  "initialize", "GCIMSChromatogram",
  function(.Object, retention_time, intensity, drift_time_idx, drift_time_ms, description, ..., baseline = NULL, peaks = NULL, peaks_debug_info = NULL) {
    dots <- list(...)
    if (length(dots) > 0) {
      if (!is.null(names(dots))) {
        invalid_names <- paste0(names(dots), collapse = ", ")
      }
      valid_names <- paste0(setdiff(names(formals()), c(".Object", "...")), collapse = ", ")
      abort(
        message = c(
          "Invalid argument to GCIMSChromatogram constructor",
          "x" = paste0("The following arguments are invalid: ", invalid_names),
          "i" = paste0("The valid argument names are: ", valid_names)
        )
      )
    }
    stopifnot(length(retention_time) == length(intensity))
    stopifnot(is.null(baseline) || length(retention_time) == length(baseline))
    .Object@retention_time <- retention_time
    .Object@drift_time_idx <- drift_time_idx
    .Object@drift_time_ms <- drift_time_ms
    .Object@intensity <- intensity
    .Object@description <- description
    .Object@baseline <- baseline
    .Object@peaks <- peaks
    .Object@peaks_debug_info <- peaks_debug_info
    .Object
  })


#' @describeIn GCIMSChromatogram-class Friendly constructor
#' @param ... See the slots section below
#' @return A GCIMSChromatogram object
GCIMSChromatogram <- function(...) {
  methods::new("GCIMSChromatogram", ...)
}

