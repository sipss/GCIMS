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
#' @slot peaks_debug_info A list with arbitrary debug information from [findPeaks()]
#'
#' @export
#' @family GCIMSChromatogram
methods::setClass(
  Class = "GCIMSChromatogram",
  slots = c(
    retention_time = "numeric",
    intensity = "numeric",
    baseline = "numericOrNULL",
    drift_time_idx = "integer",
    drift_time_ms = "numeric",
    description = "character",
    peaks = "DataFrameOrNULL",
    peaks_debug_info = "listOrNULL" # arbitrary debug info from findPeaks()
  )
)

methods::setMethod(
  "initialize", "GCIMSChromatogram",
  function(.Object, retention_time, intensity, drift_time_idx = NA_integer_,
           drift_time_ms = NA_real_, description = "", baseline = NULL,
           peaks = NULL, peaks_debug_info = NULL) {
    stopifnot(length(retention_time) == length(intensity))
    stopifnot(is.null(baseline) || length(retention_time) == length(baseline))
    .Object@retention_time <- retention_time
    .Object@drift_time_idx <- as.integer(drift_time_idx)
    .Object@drift_time_ms <- drift_time_ms
    .Object@intensity <- intensity
    .Object@description <- description
    .Object@baseline <- baseline
    .Object@peaks <- peaks
    .Object@peaks_debug_info <- peaks_debug_info
    .Object
  })


#' Create a [GCIMSChromatogram-class] object
#'
#' @param retention_time A numeric vector with retention times
#' @param intensity A numeric vector with the corresponding intensities
#' @param baseline A numeric vector of the same length as `intensity` with the corresponding baseline,
#' or `NULL` if not set. Use [estimateBaseline()] to estimate it, [baseline()] to directly access it.
#' @param drift_time_idx The index or indices used to get the intensity
#' @param drift_time_ms The drift times corresponding to `drift_time_idx`.
#' @param description A string with a description (used as plot title, useful e.g. to know the sample it came from)
#' @param peaks A data frame with peaks, use [findPeaks()] to fill it, or [peaks()] to set/get it
#' @param peaks_debug_info A list with arbitrary debug information from [findPeaks()]
#' @return A [GCIMSChromatogram-class] object
#' @examples
#' GCIMSChromatogram(
#'   retention_time = seq(from = 0, to = 10, length.out = 200),
#'   intensity = 1:200
#' )
#' @export
#' @family GCIMSChromatogram
GCIMSChromatogram <- function(
    retention_time, intensity, drift_time_idx = NA_integer_,
    drift_time_ms = NA_real_, description = "",
    baseline = NULL, peaks = NULL, peaks_debug_info = NULL) {
  methods::new("GCIMSChromatogram", retention_time, intensity, drift_time_idx,
               drift_time_ms, description, baseline = NULL, peaks = NULL,
               peaks_debug_info = NULL)
}

