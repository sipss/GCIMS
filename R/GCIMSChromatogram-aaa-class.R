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
    description = "character"
  )
)

methods::setMethod(
  "initialize", "GCIMSChromatogram",
  function(.Object, retention_time, intensity, drift_time_idx, drift_time_ms, description, ..., baseline = NULL) {
    dots <- list(...)
    if (length(dots) > 0) {
      if (!is.null(names(dots))) {
        invalid_names <- paste0(names(dots), collapse = ", ")
      }
      valid_names <- paste0(setdiff(names(formals()), c(".Object", "...")), collapse = ", ")
      rlang::abort(
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
    .Object
  })


#' @describeIn GCIMSChromatogram-class Friendly constructor
#' @param ... See the slots section below
#' @return A GCIMSChromatogram object
GCIMSChromatogram <- function(...) {
  methods::new("GCIMSChromatogram", ...)
}


#' @export
as.data.frame.GCIMSChromatogram <- function(x, ...) {
  out <- data.frame(
    retention_time_s = rtime(x),
    intensity = intensity(x)
  )
  if (!is.null(x@baseline)) {
    out$baseline <- baseline(x)
  }
  out
}


#' @export
plot.GCIMSChromatogram <- function(x, ...) {
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
