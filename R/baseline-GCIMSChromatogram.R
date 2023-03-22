#' Estimate the baseline of a GCIMS Chromatogram using a connect local minima algorithm
#'
#' The baseline is estimated by connecting local minima and interpolating from those.
#' The local minima are identified as "the minima in each region of length x"
#' The length of the regions are given in seconds in the `region_s` parameter.
#'
#' @param object A [GCIMSChromatogram] object
#' @param rt_length_s The length of the baseline region. It should be
#' comparable or longer than the peak width
#' @return The modified [GCIMSChromatogram]
#' @family GCIMSChromatogram
#' @export
methods::setMethod(
  "estimateBaseline", "GCIMSChromatogram",
  function(object, rt_length_s) {
    rt <- rtime(object)
    int <- intensity(object)
    rt_step_s <- rt[2] - rt[1]
    region_size_pts <- units_to_points(rt_length_s, rt_step_s)
    # estimate_baseline_td expects a matrix:
    int <- matrix(int, nrow = 1L)
    basel <- estimate_baseline_tr(int, region_size = region_size_pts)
    baseline(object) <- as.numeric(basel)
    object
  }
)

#' @describeIn estimateBaseline-GCIMSChromatogram-method Get the baseline
#' @inheritParams dt_rt_range_normalization
#' @param .error_if_missing A logical. If `TRUE` (default) give an error if baseline is not estimated. Returns `NULL` otherwise.
#' @family GCIMSChromatogram
#' @export
methods::setMethod(
  "baseline", "GCIMSChromatogram",
  function(object, rt_range = NULL, rt_idx = NULL, .error_if_missing = TRUE) {
    if (is.null(object@baseline)) {
      if (.error_if_missing) {
        cli_abort("Please use {.code estimateBaseline()} first")
      }
      return(NULL)
    }
    rt <- rtime(object)
    idx <- dt_rt_range_normalization(rt = rt, rt_range = rt_range, rt_idx = rt_idx)
    out <- object@baseline[idx$rt_logical]
    names(out) <- rt[idx$rt_logical]
    out
  }
)


#' @describeIn estimateBaseline-GCIMSChromatogram-method Set the baseline
#' @param value A vector with a baseline of the same length as `intensity(object)`
#' @family GCIMSChromatogram
#' @export
methods::setMethod(
  "baseline<-", "GCIMSChromatogram",
  function(object, value) {
    if (length(value) != length(object@intensity)) {
      cli_abort("The baseline (length {length(value)}) should be of the same length as the intensity (length {length(object@intensity)})")
    }
    object@baseline <- value
    object
  }
)
