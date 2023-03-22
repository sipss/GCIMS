#' Estimate the baseline of a GCIMS Spectrum using a connect local minima algorithm
#'
#' The baseline is estimated by connecting local minima and interpolating from those.
#' The local minima are identified as "the minima in each region of length x"
#' The length of the regions are estimated as `fwhm * a multiplier / 2.3482`. This
#' assumes it's several times
#'
#' @param object A [GCIMSSpectrum] object
#' @param dt_peak_fwhm_ms Full Width at Half Maximum in milliseconds. Used to
#' determine the length of the regions where local minima are searched.
#' @param dt_region_multiplier A multiplier to calculate the region
#' @return The modified [GCIMSSpectrum]
#' @export
methods::setMethod(
  "estimateBaseline", "GCIMSSpectrum",
  function(object, dt_peak_fwhm_ms, dt_region_multiplier = 12) {
    dt <- dtime(object)
    int <- intensity(object)
    dt_step_ms <- dt[2] - dt[1]
    dt_peak_fwhm_pts <- units_to_points(dt_peak_fwhm_ms, dt_step_ms)
    # estimate_baseline_td expects a matrix:
    int <- matrix(int, ncol = 1L)
    basel <- estimate_baseline_td(int, sig_mult = dt_region_multiplier, peak_fwhm_pts = dt_peak_fwhm_pts)
    baseline(object) <- as.numeric(basel)
    object
  }
)

#' @describeIn estimateBaseline-GCIMSSpectrum-method Get the baseline
#' @inheritParams dt_rt_range_normalization
#' @param .error_if_missing A logical. If `TRUE` (default) give an error if baseline is not estimated. Returns `NULL` otherwise.
#' @export
methods::setMethod(
  "baseline", "GCIMSSpectrum",
  function(object, dt_range = NULL, dt_idx = NULL, .error_if_missing = TRUE) {
    if (is.null(object@baseline)) {
      if (.error_if_missing) {
        cli_abort("Please use estimateBaseline() first")
      }
      return(NULL)
    }
    dt <- dtime(object)
    idx <- dt_rt_range_normalization(dt = dt, dt_range = dt_range, dt_idx = dt_idx)
    out <- object@baseline[idx$dt_logical]
    names(out) <- dt[idx$dt_logical]
    out
  }
)

#' @describeIn estimateBaseline-GCIMSSpectrum-method Set the baseline
#' @param value A vector with a baseline of the same length as `intensity(object)`
#' @export
methods::setMethod(
  "baseline<-", "GCIMSSpectrum",
  function(object, value) {
    if (length(value) != length(object@intensity)) {
      cli_abort("The baseline should be of the same length as the intensity")
    }
    object@baseline <- value
    object
  }
)
