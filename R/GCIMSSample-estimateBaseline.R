#' Estimate the baseline of a GCIMS Sample using a connect local minima algorithm
#'
#' The baseline is estimated by connecting local minima and interpolating from those.
#' The local minima are identified as "the minima in each region of length x"
#' The length of the regions are estimated as `fwhm * a multiplier / 2.3482`. This
#' assumes it's several times
#'
#' @param object A [GCIMSSample] object
#' @param rt_length_s The length of the baseline region. It should be
#' comparable or longer than the peak width
#' @param dt_peak_fwhm_ms Full Width at Half Maximum in milliseconds. Used to
#' determine the length of the regions where local minima are searched.
#' @param dt_region_multiplier A multiplier to calculate the region
#' @return The modified [GCIMSSample]
#' @export
methods::setMethod(
  "estimateBaseline", "GCIMSSample",
  function(object, dt_peak_fwhm_ms, dt_region_multiplier, rt_length_s) {
    rt <- rtime(object)
    dt <- dtime(object)
    int <- intensity(object)
    dt_step_ms <- dt[2] - dt[1]
    dt_peak_fwhm_pts <- units_to_points(dt_peak_fwhm_ms, dt_step_ms)
    dt_basel <- estimate_baseline_td(int, sig_mult = dt_region_multiplier, peak_fwhm_pts = dt_peak_fwhm_pts)
    rt_step_s <- rt[2] - rt[1]
    rt_region_size_pts <- units_to_points(rt_length_s, rt_step_s)

    rt_basel <- estimate_baseline_tr(int - dt_basel, region_size = rt_region_size_pts)
    full_basel <- rt_basel + dt_basel
    baseline(object) <- full_basel
    object
  }
)

#' @describeIn estimateBaseline-GCIMSSample-method Get the baseline
#' @export
methods::setMethod(
  "baseline", "GCIMSSample",
  function(object) {
    if (is.null(object@baseline)) {
      rlang::abort("Please use estimateBaseline() first")
    }
    object@baseline
  }
)

#' @describeIn estimateBaseline-GCIMSSample-method Set the baseline
#' @param value A matrix with the sample baseline of the same dimensions as `dim(object)`
#' @export
methods::setMethod(
  "baseline<-", "GCIMSSample",
  function(object, value) {
    if (any(dim(value) != dim(object@data))) {
      rlang::abort("The baseline should be of the same length as the intensity")
    }
    object@baseline <- value
    object
  }
)
