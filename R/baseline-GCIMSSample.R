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
#' @param remove A boolean, if TRUE it removes the baseline from the intensity
#' @return The modified [GCIMSSample]
#' @export
methods::setMethod(
  "estimateBaseline", "GCIMSSample",
  function(object, dt_peak_fwhm_ms, dt_region_multiplier, rt_length_s, remove = TRUE) {
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
    if (remove){
      object@data <- int - full_basel
    }
    object
  }
)

#' @describeIn estimateBaseline-GCIMSSample-method Get the baseline
#' @inheritParams dt_rt_range_normalization
#' @param .error_if_missing A logical. If `TRUE`, raise error if baseline has not been estimated. If `FALSE` returns `NULL` instead.
#' @export
methods::setMethod(
  "baseline", "GCIMSSample",
  function(object, dt_range = NULL, rt_range = NULL, dt_idx = NULL, rt_idx = NULL, .error_if_missing = TRUE) {
    if (is.null(object@baseline)) {
      if (.error_if_missing) {
        cli_abort("Please use estimateBaseline() first")
      }
      return(NULL)
    }
    dt <- dtime(object)
    rt <- rtime(object)
    if (inherits(dt_range, "dt_rt_range_normalization")) {
      idx <- dt_range
    } else {
      idx <- dt_rt_range_normalization(dt, rt, dt_range, rt_range, dt_idx, rt_idx)
    }
    dt_idx <- idx[["dt_logical"]]
    rt_idx <- idx[["rt_logical"]]
    out <- object@baseline[dt_idx, rt_idx, drop = FALSE]
    dimnames(out) <- list(dt_ms = dt[dt_idx], rt_s = rt[rt_idx])
    out
  }
)

#' @describeIn estimateBaseline-GCIMSSample-method Set the baseline
#' @param value A matrix with the sample baseline of the same dimensions as `dim(object)`
#' @export
methods::setMethod(
  "baseline<-", "GCIMSSample",
  function(object, value) {
    if (any(dim(value) != dim(object@data))) {
      cli_abort("The baseline should be of the same length as the intensity")
    }
    object@baseline <- value
    object
  }
)
