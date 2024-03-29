#' Estimate the baseline of a GCIMS Sample using a connect local minima algorithm
#'
#' @param object A [GCIMSDataset] object, modified in-place
#' @param rt_length_s The length of the baseline region. It should be
#' comparable or longer than the peak width
#' @param dt_peak_fwhm_ms Full Width at Half Maximum in milliseconds. Used to
#' determine the length of the regions where local minima are searched.
#' @param dt_region_multiplier A multiplier to calculate the region
#' @param remove A boolean, if TRUE it removes the baseline from the intensity
#' @return The modified [GCIMSDataset]
#' @export
setMethod(
  "estimateBaseline",
  "GCIMSDataset",
  function(object, dt_peak_fwhm_ms, dt_region_multiplier, rt_length_s, remove =TRUE) {

    delayed_op <- DelayedOperation(
      name = "estimateBaseline",
      fun = estimateBaseline,
      params = list(dt_peak_fwhm_ms = dt_peak_fwhm_ms, dt_region_multiplier = dt_region_multiplier, rt_length_s = rt_length_s, remove = remove)
    )
    object$appendDelayedOp(delayed_op)
    # We recompute these, but  maybe we could just reset them to zero...
    object <- extract_RIC_and_TIS(object)
    invisible(object)
  }
)
