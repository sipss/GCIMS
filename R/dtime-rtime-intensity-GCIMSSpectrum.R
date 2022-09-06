#' @describeIn GCIMSSpectrum-class Get the drift time vector
#' @export
setMethod("dtime", "GCIMSSpectrum", function(object) object@drift_time)

#' @describeIn GCIMSSpectrum-class Get the retention time where this spectrum was extracted
#' @export
setMethod("rtime", "GCIMSSpectrum", function(object) object@retention_time_s)

#' @describeIn GCIMSSpectrum-class Get the intensity matrix
#' @param object A [GCIMSSpectrum] object
#' @inheritParams dt_rt_range_normalization
setMethod("intensity", "GCIMSSpectrum", function(object, dt_range = NULL, dt_idx = NULL) {
  dt <- dtime(object)
  idx <- dt_rt_range_normalization(dt = dt, dt_range = dt_range, dt_idx = dt_idx)
  out <- object@intensity[idx$dt_logical]
  names(out) <- dt[idx$dt_logical]
  out
})
