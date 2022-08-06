#' @describeIn GCIMSChromatogram-class Get the drift time vector
#' @export
setMethod("dtime", "GCIMSChromatogram", function(object) object@drift_time_ms)

#' @describeIn GCIMSChromatogram-class Get the retention time where this spectrum was extracted
#' @import ProtGenerics
#' @export
setMethod("rtime", "GCIMSChromatogram", function(object) object@retention_time)

#' @describeIn GCIMSChromatogram-class Get the intensity matrix
#' @import ProtGenerics
#' @param rt_range The minimum and maximum retention times to extract (length 2 vector)
#' @param rt_idx A numeric vector with the retention time indices to extract (or a logical vector of the length of retention time)
setMethod("intensity", "GCIMSChromatogram", function(object, rt_range = NULL, rt_idx = NULL) {
  rt <- rtime(object)
  idx <- dt_rt_range_normalization(rt = rt, rt_range = rt_range, rt_idx = rt_idx)
  out <- object@intensity[idx$rt_logical]
  names(out) <- rt[idx$rt_logical]
  out
})
