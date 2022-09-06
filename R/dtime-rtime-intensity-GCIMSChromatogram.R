#' Get the drift time of the chromatogram
#'
#' @param object A GCIMSChromatogram
#' @return The drift time where this chromatogram was extracted from (in ms)
#' @export
setMethod("dtime", "GCIMSChromatogram", function(object) object@drift_time_ms)

#' Get the retention time vector
#'
#' @param object A GCIMSChromatogram
#' @return The retention time vector (in s)
#' @export
#' @importMethodsFrom ProtGenerics rtime
#' @export
setMethod("rtime", "GCIMSChromatogram", function(object) object@retention_time)

#' Get the intensity vector
#' @param object A GCIMSChromatogram object
#' @inheritParams dt_rt_range_normalization
#' @return The retention intesity vector
#' @importMethodsFrom ProtGenerics intensity
#' @export
setMethod("intensity", "GCIMSChromatogram", function(object, rt_range = NULL, rt_idx = NULL) {
  rt <- rtime(object)
  idx <- dt_rt_range_normalization(rt = rt, rt_range = rt_range, rt_idx = rt_idx)
  out <- object@intensity[idx$rt_logical]
  names(out) <- rt[idx$rt_logical]
  out
})
