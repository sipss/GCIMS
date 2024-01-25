#' Align a GCIMS Sample to the injection point
#'
#' Aligns a sample to the reference injection point
#'
#' @param object A [GCIMSSample] object, modified in-place
#' @param min_start The points that will be used before the injection point
#' @param rt_ref Retention time vector to be used
#' @return The modified [GCIMSSample]
#' @export
methods::setMethod(
  "align_ip",
  "GCIMSSample",
  function(object, min_start, rt_ref){
    ric <- GCIMS::getRIC(object)
    injection_point <- which.min(ric)
    object@retention_time <- rt_ref
    object@data <- object@data[, (injection_point - min_start):((injection_point - min_start)+length(rt_ref)-1)]
    return(object)
  }
)
