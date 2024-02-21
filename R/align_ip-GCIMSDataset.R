#' Align a GCIMS dataset to the injection point
#'
#' Aligns all the samples to their injection point in retention time
#'
#' @param object A [GCIMSDataset] object, modified in-place
#' @return The modified [GCIMSDataset]
#' @export
methods::setMethod(
  "align_ip",
  "GCIMSDataset",
  function(object){
    rics <- GCIMS::getRIC(object)
    mins <- apply(rics, 1,which.min)
    rt_ref <- GCIMS::rtime(object, which.min(mins))
    rt_ref <- rt_ref[1 : (length(rt_ref) - (max(mins) - min(mins)))]
    min_start <- min(mins) - 1
    delayed_op <- DelayedOperation(
      name = "align_injection_point",
      fun = align_ip,
      params = list(min_start = min_start, rt_ref = rt_ref)
    )
    object$appendDelayedOp(delayed_op)
    object$extract_dtime_rtime()
    object <- GCIMS:::extract_RIC_and_TIS(object)
    invisible(object)
  }
)




