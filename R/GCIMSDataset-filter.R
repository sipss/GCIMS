#' @importMethodsFrom ProtGenerics filterRt
#' @export
setMethod(
  "filterRt",
  "GCIMSDataset",
  function(object, rt_range) {
    delayed_op <- GCIMSDelayedOp(
      name = "filterRt",
      fun = function(x, rt_range) {
        filterRt(x, rt_range = rt_range)
      },
      params = list(rt_range = rt_range)
    )
    object <- appendDelayedOp(object, delayed_op)
    # We recompute these, but  maybe we could just reset them to zero...
    object <- extract_dtime_rtime(object)
    object <- extract_RIC_and_TIS(object)
    invisible(object)
  }
)

#' @export
setMethod(
  "filterDt",
  "GCIMSDataset",
  function(object, dt_range) {
    delayed_op <- GCIMSDelayedOp(
      name = "filterDt",
      fun = function(x, dt_range) {
        filterDt(x, dt_range = dt_range)
      },
      params = list(dt_range = dt_range)
    )
    object <- appendDelayedOp(object, delayed_op)
    # We recompute these, but  maybe we could just reset them to zero...
    object <- extract_dtime_rtime(object)
    object <- extract_RIC_and_TIS(object)
    invisible(object)
  }
)

