#' @importMethodsFrom ProtGenerics smooth
#' @export
setMethod(
  "smooth",
  "GCIMSDataset",
  function(x, method = "savgol", dt_length_ms = 0.14, rt_length_s = 3, dt_order = 2, rt_order = 2) {

    delayed_op <- GCIMSDelayedOp(
      name = "smooth",
      fun = smooth,
      params = list(method = method, dt_length_ms = dt_length_ms, rt_length_s = rt_length_s, dt_order = dt_order, rt_order = rt_order)
    )
    x <- appendDelayedOp(x, delayed_op)
    # We recompute these, but  maybe we could just reset them to zero...
    x <- extract_RIC_and_TIS(x)
    invisible(x)
  }
)
