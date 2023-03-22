#' Decimate a GCIMS dataset keeping 1 out of n points
#' @param object A [GCIMSDataset] object, modified in-place
#' @param dt_factor Keep one every `dt_factor` measurement points in drift time
#' @param rt_factor Keep one every `rt_factor` measurement points in retention time
#' @return The modified [GCIMSDataset]
#' @export
setMethod(
  "decimate",
  "GCIMSDataset",
  function(object, dt_factor = 1L, rt_factor = 1L) {
    rt_factor <- as.integer(rt_factor)
    dt_factor <- as.integer(dt_factor)
    if (rt_factor <= 0L || dt_factor <= 0L) {
      cli_abort("decimation factors must be positive integers")
    }
    if (rt_factor == 1L && dt_factor == 1L) {
      return(object)
    }
    delayed_op <- GCIMSDelayedOp(
      name = "decimate",
      fun = decimate,
      params = list(dt_factor = dt_factor, rt_factor = rt_factor)
    )
    object <- appendDelayedOp(object, delayed_op)
    # We recompute these, but  maybe we could just reset them to zero...
    object <- extract_dtime_rtime(object)
    object <- extract_RIC_and_TIS(object)
    invisible(object)
  }
)

