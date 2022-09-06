#' Decimates a GCIMS sample
#'
#' This method assumes that the sample has been low-pass filtered to avoid aliasing issues
#'
#' @param object A [GCIMSSample] object
#' @param rt_factor Keep one every `rt_factor` measurement points in retention time
#' @param dt_factor Keep one every `dt_factor` measurement points in drift time
#' @return The modified [GCIMSSample]
#' @export
methods::setMethod(
  "decimate", "GCIMSSample",
  function(object, dt_factor = 1L, rt_factor = 1L) {
    # FIXME: Don't we need some lowpass filtering to prevent aliasing issues?
    rt_factor <- as.integer(rt_factor)
    dt_factor <- as.integer(dt_factor)
    if (rt_factor <= 0L || dt_factor <= 0L) {
      abort("decimation factors must be positive integers")
    }
    if (rt_factor == 1L && dt_factor == 1L) {
      return(object)
    }
    dt <- dtime(object)
    rt <- rtime(object)
    dt_idx <- seq(from = 1L, to = length(dt), by = dt_factor)
    rt_idx <- seq(from = 1L, to = length(rt), by = rt_factor)
    subset(object, dt_idx = dt_idx, rt_idx = rt_idx)
  })
