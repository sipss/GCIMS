#' Smoothing a GCIMS sample using a Savitzky-Golay filter
#' @param x A [GCIMSSample] object
#' @param dt_length_ms the length of the filter in drift time (in ms)
#' @param dt_order The order of the filter in drift time
#' @param rt_length_s The length of the filter in retention time (in s)
#' @param rt_order The order of the filter in retention time
#' @return The modified [GCIMSSample]
#' @importMethodsFrom ProtGenerics smooth
#' @export
methods::setMethod(
  "smooth", "GCIMSSample",
  function(x, dt_length_ms, rt_length_s, dt_order = 2, rt_order = 2){
    dt <- dtime(x)
    rt <- rtime(x)

    dt_length_pts <- units_to_points(dt_length_ms, dt[2] - dt[1], must_odd = TRUE)
    rt_length_pts <- units_to_points(rt_length_s, rt[2] - rt[1], must_odd = TRUE)

    filter_rt <- signal::sgolay(p = rt_order, n = rt_length_pts)
    filter_dt <- signal::sgolay(p = dt_order, n = dt_length_pts)

    mat <- x@data
    for (i in seq_len(ncol(mat))) {
      mat[, i] <- signal::sgolayfilt(mat[, i], filter_rt)
    }
    for (i in seq_len(nrow(mat))) {
      mat[i, ] <- signal::sgolayfilt(mat[i, ], filter_dt)
    }
    x@data <- mat
    x
  })
