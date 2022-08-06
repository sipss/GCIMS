methods::setMethod(
  "smooth", "GCIMSSample",
  function(x, method = "savgol", dt_length_ms, rt_length_s, dt_order = 2, rt_order = 2){
    if (method != "savgol") {
      stop("Smoothing only implemented with the Savitzky-Golay filter")
    }

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
