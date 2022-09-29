#' Peak detection for a GCIMSSample
#' @param object A [GCIMSSample] object
#' @param dt_length_ms Length (in ms) used in the Savitzky-Golay filter to compute the second derivative with respect to drift time.
#' @param rt_length_s Length (in s) used in the Savitzky-Golay filter to compute the second derivative with respect to retention time.
#' @inheritParams findPeaksImpl
#' @return The modified [GCIMSSample], with a peak list
#' @export
setMethod(
  "findPeaks",
  "GCIMSSample",
  function(object, noise_level = 3, verbose = FALSE, dt_length_ms = 0.14, rt_length_s = 3, iou_overlap_threshold = 0.2) {
    dt <- dtime(object)
    rt <- rtime(object)
    int_mat <- intensity(object)
    dt_length_pts <- units_to_points(dt_length_ms, dt[2] - dt[1], must_odd = TRUE)
    rt_length_pts <- units_to_points(rt_length_s, rt[2] - rt[1], must_odd = TRUE)
    peak_list <- findPeaksImpl(
      drift_time = dt,
      retention_time = rt,
      int_mat = int_mat,
      noise_level = noise_level,
      verbose = verbose,
      dt_length_pts = dt_length_pts,
      rt_length_pts = rt_length_pts,
      iou_overlap_threshold = iou_overlap_threshold
    )
    peaks(object) <- peak_list
    object
  }
)
