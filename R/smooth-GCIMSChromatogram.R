#' Smoothing a GCIMS chromatogram using a Savitzky-Golay filter
#' @param x A [GCIMSChromatogram] object
#' @param rt_length_s The length of the filter in retention time (in s)
#' @param rt_order The order of the filter in retention time
#' @return The modified [GCIMSChromatogram]
#' @family GCIMSChromatogram
#' @export
methods::setMethod(
  "smooth", "GCIMSChromatogram",
  function(x, rt_length_s = 3, rt_order = 2L) {
    rt <- rtime(x)
    rt_length_pts <- units_to_points(rt_length_s, rt[2] - rt[1], must_odd = TRUE)
    if (rt_length_pts >= 1L) {
      x@intensity <- sgolayfilt(x@intensity, n = rt_length_pts, p = rt_order)
    }
    x
  }
)
