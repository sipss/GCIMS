#' Smoothing a GCIMS Spectrum using a Savitzky-Golay filter
#' @param x A [GCIMSSpectrum] object
#' @param dt_length_ms the length of the filter in drift time (in ms)
#' @param dt_order The order of the filter in drift time
#' @return The modified [GCIMSSpectrum]
#' @importMethodsFrom ProtGenerics smooth
#' @export
methods::setMethod(
  "smooth", "GCIMSSpectrum",
  function(x, dt_length_ms, dt_order = 2){
    dt <- dtime(x)
    dt_length_pts <- units_to_points(dt_length_ms, dt[2] - dt[1], must_odd = TRUE)
    x@intensity <- signal::sgolayfilt(x@intensity, n = dt_length_pts, p = dt_order)
    x
  })
