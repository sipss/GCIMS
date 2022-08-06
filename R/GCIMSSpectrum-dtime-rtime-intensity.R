#' @describeIn GCIMSSpectrum-class Get the drift time vector
#' @export
setMethod("dtime", "GCIMSSpectrum", function(object) object@drift_time)

#' @describeIn GCIMSSpectrum-class Get the retention time where this spectrum was extracted
#' @import ProtGenerics
#' @export
setMethod("rtime", "GCIMSSpectrum", function(object) object@retention_time_s)

#' @describeIn GCIMSSpectrum-class Get the intensity matrix
#' @import ProtGenerics
#' @param dt_range The minimum and maximum drift times to extract (length 2 vector)
#' @param dt_idx A numeric vector with the drift time indices to extract (or a logical vector of the length of drift time)
setMethod("intensity", "GCIMSSpectrum", function(object, dt_range = NULL, dt_idx = NULL) {
  dt <- dtime(object)
  if (!is.null(dt_range) && !is.null(dt_idx)) {
    rlang::abort("Please provide either dt_range or dt_idx, but not both in the same call")
  }
  if (!is.null(dt_range)) {
    dtmin <- min(dt_range)
    dtmax <- max(dt_range)
    dt_idx <- dt >= dtmin & dt <= dtmax
  } else if (!is.null(dt_idx)) {
    if (is.logical(dt_idx) && length(dt_idx) != length(dt)) {
      rlang::abort(
        sprintf(
          "dt_idx is a logical vector of length %d and it should be of length %d",
          length(dt_idx),
          length(dt)
        )
      )
    } else if (is.numeric(dt_idx) && (min(dt_idx) < 1 || max(dt_idx) > length(dt))) {
      rlang::abort("dt_idx out of range")
    }
  } else {
    dt_idx <- seq_along(dt)
  }


  out <- object@intensity[dt_idx]
  names(out) <- dt[dt_idx]
  out
})
