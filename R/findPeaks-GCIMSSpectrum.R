#' Peak detection for a GCIMSSpectrum
#' @param object A [GCIMSSpectrum] object
#' @inheritDotParams findPeaksImpl1D -x -y
#' @return The modified [GCIMSSpectrum], with a peak list
#' @family GCIMSSpectrum
#' @export
setMethod(
  "findPeaks",
  "GCIMSSpectrum",
  function(object, ...) {
    dt <- dtime(object)
    intens <- intensity(object)
    peak_list_and_debug_info <- findPeaksImpl1D(
      x = dt,
      y = intens,
      ...
    )
    object@peaks_debug_info <- peak_list_and_debug_info$debug_info
    peaks(object) <- peak_list_and_debug_info$peak_list
    object
  }
)
