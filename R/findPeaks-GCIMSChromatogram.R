#' Peak detection for a GCIMSChromatogram
#' @param object A [GCIMSChromatogram] object
#' @inheritDotParams findPeaksImpl1D -x -y
#' @return The modified [GCIMSChromatogram], with a peak list
#' @export
setMethod(
  "findPeaks",
  "GCIMSChromatogram",
  function(object, ...) {
    rt <- rtime(object)
    intens <- intensity(object)
    peak_list_and_debug_info <- findPeaksImpl1D(
      x = rt,
      y = intens,
      ...
    )
    object@peaks_debug_info <- peak_list_and_debug_info$debug_info
    peaks(object) <- peak_list_and_debug_info$peak_list
    object
  }
)
