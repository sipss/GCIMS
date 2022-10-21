#' Peak detection for a GCIMSSample
#' @param object A [GCIMSSample] object
#' @inheritDotParams findPeaksImpl -drift_time -retention_time -int_mat
#' @return The modified [GCIMSSample], with a peak list
#' @export
setMethod(
  "findPeaks",
  "GCIMSSample",
  function(object, ...) {
    dt <- dtime(object)
    rt <- rtime(object)
    int_mat <- intensity(object)
    peak_list_and_debug_info <- findPeaksImpl(
      drift_time = dt,
      retention_time = rt,
      int_mat = int_mat,
      ...
    )
    object@peaks_debug_info <- peak_list_and_debug_info$debug_info
    peaks(object) <- peak_list_and_debug_info$peak_list
    object
  }
)
