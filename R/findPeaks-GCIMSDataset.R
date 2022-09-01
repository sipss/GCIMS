#' Peak detection on the GCIMS dataset
#' @param object A [GCIMSDataset] object, modified in-place
#' @inheritParams findPeaks,GCIMSSample-method
#' @return The modified [GCIMSDataset], with a peak list
#' @export
setMethod(
  "findPeaks",
  "GCIMSDataset",
  function(object, noise_level = 3, verbose = FALSE, dt_length_ms = 0.14, rt_length_s = 3, iou_overlap_threshold = 0.2) {
    delayed_op <- GCIMSDelayedOp(
      name = "findPeaks",
      fun = findPeaks,
      params = list(noise_level = noise_level, verbose = verbose, dt_length_ms = dt_length_ms, rt_length_s = rt_length_s, iou_overlap_threshold = iou_overlap_threshold),
      fun_extract = peaks,
      fun_aggregate = .findPeaks_fun_aggregate
    )
    object <- appendDelayedOp(object, delayed_op)
    invisible(object)
  }
)

.findPeaks_fun_aggregate <- function(ds, objs) {
  p <- dplyr::bind_rows(purrr::map(objs, as.data.frame))
  peaks(ds) <- p
  ds
}
