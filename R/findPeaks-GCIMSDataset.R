#' Peak detection on the GCIMS dataset
#' @param object A [GCIMSDataset] object, modified in-place
#' @inheritDotParams findPeaksImpl -drift_time -retention_time -int_mat
#' @return The modified [GCIMSDataset], with a peak list
#' @export
setMethod(
  "findPeaks",
  "GCIMSDataset",
  function(object, ...) {
    params <- list(...)
    delayed_op <- GCIMSDelayedOp(
      name = "findPeaks",
      fun = findPeaks,
      params = params,
      fun_extract = peaks,
      fun_aggregate = .findPeaks_fun_aggregate
    )
    object$appendDelayedOp(delayed_op)
    invisible(object)
  }
)

.findPeaks_fun_aggregate <- function(ds, objs) {
  p <- dplyr::bind_rows(purrr::map(objs, as.data.frame))
  peaks(ds) <- p
  ds
}
