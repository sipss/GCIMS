#' Integrate peaks in a GCIMSDataset
#'
#' @param object The [GCIMSDataset] object, modified inline
#' @param peak_list A data frame with peak lists
#' @param integration_size_method Either fixed_size or free_size
#' @param rip_saturation_threshold The threshold
#'
#' @return A modified [GCIMSDataset] object
#' @export
setMethod(
  "integratePeaks",
  "GCIMSDataset",
  function(object, peak_list, integration_size_method = c("fixed_size", "free_size"), rip_saturation_threshold = 0.1) {
    integration_size_method <- match.arg(integration_size_method)
    delayed_op <- GCIMSDelayedOp(
      name = "integratePeaks",
      fun = integratePeaks,
      params = list(peak_list = peak_list, integration_size_method = integration_size_method, rip_saturation_threshold = rip_saturation_threshold),
      fun_extract = peaks,
      fun_aggregate = .findPeaks_fun_aggregate
    )
    object$appendDelayedOp(delayed_op)
    invisible(object)
  }
)
