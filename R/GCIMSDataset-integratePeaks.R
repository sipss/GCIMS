#' Title
#'
#' @param object The [GCIMSDataset] object
#' @param peak_list
#' @param integration_size_method
#' @param rip_saturation_threshold
#'
#' @return
#' @export
#'
#' @examples
integratePeaks <- function(object, peak_list, integration_size_method = c("fixed_size", "free_size"), rip_saturation_threshold = 0.1) {
  integration_size_method <- match.arg(integration_size_method)
  GCIMS:::GCIMSDelayedOp(
    name = "integratePeaks",
    fun = integratePeaks,
    params = list(peak_list = peak_list, integration_size_method = integration_size_method, rip_saturation_threshold),
    fun_extract = peaks,
    fun_aggregate = .findPeaks_fun_aggregate
  )
}
