#' Topographical plot of a GCIMSSample object
#'
#' @param object A [GCIMSSample] object
#' @param dt_range A numeric vector of length 2 with the drift time range to plot (in milliseconds)
#' @param rt_range A numeric vector of length 2 with the retention time range to plot (in seconds)
#' @param ... Ignored
#' @param remove_baseline Set to `TRUE` to subtract the estimated baseline first
#' @keywords internal
#' @return A plot of the GCIMSSample
#' @examples
#' dummy_obj <-GCIMSSample(
#'   drift_time = 1:2,
#'   retention_time = 1:3,
#'   data = matrix(1:6, nrow = 2, ncol = 3),
#'   gc_column = "Optional column name",
#'   drift_gas = "nitrogen",
#'   drift_tube_length = 98.0 # in mm
#' )
#' plotRaw(dummy_obj)
#' @export
plotRaw <- function(object, dt_range = NULL, rt_range = NULL, ..., remove_baseline = FALSE) {
  # FIXME: Before releasing remove plotRaw
  #rlang::warn(
  #  message = "plotRaw(...) is deprecated. Please use plot(..) instead (rename)"
  #)
  plot(object, dt_range = dt_range, rt_range = rt_range, ..., remove_baseline = remove_baseline)
}
