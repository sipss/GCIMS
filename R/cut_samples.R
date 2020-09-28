#' Cuts retention time and drift time
#'
#' @export
#'
cut_samples <- function(loaded_samples, drift_cut=NULL, ret_cut=NULL) {
  loaded_samples_cut <- loaded_samples
  if (is.null(drift_cut)) {
    drift_cut <- 1:dim(loaded_samples$data)[1]
  }
  if (is.null(ret_cut)) {
    ret_cut <- 1:dim(loaded_samples$data)[2]
  }
  loaded_samples_cut$data <- loaded_samples_cut$data[drift_cut, ret_cut, ]
  loaded_samples_cut$drift_time <- loaded_samples_cut$drift_time[drift_cut]
  loaded_samples_cut$ret_time <- loaded_samples_cut$ret_time[ret_cut]
  return(loaded_samples_cut)
}
