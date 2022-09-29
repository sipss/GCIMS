#' Select Regions of Interest
#'
#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param dir_out         Output directory. Where ouput data files are stored.
#' @param samples         Numeric vector of integers. Identifies the set of
#'   samples to which peaks must be detected.
#' @param noise_level     Scalar number. The number of times the standard deviation
#'   above the noise level needed to detect a peak. IUPAC recommends `noise_level = 3` for detection.
#' @return A set of S3 objects. Additionally, a peak/roi list.
#' @family Peak detection
#' @description  `gcims_select_rois()` function looks for 2-dimensional peaks in data and
#' enclose each of them in a region or interest or ROI.
#'
#' @export
#' @examples
#' \donttest{
#' # Use BiocParallel library for parallelization
#' library(BiocParallel)
#' register(SnowParam(workers = 3, progressbar = TRUE, exportglobals = FALSE), default = TRUE)
#'
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' samples <- c(3,7)
#'
#' # Example of Peak Picking
#' peak_list <- gcims_select_rois(dir_in, dir_out, samples , noise_level = 3)
#' head(peak_list)
#'
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files))
#' }
#'
gcims_select_rois <- function(dir_in, dir_out, samples, noise_level){
  peak_lists <- gcims_batch_run(
    dir_in,
    dir_out,
    gcims_rois_selection_one,
    noise_level = noise_level,
    .batch_samples = samples,
    .batch_returns = function(x) {x$data$ROIs}
  )
  out <- dplyr::bind_rows(!!!peak_lists, .id = "SampleID")
  out$UniqueID <- sprintf("%s/%s", out$SampleID, out$PeakID)
  out
}



gcims_rois_selection_one <- function(x, noise_level, verbose = FALSE) {
  x$data$ROIs <- findPeaksImpl(
    drift_time = x$data$drift_time,
    retention_time = x$data$retention_time,
    int_mat = as.matrix(x$data$data_df),
    noise_level = noise_level,
    verbose = verbose
  )
  x
}


