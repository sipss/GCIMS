#' Remove baseline in data
#'
#' @param dir_in          Input directory. Where input data files are
#'   loaded from.
#' @param dir_out         Output directory. Where data files after baseline
#'   correction are stored.
#' @param samples         A vector. Set of samples to which remove the baseline
#'   (e.g.: c(1, 2, 3)).
#' @param peak_list       A peak list after containing information regarding the 2d peaks found in data.
#' @description Performs baseline correction on sample spectra (in drift time) and chromatograms (in retention time)
#' @details `gcims_removel_baseine()`corrects the baseline in drift time and afterwards, in retention time.
#' To do it, it cuts the axis in regions of equal length and looks for the local minima of each region. A simple
#' linear interpolation is performed to create a baseline of the same length as the number
#' of data points per axis.The region width is estimated differently in each of the the axes:
#' In drift time, it is computed as 12 times the standard deviation of the RIP width, while
#' in retention time as the third quantile of peak width distributions in chromatograms (obtained
#' from the peak table). Take into account that when removing the baseline in retention
#' time, the shape of the Reactant Ion Peak is severely affected.
#'
#'
#' @family Baseline removal functions
#' @return A set of S3 objects.
#' @export
#' @examples
#' \donttest{
#' dir_in <- system.file("extdata", package = "GCIMS")
#' peak_list <- readRDS(file.path(dir_in, "peak_list.rds"))
#' dir_out <- tempdir()
#' samples <- sample_num <- 3
#'
#' # Before:
#' gcims_view_sample(
#'   dir_in,
#'   sample_num = samples,
#'   rt_range = c(155, 175),
#'   dt_range = c(8.75, 9.1),
#'   transform = FALSE
#' )
#'
#' # Example of baseline removal:
#' gcims_remove_baseline(dir_in, dir_out, samples, peak_list)
#'
#' # After:
#' gcims_view_sample(
#'   dir_out,
#'   sample_num = samples,
#'   rt_range = c(155, 175),
#'   dt_range = c(8.75, 9.1),
#'   transform = FALSE
#' )
#'
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files))
#'}
#'

gcims_remove_baseline <- function(dir_in, dir_out, samples,
                                  peak_list){


  #-------------#
  #     MAIN    #
  #-------------#


  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

  for (i in samples) {
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(file.path(dir_in, aux_string)) #new
    aux <- as.matrix(aux_list$data$data_df)
    # correct drift time
    aux <- remove_bsln_td(aux)
    # correct retention time
    aux <- remove_bsln_tr(aux, peak_list)
    aux_list$data$data_df <- round(aux)
    saveRDS(aux_list, file = file.path(dir_out, paste0("M", i, ".rds")))
  }
}



remove_bsln_td <- function(aux, sig_mult = 12) {
  baseline <- estimate_baseline_td(aux, sig_mult)
  aux - baseline
}

remove_bsln_tr <- function(aux, peak_list){
  peak_rt_lengths <- peak_list$rt_max_idx - peak_list$rt_min_idx
  # Estimate region
  region_size <- stats::quantile(peak_rt_lengths, 0.75)
  baseline <- estimate_baseline_tr(aux, region_size)
  aux_corr <-  aux - baseline
  return(aux_corr)
}
