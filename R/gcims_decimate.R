#' Data Smoothing using Savitzky-Golay filters


#' @param dir_in           The input directory.
#' @param dir_out          The output directory.
#' @param samples          The set of samples to be processed.
#' @param filter_length    Numerical. Length of the filter.
#' @param polynomial_order Numerical. Order of the polynomial.
#' @return A filtered  gcims dataset.
#' @family Smoothing functions
#' @export
#' @importFrom signal resample
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }


gcims_decimate <- function(dir_in, dir_out, samples){
  print(" ")
  print("  ///////////////////////")
  print("  / Samples Decimation /")
  print("  //////////////////////")
  print(" ")

  setwd(dir_in)
  # bin_ratiort <- 5
  # bin_ratiodt <- 2
  m = 0
  for (i in samples){
    m <- m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(aux_string)
    aux <- aux_list$data$data_df
    dtreduced <- apply(aux, 1, function(x) resample(x, p = 12, q = 32))
    rtreduced <- apply(dtreduced, 2, function(x) resample(x, p = 6, q = 32))
    resampleddt <- resample(as.matrix(aux_list$data$drift_time), p = 12, q = 32)
    resampledrt <- resample(as.matrix(aux_list$data$retention_time), p = 6, q = 32)
    aux_list$data$data_df <- rtreduced
    aux_list$data$drift_time <- resampleddt
    aux_list$data$retention_time <- resampledrt
    M <- aux_list
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }
}



