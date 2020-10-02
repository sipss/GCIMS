
#' Cut samples in a retention time - drift time rectangle


#' @param dir_in          The input directory.
#' @param dir_out         The output directory.
#' @param samples         The set of samples to be processed.
#' @param ret_cut         A  vector of two components. Beginning and end
#'                        of the retencion time cut.
#' @param drift_cut       A  vector of two components. Beginning and end
#'                        of the drift time cut.
#' @return An cut gcims dataset.
#' @family Alignment functions
#' @export
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }

gcims_cut_samples <- function(dir_in, dir_out, samples, ret_cut, dritf_cut){


  print(" ")
  print("  /////////////////////")
  print(" /    Cut Samples    /")
  print("/////////////////////")
  print(" ")


  setwd(dir_in)
  m <- -1
  for (i in c(0,samples)){
    m <- m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux <- readRDS(aux_string)
    M <- cut_samples(aux, ret_cut, drift_cut)
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)

  }

}



#' Cuts retention time and drift time
#'
#' @param loaded_sample   Current sample to be cut.
#' @param ret_cut         A  vector of two components. Beginning and end
#'                        of the retencion time cut.
#' @param drift_cut       A  vector of two components. Beginning and end
#'                        of the drift time cut.
#' @return An cut matrix.
#' @family Alignment functions
#' @export
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }
#'
#'
cut_samples <- function(loaded_sample, ret_cut = NULL, drift_cut = NULL) {
  loaded_sample_cut <- loaded_sample
  if (is.null(drift_cut)) {
    drift_cut <- 1:dim(loaded_sample)[2]
  }
  if (is.null(ret_cut)) {
    ret_cut <- 1:dim(loaded_sample)[1]
  }
  loaded_samples_cut <- loaded_sample_cut[ret_cut, drift_cut]
  #loaded_samples_cut$data <- loaded_sample_cut$data[ret_cut, drift_cut]
  # loaded_samples_cut$drift_time <- loaded_samples_cut$drift_time[drift_cut]
  # loaded_samples_cut$ret_time <- loaded_samples_cut$ret_time[ret_cut]
  return(loaded_samples_cut)
}
