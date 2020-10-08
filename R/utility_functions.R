#' Converts Sample Matrices to Feature Vectors


#' @param dir_in          The input directory.
#' @param dir_out         The output directory.
#' @param samples         The set of samples to be processed.
#' @return A set of feature vectors
#' @family Utility functions
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

gcims_unfold <- function(dir_in, dir_out, samples){
  print(" ")
  print("  //////////////////////////////////")
  print("  / Creation of the Features Vector/")
  print("  //////////////////////////////////")
  print(" ")

  setwd(dir_in)
  m = 0
  for (i in samples){
    m <- m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux <- readRDS(aux_string)
    M <- as.vector(t(aux))
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }
}
