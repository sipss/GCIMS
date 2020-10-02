#' Data Smoothing using Savitzky-Golay filters


#' @param dir_in           The input directory.
#' @param dir_out          The output directory.
#' @param samples          The set of samples to be processed.
#' @param by_rows          Logical. Direction to apply the function. If TRUE it by rows (drift time direction).
#'                         If FALSE, applied by columns (that is the retention time direction).
#' @param filter_length    Numerical. Length of the filter.
#' @param polynomial_order Numerical. Order of the polynomial.
#' @return A filtered  gcims dataset.
#' @family Smoothing functions
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


gcims_smoothing <- function (dir_in, dir_out, samples, by_rows,
                             filter_length, polynomial_order){

  print(" ")
  print("  /////////////////////////")
  print(" /   Filtering the data /")
  print("/////////////////////////")
  print(" ")

  setwd(dir_in)
  m = 0
  for (i in samples){
    m = m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".mat")
    aux <- readMat(aux_string)[[1]]
    if (by_rows == TRUE){
      dimension <- 1
    } else if (by_rows == FALSE){
      dimension <- 2
    }
    M <- t(apply(aux, dimension, function(x) sgolayfilt(x, p = polynomial_order, n = filter_length)))
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }
}
