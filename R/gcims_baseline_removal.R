#' Baseline Removal using  Using Psalsa algorithm


#' @param dir_in The input directory.
#' @param dir_out The output directory.
#' @param samples The set of samples to be processed.
#' @param by_rows Logical. Direction to apply the function. If TRUE it by rows (retention time direction). If FALSE, applied by columns (that is the drift time direction).
#' @param spectra  matrix with one spectrum per row
#' @param lambda smoothing parameter (generally 1e5 - 1e8)
#' @param p      asymmetry parameter
#' @param k      peak height parameter (usually 5\% of maximum intensity)
#' @return A baseline removed gcims dataset.
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




gcims_baseline_removal <- function(dir_in, dir_out, samples, by_rows,lambda, p, k){


  print(" ")
  print("  /////////////////////////")
  print(" /    Baseline removal   /")
  print("/////////////////////////")
  print(" ")

  setwd(dir_in)
  m = 0;
  for (i in (samples)){
    m = m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux <- readRDS(aux_string)
    M <- NULL
    setwd(wd)
    print(getwd())
    if (by_rows == TRUE){
      dimension <- 1
    } else if (by_rows == FALSE){
      dimension <- 2
    }
    M <- t(apply(aux, dimension, function(x) psalsa(x, lambda, p, k)))
    M <- aux - M
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }
}
