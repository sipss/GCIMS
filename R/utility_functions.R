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

#' Interpolates data according to retention and drift time sampling frequencies


#' @param dir_in          The input directory.
#' @param dir_out         The output directory.
#' @param samples         The set of samples to be processed.
#' @param by_rows         Logical. Direction to apply the function. If TRUE it
#'                        is applied by rows (drift time direction).
#'                        If FALSE, applied by columns
#'                        (that is the retention time direction).
#'
#' @return An interpolated dataset
#' @family Utility functions
#' @export
#' @importFrom signal interp1
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }

gcims_interpolate <- function(dir_in, dir_out, samples, by_rows){
  print(" ")
  print("  ////////////////////////")
  print(" /   Interpolating data /")
  print("////////////////////////")
  print(" ")

  setwd(dir_in)
  m = 0
  for (i in samples){
    m = m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(aux_string) #new
    aux <- as.matrix(aux_list$data$data_df)

    if (by_rows == TRUE){
      aux <- t(aux)
      x <- aux_list$data$drift_time
      step_x <- (x[length(x)]- x[1]) / (length(x) - 1)
      xi <- seq(from = x[1],
                by = step_x,
                length.out = length(x))
    } else if (by_rows == FALSE){
      x <- aux_list$data$retention_time
      step_x <- (x[length(x)]- x[1]) / (length(x) - 1)
      xi <- seq(from = x[1],
                by = step_x,
                length.out = length(x))
    }

    n <- dim(aux)[1]

    for (j in (1:n)){
      aux[j, ] <- interp1(x, aux[j, ], xi, method = "linear")
    }

    if (by_rows == TRUE){
      aux <- t(aux)
      aux_list$data$drift_time <- xi
    } else if (by_rows == FALSE){
      aux_list$data$retention_time <- xi
    }

    aux_list$data$data_df <- aux
    M <- aux_list
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }


}
