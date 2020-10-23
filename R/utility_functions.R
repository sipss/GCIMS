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


#' Removes RIP from all samples


#' @param dir_in          The input directory.
#' @param dir_out         The output directory.
#' @param samples         The set of samples to be processed.
#' @return A RIP removed dataset
#' @family Utility functions
#' @export
#' @importFrom pracma findpeaks
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#'}
gcims_remove_rip <- function(dir_in, dir_out, samples){
  print(" ")
  print("  ////////////////////////")
  print(" /    Removing the RIP  /")
  print("////////////////////////")
  print(" ")

  setwd(dir_in)
  m = 0
  for (i in samples){
    m = m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(aux_string) #new
    aux <- t(as.matrix(aux_list$data$data_df))

    # Compute the total ion spectra
    aux_2 <- colSums(aux)
    peaks_info <- findpeaks(aux_2)

    # Look for rip position
    rip_pos_ind <- which.max(peaks_info[ , 1])
    rip_pos <- peaks_info[rip_pos_ind, 2]

    # Look for the beginning and ending of the RIP (searching the closest minima to it)
    valleys_info <- findpeaks(-aux_2)
    valleys_pos <- valleys_info[ , 2]
    closest_valley_ind <- which.min(abs(valleys_pos - rip_pos))

    # Select the RIP region
    if(valleys_pos[closest_valley_ind] < rip_pos){
      rip_bounds <- valleys_pos[c(closest_valley_ind,closest_valley_ind + 1)]
    } else if (valleys_pos[closest_valley_ind] > rip_pos){
      rip_bounds <- valleys_pos[c(closest_valley_ind - 1,closest_valley_ind)]
    }

    # Erase the peak (with style...)
    for (j in (1:dim(aux)[1])){
      aux[j, rip_bounds[1]: rip_bounds[2]] <- seq(from = aux[j, rip_bounds[1]],
                                                  to = aux[j, rip_bounds[2]],
                                                  length.out = length(rip_bounds[1]: rip_bounds[2]))
    }

    aux_list$data$data_df <- t(aux)
    M <- aux_list
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }



}
