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
  print("  / Creation of the Feature Vector /")
  print("  //////////////////////////////////")
  print(" ")

  setwd(dir_in)
  m = 0
  for (i in samples){
    m <- m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(aux_string)
    aux <- aux_list$data$data_df
    aux <- as.vector(aux)
    aux_list$data$data_df <- aux
    M <-aux_list
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }
}

#' Interpolates data according to retention and drift time sampling frequencies


#' @param dir_in          The input directory.
#' @param dir_out         The output directory.
#' @param samples         The set of samples to be processed.
#' @param time            It indicates if the correction is going to be in the
#'                        drift time or in the retention time. It should be
#'                        introduce "Retention" for correcting the retention
#'                        time; or "Drift" for the drift time.
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

gcims_interpolate <- function(dir_in, dir_out, samples, time){
  print(" ")
  print("  ////////////////////////")
  print(" /   Interpolating data /")
  print("////////////////////////")
  print(" ")

  setwd(dir_in)
  m = -1
  for (i in c(0, samples)){
    m = m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(aux_string) #new
    aux <- as.matrix(aux_list$data$data_df)

    if (time == "Retention"){
      x <- aux_list$data$retention_time
      step_x <- (x[length(x)]- x[1]) / (length(x) - 1)
      xi <- seq(from = x[1],
                by = step_x,
                length.out = length(x))
    } else if (time == "Drift"){
      aux <- t(aux)
      x <- aux_list$data$drift_time
      step_x <- (x[length(x)]- x[1]) / (length(x) - 1)
      xi <- seq(from = x[1],
                by = step_x,
                length.out = length(x))
    }

    n <- dim(aux)[1]

    for (j in (1:n)){
      aux[j, ] <- interp1(x, aux[j, ], xi, method = "linear", extrap = TRUE)
    }

    if (time == "Retention"){
      aux_list$data$retention_time <- xi
    } else if (time == "Drift"){
      aux <- t(aux)
      aux_list$data$drift_time <- xi
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
  m = -1
  for (i in c(0, samples)){
    m = m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(aux_string) #new
    aux <- (as.matrix(aux_list$data$data_df))

    # Compute the total ion spectra
    aux_2 <- rowSums(aux)
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
    for (j in (1:dim(aux)[2])){
      aux[rip_bounds[1]: rip_bounds[2], j] <- seq(from = aux[rip_bounds[1], j],
                                                  to = aux[rip_bounds[2], j],
                                                  length.out = length(rip_bounds[1]: rip_bounds[2]))
    }

    aux_list$data$data_df <- (aux)
    M <- aux_list
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }
}

#' Reshape samples
#'
#' @param dir_in          The input directory.
#' @param dir_out         The output directory.
#' @param samples         The set of samples to be processed.
#' @return An cut matrix.
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
#'
#'

reshape_samples <- function(dir_in, dir_out, samples) {

  print(" ")
  print("  ///////////////////////////////")
  print(" /    Sample Matrix Reshape    /")
  print("///////////////////////////////")
  print(" ")
  m <- 0
  dimensions <- list(NULL)
  for (i in samples){
    m <- m + 1
    setwd(dir_in)
    aux_string <- paste0("M", samples[m], ".rds")
    aux_list <- readRDS(aux_string) #new
    aux <- (as.matrix(aux_list$data$data_df)) #new
    dimensions[[m]] <- dim(aux)
  }

  rts <- NULL
  dts <- NULL
  m <-0
  for (i in samples){
    m <- m + 1
    dts <- c(dts, dimensions[[m]][1])
    rts <- c(rts, dimensions[[m]][2])
  }
  rts <- min(rts)
  dts <- min(dts)

  m <- 0
  for (i in samples){
    m <- m + 1
    print(paste0("Sample ", i, " of ", length(samples)))
    aux_string <- paste0("M", samples[m], ".rds")
    aux_list <- readRDS(aux_string) #new
    aux <- (as.matrix(aux_list$data$data_df))
    aux <- aux[1:dts, 1:rts]
    aux_list$data$data_df <- aux
    aux_list$data$retention_time <- aux_list$data$retention_time[1:rts]
    aux_list$data$drift_time <- aux_list$data$drift_time[1:dts]
    M <- aux_list
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }
}

#' Data Dimentions Reduction to the Decimation

#' @param dir_in           The input directory.
#' @param dir_out          The output directory.
#' @param samples          The set of samples to be processed..
#' @return A reduced  gcims dataset.
#' @family Utility functions
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
    dtreduced <- apply(t(aux), 1, function(x) resample(x, p = 6, q = 32))
    resampleddt <- resample(aux_list$data$drift_time, p = 6, q = 32)
    rtreduced <- apply(t(dtreduced), 2, function(x) resample(x, p = 6, q = 32))
    resampledrt <- resample(aux_list$data$retention_time, p = 6, q = 32)
    aux_list$data$data_df <- t(rtreduced)
    aux_list$data$drift_time <- resampleddt
    aux_list$data$retention_time <- resampledrt
    M <- aux_list
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }
}
