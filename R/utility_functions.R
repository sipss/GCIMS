#' Converts Sample Matrices to Feature Vectors


#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param dir_out         Output directory. Where unfold data files are
#'   stored.
#' @param samples         Numeric vector of integers. Identifies the set of
#'   samples to which their data matrices need to be unfolded.
#' @return A set of S3 objects
#' @details \code{gcims_unfold} reshapes \emph{n} x \emph{m} GCIMS data into
#'   \emph{1} x (\emph{n} x \emph{m}) vectors. This is done for all samples in
#'   \code{samples}. \code{gcims_unfold} allows the latter use of pattern
#'   recognition techniques (e.g. Principal Component Analysis) that need to be
#'   fed with tabular data.
#' @note It is recommended to perform a previous reduction of data
#'   dimensionality using the functions \code{gcims_decimate} and
#'   \code{gcims_cut_samples} before applying \code{gcims_unfold} to the data.
#' @family Utility functions
#' @export
#' @examples
#' current_dir <- getwd()
#' dir_in <- system.file("extdata", "to_interpolate", package = "GCIMS")
#' dir_out <- tempdir()
#' samples <- 3
#'
#' # Example of Sample Data unfolding
#' # Before:
#' setwd(dir_in)
#' M3 <- readRDS("M3.rds")
#' data <- M3$data$data_df
#' data_dimensions <- dim(data)
#'
#' message("Before unfolding data has ", data_dimensions[1],
#' " rows and ", data_dimensions[2], " columns.")
#'
#' # After:
#' gcims_unfold(dir_in, dir_out, samples)
#' setwd(dir_out)
#' M3 <- readRDS("M3.rds")
#' data <- M3$data$data_df
#' data_length <- length(data)
#'
#' message("After unfolding data is a vector of length ", data_length,".
#' This value corresponds to product of the number
#' of columns by the number of rows of the original data.")
#'
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files))
#' setwd(current_dir)
#'
#'
gcims_unfold <- function(dir_in, dir_out, samples){
  print(" ")
  print("  //////////////////////////////////")
  print("  / Creation of the Feature Vector /")
  print("  //////////////////////////////////")
  print(" ")

  setwd(dir_in)
  m = -1
  for (i in c(0, samples)){
    m = m + 1
    if (m != 0){
      print(paste0("Sample ", m, " of ", length(samples)))
    }
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(aux_string)
    aux <- aux_list$data$data_df
    aux <- c(t(as.matrix(aux)))
    aux_list$data$data_df <- aux
    M <-aux_list
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }
}

#' Interpolates data to obtain constant sampling frequencies in retention and drift times

#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param dir_out         Output directory. Where interpolated data files are
#'   stored.
#' @param samples         Numeric vector of integers. Identifies the set of
#'   samples to which their data matrices need to be interpolated.
#' @param time            Sets the dimension to be corrected: drift time or
#'   retention time. Introduce "Retention" for retention time; or "Drift" for
#'   drift time.
#' @details \code{gcims_interpolate} performs a linear interpolation on gcims
#' data along the time axis selected in \code{gcims_interpolate}, and for all
#' samples in \code{samples}.
#' @return A set of S3 objects.
#' @references { Oppenheim, Alan V.; Schafer, Ronald W.; Buck, John R. (1999).
#'   "4". Discrete-Time Signal Processing (2nd ed.). Upper Saddle River, N.J.:
#'   Prentice Hall. p. 168. ISBN 0-13-754920-2. }
#' @family Utility functions
#' @export
#' @importFrom signal interp1
#' @examples
#' current_dir <- getwd()
#' dir_in <- system.file("extdata", "to_interpolate", package = "GCIMS")
#' dir_out <- tempdir()
#' samples <- 3
#'
#' # Example of Drift time Interpolation
#' # Before:
#' gcims_view_sample(dir_in, sample_num = samples, rt_range = NULL, dt_range = NULL)
#'
#' # After:
#' time <- "Drift"
#' gcims_interpolate(dir_in, dir_out, samples, time)
#' gcims_view_sample(dir_out, sample_num = samples, rt_range = NULL, dt_range = NULL)
#'
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files))
#' setwd(current_dir)
#'
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
    if (m != 0){
      print(paste0("Sample ", m, " of ", length(samples)))
    }
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

    aux_list$data$data_df <- round(aux)
    M <- aux_list
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }
}


#' Removes the Reactant Ion Peak (RIP) from samples

#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param dir_out         Output directory. Where RIP removed data files are
#'   stored.
#' @param samples         Numeric vector of integers. Identifies the set of
#'   samples to which their RIP has to be removed.
#' @details \code{gcims_remove_rip} substitutes the RIP by its corresponding
#'   linear approximation to the RIP baseline, for every spectrum in a sample.
#'   This process is repeated for all samples in \code{samples}. Use this
#'   function if you are interested in enhancing the contrast of peaks of sample
#'   images / chromatograms / spectra to be obtained from
#'   \code{gcims_view_sample} / \code{gcims_plot_chrom} /
#'   \code{gcims_plot_spec}.
#' @return A Set of S3 objects.
#' @family Utility functions
#' @export
#' @importFrom pracma findpeaks
#' @examples
#' current_dir <- getwd()
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' samples <- 3
#'
#' # Example of Reactant Ion Peak removal
#' # Before:
#' gcims_view_sample(dir_in, sample_num = samples, rt_range = NULL, dt_range = NULL)
#'
#' # After:
#' gcims_remove_rip(dir_in, dir_out, samples)
#' gcims_view_sample(dir_out, sample_num = samples, rt_range = NULL, dt_range = NULL)
#'
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files))
#' setwd(current_dir)
#'
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
    if (m != 0){
      print(paste0("Sample ", m, " of ", length(samples)))
    }
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
      if(closest_valley_ind > 1){
        rip_bounds <- valleys_pos[c(closest_valley_ind - 1,closest_valley_ind)]
      } else {
        rip_bounds <- valleys_pos[c(1, closest_valley_ind)]
      }

    }


    # Erase the peak (with style...)
    for (j in (1:dim(aux)[2])){
      aux[rip_bounds[1]: rip_bounds[2], j] <- seq(from = aux[rip_bounds[1], j],
                                                  to = aux[rip_bounds[2], j],
                                                  length.out = length(rip_bounds[1]: rip_bounds[2]))
    }

    aux_list$data$data_df <- round(aux)
    M <- aux_list
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }
}

#' Reshape samples
#'
#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param dir_out         Output directory. Where reshaped samples are stored.
#' @param samples         Numeric vector of integers. Identifies the set of
#'   samples to be reshaped.
#' @details \code{gcims_reshape_samples} ensures that all samples in a dataset
#'   have the same dimensions (number of data points) in retention and drift
#'   times. \code{gcims_reshape_samples} checks what are the minimum retention /
#'   drift time ranges a cuts all samples according to these ranges.
#'
#' @return A set of S3 objects.
#' @family Utility functions
#' @export
#' @examples
#' current_dir <- getwd()
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' samples <- c(3, 7, 8, 14, 20, 21, 22)
#'
#' # Example of reshaping samples
#' # (all samples must have the same
#' # retention and drift time dimensions).
#'
#' # Before:
#' nrow_before <- matrix(0, nrow = length(samples), ncol = 1)
#' ncol_before <- matrix(0, nrow = length(samples), ncol = 1)
#' setwd(dir_in)
#' for (i in seq_along(samples)){
#' aux_string <- paste0("M", samples[i], ".rds")
#' aux_list <- readRDS(aux_string)
#' nrow_before[i, 1] <- nrow(as.matrix(aux_list$data$data_df))
#' ncol_before[i, 1] <- ncol(as.matrix(aux_list$data$data_df))
#' }
#'
#' # After:
#' gcims_reshape_samples(dir_in, dir_out, samples)
#' setwd(dir_out)
#' nrow_after <- matrix(0, nrow = length(samples), ncol = 1)
#' ncol_after <- matrix(0, nrow = length(samples), ncol = 1)
#' for (i in seq_along(samples)){
#' aux_string <- paste0("M", samples[i], ".rds")
#' aux_list <- readRDS(aux_string)
#' nrow_after[i, 1] <- nrow(as.matrix(aux_list$data$data_df))
#' ncol_after[i, 1] <- ncol(as.matrix(aux_list$data$data_df))
#' }
#'
#' reshaping_info <- cbind(nrow_before, ncol_before, nrow_after, ncol_after)
#' colnames(reshaping_info) <- c("rows_before", "columns_before", "rows_after", "columns_after")
#' rownames(reshaping_info) <- c("M3", "M7", "M8", "M14", "M20", "M21", "M22")
#'
#' print(reshaping_info)
#'
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files))
#' setwd(current_dir)
gcims_reshape_samples <- function(dir_in, dir_out, samples) {

  print(" ")
  print("  ///////////////////////////////")
  print(" /    Sample Matrix Reshape    /")
  print("///////////////////////////////")
  print(" ")

  dimensions <- list(NULL)
  setwd(dir_in)
  m = 0
  for (i in samples){
    m = m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
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

#' Decimation

#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param dir_out         Output directory. Where decimated data files are
#'   stored.
#' @param samples         Numeric vector of integers. Identifies the set of
#'   samples to be decimated.
#' @param q_rt            Numeric. Binning factor in retention time.
#' @param q_dt            Numeric. Binning factor in drift time.
#' @details \code{gcims_decimate} performs downsampling in retention and drift
#'   time axes of GCIMS data. In particular, it reduces sampling frequency in
#'   retention and dritf times respectively, by the  factors \code{q_rt} and
#'   \code{q_dt}. Use this function if you are interested in both increasing the
#'   signal to noise ratio of data and compress it. Please take in to account
#'   that decimation also reduces data resolution.
#' @note \code{gcims_decimate} introduces a delay in retention and drift time
#'   axes.
#' @return A set of S3 objets.
#' @family Utility functions
#' @export
#' @references { Oppenheim, Alan V.; Schafer, Ronald W.; Buck, John R. (1999).
#'   "4". Discrete-Time Signal Processing (2nd ed.). Upper Saddle River, N.J.:
#'   Prentice Hall. p. 168. ISBN 0-13-754920-2. }
#' @importFrom signal resample
#' @examples
#' current_dir <- getwd()
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#'
#' # Before:
#' setwd(dir_in)
#' samples <- c(3, 7)
#' gcims_plot_chrom(dir_in, samples, dt_value = NULL,  rt_range = NULL, colorby = "Class")
#' gcims_plot_spec(dir_in, samples, rt_value = NULL,  dt_range = NULL, colorby = "Class")
#'
#' # After:
#' q_rt <- 5
#' q_dt <- 2
#' gcims_decimate(dir_in, dir_out, samples, q_rt, q_dt)
#' setwd(dir_out)
#' gcims_plot_chrom(dir_out, samples, dt_value = NULL,  rt_range = c(50, 226) , colorby = "Class")
#' gcims_plot_spec(dir_out, samples, rt_value = NULL,  dt_range = c(7.75, 9.6), colorby = "Class")
#'
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files))
#' setwd(current_dir)
#'

gcims_decimate <- function(dir_in, dir_out, samples, q_rt, q_dt){
  print(" ")
  print("  ///////////////////////")
  print("  / Samples Decimation /")
  print("  //////////////////////")
  print(" ")

  p_rt <- 1
  p_dt <- 1
  setwd(dir_in)
  m = -1
  for (i in c(0, samples)){
    m = m + 1
    if (m != 0){
      print(paste0("Sample ", m, " of ", length(samples)))
    }
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(aux_string)
    aux <- aux_list$data$data_df
    dtreduced <- apply(t(aux), 1, function(x) resample(x, p = p_dt, q = q_dt))
    resampleddt <- resample(aux_list$data$drift_time, p = p_dt, q = q_dt)
    rtreduced <- apply(t(dtreduced), 2, function(x) resample(x, p = p_rt, q = q_rt))
    resampledrt <- resample(aux_list$data$retention_time, p = p_rt, q = q_rt)
    aux_list$data$data_df <- t(rtreduced)
    aux_list$data$drift_time <- resampleddt
    aux_list$data$retention_time <- resampledrt
    M <- aux_list
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }
}


#' Cut samples in a retention time - drift time rectangle


#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param dir_out         Output directory. Where cut data files are
#'   stored.
#' @param samples         A vector. Set of samples to be aligned(e.g.: c(1, 2,
#'   3)).
#' @param rt_range        A  vector of two components. Beginning and end of the
#'   retention time cut.If NULL the complete retention time range is used.
#' @param dt_range        A  vector of two components. Beginning and end of the
#'   drift time cut. If NULL the complete drift time range is used.
#' @return A set of S3 objects.
#' @details \code{cut_samples} cuts a sample in a retention time - drift time
#'   rectangle according to the retention time / drift time ranges given by
#'   function arguments \code{rt_range} / \code{dt_range}. Use this function to
#'   focus on the retention time - drift time region where chemical information
#'   is more abundant, that is, where you can find a high peak densities by
#'   visual inspection.
#' @note By reducing the size of data, the computational time of signal
#'   pre-processing stage can reduced substantially.
#' @family Utility functions
#' @export
#' @examples
#' current_dir <- getwd()
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' samples <- 3
#'
#' # Example Sample data cutting:
#' # Before:
#' gcims_view_sample(dir_in, sample_num = samples, rt_range = NULL, dt_range = NULL)
#'
#' # After:
#' rt_range <-c(70, 125)
#' dt_range <- c(8, 9.25)
#' gcims_cut_samples(dir_in, dir_out, samples, rt_range, dt_range)
#' gcims_view_sample(dir_out, sample_num = samples, rt_range = NULL, dt_range = NULL)
#'
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files))
#' setwd(current_dir)
#
gcims_cut_samples <- function(dir_in, dir_out, samples, rt_range, dt_range){



  #---------------#
  #   FUNCTIONS   #
  #---------------#


  #-----------------#
  #   cut_samples   #
  #-----------------#


  cut_samples <- function(loaded_sample, rt_range = NULL, dt_range = NULL) {

    aux_list <- loaded_sample
    aux <- (as.matrix(aux_list$data$data_df))

    #SOME CHECKS
    retention_time <- aux_list$data$retention_time
    drift_time <- aux_list$data$drift_time
    cond_1_rt <- (rt_range[1] - retention_time[1]) < 0
    cond_2_rt <- (rt_range[2] - retention_time[length(retention_time)]) > 0
    cond_1_dt <-(dt_range[1] - drift_time[1]) < 0
    cond_2_dt <-(dt_range[2] - drift_time[length(drift_time)]) > 0


    if(is.null(rt_range)){# old
      rt_ind <- c(1, dim(aux)[2])

    } else{
      if(cond_1_rt | cond_2_rt){
        stop("Retention time range out of bounds.")
      }
      rt_ind  <- c(which.min(abs(retention_time - rt_range[1])), which.min(abs(retention_time - rt_range[2])))
      if( rt_ind[1] == rt_ind[2]){
        stop("Initial and Final retention time values can't be equal in the variable rt_range.")
      }
    }



    if(is.null(dt_range)){# old
      dt_ind <- c(1, dim(aux)[1])
    } else{
      if(cond_1_dt | cond_2_dt){
        stop("Drift time range out of bounds.")
      }
      dt_ind  <- c(which.min(abs(drift_time - dt_range[1])), which.min(abs(drift_time - dt_range[2])))
      if( dt_ind[1] == dt_ind[2]){
        stop("Initial and Final drift time values can't be equal in the variable dt_range.")
      }
    }

    sel_index_rt <- rt_ind[1]: rt_ind[2]
    sel_index_dt <- dt_ind[1]: dt_ind[2]

    if(is.null(rt_range)){

    } else if((class(sel_index_rt) == "integer") & (sel_index_rt[2] > sel_index_rt[1])){
    } else {
      stop("Possible errors: 1) The selected vector of indexes corresponding to the provided retention time range is not an integer vector, 2) or rt_range[2] <= rt_range[1])")
    }

    if(is.null(dt_range)){

    } else if((class(sel_index_dt) == "integer") & (sel_index_dt[2] > sel_index_dt[1])){
    } else {
      stop("Possible errors: 1) The selected vector of indexes corresponding to the provided drift time range is not an integer vector, 2) or dt_range[2] <= dt_range[1])")
    }

    aux_list$data$retention_time <- retention_time[sel_index_rt]
    aux_list$data$drift_time  <- drift_time[sel_index_dt]
    aux_list$data$data_df <- round(aux_list$data$data_df[sel_index_dt, sel_index_rt]) #transposed when matrix!

    return(aux_list)
  }



  print(" ")
  print("  /////////////////////")
  print(" /    Cut Samples    /")
  print("/////////////////////")
  print(" ")


  #-------------#
  #     MAIN    #
  #-------------#

  setwd(dir_in)
  m = -1;
  for (i in c(1, samples)){
    m = m + 1
    if (m != 0){
      print(paste0("Sample ", m, " of ", length(samples)))
    }
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(aux_string) #new
    M <- cut_samples(aux_list, rt_range, dt_range)
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)

  }

}




#' Shifting

#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param dir_out         Output directory. Where decimated data files are
#'   stored.
#' @param samples         Numeric vector of integers. Identifies the set of
#'   samples to be decimated.
#' @details \code{gcims_shifting} performs peak shifting in retention time axes
#'   of GCIMS data. In particular, it shifts all the peaks along the retention
#'   time axis in order to have the RIP peaks at the exact same time.
#' @note \code{gcims_shifting} reduce the retention time axis of some of
#'   the samples.
#' @return A set of S3 objets.
#' @family Utility functions
#' @export
#' @references { Oppenheim, Alan V.; Schafer, Ronald W.; Buck, John R. (1999).
#'   "4". Discrete-Time Signal Processing (2nd ed.). Upper Saddle River, N.J.:
#'   Prentice Hall. p. 168. ISBN 0-13-754920-2. }
#' @examples
#' current_dir <- getwd()
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#'
#' # Before:
#' setwd(dir_in)
#' samples <- c(3, 7)
#' gcims_plot_chrom(dir_in, samples, dt_value = NULL,  rt_range = NULL, colorby = "Class")
#' gcims_plot_spec(dir_in, samples, rt_value = NULL,  dt_range = NULL, colorby = "Class")
#'
#' # After:
#' gcims_shifting(dir_in, dir_out, samples)
#' setwd(dir_out)
#' gcims_plot_chrom(dir_out, samples, dt_value = NULL,  rt_range = c(50, 226) , colorby = "Class")
#' gcims_plot_spec(dir_out, samples, rt_value = NULL,  dt_range = c(7.75, 9.6), colorby = "Class")
#'
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files))
#' setwd(current_dir)
#'

gcims_shifting <- function(dir_in, dir_out, samples){


  print(" ")
  print("  ////////////////////////")
  print(" /      Peaks Sifting   /")
  print("////////////////////////")
  print(" ")


  setwd(dir_in)

  m <- -1
  tics <- NULL
  for (i in c(1,samples)){
    m <- m + 1
    if (m != 0){
      print(paste0("Sample ", m, " of ", length(samples)))
    }

    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(aux_string) #new
    aux <- as.matrix(aux_list$data$data_df)
    tic <- colSums(aux)
    tics <- rbind(tics, tic)
  }

  mins <- apply(tics, 1, which.min)
  # plot(x = c(1:64), y = mins)
  # intensities <- apply(tics, 1, min)
  # plot(x = c(1:64), y = intensities)
  referenceindex <- min(mins)
  referencetime <- which.max(mins)

  m <- 0
  for (i in c(1,samples)){
    m <- m + 1
    if (m != 0){
      print(paste0("Sample ", m, " of ", length(samples)))
    }

    reference <- readRDS(paste0("M", referencetime, ".rds"))
    reference <- reference$data$data_df
    referencetic <- colSums(reference)
    reference <- reference[,-c(1:abs(referenceindex - which.min(referencetic)))]
    dimension <- dim(reference)[2]
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(aux_string) #new
    aux <- as.matrix(aux_list$data$data_df)
    auxtic <- colSums(aux)
    if (abs(referenceindex - which.min(auxtic)) > 0){
      aux <- aux[,-c(1:abs(referenceindex - which.min(auxtic)))]
    } else {
      aux <- aux
    }
    aux_list$data$retention_time <- aux_list$data$retention_time[1:dimension]
    aux_list$data$data_df <- aux[,c(1:dimension)]
    M <- aux_list
    setwd(dir_out)
    saveRDS(M, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }
}




