#' Create a modified Gaussian (non symmetric) from two Gaussians

#' @param x              A vector of indexes.
#' @param mu             Integer.Position of the modified Gaussian maximum.
#'                       It corresponds to the mean of the left and right Gaussian
#'                       that generate the modified Gaussian
#' @param sigma1         Numeric. Standard deviation of the left Gaussian.
#' @param sigma2         Numeric. Standard deviation of the right Gaussian.
#' @return A vector of intensities of the modified Gaussian.
#' @export
#' @family Dummy data functions
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }


modgauss <- function(x, mu, sigma1, sigma2){

  gauss <- function(x, mu, sigma){
    gaussian <- (1 / (sqrt(2 * pi) * sigma)) * exp(-(1/2) * ((x - mu) / sigma) ^ 2)
  }

  gaussian_1 <- gauss(x, mu, sigma1)
  gaussian_2 <- gauss(x, mu, sigma2)

  max_pos <- which.max(gaussian_1)
  modgaussian <- as.vector(matrix(0, 1, length(x)))
  modgaussian[1:max_pos] <- gaussian_1[1:max_pos]
  modgaussian[(max_pos + 1): length(x)] <- gaussian_2[(max_pos + 1): length(x)]

  return(modgaussian)
}


#' Generate a random shift in peak positions (in retention or drift axes)

#' @param min_shift      Integer. Minimum shift in indexes allowed.
#' @param max_shift      Integer. Maximum shift in indexes allowed.
#' @param npeaks         Integer. Number of peaks to be shifted.
#' @return A vector of shifts (one per peak to be shifted).
#' @importFrom pracma randi
#' @family Dummy data functions
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }

shift <- function(min_shift, max_shift, npeaks){
  values <- randi(c(min_shift,  max_shift), n = 1, m = npeaks)
}



#' Generate a baseline matrix

#' @param max_a          Numeric. Maximum allowed baseline level.
#' @param nrow           Integer. Number of row of the output matrix.
#' @param ncol           Integer. Number of columns of the output matrix.
#' @return A matrix with a constant value (the baseline)
#' @importFrom pracma rand
#' @family Dummy data functions
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }

add_baseline <- function(max_a, nrow, ncol){


  a <- as.vector(rand(1,1) * max_a)
  baseline_mat <-matrix(a, nrow, ncol)
  return(baseline_mat)

}


#' Simulate peak tailing
#'
#' @param max_a          Numeric. Maximum allowed offset.
#' @param max_b          Numeric. Maximum allowed multiplicative constant.
#' @param max_c          Numeric. Maximum allowed exponent.
#' @param nrow           Integer. The output vector length.
#' @return A vector that simulates peak tailing as an offset plus an exponential decay
#' @importFrom pracma rand
#' @family Dummy data functions
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }

add_tail <- function(max_a, max_b, max_c, nrow){
  a <- as.vector(pracma::rand(1,1) * max_a)
  b <- as.vector(pracma::rand(1,1) * max_b)
  c <- as.vector(pracma::rand(1,1) * max_c)
  x <- 1:nrow
  tail <- a + b * exp (- c * x)
  return(tail)

}


#' Generate profiles of concentration including non-idealities
#'
#' @param C              Matrix. A matrix of concentrations to be filled
#' @param column         Column matrix to be filled with its concentration profile.
#' @param max_pos        Integer. The postion for the maximum of concentration
#'                       in the profile (in indexes)
#' @param max_value      Numeric. Maximum concentration of the profile.
#' @param rise_time      Integer. Number of indexes needed to reach the
#'                       maximum value of concentration from the baseline.
#' @param fall_time      Integer. Number of indexes needed to reach the baseline
#'                       from the maximum value of concentration.
#' @param rt_shift       Vector of integers. The shift in retention time indexes for
#'                       each of the peaks.
#' @param nscans         Integer. Length of the concentration profile (in indexes). It
#'                       its equal to the number of spectra or scans.
#' @return A matrix of Concentration profiles.
#' @family Dummy data functions
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }
fill_concentration <- function(C,column, max_pos, max_value, rise_time, fall_time, rt_shift, nscans){

  peak_rt_left_pos <- (max_pos - (rise_time - 1) + rt_shift): (max_pos + rt_shift) #(max_pos - (rise_time - 1) + rt_shift[1, 1]): (max_pos + rt_shift[1, 1])
  peak_rt_right_pos <- (max_pos + 1 + rt_shift): (max_pos + fall_time + rt_shift) #(max_pos + 1 + rt_shift[1, 1]): (max_pos + fall_time + rt_shift[1, 1])
  peak_rt_left_value <- seq(from = 0, to = max_value, length.out = rise_time)
  peak_rt_right_value <- seq(from = max_value, to = 0, length.out = fall_time)

  if (peak_rt_left_pos[1] > nscans) {
    #don't include peak
  } else if ((peak_rt_left_pos[1] < nscans) & ((peak_rt_left_pos[1] + rise_time) > nscans)) {
    peak1_rt_left_pos <- (max_pos - (rise_time - 1) + rt_shift): nscans #    peak1_rt_left_pos <- (max_pos - (rise_time - 1) + rt_shift[1, 1]): nscans
    C[peak_rt_left_pos, column] <- peak_rt_left_value[1:length(peak_rt_left_pos)]
  } else if ((peak_rt_right_pos[1] < nscans) & ((peak_rt_right_pos[1] + fall_time) > nscans)) {
    peak_rt_right_pos <- (max_pos + 1 + rt_shift): nscans #    peak_rt_right_pos <- (max_pos + 1 + rt_shift[1, 1]): nscans
    C[peak_rt_left_pos, column] <- peak_rt_left_value
    C[peak_rt_right_pos, column] <- peak_rt_right_value[1:length(peak_rt_right_pos)]
  } else if((peak_rt_right_pos[length(peak_rt_right_pos)]) < 1){
    #don't include peak
  } else if((peak_rt_right_pos[length(peak_rt_right_pos)] > 1) & ((peak_rt_right_pos[length(peak_rt_right_pos)] - fall_time) < 1)){
    peak_rt_right_pos <-  1: (max_pos + fall_time + rt_shift) # peak_rt_right_pos <-  1: (max_pos + fall_time + rt_shift[1, 1])
    C[peak_rt_right_pos, column] <- peak_rt_right_value[(length(peak_rt_right_value) - length(peak_rt_right_pos) + 1):length(peak_rt_right_value)]
  } else if ((peak_rt_left_pos[1] < 1) & ((peak_rt_left_pos[1] + rise_time) > 1)) {
    peak1_rt_left_pos <- 1: (max_pos + rt_shift) #peak1_rt_left_pos <- 1: (max_pos + rt_shift[1, 1])
    C[peak_rt_left_pos, column] <- peak_rt_left_value[(length(peak_rt_left_value) - length(peak_rt_left_pos) + 1):length(peak_rt_left_value)]
    C[peak_rt_right_pos, column] <-   peak_rt_right_value
  } else{
    C[peak_rt_left_pos, column] <- peak_rt_left_value
    C[peak_rt_right_pos, column] <- peak_rt_right_value
  }
  return(C)
}


#' Forces that sampling frequencies in retention or drift times are not constant
#'
#' @param sampling_period   Numeric. The actual sampling period.
#' @param constant_sampling Logical. TRUE: no modification is needed. FALSE: an
#'                          error in sampling is computed.
#'                          in the profile (in indexes)
#' @param axis_length       Numeric. Maximum concentration of the profile.
#' @return A vector of errors in sampling to be added to the actual retention or
#' drift time vectors. If constant sampling is TRUE, the error vector is filled with zeros.
#' @importFrom pracma rand randn
#' @family Dummy data functions
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }

force_interpolation <- function (sampling_period, constant_sampling, axis_length){

  decimal_num_count <- function(x){
    stopifnot(class(x)=="character")
    x <-gsub("(.*)(\\.)|([0]*$)","",x)
    nchar(x)
  }

  decimals <- decimal_num_count(format(sampling_period))
  sampling_error <- round(rand(n = 1, m = axis_length) * (10 ^ (-(decimals))), decimals)
  sign_sampling_error <- randn(n = 1, m = axis_length)
  sign_sampling_error <- sign_sampling_error / abs(sign_sampling_error)
  sampling_error <- sampling_error * sign_sampling_error

  if (constant_sampling){
    sampling_error <- 0 * sampling_error
  }

  sampling_error
}

#' Create sample params and labels for two different sample classes
#'
#' @return A list of with parameters for class1 (list), parameters for class (list)
#' and their corresponding labels.
#' @export
#' @family Dummy data functions
#' @examples
#' params <- gcims_create_dummy_params()
#'
gcims_create_dummy_params <- function(){

  #EXAMPLE OF PARAMETER CHOICE:

  # matrix parameters
  matrix_params <- list(npoints = 200, nscans = 50)

  # peak position parameters
  peak_positions <- list(ret_time = list(pr2 = 13, pr3 = 35),
                          drift_time = list(pd1 = 50, pd2 = 90, pd3 = 130))

  # peak height parameters
  peak_heights <- list(hp2 = 1, hp3 = 0.7)

  # peak asymmetry parameters
  peak_asym <- list(ret_time = list(rise_time = 7, fall_time = 10), drift_time = list(sigma_left = 2, sigma_right = 2.2))

  # peak asymmetry parameters
  axes_limits <- list(ret_time = list(sampling_period = 0.21, offset = 3.36),
                      drift_time = list(sampling_period = 0.0065, offset = 7.6557))

  non_idealities <- list(max_noise = 0.002, # noise parameter
                         baseline = 0.02,   # baseline parameter
                         ret_time_tail = list(max_a = 0.03, max_b = 0.05, max_c = 0.5), # tailing parameters
                         shift = list(ret_time = list(min_shift = -3, max_shift = 3),   # time shifting parameters
                                      drift_time =list(min_shift = -5, max_shift = 5)),
                         constant_sampling = list(ret_time = TRUE, drift_time = TRUE)  # Sampling shifting parameters
  )


  # CLASS 1 parameters:
  params_class1 <- list(matrix_params = matrix_params,
                        peak_positions = peak_positions,
                        peak_heights = peak_heights,
                        peak_asym = peak_asym,
                        axes_limits = axes_limits,
                        non_idealities = non_idealities
                        )

  # CLASS 2 parameters:
  # Reediting CLASS 1 parameters
  params_class2 <- params_class1

  # peak position parameters
  params_class2$peak_positions <- list(ret_time = list(pr2 = 30, pr3 = 14),
                                       drift_time = list(pd1 = 50, pd2 = 60, pd3 = 170))
  # peak height parameters
  params_class2$peak_heights <- list(hp2 = 0.5, hp3 = 0.8)

  labels <-c("Case", "Control")

  params <- list(params_class1 = params_class1, params_class2 = params_class2, labels = labels)
  return(params)

}




#' Create a dummy sample taking into account non-idealities in sampling process
#'
#' @param params         A list that include sample characteristics.
#' @return A dummy sample formally equal to the samples readed from an gc-ims
#' @importFrom pracma randn
#' @export
#' @family Dummy data functions
#' @examples
#' params <- gcims_create_dummy_params()
#' M1 <- gcims_create_dummy_sample(params$params_class1)

gcims_create_dummy_sample <- function(params){

  #Extract parameters from params:
  # matrix parameters
  npures <- 2
  npoints <- params$matrix_params$npoints
  nscans <- params$matrix_params$nscans

  # peak position parameters
  pr2 <- params$peak_positions$ret_time$pr2
  pr3 <- params$peak_positions$ret_time$pr3
  pd1 <- params$peak_positions$drift_time$pd1 #corresponds to rip
  pd2 <- params$peak_positions$drift_time$pd2
  pd3 <- params$peak_positions$drift_time$pd3

  # peak height parameters
  hp2 <- params$peak_heights$hp2
  hp3 <- params$peak_heights$hp3

  # peak asymmetry parameters
  rise_time <- params$peak_asym$ret_time$rise_time
  fall_time <- params$peak_asym$ret_time$fall_time
  sigma_left <- params$peak_asym$drift_time$sigma_left
  sigma_right <- params$peak_asym$drift_time$sigma_right

  # sampling parameters
  rt_sampling_period <- params$axes_limits$ret_time$sampling_period
  rt_offset <- params$axes_limits$ret_time$offset
  dt_sampling_period <- params$axes_limits$drift_time$sampling_period
  dt_offset <- params$axes_limits$drift_time$offset
  rt_constant_sampling <- params$non_idealities$constant_sampling$ret_time
  dt_constant_sampling <- params$non_idealities$constant_sampling$drift_time

  # noise parameter
  max_noise <- params$non_idealities$max_noise

  # Constant baseline parameter
  baseline <- params$non_idealities$baseline

  # Tailing parameters
  max_a <- params$non_idealities$ret_time_tail$max_a
  max_b <- params$non_idealities$ret_time_tail$max_b
  max_c <- params$non_idealities$ret_time_tail$max_c

  # Time shifting parameters
  rt_minshift <- params$non_idealities$shift$ret_time$min_shift
  rt_maxshift <- params$non_idealities$shift$ret_time$max_shift
  dt_minshift <- params$non_idealities$shift$drift_time$min_shift
  dt_maxshift <- params$non_idealities$shift$drift_time$max_shift

  #####################################################################################
  # COMPUTATIONS (It is assumed that only 2 pure substances plus the RIP are present) #
  #####################################################################################

  # CONCENTRATION MATRIX #

  # Initializing Concentration profiles
  Creal <- matrix(0, nrow = nscans, ncol = npures + 1)
  Sreal <- matrix(0, nrow = npures + 1, ncol = npoints)

  # Shift in retention time
  rt_shift <- shift(min_shift = rt_minshift,
                    max_shift = rt_maxshift,
                    npeaks = 2) #include rip


  # Ideal profiles of concentration (only rise and fall times of profiles can be different).
  # The 2 pure substances concentration profies are computed first

  # First substance (column 2 in Creal)
  Creal <- fill_concentration(C = Creal, column = 2,
                              max_pos = pr2, max_value = hp2,
                              rise_time = rise_time, fall_time = fall_time,
                              rt_shift = rt_shift[1, 1], nscans = nscans) # rt_shift = rt_shift, nscans = nscans)

  # Second substance (column 3 in Creal)
  Creal <- fill_concentration(C = Creal, column = 3,
                              max_pos = pr3, max_value = hp3,
                              rise_time = rise_time, fall_time = fall_time,
                              rt_shift = rt_shift[1, 2], nscans = nscans)

  # Then the RIP (column 1 in Creal).
  # Note that RIP saturation effect is included within the code
  Creal[, 1] <- 1 - Creal[, 2] - Creal[, 3]
  Creal[Creal[, 1] < 0, 1] <- 0 #saturation effect
  Creal[Creal[, 2] > 1, 2] <- 1 #saturation effect
  Creal[Creal[, 3] > 1, 3] <- 1 #saturation effect
  sat_sum_indexes <- (Creal[, 2] + Creal[, 3]) > 1
  Creal[sat_sum_indexes, 3] <- 1 - Creal[sat_sum_indexes, 2]

  # Include Tailing in concentration profiles
  C_tail <- add_tail(max_a =  max_a, max_b = max_b, max_c = max_c, nrow = nscans)
  Creal[, 1] <- Creal[, 1] + C_tail
  Creal[, 2] <- Creal[, 2] + C_tail
  Creal[, 3] <- Creal[, 3] + C_tail

  # SPECTRA MATRIX #

  # Initializing Spectra profiles
  Sreal <- matrix(0, nrow = npures + 1, ncol = npoints)

  # Shift in drift time
  dt_shift <- shift(min_shift = dt_minshift,
                    max_shift = dt_maxshift,
                    npeaks = 3)
  # Create spectra profiles including drift time shifting
  # Sreal[1, ] <- gauss(1:npoints, pd1 + dt_shift[1, 1], 2)
  # Sreal[2, ] <- gauss(1:npoints, pd2 + dt_shift[1, 2], 4)
  # Sreal[3, ] <- gauss(1:npoints, pd3 + dt_shift[1, 3], 5)

  Sreal[1, ] <- modgauss(1:npoints, pd1 + dt_shift[1, 1], sigma_left, sigma_right)
  Sreal[2, ] <- modgauss(1:npoints, pd2 + dt_shift[1, 2], sigma_left, sigma_right)
  Sreal[3, ] <- modgauss(1:npoints, pd3 + dt_shift[1, 3], sigma_left, sigma_right)

  # DATA MATRIX #

  # Compute data matrix as the product of C*S and plus E
  # (E being a matrix that includes baseline and noise effects)

  # Compute E
  E <- max_noise * randn(nscans, npoints) +
    add_baseline(max_a = baseline,
                 nrow = nscans,
                 ncol = npoints)

  # Then compute D:
  D <- Creal %*% Sreal + E

  # CREATE THE R LIST #

  #Initialize metadata and data fields
  metadata <- list(condition = NULL, experiment = NULL, group = NULL,
                   index = NULL, raw_path = NULL, replicate = NULL)
  data <- list(retention_time = NULL, drift_time = NULL, data_df = NULL)

  # Include the matrix D
  data$data_df = t(D)

  # Modify time axes to force the latter interpolation
  # Retention Time
  rt_sampling_error <- force_interpolation(rt_sampling_period, rt_constant_sampling, nscans)
  data$retention_time <- as.vector((rt_sampling_period * (1:nscans)) + rt_offset + rt_sampling_error)

  # Drift Time
  dt_sampling_error <- force_interpolation(dt_sampling_period, dt_constant_sampling, npoints)
  data$drift_time <- as.vector((dt_sampling_period * (1:npoints)) + dt_offset + dt_sampling_error)

  # Join data and metadata
  dd_list <- list(metadata = metadata, data = data)

  # Return the list
  return(dd_list)
}

 # M1 <-create_dummy_sample(params_class1)
 # image((M1$data$data_df))
 #  M2 <- create_dummy_sample(params_class2)
 #  image((M2$data$data_df))
 # # M3 <- create_dummy_sample(params_class3)
 # # image((M3$data$data_df))




#' Create a dummy dataset
#'
#' @param dir_in            The input directory.
#' @param dir_out           The output directory.
#' @param params            A list that include two list of
#'                          sample characteristics (one per sample class)
#' @param samples_per_class Integer. Number of sample repetitions per class
#' @return A dummy dataset to be stored in dir_out
#' @export
#' @family Dummy data functions
#' @examples
#' wd <- getwd()
#' dir_in <- tempdir()
#' dir_out <-  file.path(dir_in,"dummy")
#' list.files(path = dir_out, pattern = NULL, all.files = FALSE, full.names = FALSE)
#' dir.create (dir_out, FALSE)
#' samples_per_class <- 10
#' params <- gcims_create_dummy_params()
#' gcims_create_dummy_set(dir_in, dir_out, samples_per_class, params)
#' list.files(path = dir_out, pattern = NULL, all.files = FALSE, full.names = FALSE)
#' unlink(dir_out, recursive = TRUE)
#' setwd(wd)
#'

gcims_create_dummy_set <- function(dir_in, dir_out, samples_per_class, params){
  labels <- params$labels
  for (i in seq(from = 1, to = 2 * samples_per_class, by = 1)){
    if (i <= samples_per_class){
      m <- 1
    } else{
      m <- 2
    }
    #m <- m + 1
    print(paste0("Sample ", i, " of ", 2*(samples_per_class)))
    aux_list <- gcims_create_dummy_sample(params[[m]])
    aux_list$metadata$condition <- labels[m]
    setwd(dir_out)
    saveRDS(aux_list, file = paste0("M", i, ".rds"))
    setwd(dir_in)
  }
}



# GENERAR DOCUMENTACION
# VER COMO GENERAR EL DATASET EN UNA CARPETA DUMMY QUE NO SE CARGUE EN EL PAQUETE
# SUBIR A GITHUB (CHECKS Y DOCUMENTACION INCLUIDAS)
# CONTARSELO A CELIA PARA QUE LO PRUEBE
