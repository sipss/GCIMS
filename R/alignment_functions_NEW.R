#' Sample alignment in drift and retention time axes
#'
#' This alignment function corrects misalignments in drift and retention
#' time axes. In drift time, the correction is multiplicative, meaning:
#' drift_time_corr = Kcorr * drift_time. The constant Kcorr is obtained
#' by correcting the position of the Reactant Ion Peak (RIP) in a sample
#' to make it coincide with the reference position of the RIP (the median position
#' of all RIPs in the dataset). Changes in pressure and temperature
#' during measurements may lead to correction factors between 0.9 and 1.1.
#  This alignment technique can be used when the pressure and temperature
#' measurements are not available. The correction in retention time performed
#' using parametric time warping. The warping function is optimized by searching
#' the polynomial degree (between 1 and 5) that maximizes the average
#' correlation between the reference RIC (Reactant Ion Chromatogram) and each of
#' the individual samples.
#'
#' @param dir_in             Input directory. Where input data files are loaded
#'  from.
#' @param dir_out            Output directory. Where aligned data files are stored.
#' @param samples            Numeric vector of integers. Identifies the set of
#'                           sample to be visualized from the dataset.
#' @param aligment_data      A list containing two matrices: tis and rics, the Total
#'                           Ion Spectra and Reactant Ion Chromatograns of data.
#' @importFrom signal interp1
#' @importFrom ptw ptw bestref
#' @importFrom stats cor median
#' @export
#' @return A set of S3 objects. Additionally, it returns a list containing the correction
#' factors in drift time (Kcorr_samples, a vector with as many components as samples), the reference
#' Reactant Ion Chromatogram (ric_ref), and correction type for retention time alignment (correction_type).
#' The output a variable correction_type is a whole number between 0 and 5. If 0, no correction is applied.
#' Between 1 and 5 it the optimal degree of the polynomial that corrects retention time axis.
#'
gcims_align_samples <- function(dir_in, dir_out, samples, aligment_data){

  # For drift time:
  # Reference position of rip in drift time (in indexes)
  tis <- aligment_data$tis
  rip_ref_idx <- round(median(apply(X = tis, MARGIN = 1, FUN = which.max), na.rm = TRUE))
  Kcorr_samples <- 0 * samples

  # For retention time:
  # Sample idx to correct for retention time using the ric to align
  rics <- aligment_data$rics
  ref_ric_sample_idx <- find_reference_ric(rics)
  # Search for the optimal polynomial degree of the warping function:
  correction_type <- gcims_optimize_polynomial(rics, ref_ric_sample_idx)
  # Select reference RIC
  ric_ref <- rics[ref_ric_sample_idx, ]

  m <- 0
  for (i in  samples){
    m <- m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(file.path(dir_in, aux_string))

    # First, we correct drift time:
    # Compute the multiplicative correction factor:
    Kcorr_samples[m] <- gcims_compute_Kcorr(aux_list, tis, rip_ref_idx)
    # Apply the correction to data:
    aux_list <- gcims_align_dt(aux_list, Kcorr_samples[m])

    # Then, we correct retention time
    # Apply the correction to data:
    aux_list <- gcims_align_rt(aux_list, ric_ref, correction_type)

    # Finally, we save data
    saveRDS(aux_list, file = file.path(dir_out, paste0("M", i, ".rds")))
  }
  aligment_info <- list(Kcorr_samples = Kcorr_samples, ric_ref = ric_ref, correction_type = correction_type)
}


# DRIFT TIME
gcims_compute_Kcorr <- function(aux_list, rip_ref_idx) {

  # Extract drift time axis
  drift_time <- aux_list$data$drift_time
  # Compute the position of reference RIP in drift time (ms)
  rip_ref_dt <- drift_time[rip_ref_idx]
  # Extract data matrix
  aux <- as.matrix(aux_list$data$data_df)
  # Compute sample Total Ion Spectrum
  tis <- rowSums(aux)
  # Compute the position of RIP in drift time (in indexes)
  rip_idx <- which.max(tis)
  # Compute the position of RIP in drift time (in ms)
  rip_dt <- drift_time[rip_idx]
  Kcorr <- rip_ref_dt/rip_dt
  return(Kcorr)
}


gcims_align_dt <- function(aux_list, Kcorr) {

  # Extract drift time axis
  drift_time <- aux_list$data$drift_time
  # Extract data matrix
  aux <- as.matrix(aux_list$data$data_df)

  if(Kcorr == 1){
    # If Kcorr == 1 the correction is not needed
    aux_corr <- aux
  } else {
    aux_corr <- 0*aux
    drift_time_corr <- Kcorr * drift_time
    for (j in seq_len(ncol(aux))) {
      aux_corr[,j] <- signal::interp1(drift_time_corr, aux[,j], drift_time, extrap = TRUE)
    }
  }
  aux_list$data$data_df <- round(aux_corr)
  return(aux_list)
}


# RETENTION TIME

#' This function provides the index corresponding to the reference Reactant Ion Chromatogram (RIC) to correct
#' misalignments in retention time.
#'
#' @param rics            A matrix. Each row correspond to a different RIC. There are as many RICs as samples.
#' @export
#' @importFrom ptw bestref
#' @return  An Integer number that indicates the reference sample.
#' @export
#'
find_reference_ric <- function(rics){
  ref_ric_sample_idx <- ptw::bestref(rics)$best.ref
  return(ref_ric_sample_idx)
}


#' This function finds the optimal polynomial degree to apply ptw to GMIMS data
#' when correcting retention time axis. The figure of merit is the mean of the
#' correlation between the reference RIC and samples RIC.
#'
#'
#' @param rics                         Matrix. One Reactant Ion Chromatogram per row.
#' @param ref_ric_sample_idx           A scalar. It indicates the sample that is choosen
#'                                     as reference for aligning in retention time.
#' @export
#' @importFrom ptw ptw
#' @importFrom stats cor
#' @return  An Integer between 0 and 5. If zero, no correction
#'          in retention time is needed. Between 1 and 5 it indicates
#'          the polynomial degree of the warping function.
#' @export

gcims_optimize_polynomial <- function(rics, ref_ric_sample_idx) {

  # Select the reference RIC
  ric_ref <- rics[ref_ric_sample_idx]
  # Type of correction
  correction_type_vector <- c(0, 1, 2, 3, 4, 5)
  # List of initial values fo the polynomial coefficients
  init_coeff_list <- list(c(0, 1), c(0, 1, 0), c(0, 1, 0, 0), c(0, 1, 0 , 0, 0), c(0, 1, 0 , 0, 0, 0))
  # Initialization of the correlation matrix
  corr <- matrix(0,nrow = length(samples), ncol = length(init_coeff_list) + 1)

  # Compute correlation between each sample RIC and reference RIC, for the different initial value coefficients
  xi <- seq_len(dim(rics)[2])
  for (i in seq_len(dim(rics)[1])){
    ric_sample <- rics[i, ]
    corr[i, 1] <- stats::cor(ric_ref, ric_sample, use ="complete.obs")
    for (j in seq_along(init_coeff_list)){
      corr[i,j + 1] <- stats::cor(ric_ref, ptw::ptw(ref = ric_ref, samp = ric_sample, init.coef = init_coeff_list[[j]])$warped.sample[, xi], use ="complete.obs")
    }
  }
  # Look for the global warping function that corrects the set of sample RICs with respect the reference
  correction_type_index <- which.max(apply(corr, MARGIN = 2, FUN = mean))
  correction_type <- correction_type_vector[correction_type_index]
  return(correction_type)
}


#' Correction of retention time axis using Parametric Time Warping (PTW)
#'
#' This chromatrogram computes the RIC of the reference samples and
#' uses it to correct each of the EIC of a sample, for the all the selected
#' samples on the dataset. The degree ot the polynomial that relates the
#' retention time axes of the reference and the sample to be corrected
#' is optimized using the average correlation between samples and reference as a figure
#' of merit.
#'
#' @param aux_list                       Object containing GCIMS data and metadata (one sample).
#' @param ric_ref                        A vector. The reference Reactant Ion Chromatogram.                                        as reference for aligning in retention time.
#' @param correction_type                Numeric. Integer between 0 and 5. If zero, no correction
#'                                       in retention time is needed. Between 1 and 5 it indicates
#'                                       the polynomial degree of the warping function.
#' @return And object containing GCIMS data and metadata (one sample).
#' @export
#' @importFrom ptw ptw
#' @importFrom signal interp1
#'
gcims_align_rt <- function(aux_list, ric_ref, correction_type) {
  init_coeff_list <- list(c(0, 1), c(0, 1, 0), c(0, 1, 0, 0), c(0, 1, 0 , 0, 0), c(0, 1, 0 , 0, 0, 0))
  aux <- as.matrix(aux_list$data$data_df)

  if(correction_type == 0){
    # The sample is already aligned, so aux is not modified

  } else{
    # Select the starting coefficients of the polynomial for the optimal case
    init_coeff <- init_coeff_list[[correction_type]]
    # Compute sample RIC
    ric_sample <- compute_ric(aux) # This function is found in data_preparation_functions.R file
    # Correct retention time axis usic the reference RIC and sample RIC.
    xi <- seq_len(dim(aux)[2])
    x <- ptw::ptw(ref = ric_ref, samp = ric_sample, init.coef = init_coeff)$warp.fun[, xi]
    # Interpolate data after the correction
    aux <- t(apply(aux, MARGIN = 1, FUN = signal::interp1, x = x, xi = xi, extrap = "extrap"))
    aux_list$data$data_df <- round(aux)
  }
}



















#' Multiplicative correction to the drift time axis
#'
#' This spectra alignment function computes the largest peak from
#' each spectrum, assuming it corresponds to the Reactant Ion Peak.
#' Then it computes the RIP position of all the samples and
#' corrects the all the spectra, for each sample.
#' The correction is multiplicative, meaning:
#' drift_time_corr = Kcorr * drift_time
#'
#' Changes in pressure and temperature during measurements may lead to
#' correction factors between 0.9 and 1.1.
#'
#' This alignment technique can be used when the pressure and temperature
#' measurements are not available.
#'
#' @param dir_in             Input directory. Where input data files are loaded
#'  from.
#' @param dir_out            Output directory. Where aligned data files are stored.
#' @param samples            Numeric vector of integers. Identifies the set of
#'                           sample to be visualized from the dataset.
#' @export
#'
align_td <- function(dir_in, dir_out,samples) {
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  aux_string <- paste0("M0.rds")
  aux_list <- readRDS(file.path(dir_in, aux_string))
  saveRDS(aux_list, file = file.path(dir_out, paste0("M0.rds")))
  aux <- as.matrix(aux_list$data$data_df)
  drift_time <- aux_list$data$drift_time
  #rip_position <- apply(aux, MARGIN = 2, which.max)

  aux2 <- rowSums(aux)
  rip_pos_ind_ref <- which.max(aux2)
  rip_pos_dt_ref <- drift_time[rip_pos_ind_ref]

  m <- 0
  for (i in  samples){
    m <- m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(file.path(dir_in, aux_string))
    aux <- as.matrix(aux_list$data$data_df)
    aux2 <- rowSums(aux)

    if (i == 0){
      drift_time <- aux_list$data$drift_time
      rip_pos_ind_ref <- which.max(aux2)
      rip_pos_dt_ref <- drift_time[rip_pos_ind_ref]
    } else {
      rip_pos_ind <- which.max(aux2)
      rip_pos_dt <- drift_time[rip_pos_ind]
      Kcorr <- rip_pos_dt_ref/rip_pos_dt

      if(Kcorr == 1){
        aux_corr <- aux
      } else {
        aux_corr <- 0*aux
        drift_time_corr <- Kcorr * drift_time
        for (j in seq_len(ncol(aux))) {
          aux_corr[,j] <- interp1(drift_time_corr, aux[,j], drift_time, extrap = TRUE)
        }

      }
      aux_list$data$data_df <- round(aux_corr)
      saveRDS(aux_list, file = file.path(dir_out, paste0("M", i, ".rds")))
    }


  }
}



#' Correction of retention time axis using Parametric Time Warping (PTW)
#'
#' This chromatrogram computes the RIC of the reference samples and
#' uses it to correct each of the EIC of a sample, for the all the selected
#' samples on the dataset. The degree ot the polynomial that relates the
#' retention time axes of the reference and the sample to be corrected
#' is optimized using the average correlation between samples and reference as a figure
#' of merit.
#'
#' @param dir_in             Input directory. Where input data files are loaded
#'  from.
#' @param dir_out            Output directory. Where aligned data files are stored.
#' @param samples            Numeric vector of integers. Identifies the set of
#'                           sample to be visualized from the dataset.
#' @param correction_type    Numeric. Integer between 0 and 5. If zero, no correction
#'                           in retention time is needed. Between 1 and 5 it indicates
#'                           the polynomial degree of the warping function.
#'
#' @export
#' @importFrom ptw ptw
#' @importFrom signal interp1
#' @export
align_tr <- function(dir_in, dir_out, samples, correction_type) {
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  init_coeff_list <- list(c(0, 1), c(0, 1, 0), c(0, 1, 0, 0), c(0, 1, 0 , 0, 0), c(0, 1, 0 , 0, 0, 0))

  if(correction_type == 0){
    stop("You don't need to align samples in retention time since they already are")
  }else{
    init_coeff <- init_coeff_list[[correction_type]]
  }

  aux_string <- paste0("M0.rds")
  aux_list <- readRDS(file.path(dir_in, aux_string))
  aux <- as.matrix(aux_list$data$data_df)
  saveRDS(aux_list, file = file.path(dir_out, paste0("M0.rds")))

  compute_ric <- function(x){
    ric_pos <- which.max(rowSums(x))
    ric <- x[ric_pos, ]
    ric <- max(ric) - ric
    ric <- ric/sum(ric)
    return(ric)
  }

  ric_ref <- compute_ric(aux)

  m <- 0
  for (i in samples){
    m <- m + 1
    print(paste0("Sample ", m, " of ", length(samples)))
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(file.path(dir_in, aux_string))
    aux <- as.matrix(aux_list$data$data_df)
    ric_sample <- compute_ric(aux)
    xi <- seq_len(dim(aux)[2])
    x <- ptw::ptw(ref = ric_ref, samp = ric_sample, init.coef = init_coeff)$warp.fun[, xi]
    # ric_corr <- ptw::ptw(ref = ric_ref, samp = ric_sample, init.coef =c(0, 1))$warped.sample[, xi]
    # plot(ric_ref, col = "black", type = "l")
    # lines(ric_sample, col = "red")
    # lines(ric_corr, col = "blue")
    # Sys.sleep(2)
    aux <- t(apply(aux, MARGIN = 1, FUN = signal::interp1, x = x, xi = xi, extrap = "extrap"))
    aux_list$data$data_df <- round(aux)
    saveRDS(aux_list, file = file.path(dir_out, paste0("M", i, ".rds")))
  }
}



#' This function finds the optimal polynomial degree to apply ptw to GMIMS data
#' when correcting retention time axis. The figure of merit is the mean of the
#' correlation between the reference RIC and samples RIC.
#'
#'
#' @param dir_in             Input directory. Where input data files are loaded
#'  from.
#' @param samples            Numeric vector of integers. Identifies the set of
#'                           sample to be visualized from the dataset.
#' @export
#' @importFrom ptw ptw
#' @importFrom stats cor
#' @return  An Integer between 0 and 5. If zero, no correction
#'          in retention time is needed. Between 1 and 5 it indicates
#'          the polynomial degree of the warping function.
#' @export

optimize_align_tr <- function(dir_in, samples) {
  correction_type_vector <- c(0, 1, 2, 3, 4, 5)
  aux_string <- paste0("M0.rds")
  aux_list <- readRDS(file.path(dir_in, aux_string))
  aux <- as.matrix(aux_list$data$data_df)

  compute_ric <- function(x){
    ric_pos <- which.max(rowSums(x))
    ric <- x[ric_pos, ]
    ric <- max(ric) - ric
    ric <- ric/sum(ric)
    return(ric)
  }

  ric_ref <- compute_ric(aux)
  init_coeff_list <- list(c(0, 1), c(0, 1, 0), c(0, 1, 0, 0), c(0, 1, 0 , 0, 0), c(0, 1, 0 , 0, 0, 0))
  ##wcc <- matrix(0,nrow = length(samples), n = length(init_coeff_list))
  corr <- matrix(0,nrow = length(samples), ncol = length(init_coeff_list) + 1)



  m <- 0
  for (i in  samples){
    m <- m + 1
    if (m != 0){
      print(paste0("Sample ", m, " of ", length(samples)))
    }
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(file.path(dir_in, aux_string))
    aux <- as.matrix(aux_list$data$data_df)
    ric_sample <- compute_ric(aux)
    xi <- seq_len(dim(aux)[2])
    corr[m, 1] <- stats::cor(ric_ref, ric_sample, use ="complete.obs")
    for (j in seq_along(init_coeff_list)){
     corr[m,j + 1] <- stats::cor(ric_ref, ptw::ptw(ref = ric_ref, samp = ric_sample, init.coef = init_coeff_list[[j]])$warped.sample[, xi], use ="complete.obs")
    }
  }

  #print(wcc)

  #plot(apply(wcc, MARGIN = 2, FUN = mean), type = "l")
   correction_type_index <- which.max(apply(corr, MARGIN = 2, FUN = mean))
   correction_type <- correction_type_vector[correction_type_index]
   return(correction_type)


}


