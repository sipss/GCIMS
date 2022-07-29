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
#' @param alignment_data     A list containing two matrices: tis and rics, the Total
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
gcims_align_data <- function(dir_in, dir_out, samples, alignment_data){

  print(" ")
  print("  ////////////////////////")
  print(" /    Aligning data    /")
  print("////////////////////////")
  print(" ")

  # Crete folder to store data
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  # For drift time:
  # Reference position of rip in drift time (in indexes)
  tis <- alignment_data$tis
  rip_ref_idx <- round(stats::median(apply(X = tis, MARGIN = 1, FUN = which.max), na.rm = TRUE))
  Kcorr_samples <- 0 * samples

  # For retention time:
  # Sample idx to correct for retention time using the ric to align
  rics <- alignment_data$rics
  ref_ric_sample_idx <- find_reference_ric(rics)
  # Search for the optimal polynomial degree of the warping function:
  correction_type_vector <- gcims_optimize_polynomial(rics, ref_ric_sample_idx)
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

    Kcorr_samples[m] <- gcims_compute_Kcorr(aux_list, rip_ref_idx)
    # Apply the correction to data:
    aux_list <- gcims_align_dt(aux_list, Kcorr_samples[m])

    # Then, we correct retention time
    # Apply the correction to data:
    aux_list <- gcims_align_rt(aux_list, ric_ref, correction_type_vector[m])

    # Finally, we save data
    saveRDS(aux_list, file = file.path(dir_out, paste0("M", i, ".rds")))
  }
  alignment_info <- list(rip_ref_idx = rip_ref_idx, Kcorr_samples = Kcorr_samples, ref_ric_sample_idx = ref_ric_sample_idx, correction_type_vector = correction_type_vector)
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
#' @return  A vector of integers with as many components as samples to be corrected, with values
#'          between 0 and 5. If zero, no correction
#'          in retention time is needed. Between 1 and 5, it indicates
#'          the polynomial degree of the warping function.
#' @export

gcims_optimize_polynomial <- function(rics, ref_ric_sample_idx) {

  # Select the reference RIC
  ric_ref <- rics[ref_ric_sample_idx, ]
  # Type of correction
  correction_type_options <- c(0, 1, 2, 3, 4, 5)
  # List of initial values fo the polynomial coefficients
  init_coeff_list <- list(c(0, 1), c(0, 1, 0), c(0, 1, 0, 0), c(0, 1, 0 , 0, 0), c(0, 1, 0 , 0, 0, 0))
  # Initialization of the correlation matrix
  corr <- matrix(1, nrow = length(samples), ncol = length(init_coeff_list) + 1)
  correction_type_vector <- 0 * seq_along(samples)

  # Compute correlation between each sample RIC and reference RIC, for the different initial value coefficients
  xi <- seq_len(dim(rics)[2])
  # This is done to avoid computing the correlation of the reference againt itself
  samples_to_correct <- setdiff(seq_len(dim(rics)[1]),ref_ric_sample_idx)
  for (i in samples_to_correct){
    ric_sample <- rics[i, ]
    corr[i, 1] <- stats::cor(ric_ref, ric_sample, use ="complete.obs")
    for (j in seq_along(init_coeff_list)){
      corr[i,j + 1] <- stats::cor(ric_ref, ptw::ptw(ref = ric_ref, samp = ric_sample, init.coef = init_coeff_list[[j]])$warped.sample[, xi], use ="complete.obs")
    }
    # Initialize index:
    idx_max <- idx_sel <- idx_zero <- idx_sign <- (length(init_coeff_list)) + 1
    # Check when correlation decreases for the first time or does not change when increasing degree of the polynomial.
    diff_corr_i <- diff(corr[i, ])
    if (all(sign(diff_corr_i) == 1)){
      # Correlation is always increasing (don't do anything)
    } else if (any(diff_corr_i == 0)){
      # Correlation is equal for at least for  different degrees of the polynomial.
      idx_zero <- which(diff_corr_i == 0)[1]
    } else if (any(sign(diff_corr_i) == -1)){
      # Correlation decreases at least for two consecutive degrees of the polynomial.
      idx_sign <- which(sign(diff_corr_i) == -1)[1]
    }
    # Combine indexes and look for the minimum
    idx_combine <- c(idx_sign, idx_zero)
    if(any(idx_combine < idx_max)){
      # Select the minimum index in which the correlation is still increasing
      idx_sel <- min(idx_combine)
    }
    correction_type_vector[i] <- correction_type_options[idx_sel]

  }
  return(correction_type_vector)
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
    ric_sample <- compute_ric(aux_list) # This function is found in data_preparation_functions.R file
    # Correct retention time axis usic the reference RIC and sample RIC.
    xi <- seq_len(dim(aux)[2])
    x <- ptw::ptw(ref = ric_ref, samp = ric_sample, init.coef = init_coeff)$warp.fun[, xi]
    # Interpolate data after the correction
    aux <- t(apply(aux, MARGIN = 1, FUN = signal::interp1, x = x, xi = xi, extrap = "extrap"))
    aux_list$data$data_df <- round(aux)
  }
  return(aux_list)
}




