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


