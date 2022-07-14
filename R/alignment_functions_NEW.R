#' Multiplicative correction to the drift time axis
#'
#' This spectra alignment function computes the largest peak from
#' each spectrum, assuming it corresponds to the Reactant Ion Peak.
#' Then it computes the median RIP position of all the samples and
#' corrects the drift time for each sample.
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
#' @importFrom stats median
#' @importFrom signal interp1
#' @importFrom pracma findpeaks
#'
gcims_align <- function(dir_in, dir_out,samples) {
  aux_string <- paste0("M0.rds")
  aux_list <- readRDS(file.path(dir_in, aux_string))
  aux <- as.matrix(aux_list$data$data_df)
  drift_time <- aux_list$data$drift_time
  print(length(drift_time))
  rip_position <- apply(aux, MARGIN = 2, which.max)


  dt_rip_position <- drift_time[rip_position]
  drift_time_rip_ref <- median(dt_rip_position, na.rm = TRUE)
  drift_time_rip_ref_pos <- which.min(abs(drift_time - drift_time_rip_ref))
  print(drift_time_rip_ref_pos)

  m <- -1
  for (i in  samples){
    m <- m + 1
    if (m != 0){
      print(paste0("Sample ", m, " of ", length(samples)))
    }
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(file.path(dir_in, aux_string))
    aux <- as.matrix(aux_list$data$data_df)


    # Compute the total ion spectrum
    aux_2 <- rowSums(aux)
    peaks_info <- findpeaks(aux_2)

    # look for rip_position of total ion spectrum
    rip_pos_ind_sum <- which.max(peaks_info[ , 1])
    rip_pos_sum <- peaks_info[rip_pos_ind_sum, 2]

    # Look for the beginning and ending of the RIP (searching the closest minima to it)
    valleys_info <- findpeaks(-aux_2)
    valleys_pos_sum <- valleys_info[ , 2]
    closest_valley_ind_sum <- which.min(abs(valleys_pos_sum - rip_pos_sum))

    # Select the RIP region
    if(valleys_pos_sum[closest_valley_ind_sum] < rip_pos_sum){
      rip_bounds <- valleys_pos_sum[c(closest_valley_ind_sum,closest_valley_ind_sum + 1)]
    } else if (valleys_pos_sum[closest_valley_ind_sum] > rip_pos_sum){
      if(closest_valley_ind_sum > 1){
        rip_bounds <- valleys_pos_sum[c(closest_valley_ind_sum - 1,closest_valley_ind_sum)]
      } else {
        rip_bounds <- valleys_pos_sum[c(1, closest_valley_ind_sum)]
      }
    }

    # RIP position (in indexes) for all spectra in a sample
    rip_position <- apply(aux[rip_bounds[1]:rip_bounds[2],], MARGIN = 2, which.max) + rip_bounds[1] - 1
    Kcorr <- drift_time[drift_time_rip_ref_pos]/drift_time[rip_position]
    aux_corr <- 0*aux
    for (j in seq_len(ncol(aux))) {
      drift_time_corr <- Kcorr[j] * drift_time
      aux_corr[,j] <- interp1(drift_time_corr, aux[,j], drift_time, extrap = TRUE)
    }
    aux_list$data$data_df <- round(aux_corr)
    saveRDS(aux_list, file = file.path(dir_out, paste0("M", i, ".rds")))
  }
}
