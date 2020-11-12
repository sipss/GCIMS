#' Compute the signal to noise ratio
#'
#' @param dir_in            The input directory.
#' @param samples           Samples to which the SNR is computed
#'                          sample characteristics (one per sample class)
#' @return A dummy dataset to be stored in dir_out
#' @export
#' @importFrom chemometrics sd_trim
#' @importFrom stats median density
#' @family Figure of merit functions
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }

gcims_compute_snr <- function(dir_in, samples){

  #Retention_Time <- Drift_Time <- Value <- NULL

  print(" ")
  print("  /////////////////////////////////////")
  print(" /          Computing the SNR        /")
  print("/////////////////////////////////////")
  print(" ")

  signal_to_noise_ratio <-  matrix(0,nrow = length(samples), ncol = 2)
  colnames(signal_to_noise_ratio) <- c("snr", "snr_db")
  signal_to_noise_ratio <- as.data.frame(signal_to_noise_ratio)
  statistics <- matrix(0,nrow = length(samples), ncol = 2)
  colnames(statistics) <- c("mu", "sigma")
  statistics <- as.data.frame(statistics)

  noise_characterization <- list(signal_to_noise_ratio = signal_to_noise_ratio, statistics = statistics)
  rm(signal_to_noise_ratio, statistics)


  setwd(dir_in)
  m <- 0
  for (i in samples){
  m <- m + 1
  print(paste0("Computing SNR of sample ", samples[m]))
  aux_string <- paste0("M", i, ".rds")
  aux_list <- readRDS(aux_string) #new


  aux <- sort(as.vector(aux_list$data$data_df), decreasing = TRUE) #new

  #ROBUST ESTIMATION OF NOISE
  mu = mean(aux, trim = 0.2)
  sigma = sd_trim(aux, const = TRUE)

  y_limit <- mu + (3 * sigma)

  x_limit <- which.min(abs(aux - y_limit))

  # CHECKING PLOTS (commented)
  #h1 <- hist(aux, breaks = "FD")
  #h1$density=h1$counts/max(h1$counts)
  #plot(h1$mids, h1$density, type = "l", lwd = "2")
  #plot(str(h1))
  # print(h1$mids)
  # gauss <- function(x, mu, sigma){
  #   gaussian <- (1 / (sqrt(2 * pi) * sigma)) * exp(-(1/2) * ((x - mu) / sigma) ^ 2)
  # }

  # noise_gauss <- gauss(h1$mids, mu, std)
  # lines(h1$mids, noise_gauss/max(noise_gauss), type = "l", lwd = "2", col="red")
  # abline(v= mu, col = "green")
  # abline(v = (mu + 3*std), col = "orange")

  # Guess what is signal and what is noise
  aux_signal <- aux[aux > y_limit]
  aux_noise <- aux[aux <= y_limit]
  # plot(aux_signal)
  # abline(h = y_limit, col = "blue")
  # abline (v = x_limit, col = "blue")
  # plot(aux_noise)
  # plot(aux)



  #snr <- (var(aux_signal)/ var(aux_noise))
  snr <- (sum(aux_signal * aux_signal) / sum(aux_noise * aux_noise))
  snr_db <- 10*log10(snr)
  noise_characterization[[1]][m, 1:2] <- c(snr, snr_db)
  noise_characterization[[2]][m, 1:2] <- c(mu, sigma)

  #noise_est <- list(mu= mu, sigma = sigma, snr = snr, snr_db = snr_db)
  }
  saveRDS(noise_characterization, file = paste0("noise_characterization", ".rds"))
  #return(noise_est)


}


#' Compute the number of peaks in samples
#'
#' @param dir_in            The input directory.
#' @param samples           Samples to which the SNR is computed
#'                          sample characteristics (one per sample class)

#' @return A dummy dataset to be stored in dir_out
#' @export
#' @importFrom pracma findpeaks
#' @family Figure of merit functions
#' @examples
#' \dontrun{
#' dataset_2_polarities <- lcms_dataset_load(system.file("extdata",
#'                                                      "dataset_metadata.rds",
#'                                                      package = "NIHSlcms"))
#' dataset_pos <- lcms_filter_polarity(dataset_2_polarities, polarity. = 1)
#'
#' print(dataset_pos)
#' }

gcims_compute_numpeaks <- function(dir_in, samples){


  print(" ")
  print("  /////////////////////////////////////")
  print(" /   Computing the Number of Peaks   /")
  print("/////////////////////////////////////")


  setwd(dir_in)
  #m <- 0
  #for (i in samples){
  # m <- m + 1
  m <-1
  print(paste0("Computing the peaknum of sample ", samples[m]))
  aux_string <- paste0("M", samples[m], ".rds")
  aux_list <- readRDS(aux_string) #new
  aux <- t(as.matrix(aux_list$data$data_df))

  # Unfold aux
  aux_vector <- sort(as.vector(aux_list$data$data_df), decreasing = TRUE) #new

  # Compute robust estimation of noise present in a sample
  mu = mean(aux, trim = 0.2)
  sigma = sd_trim(aux_vector, const = TRUE)
  threshold <- mu + (3 * sigma)


  # Function to compute the power of a signal
  compute_power <- function (x){
    pow <- sum(x * x) / length(x)
    pow
  }

  # Compute noise power
  noise_power <- compute_power(aux_vector[aux_vector <= threshold])
  rm(aux_vector)

  # Compute a list of tentative peaks.

  # For each retention time index,
  # detect the peaks of the corresponding spectrum.

  # Create the list to contain the peaks
  rt_index_length <- dim(aux)[1]
  peak_list_td <- vector(mode = "list", length = rt_index_length)

  # Set to zero matrix values that are below the noise threshold
  aux[aux <= threshold] <- 0

  tr_index <- 1:rt_index_length
  for (j in tr_index){
    peak_list_td[[j]] <- findpeaks(as.vector(t(aux[j, ])))[, 2]
  }

  # Obtain all peaks once, that is without repetitions
  unique_peaks_td <- unique(unlist(peak_list_td)) #possible_peaks

  # Create a empty list to contain the peaks
  # (the indexes in retention time that correspond to each peak).
  grouped_peak_list_tr <- vector(mode = "list", length = length(unique_peaks_td)) #grouped_possible_peaks

  # Function to know in which retention time indexes (logical) a peak in drift time is present
  find_match <- function(x, y){
    indexes <- match(y, x, nomatch = 0) > 0
    indexes
  }
  # Fill the list of peaks applying 'find_match' to all elements in the list
  for(k in 1:length(unique_peaks_td)){
    grouped_peak_list_tr[[k]] <- which(sapply(peak_list_td, find_match, y = unique_peaks_td[k]))

  }
  # Since you can find more than one peak per retention
  # time index, it'll be necessary to split this list

  # Function to split the peak list according to which retention time indexes
  # are consecutive (they belong to a peak). Additionally, this function
  # discard peaks of peak length 0, in retention time (Note that a peak signal
  # of length 1 has only is equivalent to a pick of length 0.

  split_peaks <- function(x){
    indexes <- split(x,cumsum(c(1, diff(x) != 1)))
    flags <- lapply(indexes, function(x) length(x) > 1)
    if(any(unlist(flags)) == TRUE){
      indexes <- indexes[which(unlist(flags))]
    } else{
      indexes <- NULL
    }
    indexes
  }

  # Split the peak list (in retention time) applying split_peaks
  consecutive_grouped_peak_list_tr <- sapply(grouped_peak_list_tr, split_peaks)


#
  grouped_peak_list_td <- vector(mode = "list", length = length(consecutive_grouped_peak_list_tr))
  for (k in 1:length(consecutive_grouped_peak_list_tr)){
    condition <- length(consecutive_grouped_peak_list_tr[[k]]) == 0
    if (condition == TRUE){
      grouped_peak_list_td[[k]] <- NULL
    } else {
      grouped_peak_list_td[[k]] <- as.list(rep(unique_peaks_td[[k]], length(consecutive_grouped_peak_list_tr[[k]])))
    }

  }


  consecutive_grouped_peak_list_tr[sapply(consecutive_grouped_peak_list_tr, is.null)] <- NULL
  grouped_peak_list_td[sapply(grouped_peak_list_td, is.null)] <- NULL


  peak_length_tr <-  vector(mode = "list", length = length(consecutive_grouped_peak_list_tr))
  for (k in 1:length(consecutive_grouped_peak_list_tr)){
    peak_length_tr[[k]] <- lapply(consecutive_grouped_peak_list_tr[[k]], function(x) length(x) - 1)
  }



  peak_max_value <- peak_max_pos_abs <- peak_max_pos_rel <- peak_asymm_tr <- vector(mode = "list", length = length(consecutive_grouped_peak_list_tr))

  for (k in 1:length(consecutive_grouped_peak_list_tr)){
    length_list <- length(consecutive_grouped_peak_list_tr[[k]])
    for (l in 1: length_list){
      peak_max_pos_rel[[k]][[l]] <- which.max(aux[consecutive_grouped_peak_list_tr[[k]][[l]], grouped_peak_list_td[[k]][[l]]])
      peak_max_pos_abs[[k]][[l]] <- consecutive_grouped_peak_list_tr[[k]][[l]][peak_max_pos_rel[[k]][[l]]]
      peak_max_value[[k]][[l]]   <- aux[peak_max_pos_abs[[k]][[l]], grouped_peak_list_td [[k]][[l]]]
      left_peak_side <- abs(min(consecutive_grouped_peak_list_tr[[k]][[l]]) - peak_max_pos_abs[[k]][[l]])
      right_peak_side <- abs(max(consecutive_grouped_peak_list_tr[[k]][[l]]) - peak_max_pos_abs[[k]][[l]])
      peak_asymm_tr[[k]][[l]] <- (right_peak_side / left_peak_side )

    }
  }



peak_length_td <- peak_asymm_td <- vector(mode = "list", length = length(peak_max_pos_abs))
  for (k in 1:length(peak_max_pos_abs)){
    length_list <- length(peak_max_pos_abs[[k]])
    for (l in 1: length_list){
      zero_indexes <- which(aux[peak_max_pos_abs[[k]][[l]], ] == 0)
      left_cond <- which(zero_indexes < grouped_peak_list_td[[k]][[l]])
      left_peak_side_td <- abs((max(zero_indexes[left_cond]) + 1) - grouped_peak_list_td[[k]][[l]])
      right_cond <- which(zero_indexes > grouped_peak_list_td[[k]][[l]])
      right_peak_side_td  <- abs((min(zero_indexes[right_cond]) - 1) - grouped_peak_list_td[[k]][[l]])
      peak_length_td[[k]][[l]] <-  right_peak_side_td  + left_peak_side_td
      peak_asymm_td[[k]][[l]] <-   right_peak_side_td / left_peak_side_td
    }
  }

  m <- 0
  final_peak_list <-vector(mode = "list", length = length(unlist(peak_max_pos_abs, recursive = FALSE)))
  for (k in 1:length(peak_max_pos_abs)){
    length_list <- length(peak_max_pos_abs[[k]])
    for (l in 1: length_list){
      m <- m + 1
      final_peak_list[[m]] <- c(grouped_peak_list_td[[k]][[l]],
                                      peak_max_pos_abs[[k]][[l]],
                                      peak_max_value[[k]][[l]],
                                      peak_length_tr[[k]][[l]],
                                      peak_length_tr[[k]][[l]],
                                      peak_asymm_tr[[k]][[l]],
                                      peak_asymm_td[[k]][[l]])

    }
  }

  final_peak_list <- matrix(unlist(final_peak_list), nrow = length(unlist(peak_max_pos_abs, recursive = FALSE)),
                            ncol = length(unlist(final_peak_list))/length(unlist(peak_max_pos_abs, recursive = FALSE)),
                            byrow = TRUE)
  colnames(final_peak_list)<-c("rt_ind", "dt_ind", "value", "rt_length","dt_legth", "rt_asymm", "dt_asymm")

  final_peak_list <- as.data.frame(final_peak_list)



   return(final_peak_list)


}
