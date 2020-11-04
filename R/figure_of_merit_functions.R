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

  #@importFrom raster raster extend xyFromCell
  #@param nxn_window        Size of the nxn window
  setwd(dir_in)
  #m <- 0
  #for (i in samples){
  # m <- m + 1
  m <-1
  print(paste0("Computing the peaknum of sample ", samples[m]))
  aux_string <- paste0("M", samples[m], ".rds")
  aux_list <- readRDS(aux_string) #new
  aux <- t(as.matrix(aux_list$data$data_df))

  #by_rows
  desired_length <-50
  peak_list_td <- vector(mode = "list", length = desired_length)
  limite <- 0.0153

  for (j in 1:50){
  peak_list_td[[j]] <- list(index_tr = j, peak_info = findpeaks(as.vector(aux[, j]))[,c(1,2)])
  }


  #by_cols
  desired_length <- 200
  peak_list_tr <- vector(mode = "list", length = desired_length)
  limite <- 0.0153
  for (j in 1:200){
    peak_list_tr[[j]] <- list(index_td = j, peak_info = findpeaks(as.vector(t(aux[j, ])))[,c(1,2)])
  }

  # r <- raster(aux)
  # extent(r) <- extent(c(0, dim(aux)[1], 0, dim(aux)[2]) + 0.5)
  #
  # ## Find the maximum value within the 9-cell neighborhood of each cell
  # f <- function(X) max(X, na.rm=TRUE)
  # ww <- matrix(1, nrow= nxn_window, ncol = nxn_window) ## Weight matrix for cells in moving window
  # localmax <- focal(r, fun=f, w=ww, pad=TRUE, padValue=NA)
  #
  # ## Does each cell have the maximum value in its neighborhood?
  # r2 <- r==localmax
  #
  # ## Get x-y coordinates of those cells that are local maxima
  # maxXY <- xyFromCell(r2, Which(r2==1, cells=TRUE))
  # return(maxXY)


}
