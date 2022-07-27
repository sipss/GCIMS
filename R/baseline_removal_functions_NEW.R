#' Remove baseline of a dataset
#'
#' @param dir_in          Input directory. Where input data files are
#'   loaded from.
#' @param dir_out         Output directory. Where data files after baseline
#'   correction are stored.
#' @param samples         A vector. Set of samples to which remove the baseline
#'   (e.g.: c(1, 2, 3)).
#' @param peak_list       A peak list after containing information regarding the 2d peaks found in data.
#'
#' @details `gcims_remove_baseline` performs baseline correction on sample
#'   spectra (in drift time) and chromatograms (in retention time) .This function corrects
#'   the baseline in drift time and afterwards, in retention time. To do it, it cuts the
#'   axis in regions of equal length and looks for the local minima of each region. A simple
#'   linear interpolation is performed to create a baseline of the same length as the number
#'   of data points per axis.The region width is estimated differently in each of the the axes:
#'   In drift time, it is computed as 12 times the standard deviation of the RIP width, while
#'   in retention time as the third quantile of peak width distributions in chromatograms (obtained
#'   from the peak table). Take into account that when removing the baseline in retention
#'   time, the shape of the Reactant Ion Peak is severely affected.
#'
#'
#' @family Baseline removal functions
#' @return A set of S3 objects.
#' @export
#' @importFrom signal interp1
#' @importFrom stats  quantile
#

gcims_remove_baseline <- function(dir_in, dir_out, samples,
                                  peak_list){


  #-------------#
  #     MAIN    #
  #-------------#


  print(" ")
  print("  /////////////////////////")
  print(" /    Baseline removal   /")
  print("/////////////////////////")
  print(" ")

  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  for (i in samples) {
    aux_string <- paste0("M", i, ".rds")
    aux_list <- readRDS(file.path(dir_in, aux_string)) #new
    aux <- as.matrix(aux_list$data$data_df)
    # correct drift time
    aux <- remove_bsln_td(aux)
    # correct retention time
    aux <- remove_bsln_tr(aux, peak_list)
    aux_list$data$data_df <- round(aux)
    saveRDS(aux_list, file = file.path(dir_out, paste0("M", i, ".rds")))
  }
}



remove_bsln_td <- function(aux,sig_mult = 12){
  # Create the total ion spectrum
  y <- rowSums(aux)
  # Compute the maximum height or RIP
  y_max <- max(y)
  # Find RIP position in drift time index
  mu <- which.max(y)
  # f/2 -> FWHM = 2.3482 * sigma
  fwhm <-  2 * abs(mu - which.min(abs(y - (y_max / 2))))
  # Estimate sigma
  sig <- fwhm / 2.3482
  #Create vector of indexes
  x <- seq_along(y)
  #Selection of the region size to perform interpolation
  region_size <- sig * sig_mult
  number_of_regions <- ceiling(length(y) / region_size)
  baseline <- compute_baseline(aux, x, number_of_regions, region_size)
  aux_corr <-  aux - baseline
  return(aux_corr)
}

remove_bsln_tr <- function(aux, peak_list){
  peak_width <- peak_list$rt_max_idx - peak_list$rt_min_idx
  # Estimate region
  region_size <- stats::quantile(peak_width, 0.75)
  # Transpose matrix
  aux <- t(aux)
  # Create the total ion chromatogram
  y <- rowSums(aux)
  # Create vector of indexes
  x <- seq_along(y)
  #Selection of the region size to perform interpolation
  number_of_regions <- ceiling(length(y) / region_size)
  #compute baseline
  baseline <- compute_baseline(aux, x, number_of_regions, region_size)
  # correct baseline and transpose
  aux_corr <-  t(aux - baseline)
  return(aux_corr)
}




compute_baseline <- function(aux, x, number_of_regions, region_size){
  local_minima_idx <- numeric(number_of_regions)
  local_minima_intensities <- numeric(number_of_regions)
  baseline <- 0*aux
  for (j in seq_len(dim(aux)[2])){
    z <- aux[, j]
    for (i in 1:number_of_regions) {
      first_point_i <- (i - 1)*region_size + 1
      last_point_i <- first_point_i + region_size - 1
      last_point_i <- min(last_point_i, length(z))
      x_i <- x[first_point_i:last_point_i]
      z_i <- z[first_point_i:last_point_i]
      # We find the index of the minimum in the region:
      minimum_index_i <- which.min(z_i)
      # And we save the  value and the intensity where the local minimum is:
      local_minima_idx[i] <- x_i[minimum_index_i]
      local_minima_intensities[i] <- z_i[minimum_index_i]
    }

    # To create the full baseline, we repeat the minimum value throughout its region
    local_minima_baseline <- rep(local_minima_intensities, each = region_size)
    # The last region may be shorter, so we just cut the vector at the end:
    local_minima_baseline <- local_minima_baseline[1:length(z)]
    baseline[, j] <-signal::interp1(x = local_minima_idx, y = local_minima_intensities, xi = seq_len(dim(aux)[1]),method = "linear", extrap = TRUE)
  }
  return(baseline)
}









