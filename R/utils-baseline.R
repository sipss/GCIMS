estimate_baseline_td <- function(aux, sig_mult = 12, peak_fwhm_pts = NULL) {
  if (is.null(peak_fwhm_pts)) {
    # Create the total ion spectrum
    y <- rowSums(aux)
    # Compute the maximum height or RIP
    y_max <- max(y)
    # Find RIP position in drift time index
    mu <- which.max(y)
    # f/2 -> FWHM = 2.3482 * sigma
    peak_fwhm_pts <-  2 * abs(mu - which.min(abs(y - (y_max / 2))))
    # Estimate sigma
  }
  sig <- peak_fwhm_pts / 2.3482
  #Create vector of indexes
  x <- seq_len(nrow(aux))
  #Selection of the region size to perform interpolation
  region_size <- sig * sig_mult
  number_of_regions <- ceiling(nrow(aux) / region_size)
  baseline <- compute_baseline(aux, x, number_of_regions, region_size)
  return(baseline)
}

estimate_baseline_tr <- function(aux, region_size) {
  # Transpose matrix
  aux <- t(aux)
  # Create vector of indexes
  x <- seq_len(nrow(aux))
  #Selection of the region size to perform interpolation
  number_of_regions <- ceiling(nrow(aux)/region_size)
  #compute baseline
  baseline <- compute_baseline(aux, x, number_of_regions, region_size)
  # Transpose again
  t(baseline)
}


compute_baseline <- function(aux, x, number_of_regions, region_size){
  local_minima_idx <- numeric(number_of_regions)
  local_minima_intensities <- numeric(number_of_regions)
  baseline <- 0*aux
  for (j in seq_len(ncol(aux))) {
    z <- aux[, j]
    for (i in seq_len(number_of_regions)) {
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
    local_minima_baseline <- local_minima_baseline[seq_len(length(z))]
    baseline[, j] <- signal::interp1(x = local_minima_idx, y = local_minima_intensities, xi = seq_len(dim(aux)[1]),method = "linear", extrap = TRUE)
  }
  return(baseline)
}



