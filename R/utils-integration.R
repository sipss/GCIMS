find_regions_rip_saturated <- function(aux, rip_saturation_threshold, verbose = FALSE, retention_time = NULL, drift_time = NULL) {
  the_rip <- find_rip(aux, verbose = verbose, retention_time = retention_time, drift_time = drift_time)
  # Search saturation regions
  rip_region <- aux[the_rip$dt_idx_start:the_rip$dt_idx_end, , drop = FALSE]
  rip_chrom <- rowSums(rip_region) / nrow(rip_region)
  rt_rip_saturated_indices <- which(rip_chrom <= rip_saturation_threshold * max(rip_chrom))
  if (length(rt_rip_saturated_indices) == 0L) {
    rt_saturated_regions <- matrix(nrow = 0L, ncol = 2L)
  } else {
    saturation_list <- split(rt_rip_saturated_indices, cumsum(c(1, diff(rt_rip_saturated_indices)) != 1))
    rt_saturated_regions <- matrix(0, nrow = length(saturation_list), ncol = 2)
    for (k in seq_along(saturation_list)) {
      rt_saturated_regions[k, 1L] <- min(saturation_list[[k]])
      rt_saturated_regions[k, 2L] <- max(saturation_list[[k]])
    }
  }
  colnames(rt_saturated_regions) <- c("begin", "end")
  if (!is.null(retention_time)) {
    rt_saturated_regions_s <- matrix(0.0, nrow = length(saturation_list), ncol = 2L)
    rt_saturated_regions_s[,1L] <- retention_time[rt_saturated_regions[,1L]]
    rt_saturated_regions_s[,2L] <- retention_time[rt_saturated_regions[,2L]]
    return(rt_saturated_regions_s)
  } else {
    return(rt_saturated_regions)
  }
}

#' Find half max boundaries and FWHM metrics
#'
#' @noRd
#' @param x A numeric vector with ONE peak
#'
#' @return A list with the left and right half maximum limits of the peak,
#' the apex and widths according to the half max criteria.
find_half_max_boundaries <- function(x) {
  max_idx <- which.max(x)
  half_max <- x[max_idx]/2
  #

  left_bound_idx <- findZeroCrossings(x[seq_len(max_idx)] - half_max, direction = "up")
  if (length(left_bound_idx) == 0L) {
    left_bound_idx <- NA_integer_
  } else if (length(left_bound_idx) > 1L) {
    left_bound_idx <- left_bound_idx[length(left_bound_idx)]
  }

  right_bound_idx <- max_idx - 1L + findZeroCrossings(x[max_idx:length(x)] - half_max, direction = "down")
  if (length(right_bound_idx) == 0L) {
    right_bound_idx <- NA_integer_
  } else if (length(right_bound_idx) > 1L) {
    right_bound_idx <- right_bound_idx[1L]
  }

  half_width_left <- max_idx - left_bound_idx
  half_width_right <- right_bound_idx - max_idx
  list(
    left_bound_idx = left_bound_idx,
    apex_idx = max_idx,
    right_bound_idx = right_bound_idx,
    half_width_left = half_width_left,
    half_width_right = half_width_right,
    fwhm = half_width_left + half_width_right,
    fwhm_sym = 2*mean(c(half_width_left, half_width_right), na.rm = TRUE)
  )
}

#' Estimates RIP region
#'
#'
#' Estimate RIP region, including the indices where the RIP starts and ends in drift time and
#' the index in retention time where the RIP has a higher intensity
#' @noRd
#' @param intensity_mat A matrix of intensities of shape (drift, retention).
#' @return A named list with the positions where the RIP starts and ends in drift time indices,
#' and the positions of the RIP apex.
find_rip <- function(intensity_mat, verbose = FALSE, retention_time = NULL, drift_time = NULL) {
  # Simplified version from pracma findpeaks:
  simplest_findpeaks <- function(x)
  {
    xc <- paste(as.character(sign(diff(x))), collapse = "")
    xc <- gsub("1", "+", gsub("-1", "-", xc))
    rc <- gregexpr("[+]+[-]+", xc)[[1]]
    if (rc[1] < 0)
      return(NULL)
    x1 <- rc
    x2 <- rc + attr(rc, "match.length")
    attributes(x1) <- NULL
    attributes(x2) <- NULL
    n <- length(x1)
    xv <- xp <- numeric(n)
    for (i in 1:n) {
      xp[i] <- which.max(x[x1[i]:x2[i]]) + x1[i] - 1
      xv[i] <- x[xp[i]]
    }
    inds <- which(xv - pmax(x[x1], x[x2]) >= 0)
    X <- cbind(xv[inds], xp[inds], x1[inds], x2[inds])
    if (length(X) == 0)
      return(c())
    return(X)
  }

  total_ion_spectrum <- rowSums(intensity_mat)
  # In most rows the RIP will be the higher peak, so it will also be the maximum in the TIS
  dt_idx_rip <- which.max(total_ion_spectrum)
  # The retention time where the RIP is maximum:
  rt_idx_apex <- which.max(intensity_mat[dt_idx_rip,])
  # The valleys of the TIS. We bound them by the TIS limits:
  dt_idx_valleys <- c(1L, simplest_findpeaks(-total_ion_spectrum)[, 2L], length(total_ion_spectrum))
  # The RIP
  dt_idx_start <- dt_idx_valleys[max(which(dt_idx_rip - dt_idx_valleys > 0))] # Find starting index of RIP
  dt_idx_end <- dt_idx_valleys[min(which(dt_idx_valleys - dt_idx_rip > 0))] # Find ending index of RIP

  # The valley detection is not perfect. Sometimes the actual peak boundary is missed
  # To workaround this, we estimate the FWHM and we clip the boundaries at the apex +- 2*FWHM
  bounds_and_widths <- find_half_max_boundaries(total_ion_spectrum)
  if (is.na(bounds_and_widths$fwhm_sym)) {
    abort(
      message = c(
        "The ROI selection could not locate the RIP width",
        "i" = "Please contact the GCIMS authors at https://github.com/sipss/GCIMS/issues and if possible submit them your sample")
    )
  }
  # Clip the peak boundaries:
  dt_idx_start <- max(dt_idx_start, dt_idx_rip - 3*bounds_and_widths$fwhm_sym)
  dt_idx_end <- min(dt_idx_end, dt_idx_rip + 3*bounds_and_widths$fwhm_sym)

  if (verbose) {
    if (!is.null(retention_time)) {
      rt_units <- "s"
    } else {
      rt_units <- "pts"
      retention_time <- seq_len(ncol(intensity_mat))
    }
    if (!is.null(drift_time)) {
      dt_units <- "ms"
    } else {
      dt_units <- "pts"
      drift_time <- seq_len(nrow(intensity_mat))
    }
    dt_ms_start <- drift_time[dt_idx_start]
    dt_ms_end <- drift_time[dt_idx_end]
    dt_ms_apex <- drift_time[dt_idx_rip]
    rt_s_apex <- retention_time[rt_idx_apex]

    inform(message = c(
      "RIP was detected",
      "i" = glue(" - At drift time: [{dt_ms_start} - {dt_ms_end}] {dt_units}"),
      "i" = glue(" - Maximum RIP intensity at: (dt: {dt_ms_apex} {dt_units}, rt: {rt_s_apex} {rt_units})")
    )
    )
  }
  list(
    dt_idx_start = dt_idx_start,
    dt_idx_apex = dt_idx_rip,
    dt_idx_end = dt_idx_end,
    rt_idx_apex = rt_idx_apex,
    rip = intensity_mat[dt_idx_start:dt_idx_end, rt_idx_apex]
  )
}
