#' ROIs Selection

#' @param dir_in          Input directory. Where input data files are loaded
#'   from.
#' @param dir_out         Output directory. Where ouput data files are stored.
#' @param samples         Numeric vector of integers. Identifies the set of
#'   samples to which peaks must be detected.
#' @param noise_level     Scalar number. The number of times the standard deviation
#'   above the noise level needed to detect a peak. IUPAC recommends `noise_level = 3` for detection.
#' @return A set of S3 objects. Additionally, a peak/roi list.
#' @family Peak detection functions
#' @export
#' @importFrom signal sgolayfilt
#' @importFrom pracma findpeaks
#' @importFrom ggplot2 geom_rect
#' @importFrom rlang .data
#' @examples
#' \donttest{
#' # Use BiocParallel library for parallelization
#' library(BiocParallel)
#' register(SnowParam(workers = 3, progressbar = TRUE, exportglobals = FALSE), default = TRUE)
#'
#' dir_in <- system.file("extdata", package = "GCIMS")
#' dir_out <- tempdir()
#' samples <- c(3,7)
#'
#' # Example of Peak Picking
#' peak_list <- gcims_rois_selection(dir_in, dir_out, samples , noise_level = 3)
#' head(peak_list)
#'
#' files <- list.files(path = dir_out, pattern = ".rds", all.files = FALSE, full.names = TRUE)
#' invisible(file.remove(files))
#' }
#'
gcims_rois_selection <- function(dir_in, dir_out, samples, noise_level){
  peak_lists <- gcims_batch_run(
    dir_in,
    dir_out,
    gcims_rois_selection_one,
    noise_level = noise_level,
    .batch_samples = samples,
    .batch_returns = function(x) {x$data$ROIs}
  )
  out <- dplyr::bind_rows(!!!peak_lists, .id = "SampleID")
  out$UniqueID <- sprintf("%s/%s", out$SampleID, out$PeakID)
  out
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

  left_bound_idx <- findZeroCrossings(x[1:max_idx] - half_max, direction = "up")
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
find_rip <- function(intensity_mat) {
  total_ion_spectrum <- rowSums(intensity_mat)
  # In most rows the RIP will be the higher peak, so it will also be the maximum in the TIS
  dt_idx_rip <- which.max(total_ion_spectrum)
  # The retention time where the RIP is maximum:
  rt_idx_apex <- which.max(intensity_mat[dt_idx_rip,])
  # The valleys of the TIS. We bound them by the TIS limits:
  dt_idx_valleys <- c(1L, pracma::findpeaks(-total_ion_spectrum)[, 2L], length(total_ion_spectrum))
  # The RIP
  dt_idx_start <- dt_idx_valleys[max(which(dt_idx_rip - dt_idx_valleys > 0))] # Find starting index of RIP
  dt_idx_end <- dt_idx_valleys[min(which(dt_idx_valleys - dt_idx_rip > 0))] # Find ending index of RIP

  # The valley detection is not perfect. Sometimes the actual peak boundary is missed
  # To workaround this, we estimate the FWHM and we clip the boundaries at the apex +- 2*FWHM
  bounds_and_widths <- find_half_max_boundaries(total_ion_spectrum)
  if (is.na(bounds_and_widths$fwhm_sym)) {
    rlang::abort(
      message = c(
        "The ROI selection could not locate the RIP width",
        "i" = "Please contact the GCIMS authors at https://github.com/sipss/GCIMS/issues and if possible submit them your sample")
    )
  }
  # Clip the peak boundaries:
  dt_idx_start <- max(dt_idx_start, dt_idx_rip - 3*bounds_and_widths$fwhm_sym)
  dt_idx_end <- min(dt_idx_end, dt_idx_rip + 3*bounds_and_widths$fwhm_sym)

  list(
    dt_idx_start = dt_idx_start,
    dt_idx_apex = dt_idx_rip,
    dt_idx_end = dt_idx_end,
    rt_idx_apex = rt_idx_apex,
    rip = intensity_mat[dt_idx_start:dt_idx_end, rt_idx_apex]
  )
}

#' Estimate the minimum distance between two peaks
#'
#' The criteria to determine this threshold is as follows:
#'
#' 1. Take the RIP
#' 2. Compute its second derivative
#' 3. Estimate the standard deviation of the
#'
#' @noRd
#' @param rip A numeric vector with the RIP shape
#'
#' @return A number, measured in indices, which can be used as the minimum peak distance when finding peaks
estimate_minpeakdistance <- function(rip) {
  filter_length <- min(21L, round(length(rip)/3))
  if (filter_length %% 2 == 0) {
    filter_length <- filter_length - 1L
  }
  rip_2nd_deriv <- signal::sgolayfilt(rip, p = 2, n = filter_length, m = 2)
  bounds_and_widths <- find_half_max_boundaries(-rip_2nd_deriv)
  fwhm <- bounds_and_widths$fwhm
  # Magic:
  fwhm/sqrt(log(2))
}

gcims_rois_selection_one <- function(x, noise_level, verbose = FALSE) {
  x$data$ROIs <- peak_detection(
    drift_time = x$data$drift_time,
    retention_time = x$data$retention_time,
    int_mat = as.matrix(x$data$data_df),
    noise_level = noise_level,
    verbose = verbose
  )
  x
}

peak_detection <- function(drift_time, retention_time, int_mat, noise_level, verbose = FALSE) {

    if (anyNA(retention_time)) {
      stop("The sample has missing values in the retention_time vector")
    }
    if (anyNA(drift_time)) {
      stop("The sample has missing values in the drift_time vector")
    }
    # 1. Data load
    spec_length <- nrow(int_mat)
    num_spec <- ncol(int_mat)
    # int_mat[3, 4] # at drift time index 3, retention time index 4

    # 2. Search of RIP position
    the_rip <- find_rip(int_mat)
    minpeakdistance <- estimate_minpeakdistance(the_rip$rip)


    # Compute the 2nd derivative for both axes
    filter1 <- signal::sgolay(p = 2, n = 21, m = 2)
    filter2 <- signal::sgolay(p = 2, n = 11, m = 2)

    drt <- t(apply(int_mat, 1, function(x) -signal::sgolayfilt(x, filter1)))
    ddt <- apply(int_mat, 2, function(x) -signal::sgolayfilt(x, filter2))

    daux <- drt + ddt

    stopifnot(dim(int_mat) == dim(daux))

    # Region without peaks: PROBAR CON ORINA
    # p1 <- hist(daux)
    # plot(p1$breaks[-1], p1$counts)
    # plot(p1$breaks[-1], log(p1$counts))
    # quantile(daux)
    # quantile(daux, c(0.25, 0.75))

    region <- find_region_without_peaks(int_mat, half_min_size = c(10, 10), noise_quantile = 0.25)
    sigmaNoise <- stats::sd(daux[region$row_min:region$row_max, region$col_min:region$col_max])

    # tt <- int_mat < quantile(int_mat, 0.15)
    # indx_noise <- which(tt == TRUE, arr.ind = TRUE)
    # sd(daux[indx_noise])


    #gaussianDistr = f.a1*exp(-((tgauss-f.b1)/f.c1).^2) + f.a2*exp(-((tgauss-f.b2)/f.c2).^2) # Fitted Gaussian

    # 5. Peaks and Zero-crossings
    ## 5.a. Retention time

    # peaksrt[[j]] <- c(i) # At the ret_time[j] seconds we find a peak in drift_time[i] ms
    peaksrt <- vector(mode = "list", length = num_spec) # Initialization of vector for peaks
    zeros_rt <- vector(mode = "list", length = num_spec) # Initialization of vector for zero crossings

    # For each IMS Spectra:
    for (j in seq_len(num_spec)) {
      # Find the max (peaks)
      # cero segunda deriv gaussiana = 1/sqrt(2)*sigma = sqrt(2)/2 sigma
      #
      #locs <- pracma::findpeaks(daux[,j], minpeakheight = noise_level*sigmaNoise, minpeakdistance = 4*sqrt(2)*sigma0)[ ,2]
      locs <- pracma::findpeaks(ddt[,j], minpeakheight = noise_level*sigmaNoise, minpeakdistance = minpeakdistance)[ ,2]

      # Find the zero-crossing points
      posrt <- findZeroCrossings(ddt[,j])
      #posrt <- findZeroCrossings(daux[,j])
      tmp <- NULL
      locs_tmp <- NULL
      for (k in seq_along(locs)) {
        dist <- locs[k] - posrt
        peakaround <- findZeroCrossings(dist)
        if (length(peakaround) >= 1){
          idx1 <- posrt[peakaround]
          idx2 <- posrt[peakaround+1]
          tmp <- cbind(tmp, rbind(idx1, idx2))
          locs_tmp <- rbind(locs_tmp, locs[k])
        }
      }
      if (!is.null(locs_tmp)) {
        zeros_rt[[j]] <- tmp
        peaksrt[[j]] <- locs_tmp
      }
    }

    ## 5.b. Peaks and Zero-Crossing for Drift Time

    peaksdt <- vector(mode = "list", length = spec_length) # Initialization of vector for peaks
    zeros_dt <- vector(mode = "list", length = spec_length) # Initialization of vector for zero crossings

    # For loop that iterates through all the columns
    for (j in seq_len(spec_length)) {
      # Find the max (peaks)
      locs <- findpeaksRois(drt[j,], MinPeakHeight = noise_level*sigmaNoise, MinPeakDistance = 1)$loc
      #locs <- pracma::findpeaks(daux[j, ], minpeakheight = noise_level*sigmaNoise)[ ,2]

      # Find the zero-crossing points
      posdt <- findZeroCrossings(drt[j, ])
      #posdt <- findZeroCrossings(daux[j, ])
      tmp <- NULL
      locs_tmp <- NULL
      for (k in seq_along(locs)){
        dist <- locs[k] - posdt
        peakaround <- findZeroCrossings(dist)
        if (length(peakaround) >= 1L) {
          idx1 <- posdt[peakaround]
          idx2 <- posdt[peakaround + 1L]
          tmp <- cbind(tmp, rbind(idx1, idx2))
          locs_tmp <- rbind(locs_tmp, locs[k])
        }
      }
      if (!is.null(locs_tmp)) {
        zeros_dt[[j]] <- tmp
        peaksdt[[j]] <- locs_tmp
      }
    }


    # Compute intersection

    peaks <- NULL
    ROIs <- NULL
    for (spec_idx in seq_along(peaksrt)) {
      dt_idx_peaks <- peaksrt[[spec_idx]]
      dt_idx_zeros <- zeros_rt[[spec_idx]]
      for (dt_idx_peak in dt_idx_peaks) {
        rt_idx_peaks <- peaksdt[[dt_idx_peak]]
        rt_idx_zeros <- zeros_dt[[dt_idx_peak]]

        if ( length(intersect(spec_idx,rt_idx_peaks)) >= 1 && (dt_idx_peak > the_rip$dt_idx_end)) {
          peaks <- rbind(peaks, c(spec_idx, dt_idx_peak))
          rt_idx_min <- rt_idx_zeros[1, rt_idx_peaks == spec_idx]
          rt_idx_max <- rt_idx_zeros[2, rt_idx_peaks == spec_idx]

          dt_idx_min <- dt_idx_zeros[1, dt_idx_peaks == dt_idx_peak]
          dt_idx_max <- dt_idx_zeros[2, dt_idx_peaks == dt_idx_peak]

          dt_idx_length <- abs(dt_idx_max[1] - dt_idx_min[1])
          rt_idx_length <- abs(rt_idx_max[1] - rt_idx_min[1])
          dt_idx_half <- floor(dt_idx_length/2)
          rt_idx_half <- floor(rt_idx_length/2)

          if (dt_idx_min[1] < dt_idx_max[1] && rt_idx_min[1] < rt_idx_max[1]) {
            dt_idx_min <- max(dt_idx_min[1] - dt_idx_half, 1L)
            rt_idx_min <- max(rt_idx_min[1] - rt_idx_half, 1L)
            dt_idx_max <- min(dt_idx_max[1] + dt_idx_half, spec_length)
            rt_idx_max <- min(rt_idx_max[1] + rt_idx_half, num_spec)
            ROIs <- rbind(ROIs, c(dt_idx_min, dt_idx_max, rt_idx_min, rt_idx_max))
          }
        }
      }
    }
    colnames(ROIs) <- c("dt_idx_min", "dt_idx_max", "rt_idx_min", "rt_idx_max")


    # Merging algorithm

    thrOverlap <- 0.2
    if (nrow(ROIs) >= 1) {
      aff <- seq_len(nrow(ROIs))

      done <- NULL
      for (j in (1:nrow(ROIs))){
        done <- c(done, j)
        R1 <- ROIs[j, ]
        for (k in c((1:nrow(ROIs))[- done])){
          R2 <- ROIs[k, ]
          if (aff[k] != j){
            if (abs(overlapPercentage(R1,R2)) > thrOverlap){
              aff[aff == k] <- j
            }
          }
        }
      }
    }


    ROIs_overlap <- NULL
    peaks_overlap <- NULL
    rtmcs <- NULL
    dtmcs <- NULL

    labels <- unique(aff)
    for (n in seq_along(labels)){
      idx <- which(aff == labels[n])
      R1 <- ROIs[idx[1], ]
      if (length(idx) > 1) {
        for (m in (2:length(idx))){
          R2 <- ROIs[idx[m], ]
          R1 <- c(dt_idx_min = min(R1["dt_idx_min"], R2["dt_idx_min"]),
                  dt_idx_max = max(R1["dt_idx_max"], R2["dt_idx_max"]),
                  rt_idx_min = min(R1["rt_idx_min"], R2["rt_idx_min"]),
                  rt_idx_max = max(R1["rt_idx_max"], R2["rt_idx_max"])
                  )
        }
      }
      ROIs_overlap <- rbind(ROIs_overlap, R1)

      patch <- int_mat[
        R1["dt_idx_min"]:R1["dt_idx_max"],
        R1["rt_idx_min"]:R1["rt_idx_max"],
        drop = FALSE
      ]
      idx_mat <- arrayInd(which.max(patch), dim(patch), dimnames(patch))
      r <- idx_mat[1, 1]
      c <- idx_mat[1, 2]
      dt_idx_apex <- R1["dt_idx_min"] + r - 1L
      rt_idx_apex <- R1["rt_idx_min"] + c - 1L
      if (rt_idx_apex > length(retention_time)) {
        rlang::abort(
          glue::glue(
            "The maximum ROI is found at the retention time index {rt_idx_apex} beyond the size of the retention time {length(retention_time)}. This should not happen"
          )
        )
      }
      if (dt_idx_apex > length(drift_time)) {
        rlang::abort(
          glue::glue(
            "The maximum ROI is found at the retention time index {dt_idx_apex} beyond the size of the retention time {length(drift_time)}. This should not happen"
          )
        )
      }

      peaks_overlap <- rbind(peaks_overlap, c(dt_idx_apex, rt_idx_apex)) # Maximo del ROI

      # roi center of mass
      v <- rowSums(patch)
      dt_cm1 <- (sum(v * seq_along(v)) / sum(v)) + R1["dt_idx_min"]  - 1L
      v <- colSums(patch)
      rt_cm1 <- (sum(v * seq_along(v)) / sum(v)) + R1["rt_idx_min"] - 1L
      rtmcs <- c(rtmcs, round(rt_cm1))
      dtmcs <- c(dtmcs, round(dt_cm1))
    }

    colnames(ROIs_overlap) <- c("dt_idx_min", "dt_idx_max", "rt_idx_min", "rt_idx_max")
    colnames(peaks_overlap) <- c("dt_idx_apex", "rt_idx_apex")
    rownames(peaks_overlap) <- NULL
    rownames(ROIs_overlap) <- NULL


    fmt <- paste0("%0", nchar(as.character(nrow(ROIs_overlap))), "d")
    peaktable <- tibble::tibble(
      PeakID = sprintf(fmt, seq_len(nrow(ROIs_overlap))),
      dt_apex_ms = NA_real_,
      rt_apex_s = NA_real_,
      dt_min_ms = NA_real_,
      dt_max_ms = NA_real_,
      rt_min_s = NA_real_,
      rt_max_s = NA_real_,
      dt_cm_ms = NA_real_,
      rt_cm_s = NA_real_,
      dt_apex_idx = peaks_overlap[,"dt_idx_apex"],
      rt_apex_idx = peaks_overlap[,"rt_idx_apex"],
      dt_min_idx = ROIs_overlap[,"dt_idx_min"],
      dt_max_idx = ROIs_overlap[,"dt_idx_max"],
      rt_min_idx = ROIs_overlap[,"rt_idx_min"],
      rt_max_idx = ROIs_overlap[,"rt_idx_max"],
      dt_cm_idx = dtmcs,
      rt_cm_idx = rtmcs
    )
    peaktable <- dplyr::mutate(
      peaktable,
      dt_apex_ms = drift_time[.data$dt_apex_idx],
      rt_apex_s = retention_time[.data$rt_apex_idx],
      dt_min_ms = drift_time[.data$dt_min_idx],
      dt_max_ms = drift_time[.data$dt_max_idx],
      rt_min_s = retention_time[.data$rt_min_idx],
      rt_max_s = retention_time[.data$rt_max_idx],
      dt_cm_ms = drift_time[.data$dt_cm_idx],
      rt_cm_s = retention_time[.data$rt_cm_idx]
    )
    peaktable
}



#---------------#
#   FUNCTIONS   #
#---------------#


#' Finds a region without peaks in the matrix
#' @noRd
#'
#' @param intens The matrix where we want to find a region without peaks
#' @param half_min_size Two integers. The region will be of size `2*half_min_size+1`.
#' @param noise_quantile A percentile, intensities below this percentile are considered background
#'
#' @return A list with the indexes that bound the region without peaks (`row_min`, `row_max`, `col_min`, `col_max`)
#'
#' @examples
#' x <- matrix(0, nrow = 500, ncol = 500)
#' x[1:200, 1:300] <- 500 # you could call that a weird peak...
#' region <- find_region_without_peaks(x)
#' # The region has zero values
#' sum(x[region$row_min:region$row_max, region$col_min:region$col_max]) == 0
find_region_without_peaks <- function(intens, half_min_size = c(10, 10), noise_quantile = 0.25) {
  thr <- stats::quantile(intens, noise_quantile)
  binarized <- intens <= thr

  center_i <- seq.int(from = 1L+half_min_size[1], to = dim(binarized)[1] - half_min_size[1], by = 2L*half_min_size[1]+1L)
  center_j <- seq.int(from = 1L+half_min_size[2], to = dim(binarized)[2] - half_min_size[2], by = 2L*half_min_size[2]+1L)


  for (i in center_i) {
    for (j in center_j) {
      row_min <- i - half_min_size[1]
      row_max <- i + half_min_size[1]
      col_min <- j - half_min_size[2]
      col_max <- j + half_min_size[2]
      if (all(binarized[row_min:row_max, col_min:col_max])) {
        return(list(row_min=row_min, row_max = row_max, col_min = col_min, col_max = col_max))
      }
    }
  }
  stop("Could not find a region without peaks")
}


#----------------------#
#   findZeroCrossings  #
#----------------------#

findZeroCrossings <- function(x, direction = c("both", "up", "down")) {
  direction <- match.arg(direction)
  signs <- sign(x)
  dsigns <- diff(signs)
  if (direction == "up" || direction == "both") {
    pos_plus <- which(dsigns == 2)
  } else {
    pos_plus <- integer(0L)
  }
  if (direction == "down" || direction == "both") {
    pos_minus <- which(dsigns == -2)
  } else {
    pos_minus <- integer(0L)
  }
  pos <- sort(union(pos_plus + 1L, pos_minus))
  return(pos)
}

#----------------------#
#   overlapPercentage  #
#----------------------#


overlapPercentage <- function(ROI1, ROI2){
  if (!(ROI1[1] >= ROI2[2] | ROI1[2] <= ROI2[1] | ROI1[3] >= ROI2[4] | ROI1[4] <= ROI2[3])){
    area1 <- (ROI1[2] - ROI1[1])*(ROI1[4] - ROI1[3])
    area2 <- (ROI2[2] - ROI2[1])*(ROI2[4] - ROI2[3])
    x_left <- max(ROI1[1], ROI2[1])
    y_top <- min(ROI1[4], ROI2[4])
    x_right <- min(ROI1[2], ROI2[2])
    y_bottom <- max(ROI1[3], ROI2[3])
    overlapping_area <- (x_right - x_left)*(y_top - y_bottom)
    if (area1 == overlapping_area | area2 == overlapping_area){
      p <- 1
    } else {
      p <- overlapping_area / (area1 + area2 - overlapping_area)
    }
  } else {
    p <- 0
  }
  return(p)
}



#------------------#
#   findpeaksRois  #
#------------------#

findpeaksRois <- function(data, MinPeakDistance = 1, MinPeakHeight = 1) {
  ld <- length(data)
  tmp <- data

  df1 <- diff(data, differences = 1)[c(1, 1:(ld - 1))]
  df2 <- diff(data, differences = 2)[c(1, 1, 1:(ld - 2))]
  idx <- which(df1 * c(df1[2:length(df1)], 0) <= 0 & c(df2[2:length(df2)], 0) < 0)
  max_ind <- which(data[idx] > MinPeakHeight)
  idx <- idx[max_ind]
  if (length(idx) > 0) {
    tmp <- sort(data[idx], decreasing = TRUE, index = TRUE)
    idx_s <- idx[tmp$ix]
    D <- with(expand.grid(A = idx_s, B = t(idx_s)), abs(A - B))
    dim(D) <- c(length(idx_s), length(idx_s))
    if (dim(D)[1] > 1){
      diag(D) <- NA
    }
    if (any(D) < MinPeakDistance) {
      i <- 1
      pointsmax <- seq_along(idx_s)
      checked <- NULL
      idx_pruned <- idx_s
      while (length(pointsmax) > 0) {
        d <- D[pointsmax[1], ]
        checked <- c(checked, pointsmax[1])
        pointsmax <- pointsmax[-1]
        remains <- setdiff(which(d < MinPeakDistance), checked)
        if (length(remains) > 0) {
          idx_pruned <- setdiff(idx_pruned, idx_s[remains])
          checked <- c(checked, remains)
          pointsmax <- setdiff(pointsmax, checked)
        }
      }
      idx <- idx_pruned
    }

    idx <- sort(idx)
    idx.pruned <- idx
    n <- length(idx)
    for (i in 1:n) {
      ind <- round(max(idx[i] - MinPeakDistance/2, 1)):round(min(idx[i] + MinPeakDistance/2, ld))
      pp <- rep(0L, 3)
      H <- data[idx[i]]
      xm <- idx[i]
      pp <- rep(1L, 3)
      pp[1] <- pracma::mldivide((ind - xm)^2, (data[ind] - H))
      pp[2] <- -2 * pp[1] * xm
      pp[3] <- H + pp[1] * xm^2
      sigma <- sqrt(abs(1/(2*pp[3])))
      if (abs(pp[1]) > 0){
        width <- exp(abs(pp[1]) + abs(pp[2])*sigma^2^2/(2*sigma^2))
      } else {
        width <- exp(1 + 1*sigma^2^2/(2*sigma^2))
      }
      MinPeakWidth <- (H/2 - 15)
      MaxPeakWidth <- (H/2 + 20)
      if ((width > MaxPeakWidth || width < MinPeakWidth) ||
          pp[1] > 0 || H < MinPeakHeight || data[idx[i]] < H ||
          abs(idx[i] - xm) > MinPeakDistance/2) {
        idx.pruned <- setdiff(idx.pruned, idx[i])
      }
    }

    idx <- idx.pruned
    pks <- data[idx]
    list(pks = pks, loc = idx)
  } else {
    list(pks = NULL, loc = NULL)
  }
}


