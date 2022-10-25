merge_overlapping_rois_1D <- function(ROIs, y, iou_overlap_threshold) {
  if (iou_overlap_threshold > 1) {
    return(ROIs)
  }

  while (TRUE) {
    ROIs$group_to_merge <- seq_len(nrow(ROIs))

    ious <- tibble::as_tibble(
      expand.grid(
        ROIi = seq_len(nrow(ROIs)),
        ROIj = seq_len(nrow(ROIs))
      )
    )
    ious <- dplyr::filter(ious, .data$ROIi < .data$ROIj)
    ious$iou <- purrr::map2_dbl(ious$ROIi, ious$ROIj, function(i, j) {
      intersectionOverUnion1D(ROIs[i,], ROIs[j,])
    })

    ious <- dplyr::filter(ious, .data$iou > iou_overlap_threshold)
    if (nrow(ious) == 0) {
      break
    }
    ious <- dplyr::arrange(ious, dplyr::desc(.data$iou))
    for (i in seq_len(nrow(ious))) {
      group_to_keep <- ROIs$group_to_merge[ious$ROIi[i]]
      group_to_remove <- ROIs$group_to_merge[ious$ROIj[i]]
      ROIs$group_to_merge[ROIs$group_to_merge == group_to_remove] <- group_to_keep
    }
    ROIs <- dplyr::group_by(ROIs, .data$group_to_merge)
    ROIs <- dplyr::summarise(
      ROIs,
      idx_apex = .data$idx_apex[which.max(.data$int_apex_au)],
      idx_min = min(.data$idx_min),
      idx_max = max(.data$idx_max),
      int_apex_au = max(.data$int_apex_au),
      .groups = "drop"
    )
    ROIs <- dplyr::ungroup(ROIs)
  }
  ROIs <- ROIs[, c("idx_apex", "idx_min", "idx_max", "int_apex_au"), drop = FALSE]
  return(ROIs)
}

compute_center_of_mass_1D <- function(ROIs, y) {
  ROIs$cm_idx <- NA_integer_
  for (i in seq_len(nrow(ROIs))) {
    R1 <- ROIs[i, ]
    patch <- y[R1$idx_min:R1$idx_max]
    # roi center of mass
    ROIs$cm_idx[i] <- round((sum(patch * seq_along(patch)) / sum(patch)) + R1$idx_min  - 1L)
  }
  ROIs
}

rois_to_peaklist_1D <- function(ROIs, x, y, deriv2) {
  fmt <- paste0("%0", nchar(as.character(nrow(ROIs))), "d")
  peak_list <- tibble::tibble(
    PeakID = sprintf(fmt, seq_len(nrow(ROIs))),
    apex = NA_real_,
    int_apex_au = y[ROIs$idx_apex],
    deriv_apex_au = deriv2[ROIs$idx_apex],
    min = NA_real_,
    max = NA_real_,
    cm = NA_real_,
    apex_idx = ROIs$idx_apex,
    min_idx = ROIs$idx_min,
    max_idx = ROIs$idx_max,
    cm_idx = ROIs$cm_idx
  )
  peak_list <- dplyr::mutate(
    peak_list,
    apex = .env$x[.data$apex_idx],
    min = .env$x[.data$min_idx],
    max = .env$x[.data$max_idx],
    cm = .env$x[.data$cm_idx]
  )
  peak_list
}

#' Peak and ROI detection (1D)
#'
#' Detects regions of interest with peaks in a signal
#'
#' Peaks are detected on the intensity vector.
#'
#' In detail, the approach is as follows:
#'
#' 1. We compute the second derivative with respect to the drift and retention times.
#' 2. Based on the given peak width ranges, mexican hat wavelets are scaled
#' 3. Peaks on are detected using [MassSpecWavelet::peakDetectionCWT()].
#' 4. We merge similar ROIs using a threshold on the 1D-intersection over union
#' 5. Get some ROI metrics and return.
#'
#' For the MassSpecWavelet-based peak detection, the `scales` are computed based on
#' the requested peak widths. Besides, the scales, further tuning beyond the
#' MassSpecWavelet defaults is possible through the `peakDetectionCWTParams` argument.
#' By default, the only change we introduce is the `exclude0scaleAmpThresh = TRUE`
#' which is a reasonable peak detection setting not enabled in MassSpecWavelet
#' for backwards compatibility reasons.
#'
#' @keywords internal
#' @param x The x axis, to determine units
#' @param y The vector where to find peaks, of dimensions `length(x)`
#' @param verbose If `TRUE` information will be printed on screen
#' @param length_in_xunits Length of the filter used to compute the second derivative. See details.
#' @param peakwidth_range_xunits A vector of length 2 with the minimum and maximum peak width. See details.
#' @param peakDetectionCWTParams Additional parameters to [MassSpecWavelet::peakDetectionCWT()]. See details.
#' @param extension_factor A number to extend the ROIs beyond their default size
#' @param iou_overlap_threshold A number, between 0 and 1. Pairs of ROIs with an intersection over union larger than this threshold are merged.
#' @param debug If TRUE, return as well the debug information
#' @return A list with the `peak_list` and `debug_info` elements.
findPeaksImpl1D <- function(
    x, y,
    verbose = FALSE,
    length_in_xunits = 0.07,
    peakwidth_range_xunits = c(0.15, 0.4),
    peakDetectionCWTParams = list(exclude0scaleAmpThresh = TRUE),
    extension_factor = 0,
    iou_overlap_threshold = 0.2,
    debug = FALSE
) {
  x_step <- x[2L] - x[1L]
  x_length_pts <- units_to_points(length_in_xunits, x_step, must_odd = TRUE)
  spec_length <- length(y)
  # int_mat[3, 4] # at drift time index 3, retention time index 4

  # Compute the 2nd derivative for both axes to get the inflection points
  sg_filt <- signal::sgolay(p = 2L, n = x_length_pts, m = 2)
  deriv2 <- sgolayfilt(y, sg_filt)
  # 5. Peaks and Zero-crossings
  scales <- prep_wav_from_peakwidths(
    x,
    peakwidth_min = peakwidth_range_xunits[1L],
    peakwidth_max = peakwidth_range_xunits[2L]
  )

  if (verbose) {
    rlang::inform(
      message = c(
        "i" = paste0("Using the following scales: ", glue::glue_collapse(scales$scales, sep = ", ")),
      )
    )
  }

  peaks_and_zeros_extra <- detect_peaks_and_zeros(
    xmat = matrix(y, ncol = 1),
    rowwise = FALSE,
    scales = scales,
    peakDetectionCWTParams = peakDetectionCWTParams,
    signals_to_save_extra = if (debug) 1L else NULL,
    xmat_zeros = matrix(deriv2, ncol = 1)
  )
  peaks_and_zeros <- peaks_and_zeros_extra$peaks_and_zeros
  colnames(peaks_and_zeros) <- c("nothing", "idx_apex", "idx_min", "idx_max")
  peaks_and_zeros <- peaks_and_zeros[, c("idx_apex", "idx_min", "idx_max")]
  if (debug) {
    debug_info <- peaks_and_zeros_extra$extra
  } else {
    debug_info <- NULL
  }

  ROIs <- peaks_and_zeros

  ROIs <- dplyr::mutate(
    ROIs,
    length_pts = .data$idx_max - .data$idx_min,
    idx_min = pmax(1L,          floor(.data$idx_min - extension_factor * .data$length_pts/2)),
    idx_max = pmin(spec_length, ceiling(.data$idx_max + extension_factor * .data$length_pts/2)),
  )
  ROIs <- ROIs[,c("idx_apex", "idx_min", "idx_max"), drop = FALSE]

  # Keep only ROIs of non-zero area:
  ROIs <- dplyr::filter(
    ROIs,
    .data$idx_min < .data$idx_max
  )

  # Compute intensity:
  ROIs$int_apex_au <- y[ROIs$idx_apex]

  # Merging algorithm:
  ROIs_overlap <- merge_overlapping_rois_1D(ROIs, y, iou_overlap_threshold)

  # Compute ROI's center of mass:
  ROIs_overlap <- compute_center_of_mass_1D(ROIs = ROIs_overlap, y = y)

  # Convert into a data frame / peak list:
  peak_list <- rois_to_peaklist_1D(
    ROIs = ROIs_overlap,
    x = x,
    y = y,
    deriv2 = deriv2
  )
  list(
    peak_list = peak_list,
    debug_info = debug_info
  )
}


intersectionOverUnion1D <- function(ROI1, ROI2){
  if (
    ROI1$idx_min >= ROI2$idx_max || # ROI1 starts after ROI2 ends
    ROI1$idx_max <= ROI2$idx_min    # ROI1 ends before ROI2 starts
  ) {
    # There is no overlap between the ROIs, the intersection is zero
    return(0)
  }

  area1 <- ROI1$idx_max - ROI1$idx_min
  stopifnot(area1 >= 0) # area can't be negative
  area2 <- ROI2$idx_max - ROI2$idx_min
  stopifnot(area2 >= 0) # area can't be negative
  if (area1 == 0 || area2 == 0) {
    # The intersection of a zero area with whatever will be zero.
    # If both areas are zero, we also return zero (iou undetermined)
    return(0)
  }
  # intersection:
  ROI_int <- list(
    idx_min = max(ROI1$idx_min, ROI2$idx_min),
    idx_max = min(ROI1$idx_max, ROI2$idx_max)
  )
  area_int <- ROI_int$idx_max - ROI_int$idx_min
  stopifnot(area_int >= 0) # area can't be negative
  area_union <- area1 + area2 - area_int
  stopifnot(area_union > 0) # should be positive here
  return(area_int/area_union)
}
