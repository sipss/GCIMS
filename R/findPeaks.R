compute_second_deriv <- function(int_mat, dt_length_pts, rt_length_pts, dt_order = 2, rt_order = 2) {
  filter1 <- signal::sgolay(p = rt_order, n = rt_length_pts, m = 2)
  filter2 <- signal::sgolay(p = dt_order, n = dt_length_pts, m = 2)
  drt <- sgolayfilt(int_mat, filter1, rowwise = TRUE)
  ddt <- sgolayfilt(int_mat, filter2)
  list(
    ddt = ddt,
    drt = drt
  )
}

interleave_peaks_with_zeros <- function(peak_idx, zero_idx, signal_idx) {
  if (length(peak_idx) == 0 || length(zero_idx) == 0) {
    out <- tibble::tibble(
      signal_idx = integer(0L),
      idx_apex = integer(0L),
      idx_min = integer(0L),
      idx_max = integer(0L)
    )
    return(out)
  }
  out <- tibble::tibble(
    signal_idx = rep(signal_idx, length(peak_idx)),
    idx_apex = peak_idx,
    idx_min = NA_integer_,
    idx_max = NA_integer_
  )
  zero_idx2 <- 1L
  for (i in seq_len(nrow(out))) {
    if (out$idx_apex[i] <= zero_idx[1L]) {
      next
    }
    while (zero_idx2 <= length(zero_idx) && out$idx_apex[i] >= zero_idx[zero_idx2]) {
      zero_idx2 <- zero_idx2 + 1L
    }
    if (zero_idx2 > length(zero_idx)) {
      break
    }
    out$idx_min[i] <- zero_idx[zero_idx2 - 1L]
    out$idx_max[i] <- zero_idx[zero_idx2]
  }
  out <- dplyr::filter(out, !is.na(.data$idx_min))
  out
}

detect_peaks_and_zeros_one_signal <- function(x, x_zeros, idx, save_debugging, prep_wav, peakDetectionCWTParams) {
  peakdet <- rlang::exec(
    MassSpecWavelet::peakDetectionCWT,
    ms = x,
    scales = prep_wav,
    !!!peakDetectionCWTParams
  )
  if (is.null(x_zeros)) {
    x_zeros <- x
  }
  peak_idx <- unname(peakdet$majorPeakInfo$peakIndex)
  zero_idx <- GCIMS:::findZeroCrossings(x_zeros)
  out <- interleave_peaks_with_zeros(peak_idx, zero_idx, signal_idx = idx)
  #stopifnot(ncol(out) == 4L)
  if (save_debugging) {
    extra <- list(
      peakdet = peakdet,
      peak_idx = peak_idx,
      zero_idx = zero_idx,
      interleaved = out
    )
  } else {
    extra <- NULL
  }
  list(
    peaks_and_zeros = out,
    extra = extra
  )
}

detect_peaks_and_zeros <- function(xmat, rowwise, scales, peakDetectionCWTParams = list(exclude0scaleAmpThresh = TRUE),
                                   signals_to_save_extra = NULL, xmat_zeros = NULL) {
  if (rowwise) {
    signal_length <- ncol(xmat)
    num_signals <- nrow(xmat)
    signals <- vector("list", num_signals)
    for (i in seq_len(num_signals)) {
      signals[[i]] <- as.numeric(xmat[i, ])
    }
    signals_zeros <- vector("list", num_signals)
    if (!is.null(xmat_zeros)) {
      for (i in seq_len(num_signals)) {
        signals_zeros[[i]] <- as.numeric(xmat_zeros[i, ])
      }
    }
  } else {
    signal_length <- nrow(xmat)
    num_signals <- ncol(xmat)
    signals <- vector("list", num_signals)
    for (i in seq_len(num_signals)) {
      signals[[i]] <- as.numeric(xmat[, i])
    }
    signals_zeros <- vector("list", num_signals)
    if (!is.null(xmat_zeros)) {
      for (i in seq_len(num_signals)) {
        signals_zeros[[i]] <- as.numeric(xmat_zeros[, i])
      }
    }

  }
  prep_wav <- MassSpecWavelet::prepareWavelets(signal_length, scales = scales)
  save_debugging <- rep(FALSE, num_signals)
  if (!is.null(signals_to_save_extra)) {
    save_debugging[signals_to_save_extra] <- TRUE
  }
  peaks_and_zeros_and_extra <- mymapply(
    FUN = detect_peaks_and_zeros_one_signal,
    x = signals,
    x_zeros = signals_zeros,
    idx = seq_len(num_signals),
    save_debugging = save_debugging,
    MoreArgs = list(
      prep_wav = prep_wav,
      peakDetectionCWTParams = peakDetectionCWTParams
    ),
    SIMPLIFY = FALSE
  )
  peaks_and_zeros <- purrr::map_dfr(peaks_and_zeros_and_extra, "peaks_and_zeros")
  stopifnot(ncol(peaks_and_zeros) == 4L)
  if (!is.null(signals_to_save_extra)) {
    extra <- purrr::map(peaks_and_zeros_and_extra, "extra")
    extra <- extra[signals_to_save_extra]
    names(extra) <- as.character(signals_to_save_extra)
  } else {
    extra <- NULL
  }
  list(
    peaks_and_zeros = peaks_and_zeros,
    extra = extra
  )
}

intersect_peaks_and_zeros_in_both_axes <- function(peaks_and_zeros_drt, peaks_and_zeros_ddt) {
  ROIs <- dplyr::inner_join(
    peaks_and_zeros_drt,
    peaks_and_zeros_ddt,
    by = c(
      "dt_idx" = "dt_idx_apex",
      "rt_idx_apex" = "rt_idx"
    )
  )
  ROIs <- dplyr::select(
    ROIs,
    "dt_idx_apex" = "dt_idx",
    "rt_idx_apex" = "rt_idx_apex",
    "dt_idx_min" = "dt_idx_min",
    "dt_idx_max" = "dt_idx_max",
    "rt_idx_min" = "rt_idx_min",
    "rt_idx_max" = "rt_idx_max"
  )
  ROIs
}


merge_overlapping_rois <- function(ROIs, int_mat, iou_overlap_threshold) {
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
      intersectionOverUnion(ROIs[i,], ROIs[j,])
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
      dt_idx_apex = .data$dt_idx_apex[which.max(.data$int_apex_au)],
      rt_idx_apex = .data$rt_idx_apex[which.max(.data$int_apex_au)],
      dt_idx_min = min(.data$dt_idx_min),
      dt_idx_max = max(.data$dt_idx_max),
      rt_idx_min = min(.data$rt_idx_min),
      rt_idx_max = max(.data$rt_idx_max),
      int_apex_au = max(.data$int_apex_au),
      .groups = "drop"
    )
    ROIs <- dplyr::ungroup(ROIs)
  }
  ROIs <- ROIs[, c("dt_idx_apex", "rt_idx_apex", "dt_idx_min", "dt_idx_max", "rt_idx_min",
                   "rt_idx_max", "int_apex_au"), drop = FALSE]
  return(ROIs)
}

compute_center_of_mass <- function(ROIs, int_mat) {
  ROIs$dt_cm_idx <- NA_integer_
  ROIs$rt_cm_idx <- NA_integer_
  for (i in seq_len(nrow(ROIs))) {
    R1 <- ROIs[i, ]
    patch <- int_mat[
      R1$dt_idx_min:R1$dt_idx_max,
      R1$rt_idx_min:R1$rt_idx_max,
      drop = FALSE
    ]

    # roi center of mass
    dt_v <- rowSums(patch)
    rt_v <- colSums(patch)
    ROIs$dt_cm_idx[i] <- round((sum(dt_v * seq_along(dt_v)) / sum(dt_v)) + R1$dt_idx_min  - 1L)
    ROIs$rt_cm_idx[i] <- round((sum(rt_v * seq_along(rt_v)) / sum(rt_v)) + R1$rt_idx_min - 1L)
  }
  ROIs
}

rois_to_peaklist <- function(ROIs, drift_time, retention_time, int_mat, ddt, drt) {
  fmt <- paste0("%0", nchar(as.character(nrow(ROIs))), "d")
  peak_list <- tibble::tibble(
    PeakID = sprintf(fmt, seq_len(nrow(ROIs))),
    dt_apex_ms = NA_real_,
    rt_apex_s = NA_real_,
    int_apex_au = apply(ROIs, 1L, function(ROI) int_mat[ROI["dt_idx_apex"], ROI["rt_idx_apex"]]),
    ddt_apex_au = apply(ROIs, 1L, function(ROI) ddt[ROI["dt_idx_apex"], ROI["rt_idx_apex"]]),
    drt_apex_au = apply(ROIs, 1L, function(ROI) drt[ROI["dt_idx_apex"], ROI["rt_idx_apex"]]),
    dt_min_ms = NA_real_,
    dt_max_ms = NA_real_,
    rt_min_s = NA_real_,
    rt_max_s = NA_real_,
    dt_cm_ms = NA_real_,
    rt_cm_s = NA_real_,
    dt_apex_idx = ROIs$dt_idx_apex,
    rt_apex_idx = ROIs$rt_idx_apex,
    dt_min_idx = ROIs$dt_idx_min,
    dt_max_idx = ROIs$dt_idx_max,
    rt_min_idx = ROIs$rt_idx_min,
    rt_max_idx = ROIs$rt_idx_max,
    dt_cm_idx = ROIs$dt_cm_idx,
    rt_cm_idx = ROIs$rt_cm_idx
  )
  peak_list <- dplyr::mutate(
    peak_list,
    dt_apex_ms = .env$drift_time[.data$dt_apex_idx],
    rt_apex_s = .env$retention_time[.data$rt_apex_idx],
    dt_min_ms = .env$drift_time[.data$dt_min_idx],
    dt_max_ms = .env$drift_time[.data$dt_max_idx],
    rt_min_s = .env$retention_time[.data$rt_min_idx],
    rt_max_s = .env$retention_time[.data$rt_max_idx],
    dt_cm_ms = .env$drift_time[.data$dt_cm_idx],
    rt_cm_s = .env$retention_time[.data$rt_cm_idx]
  )
  peak_list
}

prep_wav_from_peakwidths <- function(axis, peakwidth_min, peakwidth_max) {
  step <- axis[2L] - axis[1L]
  peakwidth_min_pts <- units_to_points(peakwidth_min, step)
  peakwidth_max_pts <- units_to_points(peakwidth_max, step)
  if (peakwidth_min_pts < 1) {
    peakwidth_min_pts <- 1
  }
  # TODO: maybe warn if this happens
  if (peakwidth_max_pts > length(axis)) {
    peakwidth_max_pts <- length(axis)
  }
  base <- 1.5
  scales <- unique(
    c(
      1,
      peakwidth_min_pts,
      base^seq(from = log(peakwidth_min_pts, base), to = log(peakwidth_max_pts, base), by = 1),
      peakwidth_max_pts
    )
  )
  if (length(scales) > 12) {
    base <- 1.5
    scales <- unique(
      c(
        1,
        peakwidth_min_pts,
        base^seq(from = log(peakwidth_min_pts, base), to = log(peakwidth_max_pts, base), length.out = 12),
        peakwidth_max_pts
      )
    )
  }
  scales <- unique(round(scales))

  MassSpecWavelet::prepareWavelets(
    length(axis),
    scales = scales,
    wavelet = "mexh",
    wavelet_xlimit = 8,
    wavelet_length = 1024L,
    extendLengthScales = TRUE
  )
}

#' Peak and ROI detection
#'
#' Detects regions of interest with peaks in a sample.
#'
#' Peaks are detected on the partial second derivative of the intensity matrix.
#'
#' In detail, the approach is as follows:
#'
#' 1. We compute the second derivative with respect to the drift and retention times.
#' 2. Based on the given peak width ranges, mexican hat wavelets are scaled
#' 3. Peaks on the second derivatives are detected using [MassSpecWavelet::peakDetectionCWT()].
#' 4. We intersect the maxima to discard saddle points
#' 5. We merge similar ROIs using a threshold on the intersection over union
#' 6. Get some ROI metrics and return.
#'
#' For the MassSpecWavelet-based peak detection, the `scales` are computed based on
#' the requested peak widths. Besides, the scales, further tuning beyond the
#' MassSpecWavelet defaults is possible through the `*_peakDetectionCWTParams` argument.
#' By default, the only change we introduce is the `exclude0scaleAmpThresh = TRUE`
#' which is a reasonable peak detection setting not enabled in MassSpecWavelet
#' for backwards compatibility reasons.
#'
#' @keywords internal
#' @param drift_time The drift time vector
#' @param retention_time The retention time vector
#' @param int_mat The intensity matrix, of dimensions `length(drift_time)` rows and `length(retention_time)` columns
#' @param verbose If `TRUE` information will be printed on screen
#' @param dt_length_ms,rt_length_s Length of the filters used to compute the
#' second derivative. See details. FIXME (experimental): You can set both
#' `dt_length_ms` and `rt_length_s` to `-1` if you want to use the raw signal
#' for the peak detection instead of the second derivatives.
#' @param dt_peakwidth_range_ms,rt_peakwidth_range_s A vector of length 2 with the minimum and maximum peak width. See details
#' @param dt_peakDetectionCWTParams,rt_peakDetectinoCWTParams Additional parameters to [MassSpecWavelet::peakDetectionCWT()]. See details
#' @param exclude_rip Whether to exclude ROIs with a drift time apex smaller than the RIP drift time end.
#' @param iou_overlap_threshold A number, between 0 and 1. Pairs of ROIs with an intersection over union larger than this threshold are merged.
#' @param debug_idx A list with two numeric vectors named `dt` and `rt` each of them having a the indices to where debug info is kept
#' @return A list with the `peak_list` and `debug_info` elements.
findPeaksImpl <- function(
    drift_time, retention_time, int_mat,
    verbose = FALSE,
    dt_length_ms = 0.07,
    rt_length_s = 8,
    dt_peakwidth_range_ms = c(0.15, 0.4),
    rt_peakwidth_range_s = c(10, 25),
    dt_peakDetectionCWTParams = list(exclude0scaleAmpThresh = TRUE),
    rt_peakDetectionCWTParams = list(exclude0scaleAmpThresh = TRUE),
    dt_extension_factor = 0,
    rt_extension_factor = 0,
    exclude_rip = FALSE,
    iou_overlap_threshold = 0.2,
    debug_idx = list(dt = NULL, rt = NULL)
) {

    debug_info <- list()
    if (!is.null(debug_idx$rt) && length(debug_idx$rt) > 0) {
      debug_info$rt <- list()
    }
    if (!is.null(debug_idx$dt) && length(debug_idx$dt) > 0) {
      debug_info$dt <- list()
    }
    if (length(debug_info) == 0L) {
      debug_info <- NULL
    }

    dt_step_ms <- drift_time[2L] - drift_time[1L]
    dt_length_pts <- units_to_points(dt_length_ms, dt_step_ms, must_odd = TRUE)
    rt_step_s <- retention_time[2L] - retention_time[1L]
    rt_length_pts <- units_to_points(rt_length_s, rt_step_s, must_odd = TRUE)
    spec_length <- nrow(int_mat)
    num_spec <- ncol(int_mat)
    # int_mat[3, 4] # at drift time index 3, retention time index 4

    # Given how MassSpecWavelet works, maybe we could just do it on the original signal!
    # Compute the 2nd derivative for both axes
    deriv2 <- compute_second_deriv(
      int_mat,
      dt_length_pts = abs(dt_length_pts),
      rt_length_pts = abs(rt_length_pts),
      dt_order = 2L,
      rt_order = 2L
    )
    drt <- -deriv2$drt
    ddt <- -deriv2$ddt
    rm(deriv2)

    # 5. Peaks and Zero-crossings
    scales_rt <- prep_wav_from_peakwidths(
      retention_time,
      peakwidth_min = rt_peakwidth_range_s[1L],
      peakwidth_max = rt_peakwidth_range_s[2L]
    )
    scales_dt <- prep_wav_from_peakwidths(
      drift_time,
      peakwidth_min = dt_peakwidth_range_ms[1L],
      peakwidth_max = dt_peakwidth_range_ms[2L]
    )

    if (verbose) {
      rlang::inform(
        message = c(
          "Using the following scales",
          "i" = paste0("Retention time scales: ", glue::glue_collapse(scales_rt$scales, sep = ", ")),
          "i" = paste0("Drift time scales: ", glue::glue_collapse(scales_dt$scales, sep = ", "))
        )
      )
    }

    if (rt_length_pts > 0) {
      xmat <- drt
      xzeros <- NULL # same as xmat by default
    } else {
      xmat <- int_mat
      xzeros <- drt
    }
    peaks_and_zeros_drt_extra <- detect_peaks_and_zeros(
      xmat = xmat,
      rowwise = TRUE,
      scales = scales_rt,
      peakDetectionCWTParams = rt_peakDetectionCWTParams,
      signals_to_save_extra = debug_idx$rt,
      xmat_zeros = xzeros
    )
    peaks_and_zeros_drt <- peaks_and_zeros_drt_extra$peaks_and_zeros
    colnames(peaks_and_zeros_drt) <- c("dt_idx", "rt_idx_apex", "rt_idx_min", "rt_idx_max")
    if (!is.null(debug_idx$rt)) {
      debug_info$rt <- peaks_and_zeros_drt_extra$extra
    }

    if (rt_length_pts > 0) {
      xmat <- ddt
      xzeros <- NULL # same as xmat by default
    } else {
      xmat <- int_mat
      xzeros <- ddt
    }

    peaks_and_zeros_ddt_extra <- detect_peaks_and_zeros(
      xmat = xmat,
      rowwise = FALSE,
      scales = scales_dt,
      peakDetectionCWTParams = dt_peakDetectionCWTParams,
      signals_to_save_extra = debug_idx$dt,
      xmat_zeros = xzeros
    )
    peaks_and_zeros_ddt <- peaks_and_zeros_ddt_extra$peaks_and_zeros
    colnames(peaks_and_zeros_ddt) <- c("rt_idx", "dt_idx_apex", "dt_idx_min", "dt_idx_max")
    if (!is.null(debug_idx$dt)) {
      debug_info$dt <- peaks_and_zeros_ddt_extra$extra
    }

    ROIs <- intersect_peaks_and_zeros_in_both_axes(peaks_and_zeros_drt, peaks_and_zeros_ddt)

    ROIs <- dplyr::mutate(
      ROIs,
      dt_length_pts = .data$dt_idx_max - .data$dt_idx_min,
      rt_length_pts = .data$rt_idx_max - .data$rt_idx_min,
      dt_idx_min = pmax(1L,          floor(.data$dt_idx_min - dt_extension_factor * .data$dt_length_pts/2)),
      dt_idx_max = pmin(spec_length, ceiling(.data$dt_idx_max + dt_extension_factor * .data$dt_length_pts/2)),
      rt_idx_min = pmax(1L,          floor(.data$rt_idx_min - rt_extension_factor * .data$rt_length_pts/2)),
      rt_idx_max = pmin(num_spec,    ceiling(.data$rt_idx_max + rt_extension_factor * .data$rt_length_pts/2))
    )
    ROIs <- ROIs[,c("dt_idx_apex", "rt_idx_apex",
                    "dt_idx_min", "dt_idx_max",
                    "rt_idx_min", "rt_idx_max"), drop = FALSE]

    # Keep only ROIs of non-zero area:
    ROIs <- dplyr::filter(
      ROIs,
      .data$dt_idx_min < .data$dt_idx_max,
      .data$rt_idx_min < .data$rt_idx_max
    )

    # Compute intensity:
    ROIs$int_apex_au = purrr::map2_dbl(
      ROIs$dt_idx_apex,
      ROIs$rt_idx_apex,
      function(dt_idx, rt_idx, int_mat) {
        int_mat[dt_idx, rt_idx]
      },
      int_mat = int_mat
    )

    if (isTRUE(exclude_rip)) {
      the_rip <- find_rip(int_mat, verbose = verbose, retention_time = retention_time, drift_time = drift_time)
      ROIs <- dplyr::filter(ROIs, .data$dt_idx_apex >= !!the_rip$dt_idx_end)
    }

    # Merging algorithm:
    ROIs_overlap <- merge_overlapping_rois(ROIs, int_mat, iou_overlap_threshold)

    # Compute ROI's center of mass:
    ROIs_overlap <- compute_center_of_mass(ROIs = ROIs_overlap, int_mat = int_mat)

    # Convert into a data frame / peak list:
    peak_list <- rois_to_peaklist(
      ROIs = ROIs_overlap,
      drift_time = drift_time,
      retention_time = retention_time,
      int_mat = int_mat,
      ddt = ddt,
      drt = drt
    )
    list(
      peak_list = peak_list,
      debug_info = debug_info
    )
}


#---------------#
#   FUNCTIONS   #
#---------------#


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

#--------------------------#
#   intersectionOverUnion  #
#--------------------------#

intersectionOverUnion <- function(ROI1, ROI2){
  if (
        ROI1$dt_idx_min >= ROI2$dt_idx_max || # ROI1 starts after ROI2 ends
        ROI1$dt_idx_max <= ROI2$dt_idx_min || # ROI1 ends before ROI2 starts
        ROI1$rt_idx_min >= ROI2$rt_idx_max || # ROI1 starts after ROI2 ends
        ROI1$rt_idx_max <= ROI2$rt_idx_min    # ROI1 ends before ROI2 starts
      ) {
    # There is no overlap between the ROIs, the intersection is zero
    return(0)
  }

  area1 <- (ROI1$dt_idx_max - ROI1$dt_idx_min)*(ROI1$rt_idx_max - ROI1$rt_idx_min)
  stopifnot(area1 >= 0) # area can't be negative
  area2 <- (ROI2$dt_idx_max - ROI2$dt_idx_min)*(ROI2$rt_idx_max - ROI2$rt_idx_min)
  stopifnot(area2 >= 0) # area can't be negative
  if (area1 == 0 || area2 == 0) {
    # The intersection of a zero area with whatever will be zero.
    # If both areas are zero, we also return zero (iou undetermined)
    return(0)
  }
  # intersection:
  ROI_int <- list(
    dt_idx_min = max(ROI1$dt_idx_min, ROI2$dt_idx_min),
    dt_idx_max = min(ROI1$dt_idx_max, ROI2$dt_idx_max),
    rt_idx_min = max(ROI1$rt_idx_min, ROI2$rt_idx_min),
    rt_idx_max = min(ROI1$rt_idx_max, ROI2$rt_idx_max)
  )
  area_int <- (ROI_int$dt_idx_max - ROI_int$dt_idx_min)*(ROI_int$rt_idx_max - ROI_int$rt_idx_min)
  stopifnot(area_int >= 0) # area can't be negative
  area_union <- area1 + area2 - area_int
  stopifnot(area_union > 0) # should be positive here
  return(area_int/area_union)
}


