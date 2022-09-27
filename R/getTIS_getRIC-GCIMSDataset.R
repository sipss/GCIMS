#' @describeIn GCIMSDataset Get Total Ion Spectra matrix
#'
#' @param object A [GCIMSDataset] object
#'
#' @return A matrix with samples in rows and the drift time in columns
#' @export
setMethod("getTIS", "GCIMSDataset", function(object) {
  if (hasDelayedOps(object) || is.null(object@envir$TIS)) {
    object <- extract_RIC_and_TIS(object)
    object <- realize(object)
  }
  out <- object@envir$TIS
  dimnames(out) <- list(
    SampleID = sampleNames(object),
    drift_time_ms = object@envir$dt_ref
  )
  out
})

#' Get Reverse Ion Chromatogram
#'
#' @param object A [GCIMSDataset] object
#'
#' @return  The RIC matrix
#' @export
setMethod("getRIC", "GCIMSDataset", function(object) {
  if (hasDelayedOps(object) || is.null(object@envir$RIC)) {
    object <- extract_RIC_and_TIS(object)
    object <- realize(object)
  }
  out <- object@envir$RIC
  dimnames(out) <- list(
    SampleID = sampleNames(object),
    retention_time_s = object@envir$rt_ref
  )
  out
})


#' Plot Total Ion Spectra
#'
#' @param object A [GCIMSDataset] object
#' @inheritParams dt_rt_range_normalization
#' @param sample A number or a string with the sample index or name. If `NULL`, all samples are returned
#' @return The plot of the TIS
#' @export
setMethod(
  "plotTIS",
  "GCIMSDataset",
  function(object, dt_range = NULL, sample = NULL) {
    tis <- getTIS(object)
    dt <- dtime(object)
    sample_names <- sampleNames(object)
    if (is.null(sample)) {
      sample <- sample_names
    }
    sample_idx <- sample_name_or_number_to_both(sample, sample_names)
    idx <- dt_rt_range_normalization(dt = dt, dt_range = dt_range)
    tis_long <- reshape2::melt(tis[sample_idx$idx, idx$dt_logical, drop = FALSE], value.name = "TIS")
    gplt <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = tis_long,
        mapping = ggplot2::aes(
          x = .data$drift_time_ms,
          y = .data$TIS,
          color = .data$SampleID
        )
      ) +
      ggplot2::labs(
        x = "Drift time (ms)",
        y = "TIS Intensity (a.u.)",
        color = "SampleID"
      )
    gplt
  }
)

#' @describeIn GCIMSDataset Plot Reverse Ion Chromatograms
#'
#' @param object A [GCIMSDataset] object
#' @inheritParams dt_rt_range_normalization
#' @param sample A number or a string with the sample index or name. If `NULL`, all samples are returned
#' @return A plot
#' @export
setMethod(
  "plotRIC",
  "GCIMSDataset",
  function(object, rt_range = NULL, sample = NULL) {
    ric <- getRIC(object)
    rt <- rtime(object)
    sample_names <- sampleNames(object)
    if (is.null(sample)) {
      sample <- sample_names
    }
    sample_idx <- sample_name_or_number_to_both(sample, sample_names)
    idx <- dt_rt_range_normalization(rt = rt, rt_range = rt_range)
    ric_long <- reshape2::melt(ric[sample_idx$idx, idx$rt_logical, drop = FALSE], value.name = "RIC")

    gplt <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = ric_long,
        mapping = ggplot2::aes(
          x = .data$retention_time_s,
          y = .data$RIC,
          color = .data$SampleID
        )
      ) +
      ggplot2::labs(
        x = "Retention time (s)",
        y = "RIC Intensity (a.u.)",
        color = "SampleID"
      )
    gplt
  }
)


.extract_RIC_and_TIS_fun_extract <- function(x) {
  dt <- dtime(x)
  rt <- rtime(x)
  intmat <- intensity(x)
  tis <- rowSums(intmat)
  ric_pos <- which.max(tis)
  ric <- intmat[ric_pos, ]
  ric <- max(ric) - ric
  ric <- ric/sum(ric)
  list(ric = ric, tis = tis, rt = rt, dt = dt)
}

.extract_RIC_and_TIS_fun_aggregate <- function(ds, objs) {
  num_samples <- length(objs)
  rics <- purrr::map(objs, "ric")
  tiss <- purrr::map(objs, "tis")
  dtimes <- purrr::map(objs, "dt")
  rtimes <- purrr::map(objs, "rt")

  axes_heterogeneity <- axesHeterogeneity(ds)
  if (is.null(axes_heterogeneity)) {
    abort(c("UnexpectedError", "extract_dtime_rtime() should have been called before extract_RIC_and_TIS()"))
  }
  if (axes_heterogeneity == "needs_interpolate") {
    # Different steps, interpolate to a common step
    dt_ref <- ds@envir$dt_ref
    rt_ref <- ds@envir$rt_ref
    ds@envir$TIS <- matrix(NA, nrow = num_samples, ncol = length(dt_ref))
    ds@envir$RIC <- matrix(NA, nrow = num_samples, ncol = length(rt_ref))
    for (i in length(tiss)) {
      dt <- dtimes[[i]]
      rt <- rtimes[[i]]
      ds@envir$RIC[i,] <- signal::interp1(dt, tiss[[i]], dt_ref)
      ds@envir$TIS[i,] <- signal::interp1(rt, rics[[i]], rt_ref)
    }
  } else if (axes_heterogeneity == "needs_cutting") {
    dt_ref <- ds@envir$dt_ref
    rt_ref <- ds@envir$rt_ref
    ds@envir$TIS <- matrix(NA, nrow = num_samples, ncol = length(dt_ref))
    ds@envir$RIC <- matrix(NA, nrow = num_samples, ncol = length(rt_ref))
    for (i in seq_len(num_samples)) {
      dt <- dtimes[[i]]
      rt <- rtimes[[i]]
      dt_start_idx <- utils::head(which(dt >= dt_ref[1L]), 1L)
      dt_end_idx <- utils::tail(which(dt <= dt_ref[length(dt_ref)]), 1L)

      rt_start_idx <- utils::head(which(rt >= rt_ref[1L]), 1L)
      rt_end_idx <- utils::tail(which(rt <= rt_ref[length(rt_ref)]), 1L)
      if (any(
        length(dt_start_idx) == 0,
        length(dt_end_idx) == 0,
        length(rt_start_idx) == 0,
        length(rt_end_idx) == 0
      )) {
        abort("Unexpected error: No overlap between sample retention and drift times.")
      }
      this_ric <- rics[[i]][rt_start_idx:rt_end_idx]
      if (length(this_ric) != ncol(ds@envir$RIC)) {
        abort("Unexpected length of RIC")
      }
      ds@envir$RIC[i,] <- this_ric
      this_tiss <- tiss[[i]][dt_start_idx:dt_end_idx]
      if (length(this_tiss) != ncol(ds@envir$TIS)) {
        abort("Unexpected length of TIS")
      }
      ds@envir$TIS[i,] <- this_tiss
    }
  } else {
    ds@envir$TIS <- do.call(rbind, tiss)
    ds@envir$RIC <- do.call(rbind, rics)
  }
  stopifnot(nrow(ds@envir$TIS) == num_samples)
  stopifnot(nrow(ds@envir$RIC) == num_samples)
  ds
}

#' Extract the Reverse Ion Chromatogram and Total Ion Spectrum from the samples
#'
#' @param object A GCIMSDataset object
#'
#' @return The [GCIMSDataset] object, with the delayed operation to compute the
#' RIC and TIS matrices.
#' @noRd
extract_RIC_and_TIS <- function(object) {
  object <- extract_dtime_rtime(object)
  delayed_op <- GCIMSDelayedOp(
    name = "extract_RIC_and_TIS",
    fun = NULL,
    fun_extract = .extract_RIC_and_TIS_fun_extract,
    fun_aggregate = .extract_RIC_and_TIS_fun_aggregate
  )
  object <- appendDelayedOp(object, delayed_op)
  invisible(object)
}

