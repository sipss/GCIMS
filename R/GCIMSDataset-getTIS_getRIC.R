#' Extract the Reverse Ion Chromatogram and Total Ion Spectrum from the samples
#'
#' @param object A GCIMSDataset object
#'
#' @return The [GCIMSDataset] object, with the delayed operation to compute the
#' RIC and TIS matrices.
#' @export
extract_RIC_and_TIS <- function(object) {
  delayed_op <- GCIMSDelayedOp(
    name = "extract_RIC_and_TIS",
    fun = NULL,
    fun_extract = function(x) {
      dt <- dtime(x)
      rt <- rtime(x)
      intmat <- intensity(x)
      tis <- rowSums(intmat)
      ric_pos <- which.max(tis)
      ric <- intmat[ric_pos, ]
      ric <- max(ric) - ric
      ric <- ric/sum(ric)
      list(ric = ric, tis = tis, rt = rt, dt = dt)
    },
    fun_aggregate = function(ds, objs) {
      num_samples <- length(objs)
      rics <- purrr::map(objs, "ric")
      tiss <- purrr::map(objs, "tis")
      dtimes <- purrr::map(objs, "dt")
      rtimes <- purrr::map(objs, "rt")

      axes_heterogeneity <- ds@envir$axes_heterogeneity
      if (axes_heterogeneity == "needs_interpolate") {
        # Different steps, interpolate to a common step
        dt_ref <- object@envir$dt_ref
        rt_ref <- object@envir$rt_ref
        ds@envir$TIS <- matrix(NA, nrow = num_samples, ncol = length(dt_ref))
        ds@envir$RIC <- matrix(NA, nrow = num_samples, ncol = length(rt_ref))
        for (i in length(tiss)) {
          dt <- dtimes[[i]]
          rt <- rtimes[[i]]
          ds@envir$RIC[i,] <- signal::interp1(dt, tiss[[i]], dt_ref)
          ds@envir$TIS[i,] <- signal::interp1(rt, rics[[i]], rt_ref)
        }
      } else if (axes_heterogeneity == "needs_cutting") {
        dt_ref <- object@envir$dt_ref
        rt_ref <- object@envir$rt_ref
        ds@envir$TIS <- matrix(NA, nrow = num_samples, ncol = length(dt_ref))
        ds@envir$RIC <- matrix(NA, nrow = num_samples, ncol = length(rt_ref))
        for (i in length(tiss)) {
          dt <- dtimes[[i]]
          rt <- rtimes[[i]]
          dt_start_idx <- which(dt >= dt_ref[1])[1]
          dt_end_idx <- which(dt <= dt_ref[length(dt_ref)])[1]
          rt_start_idx <- which(rt >= rt_ref[1])[1]
          rt_end_idx <- which(rt <= rt_ref[length(rt_ref)])[1]
          ds@envir$RIC[i,] <- rics[[i]][rt_start_idx:rt_end_idx]
          ds@envir$TIS[i,] <- tiss[[i]][dt_start_idx:dt_end_idx]
        }
      } else {
        ds@envir$TIS <- do.call(rbind, tiss)
        ds@envir$RIC <- do.call(rbind, rics)
      }
      stopifnot(nrow(ds@envir$TIS) == num_samples)
      stopifnot(nrow(ds@envir$RIC) == num_samples)
      ds
    }
  )
  object <- appendDelayedOp(object, delayed_op)
  invisible(object)
}

#' @describeIn GCIMSDataset Plot Total Ion Spectra
#'
#' @param object A [GCIMSDataset] object
#'
#' @return A matrix with samples in rows and the drift time in columns
#' @export
setMethod("getTIS", "GCIMSDataset", function(object) {
  if (!hasDelayedOps(object) && !is.null(object@envir$TIS)) {
    return(object@envir$TIS)
  }
  object <- extract_RIC_and_TIS(object)
  object <- realize(object)
  out <- object@envir$TIS
  dimnames(out) <- list(
    SampleID = sampleNames(object),
    drift_time_ms = object@envir$dt_ref
  )
  out
})

#' @describeIn GCIMSDataset Plot Total Ion Spectra
#'
#' @param object A [GCIMSDataset] object
#'
#' @return  The RIC matrix
#' @export
setMethod("getRIC", "GCIMSDataset", function(object) {
  if (!hasDelayedOps(object) && !is.null(object@envir$RIC)) {
    return(object@envir$RIC)
  }
  object <- extract_RIC_and_TIS(object)
  object <- realize(object)
  out <- object@envir$RIC
  dimnames(out) <- list(
    SampleID = sampleNames(object),
    retention_time_s = object@envir$rt_ref
  )
  out
})


#' @describeIn GCIMSDataset Plot Total Ion Spectra
#'
#' @param object A [GCIMSDataset] object
#'
#' @return A plot
#' @export
setMethod(
  "plotTIS",
  "GCIMSDataset",
  function(object) {
    tis <- getTIS(object)
    matplot(
      x = object@envir$dt_ref,
      y = t(tis),
      type = "l",
      xlab = "Drift time (ms)",
      ylab = "Intensity (a.u.)"
    )
  }
)

#' @describeIn GCIMSDataset Plot Reverse Ion Chromatograms
#'
#' @param object A [GCIMSDataset] object
#'
#' @return A plot
#' @export
setMethod(
  "plotRIC",
  "GCIMSDataset",
  function(object) {
    ric <- getRIC(object)
    matplot(
      x = object@envir$rt_ref,
      y = t(ric),
      type = "l",
      xlab = "Retention time (s)",
      ylab = "Intensity (a.u.)"
    )
  }
)
