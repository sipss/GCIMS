#' Get the total ion spectrum
#'
#' @param object A [GCIMSSample] object
#'
#' @return A numeric vector with the total ion spectrum
#' @export
#'
#' @examples
#' sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
#' s <- read_mea(sample_file)
#' tis <- getTIS(s)
setMethod("getTIS", "GCIMSSample", function(object) {
  intmat <- intensity(object)
  tis <- rowSums(intmat)
  tis
})

#' Get the reverse ion chromatogram
#'
#' @param object A [GCIMSSample] object
#'
#' @return A numeric vector with the reverse ion chromatogram
#' @export
#'
#' @examples
#' sample_file <- system.file("extdata", "sample_formats", "small.mea.gz", package = "GCIMS")
#' s <- read_mea(sample_file)
#' ric <- getRIC(s)
setMethod("getRIC", "GCIMSSample", function(object) {
  intmat <- intensity(object)
  ric_pos <- rip_position(object)
  ric <- intmat[ric_pos, ]
  ric <- max(ric) - ric
  ric <- ric/sum(ric)
  ric
})

#' Drift time index of the Reactant Ion Peak (RIP)
#'
#' The RIP is the drift time with the highest total ion signal, summed
#' across all retention times (the `which.max()` of [getTIS()]).
#'
#' @param object A [GCIMSSample] object
#' @return An integer, the drift time index of the RIP
#' @noRd
rip_position <- function(object) {
  which.max(rowSums(intensity(object)))
}

#' Plot the total ion spectrum of a sample
#' @param object A [GCIMSSample] object
#' @inheritParams dt_rt_range_normalization
#' @return A ggplot2 plot object
#' @export
setMethod("plotTIS", "GCIMSSample", function(object, dt_range = NULL) {
  dt <- dtime(object)
  rt <- rtime(object)
  idx <- dt_rt_range_normalization(dt, rt, dt_range = dt_range)
  tis <- getTIS(object)
  spec <- GCIMSSpectrum(
    drift_time = dt[idx[["dt_logical"]]],
    intensity = tis[idx[["dt_logical"]]],
    retention_time_idx = unique(c(idx[["rt_idx_min"]], idx[["rt_idx_max"]])),
    retention_time_s = unique(c(idx[["rt_s_min"]], idx[["rt_s_max"]])),
    description = object@description
  )
  plot(spec) + ggplot2::labs(y = "TIS Intensity (a.u.)")
})

#' Plot the reverse ion chromatogram of a sample
#' @param object A [GCIMSSample] object
#' @inheritParams dt_rt_range_normalization
#' @return A ggplot2 plot object
#' @export
setMethod("plotRIC", "GCIMSSample", function(object, rt_range = NULL) {
  dt <- dtime(object)
  rt <- rtime(object)
  ric <- getRIC(object)
  ric_pos <- rip_position(object)
  idx <- dt_rt_range_normalization(dt, rt, dt_idx = ric_pos, rt_range = rt_range)
  chrom <- GCIMSChromatogram(
    retention_time = rt[idx[["rt_logical"]]],
    intensity = ric[idx[["rt_logical"]]],
    drift_time_idx = unique(c(idx[["dt_idx_min"]], idx[["dt_idx_max"]])),
    drift_time_ms = unique(c(idx[["dt_ms_min"]], idx[["dt_ms_max"]])),
    description = object@description
  )
  plot(chrom) + ggplot2::labs(y = "RIC Intensity (a.u.)")
})

