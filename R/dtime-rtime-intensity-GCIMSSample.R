#' Drift time, Retention time, Intensity of GCIMSSamples
#'
#' Functions to extract the drift time, the retention time and
#' the intensity.
#'
#' @name GCIMSSample-rtime-dtime-intensity
#' @param object A [GCIMSSample] object
NULL

#' @describeIn GCIMSSample-rtime-dtime-intensity Get the drift time vector
#' @return The drift time of the sample
#' @export
setMethod("dtime", "GCIMSSample", function(object) object@drift_time)

#' @describeIn GCIMSSample-rtime-dtime-intensity Get the retention time vector
#' @importMethodsFrom ProtGenerics rtime
#' @return The retention time of the sample
#' @export
setMethod("rtime", "GCIMSSample", function(object) object@retention_time)

#' @describeIn GCIMSSample-rtime-dtime-intensity Get the intensity matrix
#' @inheritParams dt_rt_range_normalization
#' @export
#' @examples
#' mea_file <- system.file("extdata",  "sample_formats", "small.mea.gz", package = "GCIMS")
#' gcims_sample <- read_mea(mea_file)
#' my_matrix <- intensity(gcims_sample, dt_range = c(7, 8), rt_range = c(1,30))
setMethod("intensity", "GCIMSSample", function(object, dt_range = NULL, rt_range = NULL, dt_idx = NULL, rt_idx = NULL) {
  dt <- dtime(object)
  rt <- rtime(object)
  if (inherits(dt_range, "dt_rt_range_normalization")) {
    idx <- dt_range
  } else {
    idx <- dt_rt_range_normalization(dt, rt, dt_range, rt_range, dt_idx, rt_idx)
  }
  dt_idx <- idx[["dt_logical"]]
  rt_idx <- idx[["rt_logical"]]
  out <- object@data[dt_idx, rt_idx, drop = FALSE]
  dimnames(out) <- list(dt_ms = dt[dt_idx], rt_s = rt[rt_idx])
  out
})


#' @describeIn GCIMSSample-rtime-dtime-intensity Set the intensity matrix
#' @importMethodsFrom ProtGenerics "intensity<-"
#' @param object A GCIMSSample object.
#' @param value A matrix of dimensions `dim(object)`
#' @export
#' @examples
#' mea_file <- system.file("extdata",  "sample_formats", "small.mea.gz", package = "GCIMS")
#' gcims_sample <- read_mea(mea_file)
#' my_matrix <- intensity(gcims_sample)
#' intensity(gcims_sample) <- my_matrix/100
setMethod("intensity<-", "GCIMSSample", function(object, value) {
  if (!is.matrix(value)) {
    cli_abort("The given intensity should be a matrix")
  }
  if (!identical(dim(object), dim(value))) {
    cli_abort("The given intensity should have dimensions {dim(object)} and it has dimensions {dim(value)}")
  }
  object@data <- value
  object
})
