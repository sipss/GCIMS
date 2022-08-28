#' Peak detection for a GCIMSSample
#' @param object A [GCIMSSample] object
#' @param noise_level     Scalar number. The number of times the standard deviation
#'   above the noise level needed to detect a peak.
#' @param verbose If `TRUE`, debug information will be printed
#' @return The modified [GCIMSSample], with a peak list
#' @export
setMethod(
  "findPeaks",
  "GCIMSSample",
  function(object, noise_level = 3, verbose = FALSE, dt_length_ms = 0.14, rt_length_s = 3, iou_overlap_threshold = 0.2) {
    dt <- dtime(object)
    rt <- rtime(object)
    int_mat <- intensity(object)
    dt_length_pts <- units_to_points(dt_length_ms, dt[2] - dt[1], must_odd = TRUE)
    rt_length_pts <- units_to_points(rt_length_s, rt[2] - rt[1], must_odd = TRUE)
    peak_list <- peak_detection(
      drift_time = dt,
      retention_time = rt,
      int_mat = int_mat,
      noise_level = noise_level,
      verbose = verbose,
      dt_length_pts = dt_length_pts,
      rt_length_pts = rt_length_pts,
      iou_overlap_threshold = iou_overlap_threshold
    )
    peaks(object) <- peak_list
    object
  }
)


#' @describeIn GCIMSSample-methods Get the peak list
#' @importMethodsFrom ProtGenerics peaks
#' @export
setMethod(
  "peaks",
  "GCIMSSample",
  function(object) {
    if (is.null(object@peaks)) {
      stop("Please run findPeaks() on the sample first")
    }
    tibble::as_tibble(cbind(SampleID = object@description, object@peaks))
  }
)

#' @describeIn GCIMSSample-methods Set the peak list
#' @param value A data frame with the peak list
#' @importMethodsFrom ProtGenerics "peaks<-"
#' @export
setMethod(
  "peaks<-",
  "GCIMSSample",
  function(object, value) {
    value <- methods::as(value, "DataFrame")
    # We don't want this column:
    value$SampleID <- NULL
    object@peaks <- value
    object
  }
)
